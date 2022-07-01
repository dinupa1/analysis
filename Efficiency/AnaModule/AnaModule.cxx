#include <iomanip>
#include <TFile.h>
#include <TTree.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <interface_main/SQHitVector_v1.h>
#include <interface_main/SQTrackVector_v1.h>
#include <interface_main/SQDimuonVector_v1.h>


#include "AnaModule.h"

AnaModule::AnaModule(const std::string& name): SubsysReco(name), p_geomSvc(GeomSvc::instance())
{}

AnaModule::~AnaModule()
{}

int AnaModule::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  eventID = 0;
  MakeTree();
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::process_event(PHCompositeNode* topNode)
{
	fill_h4();
	
	saveTree->Fill();
	
	detectorID.clear();
	elementID_exp.clear();
	elementID_closest.clear();

  ++eventID;
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::End(PHCompositeNode* topNode)
{
  saveFile->cd();
  saveTree->Write();
  saveFile->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaModule::GetNodes(PHCompositeNode* topNode)
{
  event = findNode::getClass<SQEvent>(topNode, "SQEvent");
  hitVector   = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  trackletVec = findNode::getClass<TrackletVector>(topNode, "TrackletVector");
  if(!event || !hitVector || !trackletVec)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void AnaModule::MakeTree()
{
  saveFile = new TFile(saveName, "RECREATE");

  saveTree = new TTree("save", "Efficiency tree Created by AnaModule");
  saveTree->Branch("eventID", &eventID, "eventID/I");
  saveTree->Branch("detectorID", &detectorID);
  saveTree->Branch("elementID_exp", &elementID_exp);
  saveTree->Branch("elementID_closest", &elementID_closest);
  saveTree->Branch("nHits", &nHits, "nHits/I");
  saveTree->Branch("chisq", &chisq, "chisq/D");
	//for debugging only
	saveTree->Branch("ele24", &ele24);
}

SQHit* AnaModule::findHit(int detID, int eleID)
{
  int delta = 999;
  SQHit* hit = nullptr;
  for(SQHitVector::Iter it = hitVector->begin(); it != hitVector->end(); ++it)
  {
    if((*it)->get_detector_id() != detID) continue;
    int delta_temp = abs(eleID - (*it)->get_element_id());
    if(delta > delta_temp)
    {
      delta = delta_temp;
      hit = (*it);
    }
  }

  return hit;
}

int AnaModule::fit_prop(int det_id, Tracklet* tracklet)
{
  std::vector<int> track3 = {19, 21, 22, 51, 52, 53, 54};

  TGraphErrors* gx = new TGraphErrors();
  TGraphErrors* gy = new TGraphErrors();

  int ndet = track3.size();

  
  double zz0;
  double xx0;
  double yy0;
  double exx0;
  double eyy0;

  for(int i = 0; i < ndet; i++)
  {
    zz0 = p_geomSvc->getPlanePosition(track3.at(i));
    xx0 = tracklet->getExpPositionX(zz0);
    yy0 = tracklet->getExpPositionY(zz0);
    exx0 = tracklet->getExpPosErrorX(zz0);
    eyy0 = tracklet->getExpPosErrorY(zz0);
    
    if(i > 50)
    {
      int nhits = hitVector->size();
      
      for(int j = 0; j < nhits; j++)
      {
        
        SQHit* hit = hitVector->at(j);
        
        if(hit->get_detector_id() == i)
        {
          zz0 = hit->get_truth_z();
          xx0 = hit->get_truth_y();
          yy0 = hit->get_truth_x();
          exx0 = 0.;
          eyy0 = 0.;
        }
      
      }
    }
    
    // set x points
    gx->SetPoint(i, zz0, xx0);
    gx->SetPointError(i, 0., exx0);

    // set y points
    gy->SetPoint(i, zz0, yy0);
    gy->SetPointError(i, 0., eyy0);

    //std::cout << "det : " << track3.at(i) << " x : " << xx0 << " y : " << yy0 << " z : " << zz0 << " ex :" << exx0 << " ey : " << eyy0 << std::endl;
  }

  // fit functions
  TF1* fx = new TF1("fx", "[0]* x + [1]", 1900., 2400.);
  TF1* fy = new TF1("fy", "[0]* x + [1]", 1900., 2400.);

  gx->Fit("fx");
  gy->Fit("fy");

  double axx = fx->GetParameter(0);
  double cxx = fx->GetParameter(1);

  double ayy = fy->GetParameter(0);
  double cyy = fy->GetParameter(1);

  double zz1 = p_geomSvc->getPlanePosition(det_id);
  double xx1 = axx* zz1 + cxx;
  double yy1 = ayy* zz1 + cyy;

  /*std::cout << "***     ***" << std::endl;
  std::cout << "hodo : " << det_id << " x : " << xx1 << " y :" << yy1 << " z : " << zz1 << std::endl;
  std::cout << "***     ***" << std::endl;*/

  if(p_geomSvc->isInPlane(det_id, xx1, yy1))
  {
    double pos = p_geomSvc->getCostheta(det_id)*xx1 + p_geomSvc->getSintheta(det_id)*yy1;
    return p_geomSvc->getExpElementID(det_id, pos);
  }
  
  return -1;
}

int AnaModule::fit2d_prop(int det_id, Tracklet* tracklet)
{
  std::vector<int> track3 = {19, 21, 22, 51, 52, 53, 54};
  
  TGraph2DErrors* gg = new TGraph2DErrors();

  int ndet = track3.size();
  
  double zz0;
  double xx0;
  double yy0;
  double exx0;
  double eyy0;

  for(int i = 0; i < ndet; i++)
  {
    zz0 = p_geomSvc->getPlanePosition(track3.at(i));
    xx0 = tracklet->getExpPositionX(zz0);
    yy0 = tracklet->getExpPositionY(zz0);
    exx0 = tracklet->getExpPosErrorX(zz0);
    eyy0 = tracklet->getExpPosErrorY(zz0);
    
    if(i > 50)
    {
      int nhits = hitVector->size();
      
      for(int j = 0; j < nhits; j++)
      {
        
        SQHit* hit = hitVector->at(j);
        
        if(hit->get_detector_id() == i)
        {
          zz0 = hit->get_truth_z();
          xx0 = hit->get_truth_y();
          yy0 = hit->get_truth_x();
          exx0 = 0.;
          eyy0 = 0.;
        }
      
      }
    }
    
    // set points
    gg->SetPoint(i, zz0, xx0, yy0);
    gg->SetPointError(i, 0.0, exx0, eyy0);

    //std::cout << "det : " << track3.at(i) << " x : " << xx0 << " y : " << yy0 << " z : " << zz0 << " ex :" << exx0 << " ey : " << eyy0 << std::endl;
  }
  
  // fit function
  
  TF2* ff = new TF2("ff", "[0]* x +[1]* y + [2]");
  
  gg->Fit("ff");
  
  double axx = ff->GetParameter(0);
  double ayy = ff->GetParameter(1);
  
  double zz1 = p_geomSvc->getPlanePosition(det_id);
  double zzp = p_geomSvc->getPlanePosition(21);
  double xxp = tracklet->getExpPositionX(zzp);
  double yyp = tracklet->getExpPositionX(zzp);
  
  double xx1 = xxp + axx* (zz1 - zzp);
  double yy1 = yyp + ayy* (zz1 - zzp);
  

  if(p_geomSvc->isInPlane(det_id, xx1, yy1))
  {
    double pos = p_geomSvc->getCostheta(det_id)*xx1 + p_geomSvc->getSintheta(det_id)*yy1;
    return p_geomSvc->getExpElementID(det_id, pos);
  }
  
  return -1;
}


// accept event for H4 : an event should have hits in all 4 hodo. planes
bool AnaModule::acc_h4(Tracklet* tracklet)
{
  int hodoId[6] = {41, 42, 43, 44, 45, 46};
  
  std::vector<int> acc_mask;
  
  int nhits = hitVector->size();
  
  for(int i = 0; i < nhits; i++)
  {
    for(int j = 0; j < 6; j++)
    {
      SQHit* hit = hitVector->at(i);
      
      //if(!(tracklet->get_track_id() == hit->get_track_id())) continue;
      
      if(hit->get_detector_id() == hodoId[j])
      {
        acc_mask.push_back(hit->get_detector_id());
      }
    }
  }
  
  int mask_hits = acc_mask.size();
  
  if(mask_hits != 6){return false;}
  return true;
}


void AnaModule::effi_h4(Tracklet* tracklet)
{
  // only NIM4 events are considered
  std::vector<int> hodo4 = {41, 42, 43, 44, 45, 46};
  int nhodo = hodo4.size();
  for(int i = 0; i < nhodo; i++)
  {
    int det_id = hodo4.at(i);
    int exp_id = fit_prop(det_id, tracklet);
		//int exp_id = fit2d_prop(det_id, tracklet);
    
    //double z_exp = p_geomSvc->getPlanePosition(det_id);
    //double x_exp = tracklet->getExpPositionX(z_exp);
    //double y_exp = tracklet->getExpPositionY(z_exp);
    
    //if(!p_geomSvc->isInPlane(det_id, x_exp, y_exp)) continue;
    
    //int exp_id = p_geomSvc->getExpElementID(det_id, tracklet->getExpPositionW(det_id));
    
    if(exp_id < 1 || exp_id > p_geomSvc->getPlaneNElements(det_id)) continue;
    
    SQHit* hit = findHit(det_id, exp_id);
    int close_id = hit == nullptr ? -1 : hit->get_element_id();
    
    detectorID.push_back(det_id);
    elementID_exp.push_back(exp_id);
    elementID_closest.push_back(close_id);
  }
}

void AnaModule::hodo24(Tracklet* tracklet)
{
	std::vector<int> vec24 = {37, 46};
	int nhodo24 = vec24.size();
	for(int i = 0; i <  nhodo24; i++)
	{
		int hodoid = vec24.at(i);
		double z_exp = p_geomSvc->getPlanePosition(hodoid);
		double x_exp = tracklet->getExpPositionX(z_exp);
		double y_exp = tracklet->getExpPositionY(z_exp);
		if(!p_geomSvc->isInPlane(hodoid, x_exp, y_exp)) continue;
		int id24 = p_geomSvc->getExpElementID(hodoid, tracklet->getExpPositionW(hodoid));
		ele24.push_back(id24);
	}
}

void AnaModule::fill_h4()
{
	int nTracklets = trackletVec->size();
	for(int i = 0; i < nTracklets; ++i)
	{
		Tracklet* tracklet = trackletVec->at(i);
    nHits = tracklet->getNHits();
    chisq = tracklet->getChisq();
    
    // get only acctepted NIM4 || MATRIX5 events
    // beam like and reverse beam like
    //if(!(event->get_trigger(SQEvent::NIM4)|| event->get_trigger(SQEvent::MATRIX5))) continue;
    //if(!event->get_trigger(SQEvent::NIM4)) continue;
		if(!event->get_trigger(SQEvent::MATRIX5)) continue;
    if(!acc_h4(tracklet)) continue;

    //very loose cuts here
    if(nHits < 9) continue;
    if(chisq > 10.) continue;

    effi_h4(tracklet);
		hodo24(tracklet);
  }
}
