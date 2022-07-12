#include <iomanip>
#include <TFile.h>
#include <TTree.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <interface_main/SQHitVector_v1.h>
#include <interface_main/SQTrackVector_v1.h>
#include <interface_main/SQDimuonVector_v1.h>

#include <string>
#include <vector>

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
	int nTracklets = trackletVec->size();
	for(int i = 0; i < nTracklets; ++i)
	{
		Tracklet* tracklet = trackletVec->at(i);
		nHits = tracklet->getNHits();
		chisq = tracklet->getChisq();
		
		//very loose cuts here
		if(nHits < 9) continue;
		if(chisq > 10.) continue;
		
		effi_h4(tracklet);
		hodo24(tracklet);
		hodo42(tracklet);
		
		saveTree->Fill();
		
		detectorID.clear();
		elementID_exp.clear();
		elementID_closest.clear();
		ele24.clear();
		ele42.clear();
	}
	
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
	saveTree->Branch("ele42", &ele42);
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

int AnaModule::acc_plane(Tracklet* tracklet, std::vector<std::string> &vec)
{
	int nvec = vec.size();
	//int nhits = hitVector->size();
	std::vector<int> acc_mask;

	for(int i = 0; i < nvec; i++)
	{
		int id_vec = p_geomSvc->getDetectorID(vec.at(i));
		double z_exp = p_geomSvc->getPlanePosition(id_vec);
		double x_exp = tracklet->getExpPositionX(z_exp);
		double y_exp = tracklet->getExpPositionY(z_exp);
		
		if(!p_geomSvc->isInPlane(id_vec, x_exp, y_exp)) continue;
		
		int id_acc = p_geomSvc->getExpElementID(id_vec, tracklet->getExpPositionW(id_vec));
		
		if(id_acc < 1 || id_acc > p_geomSvc->getPlaneNElements(id_vec)) continue;

		SQHit* hit = findHit(id_vec, id_acc);
		int clo_id = hit == nullptr ? -1 : hit->get_element_id();
		if(clo_id < 0) continue;
		std::cout <<"vec_id : " << id_vec <<" acc_id : " << id_acc << " clo_id : " << clo_id << std::endl;
		if(abs(id_acc - clo_id) < 2){acc_mask.push_back(clo_id);}
	}

	/*for(int i = 0; i < nhits; i++)
	{
		for(int j = 0; j < nvec; j++)
		{
			SQHit* hit = hitVector->at(i);
			if(hit->get_detector_id() == vec.at(j)){acc_mask.push_back(hit->get_detector_id());}
		}
	}*/

	int mask_hits = acc_mask.size();
	acc_mask.clear();

	return mask_hits;
}

// use only good tracks in the denominator
int AnaModule::acc_h4(Tracklet* tracklet, int id)
{
	int nacc = -1;

	// H4Y1L -> H3B, H3T, H4B, H4T
	if(id == p_geomSvc->getDetectorID("H4Y1L"))
	{
		std::vector<std::string> acc41 = {"H3B", "H3T", "H4B", "H4T"};
		nacc = acc_plane(tracklet, acc41);
		std::cout << "det_id : " << id << " nacc : " << nacc << std::endl;
	}

	if(id == p_geomSvc->getDetectorID("H4Y1R"))
	{
		std::vector<std::string> acc42 = {"H3B", "H3T", "H4B", "H4T"};
		nacc = acc_plane(tracklet, acc42);
		std::cout << "det_id : " << id << " nacc : " << nacc << std::endl;
	}

	if(id == p_geomSvc->getDetectorID("H4Y2L"))
	{
		std::vector<std::string> acc43 = {"H3B", "H3T", "H4B", "H4T"};
		nacc = acc_plane(tracklet, acc43);
		std::cout << "det_id : " << id << " nacc : " << nacc << std::endl;
	}

	if(id == p_geomSvc->getDetectorID("H4Y2R"))
	{
		std::vector<std::string> acc44 = {"H3B", "H3T", "H4B", "H4T"};
		nacc = acc_plane(tracklet, acc44);
		std::cout << "det_id : " << id << " nacc : " << nacc << std::endl;
	}

	// H4B -> P1Y1, P1Y2, P2Y1, P2Y2
	if(id == p_geomSvc->getDetectorID("H4B"))
	{
		std::vector<std::string> acc45 = {"P1Y1", "P1Y2", "P2Y1", "P2Y2"};
		nacc = acc_plane(tracklet, acc45);
		std::cout << "det_id : " << id << " nacc : " << nacc << std::endl;
	}

	if(id ==  p_geomSvc->getDetectorID("H4T"))
	{
		std::vector<std::string> acc46 = {"P1Y1", "P1Y2", "P2Y1", "P2Y2"};
		nacc = acc_plane(tracklet, acc46);
		std::cout << "det_id : " << id << " nacc : " << nacc << std::endl;
	}

	return nacc;
}


void AnaModule::effi_h4(Tracklet* tracklet)
{
  std::vector<std::string> hodo4 = {"H4Y1L", "H4Y1R", "H4Y2L", "H4Y2R", "H4B", "H4T"};
  int nhodo = hodo4.size();
		
	for(int j = 0; j < nhodo; j++)
	{
		// get only acctepted NIM4 || MATRIX5 events
		// beam like and reverse beam like
		if(!(event->get_trigger(SQEvent::NIM4)|| event->get_trigger(SQEvent::MATRIX5))) continue;
		//if(!event->get_trigger(SQEvent::NIM4)) continue;
		//if(!event->get_trigger(SQEvent::MATRIX5)) continue;

		int det_id = p_geomSvc->getDetectorID(hodo4.at(j));

		int nacc_hits = acc_h4(tracklet, det_id);
		if(nacc_hits != 2) continue;

		std::cout << "*** accepted track ***" << std::endl;
		
		//std::cout << "det_id : " << det_id << std::endl;
		//int exp_id = fit_prop(det_id, tracklet);
		//int exp_id = fit2d_prop(det_id, tracklet);
		

		double z_exp = p_geomSvc->getPlanePosition(det_id);
		double x_exp = tracklet->getExpPositionX(z_exp);
		double y_exp = tracklet->getExpPositionY(z_exp);
		
		if(!p_geomSvc->isInPlane(det_id, x_exp, y_exp)) continue;


		int exp_id = p_geomSvc->getExpElementID(det_id, tracklet->getExpPositionW(det_id));

		if(exp_id < 1 || exp_id > p_geomSvc->getPlaneNElements(det_id)) continue;
		
		SQHit* hit = findHit(det_id, exp_id);
		int close_id = hit == nullptr ? -1 : hit->get_element_id();
		detectorID.push_back(det_id);
		elementID_exp.push_back(exp_id);
		elementID_closest.push_back(close_id);
	}
}


void AnaModule::hodo42(Tracklet* tracklet)
{
	std::vector<std::string> vec42 = {"H2L", "H2R", "H4Y2L", "H4Y2R"};
	int nhodo42 = vec42.size();
		
	for(int j = 0; j <  nhodo42; j++)
	{
		// get only acctepted NIM4 || MATRIX5 events
		// beam like and reverse beam like
		//if(!(event->get_trigger(SQEvent::NIM4)|| event->get_trigger(SQEvent::MATRIX5))) continue;
		//if(!event->get_trigger(SQEvent::NIM4)) continue;
		if(!event->get_trigger(SQEvent::MATRIX5)) continue;

		int hodoid = p_geomSvc->getDetectorID(vec42.at(j));
		
		//int nacc_hits = acc_h4(tracklet, hodoid);
		//if(nacc_hits != 2) continue;


		double z_exp = p_geomSvc->getPlanePosition(hodoid);
		double x_exp = tracklet->getExpPositionX(z_exp);
		double y_exp = tracklet->getExpPositionY(z_exp);
		
		if(!p_geomSvc->isInPlane(hodoid, x_exp, y_exp)) continue;
		
		int id42 = p_geomSvc->getExpElementID(hodoid, tracklet->getExpPositionW(hodoid));

		if(id42 < 1 || id42 > p_geomSvc->getPlaneNElements(hodoid)) continue;

		ele42.push_back(id42);
	}
}


void AnaModule::hodo24(Tracklet* tracklet)
{
	std::vector<std::string> vec24 = {"H2L", "H2R", "H4Y2L", "H4Y2R"};
	int nhodo24 = vec24.size();
		
	for(int j = 0; j <  nhodo24; j++)
	{
		// get only acctepted NIM4 || MATRIX5 events
		// beam like and reverse beam like
		//if(!(event->get_trigger(SQEvent::NIM4)|| event->get_trigger(SQEvent::MATRIX5))) continue;
		if(!event->get_trigger(SQEvent::NIM4)) continue;
		//if(!event->get_trigger(SQEvent::MATRIX5)) continue;
	
		int hodoid = p_geomSvc->getDetectorID(vec24.at(j));
		
		//int nacc_hits = acc_h4(tracklet, hodoid);
		//if(nacc_hits != 2) continue;

		double z_exp = p_geomSvc->getPlanePosition(hodoid);
		double x_exp = tracklet->getExpPositionX(z_exp);
		double y_exp = tracklet->getExpPositionY(z_exp);

		if(!p_geomSvc->isInPlane(hodoid, x_exp, y_exp)) continue;

		int id24 = p_geomSvc->getExpElementID(hodoid, tracklet->getExpPositionW(hodoid));

		if(id24 < 1 || id24 > p_geomSvc->getPlaneNElements(hodoid)) continue;

		ele24.push_back(id24);
	}
}
