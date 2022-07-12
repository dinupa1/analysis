#ifndef _ANA_Module__H_
#define _ANA_Module__H_

#include <map>
#include <set>
#include <fun4all/SubsysReco.h>
#include <TString.h>
#include <TVector3.h>
#include <interface_main/SQEvent.h>
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TF1.h>
#include <TF2.h>
#include <ktracker/SRecEvent.h>
#include <ktracker/FastTracklet.h>
#include <geom_svc/GeomSvc.h>
#include <interface_main/SQHit_v1.h>
#include <vector>
#include <string>

class TFile;
class TTree;
class SQHitVector;
class SQTrackVector;
class SQDimuonVector;

class AnaModule: public SubsysReco 
{
public:
  AnaModule(const std::string& name = "AnaModule");
  virtual ~AnaModule();

  int Init(PHCompositeNode* topNode);
  int InitRun(PHCompositeNode* topNode);
  int process_event(PHCompositeNode* topNode);
  int End(PHCompositeNode* topNode);

  void set_output_filename(const TString& n) { saveName = n; }

private:
  int GetNodes(PHCompositeNode* topNode);
  void MakeTree();
  int fit_prop(int det_id, Tracklet* tracklet);
  int fit2d_prop(int det_id, Tracklet* tracklet);
  void effi_h4(Tracklet* tracklet);
  int acc_h4(Tracklet* tracklet, int id);
	int acc_plane(Tracklet* tracklet, std::vector<int> &vec);
//	void fill_h4();
	void hodo24(Tracklet* tracklet);
	void hodo42(Tracklet* tracklet);

  SQHit* findHit(int detectorID, int elementID);
  std::set<int> detectorIDs;

  GeomSvc* p_geomSvc;

  // Input
  SQEvent* event;
  SQHitVector*    hitVector;
  TrackletVector* trackletVec;

  // Output
  TString saveName;
  TFile*  saveFile;
  TTree*  saveTree;

  int eventID;
  std::vector<int> detectorID;
  std::vector<int> elementID_exp;
  std::vector<int> elementID_closest;
  int nHits;
  double chisq;

	// fro debugging
	std::vector<int> ele24;
	std::vector<int> ele42;
};

#endif
