/*
 *  embedAnalMaker.h
 *  2012Embed
 *
 *  Created by Keith Landry on 2/3/16.
 *  UCLA
 *
 */


#ifndef EMBEDANALMAKER_H
#define EMBEDANALMAKER_H

//root
#include "TTree.h"
#include "TFile.h"
#include "TH3D.h"
#include "TH1D.h"
#include "TParticle.h"

#include "StMaker.h"

//StMuDstMaker
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
//#include "StMuDSTMaker/COMMON/StMuEmcCollection.h"
//#include "StMuDSTMaker/COMMON/StMuEmcUtil.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"

/*//StEmc
#include "StEmcClusterCollection.h"
#include "StEmcPoint.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/others/emcDetectorName.h"
#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "StEmcRawMaker/StEmcRawMaker.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcUtil/database/StBemcTables.h"
#include "StEvent/StEmcCollection.h"
 */

//StSpinDb
#include "StSpinPool/StSpinDbMaker/StSpinDbMaker.h"

//Trigger
#include "StMuDSTMaker/COMMON/StMuTriggerIdCollection.h"
#include "StEvent/StTriggerId.h"
#include "StTriggerUtilities/StTriggerSimuMaker.h"

//MC
//#include "StMcEventMaker/StMcEventMaker.h"
//#include "StMcEventTypes.hh"

//Pythia
#include "StSpinPool/StJetSkimEvent/StPythiaEvent.h"


#include <string>
#include <sstream>
#include <vector>

//
struct binInfo{
  //Int_t nx;
  Int_t ny;
  Int_t nz;
  
  //Double_t* xlims;
  Double_t* ylims;
  Double_t* zlims;
};



class embedAnalMaker : public StMaker
{
private:
  
  StMuDstMaker*  mMuDstMaker;
  TFile*         mOutputFile;
  string         mOutFileName;
  
  TTree*         mPythTree;
  StPythiaEvent* mPythEvent;
  StTriggerSimuMaker* mTrgSimMkr;

  
  //bin information//
  double mBinWeight;
  double* mPtBinLiimits;
  double* mMBinLiimits;
  double* mEtaBinLiimits;
  
  int nBinsPt;
  int nBinsM;
  int nBinsEta;
  
  bool mBinLimits1Set;
  bool mBinLimits2Set;
  bool mBinLimits3Set;
  bool mPythiaTreeSet;
  
  binInfo mBinsPtM;
  binInfo mBinsPtEta;
  binInfo mBinsMEta;
    
  //Histograms//
  std::vector<std::vector<TH1D*> > mH_ProcIds_PtM;
  std::vector<std::vector<TH1D*> > mH_ProcIds_MEta;
  std::vector<std::vector<TH1D*> > mH_ProcIds_PtEta;
  std::vector<std::vector<TH1D*> > mH_X1_PtM;
  std::vector<std::vector<TH1D*> > mH_X1_MEta;
  std::vector<std::vector<TH1D*> > mH_X1_PtEta;
  std::vector<std::vector<TH1D*> > mH_X2_PtM;
  std::vector<std::vector<TH1D*> > mH_X2_MEta;
  std::vector<std::vector<TH1D*> > mH_X2_PtEta;
  std::vector<std::vector<TH1D*> > mH_Z_PtM;
  std::vector<std::vector<TH1D*> > mH_Z_MEta;
  std::vector<std::vector<TH1D*> > mH_Z_PtEta;  
  std::vector<std::vector<TH1D*> > mH_matchedPdgCode_PtM;
  std::vector<std::vector<TH1D*> > mH_matchedPdgCode_MEta;
  std::vector<std::vector<TH1D*> > mH_matchedPdgCode_PtEta;
  
  
  std::vector<std::vector<TH1D*> > mH_ProcIds_PtM_pyth;
  std::vector<std::vector<TH1D*> > mH_ProcIds_MEta_pyth;
  std::vector<std::vector<TH1D*> > mH_ProcIds_PtEta_pyth;
  std::vector<std::vector<TH1D*> > mH_X1_PtM_pyth;
  std::vector<std::vector<TH1D*> > mH_X1_MEta_pyth;
  std::vector<std::vector<TH1D*> > mH_X1_PtEta_pyth;
  std::vector<std::vector<TH1D*> > mH_X2_PtM_pyth;
  std::vector<std::vector<TH1D*> > mH_X2_MEta_pyth;
  std::vector<std::vector<TH1D*> > mH_X2_PtEta_pyth;
  std::vector<std::vector<TH1D*> > mH_matchedPdgCode_PtM_pyth;
  std::vector<std::vector<TH1D*> > mH_matchedPdgCode_MEta_pyth;
  std::vector<std::vector<TH1D*> > mH_matchedPdgCode_PtEta_pyth;
  std::vector<std::vector<TH1D*> > mH_Z_PtM_pyth;
  std::vector<std::vector<TH1D*> > mH_Z_MEta_pyth;
  std::vector<std::vector<TH1D*> > mH_Z_PtEta_pyth;  
  TH1D* mMatchedPdgCode_total;
  TH1D* mMatchedPdgCode_total_pyth;
  
  
  //QA Histograms//
  TH1D* mH_PairPt;
  TH1D* mH_PairPhi;
  TH1D* mH_PairEta;
  TH1D* mH_PairMass;
  TH1D* mH_PairPt_pyth;
  TH1D* mH_PairPhi_pyth;
  TH1D* mH_PairEta_pyth;
  TH1D* mH_PairMass_pyth;
  
  const TParticle* mParton1;
  const TParticle* mParton2;
  const TParticle* mParton3;
  const TParticle* mParton4;
  
  

  int iEntry; // number of event for pythia tree
  
  
public:
  
  embedAnalMaker(StMuDstMaker* muDstMaker);
  //~embedAnalMaker();
  
  Int_t          Init();
  Int_t          Make();
  Int_t          Finish();
  
  void           initializeHistograms();
  void           findWeight();
  void           setBinInfo(int axis, int nBins, double* limits);
  void           setPythiaTree(TTree* tree);
  void           setOutputFileName(string name) {mOutFileName = name;}
  
  void           createNMbins(binInfo bins, vector<vector<TH1D*> >& histMatrix, int histbins, double histmin, double histmax, string name);
  void           writeNMbins(binInfo bins, vector<vector<TH1D*> >& histMatrix);
  void           fillNMbin(binInfo bins, vector<vector<TH1D*> >& histMatrix, double yval, double zval, double fillvall);  

  int            binarySearch   (double val, int nBins, double* binLimits               );
  int            findValueInHalf(double val, double* binLimits, int lowLimit, int highLimit);
  void           moveSearchRange(int half, int& lowLimit, int& highLimit                );
  bool           lowestPossibleDivision(int lowLimit, int highLimit);
  
  bool             checkForTriggers(StMuDst* muDst);
  bool             checkTracksForPair(StMuTrack* track1, StMuTrack* track2);
  bool             checkParticlesForPair(const TParticle* particle1, const TParticle* particle2);
  int              findHighestRankingVertex(StMuDst* muDst);
  StMuTrack*       findTrack(StMuDst* muDst, int trackId, int vertexId);
  const TParticle* findPythParticle(StPythiaEvent* pythEvent,int particleId);
  bool             usePythParticle(const TParticle* particle, int pdgCode);

  double         getRadius(StMuTrack* track1, StMuTrack* track2);
  double         getRadius(const TParticle* particle1, const TParticle* particle2);
  double         getXvalue(int partonNumber, int eventNumber, StPythiaEvent* pythEvent);
  int            getProccessId(int eventNumber, StPythiaEvent* pythEvent);
  
  TLorentzVector getSumLV(StMuTrack* track1, StMuTrack* track2);
  TLorentzVector getSumLV(const TParticle* particle1,const TParticle* particle2);
  double         convertPhiScheme(double phi);
  
  
  
  
  
  
  
  //newer ones
  double getZValue(TLorentzVector pairLV, const TParticle* matchedParton);
  TVector3 TVector3FromTParticle(const TParticle* particle);
  const TParticle* matchPairToParton(TLorentzVector pairLV);
  void getPartonMomenta();
  void fillQAHistograms(TLorentzVector pairLV, string level);
  void fillAllHistograms(TLorentzVector pairLV, double fillVal, string fillType, string level);
  void writeAllHistograms();
  
  double getSpinTransFactor(const TParticle* matchedParton, TLorentzVector pairLV);
  bool fromPolarizedBeam(const TParticle* matchedParton);
  bool samePdgCode(const TParticle* part1, const TParticle* part2);
  bool particleAntiPariclePair(const TParticle* part1, const TParticle* part2);
  bool allSameFlavors();
  bool quarkGluonScat();
  bool quarkAntiquarkScat();
  bool quarkQuarkScat();
  bool sameParticle(const TParticle* part1, const TParticle* part2);
  double firstSpinTransFormula(double s, double t, double u);
  double secondSpinTransFormula(double s, double t, double u);
  double thirdSpinTransFormula(double s, double t, double u);
  void printHardScatScenario();
  
  
  
  
  
  void setTriggerSimuMaker(StTriggerSimuMaker* sim) {mTrgSimMkr = sim;}

  
  ClassDef(embedAnalMaker,1); //TCollectionProxyInfo errors when uncommented - will compile with no errors in starver new (root 5.34) but not in SL12d (root 5.22)
};

#endif
