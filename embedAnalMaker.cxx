/*
 *  embedAnalMaker.cxx
 *  2012Embed
 *
 *  Created by Keith Landry on 2/3/16.
 *  UCLA
 *
 */

#include "embedAnalMaker.h"

ClassImp(embedAnalMaker) //TCollectionProxyInfo errors when uncommented

embedAnalMaker::embedAnalMaker(StMuDstMaker* muDstMaker) : StMaker() {
   
  mMuDstMaker = muDstMaker;

  string defaultName = "/star/u/klandry/gpfs01/2012Embed/output/defaultOutput.root";
  mOutFileName = defaultName;  
}

void embedAnalMaker::setBinInfo(int axis, int nBins, double* limits){
 
  //pt = 1, M = 2, eta = 3
  
    switch (axis)
    {
      case 1:
        mBinsPtM.ny      = nBins;
        mBinsPtM.ylims   = limits;
        mBinsPtEta.ny    = nBins;
        mBinsPtEta.ylims = limits;
        mBinLimits1Set   = true;
        break;
      case 2:
        mBinsPtM.nz      = nBins;
        mBinsPtM.zlims   = limits;
        mBinsMEta.ny     = nBins;
        mBinsMEta.ylims  = limits;
        mBinLimits2Set   = true;
        break;
      case 3:
        mBinsPtEta.nz    = nBins;
        mBinsPtEta.zlims = limits;
        mBinsMEta.nz     = nBins;
        mBinsMEta.zlims  = limits;
        mBinLimits3Set   = true;
        break;
      default:
        break;
    }
}

void embedAnalMaker::setPythiaTree(TTree* tree){
 
  mPythTree = (TTree*) tree;
  assert(mPythTree);
  mPythEvent = new StPythiaEvent;
  mPythTree->SetBranchAddress("PythiaBranch", &mPythEvent);
  mPythiaTreeSet = true;
}

void embedAnalMaker::fillAllHistograms(TLorentzVector sum, double fillVal, string fillType, string level){
  
  fillQAHistograms(sum, level);
  
  map <string, int> flags;
  flags["detector"]       = 10;
  flags["pythia"]         = 20;
  flags["procId"]         =  1;
  flags["X1"]             =  2;
  flags["X2"]             =  3;
  flags["Z"]              =  4;
  flags["matchedPdgCode"] =  5;
    
  switch (flags[fillType] + flags[level])
  {
    case 11:
      fillNMbin(mBinsPtM, mH_ProcIds_PtM, sum.Pt(), sum.M(), fillVal);
      fillNMbin(mBinsPtEta, mH_ProcIds_PtEta, sum.Pt(), sum.Eta(), fillVal);
      fillNMbin(mBinsMEta, mH_ProcIds_MEta, sum.M(), sum.Eta(), fillVal);
      break;
    case 12:
      fillNMbin(mBinsPtM, mH_X1_PtM, sum.Pt(), sum.M(), fillVal);
      fillNMbin(mBinsPtEta, mH_X1_PtEta, sum.Pt(), sum.Eta(), fillVal);
      fillNMbin(mBinsMEta, mH_X1_MEta, sum.M(), sum.Eta(), fillVal);      
      break;
    case 13:
      fillNMbin(mBinsPtM, mH_X2_PtM_pyth, sum.Pt(), sum.M(), fillVal);
      fillNMbin(mBinsPtEta, mH_X2_PtEta_pyth, sum.Pt(), sum.Eta(), fillVal);
      fillNMbin(mBinsMEta, mH_X2_MEta_pyth, sum.M(), sum.Eta(), fillVal);
      break;
    case 14:
      fillNMbin(mBinsPtM, mH_Z_PtM_pyth, sum.Pt(), sum.M(), fillVal);
      fillNMbin(mBinsPtEta, mH_Z_PtEta_pyth, sum.Pt(), sum.Eta(), fillVal);
      fillNMbin(mBinsMEta, mH_Z_MEta_pyth, sum.M(), sum.Eta(), fillVal);
      break;      
    case 15:
      fillNMbin(mBinsPtM, mH_matchedPdgCode_PtM, sum.Pt(), sum.M(), fillVal);
      fillNMbin(mBinsPtEta, mH_matchedPdgCode_PtEta, sum.Pt(), sum.Eta(), fillVal);
      fillNMbin(mBinsMEta, mH_matchedPdgCode_MEta, sum.M(), sum.Eta(), fillVal);
      break;
      
    case 21:
      fillNMbin(mBinsPtM, mH_ProcIds_PtM_pyth, sum.Pt(), sum.M(), fillVal);
      fillNMbin(mBinsPtEta, mH_ProcIds_PtEta_pyth, sum.Pt(), sum.Eta(), fillVal);
      fillNMbin(mBinsMEta, mH_ProcIds_MEta_pyth, sum.M(), sum.Eta(), fillVal);
      break;
    case 22:
      fillNMbin(mBinsPtM, mH_X1_PtM_pyth, sum.Pt(), sum.M(), fillVal);
      fillNMbin(mBinsPtEta, mH_X1_PtEta_pyth, sum.Pt(), sum.Eta(), fillVal);
      fillNMbin(mBinsMEta, mH_X1_MEta_pyth, sum.M(), sum.Eta(), fillVal);
      break;
    case 23:
      fillNMbin(mBinsPtM, mH_X2_PtM_pyth, sum.Pt(), sum.M(), fillVal);
      fillNMbin(mBinsPtEta, mH_X2_PtEta_pyth, sum.Pt(), sum.Eta(), fillVal);
      fillNMbin(mBinsMEta, mH_X2_MEta_pyth, sum.M(), sum.Eta(), fillVal);
      break;
    case 24:
      fillNMbin(mBinsPtM, mH_Z_PtM_pyth, sum.Pt(), sum.M(), fillVal);
      fillNMbin(mBinsPtEta, mH_Z_PtEta_pyth, sum.Pt(), sum.Eta(), fillVal);
      fillNMbin(mBinsMEta, mH_Z_MEta_pyth, sum.M(), sum.Eta(), fillVal);
      break;
    case 25:
      fillNMbin(mBinsPtM, mH_matchedPdgCode_PtM_pyth, sum.Pt(), sum.M(), fillVal);
      fillNMbin(mBinsPtEta, mH_matchedPdgCode_PtEta_pyth, sum.Pt(), sum.Eta(), fillVal);
      fillNMbin(mBinsMEta, mH_matchedPdgCode_MEta_pyth, sum.M(), sum.Eta(), fillVal);
      break;      
    default:
      cout << "filling histograms: flags not correct" << endl;
      break;
  }
}

void embedAnalMaker::fillQAHistograms(TLorentzVector pairLV, string level){
  
  map <string, int> flags;
  flags["detector"] = 1;
  flags["pythia"]   = 2;
  
  switch (flags[level])
  {
    case 1:
      mH_PairPt->Fill(pairLV.Pt(),mBinWeight);
      mH_PairEta->Fill(pairLV.Eta(),mBinWeight);
      mH_PairPhi->Fill(pairLV.Phi(),mBinWeight);
      mH_PairMass->Fill(pairLV.M(),mBinWeight);
      break;
    case 2:
      mH_PairPt_pyth->Fill(pairLV.Pt(),mBinWeight);
      mH_PairEta_pyth->Fill(pairLV.Eta(),mBinWeight);
      mH_PairPhi_pyth->Fill(pairLV.Phi(),mBinWeight);
      mH_PairMass_pyth->Fill(pairLV.M(),mBinWeight);   
      break;
    default:
      cout << "filling QA Histograms: flags not correct";
      break;
  }
}

void embedAnalMaker::writeAllHistograms(){
  
  writeNMbins(mBinsPtM, mH_ProcIds_PtM);
  writeNMbins(mBinsMEta, mH_ProcIds_MEta);
  writeNMbins(mBinsPtEta, mH_ProcIds_PtEta);
  writeNMbins(mBinsPtM, mH_X1_PtM);
  writeNMbins(mBinsMEta, mH_X1_MEta);
  writeNMbins(mBinsPtEta, mH_X1_PtEta);
  writeNMbins(mBinsPtM, mH_X2_PtM);
  writeNMbins(mBinsMEta, mH_X2_MEta);
  writeNMbins(mBinsPtEta, mH_X2_PtEta); 
  writeNMbins(mBinsPtM, mH_Z_PtM);
  writeNMbins(mBinsMEta, mH_Z_MEta);
  writeNMbins(mBinsPtEta, mH_Z_PtEta); 
  writeNMbins(mBinsPtM, mH_matchedPdgCode_PtM);
  writeNMbins(mBinsMEta, mH_matchedPdgCode_MEta);
  writeNMbins(mBinsPtEta, mH_matchedPdgCode_PtEta);  
  
  writeNMbins(mBinsPtM, mH_ProcIds_PtM_pyth);
  writeNMbins(mBinsMEta, mH_ProcIds_MEta_pyth);
  writeNMbins(mBinsPtEta, mH_ProcIds_PtEta_pyth);
  writeNMbins(mBinsPtM, mH_X1_PtM_pyth);
  writeNMbins(mBinsMEta, mH_X1_MEta_pyth);
  writeNMbins(mBinsPtEta, mH_X1_PtEta_pyth);
  writeNMbins(mBinsPtM, mH_X2_PtM_pyth);
  writeNMbins(mBinsMEta, mH_X2_MEta_pyth);
  writeNMbins(mBinsPtEta, mH_X2_PtEta_pyth);  
  writeNMbins(mBinsPtM, mH_Z_PtM_pyth);
  writeNMbins(mBinsMEta, mH_Z_MEta_pyth);
  writeNMbins(mBinsPtEta, mH_Z_PtEta_pyth);  
  writeNMbins(mBinsPtM, mH_matchedPdgCode_PtM_pyth);
  writeNMbins(mBinsMEta, mH_matchedPdgCode_MEta_pyth);
  writeNMbins(mBinsPtEta, mH_matchedPdgCode_PtEta_pyth);  
  
  mH_PairPt->Write();
  mH_PairPhi->Write();
  mH_PairEta->Write();
  mH_PairMass->Write();
  
  mH_PairPt_pyth->Write();
  mH_PairPhi_pyth->Write();
  mH_PairEta_pyth->Write();
  mH_PairMass_pyth->Write();
}

void embedAnalMaker::initializeHistograms(){
  
  cout << "initializing histograms" << endl;
  createNMbins(mBinsPtM, mH_ProcIds_PtM, 100, 0, 100, "hProcIds_PtMbin_");
  createNMbins(mBinsPtEta, mH_ProcIds_PtEta, 100, 0, 100, "hProcIds_PtEtabin_");
  createNMbins(mBinsMEta, mH_ProcIds_MEta, 100, 0, 100, "hProcIds_MEtabin_");
  createNMbins(mBinsPtM, mH_X1_PtM, 100, 0, 1, "hX1_PtMbin_");
  createNMbins(mBinsPtEta, mH_X1_PtEta, 100, 0, 1, "hX1_PtEtabin_");
  createNMbins(mBinsMEta, mH_X1_MEta, 100, 0, 1, "hX1_MEtabin_");
  createNMbins(mBinsPtM, mH_X2_PtM, 100, 0, 1, "hX2_PtMbin_");
  createNMbins(mBinsPtEta, mH_X2_PtEta, 100, 0, 1, "hX2_PtEtabin_");
  createNMbins(mBinsMEta, mH_X2_MEta, 100, 0, 1, "hX2_MEtabin_");  
  createNMbins(mBinsPtM, mH_Z_PtM, 100, 0, 1, "hZ_PtMbin_");
  createNMbins(mBinsPtEta, mH_Z_PtEta, 100, 0, 1, "hZ_PtEtabin_");
  createNMbins(mBinsMEta, mH_Z_MEta, 100, 0, 1, "hZ_MEtabin_");  
  createNMbins(mBinsPtM, mH_matchedPdgCode_PtM, 2000, -999, 1001, "hMatchedPdgCode_PtMbin_");
  createNMbins(mBinsPtEta, mH_matchedPdgCode_PtEta, 2000, -999, 1001, "hMatchedPdgCode_PtEtabin_");
  createNMbins(mBinsMEta, mH_matchedPdgCode_MEta, 2000, -999, 1001, "hMatchedPdgCode_MEtabin_");  
  
  createNMbins(mBinsPtM, mH_ProcIds_PtM_pyth, 100, 0, 100, "hProcIds_pyth_PtMbin_");
  createNMbins(mBinsPtEta, mH_ProcIds_PtEta_pyth, 100, 0, 100, "hProcIds_pyth_PtEtabin_");
  createNMbins(mBinsMEta, mH_ProcIds_MEta_pyth, 100, 0, 100, "hProcIds_pyth_MEtabin_");
  createNMbins(mBinsPtM, mH_X1_PtM_pyth, 100, 0, 1, "hX1_pyth_PtMbin_");
  createNMbins(mBinsPtEta, mH_X1_PtEta_pyth, 100, 0, 1, "hX1_pyth_PtEtabin_");
  createNMbins(mBinsMEta, mH_X1_MEta_pyth, 100, 0, 1, "hX1_pyth_MEtabin_");
  createNMbins(mBinsPtM, mH_X2_PtM_pyth, 100, 0, 1, "hX2_pyth_PtMbin_");
  createNMbins(mBinsPtEta, mH_X2_PtEta_pyth, 100, 0, 1, "hX2_pyth_PtEtabin_");
  createNMbins(mBinsMEta, mH_X2_MEta_pyth, 100, 0, 1, "hX2_pyth_MEtabin_");  
  createNMbins(mBinsPtM, mH_Z_PtM_pyth, 100, 0, 1, "hZ_pyth_PtMbin_");
  createNMbins(mBinsPtEta, mH_Z_PtEta_pyth, 100, 0, 1, "hZ_pyth_PtEtabin_");
  createNMbins(mBinsMEta, mH_Z_MEta_pyth, 100, 0, 1, "hZ_pyth_MEtabin_");  
  createNMbins(mBinsPtM, mH_matchedPdgCode_PtM_pyth, 2000, -999, 1001, "hMatchedPdgCode_pyth_PtMbin_");
  createNMbins(mBinsPtEta, mH_matchedPdgCode_PtEta_pyth, 2000, -999, 1001, "hMatchedPdgCode_pyth_PtEtabin_");
  createNMbins(mBinsMEta, mH_matchedPdgCode_MEta_pyth, 2000, -999, 1001, "hMatchedPdgCode_pyth_MEtabin_");  
  
    
  mH_PairPt   = new TH1D("hPairPt","hPairPt",300,0,30);
  mH_PairEta  = new TH1D("hPairEta","hPairEta",200,-2,2);
  mH_PairPhi  = new TH1D("hPairPhi","hPairPhi",200,-3.14159,3.14159);
  mH_PairMass = new TH1D("hPairMass","hPairMass",250,0,5); 
  
  mH_PairPt_pyth   = new TH1D("hPairPt_pyth","hPairPt_pyth",300,0,30);
  mH_PairEta_pyth  = new TH1D("hPairEta_pyth","hPairEta_pyth",200,-2,2);
  mH_PairPhi_pyth  = new TH1D("hPairPhi_pyth","hPairPhi_pyth",200,-3.14159,3.14159);
  mH_PairMass_pyth = new TH1D("hPairMass_pyth","hPairMass_pyth",250,0,5); 
}

void embedAnalMaker::createNMbins(binInfo bins, vector<vector<TH1D*> >& histMatrix, int histbins, double histmin, double histmax, string name){

  int nBins = bins.ny;
  int mBins = bins.nz;
  
	for (int iBin = 0; iBin < nBins; iBin++)
	{
    vector<TH1D*> tempVec;
		for (int jBin = 0; jBin < mBins; jBin++)
		{
			stringstream ss;
			ss << iBin;
			stringstream ss2;
			ss2 << jBin;
			
			string fullname = name + ss.str() + "_" + ss2.str();
			//cout << fullname << endl;
			
			TH1D* h = new TH1D(fullname.c_str(),fullname.c_str(),histbins,histmin,histmax);
      h->Sumw2();
      tempVec.push_back(h);
		}
    histMatrix.push_back(tempVec);
  }
}

void embedAnalMaker::writeNMbins(binInfo bins, vector<vector<TH1D*> >& histMatrix){
  
  int nBins = bins.ny;
  int mBins = bins.nz;

	for (int i=0; i<nBins; i++)
	{
		for (int j=0; j<mBins; j++)
		{
			histMatrix[i][j]->Write();
		}
	}	
}

void embedAnalMaker::fillNMbin(binInfo bins, std::vector<std::vector<TH1D*> >& histMatrix, double yval, double zval, double fillval){
  
  int yBin = binarySearch(yval, bins.ny, bins.ylims); 
  int zBin = binarySearch(zval, bins.nz, bins.zlims);
  
  if (yBin != -1 && zBin != -1){                             // -1 is return value of binarySearch when yval or zval is not in range
    histMatrix[yBin][zBin]->Fill(fillval, mBinWeight);
    //cout << "filled hist " << histMatrix[yBin][zBin]->GetName() << endl;
  }
}

int embedAnalMaker::binarySearch(double val, int nBins, double* binLimits){
  
	int currentLowLimit  = 0;
	int currentHighLimit = nBins;
	bool found = false;
  int binFound;
  
	while (!found){
		int halfFoundIn = findValueInHalf(val, binLimits, currentLowLimit, currentHighLimit);
    
		if (halfFoundIn == 1 && lowestPossibleDivision(currentLowLimit, currentHighLimit)){
			binFound = currentLowLimit;
			found = true;
		}
		else if (halfFoundIn == 2 && lowestPossibleDivision(currentLowLimit, currentHighLimit)){
			binFound = currentHighLimit - 1;
			found = true;
		}
		else if (halfFoundIn == 4){
			binFound = currentLowLimit;
			found = true;
		}
		else if (halfFoundIn == 5){
			binFound = currentHighLimit - 1;
			found = true;
		}
		//error case
		else if (halfFoundIn == 3){
			binFound = - 1;
			cout << "error: value not in range" << endl;
			found = true;
		}
		moveSearchRange(halfFoundIn, currentLowLimit, currentHighLimit);
	}
	return binFound;
}

int embedAnalMaker::findValueInHalf(double val, double* binLimits, int lowLimit, int highLimit){
  
	int halfwayBin = ceil((lowLimit + highLimit)/ 2);
  
	if      (val > binLimits[lowLimit] && val < binLimits[halfwayBin])   return 1;
	else if (val >= binLimits[halfwayBin] && val < binLimits[highLimit]) return 2;
	else if (val == binLimits[lowLimit])                                 return 4;
	else if (val == binLimits[highLimit])                                return 5;
  else return 3; //this handles errors resulting from value not being in search range at all
}

void embedAnalMaker::moveSearchRange(int half, int& lowLimit, int& highLimit){
  
	if (half == 1) highLimit = ceil ((lowLimit + highLimit)/2.); 
	if (half == 2) lowLimit  = floor((lowLimit + highLimit)/2.);
}

bool embedAnalMaker::lowestPossibleDivision(int lowLimit, int highLimit){
	if ( abs(highLimit - lowLimit) <=2 ) return true;
	else return false;
}

bool embedAnalMaker::checkForTriggers(StMuDst* muDst){
  
  bool triggerFired = false;
  bool printTrigs   = false;
  
  assert(mTrgSimMkr);
    
  vector<int> trigvector = mTrgSimMkr->triggerIds();
	
  if (printTrigs){
    for (int j = 0; j<trigvector.size(); j++)
      cout << "trigger keith" << trigvector[j] << endl;
  }
	
  //JP0,JP1,JP2,AJP
  if (mTrgSimMkr->isTrigger(370601) || mTrgSimMkr->isTrigger(370611) || mTrgSimMkr->isTrigger(370621) || mTrgSimMkr->isTrigger(370641))
    triggerFired = true;

  //BHT0VPD,BHT1VPD,BHT2BBC,BHT2
  if (mTrgSimMkr->isTrigger(370501) || mTrgSimMkr->isTrigger(370511) || mTrgSimMkr->isTrigger(370522) || mTrgSimMkr->isTrigger(370531))
    triggerFired = true;

  //vpd min bias?
  /*
  if (mTrgSimMkr.isTrigger(370001) || mTrgSimMkr.isTrigger(370011) || mTrgSimMkr.isTrigger(370983))
    triggerFired = true;
  //*/
  
  return triggerFired;
}

int embedAnalMaker::findHighestRankingVertex(StMuDst* muDst){
  
  int highRankId = 0;
  double highRank = 0;
  bool vertexFound = false;
  
  for (int iVert=0; iVert < muDst->numberOfPrimaryVertices(); iVert++)
  {
    StMuPrimaryVertex* vertex = muDst->primaryVertex(iVert);
    assert(vertex);
    
    double ranking = vertex->ranking();
    
    if (ranking > highRank && fabs(vertex->position().z()) < 60)
    {
      highRank   = ranking;
      highRankId = iVert;
      vertexFound = true;
    }
  }
  if (vertexFound) return highRankId;
}

StMuTrack* embedAnalMaker::findTrack(StMuDst* muDst, int trackId, int vertexId){
 
  StMuTrack* track = muDst->primaryTracks(trackId);
  if (track->vertexIndex() == vertexId                      && 
      track->nSigmaPion() > -1 && track->nSigmaPion() < 2.5 &&
      track->nHitsFit()/track->nHitsPoss() >- 0.5           &&
      track->pt() > 1.5 && fabs(track->eta()) < 2              )
  {
    return track;
  }
}

bool embedAnalMaker::checkTracksForPair(StMuTrack* track1, StMuTrack* track2){
  
  double minRad = 0.05;
  double maxRad = 0.7;
  
  if (track1->charge() != -1*track2->charge()) return false;
  else{
    double radius = getRadius(track1, track2);
    if (radius >= minRad && radius < maxRad) return true;
    else return false;
  }
}

double embedAnalMaker::getRadius(StMuTrack* track1, StMuTrack* track2){
 
  assert(track1->charge() == -1*track2->charge());
  
  StMuTrack* plus;
  StMuTrack* minus;
  
  if (track1->charge() == 1 && track2->charge() == -1){
    plus  = track1;
    minus = track2;
  }
  else{
    plus  = track2;
    minus = track1;
  }

  double pi = 3.14159265359;
  double radius;
  double radius_patch = 100; //must start with radius_patch higher than radius could be

  radius = sqrt( pow(plus->eta() - minus->eta(), 2) + pow(plus->phi() - minus->phi(), 2) );
  
  if (plus->phi() < 0 && minus->phi() > 0){
    radius_patch = sqrt( pow(plus->eta() - minus->eta(), 2) + pow(2.0*pi+plus->phi() - minus->phi(), 2) );
  }
  if (plus->phi() > 0 && minus->phi() < 0){
    radius_patch = sqrt( pow(plus->eta() - minus->eta(), 2) + pow(plus->phi() - (2.0*pi+minus->phi()) , 2) );
  }
  radius = min(radius, radius_patch);
  return radius;
}

double embedAnalMaker::getRadius(const TParticle* particle1, const TParticle* particle2){
  
  double pi = 3.14159265359;
  double radius;
  double radius_patch = 100; //must start with radius_patch higher than radius could be
  
  double phi1 = convertPhiScheme(particle1->Phi());
  double phi2 = convertPhiScheme(particle2->Phi());
  
  radius = sqrt( pow(particle1->Eta() - particle2->Eta(), 2) + pow(phi1 - phi2, 2) );

  if (phi1 < 0 && phi2 > 0){
    radius_patch = sqrt( pow(particle1->Eta() - particle2->Eta(), 2) + pow(2.0*pi+phi1 - phi2, 2) );
  }  
  if (phi1 > 0 && phi2 < 0){
    radius_patch = sqrt( pow(particle1->Eta() - particle2->Eta(), 2) + pow(phi1 - (2.0*pi+phi2) , 2) );
  }
  radius = min(radius, radius_patch);
  return radius;
}

double embedAnalMaker::convertPhiScheme(double phi){
  
  //convert from phi ranging from zero to 2*pi to phi ranging from -pi to pi
  double pi = 3.14159265359;
  double newPhi;

  if (phi > pi) newPhi = phi - 2*pi;
  else newPhi = phi;
  return newPhi;
}

double embedAnalMaker::getXvalue(int partonNumber, int eventNumber, StPythiaEvent* pythEvent){
  
  double X; 
  
  assert(pythEvent->eventId() == eventNumber);
  if (partonNumber == 1) X = pythEvent->x1();
  if (partonNumber == 2) X = pythEvent->x2();
  
  return X;
}

int embedAnalMaker::getProccessId(int eventNumber, StPythiaEvent* pythEvent){

  int proccessId;
  
  assert(pythEvent->eventId() == eventNumber);
  proccessId = pythEvent->processId();
  return proccessId;
}

TLorentzVector embedAnalMaker::getSumLV(StMuTrack* track1, StMuTrack* track2){
  
  TLorentzVector track1LV;
  TLorentzVector track2LV;
  TLorentzVector sum;
  
  double pionMass = 0.13957;
  
  track1LV.SetPtEtaPhiM(track1->pt(), track1->eta(), track1->phi(), pionMass);
  track2LV.SetPtEtaPhiM(track2->pt(), track2->eta(), track2->phi(), pionMass);
  
  sum = track1LV + track2LV;
  return sum;
}

TLorentzVector embedAnalMaker::getSumLV(const TParticle* particle1, const TParticle* particle2){
  
  TLorentzVector part1LV;
  TLorentzVector part2LV;
  TLorentzVector sum;
  
  double pionMass = 0.13957;
  
  part1LV.SetPtEtaPhiM(particle1->Pt(), particle1->Eta(), particle1->Phi(), pionMass);
  part2LV.SetPtEtaPhiM(particle2->Pt(), particle2->Eta(), particle2->Phi(), pionMass);
  
  sum = part1LV + part2LV;
  return sum;
}

void embedAnalMaker::findWeight(){
 
  Double_t crossSection[11] = {9.02,1.466,3.56E-1,1.52E-1,2.47E-2,5.8E-3,2.32E-3,3.41E-4,4.53E-5,1.029E-5,5.04E-7};
  Double_t numEvents[11]    = {54055,15453,15453,7733,7733,7733,4128,2582,2068,1965,601};
  Int_t    binEdges[12]     = {2,3,4,5,7,9,11,15,20,25,35,1000};
  
  int    bin;
  double binLumi; 
  mPythTree->GetEntry(0);
  double partonicPt = mPythEvent->pt();
  
  for (int i=0; i<11; i++){
    if (partonicPt >= binEdges[i] && partonicPt < binEdges[i+1]) bin = i;
  }
  
  binLumi    = (numEvents[bin]*crossSection[10])/(crossSection[bin]*numEvents[10]); // Weight normalized to Kevin's highest bin
  mBinWeight = 1.0/binLumi;

}

const TParticle* embedAnalMaker::findPythParticle(StPythiaEvent* pythEvent,int particleId){
 
  const TParticle* particle = pythEvent->particle(particleId);
  return particle;
}

bool embedAnalMaker::usePythParticle(const TParticle* particle, int pdgCode){
  
  if (particle->GetPdgCode() == pdgCode                 && 
      particle->Pt() > 1.5 && fabs(particle->Eta()) < 2 &&
      particle->GetNDaughters() == 0                      ) 
    return true;
  else return false;
}

bool embedAnalMaker::checkParticlesForPair(const TParticle* particle1, const TParticle* particle2){
  
  double minRad = 0.05;
  double maxRad = 0.7;
  
  double radius = getRadius(particle1, particle2);
  if (radius >= minRad && radius < maxRad) return true;
  else return false;
}

void embedAnalMaker::getPartonMomenta(){
  
  mParton1 = findPythParticle(mPythEvent, 4);
  mParton2 = findPythParticle(mPythEvent, 5);
  mParton3 = findPythParticle(mPythEvent, 6);
  mParton4 = findPythParticle(mPythEvent, 7);
}

const TParticle* embedAnalMaker::matchPairToParton(TLorentzVector pairLV){
  
  TVector3 parton3_3mom = TVector3FromTParticle(mParton3);
  TVector3 parton4_3mom = TVector3FromTParticle(mParton4);

  double dotProd3 = pairLV.Vect().Unit() * parton3_3mom.Unit();
  double dotProd4 = pairLV.Vect().Unit() * parton4_3mom.Unit();
  
  if (dotProd3 > dotProd4) return mParton3;
  else return mParton4;
}

TVector3 embedAnalMaker::TVector3FromTParticle(const TParticle* particle){
 
  TVector3 vec;
  vec.SetXYZ(particle->Px(), particle->Py(), particle->Pz());
  return vec;
}
  

double embedAnalMaker::getZValue(TLorentzVector pairLV, const TParticle* matchedParton){
  
  TVector3 partonMom = TVector3FromTParticle(matchedParton);
  double Z = pairLV.Vect().Mag()/partonMom.Mag();
  return Z;
}

double embedAnalMaker::getSpinTransFactor(const TParticle* matchedParton, TLorentzVector pairLV){
  
  double s = mPythEvent->s();
  double t = mPythEvent->t();
  double u = mPythEvent->u();
  double STF;

  if (mParton1->GetPdgCode() == 21) return 0; //gluons transfer no spin info from polarized beam
  else if (quarkGluonScat() && fromPolarizedBeam(matchedParton)){
    STF = thirdSpinTransFormula(s,t,u);
    //cout << "used 3rd form" << endl;
    return STF;
  }
  
  if (quarkQuarkScat()){
    if (allSameFlavors()){
      if ( (pairLV.Eta() > 0 && u < t) || (pairLV.Eta() < 0 && t < u) ){
        STF = firstSpinTransFormula(s,t,u);
        //cout << "used 1st form" << endl;
      }
    }
    else if (fromPolarizedBeam(matchedParton)){
      STF = thirdSpinTransFormula(s,t,u);
      //cout << "used 3rd form" << endl;
    }
    return STF;
  }
 
  if (quarkAntiquarkScat()){
    if (particleAntiPariclePair(mParton1, mParton2) && !samePdgCode(mParton1, mParton3) && !samePdgCode(mParton1, mParton4) ){
      if ( (mParton1->GetPdgCode()*mParton3->GetPdgCode() > 0 && sameParticle(mParton3, matchedParton)) ||
           (mParton1->GetPdgCode()*mParton4->GetPdgCode() > 0 && sameParticle(mParton4, matchedParton))  ){
        STF = thirdSpinTransFormula(s,t,u);
        //cout << "used 3rd form" << endl;
      }
    }
    else if (fromPolarizedBeam(matchedParton)){
      STF = secondSpinTransFormula(s,t,u);
      //cout << "used 2nd form" << endl;
    }
    return STF;
  }
  
  return 0; //if no condition is met spin transfer is zero
}
        
bool embedAnalMaker::fromPolarizedBeam(const TParticle* matchedParton){
 
  if ( (samePdgCode(mParton1, mParton3) && sameParticle(matchedParton, mParton3)) ||
       (samePdgCode(mParton1, mParton4) && sameParticle(matchedParton, mParton4))  )
    return true;
  else return false;
}
        

bool embedAnalMaker::samePdgCode(const TParticle* part1, const TParticle* part2){
 
  if (part1->GetPdgCode() == part2->GetPdgCode()) return true;  
  else return false;
}

bool embedAnalMaker::quarkQuarkScat(){

  if (!quarkGluonScat() && !quarkAntiquarkScat()) return true;
  else return false;
}

bool embedAnalMaker::particleAntiPariclePair(const TParticle* part1, const TParticle* part2){
  
  if (part1->GetPdgCode() == -1*part2->GetPdgCode()) return true;
  else return false;
}

bool embedAnalMaker::allSameFlavors(){

  if (samePdgCode(mParton1, mParton2) && 
      samePdgCode(mParton1, mParton3) && 
      samePdgCode(mParton1, mParton4))
    return true;
  else return false;
}

bool embedAnalMaker::quarkAntiquarkScat(){
 
  if (!quarkGluonScat() && mParton1->GetPdgCode()*mParton2->GetPdgCode() < 0) return true;
  else return false;
}

bool embedAnalMaker::quarkGluonScat(){

  if ( (mParton1->GetPdgCode() == 21 || mParton2->GetPdgCode() == 21) &&
       (mParton1->GetPdgCode() != 21 || mParton2->GetPdgCode() != 21)  )
    return true;
  else return false;
}

double embedAnalMaker::firstSpinTransFormula(double s, double t, double u){
  
  double STF = -2*s*u*(1-t/(3*u)) / (pow(s,2)+pow(u,2)+(pow(s,2)+pow(t,2))*pow(t,2)/pow(u,2)-2*pow(s,2)*t/(3*u));
  return STF;
}

double embedAnalMaker::secondSpinTransFormula(double s, double t, double u){
 
  double STF = -2*s*u*(1-t/(3*s)) / (pow(s,2)+pow(u,2)+(pow(u,2)+pow(t,2))*pow(t,2)/pow(s,2)-2*pow(u,2)*t/(3*s));
  return STF;
}

double embedAnalMaker::thirdSpinTransFormula(double s, double t, double u){
  double STF = -2*s*u/(pow(s,2)+pow(u,2));
  return STF;
}

bool embedAnalMaker::sameParticle(const TParticle* part1, const TParticle* part2){

  if (samePdgCode(part1, part2)){
    if (part1->Px() == part2->Px() &&
        part1->Py() == part2->Py() &&
        part1->Pz() == part2->Pz()   )
      return true;
  }
  else return false;
}

void embedAnalMaker::printHardScatScenario(){
  
  cout << "\nprocId: " << mPythEvent->processId();
  cout << "\t" << mParton1->GetPdgCode();
  cout << "\t" << mParton2->GetPdgCode();
  cout << "\t" << mParton3->GetPdgCode();
  cout << "\t" << mParton4->GetPdgCode() << endl;
  if (mParton1->GetPdgCode() == 21) cout << "parton1 is gluon disreguarding event" << endl;
  else{
    if (quarkGluonScat())     cout << "determined qg scattering" << endl;
    if (quarkQuarkScat())     cout << "determined qq scattering" << endl;
    if (quarkAntiquarkScat()) cout << "determined qqbar scattering" << endl;
    if (allSameFlavors())     cout << "determined all same flavor" << endl;
    if (particleAntiPariclePair(mParton1, mParton2)) cout << "determined part anti pair" << endl;
  }
}




//==============================================================================
Int_t embedAnalMaker::Init(){
  
  assert(mBinLimits1Set && mBinLimits2Set && mBinLimits3Set);
  assert(mPythiaTreeSet);
  
  mTrgSimMkr = dynamic_cast<StTriggerSimuMaker*>(GetMakerInheritsFrom("StTriggerSimuMaker"));

  mOutputFile = new TFile(mOutFileName.c_str(),"RECREATE");
  iEntry = 0;
  
  mMatchedPdgCode_total = new TH1D("hMatchedPdgCode_total", "hMatchedPdgCode_total", 2000, -999, 1001);
  mMatchedPdgCode_total_pyth = new TH1D("hMatchedPdgCode_total_pyth", "hMatchedPdgCode_total_pyth", 2000, -999, 1001);
  initializeHistograms();
  findWeight();
  
  return kStOK;
}

Int_t embedAnalMaker::Finish(){
  
  writeAllHistograms();
  mMatchedPdgCode_total->Write();
  mMatchedPdgCode_total_pyth->Write();
  mOutputFile->Write();
  cout << "\n";
  cout << "\n";
  cout << "finished" << endl;
  cout << "\n";
  cout << "\n";
   
	return kStOK;
}

Int_t embedAnalMaker::Make(){
  assert(mMuDstMaker);
	
	StMuDst* muDst = mMuDstMaker->muDst();	
	assert(muDst);
  
  mPythTree->GetEntry(iEntry);
  assert(mPythEvent);

  //find event wide values
  int eventId    = muDst->event()->eventId();
  double X1      = getXvalue(1, eventId, mPythEvent);
  double X2      = getXvalue(2, eventId, mPythEvent);
  int proccessId = getProccessId(eventId, mPythEvent);
  getPartonMomenta();
  
  //dectector level
  if (checkForTriggers(muDst))
  {
    int vertexId = -100;
    vertexId = findHighestRankingVertex(muDst);
    if (vertexId < 0) return kStOK; // no viable vertex in this event

    for (int i = 0; i<muDst->numberOfPrimaryTracks() - 1; i++)
    {
      StMuTrack* track1 = findTrack(muDst, i, vertexId);
      if (!track1) continue;
      
      for (int j = i+1; j<muDst->numberOfPrimaryTracks(); j++) //don't double count
      {
        StMuTrack* track2 = findTrack(muDst, j, vertexId); 
        if (!track2) continue;
        
        if (checkTracksForPair(track1, track2))
        {
          TLorentzVector sum = getSumLV(track1, track2);
          
          const TParticle* matchedParton = matchPairToParton(sum);
          double Z = getZValue(sum, matchedParton);
          //printHardScatScenario();
          double d_nn = getSpinTransFactor(matchedParton, sum);
          
          mMatchedPdgCode_total->Fill(matchedParton->GetPdgCode());
          fillAllHistograms(sum, proccessId, "procId", "detector");
          fillAllHistograms(sum, X1, "X1", "detector");
          fillAllHistograms(sum, X2, "X2", "detector");
          fillAllHistograms(sum, Z, "Z", "detector");
          fillAllHistograms(sum, matchedParton->GetPdgCode(), "matchedPdgCode", "detector");
        }
      }//end j
    }//end i
  }//end detector level

  //pythia level
  for (int i=0; i<mPythEvent->particles()->GetEntries(); i++)
  {
    const TParticle* piPlus = findPythParticle(mPythEvent,i);
    if (!usePythParticle(piPlus, 211)) continue;
    
    for (int j=0; j<mPythEvent->particles()->GetEntries(); j++)
    {
      const TParticle* piMinus = findPythParticle(mPythEvent,j);
      if (!usePythParticle(piMinus, -211)) continue;

      if (checkParticlesForPair(piPlus, piMinus))
      {
        TLorentzVector sum = getSumLV(piPlus, piMinus);
                
        const TParticle* matchedParton = matchPairToParton(sum);
        double Z = getZValue(sum, matchedParton);
        //printHardScatScenario();
        double d_nn = getSpinTransFactor(matchedParton, sum);

        mMatchedPdgCode_total_pyth->Fill(matchedParton->GetPdgCode());
        fillAllHistograms(sum, proccessId, "procId", "pythia");
        fillAllHistograms(sum, X1, "X1", "pythia");
        fillAllHistograms(sum, X2, "X2", "pythia");  
        fillAllHistograms(sum, Z, "Z", "pythia");
        fillAllHistograms(sum, matchedParton->GetPdgCode(), "matchedPdgCode", "pythia");
      }
    }//end j
  }//end i
  iEntry++;
  return kStOK;

}


