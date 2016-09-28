/*
 *  runEmbedAnal.C
 *  2012Embed
 *
 *  Created by Keith Landry on 2/5/16.
 *  UCLA
 *
 */

//#include "StTriggerSimuMaker.h"


void loadLibraries(){
  cout << "\n";
  //gROOT->Macro("/star/u/klandry/gpfs01/2012Embed/code/LoadLibs.C");
  gROOT->Macro("/global/homes/k/klandry/2012Embed/code/LoadLibs.C");
  gSystem->Load("StEEmcA2EMaker");
  gSystem->Load("embedAnalMaker");
  cout << " loading of embedAnalMaker library done" << endl;
}

string getPythiaFileNameFromMudstName(string muDstFileName){
  string pythFileName = muDstFileName;
  string searchStr  = "MuDst";
  string replaceStr = "pythia";
  
  pythFileName.replace(pythFileName.find(searchStr), searchStr.length(), replaceStr);
  return pythFileName;
}

string createOutputFileName(string muDstFileName){
  
  string searchStr  = "pt";  
  std::size_t found = muDstFileName.rfind(searchStr);
  string outFileName = muDstFileName.substr(found);

  string searchStr2 = "MuDst";
  string replaceStr = "EmbedAnalysis";
  outFileName.replace(outFileName.find(searchStr2), searchStr2.length(), replaceStr);

  string path = "/global/homes/k/klandry/2012Embed/outputfull/";
  //string path = "/star/u/klandry/gpfs01/2012Embed/output/";

  //outFileName = path + outFileName;
  return outFileName;
  
}

//==============================================================================
void runEmbedAnal(string muDstFileName){
  
  //string muDstFileName = "/star/u/klandry/gpfs01/2012Embed/rootFilesForTest/pt9_11_13044126_7.MuDst.root";
  //string muDstFileName = "/global/homes/k/klandry/2012Embed/rootFilesForTest/pt25_35_13049042_1.MuDst.root";
  loadLibraries();
  
  string pythFileName = getPythiaFileNameFromMudstName(muDstFileName);  
  cout << "\npythia file: " << pythFileName << endl;
  TFile* pythFile = new TFile(pythFileName.c_str());
  TTree* pythTree = (TTree*)pythFile->Get("PythiaTree");
  assert(pythTree);

  string outFileName = createOutputFileName(muDstFileName);
  //string outFileName = "test123.root";

  cout << "\noutput file: " << outFileName << endl;
  
  //SET UP CHAIN OF MAKERS
  cout << "\n";
	cout << "Creating chain ......." << endl;
  StChain* chain = new StChain; 
	
  //MUDST READER
  cout << "\n";
	cout << "adding StMuDstMaker ......." << endl;
  StMuDstMaker* muDstMaker = new StMuDstMaker(0,0,"",muDstFileName.c_str(),"",100000,"MuDst");

  //StarDbMaker
  cout << "\n";
	cout << "adding St_db_Maker ......." << endl;
  St_db_Maker* dbMaker = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb");
  
  //ACCESS ENDCAP DATABASE
  StEEmcDbMaker* eemcDb = new StEEmcDbMaker;

  bool realdata   = true;  //for data or embedding 
	bool simulation = false;  
  
  if (simulation)
	{
		//EMC SIMULATOR MAKER
		//for simulated data for year > 2009 must have:
		StEmcSimulatorMaker* ecmSimu = new StEmcSimulatorMaker();                   // needed by StTriggerSimuMaker for use with simulated data
		
	}	
  if (realdata)
  {
    //EmcADCtoEMaker (needed for trigger sim maker)
    StEmcADCtoEMaker* adc = new StEmcADCtoEMaker;
    adc->saveAllStEvent(true);
    
    //EEMC ADC TO ENERGY MAKER - don't actually need to include this. 
    //also for real data - used in conjuntion with StEmcADCtoEMaker instead of StEmcSimulatorMaker
    //StEEmcA2EMaker *a2EMakerPtr = new StEEmcA2EMaker("EEmcA2EMaker");
    //a2EMakerPtr->database("eemcDb");          // sets db connection
    //a2EMakerPtr->source("MuDst",1);           // sets mudst as input
    
  }
    
  //TriggerSimuMaker
  bool useOnline = false;
  StTriggerSimuMaker* simuTrig = new StTriggerSimuMaker;
  simuTrig->setMC(2); //0=data, 1=simulation, 2=embedding
  simuTrig->useBemc();
  simuTrig->useEemc();
  
  if (useOnline)
  {
    simuTrig->bemc->setConfig(StBemcTriggerSimu::kOnline);
    simuTrig->useOnlineDB();
  }
  else
  {
    simuTrig->bemc->setConfig(StBemcTriggerSimu::kOffline);
    simuTrig->useOfflineDB();
    
    //Taken from kevin's example
    // Set trigger thresholds
    simuTrig->setBarrelJetPatchTh(0,20);
    simuTrig->setBarrelJetPatchTh(1,28);
    simuTrig->setBarrelJetPatchTh(2,36);
    
    simuTrig->setOverlapJetPatchTh(0,20);
    simuTrig->setOverlapJetPatchTh(1,28);
    simuTrig->setOverlapJetPatchTh(2,36);
    
    simuTrig->setEndcapJetPatchTh(0,20);
    simuTrig->setEndcapJetPatchTh(1,28);
    simuTrig->setEndcapJetPatchTh(2,36);
    
    simuTrig->setBarrelHighTowerTh(0,11);
    simuTrig->setBarrelHighTowerTh(1,15);
    simuTrig->setBarrelHighTowerTh(2,18);
    
    simuTrig->setEndcapHighTowerTh(0,18);
    simuTrig->setEndcapHighTowerTh(1,25);
    
    // Define triggers using run 13072020 /star/u/jkadkins/KeithTriggerDef.txt
    // triggerIndex,name,triggerId,onbits,offbits,onbits1,onbits2,onbits3,offbits1,offbits2,offbits3
    
    simuTrig->emc->defineTrigger(8, "JP0"          ,370601,0x00008000,0x00000000,0x80000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000001);
    simuTrig->emc->defineTrigger(9, "JP1"          ,370611,0x00000040,0x00000000,0x00400000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000001);
    simuTrig->emc->defineTrigger(23,"JP2"          ,370621,0x00000080,0x00000000,0x00800000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000001);
    simuTrig->emc->defineTrigger(18,"AJP"          ,370641,0x00001000,0x00000000,0x10000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000001);
    
    simuTrig->emc->defineTrigger(15,"BHT0*VPDMB"   ,370501,0x00000001,0x00000000,0x00010000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000001);
    simuTrig->emc->defineTrigger(17,"BHT1*VPDMB"   ,370511,0x00000002,0x00000000,0x00020000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000001);
    simuTrig->emc->defineTrigger(5, "BHT2*BBCMB"   ,370522,0x00000004,0x00000000,0x00040000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000001);
    simuTrig->emc->defineTrigger(4, "BHT2"         ,370531,0x00000004,0x00000000,0x00040000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000001);

    /* triggers that kevin added that I don't use in my analysis - figured I'd keep them here just in case
    simuTrig->emc->defineTrigger(10,"BBCMB"        ,370022,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000001);
    simuTrig->emc->defineTrigger(11,"VPDMB"        ,370001,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000001);
    simuTrig->emc->defineTrigger(12,"VPDMB-nobsmd" ,370011,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000001);
    simuTrig->emc->defineTrigger(54,"JP2*L2JetHigh",370982,0x00000080,0x00000000,0x00800000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000001);
    simuTrig->emc->defineTrigger(55,"VPDMB-peds"   ,370983,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000001);
    */
    
  }

  //embedAnalMaker
  cout << "\n";
	cout << "adding embedAnalMaker ......." << endl;
  embedAnalMaker* eaMaker = new embedAnalMaker(muDstMaker);
 
  cout << "\n";
	cout << "setting pythia tree in embedAnalMaker ......." << endl;
  eaMaker->setPythiaTree(pythTree);
  
  const int nPtBins    = 5;
  const int nEtaBins   = 4;
  const int nMassBins  = 5;
  
  double    ptBinLimits[nPtBins+1]   = {3.0, 4.0, 5.0, 6.5, 8.0, 50};
  double    mBinLimits[nMassBins+1]  = {0.0, 0.4, 0.6, 0.8, 1.0, 100.5};
  double    etaBinLimits[nEtaBins+1] = {-2,-0.5,0.0,0.5, 2};
  
  cout << "\n";
	cout << "setting hist bin limits in embedAnalMaker ......." << endl;
  eaMaker->setBinInfo(1, nPtBins,   ptBinLimits);
  eaMaker->setBinInfo(2, nMassBins, mBinLimits);
  eaMaker->setBinInfo(3, nEtaBins,  etaBinLimits);
  
  eaMaker->setOutputFileName(outFileName.c_str());

  cout << "\n";
	cout << "\n";
	cout << "Init() ......." << endl;
	cout << "\n";
	cout << "\n";
  
	chain->Init();
  
	cout << "\n";
	cout << "\n";
	cout << "Make() ......." << endl;
	cout << "\n";
	cout << "\n";
  
  for(int iEvent = 0; iEvent < 10000; iEvent++)
  {
    chain->Clear();
    int status = chain->Make(iEvent);
    if(status == kStSkip) continue;
    if(status % 10 == kStEOF || status % 10 == kStFatal) break;
  }
  
  
  
	
  cout << "\n";
	cout << "\n";
	cout << "Finish() ......." << endl;
	cout << "\n";
	cout << "\n";
  
	chain->Finish();
	
	
	delete chain;
  //*/

}