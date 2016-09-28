/*
 *  makeEmbedPlots.C
 *  2012Embed
 *
 *  Created by Keith Landry on 3/31/16.
 *  UCLA
 *
 */

#include <vector>
#include "TH1D.h"
//#include "TGraphErrors.h"
//#include "TMultiGraph.h"

void getHistogramsFromFile(int nBins, int mBins, vector<vector<TH1D*> >& histMatrix, string name, TFile* file){

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
			
			TH1D* h = file->Get(fullname.c_str());

      tempVec.push_back(h);
		}
    histMatrix.push_back(tempVec);
  }
}

double** numberMatchedWithParton(string type, vector<vector<TH1D*> >& histMatrix){
 
  int binPdgCodeOffset = 1000;
  int lastQuarkCode    = 8;
  int nPairs = 0;
  int gluonPdgCode = 21;
  
  int a = histMatrix.size();
  int b = histMatrix[0].size();
  
  const int n = a;
  const int m = b;
  
  double** x = 0;
  x = new double*[n];
  
  if (type == "quark"){
    for (int i = 0; i < n; i++){
      x[i] = new double[m];
      for (int j = 0; j < m; j++){
        for (int bin = binPdgCodeOffset - lastQuarkCode; bin < binPdgCodeOffset + lastQuarkCode; bin++){
          x[i][j] += histMatrix[i][j]->GetBinContent(bin);
        }
      }
    }
  }
  else if (type == "gluon"){
    for (int i = 0; i < n; i++){
      x[i] = new double[m];
      for (int j = 0; j < m; j++){
        x[i][j] = histMatrix[i][j]->GetBinContent(binPdgCodeOffset + gluonPdgCode);
      }
    }
  }
  else 
    cout << "wrong parton type" << endl;
  
  return x;
}

double** errorMatchedWithParton(string type, vector<vector<TH1D*> >& histMatrix){
  
  int binPdgCodeOffset = 1000;
  int lastQuarkCode    = 8;
  int nPairs = 0;
  int gluonPdgCode = 21;
  
  int a = histMatrix.size();
  int b = histMatrix[0].size();
  
  const int n = a;
  const int m = b;
  
  double** x = 0;
  x = new double*[n];
  
  if (type == "quark"){
    for (int i = 0; i < n; i++){
      int err = 0;
      x[i] = new double[m];
      for (int j = 0; j < m; j++){
        for (int bin = binPdgCodeOffset - lastQuarkCode; bin < binPdgCodeOffset + lastQuarkCode; bin++){
          err += pow(histMatrix[i][j]->GetBinError(bin),2);
        }
        x[i][j] = sqrt(err);
      }
    }
  }
  else if (type == "gluon"){
    for (int i = 0; i < n; i++){
      x[i] = new double[m];
      for (int j = 0; j < m; j++){
        x[i][j] = histMatrix[i][j]->GetBinContent(binPdgCodeOffset + gluonPdgCode);
      }
    }
  }
  else 
    cout << "wrong parton type" << endl;
  
  return x;
}

double getPercQuark(double q, double g){
  double percQuark = q/(q+g); 
  return percQuark;
}
 
double getErrPercQuark(double q, double g, double Eq, double Eg){
  double EpercQuark  = sqrt(g**2 * Eq**2/pow(q+g, 4) + q**2 * Eg**2/pow(q+g, 4));
  return EpercQuark;
}

TGraphErrors** createGraph(double** ratio, double** Eratio, const int m, int markerStyle){
  
  TGraphErrors* graph[m];
  for (int i = 0; i < m; i++){
    graph[i] = new TGraphErrors(m,ratio[i],Eratio[i]);
    graph[i]->SetMarkerStyle(markerStyle);
  }
  return graph;
}

TMultiGraph* createMultiGraph(double* plottingPoints, const int n, const int m,
                              double** quark, double** gluon, double** Equark, double** Egluon, 
                              double** quarkPyth, double** gluonPyth, double** EquarkPyth, double** EgluonPyth){
  
  double percQuarkDet;
  double percQuarkPyth;
  double EpercQuarkDet;
  double EpercQuarkPyth;
  double ratio[n][m];
  double Eratio[n][m];

  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++){
            
      percQuarkDet   = getPercQuark(quark[i][j], gluon[i][j]);
      percQuarkPyth  = getPercQuark(quarkPyth[i][j], gluonPyth[i][j]);
      EpercQuarkDet  = getErrPercQuark(quark[i][j], gluon[i][j], Equark[i][j], Egluon[i][j]);
      EpercQuarkPyth = getErrPercQuark(quarkPyth[i][j], gluonPyth[i][j], EquarkPyth[i][j], EgluonPyth[i][j]);
      
      ratio[i][j]  = percQuarkDet/percQuarkPyth;
      Eratio[i][j] = sqrt(EpercQuarkDet**2/percQuarkPyth**2 + percQuarkDet**2*EpercQuarkPyth**2/pow(percQuarkPyth, 4));
      
    }
  }
  TMultiGraph* multiG = new TMultiGraph();
  addToMulti(n, m, plottingPoints, ratio, Eratio, multiG, 22);
  return multiG;  
}



void addToMulti(const int n, const int m, double* x, double** y, double** Ey, TMultiGraph* multiG, int markerStyle){
  
  //must reverse the order to plot in the wanted variable (first index)
  double yPlot[n][m];
  double EyPlot[n][m];
  
  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++){
      yPlot[j][i] = y[i][j];
      EyPlot[j][i] = Ey[i][j];
    }
  }
  
  for (int i = 0; i < m; i++)
  {
    TGraphErrors* g = new TGraphErrors(m, x, yPlot[i], 0, EyPlot[i]); 
    g->SetMarkerColor(i+1);
    g->SetMarkerStyle(markerStyle);
    g->SetLineColor(i+1);
    multiG->Add(g);  
  }
}



void drawNMbins_samePlot(vector<vector<TH1D*> >& histMatrix, vector<vector<TH1D*> >& histMatrix2,  TCanvas* can){
	
  int a = histMatrix.size();
  int b = histMatrix[0].size();
  const int nBins = a;
  const int mBins = b;
  
	int xSide;
	int ySide;
	int maxSide;
  int numHists = nBins*mBins;
	
	
	if (sqrt(numHists) != floor(sqrt(numHists)))
	{
		xSide = floor(sqrt(numHists)) + 1;
		ySide = numHists/xSide + 1;	
	}
	else
	{
		xSide = sqrt(numHists);
		ySide = xSide;
	}
	
	can->Divide(xSide, ySide);
	
  int k = 0;
	for (int i = 0; i < nBins; i++)
	{
    for (int j = 0; j < mBins; j++)
    {
      
      can->cd(k+1);
      histMatrix[i][j]->Draw();
      histMatrix2[i][j]->SetLineColor(2);
      histMatrix2[i][j]->Draw("SAME");
      k++;
    }
	}
	
}

void drawNMbins(vector<vector<TH1D*> >& histMatrix, TCanvas* can){
	
  int a = histMatrix.size();
  int b = histMatrix[0].size();
  const int nBins = a;
  const int mBins = b;
  
	int xSide;
	int ySide;
	int maxSide;
  int numHists = nBins*mBins;
	
	
	if (sqrt(numHists) != floor(sqrt(numHists)))
	{
		xSide = floor(sqrt(numHists)) + 1;
		ySide = numHists/xSide + 1;	
	}
	else
	{
		xSide = sqrt(numHists);
		ySide = xSide;
	}
	
	can->Divide(xSide, ySide);
	
  int k = 0;
	for (int i = 0; i < nBins; i++)
	{
    for (int j = 0; j < mBins; j++)
    {
      
      can->cd(k+1);
      histMatrix[i][j]->Draw();
      k++;
    }
	}
	
}


double ** getAvgKinematicValues(vector<vector<TH1D*> >& histMatrix){
  
  int a = histMatrix.size();
  int b = histMatrix[0].size();
  const int n = a;
  const int m = b;
  
  double ** points = 0;
  points = new double*[n];
  
  for (int i = 0; i < n; i++){
    points[i] = new double[m];
    for (int j = 0; j < m; j++){
      points[i][j] = histMatrix[i][j]->GetMean();
    }
  }
  return points;
}

double ** getErrKinematicValues(vector<vector<TH1D*> >& histMatrix){
  
  int a = histMatrix.size();
  int b = histMatrix[0].size();
  const int n = a;
  const int m = b;
  
  double ** points = 0;
  points = new double*[n];
  
  for (int i = 0; i < n; i++){
    points[i] = new double[m];
    for (int j = 0; j < m; j++){
      points[i][j] = histMatrix[i][j]->GetMeanError();
    }
  }
  return points;
}





void makeEmbedPlots(){
  
 // TFile* infile = new TFile("/Users/keithlandry/Desktop/Research/run12Embed/rootFiles/pt9_11_13044126_7.EmbedAnalysis.root");
  TFile* infile = new TFile("/Users/keithlandry/Desktop/Research/run12Embed/rootFiles/all5.root");
  
  //Histograms -- copied from embedAnalMaker.h
  vector<vector<TH1D*> > mH_ProcIds_PtM; 
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
  
  std::vector<std::vector<TH1D*> > mH_Pt_PtM;
  std::vector<std::vector<TH1D*> > mH_Pt_PtEta;
  std::vector<std::vector<TH1D*> > mH_Pt_MEta;
  std::vector<std::vector<TH1D*> > mH_M_PtM;
  std::vector<std::vector<TH1D*> > mH_M_PtEta;
  std::vector<std::vector<TH1D*> > mH_M_MEta;
  std::vector<std::vector<TH1D*> > mH_Eta_PtM;
  std::vector<std::vector<TH1D*> > mH_Eta_PtEta;
  std::vector<std::vector<TH1D*> > mH_Eta_MEta;
  
  std::vector<std::vector<TH1D*> > mH_Pt_PtM_pyth;
  std::vector<std::vector<TH1D*> > mH_Pt_PtEta_pyth;
  std::vector<std::vector<TH1D*> > mH_Pt_MEta_pyth;
  std::vector<std::vector<TH1D*> > mH_M_PtM_pyth;
  std::vector<std::vector<TH1D*> > mH_M_PtEta_pyth;
  std::vector<std::vector<TH1D*> > mH_M_MEta_pyth;
  std::vector<std::vector<TH1D*> > mH_Eta_PtM_pyth;
  std::vector<std::vector<TH1D*> > mH_Eta_PtEta_pyth;
  std::vector<std::vector<TH1D*> > mH_Eta_MEta_pyth;
  
  
  
  //QA Histograms//
  TH1D* mH_PairPt;
  TH1D* mH_PairPhi;
  TH1D* mH_PairEta;
  TH1D* mH_PairMass;
  TH1D* mH_PairPt_pyth;
  TH1D* mH_PairPhi_pyth;
  TH1D* mH_PairEta_pyth;
  TH1D* mH_PairMass_pyth;
  
  
  const int nPtBins  = 5;
  const int nMBins   = 5;
  const int nEtaBins = 4;
    
  
  //for testing 
  double ptPoints[5] = {3.5,4.5,5.87.3,9};
  
  
  getHistogramsFromFile(nPtBins, nMBins,   mH_ProcIds_PtM,   "hProcIds_PtMbin_", infile);
  getHistogramsFromFile(nPtBins, nEtaBins, mH_ProcIds_PtEta, "hProcIds_PtEtabin_", infile);
  getHistogramsFromFile(nMBins,  nEtaBins, mH_ProcIds_MEta,  "hProcIds_MEtabin_", infile);
  getHistogramsFromFile(nPtBins, nMBins,   mH_X1_PtM, "hX1_PtMbin_", infile);
  getHistogramsFromFile(nPtBins, nEtaBins, mH_X1_PtEta, "hX1_PtEtabin_", infile);
  getHistogramsFromFile(nMBins,  nEtaBins, mH_X1_MEta, "hX1_MEtabin_", infile);
  getHistogramsFromFile(nPtBins, nMBins,   mH_X2_PtM, "hX2_PtMbin_", infile);
  getHistogramsFromFile(nPtBins, nEtaBins, mH_X2_PtEta,  "hX2_PtEtabin_", infile);
  getHistogramsFromFile(nMBins,  nEtaBins, mH_X2_MEta,  "hX2_MEtabin_", infile);  
  getHistogramsFromFile(nPtBins, nMBins,   mH_Z_PtM, "hZ_PtMbin_", infile);
  getHistogramsFromFile(nPtBins, nEtaBins, mH_Z_PtEta,  "hZ_PtEtabin_", infile);
  getHistogramsFromFile(nMBins,  nEtaBins, mH_Z_MEta,  "hZ_MEtabin_", infile);  
  getHistogramsFromFile(nPtBins, nMBins,   mH_matchedPdgCode_PtM, "hMatchedPdgCode_PtMbin_", infile);
  getHistogramsFromFile(nPtBins, nEtaBins, mH_matchedPdgCode_PtEta,  "hMatchedPdgCode_PtEtabin_", infile);
  getHistogramsFromFile(nMBins,  nEtaBins, mH_matchedPdgCode_MEta, "hMatchedPdgCode_MEtabin_", infile);  
  getHistogramsFromFile(nPtBins, nMBins, mH_Pt_PtM, "hPt_PtMbin_", infile);
  getHistogramsFromFile(nPtBins, nEtaBins, mH_Pt_PtEta, "hPt_PtEtabin_", infile);
  getHistogramsFromFile(nMBins, nEtaBins, mH_Pt_MEta, "hPt_MEtabin_", infile);
  getHistogramsFromFile(nPtBins, nMBins, mH_M_PtM, "hM_PtMbin_", infile);
  getHistogramsFromFile(nPtBins, nEtaBins, mH_M_PtEta, "hM_PtEtabin_", infile);
  getHistogramsFromFile(nMBins, nEtaBins, mH_M_MEta, "hM_MEtabin_", infile);
  getHistogramsFromFile(nPtBins, nMBins, mH_Eta_PtM, "hEta_PtMbin_", infile);
  getHistogramsFromFile(nPtBins, nEtaBins, mH_Eta_PtEta, "hEta_PtEtabin_", infile);
  getHistogramsFromFile(nMBins, nEtaBins, mH_Eta_MEta, "hEta_MEtabin_", infile);
  
  
  getHistogramsFromFile(nPtBins, nMBins,   mH_ProcIds_PtM_pyth,  "hProcIds_pyth_PtMbin_", infile);
  getHistogramsFromFile(nPtBins, nEtaBins, mH_ProcIds_PtEta_pyth, "hProcIds_pyth_PtEtabin_", infile);
  getHistogramsFromFile(nMBins,  nEtaBins, mH_ProcIds_MEta_pyth, "hProcIds_pyth_MEtabin_", infile);
  getHistogramsFromFile(nPtBins, nMBins,   mH_X1_PtM_pyth, "hX1_pyth_PtMbin_", infile);
  getHistogramsFromFile(nPtBins, nEtaBins, mH_X1_PtEta_pyth, "hX1_pyth_PtEtabin_", infile);
  getHistogramsFromFile(nMBins,  nEtaBins, mH_X1_MEta_pyth, "hX1_pyth_MEtabin_", infile);
  getHistogramsFromFile(nPtBins, nMBins,   mH_X2_PtM_pyth, "hX2_pyth_PtMbin_", infile);
  getHistogramsFromFile(nPtBins, nEtaBins, mH_X2_PtEta_pyth, "hX2_pyth_PtEtabin_", infile);
  getHistogramsFromFile(nMBins,  nEtaBins, mH_X2_MEta_pyth, "hX2_pyth_MEtabin_", infile);  
  getHistogramsFromFile(nPtBins, nMBins,   mH_Z_PtM_pyth, "hZ_pyth_PtMbin_", infile);
  getHistogramsFromFile(nPtBins, nEtaBins, mH_Z_PtEta_pyth, "hZ_pyth_PtEtabin_", infile);
  getHistogramsFromFile(nMBins,  nEtaBins, mH_Z_MEta_pyth, "hZ_pyth_MEtabin_", infile);  
  getHistogramsFromFile(nPtBins, nMBins,   mH_matchedPdgCode_PtM_pyth, "hMatchedPdgCode_pyth_PtMbin_", infile);
  getHistogramsFromFile(nPtBins, nEtaBins, mH_matchedPdgCode_PtEta_pyth, "hMatchedPdgCode_pyth_PtEtabin_", infile);
  getHistogramsFromFile(nMBins,  nEtaBins, mH_matchedPdgCode_MEta_pyth, "hMatchedPdgCode_pyth_MEtabin_", infile);  
  getHistogramsFromFile(nPtBins, nMBins, mH_Pt_PtM_pyth, "hPt_pyth_PtMbin_", infile);
  getHistogramsFromFile(nPtBins, nEtaBins, mH_Pt_PtEta_pyth, "hPt_pyth_PtEtabin_", infile);
  getHistogramsFromFile(nMBins, nEtaBins, mH_Pt_MEta_pyth, "hPt_pyth_MEtabin_", infile);
  getHistogramsFromFile(nPtBins, nMBins, mH_M_PtM_pyth, "hM_pyth_PtMbin_", infile);
  getHistogramsFromFile(nPtBins, nEtaBins, mH_M_PtEta_pyth, "hM_pyth_PtEtabin_", infile);
  getHistogramsFromFile(nMBins, nEtaBins, mH_M_MEta_pyth, "hM_pyth_MEtabin_", infile);
  getHistogramsFromFile(nPtBins, nMBins, mH_Eta_PtM_pyth, "hEta_pyth_PtMbin_", infile);
  getHistogramsFromFile(nPtBins, nEtaBins, mH_Eta_PtEta_pyth, "hEta_pyth_PtEtabin_", infile);
  getHistogramsFromFile(nMBins, nEtaBins, mH_Eta_MEta_pyth, "hEta_pyth_MEtabin_", infile);
  
  
    
  double** matchedQuarkPtM = numberMatchedWithParton("quark", mH_matchedPdgCode_PtM);
  double** matchedGluonPtM = numberMatchedWithParton("gluon", mH_matchedPdgCode_PtM);
  double** errorQuarkPtM   = errorMatchedWithParton("quark", mH_matchedPdgCode_PtM);
  double** errorGluonPtM   = errorMatchedWithParton("gluon", mH_matchedPdgCode_PtM);

  
  double** matchedQuarkPythPtM = numberMatchedWithParton("quark", mH_matchedPdgCode_PtM_pyth);
  double** matchedGluonPythPtM = numberMatchedWithParton("gluon", mH_matchedPdgCode_PtM_pyth);
  double** errorQuarkPythPtM   = errorMatchedWithParton("quark", mH_matchedPdgCode_PtM_pyth);
  double** errorGluonPythPtM   = errorMatchedWithParton("gluon", mH_matchedPdgCode_PtM_pyth);
  
  //double** ptPoints_PtM = getAvgKinematicValues(mH_Pt_PtM); //these hists are empty still
  
  //TMultiGraph* multiPtM = createMultiGraph(ptPoints, nPtBins, nMBins, matchedQuarkPtM, matchedGluonPtM, errorQuarkPtM, errorGluonPtM, matchedQuarkPythPtM, matchedGluonPythPtM, errorQuarkPythPtM, errorGluonPythPtM);

  //multiPtM->Draw("AP");
  //multiPtM->GetXaxis()->SetTitle("Pt");
  //multiPtM->GetYaxis()->SetTitle("%pairs from q recon/%pairs from q pythia"); //problem getting axis before call to draw()
  //multiPtM->Draw("AP");

  
  
  TCanvas* cX12PtM = new TCanvas("cX12PtM");
  drawNMbins_samePlot(mH_X1_PtM,mH_X1_PtM_pyth,cX12PtM);
  
  TCanvas* CZPtM = new TCanvas("cZPtM");
  drawNMbins(mH_Z_PtM, CZPtM);
  
  
  double** x1PtM = getAvgKinematicValues(mH_X1_PtM);
  double** errX1PtM = getErrKinematicValues(mH_X1_PtM);
  double** x2PtM = getAvgKinematicValues(mH_X2_PtM_pyth);
  double** errX2PtM = getErrKinematicValues(mH_X2_PtM_pyth);
  
  
  TCanvas* cAvgX12PtM = new TCanvas("cAvgX1X2PtM");
  cAvgX12PtM->cd();
  TMultiGraph* multiX12PtM = new TMultiGraph();
  addToMulti(nPtBins, nMBins, ptPoints, x1PtM, errX1PtM, multiX12PtM, 21);
  addToMulti(nPtBins, nMBins, ptPoints, x2PtM, errX2PtM, multiX12PtM, 22);
  multiX12PtM->Draw("AP");


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  
  
  
  
}

