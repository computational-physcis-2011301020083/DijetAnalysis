///////////////////////////////////////////////////////////////
//
// Fit MC with an analytic function 
// + draw pseudo-experiments from the fit 
//
///////////////////////////////////////////////////////////////

#include "../ATLASStyle/AtlasStyle.C"                                                                                                                                      
#include "../ATLASStyle/AtlasUtils.C"
#include "../code/helperFunctions/histPlotFormatFunctions.C"
#include "../code/helperFunctions/statFitFunctions.C"
#include "../code/helperFunctions/fitFunctions.C" 
#include "../code/helperFunctions/swiftBkgGenerator.C"

using namespace std; 
using namespace TMath;


//====================
//
//  Inputs and outputs  
//
//====================

TString dataInput   = "data/dijets/dataLikeHistograms.2016.37fb-Resonance.root"; 
TString dataDir     = "Nominal"; 
TString dataHist    = "mjj_Data_2016_37p4fb"; 

TString mcInput     = "MC/mc_moriond_NLO_EW.root"; 
TString mcHistName  = "mjj_mc15_13TeV_361023_Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3W_total_final"; 

TString pdf         = "smoothMC/mjj_pythia_NLO_EW_37fb_6paramFit.pdf"; 
TString outFile     = "smoothMC/mjj_pythia_NLO_EW_37fb_6paramFit.root"; 

TString outPEPath   = "smoothMC/PE/";

//============
//
// Set options
//
//============

double xminfit      = 1100  ; // GeV
double xmaxfit      = 9000 ; // GeV 
double lumi         = 37.0 ; // fb-1 
TString bkgFunc     = "dijet6Param" ; 
double par[6]       = {0.06, 6.2, -3.5, 2.6, 1.3, 0.2}; 

bool doAndSavePEs   = 1       ; 
int numPE           = 1000    ; 


//////////////////////////////////////////////////////////////////
//
// Main code starts here
//
//////////////////////////////////////////////////////////////////

void smoothMC()
{ 

  SetAtlasStyle (); 

  //==// Make TCanvas to store plots 
  TCanvas *canvas = new TCanvas();
  canvas -> Update();
  canvas -> Print(pdf+"[");
  gPad   -> SetLeftMargin(0.12);
  gPad   -> SetBottomMargin(0.12);
  canvas -> SetLogy();
  //canvas -> SetLogx();  

  //==============================
  //
  // Open & pull in data hist
  //
  //==============================  
  TFile *dFile    = TFile::Open( dataInput, "r" );
  dFile -> cd();
  TDirectoryFile *d = (TDirectoryFile*)dFile->Get(dataDir);
  TH1D *dHist = (TH1D*)d->Get(dataHist);
  formatHist(dHist);
  dHist -> SetMarkerSize(0.5); 
  int binMin = dHist -> GetXaxis()->FindBin(xminfit); 
  int binMax = dHist -> GetXaxis()->FindBin(xmaxfit); 
  dHist -> GetXaxis() -> SetRangeUser(xminfit, xmaxfit);

  //===================================
  //
  // Open MC file and pull in histogram  
  //
  //===================================
  TFile *mFile    = TFile::Open( mcInput, "r" );
  mFile -> cd(); 
  /*
  //==// fill vector with histograms that contain "mjj" in their names in the TFile 
  vector<TH1*> mVec;
  mVec = getHist(mFile, mVec, "mjj");

  //==// add all the mass histograms into one
  TH1D *mHist=(TH1D*)mVec[1]->Clone();
  for(int i= 2; i < mVec.size(); i++)
    mHist->Add(mVec[i]);         
  */  

  TH1D *mHist = (TH1D*)mFile->Get(mcHistName);
  formatHist(mHist, "Events", "m_{jj} [GeV]");
  mHist -> SetMarkerSize(0.1); 
  mHist -> SetLineColor(kBlue);  
  mHist -> SetMarkerColor(kBlue); 
  mHist -> Scale(dHist->Integral(binMin, binMax)/mHist->Integral(binMin, binMax) ); 
  mHist -> SetLineWidth(2); 

  //==// set sqrtN error 
  for (int i = 1; i <= mHist->GetNbinsX(); i++) 
    mHist -> SetBinError( i, sqrt(mHist->GetBinContent(i)) );

  //==// draw to canvas 
  mHist -> Draw("hist ][");
  dHist -> Draw("same ]["); 
  mHist -> SetMinimum(3e-3);
  mHist -> GetXaxis() -> SetRangeUser(xminfit, xmaxfit);
  ATLAS_LABEL(0.70, 0.85 ,kBlack, 0.045); 
  myText(0.70, 0.80, kBlack, (char*)"#sqrt{s} = 13TeV, 37fb^{-1}", 0.035);
  myText(0.70, 0.76, kBlack, (char*)"Data", 0.035);
  myText(0.70, 0.72, kBlue, (char*)"NLO+LO MC", 0.035);
  canvas -> Print(pdf); canvas -> Clear(); 

  //==// create bin-width divided hist 
  TH1D *mHistBinWidDiv=(TH1D*)mHist->Clone();
  for (int i = 1; i <= mHist->GetNbinsX(); i++) 
  {    
    mHistBinWidDiv -> SetBinContent(i, mHist->GetBinContent(i)/mHist->GetBinWidth(i) );                                                                     
    mHistBinWidDiv -> SetBinError(i, mHist->GetBinError(i)/mHist->GetBinWidth(i) );
  }  
  
  //==// Get difference in data and MC  
  canvas -> SetLogy(0);
  TH1D *diff = (TH1D*)dHist->Clone("Diff");
  diff -> Add(mHist, -1);
  for (int i = 1; i <= diff->GetNbinsX(); i++) 
  { 
    diff -> SetBinContent(i, diff->GetBinContent(i)/mHist->GetBinContent(i)) ;    
    diff -> SetBinError(i, diff->GetBinError(i)/mHist->GetBinContent(i)) ;    
  }
  diff -> SetMaximum(0.5); 
  diff -> SetMinimum(-0.5); 
  //divideUncert(diff, dHist);
  diff -> Draw("e1 ][");  
  diff -> GetYaxis() -> SetTitle("#frac{Data-MC}{MC}");
  diff -> GetXaxis() -> SetTitle("m_{jj} [GeV]");
  makeTLine(xminfit, 0, xmaxfit,0, kBlack, "hist same", 1, 2);  
  ATLAS_LABEL(0.20, 0.85 ,kBlack, 0.045);
  canvas -> Print(pdf); canvas -> Clear();   
  
  //===================================
  //
  // Fit MC with function and plot  
  //
  //===================================
  gErrorIgnoreLevel = kFatal;

  //==// set up global function 
  int numParam = sizeof(par)/sizeof(par[0]); 
  vector<double> nomGlob_fitQuality; 
  TF1 *BFuncGlobal = new TF1("BFuncGlobal", getFunction(bkgFunc) , xminfit, xmaxfit);  
  formatFunc(BFuncGlobal, kRed, 2);
  for (int i = 0; i < numParam; i++) 
    BFuncGlobal -> SetParameter(i, par[i]);

  //==// fit 
  double BFuncGlobalLLH = minuitFit ( mHist, BFuncGlobal, 0, xminfit, xmaxfit, 0, &nomGlob_fitQuality);                          
  if ( nomGlob_fitQuality[0] || nomGlob_fitQuality[1] != 0 )  
  { 
    cout << endl; 
    printStuff("../ascii/crab.txt");
    cout << "NOMINAL global fit is failing ... retry with different initial conditions in the config. " << endl; 
    cout << endl; 
    exit(0); 
  }

  //==// convert fit to histogram 
  TH1D *BHistGlobal   = (TH1D*) mHist -> Clone("BHistGlobal");
  fit_to_hist(BHistGlobal, BFuncGlobal, xminfit, xmaxfit); 

  //==// plot MC and fit
  plotBkg (canvas, pdf, mHist, BHistGlobal, "NLO+EW MC", to_string(numParam)+" Paramter Fit" );
 
  //==// plot MC and data 
  plotBkg (canvas, pdf, dHist, mHist, "Data", "NLO+EW MC");
  
  //===================================
  //
  // Generate PE from fit to MC
  //
  //===================================
  canvas -> SetLogy();
  if (doAndSavePEs)
  { 
    for (int i = 1; i<=numPE; i++ )
    {
      TFile *outPEFile; 
      outPEFile = new TFile( outPEPath+"PE"+to_string(i)+".root"  , "RECREATE");     
      TH1D *PEHist = doPseudoExp(BHistGlobal, i);
      PEHist -> Draw("E1"); 
      PEHist -> SetLineColor(kBlack);
      PEHist -> Write("myy");
      canvas -> Print(pdf); canvas -> Clear();     
      outPEFile  -> Close();
    }
  }  
  canvas->Print(pdf+"]");


  //===================================
  //
  // Save to root file 
  //
  //===================================
  TFile *outRootFile; 
  outRootFile = new TFile(outFile, "RECREATE");     
  BHistGlobal  -> Write("BHistGlobal");
  mHist        -> Write("mHist");  
  outRootFile  -> Close();

} 


