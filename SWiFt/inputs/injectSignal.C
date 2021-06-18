///////////////////////////////////////////////////////////////
//
// This code injects signal into a distribution. 
// It also generates PEs from the signal injected distribution 
//
///////////////////////////////////////////////////////////////

#include "../ATLASStyle/AtlasStyle.C"                                                                                             
#include "../ATLASStyle/AtlasUtils.C"
#include "../code/helperFunctions/histPlotFormatFunctions.C"
#include "../code/helperFunctions/statFitFunctions.C"
#include "../code/helperFunctions/fitFunctions.C" 

using namespace std; 
using namespace TMath;


//====================
//
//  Inputs and outputs  
//
//====================

TString inputHist      = "smoothMC/mjj_pythia_NLO_EW_37fb_6paramFit.root";
TString inputHistName  = "BHistGlobal"; 
TString sigFile        = "../morphing/output/gaus_dijet/gaus_res.root"; 

TString pdf            = "smoothMC_signal/smoothMCNLOEW_37fb_6paramFit_2016GeV_gausRes_5000events.pdf"; 
TString outFile        = "smoothMC_signal/smoothMCNLOEW_37fb_6paramFit_2016GeV_gausRes_5000events.root"; 
TString outPEPath      = "smoothMC_signal/PE_6paramFit_2016GeV_5000events_gausRes/";

//============
//
// Set options
//
//============

double xminfit      = 1100.0  ; //GeV
double xmaxfit      = 9000.0  ; //GeV
double lumi         = 37.0    ; // fb-1 

TString SonlyFunc   = "GFunc" ; 
double sigMass      = 2016.0  ; // GeV 
double scale        = 1       ; 
int numSignalEvents = 5000    ; 

bool doAndSavePEs   = 0       ;
int numPE           = 1000    ; 


void injectSignal()
{ 

  //===================================
  //
  // Open input file and read histogram 
  //
  //===================================
  TFile *inputFile    = TFile::Open( inputHist, "r" );
  inputFile -> cd(); 
  TH1D *inputHist = (TH1D*)inputFile->Get(inputHistName);
  
  //===================================
  //
  // Open signal spline and generate signal shape  
  //
  //===================================

  //==// get signal TF1 
  TFile *signalFile       = openFile(sigFile);
  TF1 *sigFunc            = new TF1("sigFunc", getFunction(SonlyFunc), xminfit, xmaxfit);
  vector<double> SigParam = generateSignal(sigMass, signalFile, sigFunc);
 
  //==// convert TF1 to histogram 
  TH1D *sigHist  = (TH1D*) inputHist -> Clone("sigHist");
  fit_to_hist(sigHist, sigFunc, xminfit, xmaxfit);

  sigHist -> Scale(numSignalEvents);
  cout << "Signal Events injected: " << sigHist -> Integral() << endl;  

  //===================================
  //
  // Inject and set errors 
  //
  //===================================
  
  //==// inject signal into input hist 
  inputHist -> Add(sigHist);
  
  //==// set sqrtN error 
  for (int i = 1; i <= inputHist->GetNbinsX(); i++) 
    inputHist -> SetBinError( i, sqrt(inputHist->GetBinContent(i)) );

 
  //===================================
  //
  // Inject and set errors 
  //
  //===================================
  
  if (doAndSavePEs)
  { 
    for (int i = 1; i<=numPE; i++ )
    {
      TFile *outPEFile; 
      outPEFile = new TFile( outPEPath+"PE"+to_string(i)+".root"  , "RECREATE");     
      TH1D *PEHist = doPseudoExp(inputHist, i);
      PEHist -> SetLineColor(kBlack);
      PEHist -> Write("myy");
      outPEFile  -> Close();
    }
  } 
  
  
  //===================================
  //
  // Save output
  //
  //===================================

  //==// Save final output as root file 
  TFile *outRootFile; 
  outRootFile = new TFile(outFile, "RECREATE");     
  sigHist     -> Write("injectedSignal"); 
  inputHist   -> Write("histWithInjectedSignal");
  outRootFile -> Close();

} 


