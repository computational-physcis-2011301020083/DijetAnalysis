///////////////////////////////////////////////////////////////
//
// This code injects signal into a distribution. 
// It also generates PEs from the signal injected distribution 
//
///////////////////////////////////////////////////////////////

#include "../../code/helperFunctions/histPlotFormatFunctions.C"
#include "../../code/helperFunctions/statFitFunctions.C"

using namespace std; 
using namespace TMath;


void doPE(int numPE, TString inputFile, TString inputHistName, TString outPEPath)
{ 

  //==// Open input file and read histogram 
  TFile *file    = TFile::Open( inputFile, "r" );
  file -> cd(); 
  TH1D *inputHist = (TH1D*)file->Get(inputHistName);
  
  //==// Do PE and save  
  for (int i = 1; i<=numPE; i++ )
  {
    TFile *outPEFile; 
    outPEFile = new TFile( outPEPath+"/PE"+to_string(i)+".root"  , "RECREATE");     
    TH1D *PEHist = doPseudoExp(inputHist, i);
    PEHist -> SetLineColor(kBlack);
    PEHist -> Write("mjj");
    outPEFile  -> Close();
  }
  
} 


