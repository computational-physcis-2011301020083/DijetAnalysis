#include "../../ATLASStyle/AtlasStyle.C"        
#include "../../ATLASStyle/AtlasUtils.C"
#include "../../code/helperFunctions/histPlotFormatFunctions.C"
#include "../../code/helperFunctions/statFitFunctions.C"                                                                                     

using namespace std; 
using namespace TMath;


TString DataPath = "../../output/root/data/dijet_gausRes_LLHScan.root"; 
TString PEPath   = "../../condor/gausRes/root"; 
int numPE        = 100 ; 

double lumi = 37.0   ; 
TString pdf = "./pdf/dijets_2015_2016_37fb_pValuesPE_gausRes.pdf"; 

void pValuesFromPE() 
{ 

  SetAtlasStyle ();

  //==// Make canvas 
  TCanvas *canvas = new TCanvas();
  canvas -> Update();
  canvas -> Print(pdf+"[");
  gPad   -> SetLeftMargin(0.12);
  gPad   -> SetBottomMargin(0.12);

  //==// Read data file and get TGraphs/TH1Ds
  TFile *dataFile = TFile::Open( DataPath, "r" ); 
  dataFile -> cd();
  TGraph *TG_data_LLHR          = (TGraph*)dataFile -> Get("LLHR_positive");  // LLH ratio for positive signals 
  TGraph *TG_data_LLHRpVal      = (TGraph*)dataFile -> Get("LLHRpVal");       // LLHR local p-values 
  TH1D   *TH1_data_swiftBkgBest = (TH1D*)dataFile   -> Get("swiftBkgBest");   // swift background from data 
  int numPoints                 = TG_data_LLHRpVal   -> GetN(); 

  //==// Find largest LLHR value in data 
  double maxLLHR_data      = TMath::MaxElement(TG_data_LLHR->GetN(), TG_data_LLHR->GetY()); 
  double minLLHRpVal_data  = TMath::MinElement(TG_data_LLHRpVal->GetN(), TG_data_LLHRpVal->GetY());
  double localpValSig_data = ROOT::Math::gaussian_quantile_c(minLLHRpVal_data, 1); // convert local p-value to sigma using quantile function 

  //==// Find mass that has the lowest p-value 
  double bestMass    = 0.0; 
  double x, y;  
  for (int i = 0; i < numPoints; i++) 
  {  
    int k = TG_data_LLHR -> GetPoint(i, x, y);
    if (y == maxLLHR_data) 
    { 
      bestMass = x; 
      break; 
    }
  }  
 
  //==// Now loop over all pseudo-experiments and save the largest LLHR value from each   
  vector<double> vec_PE_LLHR; 
  TH1D *TH1_PE_LLHR = new TH1D("PE_LLHR","", 50 , 0, 14 );
  for (int i = 0; i < numPE; i++)
  { 
    if (i%100==0) cout << "# PE: " << i << endl;

    //==// open PE file 
    TFile *PEFile       = TFile::Open( PEPath+"/PE_seed"+to_string(i+1)+"_LLHScan.root", "r" ); 
    if (PEFile==NULL) 
    {
      cout << "Skipping PE" << i+1 << endl;  
      continue; 
    }
    if (PEFile -> IsZombie())
    { 
      cout << "Skipping PE" << i+1 << endl;
      continue; 
    }
    PEFile              -> cd();  
 
    //==// get largest LLHR value 
    if ( !(PEFile->GetListOfKeys()->Contains("LLHR_positive")) ) 
    {   
      cout << "Skipping PE " << i+1 << endl; 
      continue;  
    } 

    TGraph *TG_PE_LLHR  = (TGraph*)PEFile -> Get("LLHR_positive");            // LLH ratio for positive signals  
    double maxLLHR_PE   = TMath::MaxElement(TG_PE_LLHR->GetN(), TG_PE_LLHR->GetY());  
    vec_PE_LLHR.push_back( maxLLHR_PE ); 
    TH1_PE_LLHR       -> Fill( maxLLHR_PE );
    PEFile -> Close(); 
  }    
  cout << "Using " << vec_PE_LLHR.size() << "/" << numPE << " pseudo-experiments. " << endl; 

  //==// COMPUTE GLOBAL P-VALUE
  //==// by counting how many times in the PEs we see a LLHR value larger than what we see in the data 
  double numPElargerThanData = 0.0; 
  for (int i = 0; i < vec_PE_LLHR.size(); i++)
    if ( vec_PE_LLHR[i] > maxLLHR_data ) 
      numPElargerThanData += 1; 

  double global_pValue = numPElargerThanData / vec_PE_LLHR.size();            // Compute global p-value 
  double globalpValSig = ROOT::Math::gaussian_quantile_c(global_pValue, 1);   // convert p-value to sigma using quantile function 


  //==// Print out results 
  cout << endl; 
  cout << "Results from "      << DataPath << endl; 
  cout << "------------------------------------------------------------------------------" << endl; 
  cout << "Mass           : "  << bestMass << " GeV" << endl; 
  cout << "Local p-value  : "  << minLLHRpVal_data   << " (" << localpValSig_data  << " sig)" << endl; 
  cout << "Global p-value : "  << global_pValue << " (" << globalpValSig << " sig)" << endl;  
  cout << endl;  

  //==// Format histogram and plot  
  formatHist(TH1_PE_LLHR, "# PE", "LLHR"); 
  TH1_PE_LLHR -> Draw("hist"); 
  drawText( 0.65, 0.85, "#bf{#it{ATLAS}} Internal", kBlack, 0.038);
  drawText( 0.65, 0.81, "#sqrt{s} = 13TeV, "+roundUp(lumi, 1) +"fb^{-1}", kBlack, 0.032);
  drawText( 0.65, 0.77, "Data", kRed, 0.032); 
  drawText( 0.65, 0.72, "Mass : "+roundUp(bestMass, 0)+"GeV", kBlack, 0.032);
  drawText( 0.65, 0.68, "Local p-value : "+roundUp(minLLHRpVal_data,3)+" ("+roundUp(localpValSig_data, 2)+"#sigma)", kBlack, 0.032);
  drawText( 0.65, 0.64, "Global p-value: "+roundUp(global_pValue,3)+" ("+roundUp(globalpValSig)+"#sigma)", kBlack, 0.032);
  makeTLine(maxLLHR_data, 0, maxLLHR_data, TH1_PE_LLHR->GetMaximum() , kRed, "same", 2, 1, 0); 
  canvas -> Print(pdf);  

  canvas -> Print(pdf+"]"); 




}
