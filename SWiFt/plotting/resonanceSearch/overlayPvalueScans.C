
#include "../../ATLASStyle/AtlasStyle.C"    
#include "../../ATLASStyle/AtlasUtils.C"
#include "../../code/helperFunctions/histPlotFormatFunctions.C"
#include "../../code/helperFunctions/statFitFunctions.C"

using namespace std; 
using namespace TMath;

TString filesPath = "../../output/root/data"; 
TString files[] = 
{
"dijet_gausRes_LLHScan.root", 
//"dijet_gaus5per_LLHScan.root",
}; 


TString legendLabel[] = 
{
"Res. Gaussian", 
//"5% Gaussian", 
};

int col[] = 
{
kTeal, 
//kTeal-1, 
};

double low  = 1416.0 ; 
double high = 6918.0 ; 
TString pdf = "./pdf/dijets_2015_2016_37fb.pdf"; 

void overlayPvalueScans()
{ 

  SetAtlasStyle ();

  //==// Make canvas 
  TCanvas *canvas = new TCanvas();
  canvas -> Update();
  canvas -> Print(pdf+"[");
  gPad   -> SetLeftMargin(0.12);
  gPad   -> SetBottomMargin(0.12);

  //==// Read in root files  
  int num = sizeof(files)/sizeof(files[0]); 
  TFile *file;
  TGraph *TG_signal[num];
  TGraph *TG_signifi[num];
  TGraph *TG_LLHR[num];

  for (int i = 0; i<num; i++)
  { 
    file        = TFile::Open( filesPath+"/"+files[i], "r" );                                                                                   
    file        -> cd(); 
    TG_signal[i]          = (TGraph*)file->Get("ExtractedSignal");
    TG_signifi[i]         = (TGraph*)file->Get("Significance");
    TG_LLHR[i]            = (TGraph*)file->Get("LLHRpVal");
 
  }


  ////////////////////////////////////////////////////////
  //
  // Now draw 
  //
  ////////////////////////////////////////////////////////


  //==// LLHR p-value
  canvas -> SetLogy(); 
  for (int i = 0; i<num; i++) 
  {   
    TG_LLHR[i] -> SetMarkerColor(col[i]); 
    TG_LLHR[i] -> SetMaximum(1);
    TG_LLHR[i] -> SetMinimum(9e-4);
    formatTGraph(TG_LLHR[i], "M_{jj} [GeV]", "Local p-value", col[i]); 
    if (i == 0)
      TG_LLHR[i]-> Draw("APL");
    else                                                                                                                    
      TG_LLHR[i]-> Draw("PL");
  }
  drawText(0.70, 0.45,"#bf{#it{ATLAS}} Internal", kBlack, 0.04);
  drawText(0.70, 0.41,"#sqrt{s} = 13TeV, 37fb^{-1}", kBlack, 0.032);
  for (int i = 0; i < num; i++)
    drawText(0.70, 0.36-(i*0.04), legendLabel[i], col[i], 0.032);
  
  makeTLine(low, Prob(1*1, 1)/2 , high , Prob(1*1, 1)/2, 17, "", 2, 2, 1);
  makeTLine(low, Prob(2*2, 1)/2 , high , Prob(2*2, 1)/2, 17, "", 2, 2, 1);
  makeTLine(low, Prob(3*3, 1)/2 , high , Prob(3*3, 1)/2, 17, "", 2, 2, 1); 
  makeTText(high+50, Prob(1*1, 1)/2, "1s", kBlack, 1, 0.030);
  makeTText(high+50, Prob(2*2, 1)/2, "2s", kBlack, 1, 0.030);
  makeTText(high+50, Prob(3*3, 1)/2, "3s", kBlack, 1, 0.030);

  canvas -> Print(pdf); canvas -> Clear(); 
  canvas -> SetLogy(0); 






  canvas->Print(pdf+"]");
} 
