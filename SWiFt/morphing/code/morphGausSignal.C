//////////////////////////////////////////////////////////////////////////////
//
// Code for fitting signal samples to a function: gaussian + reverse landau 
// By: Karishma Sekhon (ksekhon@umich.edu)
//
////////////////////////////////////////////////////////////////////////////////

#include "functions.C"
#include "config.C"

using namespace std; 


////////////// Main function 
void morphGausSignal(TString cfgName) 
{ 

  //==// load config settings 
  loadConfig(config, cfgName);

  //==// make canvas and set canvas style
  TCanvas *canvas = new TCanvas();
  canvas -> Update();
  gStyle -> SetOptFit();  // show stats on plots
  gStyle -> SetOptStat();
  gPad   -> SetLeftMargin(0.15);
  gPad   -> SetBottomMargin(0.15);
  TString pdfFile = config.path+"/"+config.sigName+config.pdf;
  canvas -> Print(pdfFile+"[");

  //==// make vector of masses 
  vector<double> MassNum; 
  vector<double> Norm; 
  for (double d = config.massMin ; d <= config.massMax; d=d+50)
  { 
    Norm.push_back(1.0); 
    MassNum.push_back(d); 
  } 

  //==// read in resolution function 
  TF1 *resFun = new TF1("Resolution", config.resFunc ,config.massMin, config.massMax);                                                                   
  resFun -> SetNpx(500);
  for (int i = 0; i<config.par.size(); i++) 
    resFun    -> FixParameter(i, config.par[i]); 
  resFun -> Draw(); 
  canvas -> Print(pdfFile);

  //==// calculate width 
  vector<double> Width; 
  for (int i = 0; i < MassNum.size(); i++) 
    Width.push_back( MassNum[i] * resFun->Eval(MassNum[i]) ); 
 
  //==// display fit paramters on TGraphs
  TGraph* TG_par0 = makeTGraph(MassNum, Norm, "Normalization");  
  TG_par0 -> Draw("AP");
  TSpline3 *s0 = new TSpline3("",TG_par0);                                                                                                               
  s0      -> SetLineColor(kRed);
  s0      -> Draw("c same");
  canvas  -> Print(pdfFile);

  TGraph* TG_par1 = makeTGraph(MassNum, MassNum, "Mean");  
  TG_par1 -> Draw("AP");
  TSpline3 *s1 = new TSpline3("",TG_par1);                                                                                                               
  s1      -> SetLineColor(kRed);
  s1      -> Draw("c same");
  canvas  -> Print(pdfFile);

  TGraph* TG_par2 = makeTGraph(MassNum, Width, "Width");  
  TG_par2 -> Draw("AP");
  TSpline3 *s2 = new TSpline3("",TG_par2);                                                                                                               
  s2      -> SetLineColor(kRed);
  s2      -> Draw("c same");
  canvas  -> Print(pdfFile);

  //==// close pdf 
  canvas  -> Print(pdfFile+"]");

  //==// save interpolations
  TFile *oRootFile = new TFile(config.path+"/"+config.sigName+config.oRoot, "RECREATE");
  s0 -> Write("s0");
  s1 -> Write("s1");
  s2 -> Write("s2"); 
  TG_par0 -> Write("TG_par0");
  TG_par1 -> Write("TG_par1");
  TG_par2 -> Write("TG_par2");
  oRootFile -> Close(); 
  

}
