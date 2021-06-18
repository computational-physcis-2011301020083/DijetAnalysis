//////////////////////////////////////////////////////////////////////////////
//
// Code for fitting signal samples to a function: gaussian + reverse landau 
// By: Karishma Sekhon (ksekhon@umich.edu)
//
////////////////////////////////////////////////////////////////////////////////

#include "functions.C"
#include "sigFunctions.C"
#include "config.C"

using namespace std; 


////////////// Main function 
void morphMCSignal(TString cfgName) 
{ 

  TVirtualFitter::SetMaxIterations(1000000);

  //====// load config settings 
  loadConfig(config, cfgName);

  //====// counter for number of signal masses to morph 
  int numMC = config.HNAMES.size(); 

  //====// make arrays for looping over MC signal shapes  
  vector<TFile*> FILES(numMC);
  vector<TH1D*>  HISTS(numMC); 
  vector<TFitResultPtr> FITS(numMC);  

  //====// make canvas and set canvas style
  TCanvas *canvas = new TCanvas();
  canvas -> Update();
  gStyle -> SetOptFit();  // show stats on plots
  gStyle -> SetOptStat();
  gPad   -> SetLeftMargin(0.15);
  gPad   -> SetBottomMargin(0.15);
  TString pdfFile = config.path+"/"+config.sigName+config.pdf; 
  canvas -> Print(pdfFile+"[");
  canvas -> SetLogy(); 

  //====// vectors to store fit parameters in 
  vector<double> Param0;
  vector<double> Param1;
  vector<double> Param2;
  vector<double> Param3;
  vector<double> Param4;
  vector<double> Param5;
  vector<double> Param6;
  vector<double> Param7;
  vector<double> Param8;
  vector<double> MassNum; 
  vector< vector<double> > Parameters; 
  vector<TF1*> Fits; 
  double MinCSq = 0.0 ; 

  //====// loop over all the MC signal masses
  for (int i = 0; i<numMC; i++)
  //for (int i = 8; i < 10; i++)
  {
    cout << "Runing over histogram: " << config.HNAMES[i] << endl; 

    //====// read from individual files 
    /*
    FILES[i] =  openFile(config.path+"/"+config.FNAMES[i]); 
    TDirectoryFile *d = NULL; 
    d = (TDirectoryFile*)FILES[i]->Get("Nominal");                                                                         
    HISTS[i] = (TH1D*)d->Get(config.HNAMES[i]);
    */    

    //====// read from common file 
    if (config.oneInFile)
      FILES[i] =  openFile(config.fileName);
    else 
      FILES[i] =  openFile(config.FNAMES[i]);
    HISTS[i] =  getHisto(FILES[i], config.HNAMES[i] ); 
    

    if (!config.quietMode) HISTS[i] -> ResetBit(TH1::kNoStats);
    HISTS[i] -> GetYaxis() -> SetTitle("Events"); 
    HISTS[i] -> Rebin(config.rebin);
    HISTS[i] -> GetXaxis()->SetRangeUser(config.histMin, config.histMax);;    


    //====// take care of variable bin-widths 
    for (int j = 1; j <= HISTS[i]->GetNbinsX(); j++)
    { 
      HISTS[i]->SetBinContent(j, HISTS[i]->GetBinContent(j)/HISTS[i]->GetBinWidth(j) ); 
      HISTS[i]->SetBinError(j, HISTS[i]->GetBinError(j)/HISTS[i]->GetBinWidth(j) );
    } 

    int binmax = HISTS[i]->GetMaximumBin();
    double x = HISTS[i]->GetXaxis()->GetBinCenter(binmax);


    //====// define function to fit
    double fitLow = config.fitMin;  
    TF1 *f1 = new TF1((to_string((int)config.MASS[i])+"fit").c_str(), getFunction(config.sigFunc) ,fitLow,config.fitMax);
    

    //====// initialize fit parameters 
    if (config.paramOpt[i]) // == 1: use generalized parameters 
    {  
      f1       -> SetParameter(0, HISTS[i] -> Integral()*config.rebin); 
      f1       -> SetParameter(1, x - 150.);
      f1       -> SetParameter(2, HISTS[i] -> GetRMS());
      f1       -> SetParameter(3, 0.4);
      f1       -> SetParameter(4, -x); 
      f1       -> SetParameter(5, HISTS[i] -> GetRMS()/3.); 
    } 
    else //====// == 0: use hand tuned parameters 
    { 
      f1       -> SetParameter(0, HISTS[i] -> Integral()*config.rebin );   
      f1       -> SetParameter(1, config.PAR1[i]  );
      f1       -> SetParameter(2, config.PAR2[i]  );
      f1       -> SetParameter(3, config.PAR3[i]  );
      f1       -> SetParameter(4, config.PAR4[i]  );  
      f1       -> SetParameter(5, config.PAR5[i]  ); 
    }
    

    //====// fit to the histograms and store fit results 
    vector<double> fitQuality;
    FITS[i] = HISTS[i] -> Fit(f1, "SRQL"); // Do the first fit to signal model

    //====// find better fits with autoFudge 
    if (config.autoFudge)
      autoFudge(FITS[i], HISTS[i], f1, config.refit, config.delta );
    FITS[i] = HISTS[i] -> Fit(f1, "SRQ");

    //====// print out final parameters 
    cout << config.MASS[i] << endl; 
    cout << f1->GetParameter(0) << endl; 
    cout << f1->GetParameter(1) << endl; 
    cout << f1->GetParameter(2) << endl; 
    cout << f1->GetParameter(3) << endl; 
    cout << f1->GetParameter(4) << endl; 
    cout << f1->GetParameter(5) << endl; 
    cout << endl;

    //====// draw fits to pdf 
    HISTS[i] -> Draw();   
    drawText(0.80,0.20,to_string((int)config.MASS[i])+"GeV"); // 0.80, 0.40
    canvas   -> Print(pdfFile); 

    //====// save parameters in vectors 
    MassNum.push_back(config.MASS[i]);
    Param0.push_back(f1->GetParameter(0));
    Param1.push_back(f1->GetParameter(1));
    Param2.push_back(abs(f1->GetParameter(2)));
    Param3.push_back(f1->GetParameter(3));
    Param4.push_back(f1->GetParameter(4));
    Param5.push_back(f1->GetParameter(5));   
    
    Fits.push_back(f1); 
  }

  ////////////////////////////// display fit paramters on TGraphs

  TGraph* TG_par0 = makeTGraph(MassNum, Param0, "Normalization");  
  canvas  -> SetLogy(); 
  TG_par0 -> Draw("AP");
  TSpline3 *s0 = new TSpline3("",TG_par0);                                                                                                               
  s0      -> SetLineColor(kRed);
  s0      -> Draw("c same");
  canvas  -> Print(pdfFile);
  canvas  -> SetLogy(0);

  TGraph* TG_par1 = makeTGraph(MassNum, Param1, "Gaussian Mean");  
  TG_par1 -> Draw("AP");
  //TG_par1 -> SetMinimum(0); 
  //TG_par1 -> SetMaximum(6000); 
  TSpline3 *s1 = new TSpline3("",TG_par1);                                                                                                               
  s1      ->SetLineColor(kRed);
  s1      ->Draw("c same");
  canvas  -> Print(pdfFile);

  TGraph* TG_par2 = makeTGraph(MassNum, Param2, "Gaussian Width");  
  TG_par2 -> Draw("AP");
  //TG_par2 -> SetMinimum(-1000);
  //TG_par2 -> SetMaximum(2000); 
  TSpline3 *s2 = new TSpline3 ("",TG_par2);                                                                                                               
  s2      -> SetLineColor(kRed);
  s2      -> Draw("c same");
  canvas  -> Print(pdfFile);

  TGraph* TG_par3 = makeTGraph(MassNum, Param3, "G/L Fraction");  
  TG_par3 -> Draw("AP");
  TG_par3 -> SetMinimum(-1); 
  TG_par3 -> SetMaximum(2); 
  TSpline3 *s3 = new TSpline3("",TG_par3);                                                                                                               
  s3      -> SetLineColor(kRed);
  s3      -> Draw("c same");
  canvas  -> Print(pdfFile);

  TGraph* TG_par4 = makeTGraph(MassNum, Param4, "Landau Mean");  
  TG_par4 -> Draw("AP");
  TG_par4 -> SetMinimum(-7000); 
  TG_par4 -> SetMaximum(0); 
  TSpline3 *s4 = new TSpline3("",TG_par4);                                                                                                               
  s4      -> SetLineColor(kRed);
  s4      -> Draw("c same");
  canvas  -> Print(pdfFile);

  TGraph* TG_par5 = makeTGraph(MassNum, Param5, "Landau Width");  
  TG_par5 -> Draw("AP");
  TG_par5 -> SetMinimum(0); 
  TG_par5 -> SetMaximum(400); 
  TSpline3 *s5 = new TSpline3("",TG_par5);                                                                                                               
  s5      -> SetLineColor(kRed);
  s5      -> Draw("c same");
  canvas  -> Print(pdfFile);

  //====// overlay fits 
  canvas -> SetLogy(); 
  for (int i = 0; i<Fits.size(); i++)
  { 
    //Fits[i] -> FixParameter(0, 1.0);
    if (i==0) { 
      Fits[i] -> Draw();
      Fits[i] -> SetMinimum(1e-9);  
    }
    else 
      Fits[i] -> Draw("same");
  }
  canvas -> Print(pdfFile);  
  canvas -> SetLogy(0); 

  canvas   -> Print(pdfFile+"]");


  ////////////////////////////// save interpolations
  /*
  TFile *oRootFile = new TFile(config.path+"/"+config.sigName+config.oRoot, "RECREATE");
  s0 -> Write("s0"); 
  s1 -> Write("s1"); 
  s2 -> Write("s2");
  s3 -> Write("s3");
  s4 -> Write("s4");
  s5 -> Write("s5"); 
  TG_par0 -> Write("TG_par0");
  TG_par1 -> Write("TG_par1");
  TG_par2 -> Write("TG_par2");
  TG_par3 -> Write("TG_par3");
  TG_par4 -> Write("TG_par4");
  TG_par5 -> Write("TG_par5");
  oRootFile -> Close(); 
  */
  

}
