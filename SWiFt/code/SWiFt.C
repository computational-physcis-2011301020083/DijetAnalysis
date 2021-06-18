////////////////////////////////////////////////////////////////////////////////                                                                        
// 
// Sliding Window Fit (SWiFt) 
//   
////////////////////////////////////////////////////////////////////////////////

#include "TError.h"
#include "helperFunctions/fitFunctions.C"
#include "helperFunctions/histPlotFormatFunctions.C"
#include "helperFunctions/statFitFunctions.C"
#include "helperFunctions/windowSelectorBins.C"
#include "helperFunctions/swiftBkgGenerator.C"
#include "config.C"
#include "../ATLASStyle/AtlasStyle.C"                                                                                                                  
#include "../ATLASStyle/AtlasUtils.C"

using namespace std; 
using namespace TMath; 


void SWiFt(TString cfgName = "") 
{ 

  //======================================================== 
  //========================================================                                                                                                   
  //  
  // PART O: The preperation for SWiFt  
  //   -> Load input files (input histogram, signal files)
  //   -> Prepare output pdf, root file
  //   -> Perform two global fits 
  //   -> Initialize quantities for SWiFt's slide 
  //========================================================
  //========================================================
    
  SetAtlasStyle ();

  cout << endl; 
  cout << "-------------------------------------------------------------------" << endl;
  cout << "------------> Welcome to SWiFt ... " << endl; 
  cout << "-------------------------------------------------------------------" << endl;
  cout << endl; 
  double value1,value2,value3;
  printStuff("../ascii/dragon.txt"); 

  //==// Suppress all the annoying root messages
  gErrorIgnoreLevel = kFatal;

  //==// Load config
  loadConfig(config, cfgName);
  cout << "Config File  : " << cfgName << endl;

  //==// Set up output files 
  TString rootFileName  = config.path+"/"+config.oRootName+".root";    
  TString pdfFileName   = config.path+"/"+config.pdfName+".pdf"; 

  cout << "Output files : " << endl; 
  cout << "          PDF: " << pdfFileName  << endl; 
  cout << "         ROOT: " << rootFileName << endl; 

  TFile* rootFile  = new TFile( rootFileName, "recreate" );  
  TString pdf      = pdfFileName;                                                

  //==// Prepare TCanvas for printing out pdfs 
  TCanvas *canvas = new TCanvas();
  canvas -> Update();
  gStyle -> SetOptFit();
  canvas -> Print(pdf+"[");
  gPad   -> SetLeftMargin(0.12);
  gPad   -> SetBottomMargin(0.12);
  
  //======================================================                                                                                                    
  //  
  // Load input histogram 
  //  
  //======================================================

  cout << endl; 
  cout << "-------------------------------------------------------------------" << endl;
  cout << "------------> Opening input histogram" << endl; 
  cout << "-------------------------------------------------------------------" << endl;
  cout << endl; 
  
  //==// open input file
  TFile *dFile      = openFile(config.inputFile); 
  dFile -> cd();
  
  //==// read data histogram from root file 
  TH1D *histRaw    = NULL;
  TDirectoryFile *d = NULL;

  if ( config.inputDir == "x") 
    histRaw = getHisto(dFile, config.inputHist); 
  else  
    histRaw = getHisto(dFile, d, config.inputDir, config.inputHist); 

  //==// format and draw histogram
  histRaw -> ResetBit(TH1::kNoStats); // show fit parameters
  formatHist(histRaw, "Events", "m_{jj} [GeV]");
  histRaw -> Rebin(config.rebin);  
  histRaw -> GetXaxis() -> SetRangeUser(config.xminfit, config.xmaxfit);
  cout << "Data Hist Integral: "<< histRaw -> Integral() << endl;
  canvas -> SetLogy();
  histRaw -> Draw("e1"); 
  drawText( 0.70, 0.85, "#bf{#it{ATLAS}} Internal", kBlack, 0.038);
  drawText( 0.70, 0.80, "#sqrt{s} = 13TeV, "+roundUp(config.lumi, 1) +"fb^{-1}", histRaw->GetLineColor(), 0.032);
  canvas -> Print(pdf); canvas -> Clear();  
  canvas -> SetLogy(0); 

  //======================================================                                                                                                    
  //  
  // Record bin edges and create bin-width divided histogram
  //  
  //======================================================
  cout << endl; 
  cout << "-------------------------------------------------------------------" << endl;
  cout << "------------> Record bin edges and divide bin-width from input"   << endl; 
  cout << "-------------------------------------------------------------------" << endl;
  cout << endl; 

  //==// store bin edges for window centers and signal points
  vector<double> v_binLowEdge; 
  for (int i = 0; i<histRaw->GetNbinsX(); i++ ) {   
    v_binLowEdge.push_back( histRaw->GetBinLowEdge(i) );
    cout << histRaw->GetBinLowEdge(i) << endl; 
  }
  cout << "Number of bins: " << v_binLowEdge.size() << " [ "<< v_binLowEdge[0] << ", " << v_binLowEdge[v_binLowEdge.size()-1] << " ]"<<  endl; 

  //==// clone input histogram and divide out bin-width 
  TH1D *histDivBinWidth = (TH1D*)histRaw -> Clone();
  histDivBinWidth -> GetYaxis() -> SetTitle("Events/GeV"); 

  for (int i = 1; i <= histDivBinWidth->GetNbinsX(); i++)
  {   
    histDivBinWidth -> SetBinContent(i, histRaw -> GetBinContent(i) / histRaw -> GetBinWidth(i) );  
    histDivBinWidth -> SetBinError(i, histRaw -> GetBinError(i) / histRaw -> GetBinWidth(i) );
  }  
  cout << "Divided out bin-width from clone of input histogram" << endl;  



  //======================================================                                                                                                    
  //  
  // Create histograms that will store the SWiFt bkg
  //  
  //======================================================
  cout << endl; 
  cout << "-------------------------------------------------------------------" << endl;
  cout << "------------> Make histograms to store nominal/best SWiFt bkg" << endl; 
  cout << "-------------------------------------------------------------------" << endl;
  cout << endl;

  TH1D *bkgNomHist  = (TH1D*)histRaw->Clone(); 
  TH1D *bkgBestHist = (TH1D*)histRaw->Clone();
  TH1D *bkgSSHist   = (TH1D*)histRaw->Clone();
  bkgNomHist        -> SetLineColor(kBlue);
  bkgBestHist       -> SetLineColor(kGreen+2);
  bkgSSHist         -> SetLineColor(kGreen); 

  //======================================================                                                                                                    
  //  
  // Load signal files - they should contain the spline interpolation results
  //  
  //======================================================
  cout << endl; 
  cout << "-------------------------------------------------------------------" << endl;
  cout << "------------> Opening signal spline interpolation file" << endl; 
  cout << "-------------------------------------------------------------------" << endl;
  cout << endl;

  //==// read in nominal signal path 
  cout << "Nominal signal splines ... " << endl; 
  TFile *nomSigFile = openFile(config.nomSigFile);

  //==// read in JES varied signal paths 
  vector<TFile* > JES1upSigFile; 
  vector<TFile* > JES1downSigFile;  

  if (config.addSys)
  {
    cout << endl;  
    cout << "JES signal systematics are requested. Loading TFiles ...  " << endl; 
    for (int i = 0; i < config.JES1upSigFile.size(); i++) 
    { 
      JES1upSigFile.push_back( openFile(config.JES1upSigFile[i]) );
      JES1downSigFile.push_back( openFile(config.JES1downSigFile[i]) );
    } 
  }     
 
  //======================================================                                                                                                    
  //  
  // Perform global fits - these will initialize the 1st window fits
  //  
  //======================================================
  
  //*
  //* NOMINAL FUNCTION FIT 
  //*
  cout << endl; 
  cout << "-------------------------------------------------------------------" << endl;  
  cout << "------------> Performing global background-only fit: NOMINAL "   << endl; 
  cout << "-------------------------------------------------------------------" << endl;

  //==// set up global function 
  vector<double> nomGlob_fitQuality; 
  TF1 *BFuncGlobal = new TF1("BFuncGlobal", getFunction(config.bkgFunc) ,config.xminfit, config.xmaxfit);  
  formatFunc(BFuncGlobal, kBlue-2, 2);
  for (int i = 0; i<config.par.size(); i++) 
    BFuncGlobal    -> SetParameter(i, config.par[i]);

  //==// fit 
  double BFuncGlobalLLH = minuitFit ( histRaw, BFuncGlobal, 0, config.xminfit, config.xmaxfit, config.quietFitInfo, &nomGlob_fitQuality); 
  if ( nomGlob_fitQuality[0] || nomGlob_fitQuality[1] != 0 ) 
  { 
    cout << endl; 
    printStuff("../ascii/crab.txt");
    cout << "NOMINAL global fit is failing ... retry with different initial conditions in the config. " << endl; 
    cout << endl; 
    exit(0); 
  }

  //==// convert fit to histogram 
  TH1D *BHistGlobal   = (TH1D*) histRaw -> Clone("BHistGlobal");
  fit_to_hist(BHistGlobal, BFuncGlobal, config.xminfit, config.xmaxfit); 

  //==// plot
  double chisqpVal = ChiSqNDFTest(histRaw, BHistGlobal, config.xminfit, config.xmaxfit, BFuncGlobal->GetNDF(), 0, 0, 1);
  cout<<"CHI2 p-value of nominal global fit bkg "<<chisqpVal<<endl;
  value1=chisqpVal;
  plotBkg (canvas, pdf, histRaw, BHistGlobal, "Data", to_string(config.par.size())+" Paramter Fit", chisqpVal );

  //*
  //* ALTERNATE FUNCTION FIT
  //*
  cout << endl; 
  cout << "-------------------------------------------------------------------" << endl;  
  cout << "------------> Performing global background-only fit: ALTERNATE "   << endl; 
  cout << "-------------------------------------------------------------------" << endl;

  //==// set up global function 
  vector<double> altGlob_fitQuality;
  TF1 *BFunc2Global = new TF1("BFunc2Global", getFunction(config.bkgFunc2) ,config.xminfit, config.xmaxfit);  
  formatFunc(BFunc2Global, kBlue-5, 2);
  for (int i = 0; i<config.par2.size(); i++) 
    BFunc2Global    -> SetParameter(i, config.par2[i]);

  //==// fit 
  double BFunc2GlobalLLH = minuitFit ( histRaw, BFunc2Global, 0, config.xminfit, config.xmaxfit, config.quietFitInfo, &altGlob_fitQuality);
  if ( altGlob_fitQuality[0] || altGlob_fitQuality[1] != 0 ) 
  { 
    cout << endl; 
    printStuff("../ascii/crab.txt");
    cout << "ALTERNATE global fit is failing ... retry with different initial conditions in the config. " << endl; 
    cout << endl; 
    exit(0); 
  }

  //==// convert fit to histogram 
  TH1D *BHist2Global   = (TH1D*) histRaw -> Clone("BHist2Global");
  fit_to_hist(BHist2Global, BFunc2Global, config.xminfit, config.xmaxfit); 

  //==// plot
  chisqpVal = ChiSqNDFTest(histRaw, BHist2Global, config.xminfit, config.xmaxfit, BFunc2Global->GetNDF(), 0, 0, 1);
  cout<<"CHI2 p-value of alternate global fit bkg "<<chisqpVal<<endl;
  value2=chisqpVal;
  plotBkg (canvas, pdf, histRaw, BHist2Global, "Data", to_string(config.par2.size())+" Paramter Fit", chisqpVal );
 

  //======================================================                                                                                                    
  //  
  // Initialize parameters and settings for the SWiFT's slide 
  //  
  //======================================================
  
  cout << endl;
  cout << "-------------------------------------------------------------------" << endl;
  cout << "------------> Initialize TF1s. vectors, canvas, parameters " << endl;
  cout << "-------------------------------------------------------------------" << endl;
  cout << endl;

  cout << "... prepare vectors to store output variables" << endl; 
  
  //==// swift scan vectors
  vector<double> v_MassNum;  
  vector<double> v_Significance;
  vector<double> v_ExtractSig; 
  vector<double> v_ExtractSigError;
  vector<double> v_ExtractSigRaw; 
  vector<double> v_WindowBkg;
  vector<double> v_ExtractBkg; 
  vector<double> v_LogLikeHood;  
  vector<double> v_LogLikeHoodPvalue; 
  vector<double> v_Signal95CL; 
  vector<double> v_Limit95CL;
  vector<double> v_Chisquare_SB;  
  vector<double> v_ChisquareNdf_SB; 
  vector<double> v_ChisquarepVal_SB; 
  vector<double> v_LogLikeHood_SB;
  vector<double> v_Chisquare_B;  
  vector<double> v_ChisquareNdf_B; 
  vector<double> v_ChisquarepVal_B;
  vector<double> v_LogLikeHood_B; 
  vector<double> v_WindowCorr;
  vector<double> v_numParamChoosen; 

  //==// parameters in each window
  vector<vector<double>> vv_BkgParam_Bnom( config.par.size() );
  vector<vector<double>> vv_BkgParam_Balt( config.par2.size() );
  vector<vector<double>> vv_BkgParam_SBnom( config.par.size() );
  vector<vector<double>> vv_BkgParam_SBalt( config.par2.size() );

  //==// vectors for fit errors
  vector<vector<double>> vv_fitQuality_BnomOnly( v_binLowEdge.size() ); 
  vector<vector<double>> vv_fitQuality_Bnom( v_binLowEdge.size() ); 
  vector<vector<double>> vv_fitQuality_Balt( v_binLowEdge.size() ); 
  vector<vector<double>> vv_fitQuality_SBnom( v_binLowEdge.size() ); 
  vector<vector<double>> vv_fitQuality_SBalt( v_binLowEdge.size() );

  //==// vectors for pulls
  vector<vector<double>> vv_NP_values( config.numFlatNP+config.numJESNP ); 
  vector<vector<double>> vv_NP_errors( config.numFlatNP+config.numJESNP ); 
  vector<vector<double>> vv_NP95_values( config.numFlatNP+config.numJESNP ); 
  vector<vector<double>> vv_NP95_errors( config.numFlatNP+config.numJESNP );

  //==// nominal bkg-only fits while window selection 
  vector<double> v_bestWindow; 
  vector<double> v_windowLow; 
  vector<double> v_windowHigh;
  vector<double> v_ChisquarepVal_Bnom; 
  double currentWindowSize = 0.0;
  int currentWindowIdx     = -1; 


  //==// initialize parameters  
  cout << "... initialize 1st window fit from global fit parameters " << endl; 
  vector<double> paramNomGlob, paramNomBonly, paramNomSB; 
  for (int i = 0; i<config.par.size(); i++) // nominal fit  
  {
    paramNomGlob.push_back(BFuncGlobal->GetParameter(i));  
    paramNomBonly.push_back(paramNomGlob[i]);
  }

  vector<double> paramAltGlob, paramAltBonly;
  for (int i = 0; i<config.par2.size(); i++) // alternate fit  
  {
    paramAltGlob.push_back(BFunc2Global->GetParameter(i)); 
    paramAltBonly.push_back(paramAltGlob[i]); 
  }  

  //==// get the bins for all the important mass points  
  cout << "... set window conditions for the slide " << endl; 
  int start_bin    = histRaw -> FindBin(config.start);                                                                                                                
  int end_bin      = histRaw -> FindBin(config.end); 
  int xminfit_bin  = histRaw -> FindBin(config.xminfit);
  int xmaxfit_bin  = histRaw -> FindBin(config.xmaxfit);
  
  cout << endl; 
  cout << "WINDOW SETTINGS:::::::::::::::::::::::::::: " << endl; 
  cout << "First window center: "  << config.start   << " (bin " << start_bin   << ")" << endl; 
  cout << "Last window center : "  << config.end     << " (bin " << end_bin     << ")" << endl;
  cout << "SWiFt ultimate range: " << config.xminfit << " (bin " << xminfit_bin << ")" << " - " 
                                   << config.xmaxfit << " (bin " << xmaxfit_bin << ")" << endl; 
  cout << endl; 

  //======================================================== 
  //========================================================                                                                                                   
  //  
  // PART A: Pick window sizes 
  // This is done by trying various window sizes and picking
  // one that gives the best nominal bkg-only fit. 
  //  
  //========================================================
  //========================================================
    
  cout << endl;
  cout << "-------------------------------------------------------------------" << endl;
  cout << "------------> Pick window size "<< endl;
  cout << "-------------------------------------------------------------------" << endl;
  cout << endl;

  printStuff("../ascii/owl.txt"); 
  cout << endl; 
 
  //==// initialize window center and edges 
  double windowCenter = 0.0; 
  double lowEdge      = 0.0; 
  double highEdge     = 0.0; 

  /////////////////////////////////////////////////////// 
  //
  //== loop over bin edges 
  //
  /////////////////////////////////////////////////////// 
  for (int winIdx = start_bin; winIdx <= end_bin; winIdx++)
  {

    //==// define temporary variables needed in each window 
    TH1D *BHistBkg; 
    double bestpVal    = 0.0; 
    double bestChisq   = 0.0; 
    double bestWindow  = 0.0; 
    int windowNDF      = 0; 
    int bestWindowIdx  = 0; 
    int numGlobalFits  = 0;
    double lowEdgeTmp  = 0.0; 
    double highEdgeTmp = 0.0;
    vector<double> bbest_fitQuality; 

    //==// Set window center around which a window size will be picked 
    windowCenter  = v_binLowEdge[winIdx]; 
    v_MassNum.push_back( windowCenter );
    cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Window Center: " << windowCenter << endl;  
  
    //==// Fill vector with window sizes to test  
    vector<double> windowWidth; 
    double step = (config.windowWidthHigh - config.windowWidthLow) / config.windowNum ; 
    for (int n = 0; n <= config.windowNum; n++) 
      windowWidth.push_back(config.windowWidthLow + (n*step));

    //==// Set indexes of windows to loop over 
    // this is done such that for the 1st window center, all the window sizes in vector windowWidth are scanned 
    // for all other window centers, a total of config.windowNumAlways windows are scanned - this allows the window size to decrease or increase in a controlled manner
    // the config.windowNumAlways windows are centered around the previous window center's best window size 
    int startIdx = currentWindowIdx>0 ? currentWindowIdx-2: 0; 
    int endIdx   = currentWindowIdx>0 ? currentWindowIdx+2: windowWidth.size()-1; 

    //==// always scan config.windowNumAlways windows, even if you get to the edges of the windowWidth vector
    startIdx = 0;
    endIdx   = windowWidth.size()-1;
    /*
    if (startIdx < 0 || currentWindowIdx == 0) 
    {                       
      startIdx = 0;
      endIdx   = config.windowNumAlways - 1;
    }
  
    if (endIdx > windowWidth.size()-1)
    { 
      startIdx = (windowWidth.size()-1) - (config.windowNumAlways-1);  
      endIdx   = windowWidth.size()-1; 
    }
    */
     
    /////////////////////////////////////////////////////// 
    //
    //== Loop over the window sizes 
    //
    ///////////////////////////////////////////////////////
    for (int testWinIdx = startIdx; testWinIdx < endIdx; testWinIdx++)
    { 
    
      //==// set window size   
      setWindowSize(lowEdgeTmp, highEdgeTmp, windowCenter, v_binLowEdge, windowWidth[testWinIdx],winIdx);  
      
      //==// skip if there is more than 1 global window size 
      if ( lowEdgeTmp == config.xminfit && highEdgeTmp == config.xmaxfit) 
        numGlobalFits++; 
      if (numGlobalFits > 1)
      { 
        cout << "Global window reached! Skip " << windowWidth[testWinIdx] << "x" << endl; 
        continue; 
      }
        
      cout << endl;
      cout << "--------------------------------------------------------------------------------------------------------------" << endl;
      cout << "------------> TESTING "<< windowWidth[testWinIdx] << "x window size: [ " << lowEdgeTmp << ", " << highEdgeTmp << " ]" << endl;
      cout << "--------------------------------------------------------------------------------------------------------------" << endl;
      cout << endl;

      cout << "::::::::::::::::::::::::::::: NOMINAL BACKGROUND ONLY FIT :::::::::::::::::::::::::::::::: " << endl;

      //==// prepare nominal function 
      TF1 *BFuncTmp = new TF1("BFuncTmp", getFunction(config.bkgFunc), lowEdgeTmp, highEdgeTmp);
      formatFunc(BFuncTmp, kBlue, 2);
      armFunc(BFuncTmp, paramNomBonly);  
   
      //==// perform fit  
      vector<double> btemp_fitQuality;
      double LLH_BTmp_Nom = -minuitFit ( histRaw, BFuncTmp, 0, lowEdgeTmp, highEdgeTmp, config.quietFitInfo, &btemp_fitQuality);    
    
      //==// convert bkg fit to histogram for calculating chisquare p-value
      TH1D *BHistTmp   = (TH1D*) histRaw -> Clone("BHistTmp");
      fit_to_hist(BHistTmp, BFuncTmp, lowEdgeTmp, highEdgeTmp); 

      //======================================================                                                                                                    
      //  
      // Pick window with best bkg only chisq p-value
      //  
      //======================================================
      double chisqpValTmpBonly   = ChiSqNDFTest(histRaw, BHistTmp, lowEdgeTmp, highEdgeTmp, BFuncTmp->GetNDF(), 0, 0, 1);
      double chisqTmpBonly       = ChiSqNDFTest(histRaw, BHistTmp, lowEdgeTmp, highEdgeTmp, BFuncTmp->GetNDF(), 1, 0, 0);
      double NDF                 = BFuncTmp->GetNDF(); 
      cout << "                                                           ChiSquare p-value: " << chisqpValTmpBonly << endl;  
      if ( chisqpValTmpBonly >= bestpVal )
      { 
        BHistBkg         = (TH1D*) BHistTmp->Clone("BHistBkg");
        bestpVal         = chisqpValTmpBonly;
        windowNDF        = NDF; 
        bestChisq        = chisqTmpBonly; 
        bbest_fitQuality = btemp_fitQuality; 
        bestWindow       = windowWidth[testWinIdx];
        bestWindowIdx    = testWinIdx; 
        lowEdge          = lowEdgeTmp; 
        highEdge         = highEdgeTmp;

      }
 
      //==// delete variables created in the loop 
      delete BFuncTmp; 
      delete BHistTmp; 
    }

    cout << endl; 
    cout << "==============================================================================================================" << endl; 
    cout << "------------> PICKING "<< bestWindow << "x window size: [ " << lowEdge << ", " << highEdge << " ] for window " << windowCenter << endl;
    cout << "==============================================================================================================" << endl;
    cout << endl;   

    //==// check nominal bkg-only fit for fails 
    checkForFails( bbest_fitQuality, &vv_fitQuality_BnomOnly, winIdx, start_bin, windowCenter );
    
    //==// save best window sizes to be used later 
    v_bestWindow.push_back(bestWindow); 
    v_windowLow.push_back(lowEdge); 
    v_windowHigh.push_back(highEdge); 
    v_ChisquarepVal_Bnom.push_back(bestpVal);   

    //==// pass best window information to initialize window scan for next windows 
    currentWindowSize = bestWindow; 
    currentWindowIdx  = bestWindowIdx; 

    //==// create SWiFt bkg using nominal bkg function
    makeBkg(winIdx, bkgNomHist, BHistBkg ); 

    //==// plot fits in window 
    if (!config.doSWiFtScan)
      plotWindowFits(windowCenter, lowEdge, highEdge, histRaw, BHistBkg, canvas, pdf, bestpVal, bestChisq, windowNDF  );

    //==// delete variables created in the loop 
    delete BHistBkg; 
  }  


  //==// print out fit quality info for fits that fail 
  cout << endl; 
  cout << "Failed Nominal Background-only Fits: " << endl;                                                                                
  //cout << "Mass, Simplex error Flag, Minuit error Flag, Hesse error Flag, LLH" << endl; 
  reportFail(vv_fitQuality_BnomOnly); 
  cout << endl; 

 
  //===============================================================================================                                                                  
  //===============================================================================================
  //  
  // PART B: Perform SWiFt's LLH ratio scan for resonance search  
  //   -> Window sizes decided by part A are used here
  //   -> Nominal and alternate bkg-only, nominal and alternate S+B fits are performed 
  //   -> A choice of using the nonimal or alternate bkg function is made by picking
  //      the S+B fit with the better chiSquare p-value. This accounts for the bkg 
  //      choice uncertainty 
  //   -> the LLH ratio of S+B and b-only is converted to a local p-value using wilk's theorem 
  //   -> If requested, the 95% upper limit is calculated using a binary search method 
  //   -> If requested, signal systematics are included in the LLH ratio scan and 95% upper limit 
  //      Signal systematics are profiled by adding nuisance parameters to the likelihood 
  //  
  //===============================================================================================
  //================================================================================================

  if (config.doSWiFtScan)
  { 
    
    cout << endl;
    cout << "-------------------------------------------------------------------" << endl;
    cout << "------------> Perform the resonance scan "<< endl;
    cout << "-------------------------------------------------------------------" << endl;
    cout << endl;

    printStuff("../ascii/owl.txt"); 
    cout << endl; 

    //==// make pads for plots per window                                                                                    
    cout << "... make TPads " << endl; 
    TPad *pad1 = new TPad("pad1","pad1",0,0.51,1,1);
    pad1 -> SetBottomMargin(0);
    pad1 -> SetLeftMargin(0.12);  
    pad1 -> SetLogy(1);
    pad1 -> Draw(); 

    TPad *pad2 = new TPad("pad2","pad2",0,0.31,1,0.51);
    pad2 -> SetBottomMargin(.4);
    pad2 -> SetLeftMargin(0.12);
    pad2 -> SetGridx(1);
    pad2 -> SetGridy(1);
    pad2 -> Draw();
  
    TPad *pad3 = new TPad("pad3","pad3",0,0.18,1,0.38);
    pad3 -> SetBottomMargin(.4);
    pad3 -> SetLeftMargin(0.12);
    pad3 -> SetGridx(1);
    pad3 -> SetGridy(1);
    pad3 -> Draw();
  
    TPad *pad4 = new TPad("pad4","pad4",0,0.05,1,0.25);                                                                              
    pad4 -> SetBottomMargin(.4);                                                                                                    
    pad4 -> SetLeftMargin(0.12);
    pad4 -> SetGridx(1); 
    pad4 -> SetGridy(1);
    pad4 -> Draw();

    /////////////////////////////////////////////////////// 
    //
    //== loop over bin edges 
    //
    /////////////////////////////////////////////////////// 
    for (int winIdx = start_bin; winIdx <= end_bin; winIdx++)
    {

      //==// define variables, etc needed in the loop
      TF1 *sigFunc; 
      TH1D *sigHist;  
      TF1 *BFunc, *BFuncNom, *BFuncAlt;
      TH1D *BHist, *BHistNom, *BHistAlt;
      TF1 *SBFunc, *SBFuncNom, *SBFuncAlt;
      TH1D *SBHist, *SBHistNom, *SBHistAlt;
      int bMin, bMax; 
      double LLH_B, LLH_B_Nom, LLH_B_Alt; 
      double LLH_SB, LLH_SB_Nom, LLH_SB_Alt;

      //==// set window center and ranges 
      windowCenter  = v_binLowEdge[winIdx]; 
      lowEdge       = v_windowLow[winIdx - start_bin]; 
      highEdge      = v_windowHigh[winIdx - start_bin]; 
      bMin          = histRaw -> GetXaxis()->FindBin(lowEdge);
      bMax          = histRaw -> GetXaxis()->FindBin(highEdge); 

      cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Window Center: " << windowCenter << endl;
      cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Low Edge  : "  << lowEdge   << " (bin " << bMin << ")" << endl; 
      cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& High Edge : "  << highEdge  << " (bin " << bMax << ")" << endl;
      cout << endl;  

      //==// set the range of the raw data histogram                                                                                     
      histRaw         -> GetXaxis() -> SetRangeUser(lowEdge, highEdge);
      histDivBinWidth -> GetXaxis() -> SetRangeUser(lowEdge, highEdge);

      //======================================================                                                                                                    
      //  
      // Perform nominal bkg-only fit   
      //  
      //======================================================
      cout << endl; 
      cout << "::::::::::::::::::::::::::::: NOMINAL BACKGROUND ONLY FIT :::::::::::::::::::::::::::::::: " << endl;

      //==// prepare nominal function 
      BFuncNom = new TF1("BFuncNom", getFunction(config.bkgFunc), lowEdge, highEdge);
      formatFunc(BFuncNom, kBlue, 2);                                                                                                                           
      armFunc(BFuncNom, paramNomBonly);  
   
      //==// perform fit  
      vector<double> bNom_fitQuality;
      LLH_B_Nom = -minuitFit ( histRaw, BFuncNom, 0, lowEdge, highEdge, config.quietFitInfo, &bNom_fitQuality);    
    
      //==// store bkg parameters for plotting later 
      for (int i = 0; i < config.par.size(); i++) 
        vv_BkgParam_Bnom[i].push_back( BFuncNom->GetParameter(i) );

      //==// check nominal bkg-only fit for fails 
      checkForFails( bNom_fitQuality, &vv_fitQuality_Bnom, winIdx, start_bin, windowCenter );

      //==// convert bkg fit to histogram for each window 
      BHistNom   = (TH1D*) histRaw -> Clone("BHistNom");
      fit_to_hist(BHistNom, BFuncNom, lowEdge, highEdge);

      //==// save parameters to initialize next window 
      for (int i = 0; i < config.par.size() ; i++) 
        paramNomBonly[i] = BFuncNom->GetParameter(i);


      //======================================================                                                                                                    
      //  
      // Perform alternate bkg-only fit   
      //  
      //======================================================
      cout << endl; 
      cout << "::::::::::::::::::::::::::::: ALTERNATE BACKGROUND ONLY FIT :::::::::::::::::::::::::::::::: " << endl;

      //==// prepare alternate function 
      BFuncAlt = new TF1("BFuncAlt", getFunction(config.bkgFunc2), lowEdge, highEdge);
      formatFunc(BFuncAlt, kCyan, 2);                                                                                                                           
      armFunc(BFuncAlt, paramAltBonly);  
   
      //==// perform fit  
      vector<double> bAlt_fitQuality;
      LLH_B_Alt = -minuitFit ( histRaw, BFuncAlt, 0, lowEdge, highEdge, config.quietFitInfo, &bAlt_fitQuality);    
    
      //==// store bkg parameters for plotting later 
      for (int i = 0; i < config.par2.size(); i++) 
        vv_BkgParam_Balt[i].push_back( BFuncAlt->GetParameter(i) );
      
      //==// check nominal bkg-only fit for fails 
      checkForFails( bAlt_fitQuality, &vv_fitQuality_Balt, winIdx, start_bin, windowCenter );

      //==// convert bkg fit to histogram for each window 
      BHistAlt   = (TH1D*) histRaw -> Clone("BHistAlt");
      fit_to_hist(BHistAlt, BFuncAlt, lowEdge, highEdge);

      //==// save parameters to initialize next window 
      for (int i = 0; i < config.par2.size() ; i++) 
        paramAltBonly[i] = BFuncAlt->GetParameter(i);


      //======================================================                                                                                                    
      //  
      // Build signal shape at window center
      //  
      //======================================================
      //==// select signal shape to be define at window center
      sigFunc                 = new TF1("sigFunc", getFunction(config.SonlyFunc), config.xminfit, config.xmaxfit);
      vector<double> SigParam = generateSignal(windowCenter, nomSigFile, sigFunc);
    
      //==// convert TF1 to histogram to obtain RMS
      sigHist  = (TH1D*) histRaw -> Clone("sigHist");
      fit_to_hist(sigHist, sigFunc, config.xminfit, config.xmaxfit);
      sigHist  -> GetXaxis() -> SetRangeUser(config.xminfit, config.xmaxfit);
      sigHist  -> SetLineColor(kBlack); 
      sigHist  -> SetLineWidth(2); 
      sigHist  -> Scale(1/sigHist->Integral(xminfit_bin, xmaxfit_bin));
      int numSigParam = SigParam.size();

      //==// when including the JES signal systematics, calculate the 1sigma variations to generate the systematically varied shape  
      //==// the results are stored in vvv_diffVec
      vector<vector<vector<double>>> vvv_diffVec(2, vector<vector<double>>(JES1upSigFile.size(), vector<double>(numSigParam)));                
      // The indexes of the vector above are as follow: 
      // 1st: number of sigmas, in this case we always have 2, i.e. 1sigma up and 1sigma down
      // 2nd: number of JES variations, JES1, JES2, JES3
      // 3rd: number of parameters of the signal shape, 1 normalization + 5 shape parameters  
      if (config.addSys)
        calculateDiff(windowCenter, sigFunc, JES1upSigFile, JES1downSigFile, vvv_diffVec);


      //======================================================                                                   
      //  
      // S+B fit with Nominal bkg function 
      //  
      //======================================================
      cout << endl; 
      cout << "::::::::::::::::::::::::::::: S+B FIT (Nominal) :::::::::::::::::::::::::::::::: " << endl;
 
      //==// pick s+b fit
      TString SBFuncType = getSBFunction(config.bkgFunc, config.SonlyFunc);

      //==// prepare s+b fit 
      SBFuncNom = new TF1("SBFuncNom", getFunction(SBFuncType), lowEdge, highEdge);
      formatFunc(SBFuncNom, kRed, 2);
      armFunc(SBFuncNom, sigFunc, BFuncNom); 

      //==// perform fit
      vector<double> sbNom_fitQuality;
      vector<double> NPNom_values_window;
      vector<double> NPNom_errors_window;
      LLH_SB_Nom = -minuitFit ( histRaw, SBFuncNom, 1, lowEdge, highEdge, config.quietFitInfo, &sbNom_fitQuality,
                                &vvv_diffVec, &NPNom_values_window, &NPNom_errors_window );

      //==// store bkg parameters for plotting later 
      for (int i = numSigParam; i < numSigParam+config.par.size(); i++)  
        vv_BkgParam_SBnom[i-numSigParam].push_back( SBFuncNom->GetParameter(i) ); 
      
      //==// convert s+b bkg component to histogram for each window 
      SBHistNom = (TH1D*) histRaw -> Clone("SBHistNom");
      fit_to_hist(SBHistNom, SBFuncNom, lowEdge, highEdge); 

      //==// save details if fit fails 
      checkForFails( sbNom_fitQuality, &vv_fitQuality_SBnom, winIdx, start_bin, windowCenter );


      //======================================================                                                   
      //  
      // S+B fit with Alternate bkg function 
      //  
      //======================================================
      cout << endl; 
      cout << "::::::::::::::::::::::::::::: S+B FIT (Alternate) :::::::::::::::::::::::::::::::: " << endl;
 
      //==// pick s+b fit
      SBFuncType = getSBFunction(config.bkgFunc2, config.SonlyFunc);

      //==// prepare s+b fit 
      SBFuncAlt = new TF1("SBFuncAlt", getFunction(SBFuncType), lowEdge, highEdge);
      formatFunc(SBFuncAlt, kRed, 2);
      armFunc(SBFuncAlt, sigFunc, BFuncAlt); 

      //==// perform fit
      vector<double> sbAlt_fitQuality;
      vector<double> NPAlt_values_window;
      vector<double> NPAlt_errors_window;
      LLH_SB_Alt = -minuitFit ( histRaw, SBFuncAlt, 1, lowEdge, highEdge, config.quietFitInfo, &sbAlt_fitQuality, 
                                &vvv_diffVec, &NPAlt_values_window, &NPAlt_errors_window );    

      //==// store bkg parameters for plotting later 
      for (int i = numSigParam; i < numSigParam+config.par2.size(); i++) 
        vv_BkgParam_SBalt[i-numSigParam].push_back( SBFuncAlt->GetParameter(i) );
      
      //==// convert s+b bkg component to histogram for each window 
      SBHistAlt = (TH1D*) histRaw -> Clone("SBHistAlt");
      fit_to_hist(SBHistAlt, SBFuncAlt, lowEdge, highEdge); 

      //==// save details if fit fails 
      checkForFails( sbAlt_fitQuality, &vv_fitQuality_SBalt, winIdx, start_bin, windowCenter );


      //===================================================================================                                                          
      //  
      // Pick Bonly and S+B function, based on chisquare p-value test, btw nominal and alternate S+B fits  
      //  
      //===================================================================================
      TString bkgFuncType;
      vector<double> NP_values_window;
      vector<double> NP_errors_window;
      double ChiSqPValNom = ChiSqNDFTest(histRaw, SBHistNom, lowEdge, highEdge, SBFuncNom->GetNDF(), 0, 0, 1); 
      double ChiSqPValAlt = ChiSqNDFTest(histRaw, SBHistAlt, lowEdge, highEdge, SBFuncAlt->GetNDF(), 0, 0, 1); 
      cout << endl; 
      cout << "==============================================================================================================" << endl; 
      cout << "------------> Choosing btw nominal or alternate S+B fit"  << endl; 
      cout << "-----> ChiSq p-value of Nominal S+B   : " << ChiSqPValNom << endl; 
      cout << "-----> ChiSq p-value of Alternate S+B : " << ChiSqPValAlt << endl; 
      cout << endl;
      //if ( ChiSqPValNom > ChiSqPValAlt )
      //{
        cout << "-----> NOMINAL S+B and B-only picked" << endl;  
        bkgFuncType      = config.bkgFunc; 
        LLH_B            = LLH_B_Nom; 
        BFunc            = (TF1*)BFuncNom->Clone("BFunc"); 
        BHist            = (TH1D*)BHistNom->Clone("BHist"); 
        LLH_SB           = LLH_SB_Nom;  
        SBFunc           = (TF1*)SBFuncNom->Clone("SBFunc"); 
        SBHist           = (TH1D*)SBHistNom->Clone("SBHist");
        if (config.addSys)
        { 
          for (int n = 0; n < vv_NP_values.size(); n++)
          {                                                                                                                                                         
            vv_NP_values[n].push_back(NPNom_values_window[n]); 
            vv_NP_errors[n].push_back(NPNom_errors_window[n]); 
          }  
        }
        v_numParamChoosen.push_back(BFuncNom->GetNpar());
      //}  
      /*else 
      { 
        cout << "-----> ALTERNATE S+B and B-only picked" << endl;
        bkgFuncType      = config.bkgFunc2;
        LLH_B            = LLH_B_Alt; 
        BFunc            = (TF1*)BFuncAlt->Clone("BFunc"); 
        BHist            = (TH1D*)BHistAlt->Clone("BHist"); 
        LLH_SB           = LLH_SB_Alt;  
        SBFunc           = (TF1*)SBFuncAlt->Clone("SBFunc"); 
        SBHist           = (TH1D*)SBHistAlt->Clone("SBHist");
        if (config.addSys)
        { 
          for (int n = 0; n < vv_NP_values.size(); n++)
          {                                                                                                                                                         
            vv_NP_values[n].push_back(NPAlt_values_window[n]); 
            vv_NP_errors[n].push_back(NPAlt_errors_window[n]); 
          }  
        }
        v_numParamChoosen.push_back(BFuncAlt->GetNpar());
      }
      */
      cout << "==============================================================================================================" << endl;
      cout << endl; 
  
      //======================================================                                                                                                    
      //  
      // Make and save SWiFt bkg  
      //  
      //======================================================
      makeBkg(winIdx, bkgBestHist, BHist); 

      //======================================================                                                                                                    
      //  
      // Get components of S+B fit 
      //  
      //======================================================
      //==// get background component
      TF1 *SBFuncB = new TF1("SBFuncB", getFunction(bkgFuncType), lowEdge, highEdge);
      formatFunc(SBFuncB, kRed, 2, 2);
      armFunc(SBFuncB, SBFunc, sigFunc->GetNpar());

      //==// convert bkg component to histogram
      TH1D *SBHistB = (TH1D*) histRaw -> Clone("SBHistB");
      fit_to_hist(SBHistB, SBFuncB, lowEdge, highEdge); 

      //==// get signal componenet
      TF1 *SBFuncS = new TF1("SBFuncS", getFunction(config.SonlyFunc), lowEdge, highEdge);
      formatFunc(SBFuncS, kRed, 1);
      armFunc(SBFuncS, SBFunc, 0);

      //==// convert signal component histogram   
      TH1D *SBHistS  = (TH1D*) histRaw -> Clone("SBHistS");
      fit_to_hist(SBHistS, SBFuncS, lowEdge, highEdge);

      //======================================================                                                                                                    
      //  
      // Calculate variables from fits 
      //  
      //======================================================
      //==// extracted signal, error, bkg estimates 
      double signal    = SBFunc -> Integral(lowEdge, highEdge) - SBFuncB->Integral(lowEdge, highEdge); 
      double signalErr = abs(SBFunc -> GetParError(0) * (signal / SBFunc -> GetParameter(0)));
      v_ExtractSig.push_back( signal );                        
      v_ExtractSigError.push_back ( signalErr );
      v_ExtractSigRaw.push_back( SBFunc->GetParameter(0) );  
      v_WindowBkg.push_back ( histRaw -> Integral(bMin, bMax) );                                        
      v_ExtractBkg.push_back( SBFuncB-> Integral(lowEdge, highEdge) );       
              
      //==// significances, likelihoods and ratios 
      double sig       = SBFunc -> GetParameter(0) / SBFunc->GetParError(0);
      double LLHR      = LLH_SB - LLH_B;
      double LLHR_pval = LLHRPval(LLHR, 1)/2;
      v_Significance.push_back( sig ); 
      v_LogLikeHood.push_back( LLHR );
      if (signal > 0)
        v_LogLikeHoodPvalue.push_back(LLHR_pval);
      else 
        v_LogLikeHoodPvalue.push_back(1-LLHR_pval);
      v_LogLikeHood_SB.push_back( LLH_SB );
      v_LogLikeHood_B.push_back( LLH_B );    
  
      //==// chisquares, goodness of fit measures 
      double chisq          = ChiSqNDFTest(histRaw, SBHist, lowEdge, highEdge, SBFunc->GetNDF(), 1, 0, 0 ); 
      double chisqBonly     = ChiSqNDFTest(histRaw, BHist, lowEdge, highEdge, BFunc->GetNDF(), 1, 0, 0 );
      double chisqNDF       = ChiSqNDFTest(histRaw, SBHist, lowEdge, highEdge, SBFunc->GetNDF(), 0, 1, 0 ); 
      double chisqNDFBonly  = ChiSqNDFTest(histRaw, BHist, lowEdge, highEdge, BFunc->GetNDF(), 0, 1, 0 );
      double chisqpVal      = ChiSqNDFTest(histRaw, SBHist, lowEdge, highEdge, SBFunc->GetNDF(), 0, 0, 1 );
      double chisqpValBonly = ChiSqNDFTest(histRaw, BHist, lowEdge, highEdge, BFunc->GetNDF(), 0, 0, 1 );
      v_Chisquare_SB.push_back( chisq );  
      v_Chisquare_B.push_back( chisqBonly );
      v_ChisquareNdf_SB.push_back( chisqNDF ); 
      v_ChisquareNdf_B.push_back( chisqNDFBonly );
      v_ChisquarepVal_SB.push_back( chisqpVal );
      v_ChisquarepVal_B.push_back( chisqpValBonly );

      //==// signal window correction factor 
      double sigHistWin  = sigHist -> Integral( sigHist->FindBin(lowEdge)       , sigHist->FindBin(highEdge)       ) ;  
      double sigHistFull = sigHist -> Integral( sigHist->FindBin(config.xminfit), sigHist->FindBin(config.xmaxfit) ) ;  
      double sigLeakage  = sigHistWin / sigHistFull ;  
      //double sigLeakage = sigFunc->Integral(lowEdge, highEdge) / sigFunc->Integral(config.xminfit, config.xmaxfit); 
      v_WindowCorr.push_back ( sigLeakage );  

      //======================================================                                                   
      //  
      // Calculate 95% upper limit  
      //  
      //======================================================
      TF1* SBFunc95  = (TF1*)SBFunc->Clone();
      TF1* SBFuncS95 = (TF1*)sigFunc->Clone();
      TH1D *SBHistS95 = (TH1D*) histRaw -> Clone("SBHistS95"); 
      double Signal95 = 0; 
      if (config.do95Scan)
      {
        cout << "::::::::::::::::::::::::::::: SCAN FOR 95% UPPER LIMIT :::::::::::::::::::::::::::::::: " << endl;  
        vector<double> NP95_values_window;
        vector<double> NP95_errors_window;
        Signal95 = LLHScan(LLH_SB, SBFunc, SBFunc95, BFunc, sigFunc, signal, signalErr, histRaw, lowEdge, highEdge, windowCenter,
                           &vvv_diffVec, &NP95_values_window, &NP95_errors_window );
       
        //==// get pulls 
        if (config.addSys)
        {    
          for (int n = 0; n < vv_NP95_values.size(); n++) 
          {                                                                                                                                                              
            vv_NP95_values[n].push_back(NP95_values_window[n]); 
            vv_NP95_errors[n].push_back(NP95_errors_window[n]); 
          }    
        } 
        
        //==// save 95% upper limit on signal and cross-section 
        v_Signal95CL.push_back(Signal95);
        v_Limit95CL.push_back(Signal95/ (config.lumi*1000*sigLeakage) ); 
  
        // plot on histograms for each window 
        formatFunc(SBFunc95, kGreen+1, 2); 
    
        //==// get signal shape after 95% UL scan 
        for (int s = 0; s < sigFunc->GetNpar() ; s++)
          SBFuncS95 -> FixParameter(s, SBFunc95->GetParameter(s)); 

        //==// convert s component of s+b to histogram 
        fit_to_hist(SBHistS95, SBFuncS95, lowEdge, highEdge); 
        SBHistS95 -> SetLineColor(kGreen+1);
        SBHistS95 -> SetLineWidth(2);
    
      }
      else
      { 
        v_Signal95CL.push_back(0.0); 
        v_Limit95CL.push_back(0.0); 
      }  

      //======================================================                                                                                                    
      //  
      // Draw histogram, fits, details for each window 
      //  
      //======================================================
     
      //===================== enter pad1
      pad1 -> cd(); 
      pad1 -> SetLogy(1);   
 
      //==// print all plots to pad1 of the conavs  
      histDivBinWidth -> Draw("P E1");
      BFunc           -> Draw("same"); 
      if (config.do95Scan) SBFunc95 -> Draw("same");
      SBFunc          -> Draw("same"); 
      SBFuncB         -> Draw("same");  
      histDivBinWidth -> SetMarkerSize(0.8);
      histDivBinWidth -> SetMinimum( BFunc->Eval(highEdge) );  

      //==// draw line to mark center  
      double ymax = SBFuncB -> Eval(windowCenter); 
      makeTLine(windowCenter, 0, windowCenter, ymax, kGreen-10); 
 
      //==// print fit details to canvas 
      drawText(0.68, 0.85,"#bf{#it{ATLAS}} Internal", kBlack, 0.07);
      drawText(0.68, 0.78, "Mass: "+roundUp(windowCenter,1)+" GeV", kBlack, 0.05 );
      drawText(0.68, 0.72, "Window: ["+roundUp(lowEdge,0)+", "+roundUp(highEdge,0)+"] GeV", kBlack, 0.05 );
      drawText(0.68, 0.66, "Signal: "+roundUp(signal)+" #pm "+roundUp(signalErr), kBlack, 0.05 );
      drawText(0.68, 0.60, "LLH Ratio p-value: "+roundUp(LLHR_pval, 5, 0), kBlack, 0.05 ); 
      drawText(0.68, 0.54, "#chi^{2}/NDF (SB, B): "+roundUp(chisq)+"/"+to_string(SBFunc->GetNDF())+", "
               +roundUp(chisqBonly)+"/"+to_string(BFunc->GetNDF()), kBlack, 0.05 );
    
      //===================== enter pad2
      pad2 -> cd(); 
      TH1D *diffBkgOnly = (TH1D*) histRaw -> Clone("diffBkgOnly");
      diffBkgOnly -> Add(BHist , -1);
      diffBkgOnly -> Draw("hist E1 P");
      diffBkgOnly -> SetMarkerSize(0.8);
      formatDiff(diffBkgOnly, "Data-Bkg_{only}", "m_{jj} [GeV]", 0, 1, 1);
      makeTLine(lowEdge, 0, highEdge, 0, kBlack, "same"); 

      //===================== enter pad3
      pad3->cd();
      TH1D *diffBfromSB = (TH1D*) histRaw -> Clone("diffBfromSB");
      diffBfromSB -> Add(SBHistB , -1);
      diffBfromSB -> Draw("hist E1 P");
      diffBfromSB -> SetMarkerSize(0.8);
      formatDiff(diffBfromSB, "Data-Bkg_{SB}", "m_{jj} [GeV]", 0, 1, 1);
      SBHistS -> SetLineWidth(2); 
      SBHistS -> Draw("hist same"); 
      makeTLine(lowEdge, 0, highEdge, 0, kBlack, "same");
    
      //===================== enter pad4
      pad4->cd();
      TH1D *diffBfromSBUn = (TH1D*) histRaw -> Clone("diffBfromSBUn");
      diffBfromSBUn -> Add(SBHistB , -1);
      divideUncert(diffBfromSBUn, histRaw); 
      diffBfromSBUn -> Draw("hist P");
      diffBfromSBUn -> SetMarkerSize(0.8);
      formatDiff(diffBfromSBUn, "#frac{Data-Bkg_{SB}}{#sqrt{Data}} ", "m_{jj} [GeV]", 0, 1, 0);
      makeTLine(lowEdge, 0, highEdge, 0, kBlack, "same");
    
      //==// print all pads to canvas 
      canvas     -> Print(pdf);

      //==// if signal systematics are asked for, print out the nominal, best and 95% signal shapes 
      //==// All these shapes should change abit. 
      if (config.addSys)
      {  
        pad1 -> cd();  
        pad1 -> SetLogy(0);   
        SBHistS -> Scale(1/SBHistS->Integral(bMin, bMax));
        SBHistS -> Draw("hist");
        SBHistS -> GetYaxis() -> SetTitle("Unit Area"); 
        sigHist -> Scale(1/sigHist->Integral(bMin, bMax)); 
        sigHist -> Draw("hist same");
        SBHistS -> Draw("hist same"); 
        if (config.do95Scan) 
        {    
          SBHistS95 -> Scale(1/SBHistS95->Integral(bMin, bMax)); 
          SBHistS95 -> Draw("hist same"); 
        }    
        drawText(0.72, 0.85,"#bf{#it{ATLAS}} Internal", kBlack, 0.07);
        drawText(0.72, 0.78, "Nominal Signal", kBlack, 0.05 );                                                                      
        drawText(0.72, 0.72, "Best Fit Signal", kRed, 0.05 );
        drawText(0.72, 0.66, "95% C.L. Signal", kGreen+1, 0.05 );
        canvas  -> Print(pdf); 
      } 

      //==// delete objects 
      delete sigFunc, sigHist; 
      delete BFunc, BFuncNom, BFuncAlt; 
      delete BHist, BHistNom, BHistAlt;
      delete SBFunc, SBFuncNom, SBFuncAlt;  
      delete SBHist, SBHistNom, SBHistAlt;  
      delete SBFuncB; 
      delete SBHistB;
      delete SBFuncS; 
      delete SBHistS;
      delete SBFunc95; 
      delete SBHistS95; 
      delete SBFuncS95; 
      delete diffBkgOnly; 
      delete diffBfromSB; 
      delete diffBfromSBUn; 
    
    }   

    //==// print out fit quality info for fits that fail 
    cout << endl; 
    cout << "Failed Nominal Background-only Fits: " << endl;
    //cout << "Mass, Simplex error Flag, Minuit error Flag, Hesse error Flag, LLH" << endl; 
    reportFail(vv_fitQuality_Bnom); 
    cout << endl; 

    cout << endl; 
    cout << "Failed Alternate Background-only Fits: " << endl;
    reportFail(vv_fitQuality_Balt); 
    cout << endl; 
 
    cout << endl; 
    cout << "Failed Nominal S+B Fits: " << endl;
    reportFail(vv_fitQuality_SBnom); 
    cout << endl; 

    cout << endl; 
    cout << "Failed Alternate S+B Fits: " << endl;
    reportFail(vv_fitQuality_SBalt); 
    cout << endl; 

  } 


  //===============================================================================================                                                                  
  //  
  // PART C: Signal subtraction stage  
  //   -> If these is a significant fluctuation in the data, perform signal subtraction   
  //   -> The user gets to define what significant means here - based on local p-value
  //   -> The SWiFt bkg is re-made after signal subtraction is performed 
  //   -> This reduces the bais in the SWiFt bkg 
  //  
  //================================================================================================
 
  TH1D* histRawSS; 
  if (config.doSWiFtScan)
  { 

    histRaw -> GetXaxis() -> SetRangeUser(config.xminfit, config.xmaxfit);

    //==// Find maximum likelihood   
    double maxLLH     = 10.0; 
    double bestMass   = 0.0; 
    double bestSig    = 0.0; 
    double bestRawSig = 0.0; 
    for (int i = 0; i < v_MassNum.size(); i++)
      if (v_LogLikeHoodPvalue[i] < maxLLH) 
      {
        maxLLH     = v_LogLikeHoodPvalue[i];  
        bestMass   = v_MassNum[i];
        bestSig    = v_ExtractSig[i]; 
        bestRawSig = v_ExtractSigRaw[i]; 
      } 


    //==// If LLH < config.sigSubLLHCut, trigger signal subtraction and re-create the SWiFt bkg 
    if (maxLLH < config.sigSubLLHCut)
    { 
      canvas -> Clear(); 
  
      //==// select signal shape to be define at window center
      TF1* sigFuncSub            = new TF1("sigFuncSub", getFunction(config.SonlyFunc), config.xminfit, config.xmaxfit);
      vector<double> SigParamSub = generateSignal(bestMass, nomSigFile, sigFuncSub);
      sigFuncSub                 -> SetParameter(0, bestRawSig);   

      //==// convert TF1 to histogram
      TH1D* sigHistSub  = (TH1D*) histRaw -> Clone("sigHistSub");
      fit_to_hist(sigHistSub, sigFuncSub, config.xminfit, config.xmaxfit);
      sigHistSub  -> GetXaxis() -> SetRangeUser(config.xminfit, config.xmaxfit);

      //==// subtract signal from data 
      histRawSS = (TH1D*)histRaw -> Clone();
      histRawSS -> Add(sigHistSub, -1);  

      //==// re-generate SWiFt bkg from signal subtracted data 
      printStuff("../ascii/owl.txt"); 
      cout << endl; 
 
      /////////////////////////////////////////////////////// 
      //
      //== loop over bin edges 
      //
      /////////////////////////////////////////////////////// 
       
      vector<double> paramSSBonly(10, 0.0); 
      for (int winIdx = start_bin; winIdx <= end_bin; winIdx++)
      {

        //==// define variables, etc needed in the loop
        TF1 *BFuncSS; 
        TH1D *BHistSS;
        double LLH_B_SS; 

        //==// set window center and ranges 
        double windowCenter  = v_binLowEdge[winIdx]; 
        double lowEdge       = v_windowLow[winIdx - start_bin]; 
        double highEdge      = v_windowHigh[winIdx - start_bin]; 
        int bMin             = histRawSS -> GetXaxis()->FindBin(lowEdge);
        int bMax             = histRawSS -> GetXaxis()->FindBin(highEdge); 
        int numParamChoosen  = v_numParamChoosen[winIdx - start_bin]; 

        cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Window Center: " << windowCenter << endl;
        cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Low Edge  : "  << lowEdge   << " (bin " << bMin << ")" << endl; 
        cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& High Edge : "  << highEdge  << " (bin " << bMax << ")" << endl;
        cout << endl;  
  
        //==// set the range of the raw data histogram
        histRawSS  -> GetXaxis() -> SetRangeUser(lowEdge, highEdge);
  
        //==// perform nominal bkg-only fit   
        cout << endl; 
        cout << "::::::::::::::::::::::::::::: NOMINAL BACKGROUND ONLY FIT :::::::::::::::::::::::::::::::: " << endl;

        //==// initialize parameters 
        if (numParamChoosen == config.par.size()) 
          for (int i = 0; i < config.par.size(); i++ )
            paramSSBonly[i] = vv_BkgParam_Bnom[i][winIdx-start_bin];   

        if (numParamChoosen == config.par2.size()) 
          for (int i = 0; i < config.par2.size(); i++ )
            paramSSBonly[i] = vv_BkgParam_Balt[i][winIdx-start_bin];   


        //==// prepare bkg function 
        TString pickBkgFunc = ""; 
        if ( numParamChoosen == config.par.size() ) 
          pickBkgFunc = config.bkgFunc; 
        if ( numParamChoosen == config.par2.size() ) 
          pickBkgFunc = config.bkgFunc2; 
        
        cout << "numParamChoosen: "<< numParamChoosen << endl; 
        BFuncSS = new TF1("BFuncSS", getFunction(pickBkgFunc), lowEdge, highEdge);
        formatFunc(BFuncSS, kBlue, 2);                       
        armFunc(BFuncSS, paramSSBonly);  
          

        //==// perform fit  
        vector<double> bNom_fitQuality;
        LLH_B_SS = -minuitFit ( histRawSS, BFuncSS, 0, lowEdge, highEdge, config.quietFitInfo, &bNom_fitQuality);    
  
        //==// convert bkg fit to histogram for each window 
        BHistSS   = (TH1D*) histRawSS -> Clone("BHistNomSS");
        fit_to_hist(BHistSS, BFuncSS, lowEdge, highEdge);

        //==// make swift bkg
        makeBkg(winIdx, bkgSSHist, BHistSS );

        //==// delete 
        delete BFuncSS; 
        delete BHistSS; 
      } 

      //==// plot signal subtracted bkg   
      histRawSS  -> GetXaxis() -> SetRangeUser(config.xminfit, config.xmaxfit);
      double chisqpVal = ChiSqNDFTest(histRawSS, bkgSSHist, config.xminfit, config.xmaxfit, BFuncGlobal->GetNDF()+config.par.size(), 0, 0, 1);
      plotBkg(canvas, pdf, histRawSS, bkgSSHist, "SS Data", "SS SWiFt Bkg", chisqpVal);

      histRaw  -> GetXaxis() -> SetRangeUser(config.xminfit, config.xmaxfit);
      chisqpVal = ChiSqNDFTest(histRaw, bkgSSHist, config.xminfit, config.xmaxfit, BFuncGlobal->GetNDF()+config.par.size(), 0, 0, 1);
      plotBkg(canvas, pdf, histRaw, bkgSSHist, "Data", "SS SWiFt Bkg", chisqpVal);

      //==// delete stuff 
      delete sigFuncSub; 
      delete sigHistSub; 
    }
    else 
      bkgSSHist = (TH1D*)bkgBestHist->Clone();
  } 

  //===============================================================================================                                                                  
  //  
  // PART D: All the plotting fun!  
  //   -> Plots dedicated to the SWiFt bkg outputted
  //   -> Plots dedicated to the SWiFt scan outputted
  //   -> Plotting code can be nasty looking, apologies in advance T_T 
  //  
  //================================================================================================

  //==// Define names
  TGraph* TG_paramChoosen; 
  TGraph* TG_windowSizeX; 
  TGraph* grshade; 
  TGraph* TG_windowHigh; 
  TGraph* TG_windowCenter; 
  TGraph* TG_windowLow;  
  TGraph* TG_ChisquarepVal_Bnom; 
  TGraph* TG_Signal; 
  TGraph* TG_SignalError; 
  TGraph* TG_WindowCorr;
  TGraph* TG_Signal95; 
  TGraph* TG_Limit95; 
  TGraph* TG_Signifi; 
  TGraph* TG_Chi_SB; 
  TGraph* TG_Chi_B; 
  TGraph* TG_ChiNdf_SB; 
  TGraph* TG_ChiNdf_B; 
  TGraph* TG_ChiNdfpVal_SB; 
  TGraph* TG_ChiNdfpVal_B;
  TGraph* TG_LLH_SB;  
  TGraph* TG_LLH_B; 
  TGraph* TG_LLH; 
  TGraph* TG_LLHpval; 
  TGraph* TG_LLH_positive;  
  TGraph* TG_LLH_negative; 
  TGraph* TG_LLHpval_positive; 
  TGraph* TG_LLHpval_negative; 
  TGraph* TG_2sig; 
  TGraph* TG_1sig; 
  TGraphErrors* TG_NP; 
  vector<TGraph*> TG_bkgParam_Bnom(config.par.size(), NULL); 
  vector<TGraph*> TG_bkgParam_SBnom(config.par.size(), NULL); 
  vector<TGraph*> TG_bkgParam_Balt(config.par2.size(), NULL); 
  vector<TGraph*> TG_bkgParam_SBalt(config.par2.size(), NULL); 


  //======================================================   
  //                                                                                                 
  // Draw & save TGraphs related to bkgOnly mode   
  // 
  //======================================================
  if (!config.doSWiFtScan)
  { 
    canvas -> Clear(); 

    //==// Draw nominal SWiFt bkg  
    histRaw  -> GetXaxis() -> SetRangeUser(config.xminfit, config.xmaxfit);
    double chisqpVal = ChiSqNDFTest(histRaw, bkgNomHist, config.xminfit, config.xmaxfit, BFuncGlobal->GetNDF() + config.par.size(), 0, 0, 1);
    cout<<"CHI2 p-value of SWift bkg "<<chisqpVal<<endl;   //WEI
    value3=chisqpVal;
    plotBkg(canvas, pdf, histRaw, bkgNomHist, "Data", "SWiFt Bkg", chisqpVal);
    
    //==// window size
    TG_windowSizeX = makeTGraph(v_MassNum, v_bestWindow, "Window Size"); 
    TG_windowSizeX -> Draw("APC"); 
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    canvas -> Print(pdf); 

    //==// window size
    int n = v_MassNum.size(); 
    grshade = new TGraph(2*n);
    formatTGraph(grshade, "m_{jj} [GeV]", "Window [GeV]"); 
    for (int i=0;i<n;i++) 
    {
       grshade->SetPoint(i,v_MassNum[i],v_windowHigh[i]);
       grshade->SetPoint(n+i,v_MassNum[n-i-1],v_windowLow[n-i-1]);
    }   
    grshade -> SetFillColor(18);
    grshade -> Draw("AF"); 
    grshade -> SetMaximum(config.xmaxfit+100); 
    grshade -> SetMinimum(config.xminfit-100);

    TG_windowHigh = makeTGraph(v_MassNum, v_windowHigh, "Window [GeV]", kRed); 
    TG_windowHigh -> Draw("PL"); 
    TG_windowCenter = makeTGraph(v_MassNum, v_MassNum, "Window [GeV]"); 
    TG_windowCenter -> Draw("PL"); 
    TG_windowLow = makeTGraph(v_MassNum, v_windowLow, "Window Size", kRed); 
    TG_windowLow -> Draw("PL");
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    canvas -> Print(pdf); 

    //==// chisq p-value of best window size 
    canvas -> SetLogy(); 
    TG_ChisquarepVal_Bnom = makeTGraph(v_MassNum, v_ChisquarepVal_Bnom, "ChiSquare p-value"); 
    TG_ChisquarepVal_Bnom -> Draw("APC"); 
    TG_ChisquarepVal_Bnom -> SetMaximum(1.0);   
    canvas -> Print(pdf);
    canvas -> SetLogy(0); 
     
    //==// Close canvas 
    canvas->Print(pdf+"]");

    //==// save TGraphs to root file 
    rootFile -> cd();    
    bkgNomHist            -> Write("swiftBkgNominal");     
    BHistGlobal           -> Write("globalFitNominal");
    BHist2Global          -> Write("globalFitAlternate");
    TG_windowSizeX        -> Write("windowSizeX"); 
    TG_windowLow          -> Write("windowLow"); 
    TG_windowHigh         -> Write("windowHigh"); 
    TG_ChisquarepVal_Bnom -> Write("chiSqPvalBnom"); 
    rootFile -> Close(); 
    delete rootFile; 

  }

  //======================================================     
  //                                                                                               
  // Draw TGraphs related to SWiFt scan   
  // 
  //======================================================
  if (config.doSWiFtScan)
  { 
    canvas -> Clear(); 
 
    //==// Draw best SWiFt bkg  
    histRaw  -> GetXaxis() -> SetRangeUser(config.xminfit, config.xmaxfit);
    double chisqpVal = ChiSqNDFTest(histRaw, bkgBestHist, config.xminfit, config.xmaxfit, BFuncGlobal->GetNDF()+config.par.size(), 0, 0, 1);
    cout<<"CHI2 p-value of SWift bkg "<<chisqpVal<<endl; 
    value3=chisqpVal;
    plotBkg(canvas, pdf, histRaw, bkgBestHist, "Data", "SWiFt Bkg", chisqpVal);
    
    //==// parameters choosen 
    TG_paramChoosen = makeTGraph(v_MassNum, v_numParamChoosen, "Parameters Choosen"); 
    TG_paramChoosen -> Draw("APC"); 
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    canvas    -> Print(pdf); 

    //==// window size
    TG_windowSizeX = makeTGraph(v_MassNum, v_bestWindow, "Window Size"); 
    TG_windowSizeX -> Draw("APC"); 
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    canvas    -> Print(pdf); 

    //==// window size
    int n = v_MassNum.size(); 
    grshade = new TGraph(2*n);
    formatTGraph(grshade, "m_{jj} [GeV]", "Window [GeV]"); 
    for (int i=0;i<n;i++) 
    {
       grshade->SetPoint(i,v_MassNum[i],v_windowHigh[i]);
       grshade->SetPoint(n+i,v_MassNum[n-i-1],v_windowLow[n-i-1]);
    }   
    grshade -> SetFillColor(18);
    grshade -> Draw("AF"); 
    grshade -> SetMaximum(config.xmaxfit+100); 
    grshade -> SetMinimum(config.xminfit-100);

    TG_windowHigh = makeTGraph(v_MassNum, v_windowHigh, "Window [GeV]", kRed); 
    TG_windowHigh -> Draw("PL"); 
    TG_windowCenter = makeTGraph(v_MassNum, v_MassNum, "Window [GeV]"); 
    TG_windowCenter -> Draw("PL"); 
    TG_windowLow = makeTGraph(v_MassNum, v_windowLow, "Window Size", kRed); 
    TG_windowLow -> Draw("PL");
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    canvas    -> Print(pdf); 

    //==// Window correction 
    TG_WindowCorr = makeTGraph(v_MassNum, v_WindowCorr, "Window Correction");
    TG_WindowCorr -> Draw("APC");
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    canvas    -> Print(pdf);

    //==// Extracted signal 
    TG_Signal = makeTGraph(v_MassNum, v_ExtractSig, "Extracted Signal"); 
    TG_Signal -> Draw("APC"); 
    makeTLine(config.start, 0, config.end, 0, kBlack, "", 1, 2); 
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    canvas    -> Print(pdf);

    TG_SignalError = makeTGraph(v_MassNum, v_ExtractSigError, "Extracted Signal Error");
    TG_SignalError -> Draw("APC");
    makeTLine(config.start, 0, config.end, 0, kBlack, "", 1, 2);
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    canvas    -> Print(pdf);
 

    //==// 95% U.L. on extracted signal 
    canvas -> SetLogy();
    TG_Signal95 = makeTGraph(v_MassNum, v_Signal95CL, "Extracted 95% Signal");
    TG_Signal95 -> Draw("APC"); 
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    canvas    -> Print(pdf); 

    //==// limit on cross-section 
    TG_Limit95 = makeTGraph(v_MassNum, v_Limit95CL, "#sigma x A x BR");
    TG_Limit95 -> Draw("APC"); 
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    canvas    -> Print(pdf); 
    canvas -> SetLogy(0);

    //==// Significance 
    TG_Signifi = makeTGraph(v_MassNum, v_Significance, "Significance [#Signal / Error in #Signal]");
    TG_Signifi -> Draw("APC"); 
    makeTLine(config.start, 0, config.end, 0, kBlack, "", 1, 2); 
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    //makeTLine(config.start, -0.4, config.end, -0.4, kRed,"", 1, 2);
    //makeTLine(config.start, 0.4, config.end, 0.4, kRed,"", 1, 2);
    canvas -> Print(pdf); 

    //==// Chisquare
    TG_Chi_SB = makeTGraph(v_MassNum, v_Chisquare_SB, "#chi^{2}", kRed); 
    TG_Chi_SB -> Draw("APC");
    TG_Chi_B  = makeTGraph(v_MassNum, v_Chisquare_B, "#chi^{2}");
    TG_Chi_B  -> Draw("PC"); 
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    drawText(0.70, 0.83,"S+B Fit", kRed, 0.032);
    drawText(0.70, 0.79,"B-only Fit", kBlack, 0.032);
    canvas -> Print(pdf); 

    //==// Chisquare / NDF 
    TG_ChiNdf_SB = makeTGraph(v_MassNum, v_ChisquareNdf_SB, "#chi^{2}/NDF", kRed);  
    TG_ChiNdf_SB -> Draw("APC"); 
    TG_ChiNdf_B  = makeTGraph(v_MassNum, v_ChisquareNdf_B, "#chi^{2}/NDF"); 
    TG_ChiNdf_B  -> Draw("PC"); 
    makeTLine(config.start, 1, config.end, 1, kBlack, "", 1, 2);
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    drawText(0.70, 0.83,"S+B Fit", kRed, 0.032);
    drawText(0.70, 0.79,"B-only Fit", kBlack, 0.032);
    canvas -> Print(pdf); 

    //==// Chisquare p-value
    canvas -> SetLogy();
    TG_ChiNdfpVal_SB = makeTGraph(v_MassNum, v_ChisquarepVal_SB, "#chi^{2} p-value", kRed);
    TG_ChiNdfpVal_SB -> Draw("APC"); 
    TG_ChiNdfpVal_SB -> SetMaximum(1); 
    TG_ChiNdfpVal_B  = makeTGraph(v_MassNum, v_ChisquarepVal_B, "#chi^{2} p-value");
    TG_ChiNdfpVal_B  -> Draw("PC"); 
    drawText(0.70, 0.30,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    drawText(0.70, 0.26,"S+B Fit", kRed, 0.032);
    drawText(0.70, 0.22,"B-only Fit", kBlack, 0.032);
    canvas -> Print(pdf); 
    canvas -> SetLogy(0);

    //==// LLH SB and B-only 
    TG_LLH_SB = makeTGraph(v_MassNum, v_LogLikeHood_SB, "LogLikeHood", kRed);
    TG_LLH_SB -> Draw("APC");
    TG_LLH_B = makeTGraph(v_MassNum, v_LogLikeHood_B, "LogLikeHood");
    TG_LLH_B -> Draw("PC"); 
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    drawText(0.70, 0.83,"S+B Fit", kRed, 0.032);
    drawText(0.70, 0.79,"B-only Fit", kBlack, 0.032);
    canvas -> Print(pdf); 

    //==// LLHR
    TG_LLH = makeTGraph(v_MassNum, v_LogLikeHood, "log [ L(s+b)/L(b) ]");
    TG_LLH -> Draw("APC");
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    canvas -> Print(pdf); 

    //==// LLHR p-value
    canvas -> SetLogy();
    TG_LLHpval = makeTGraph(v_MassNum, v_LogLikeHoodPvalue, "Local p-value"); 
    TG_LLHpval -> SetMaximum(1);
    TG_LLHpval -> Draw("APC");
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    canvas -> Print(pdf); 
    canvas -> SetLogy(0);  

    //==// get LLHR and p-value for positive and negative signal 
    vector<double> v_LogLikeHood_p, v_LogLikeHood_n; 
    vector<double> v_LogLikeHoodpval_p, v_LogLikeHoodpval_n; 
  
    for(int i = 0; i < v_ExtractSig.size(); i++)
    { 
      if (v_ExtractSig[i] < 0) 
      {
        v_LogLikeHood_n.push_back( v_LogLikeHood[i] ); 
        v_LogLikeHood_p.push_back( 0 );
        v_LogLikeHoodpval_n.push_back( v_LogLikeHoodPvalue[i] );
        v_LogLikeHoodpval_p.push_back( 1 );
        //LogLikeHood[i] = 0; 
      } 
      else 
      { 
        v_LogLikeHood_n.push_back( 0 ); 
        v_LogLikeHood_p.push_back( v_LogLikeHood[i] );
        v_LogLikeHoodpval_n.push_back( 1 );
        v_LogLikeHoodpval_p.push_back( v_LogLikeHoodPvalue[i] );
      }
    }

    //==// LLHR 
    TG_LLH_positive = makeTGraph(v_MassNum, v_LogLikeHood_p, "log [ L(s+b)/L(b) ]", kRed);
    TG_LLH_positive -> Draw("APC");
    TG_LLH_negative = makeTGraph(v_MassNum, v_LogLikeHood_n, "log [ L(s+b)/L(b) ]"); 
    TG_LLH_negative -> Draw("PC"); 
    drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    drawText(0.70, 0.83,"Positive Signal", kRed, 0.032);
    drawText(0.70, 0.79,"Negative Signal", kBlack, 0.032);
    canvas -> Print(pdf); 

    //==// LLHR p-value 
    canvas -> SetLogy(); 
    TG_LLHpval_positive = makeTGraph(v_MassNum, v_LogLikeHoodpval_p, "Local p-value", kRed); // p-value 
    TG_LLHpval_positive -> Draw("APC"); 
    TG_LLHpval_negative = makeTGraph(v_MassNum, v_LogLikeHoodpval_n, "Local p-value"); // p-value 
    TG_LLHpval_negative -> Draw("PC"); 
    drawText(0.70, 0.30,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
    drawText(0.70, 0.26,"Positive Signal", kRed, 0.032);
    drawText(0.70, 0.22,"Negative Signal", kBlack, 0.032);
    canvas -> Print(pdf); 
    canvas -> SetLogy(0); 

  
    //==// nominal fit bkg parameters
    for (int i = 0; i < config.par.size(); i++)
    { 
      TG_bkgParam_Bnom[i] = makeTGraph(v_MassNum, vv_BkgParam_Bnom[i] , to_string(config.par.size())+" Param Function: p"+to_string(i), kBlack); 
      TG_bkgParam_Bnom[i] -> Draw("APC");
      TG_bkgParam_SBnom[i] = makeTGraph(v_MassNum, vv_BkgParam_SBnom[i] , to_string(config.par.size())+" Param Function: p"+to_string(i), kRed); 
      TG_bkgParam_SBnom[i] -> Draw("PC");
      setBestRangeTGraph(TG_bkgParam_Bnom[i], TG_bkgParam_SBnom[i]);
      makeTLine(config.start, paramNomGlob[i], config.end, paramNomGlob[i], kBlack, "", 1, 2);
      drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
      drawText(0.70, 0.83,"Nominal B-only", kBlack, 0.032);
      drawText(0.70, 0.79,"Nominal S+B", kRed, 0.032);
      canvas -> Print(pdf);
    }

    //==// alternate fit bkg parameters
    for (int i = 0; i < config.par2.size(); i++)
    { 
      TG_bkgParam_Balt[i] = makeTGraph(v_MassNum, vv_BkgParam_Balt[i] , to_string(config.par2.size())+" Param Function: p"+to_string(i), kBlack); 
      TG_bkgParam_Balt[i] -> Draw("APC");
      TG_bkgParam_SBalt[i] = makeTGraph(v_MassNum, vv_BkgParam_SBalt[i] , to_string(config.par2.size())+" Param Function: p"+to_string(i), kRed); 
      TG_bkgParam_SBalt[i] -> Draw("PC");
      setBestRangeTGraph(TG_bkgParam_Balt[i], TG_bkgParam_SBalt[i]);
      makeTLine(config.start, paramAltGlob[i], config.end, paramAltGlob[i], kBlack, "", 1, 2);
      drawText(0.70, 0.87,"#bf{#it{ATLAS}} Internal", kBlack, 0.035);
      drawText(0.70, 0.83,"Alternate B-only", kBlack, 0.032);
      drawText(0.70, 0.79,"Alternate S+B", kRed, 0.032);
      canvas -> Print(pdf);
    }


    //==// NP pull plots
    if (config.addSys)
    { 
      vector<double> zeros(v_MassNum.size()   ,  0.0); 
      vector<double> sig1Up(v_MassNum.size()  ,  1.0);
      vector<double> sig1Down(v_MassNum.size(), -1.0); 
      vector<double> sig2Up(v_MassNum.size()  ,  2.0);
      vector<double> sig2Down(v_MassNum.size(), -2.0); 
   
      //==// from best-fit 
      for (int i = 0; i < vv_NP_values.size(); i++ ) // loop over NPs
      {
        TG_2sig  = makeBandTGraph(v_MassNum, sig2Up, sig2Down, kYellow, " NP "+to_string(i)); 
        TG_2sig -> Draw("ALF"); 
        TG_1sig  = makeBandTGraph(v_MassNum, sig1Up, sig1Down, kGreen+1, " NP "+to_string(i)); 
        TG_1sig -> Draw("LF"); 
        TG_NP    = makeTGraphErrors(v_MassNum, vv_NP_values[i], zeros, vv_NP_errors[i], " NP "+to_string(i)); 
        TG_NP   -> Draw("P");
        makeTLine(config.start, 0.0, config.end, 0.0, kBlack, "", 1, 1); 
        canvas  -> Print(pdf); 
      }

      //==// from limits 
      for (int i = 0; i < vv_NP95_values.size(); i++ ) // loop over NPs
      {    
        TG_2sig  = makeBandTGraph(v_MassNum, sig2Up, sig2Down, kYellow, " NP "+to_string(i)); 
        TG_2sig -> Draw("ALF"); 
        TG_1sig  = makeBandTGraph(v_MassNum, sig1Up, sig1Down, kGreen+1, " NP "+to_string(i)); 
        TG_1sig -> Draw("LF"); 
        TG_NP    = makeTGraphErrors(v_MassNum, vv_NP95_values[i], zeros, vv_NP95_errors[i], " NP "+to_string(i)); 
        TG_NP   -> Draw("P");
        makeTLine(config.start, 0.0, config.end, 0.0, kBlack, "", 1, 1);  
        canvas  -> Print(pdf); 
      } 
    } 

    //==// close canvas 
    canvas->Print(pdf+"]");

    histRawSS -> GetXaxis() -> SetRangeUser(config.xminfit, config.xmaxfit);

    //==// print out histograms to root file 
    rootFile -> cd();   
    histRawSS           -> Write("sigSubData"); 
    bkgBestHist         -> Write("swiftBkgBest");  
    bkgSSHist           -> Write("swiftSSBkg");                                                                                           
    BHistGlobal         -> Write("globalFitNominal");
    BHist2Global        -> Write("globalFitAlternate");
    TG_paramChoosen     -> Write("parameterChoosen"); 
    TG_windowSizeX      -> Write("windowSizeX"); 
    TG_windowCenter     -> Write("windowCenter"); 
    TG_windowLow        -> Write("windowLow"); 
    TG_windowHigh       -> Write("windowHigh"); 
    TG_WindowCorr       -> Write("windowCorr");
    TG_Signal           -> Write("ExtractedSignal");     
    TG_SignalError      -> Write("ExtractedSignalError");
    TG_Signal95         -> Write("Extracted95Signal");     
    TG_Limit95          -> Write("Limit95"); 
    TG_Signifi          -> Write("Significance");     
    TG_Chi_SB           -> Write("ChiSq_SB");
    TG_Chi_B            -> Write("ChiSq_B");     
    TG_ChiNdf_SB        -> Write("ChiSqNDF_SB");     
    TG_ChiNdf_B         -> Write("ChiSqNDF_B");
    TG_ChiNdfpVal_SB    -> Write("ChiSquaredPvalue_SB");
    TG_ChiNdfpVal_B     -> Write("ChiSquaredPvalue_B");
    TG_LLH_SB           -> Write("LLH_SB");     
    TG_LLH_B            -> Write("LLH_B");    
    TG_LLH              -> Write("LLHR");
    TG_LLHpval          -> Write("LLHRpVal");
    TG_LLH_positive     -> Write("LLHR_positive");  
    TG_LLH_negative     -> Write("LLHR_negative"); 
    TG_LLHpval_positive -> Write("LLHRpVal_positive");
    TG_LLHpval_negative -> Write("LLHRpVal_negative");
    for (int i = 0; i < config.par.size(); i++)
    { 
      TG_bkgParam_Bnom[i]  -> Write(("bkgParam"+to_string(i)+"_Bnom").c_str()); 
      TG_bkgParam_SBnom[i] -> Write(("bkgParam"+to_string(i)+"_SBnom").c_str()); 
    }  
    for (int i = 0; i < config.par2.size(); i++)
    { 
      TG_bkgParam_Balt[i]  -> Write(("bkgParam"+to_string(i)+"_Balt").c_str()); 
      TG_bkgParam_SBalt[i] -> Write(("bkgParam"+to_string(i)+"_SBalt").c_str()); 
    }
    rootFile -> Close(); 
    delete rootFile;    

  }
    
  //==// End of SWiFt 
  cout << endl; 
  cout << "Output files : " << endl; 
  cout << "          PDF: " << pdfFileName  << endl; 
  cout << "         ROOT: " << rootFileName << endl; 
  printStuff("../ascii/zap.txt"); 
  cout<<"CHI2 p-value of nominal global fit bkg "<<value1<<endl;
  cout<<"CHI2 p-value of alternate global fit bkg "<<value2<<endl;
  cout<<"CHI2 p-value of SWift bkg "<<value3<<endl;

 
}



