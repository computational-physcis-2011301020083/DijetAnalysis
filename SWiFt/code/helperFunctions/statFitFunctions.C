////////////////////////////////////////////////////////////////////////////////                                                                       
// 
// Collection of statistical, goodness of fit, etc functions used by SWiFt
//   
////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace TMath;

#include "../config.C"
#include "signalGenerator.C"

bool isParamFixed (TF1* function, int i); 

double minuitFit ( TH1D* hist, TF1 *function, int funcType, double min, double max, bool quietFitInfo,  
                 vector<double> *fitQuality, vector<vector<vector<double>>> *diffVec, vector<double> *NP_values, vector<double> *NP_errors);

// global TH1 and TF1 for Minuit
TH1D *theHist;
TF1 *theFunc;
TF1 *theFuncNom; 
TF1 *theBestFunc; 
int theFuncType; 
double theXmin;
double theXmax;
double theLLH; 
vector<vector<vector<double>>> *theDiffVec; 
double width[15] = {0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17};

////////////////////// refit function 
TFitResultPtr reFit(TFitResultPtr r, TH1D *hist, TF1 *func, TString label)
{
  //cout << "============================================= "+label+" fit success before: " << r->IsValid() << endl;
  for (int n = 0; n < 50; n++)
  {
    //r = hist -> Fit(func, "SRLIE");
    int status = int ( r );
    //cout << "TOOT:::" << status << endl;
    //cout << "TOOT2:::" << r->MinFcnValue() << endl;
    if (r->IsValid() == 1 && r->MinFcnValue() < 1000000 ) break;
    r = hist -> Fit(func, "SRLI");
    cout << "============================================= "+label+" fit success after: " << r->IsValid() << endl;
  }
  return r;
}

////////////////////// auto refit                                                                                                           
void autoFudge(TFitResultPtr r, TH1D *hist, TF1 *func, double refit, double delta)
{ 

  cout << "Min ChiSq initial: " << r->MinFcnValue() << endl;
  double bestFit = r->MinFcnValue();     

  const int numPar = 6; // 9

  // fit parameters arrays
  double param[numPar] = 
  {
    func->GetParameter(0), // normalization 
    func->GetParameter(1), // gaussian mean
    func->GetParameter(2), // gaussian width 
    func->GetParameter(3), // fraction
    func->GetParameter(4), // landau mean 
    func->GetParameter(5), // landau width
    //func->GetParameter(6), // normalization 2 - W* only 
    //func->GetParameter(7), // gaussian mean 2 - W* only
    //func->GetParameter(8), // gaussian width 2 - W* only
    
  }; 
  double randomParam[numPar] = {0,0,0,0,0,0};   
  double bestParam[numPar] = {0,0,0,0,0,0}; 

  // initialize 
  for (int j = 0; j < numPar; j++) 
    bestParam[j] = param[j]; 

  // number of tries 
  for (int i = 0; i<refit; i++)
  {
  
    // get randomly twiked parameters and set parameters in function 
    TRandom3 random(i); 
    for (int j = 0; j < numPar; j++)
    { 
      // skip random twiddle if param is fixed 
      if ( isParamFixed (func,j) ) continue;  

      //// uniform random numbers in a range 
      //double start = param[j]-(delta*param[j]) ; 
      //double end   = param[j]+(delta*param[j]) ; 
      //randomParam[j] = random.Uniform(start, end);

      //// gaussian random numbers 
      randomParam[j] = random.Gaus( abs(param[j]), delta*abs(param[j]) );
      if (j==4) 
        randomParam[j] = -randomParam[j]; // negative for landau mean 

      // re-set function paraeters
      func -> SetParameter(j, randomParam[j]);   

    }

    // refit histogram 
    TFitResultPtr rNew = hist -> Fit(func, "SRQ");
    cout << "Min ChiSq final: " << rNew->MinFcnValue() << endl;
    double tempFit = rNew->MinFcnValue();

    // pass fit results out of function 
    //if (tempFit < bestFit) 
    if( (tempFit < bestFit) && abs(tempFit - bestFit) > 0.01 )
    {
      bestFit = tempFit;
      for (int j = 0; j < numPar; j++)
      {
        bestParam[j] = func->GetParameter(j);
        cout << bestParam[j] << endl;
      }
    }

  }

  // reset parameters to best 
  for (int j = 0; j < numPar; j++)
  { 
    if ( isParamFixed (func,j) ) continue;
    func -> SetParameter(j, bestParam[j]);
  }

  // fit with best parameters
  //TFitResultPtr rNeww = hist -> Fit(func, "SRIQ"); 

}

////////////////////// calculate the loglikehood of TH1 and TF1 
double calculateLogLikeHood(TH1D *hist, double xmin, double xmax, TF1 *func)
{   
  // find number of bins and lowest bin in histogram  
  int nBins = hist->GetXaxis()->FindBin(xmax) - hist->GetXaxis()->FindBin(xmin); 
  //int nBins = (xmax-xmin)/config.rebin;
  int lowBin = hist->FindBin(xmin);
  double totalLogLikeHood = 0;
  
  // loop over all bins in histogram
  for(int i = 0; i < nBins; i++)
  { 
    double Bin_xmin = hist -> GetBinLowEdge(lowBin); 
    double Bin_xmax = hist -> GetBinLowEdge(lowBin+1);
    double center   = hist -> GetBinCenter(lowBin);
   
    // for each bin, calculate the LLH
    double n           = hist->GetBinContent(lowBin); 
    double lambda      = func -> Integral(Bin_xmin, Bin_xmax);  
    //if (lambda < 1e-10) lambda = 1e-10;
    double logLikeHood = n*log(abs(lambda)) - LnGamma(n+1) - lambda;   

    double cutoff = 1e-4;                                                                                                                             
    if (lambda < cutoff)
      logLikeHood += lambda/cutoff - 1;

    //double cutOff      = 1e-12; 
    //if (lambda > cutOff) logLikeHood = n*log(lambda) - LnGamma(n+1) - lambda;
    //else 
    //{
    //  double logL        = (log(cutOff)+(lambda/cutOff-1) - 0.5*(lambda/cutOff-1)*(lambda/cutOff-1));  
    //  logLikeHood = n*logL - LnGamma(n+1) - exp(logL);
    //}

    //cout << "bin: " << Bin_xmin << ", lambda: " << lambda << ", LLH: " << logLikeHood << endl; 
    //cout << "Bin: " << Bin_xmin << ", " << n << ", " << lambda << endl; 

    // add all the LLH for all bins 
    totalLogLikeHood = totalLogLikeHood + logLikeHood;
    
    lowBin++;
  }
  
  // return LLH for histogram("a window") 
  return totalLogLikeHood;
}

////////////////////// LLH ratio to p-value 
double LLHRPval(double LLHR, int NDF)
{ 
  LLHR = 2*LLHR; 
  return TMath::Prob(LLHR, NDF);
}


////////////////////// LLH scan function 
double LLHScan(double lnLikehoodSB, TF1 *SBfunc, TF1 *SBFunc95, TF1 *BfuncLocal, TF1 *sigFunc, double signal, double signalErr, 
               TH1D *dHistRaw, double low, double high, double center, vector<vector<vector<double>>> *diffVec = NULL, 
               vector<double> *NP_values_window = NULL, vector<double> *NP_errors_window = NULL)
{
  //==// get bkg parameters 
  int sigParam = sigFunc->GetNpar(); 
  int bkgParam = BfuncLocal->GetNpar(); 
  vector<double> v_bkgParam(bkgParam);
  for (int i = 0; i < bkgParam; i++)  
    v_bkgParam[i] = SBfunc -> GetParameter(sigParam+i); 
 
  //==// define variables
  double LLH95     = -lnLikehoodSB+1.92073;      // target likelihood, i.e. ~2sigma worse than original 
  double signal95  = -123456. ;                  // signal events that correspond to LLH95
  double sigTemp   = signal+abs(4*signalErr);    // initial signal estimate for binary search  
  double unit      = abs(2*signalErr);           // step size for binary search 
  double sigRaw    = SBfunc -> GetParameter(0);  // raw signal from best s+b fit 

  //==// recalculate LLH for negative signals 
  if (signal < 0)
  {   
    LLH95 = -calculateLogLikeHood(dHistRaw, low, high, BfuncLocal)+1.92073;
    sigTemp = abs(4.5*signalErr);
    for (int i = 0; i< bkgParam; i++) 
      v_bkgParam[i] = BfuncLocal -> GetParameter(i); 
  }   
  
  //==// print out info
  cout << "Window center    : " << center << endl; 
  cout << "Best s+b fit LLH : " << -lnLikehoodSB << endl;
  cout << "Signal extracted : " << signal << " (raw: " << sigRaw << " )" << endl;  
  cout << "95% LLH to reach : " << LLH95  << endl;
  cout << "unit tolerance   : " << config.tol << endl;     

  //==// start binary search for number of events that make the best LLH worse by ~2sigma 
  int lsum = -1; 
  int rsum = -1;
  double LLHTemp = 99999999.0; 
  for (int i = 0; i < 200; i++) 
  { 

    bool left = 0; 
    bool right = 0; 
     
    //==// sigTemp always needs to be positive 
    if (sigTemp < 0) sigTemp = 1.0;

    //==// fix normalization and reset bkg parameters 
    SBFunc95 -> FixParameter(0, sigTemp*(sigRaw/signal));  // normalization
    for (int i = 1; i < sigParam; i++) 
      SBFunc95 -> FixParameter(i, sigFunc->GetParameter(i));
    for (int i = 0; i < bkgParam; i++) 
      SBFunc95 -> SetParameter(sigParam+i, v_bkgParam[i]);

    //==// perform fit 
    vector<double> Lim_fitQuality;
    LLHTemp = minuitFit ( dHistRaw, SBFunc95, 1, low, high, config.quiet95FitInfo, &Lim_fitQuality, diffVec, NP_values_window, NP_errors_window);
    cout << "Iteration " << i << ", LLHTemp: " << LLHTemp << ", LLH95: " << LLH95 << ", sigTemp: " << sigTemp << ", unit: " << unit <<  endl;

    //==// check position of LLHTemp w.r.t. LLH95
    if(LLHTemp < LLH95)
      { lsum++;  left = 1; }  
    else
      { rsum++;  right = 1; }

    //==// if unit in less than config.tol, break the loop 
    if (abs(LLHTemp - LLH95) < config.tol)
    {  
      signal95 = sigTemp;   
      break; 
    } 

    //==// decrease unit by 2 if LLHTemp crosses LLH95, i.e. when the sign changes
    // if the sign does not change, decrease unit linearly 
    // once the sign changes, reduce unit by half with each loop iteration 
    if (left) 
    { 
      if (i!=lsum) 
        unit=unit/2;
       else                                                                                                                                     
        unit = unit*1.5;
      sigTemp += unit;
    } 
    if (right) 
    { 
      if (i!=rsum) 
        unit=unit/2;
      else
        unit = unit*1.5;
      sigTemp -= unit;
    }

    //==// kick signal tmp
    if (i == 100 && abs(LLHTemp - LLH95) > config.tol)                                                                                              
      unit = abs(2*signalErr);  
  
    //==// pass bkg fit values to next fit
    for (int i = 0; i < bkgParam; i++) 
      v_bkgParam[i] = SBFunc95 -> GetParameter(sigParam+i); 

  }

  cout << "====================================================>>>>>>>>>>> 95% Signal: " << signal95 << endl; 
  cout << "====================================================>>>>>>>>>>> abs(LLHTemp - LLH95): " << abs(LLHTemp - LLH95) << endl; 

  //==// return 95% upper limit on extracted signal
  return signal95; 
}


////////////////////// Do KS test 
double KSTest(TH1D *Data, TH1D *MC, double min, double max)
{ 

  // get bins of min and max 
  int bMin = Data -> FindBin(min);
  int bMax = Data -> FindBin(max); 

  // re-create TH1Ds for given ranges: min and max 
  TH1D *h_data = new TH1D("h_data", "Data", bMax - bMin , min, max);
  TH1D *h_mc   = new TH1D("h_mc"  , "MC"  , bMax - bMin , min, max);

  // fill re-created histograms with incoming histograms 
  for (int i = 1; i <= h_data->GetNbinsX(); i++) 
  { 
    if ( Data->GetBinContent(bMin-1+i) < 1 )
      h_data -> SetBinContent(i, 0);
    else 
      h_data -> SetBinContent(i, Data->GetBinContent(bMin-1+i)); 
    //cout << "Data entry bin: " << bMin-1+i << " , BinContent: " << h_data->GetBinContent(i) << endl; 
  } 

  for (int i = 1; i <= h_mc->GetNbinsX(); i++) 
  { 
    h_mc -> SetBinContent(i, MC->GetBinContent(bMin-1+i)); 
    //cout << "MC entry bin: " << bMin-1+i << " , BinContent: " << MC->GetBinContent(bMin-1+i) << endl; 
  } 


  double KS = h_mc -> KolmogorovTest(h_data, "X"); 

  h_data -> Delete(); 
  h_mc   -> Delete(); 

  return KS; 
} 

////////////////////// Compute Mahalanobis distance 
double MahalDist(TMatrixDSym S, vector<double> globalParam, vector<double> localParam)
{ 
  // invert the sum of the matrices    
  S.Invert();

  // generate difference between the fit parameters of local and global bkgs 
  vector<double> diff; 
  for (int i = 0; i < globalParam.size(); i++)
    diff.push_back(globalParam[i]-localParam[i]);  

  // calulate the Mahalanobis distance with some matrix maniplulation 
  double MDist = 0; 
  for (int i = 0; i < globalParam.size(); i++)
  { 
    double temp = 0; 
    for (int j = 0; j < globalParam.size(); j++)
      temp = temp + S[i][j] * diff[j]; 
    MDist = MDist + temp*diff[i]; 
  } 

  // output the distance
  return sqrt(MDist); 
}

////////////////////// Do Chisquare/NDF test 
double ChiSqNDFTest(TH1D *data, TH1D *hist, double min, double max, double NDF=-1, bool getChiSq = 0, bool getChiSqNDF = 0, bool getpValue = 1, double offset = 0.0 )
{ 
  double ChiSq = 0.0; 

  // get bins of min and max 
  int bMin = data -> FindBin(min);
  int bMax = data -> FindBin(max); 
  int bMinH = hist -> FindBin(min);
  int bMaxH = hist -> FindBin(max); 

  /*
  cout << "min, max: " << min << ", " << max << endl; 
  cout << "bMin, bMax: " << bMin << ", " << bMax << endl; 
  cout << "bMinH, bMaxH: " << bMinH << ", " << bMaxH << endl;
  cout << "bMax-bMin: " << bMax-bMin << endl; 
  */

  // calculate chisquare
  for (int i = 0; i < (bMax-bMin)+1; i++) 
  { 
    double ChiTmp = (data->GetBinContent(i+bMin) - hist->GetBinContent(i+bMinH)) * (data->GetBinContent(i+bMin) - hist->GetBinContent(i+bMinH)); 
    ChiTmp = ChiTmp/hist->GetBinContent(i+bMin);
    //cout <<"i: "<< i << ", Mass: " << data->GetBinLowEdge(i+bMin) << ", ChiTmp: " << ChiTmp << endl; 
    //if (data->GetBinContent(i+bMin) < 0) continue;   
    ChiSq = ChiSq + ChiTmp;        
  }

  // penalize chiSq by offset 
  ChiSq = ChiSq + offset; 
  //ChiSq = ChiSq + (bMax-bMin);  

  double ChisqNDF = ChiSq/NDF; 
  double ChipVal  = TMath::Prob(ChiSq, NDF);
  //cout << "ChiSq: " << ChiSq << ", NDF: " << NDF << ", p-value: " << ChipVal << endl; 

  if (getChiSq) return ChiSq; 
  if (getChiSqNDF) return ChisqNDF; 
  if (getpValue) return ChipVal; 
}

////////////////////// Wilks test 
double WilksTest( TH1D* hist, TF1* nominalF, TF1* alternateF, double min, double max)
{ 

  double LLH_nominal   = calculateLogLikeHood(hist, min, max, nominalF); 
  double LLH_alternate = calculateLogLikeHood(hist, min, max, alternateF);
  int NDF              = alternateF->GetNumberFreeParameters() - nominalF->GetNumberFreeParameters();  

  double LLH           = 2*(LLH_alternate-LLH_nominal); 
  double wilks         = 1 - ( Gamma(0.5, LLH/2.0)/Gamma(0.5, 1e20)); 
  //double wilks         = TMath::Prob(LLH, NDF);

  return wilks;  

}

////////////////////// do pseudo-experiment 
TH1D *doPseudoExp(TH1D* hist, int seed = 1)
{ 
  TRandom3 s(seed);
  TH1D* tmpHist = (TH1D*) hist -> Clone("tmpHist");                                                                                                              
  for (int i = 0; i < hist -> GetSize(); i++)
  {   
    double bin_cont = hist -> GetBinContent(i);
    double random   = s.Poisson(bin_cont); 
    tmpHist -> SetBinContent(i,random);                                                                                                                        
    tmpHist -> SetBinError(i, sqrt(random)); 
  }   

  return tmpHist; 

} 


/////////////////////// solve quadratic equation 
double solveQuad( double a, double b, double c )
{ 
  a = a/c; 
  b = b/c; 
  c = 1;  
  double sq = sqrt (b*b - (4*(a-0.5)*c)); 
  double pos = (-b + sq)/2*c; 
  double neg = (-b - sq)/2*c;
 
  return (pos-neg)/2;  
}

//////////////////////// print out fails 
void reportFail( vector<vector<double>> vv_fitQuality )
{
  for (int i = 0; i < vv_fitQuality.size(); i++)
  { 
    for (int j = 0; j < vv_fitQuality[i].size(); j++)
    {
      if (j == (vv_fitQuality[i].size())-1)
        cout << vv_fitQuality[i][j] << endl; 
      else 
        cout << vv_fitQuality[i][j] << ", " ;    
    } 
  } 

}

//////////////////////// get number of free parameters - redefinition since there is a bug in the official ROOT function 
int GetNumberFreeParameters(TF1* func) 
{ 
  double al, bl;
  int nfree = func -> GetNpar();                                                                                                                           
  int ntot = func -> GetNpar();
  for (Int_t f = 0; f < ntot; f++) 
  { 
    func->GetParLimits(f, al, bl); 
    if (al * bl != 0 && al >= bl) nfree--;
  }    

  return nfree; 
}

////////////////////// -LLH minimization 
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, int iflag) 
{
 
  int numTotFitParam  = theFunc->GetNpar();               // total number of parameters in TF1 
  int numFreeParam    = GetNumberFreeParameters(theFunc); // number of free parameters in TF1
  int numNPParam      = npar-numFreeParam;                // number of nuisance parameters - 1 for each uncertainty
  double LLH          = 0.0;   
  
  //======================================================
  //
  // NO signal systematics : always true for bkg only fits 
  // Pass minuit parameters to TF1 
  //
  //======================================================
  for (int i = 0; i<numTotFitParam; i++) 
    theFunc  -> SetParameter(i, par[i]);

  //======================================================
  //
  // WITH signal systematics: only for the s+b fit  
  // Modify parameters, i.e., signal normalization and 
  // signal shape parameters to take into account uncertainties 
  // Add NPs to LLH 
  //
  //======================================================
  if (config.addSys == 1 && theFuncType == 1)
  { 
  
    //==// combine flat uncertaities in quadrature 
    double flatUncrt = config.scaleUn;
    if (config.do95Scan) flatUncrt = sqrt( (config.scaleUn*config.scaleUn) + (config.lumiUn*config.lumiUn) );  
    
    //==// modify signal normalization as per the flat uncertainty
    //* par[0] = signal normalization 
    //* par[numTotFitParam] = NP for flat uncertainty
    double modifiedNorm = par[0] + (par[0] * par[numTotFitParam] * (flatUncrt/100.)); 
    theFunc -> SetParameter( 0, modifiedNorm );

    //==// modify signal shape parameters as per the JES uncertainties
    //* function generateSysSignal()  modifies the signal parameters of theFunc
    generateSysSignal(theDiffVec, theFunc, theFuncNom, par);
   
    //==// add gaussian penalty terms to likelihood for each uncertainty
    //* gaussian penelty term is of form: (NP*NP)/2
    for (int i = numTotFitParam; i<numTotFitParam+numNPParam; i++)  
      LLH = LLH + (par[i]*par[i])/2 ;
   
  } 

  //======================================================
  //
  // Calculate final likelihood 
  // This is what minuit minimizes (f)
  //
  //======================================================
  LLH = LLH + (-calculateLogLikeHood(theHist, theXmin, theXmax, theFunc));
  f = LLH; 
  if (LLH < theLLH) // same best likelihood and function parameters 
  {
    theLLH = LLH;
    delete theBestFunc; 
    theBestFunc = (TF1*)theFunc->Clone();
    //theBestFunc = theFunc; 
    //cout << "=====================================================================This is the BEST" << endl; 
  }

  //if (config.addSys == 1 && theFuncType == 1) 
  //{ 
    //cout << "LLH:::::::::::::::::::: " << LLH << endl; 
    //for ( int i = 0; i<numTotFitParam+numNPParam; i++ ) { 
    //  cout << "Best fit Param " << par[i] << endl; 
    //  //cout << "Fixed Param    " << theBestFunc -> GetParameter(i) << endl; 
    //} 
    //cout << endl;  
  //}

}


////////////////////// check if parameter is fixed 
bool isParamFixed (TF1* function, int i)
{
  double al, bl; 
  function -> GetParLimits(i,al,bl);
  return (al*bl != 0 && al >= bl); 
}


////////////////////// initialize parameters 
void initialize(TMinuit *gMinuit, TF1 *function, int numTotFitParam = 0, int ierflg = 0, int numNP = 0)
{
  //==// tell minuit about all the parameters in the TF1  
  for (int i = 0; i<numTotFitParam; i++)
  { 
    // get parameter value and set step size 
    double ini     = function->GetParameter(i);
    double stepval = (ini != 0) ? abs(ini) * 1e-6 : 1.0;  

    // get limits on parameters from input TF1 
    double a, b; 
    function->GetParLimits(i,a,b); 

    // initialize fixed and free paramters 
    if ( isParamFixed(function,i) ) // fixed 
    { 
      gMinuit->mnparm(i, "p"+to_string(i), ini, stepval, 0, 0, ierflg);
      gMinuit->FixParameter(i);
    }
    else // freeeeeee  
    {
      gMinuit->mnparm(i, "p"+to_string(i), ini, stepval, a, b, ierflg); 
      gMinuit->Release(i);
    }
  }

  //==// tell minuit about the nuisance parameters for the systematics 
  if (config.addSys)
  { 
    // loop over the number of nuisance parameters 
    for (int i = numTotFitParam; i<numTotFitParam+numNP; i++)
    { 
      gMinuit->mnparm(i, "np"+to_string(i-numTotFitParam), 0.0, 1e-6 , 0.0 , 0.0, ierflg);
      gMinuit->Release(i);  
    } 
  } 

 
}


////////////////////// fudge 
void fudge(TF1 *function, TF1 *functionLim, TRandom3 *random, int numTotFitParam = 0, int funcType = 0, double wid = 0.02)
{ 
  
  //cout << "Parameters that are fudged: " << endl; 
  //for (int i = 0; i < numTotFitParam; i++)
  //  cout << function -> GetParameter(i) << endl; 
 

  //==// loop over all parameters 
  for (int i = 0; i < numTotFitParam; i++)
  { 
    //==// skip fudging fixed parameter 
    if ( isParamFixed (function,i) ) continue;

    //==// get limits  
    double start = -999; 
    double end   = -999; 
    functionLim  -> GetParLimits(i, start, end); 
    //cout << "Parameter " << i << " [ " << start << ", " << end << " ]" << endl; 

    //======================================================
    //
    // For s+b functions
    //
    //======================================================
    /*
    if (funcType == 1)  
    { 
      //==// set the normalization to what is in the config file 
      if (i == 0)
        function -> SetParameter(i, config.norm);
      else //==// fudge the bkg components using a gaussian distribution with center = parameter, wid = 15%*parameter 
      {  
        double param      = function -> GetParameter(i); 
        double paramPrime = random   -> Gaus( param, 0.15*abs(param) );
        function -> SetParameter(i, paramPrime);
      }

      //==// do some fixed width gaussians 
      //if (config.freeFloat) 
      //  function -> FixParameter(3, wid*function->GetParameter(2));
    }
    */

    //======================================================
    //
    // For bkg only functions    
    //
    //======================================================
    // make sure that the parameter limits in the config are not set to zero 
    // then initialize fits from a uniform distribution drawn from the parameter limits   
    // and scale the normalization based on the number of events in the sliding window 
    if (funcType == 0) 
    {
      //==// check the parameter limits in the config. They should not be zero
      if (start == 0 && end == 0) 
      {
        cout << endl; 
        printStuff("../ascii/crab.txt");
        cout << "Parameter limits in the config are set to zero. Please change. " << endl; 
        cout << endl; 
        exit(0); 
      } 

      //==// initialize parameters by selecting random numbers from a uniform range specified by the parameter limits 
      double paramPrime = random -> Uniform(start, end);
      function -> SetParameter( i, paramPrime); 

      //==// scale normalization based on the number of events in the sliding window
      double funcIntegral = function -> Integral(theXmin, theXmax); 
      double histIntegral = theHist  -> Integral(theHist->FindBin(theXmin), theHist->FindBin(theXmax));
      double scaleInt     = histIntegral/funcIntegral; 
      function -> SetParameter(0, function->GetParameter(0)*scaleInt); 
    } 
  }  

  
  //cout << "Fudged parameters: " << endl; 
  //for (int i = 0; i < numTotFitParam; i++)
  //  cout << function -> GetParameter(i) << endl; 
  //cout << endl; 
    

}


////////////////////// Minuit fit
double minuitFit ( TH1D* hist, TF1 *function, int funcType, double min, double max, bool quietFitInfo = 0,   
                   vector<double> *fitQuality = NULL, vector<vector<vector<double>>> *diffVec = NULL, 
                   vector<double> *NP_values = NULL, vector<double> *NP_errors = NULL   )
{ 

  // NOTE: parmaeter "funcType" defines the kind of fit. 0 = bkgOnly fit, 1 = s+b fit

  //======================================================
  //
  // Set up inputs and conditions
  //
  //======================================================

  //==// pass function, hist, max, min to global parameters so that fcn() can see them 
  theHist     = (TH1D*)hist->Clone(); 
  theFunc     = (TF1*)function->Clone();
  theFuncNom  = (TF1*)function->Clone();
  theBestFunc = (TF1*)function->Clone();
  theFuncType = funcType; 
  theXmin     = min; 
  theXmax     = max; 
  theDiffVec  = diffVec; 
  theLLH      = 1e9;

  //==// calculate and store number of parameters  
  int numFreeParam   = GetNumberFreeParameters(function);  // number of free parameters to fit in TF1 
  int numTotFitParam = function->GetNpar();                // nummber of total parameters in TF1 (free + fixed)
  int numFixedParam  = numTotFitParam - numFreeParam;      // number of fixed parameters 
  int numNP          = 0;                                  // number of nuisance parameters (NPs) to add to S+B fit to incorporate systematics
  if (config.addSys == 1 && funcType == 1)                 // NP only added if incorporating signal systematics && is (s+b) function 
    numNP = config.numFlatNP+config.numJESNP;    
  int numTotParam    = numTotFitParam + numNP;             // number of total parmaeters: all TF1 parameters + NPs

  //==// print out details
  if (!quietFitInfo) 
  { 
    //cout << endl; 
    cout << "###### Minuit called to do the fit ######" << endl; 
    cout << "Total parameters in function: " << numTotParam    << endl; 
    cout << "# free parameters           : " << numFreeParam   << endl; 
    cout << "# fixed parameters          : " << numFixedParam  << endl; 
    cout << "# nuisance parameters       : " << numNP          << endl;  
    cout << endl; 
  }  

  //==// clone input function and use it for the 1st fit with no parameter limits 
  TF1* functionNoLim = (TF1*)function->Clone();
  if (funcType == 0) // if bkg-only fit
    for (int i = 0; i < numTotFitParam; i++)
      functionNoLim -> SetParLimits(i, 0, 0); 
 
  //==// print out initial parameters and limits  
  if (!quietFitInfo)
  {
    cout << "Initialize fits with these parameters (NPs are always initialized to zero): " << endl; 
    cout << "-------------------------------------- " << endl; 
    double limLow  = -999; 
    double limHigh = -999;
    for (int i = 0; i<numTotFitParam; i++)
    { 
      functionNoLim->GetParLimits(i,limLow,limHigh);  
      cout << "Parameter " << i << ": " << functionNoLim->GetParameter(i) << ", [" << limLow << ", " << limHigh << "]" << endl;     
    }
    cout << endl; 
  }


  //======================================================
  //
  // Set up the minimizer 
  //
  //======================================================
  
  //==// set intergrator type 
  ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator(config.integrator); 
  
  //==// initialize minuit with 12 parameters and pass fcn() 
  TMinuit *gMinuit = new TMinuit(12); 
  gMinuit -> SetFCN(fcn); 

  //==// silence minuit output 
  if (config.quietMinuit)
    gMinuit -> SetPrintLevel(-1); 
  
  //==// arglist for error,migrand,hesse parameters 
  Double_t arglist[10];
  Int_t ierflg  = 0;
  Int_t ierflgS = 0; 
  Int_t ierflgM = 0;
  Int_t ierflgH = 0;   
  Int_t ierflgI = 0; 

  //==// minuit settings 
  arglist[0] = config.errorType;                      // 0.5 errors since we are minimizing the LLH 
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  arglist[0] = config.fitStrategy;                    // set fit strategy to 2 for better fitting 
  gMinuit->mnexcm("SET STR", arglist ,1,ierflg);
  
  arglist[0] = config.maxIterations;                  // number of maximum iterations
  arglist[1] = config.tolSimplex ;                    // simplex tolerance 

  //==// structures to store fit results 
  vector<double> parFinal(numTotParam, 0.0); 
  vector<double> parerrFinal(numTotParam, 0.0); 
  double NDF; 

  //======================================================
  //
  // Minimize !! 
  //
  //======================================================

  //==// initialize parameters for fit 
  initialize(gMinuit, functionNoLim, numTotFitParam, ierflg, numNP);
  
  if (!quietFitInfo) cout << "Ok, Minimize: " << endl; 
  if (!quietFitInfo) cout << "--------------" << endl; 
  double bestLLH = 9999999999.0;

  //==// Use SIMPLEX 
  gMinuit->mnexcm("SIMPLEX", arglist, 2, ierflgS);
  bestLLH = theLLH;
  cout << "                                                           Simplex Status: " << ierflgS << endl; 
  cout << "                                                           LLH: "            << bestLLH << endl;

  //==// Use MIGRAD
  arglist[1] = config.tolMigrad ;
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflgM);
  bestLLH = theLLH; 
  cout << "                                                           Migrad Status: " << ierflgM << endl; 
  cout << "                                                           LLH: "           << bestLLH << endl;

  //==// Use SIMPLEX & MIGRAD again if MIGRAD fails
  if (ierflgM != 0 ) 
  { 
    gMinuit->mnexcm("SIMPLEX", arglist, 2, ierflgS);
    bestLLH = theLLH;
    cout << "                                                           2nd Simplex Status: " << ierflgS << endl; 
    cout << "                                                           LLH: "                << bestLLH << endl;

    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflgM);
    bestLLH = theLLH;
    cout << "                                                           2nd Migrad Status: " << ierflgM << endl; 
    cout << "                                                           LLH: "               << bestLLH << endl;

    if (ierflgM != 0)
    {  

      gMinuit->mnexcm("SIMPLEX", arglist, 2, ierflgS);
      bestLLH = theLLH;
      cout << "                                                           3nd Simplex Status: " << ierflgS << endl; 
      cout << "                                                           LLH: "                << bestLLH << endl;

      gMinuit->mnexcm("MIGRAD", arglist, 2, ierflgM);
      bestLLH = theLLH;
      cout << "                                                           3nd Migrad Status: " << ierflgM << endl; 
      cout << "                                                           LLH: "               << bestLLH << endl;
    }  
      
  }  

  ////==// Use IMPROVE
  //gMinuit->mnexcm("IMP", arglist, 1, ierflgI);
  //bestLLH = theLLH; 
  //cout << "                                                           Improve Status: " << ierflgI << endl;  
  //cout << "                                                           LLH: "            << bestLLH << endl; 

  //==// Use HESSE
  gMinuit->mnexcm("HESSE", arglist, 1, ierflgH);
  bestLLH = theLLH; 
  cout << "                                                           HESSE Status: "   << ierflgH << endl;
  cout << "                                                           LLH: "            << bestLLH << endl; 
  cout << "                                                           # total fcn calls: "         << gMinuit->fNfcn << endl;  
  if (!quietFitInfo) cout << endl; 

  //==// push_back fit quality  
  (*fitQuality).push_back(ierflgS);  
  (*fitQuality).push_back(ierflgM);
  (*fitQuality).push_back(ierflgH);
  (*fitQuality).push_back(bestLLH); // LLH after HESSE 

  //==// get fit parameters and error on them 
  if (!quietFitInfo) cout << "Parameters after fit: " << endl; 
  if (!quietFitInfo) cout << "--------------------------" << endl; 
  int endNum = numTotFitParam; 
  if (config.addSys && theFuncType == 1) 
    endNum = numTotParam ; 

  for (int i = 0; i < endNum; i++)
  { 
    //==// get parameters and uncertainties from minuit
    gMinuit->GetParameter(i, parFinal[i], parerrFinal[i]);

    //==// get best signal shape when systematic uncertainties are added
    if (config.addSys && theFuncType == 1)
      if ( isParamFixed (theBestFunc, i) ) 
        parFinal[i] = theBestFunc -> GetParameter(i); 

    //==// print out 
    if (!quietFitInfo) cout << "Parameter " << i << ": " << parFinal[i] << " +/- " << parerrFinal[i] << endl; 
  }
  if (!quietFitInfo) cout << endl;

  if (config.addSys && theFuncType == 1)
  { 
    // NP values 
    for (int i = 0; i < numNP; i++ ) 
    {  
      (*NP_values).push_back(parFinal[numTotFitParam+i]); 
      (*NP_errors).push_back(parerrFinal[numTotFitParam+i]);
    }
  }
   
  ////////////////////////////////////////////////////////////////////
  //
  // Set parameter values, errors, ndfs for output TF1
  //
  ////////////////////////////////////////////////////////////////////

  //==// calculate NDF for the fit 
  int bins = theHist->FindBin(theXmax)-theHist->FindBin(theXmin); 
  int ndf  = bins - numFreeParam;  

  //==// pass parameters and errors to output TF1
  for (int i = 0; i<numTotFitParam; i++) 
  { 
    function  -> SetParameter(i, parFinal[i]);
    function  -> SetParError(i, parerrFinal[i]);
    function  -> SetNDF(ndf); 
  } 

  //==// delete all the clones 
  delete functionNoLim; 
  delete theHist; 
  delete theFunc; 
  delete theFuncNom;  
  delete theBestFunc;
 
  //==// return best LLH (this included the penelty terms when systematics are on) 
  return bestLLH; 
  
}







