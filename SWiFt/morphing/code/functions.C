////////////////////////////////////////////////////////////////////////////////                                                                       
// 
// Code that has all the functions used by SWiFt
// By: Karishma Sekhon (ksekhon@umich.edu)
//   
////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace TMath;

#include "config.C"

bool isParamFixed (TF1* function, int i); 

//void minuitFit ( TH1D* hist, TF1 *function, double min, double max, int numFits = 20, double center = 1000.);

// global TH1 and TF1 for Minuit
//TH1D *theHist;
//TF1 *theFunc;
//double theXmin;
//double theXmax;
//double theLLH; 
//double width[15] = {0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17};

double lineY; 

///////////////////// print stuff
void printStuff(TString file)
{ 
  string l; 
  ifstream file0(file); 
  while ( getline (file0,l) ) 
    cout << l << endl;  
  file0.close();
}

////////////////////// function to write text onto canvas                                                             
void fatal(TString msg) { printf("\nFATAL\n  %s\n\n",msg.Data()); abort(); }

////////////////////// open file 
TFile *openFile(TString fn) 
{
  cout << "File Name: " << fn << endl;
  TFile *f = TFile::Open(fn, "r");
  if (f==NULL) 
  {
    cout << endl; 
    printStuff("../ascii/cat.txt");
    fatal("YOUR FILE IS WRONG/MISSING/OUT-for-COFFEE. Cannot open file: "+fn);
  }
  return f;
}

////////////////////// get histo from TFile
TH1D* getHisto(TFile *f, TString hName)                                                            
{
  cout << "Hist Name: " << hName << endl;
  TH1D *h = (TH1D*)f->Get(hName);
  if (h==NULL) 
  {
    cout << endl; 
    printStuff("../ascii/cat.txt"); 
    fatal("Cannot access histogram from file. Does the histogram exist?");
  }
  return h; 
}

///////////////////// get histo from TDir 
TH1D* getHisto(TFile *f, TDirectoryFile* d, TString dName, TString hName)                                                            
{
  cout << "Dir Name : " << dName << endl; 
  cout << "Hist Name: " << hName << endl; 
  d = (TDirectoryFile*)f->Get(dName);
  if (d==NULL) 
  { 
    cout << endl; 
    printStuff("../ascii/cat.txt");
    fatal("Do not recognize TDir name. Check spellings?");
  }
    
 
  TH1D* h = (TH1D*)d->Get(hName);
  if (h==NULL) 
  { 
    cout << endl; 
    printStuff("../ascii/cat.txt");
    fatal("Do not recognize histogram. Check spellings?");
  }
  return h; 
}



////////////////////// draw text 
void drawText(double x, double y, TString txt, int col=kBlack, double size = 0.032, bool NDC = 1)
{
  static TLatex *tex = new TLatex(); 
  if (NDC) tex->SetNDC();
  tex->SetTextFont(42); 
  tex->SetTextSize(size); // 0.032 
  tex->SetTextColor(col);
  tex->DrawLatex(x,y,txt);
}

////////////////////// labels for histogram
void label(TString label)
{ 
  drawText(0.70, 0.85,"#bf{#it{ATLAS}} Internal");
  drawText(0.70, 0.80, label);
} 


////////////////////// make TGraph
TGraph *makeTGraph(vector<double> MassNum, vector<double> Vect, TString ylabel, int col=kBlack, TString xlabel = "m_{jj} [GeV]")
{ 
  TGraph * graph  = new TGraph(MassNum.size(),&MassNum[0],&Vect[0]);
  graph->SetMarkerStyle(20);
  graph->SetLineColor(col);
  graph -> SetMarkerSize(0.25);
  graph -> SetTitle("");

  // x-axis 
  graph -> GetXaxis() -> SetTitle(xlabel);
  graph -> GetXaxis() -> SetTitleFont(43);
  graph -> GetXaxis() -> SetTitleSize(14);
  graph -> GetXaxis() -> SetLabelFont(43);
  graph -> GetXaxis() -> SetLabelSize(14); 
  graph -> GetXaxis() -> SetMoreLogLabels(); 
  graph -> GetXaxis() -> SetTitleOffset(1.2);
 
  // y-axis
  graph -> GetYaxis() -> SetTitle(ylabel);
  graph -> GetYaxis() -> SetTitleFont(43);
  graph -> GetYaxis() -> SetTitleSize(14);
  graph -> GetYaxis() -> SetLabelFont(43);
  graph -> GetYaxis() -> SetLabelSize(14); 
  graph -> GetYaxis() -> SetTitleOffset(1.2);

  return graph; 
}

////////////////////// format TGraph
void formatTGraph(TGraph* graph, TString xlabel, TString ylabel, int col=kBlack)
{
  //graph->GetXaxis()->SetRangeUser(config.xminfit,config.xmaxfit);
  graph -> GetXaxis()->SetTitle(xlabel);
  graph -> GetXaxis() -> SetTitleSize(0.06);
  //graph -> GetXaxis() -> SetLabelSize(0.03);
  //graph -> GetXaxis() -> SetTitleOffset(0.9);
  
  graph -> GetYaxis()->SetTitle(ylabel);
  graph -> GetYaxis() -> SetTitleSize(0.06);
  //graph -> GetYaxis() -> SetLabelSize(0.03);
  graph -> GetYaxis() -> SetTitleOffset(0.7);
  
  graph->SetLineColor(col);
  graph -> SetMarkerStyle(20);
  graph -> SetMarkerSize(0.5);

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
    float Bin_xmin = hist -> GetBinLowEdge(lowBin); 
    float Bin_xmax = hist -> GetBinLowEdge(lowBin+1);
    float center   = hist -> GetBinCenter(lowBin);
   
    // for each bin, calculate the LLH
    double n           = hist->GetBinContent(lowBin); 
    double lambda      = func -> Integral(Bin_xmin, Bin_xmax);  
    if (lambda < 1e-10) lambda = 1e-10;
    double logLikeHood = n*log(lambda) - LnGamma(n+1) - lambda; 
   
    //if (lambda < 1e-10) lambda = 1e-10; 

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


///////////////////// create a TLine 
void makeTLine(double x1, double y1, double x2, double y2, int col = kBlack, TString drawOpt = " ", int lineWidth = 1, int lineStyle = 1, int trans = 0)
{ 
  TLine *line = new TLine(x1,y1,x2,y2);                                                                                       
  line -> SetLineColor(col); //Alpha(col, 0.5);
  if (trans) line -> SetLineColorAlpha(col, 0.3);
  line -> SetLineWidth(lineWidth);
  line -> SetLineStyle(lineStyle);
  line -> Draw(drawOpt); 

  lineY = line->GetY1(); 
  
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

////////////////////// check if parameter is fixed 
bool isParamFixed (TF1* function, int i)                                                                                                                                          
{
  double al, bl; 
  function -> GetParLimits(i,al,bl);
  return (al*bl != 0 && al >= bl); 
}




