#ifndef __CONFIG
#define __CONFIG

using namespace std;

struct Config 
{

  //==// swift options 
  int doSWiFtScan; 

  //==// window details 
  int    rebin;
  double windowWidthLow; 
  double windowWidthHigh; 
  double windowNum; 
  double windowNumAlways; 
  double start; 
  double end; 
  double xminfit; 
  double xmaxfit;
  double tol; 

  //==// input, output files 
  TString inputFile; 
  TString inputDir; 
  TString inputHist; 
  TString nomSigFile; 
  vector<TString> JES1upSigFile; 
  vector<TString> JES1downSigFile; 
  TString path; 
  TString pdfName;
  TString oRootName; 
  double lumi;

  //==// fitting and plotting 
  bool do95Scan; 
  bool quietMinuit; 
  bool quietFitInfo; 
  bool quiet95FitInfo; 
  TString bkgFunc;
  TString bkgFunc2; 
  double norm; 
  vector<double> par;
  vector<double> par2; 
  TString SonlyFunc; 

  //==// signal subtraction 
  double sigSubLLHCut; 

  //==// systematics 
  bool addSys;
  int numFlatNP; 
  int numJESNP; 
  double scaleUn; 
  double lumiUn;  

  //==// minuit fitting parameters 
  TString integrator = "Gauss"; 
  double errorType   = 0.5; 
  int fitStrategy    = 2; 
  int maxIterations  = 500000; 
  double tolSimplex  = 0.05;
  double tolMigrad   = 0.05;  
 

} config; 


bool check(string line, int num, string name)
{
  string sub = line.substr(0,num); 
  sub.erase(remove(sub.begin(), sub.end(), ' '), sub.end()); 
  return sub==name ; 
}    

void loadConfig(Config& config, TString cfgName) 
{
  //ifstream fin(cfgName);

  istream *fin;
  if (cfgName == "-")
    fin = &cin;
  else
    fin = new ifstream(cfgName);

  if( fin->fail() )
  { 
    cout << endl; 
    cout << "Config file does not exist!" << endl; 
    cout << "Perhaps the path to the config is incorrect? : " << cfgName << endl; 
    cout << endl; 
    exit(0); 
  }

  string line;
  while (getline(*fin, line)) 
  {
    int num = line.find(":"); 
    int cmd = line.find(";"); 
    istringstream sin(line.substr(num + 1 , cmd - num -1 )); 

    // ignore lines starting with #: comments 
    if (line.substr(0,1).find("#") != -1) continue; 
    if ( check(line, num, "rebin") ) 
      sin >> config.rebin;
    else if ( check(line, num, "doSWiFtScan") ) 
      sin >> config.doSWiFtScan;
    else if ( check(line, num, "windowWidthLow") ) 
      sin >> config.windowWidthLow;
    else if ( check(line, num, "windowWidthHigh") ) 
      sin >> config.windowWidthHigh;
    else if ( check(line, num, "windowNum") ) 
      sin >> config.windowNum;
    else if ( check(line, num, "windowNumAlways") ) 
      sin >> config.windowNumAlways;
    else if ( check(line, num, "start") ) 
      sin >> config.start;
    else if ( check(line, num, "end") )
      sin >> config.end;
    else if ( check(line, num, "xminfit") )
      sin >> config.xminfit;
    else if ( check(line, num, "xmaxfit") )
      sin >> config.xmaxfit;
    else if ( check(line, num, "tol") )
      sin >> config.tol;
    else if ( check(line, num, "inputFile") )
      sin >> config.inputFile; 
    else if ( check(line, num, "inputDir") )
      sin >> config.inputDir; 
    else if ( check(line, num, "inputHist") )
      sin >> config.inputHist; 
    else if ( check(line, num, "nomSigFile") )
      sin >> config.nomSigFile; 
    else if (check(line, num, "JES1upSigFile"))  
      { TString str; sin >> str; config.JES1upSigFile.push_back( str ); } 
    else if (check(line, num, "JES1downSigFile"))  
      { TString str; sin >> str; config.JES1downSigFile.push_back( str ); } 
    else if ( check(line, num, "path") )
      sin >> config.path; 
    else if ( check(line, num, "pdfName") )
      sin >> config.pdfName; 
    else if ( check(line, num, "oRootName") )
      sin >> config.oRootName;
    else if ( check(line, num, "lumi") )
      sin >> config.lumi; 
    else if ( check(line, num, "do95Scan") )
      sin >> config.do95Scan; 
    else if ( check(line, num, "bkgFunc") )
      sin >> config.bkgFunc;
    else if ( check(line, num, "bkgFunc2") )
      sin >> config.bkgFunc2;
    else if ( check(line, num, "SonlyFunc") )
      sin >> config.SonlyFunc;
    else if ( check(line, num, "norm") )
      sin >> config.norm;
    else if ( check(line, num, "quietMinuit") )
      sin >> config.quietMinuit;
    else if ( check(line, num, "quietFitInfo") )
      sin >> config.quietFitInfo;
    else if ( check(line, num, "quiet95FitInfo") )
      sin >> config.quiet95FitInfo;
    else if (check(line, num, "par"))  
      { double d; sin >> d; config.par.push_back( d ); } 
    else if (check(line, num, "par2"))  
      { double d; sin >> d; config.par2.push_back( d ); } 
    else if ( check(line, num, "sigSubLLHCut") )
      sin >> config.sigSubLLHCut;

    //==// systematics 
    else if ( check(line, num, "addSys") )
      sin >> config.addSys; 
    else if ( check(line, num, "numFlatNP") )
      sin >> config.numFlatNP; 
    else if ( check(line, num, "numJESNP") )
      sin >> config.numJESNP; 
    else if ( check(line, num, "scaleUn") )
      sin >> config.scaleUn; 
    else if ( check(line, num, "lumiUn") )
      sin >> config.lumiUn; 
    
    //==// minuit integration settings
    else if ( check(line, num, "integrator") )
      sin >> config.integrator;
    else if ( check(line, num, "errorType") )
      sin >> config.errorType;
    else if ( check(line, num, "fitStrategy") )
      sin >> config.fitStrategy;
    else if ( check(line, num, "maxIterations") )
      sin >> config.maxIterations;
    else if ( check(line, num, "tolSimplex") )
      sin >> config.tolSimplex;
    else if ( check(line, num, "tolMigrad") )
      sin >> config.tolMigrad;
     
  }
}  





#endif
