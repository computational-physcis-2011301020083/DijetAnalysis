#ifndef __CONFIG
#define __CONFIG

using namespace std;

string nukeSpaces(string s)
{
  s.erase(remove(s.begin(),s.end(),' '),s.end());
  return s;
}

struct Config 
{
  bool getHistBin; 
  TString dataInput; 
  TString dataDir; 
  TString dataHist;  
  int    rebin;
  int    intBin; 
  TString sigFunc; 
  string sigName;
  TString resFunc; 
  double fitMin; 
  double fitMax; 
  double histMin; 
  double histMax; 
  double massMin; 
  double massMax; 
  double refit; 
  double delta; 
  bool goofyMode; 
  bool quietMode; 
  bool autoFudge;
  bool interpolate;   
  TString pathIn; 
  TString path; 
  TString pdf; 
  TString oFile;
  TString oRoot; 
  bool oneInFile; 
  TString fileName; 
  vector<bool> paramOpt; 
  vector<TString> FNAMES; 
  vector<TString> HNAMES; 
  vector<double> par; 
  vector<double> MASS; 
  vector<double> PAR1; 
  vector<double> PAR2; 
  vector<double> PAR3; 
  vector<double> PAR4; 
  vector<double> PAR5; 
  vector<double> PAR6; 
  vector<double> PAR7;                                                                                                                                      
  vector<double> PAR8; 
  vector<double> PAR9;                                                                                                                                      
  vector<double> PAR10; 

} config;

void pad(vector<double> &MassVec, vector<double> &vec, double val)
{ 
  while (vec.size() < MassVec.size() - 1)
    vec.push_back(0);
  vec.push_back(val);
}

bool check(string line, int num, string name)
{
  string sub = line.substr(0,num); 
  sub.erase(remove(sub.begin(), sub.end(), ' '), sub.end()); 
  return sub==name ; 
} 

void loadConfig(Config& config, TString cfgName) 
{
  ifstream fin(cfgName);

  if( fin.fail() )
  {   
    cout << endl; 
    cout << "Config file does not exist!" << endl; 
    cout << "Perhaps the path to the config is incorrect? : " << cfgName << endl; 
    cout << endl; 
    exit(0); 
  }

  string line;
  while (getline(fin, line)) 
  {
    int num = line.find("="); 
    int cmd = line.find(";"); 
    istringstream sin(line.substr(num + 1 , cmd - num -1 ));

    // ignore lines starting with #: comments 
    if (line.substr(0,1).find("#") != -1) continue; 
    
    if ( check(line, num, "rebin") ) 
      sin >> config.rebin;
    else if ( check(line, num, "getHistBin") ) 
        sin >> config.getHistBin; 
    else if ( check(line, num, "dataInput") )
        sin >> config.dataInput;
    else if ( check(line, num, "dataDir") )
        sin >> config.dataDir;
    else if ( check(line, num, "dataHist") )
        sin >> config.dataHist;
    else if ( check(line, num, "sigFunc") )
        sin >> config.sigFunc;
    else if ( check(line, num, "sigName") )
        sin >> config.sigName;
    else if ( check(line, num, "resFunc") )
        sin >> config.resFunc;
    else if ( check(line, num, "fitMin") )
        sin >> config.fitMin;
    else if ( check(line, num, "fitMax") )
        sin >> config.fitMax;
    else if ( check(line, num, "histMin") )
        sin >> config.histMin;
    else if ( check(line, num, "histMax") )
        sin >> config.histMax;
    else if ( check(line, num, "intBin") )
        sin >> config.intBin;
    else if ( check(line, num, "massMin") )
        sin >> config.massMin;
    else if ( check(line, num, "massMax") )
        sin >> config.massMax;
    else if ( check(line, num, "refit") )
        sin >> config.refit;
    else if ( check(line, num, "delta") )
        sin >> config.delta;
    else if ( check(line, num, "goofyMode") )
        sin >> config.goofyMode;
    else if ( check(line, num, "quietMode") )
        sin >> config.quietMode;
    else if ( check(line, num, "autoFudge") )
        sin >> config.autoFudge;
    else if ( check(line, num, "interpolate") )
        sin >> config.interpolate;
    else if ( check(line, num, "pdf") )
        sin >> config.pdf; 
    else if ( check(line, num, "oFile") )
        sin >> config.oFile; 
    else if ( check(line, num, "oRoot") )
        sin >> config.oRoot; 
    else if ( check(line, num, "pathIn") )
        sin >> config.pathIn;
    else if ( check(line, num, "path") )
        sin >> config.path;
    else if ( check(line, num, "sigFunc") )
        sin >> config.sigFunc;
    else if ( check(line, num, "fileName") )
        sin >> config.fileName;
    else if ( check(line, num, "oneInFile") )
        sin >> config.oneInFile;
    else if ( check(line, num, "paramOpt") )
        { bool f; sin >> f; config.paramOpt.push_back( f ); }        
    else if ( check(line, num, "FNAMES") )
        config.FNAMES.push_back( nukeSpaces(sin.str()) ); 
    else if ( check(line, num, "HNAMES") )
        config.HNAMES.push_back( nukeSpaces(sin.str()) ); 
    else if ( check(line, num, "par") )
        { double d; sin >> d; config.par.push_back( d ); }        
    else if ( check(line, num, "MASS") )
        { double d; sin >> d; config.MASS.push_back( d ); }        
    else if ( check(line, num, "PAR1") )  
        { double d; sin >> d; pad(config.MASS, config.PAR1, d); } 
    else if ( check(line, num, "PAR2") )
        { double d; sin >> d; pad(config.MASS, config.PAR2, d); } 
    else if ( check(line, num, "PAR3") )
        { double d; sin >> d; pad(config.MASS, config.PAR3, d); } 
    else if ( check(line, num, "PAR4") )
        { double d; sin >> d; pad(config.MASS, config.PAR4, d); } 
    else if ( check(line, num, "PAR5") )
        { double d; sin >> d; pad(config.MASS, config.PAR5, d); } 
    else if ( check(line, num, "PAR6") )  
        { double d; sin >> d; pad(config.MASS, config.PAR6, d); } 
    else if ( check(line, num, "PAR7") )  
        { double d; sin >> d; pad(config.MASS, config.PAR7, d); } 
    else if ( check(line, num, "PAR8") )  
        { double d; sin >> d; pad(config.MASS, config.PAR8, d); } 
    else if ( check(line, num, "PAR9") )  
        { double d; sin >> d; pad(config.MASS, config.PAR9, d); } 

 
  }
}

TH1D* makeVarBinHist( double bin )
{ 
  // TLA binning 
  double bins[] = { 362, 381, 400, 420, 441, 462, 484, 507, 531, 555, 580, 606, 633, 661, 690, 720, 
                    751, 784, 818, 853, 889, 927, 966, 1007, 1049, 1093, 1139, 1186, 1235, 1286, 1339, 
                    1394, 1451, 1511, 1573, 1637, 1704, 1773, 1845, 1920, 1998, 2079, 2163, 2251, 2342, 
                    2437, 2536, 2639, 2746, 2857, 2973, 3094, 3220, 3351, 3487, 3629, 3776, 3929, 4088, 
                    4254, 4427, 4607, 4794, 4989, 5191, 5402, 5621, 5849, 6086, 6333, 6590, 6857, 7135, 
                    7425, 7726, 8040, 8366, 8706, 9059, 9427, 9810, 10208, 10622, 11053, 11502, 11969, 
                    12455, 12960, 13486, 
                  };

  int binNum = sizeof(bins)/sizeof(double) -1 ;
  TH1D* hist = new TH1D( to_string((int)bin).c_str(),"", binNum, bins);
 
  return hist;  

}


#endif
