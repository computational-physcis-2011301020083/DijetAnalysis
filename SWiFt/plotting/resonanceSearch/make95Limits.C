#include "../../ATLASStyle/AtlasStyle.C"                                                                                                          
#include "../../ATLASStyle/AtlasUtils.C"
#include "../../code/helperFunctions/histPlotFormatFunctions.C"
#include "../../code/helperFunctions/statFitFunctions.C"                                                                                         

using namespace std; 
using namespace TMath;

TString DataPath  = "../../output/root/data/dijet_gausRes_LLHScan.root"; 
TString PEPath    = "../../condor/gausRes/root"; 
int numPE         = 100 ; 

// NOTE: to draw theory curve, uncomment lines 304 & 316
TString theoryTxt = "./theory/qStar.txt"; 
TString theory    = "q*"; 

double lumi       = 37.0   ;   
double minMass    = 1300.0 ;  
double maxMass    = 7100.0;
TString pdf = "./pdf/dijets_2015_2016_37fb_95CLs_gausRes.pdf"; 


//////////////////////////////////////////////
//
// Calculate meadian and uncertainty for expected limit 
//
//////////////////////////////////////////////
double findUncert( vector<double> vec_limit_PE, double value = 34.13)
{
  //==// get total number of PEs and find uncertainty bin index
  int numTotPE = vec_limit_PE.size(); 
  double uncert = 0.0;    
  int indexU = int(round( numTotPE * (value)/100  ));
 
  int index = numTotPE/2 + indexU; 
  if (numTotPE%2==0)                          
    uncert = (vec_limit_PE[index-1]+vec_limit_PE[index])/2; 
  else                                     
    uncert = vec_limit_PE[index];

  return uncert;   
}


//////////////////////////////////////////////
//
// Read text from csv file  
//
//////////////////////////////////////////////
int getInfo ( vector< vector<double> > &parameters, TString inputFile )
{
  // read entries 
  int num = 0; 
  string line; 
  ifstream file (inputFile);

  // skip file if it doesnt exist 
  bool toot = file.good();
  if (!toot) 
  {   
    cout << "File does not exist: " << inputFile << endl;  
    return true; 
  } 

  while ( getline (file,line) )
  {
    stringstream ss;                                                                                                                                    
    ss << line; 

    string each; 
    while ( getline(ss, each, ','))
      parameters[num].push_back(atof(each.c_str())); 
 
    num++; 
  }
  file.close();

  return false; 

}


///////////////////////////////////////////////
//
// Calculate 95% CLs and plot  
//
///////////////////////////////////////////////

void make95Limits()
{ 

  SetAtlasStyle ();

  //==// Make canvas 
  TCanvas *canvas = new TCanvas();
  canvas -> Update();
  canvas -> Print(pdf+"[");
  gPad   -> SetLeftMargin(0.12);
  gPad   -> SetBottomMargin(0.12);

  cout << "------------------------------------" << endl;
  cout << "Getting Observed limit from data " << endl; 
  cout << "------------------------------------" << endl; 

  //==// Read data file and get observed limit 
  TFile *dataFile = TFile::Open( DataPath, "r" ); 
  dataFile -> cd();
  TGraph *TG_obs_95CL          = (TGraph*)dataFile -> Get("Limit95");        // 95% CL on cross-section 
  formatTGraph(TG_obs_95CL, "m_{jj} [GeV]", "#sigma x A x BR [pb]", kBlack); 
  TG_obs_95CL                  -> SetLineWidth(2);
  int numPoints                 = TG_obs_95CL      -> GetN(); 


  cout << "------------------------------------" << endl;
  cout << "Getting Limits from pseudo-experiments" << endl; 
  cout << "------------------------------------" << endl; 

  //==// Now loop over all pseudo-experiments and get limits for each mass   
  vector<vector<double>> vec_limit_PE(numPoints); // numPoints x numPE vector 
  vector<double> MassNum; 

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
 
    //==// get limit for each mass 
    if ( !(PEFile->GetListOfKeys()->Contains("Limit95")) ) 
    {
      cout << "Skipping PE " << i+1 << endl; 
      continue;  
    }

    TGraph *TG_PE_95CL  = (TGraph*)PEFile -> Get("Limit95");
    double x, y;  
    for (int j = 0; j < numPoints; j++) 
    {  
      int k = TG_PE_95CL -> GetPoint(j, x, y);
      vec_limit_PE[j].push_back(y);  
      if (i ==0) MassNum.push_back(x);
    }   

    PEFile -> Close(); 
  }    
  cout << "Using " << vec_limit_PE[0].size() << "/" << numPE << " pseudo-experiments. " << endl; 
  cout << endl; 


  cout << "------------------------------------------" << endl;
  cout << "Calculating expected limit + uncertainties" << endl; 
  cout << "------------------------------------------" << endl; 

  //==// Calculate expected limit
  vector<double> exp_limit; 
  vector<double> exp_limit_1up; 
  vector<double> exp_limit_1down;
  vector<double> exp_limit_2up; 
  vector<double> exp_limit_2down;

  for (int i = 0; i < numPoints; i++ )
  { 
    //==// sort vector - increasing order
    sort(vec_limit_PE[i].begin(), vec_limit_PE[i].end() );

    //==// remove failed 95% fits  
    int numNegative = 0;  
    for (int j = 0; j < vec_limit_PE[i].size(); j++)
      if (vec_limit_PE[i][j] < 0) numNegative++; 
    vec_limit_PE[i].erase(vec_limit_PE[i].begin(), vec_limit_PE[i].begin() + numNegative);
    cout << "Mass " << MassNum[i] << " : " <<  vec_limit_PE[i].size() << "/" << numPE << endl;

    //==// calculate expected and uncertainies  
    exp_limit.push_back( findUncert(vec_limit_PE[i], 0.0) );
    exp_limit_1up.push_back( findUncert(vec_limit_PE[i], 34.13) ); 
    exp_limit_1down.push_back( findUncert(vec_limit_PE[i], -34.13) ); 
    exp_limit_2up.push_back( findUncert(vec_limit_PE[i], 47.72) ); 
    exp_limit_2down.push_back( findUncert(vec_limit_PE[i], -47.42) ); 

    //==// fill histograms 
    /*
    canvas -> SetLogy(0);  
    TH1D* hist = new TH1D("hist", "", 100, vec_limit_PE[i][0], vec_limit_PE[i][vec_limit_PE[i].size()-1]);
    for (int j = 0; j < vec_limit_PE[i].size(); j++) 
      hist -> Fill(vec_limit_PE[i][j]);
    hist -> Draw(); 
    makeTLine(exp_limit[i], 0, exp_limit[i], hist->GetMaximum(), kRed, "same", 2);
    makeTLine(exp_limit_1up[i], 0, exp_limit_1up[i], hist->GetMaximum(), kGreen+1, "same", 2);
    makeTLine(exp_limit_1down[i], 0, exp_limit_1down[i], hist->GetMaximum(), kGreen+1, "same", 2);
    makeTLine(exp_limit_2up[i], 0, exp_limit_2up[i], hist->GetMaximum(), kYellow, "same", 2);
    makeTLine(exp_limit_2down[i], 0, exp_limit_2down[i], hist->GetMaximum(), kYellow, "same", 2);
    canvas -> Print(pdf);
    delete hist;    
    */
            
  } 

  //==// Prepare vectors to store uncertainties - for drawing filled bands  
  vector<double> Mass2x; 
  vector<double> sig1; 
  vector<double> sig2; 
  for (int i = 0; i < numPoints; i++)
  {
    Mass2x.push_back(MassNum[i]);                                                                                                            
    sig1.push_back(exp_limit_1up[i]);
    sig2.push_back(exp_limit_2up[i]);
  }
  for (int i = numPoints-1; i >= 0; i--)
  {   
   Mass2x.push_back(MassNum[i]);
   sig1.push_back(exp_limit_1down[i]); 
   sig2.push_back(exp_limit_2down[i]); 
  } 


  cout << "------------------------------------" << endl;
  cout << "Get theory limits" << endl; 
  cout << "------------------------------------" << endl; 

  // count lines
  string l;
  int lines = 0 ;                                                                                                                 
  ifstream file (theoryTxt);
  while ( getline (file,l) ) lines++; 
  file.close();

  // read theory text file 
  vector< vector<double> > theoryPoints(lines);
  getInfo(theoryPoints, theoryTxt);  

  vector<double> mass_theory; 
  vector<double> xsec_theory;
  vector<double> accp_theory; 
  for (int i = 0; i < lines; i++)
  {   
    mass_theory.push_back(theoryPoints[i][1]); 
    accp_theory.push_back(theoryPoints[i][4]/theoryPoints[i][5]);
    xsec_theory.push_back((theoryPoints[i][2])*0.001*accp_theory[i]); 
  }
  

  cout << "------------------------------------" << endl;
  cout << "Make all the TGraphs" << endl; 
  cout << "------------------------------------" << endl; 

  //==// Make expected limit, theory TGraphs  
  TGraph* TG_exp_95CL = makeTGraph(MassNum, exp_limit, "#sigma x A x BR [pb]");
  TG_exp_95CL -> SetLineStyle(2);
  TG_exp_95CL -> SetLineWidth(2); 

  TGraph* TG_exp1sig_95CL = makeTGraph(Mass2x, sig1, "#sigma x A x BR [pb]");
  TG_exp1sig_95CL -> SetLineWidth(2); 
  TG_exp1sig_95CL -> SetFillColor(kGreen+1);

  TGraph* TG_exp2sig_95CL = makeTGraph(Mass2x, sig2, "#sigma x A x BR [pb]");
  TG_exp2sig_95CL -> SetLineWidth(2); 
  TG_exp2sig_95CL -> SetFillColor(kYellow);

  TGraph* TG_theory = makeTGraph(mass_theory, xsec_theory, " sigma x acc x BR [pb]"); 
  TG_theory -> SetLineStyle(9); 
  TG_theory -> SetLineWidth(2); 
  TG_theory -> SetLineColor(kBlue); 


  cout << "------------------------------------" << endl;
  cout << "Plot 95% CL plot" << endl; 
  cout << "------------------------------------" << endl; 

  //==// Draw limits plot
  canvas -> SetLogy(); 
  TG_exp2sig_95CL   -> Draw("AFL"); 
  TG_exp1sig_95CL   -> Draw("FL"); 
  TG_exp_95CL       -> Draw("L");  
  TG_obs_95CL       -> Draw("LP");
  //TG_theory         -> Draw("L");

  TG_exp2sig_95CL   -> GetXaxis() -> SetRangeUser(minMass, maxMass); 
  TG_exp2sig_95CL   -> SetMaximum(1); //5
  TG_exp2sig_95CL   -> SetMinimum(1e-5); //1e-4
 
  drawText(0.70, 0.85,"#bf{#it{ATLAS}} Internal", kBlack, 0.043);
  drawText(0.70, 0.79, "#sqrt{s}=13 TeV, "+roundUp(lumi,1)+"fb^{-1}", kBlack, 0.035);
  drawText(0.70, 0.74, "|y*| < 0.6 ", kBlack, 0.035);

  TLegend *leg = new TLegend(0.16,0.16,0.51,0.38);
  leg->SetBorderSize(0);
  //leg->AddEntry(TG_theory,theory,"l");
  leg->AddEntry(TG_obs_95CL,"Observed 95% CL upper limit","lp");
  leg->AddEntry(TG_exp_95CL,"Expected 95% CL upper limit","l");
  leg->AddEntry(TG_exp1sig_95CL,"Expected #pm 1#sigma","f");
  leg->AddEntry(TG_exp2sig_95CL,"Expected #pm 2#sigma","f");
  leg->Draw();

  canvas -> Print(pdf); 
  canvas -> Print(pdf+"]");   


}




