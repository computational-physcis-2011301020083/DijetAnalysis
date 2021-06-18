////////////////////////////////////////////////////////////////////////////////                                                                      
// 
// Generate signal for a mass by evaluating spline fits from morphing 
//   
////////////////////////////////////////////////////////////////////////////////

#include "histPlotFormatFunctions.C"

//* 
//* Get TSpline
//*
TSpline3* getTSpline(TFile *f, TString name)
{ 
  TSpline3 *s = (TSpline3*)f->Get(name);
  if (s==NULL)
  {   
    cout << endl; 
    printStuff("../ascii/cat.txt");
    fatal("Cannot access TSpline from file. Does it exist?");
  }   
  return s;                                                                                                                                           
} 

//* 
//* Get TGraph 
//*
TGraph* getTGraph(TFile *f, TString name)
{ 
  TGraph *g = (TGraph*)f->Get(name); 
  if (g==NULL)
  { 
    cout << endl; 
    printStuff("../ascii/cat.txt");
    fatal("Cannot access TGraph from file. Does it exist?");
  }
  return g;  
}

//*
//* Read in spline interpolations and return TF1 
//*
vector<double> generateSignal(double windowCenter, TFile* file, TF1 *SigFunc)
{
  int numPar = SigFunc->GetNpar(); 
  vector<double> SigParam;  
  TSpline3 *splines[numPar];
  TGraph   *tgraphs[numPar];  
  TString splineNames[] = {"s0", "s1", "s2", "s3", "s4", "s5"}; 
  TString tgraphNames[] = {"TG_par0", "TG_par1", "TG_par2", "TG_par3", "TG_par4", "TG_par5"}; 
 

  //==// get all the tsplines and tgraphs from input TFile 
  for (int i = 0; i< numPar; i++)
  { 
    splines[i] = getTSpline(file, splineNames[i]); 
    tgraphs[i] = getTGraph(file, tgraphNames[i]); 
  } 

  // evaluate splines at input mass and get signal shape 
  for (int i = 0; i< numPar; i++)
  { 
    SigParam.push_back(tgraphs[i]->Eval(windowCenter, splines[i]));
    SigFunc -> FixParameter( i, tgraphs[i]->Eval(windowCenter, splines[i]) );
  }

  return SigParam; 
} 

//*
//* Calculate variation between nominal and systematically varied shape
//*
void calculateDiff(double windowCenter, TF1* SigFunc, vector<TFile*> JES1upSigFile, vector<TFile*> JES1downSigFile, vector<vector<vector<double>>> &diffVec )
{ 

  //==// clone input TF1 so as to not modify it by mistake
  TF1* SigFuncCopy = (TF1*)SigFunc->Clone();

  //==// get nominal signal shape parameters from input TF1 
  vector<double> nomSigParam; 
  for (int i = 0; i< SigFuncCopy -> GetNpar(); i++) 
    nomSigParam.push_back(SigFuncCopy -> GetParameter(i)); 

  //==// get systematically varied shapes from input TFiles 
  vector< vector<double> > JES1upSigParam(JES1upSigFile.size()); 
  vector< vector<double> > JES1downSigParam(JES1downSigFile.size()); 
  
  for (int i = 0; i< JES1upSigFile.size(); i++)
  { 
    JES1upSigParam[i] = generateSignal(windowCenter, JES1upSigFile[i], SigFuncCopy);
    JES1downSigParam[i] = generateSignal(windowCenter, JES1downSigFile[i], SigFuncCopy);
  }  

  //==// calculate difference in parameters btw nominal and varied shapes and store results in vector 
  for (int j = 0; j < JES1upSigFile.size(); j++) // loop over 2nd index: number of JES kinds 
  { 
    for (int k = 0; k < nomSigParam.size(); k++) // loop over 3rd index: number of parameters 
    { 
      diffVec[0][j][k] = nomSigParam[k] - JES1upSigParam[j][k]; 
      diffVec[1][j][k] = nomSigParam[k] - JES1downSigParam[j][k];
    }
  }

  delete SigFuncCopy; 
}   

//*
//* Generate systematically varied shape 
//*
void generateSysSignal(vector<vector<vector<double>>> *diffVec, TF1* theFunc, TF1* theFuncRef, double *par)
{

  
  int sigmaNum = (*diffVec).size();  //2
  int JESNum   = (*diffVec)[0].size(); //3
  int paramNum = (*diffVec)[0][0].size(); //6
  int numTotParam = theFunc->GetNpar(); //9 
  vector<double> newParam(paramNum); 
  TF1* tempFunc    = (TF1*)theFuncRef->Clone();

  //==// loop over JES types and parameters and calculate new signal shape parameters 
  for (int i = 0; i < JESNum; i++)
  { 
    //cout << "JES" << i << " : " << par[numTotParam+1+i] << endl;
    //==// loop over all the parameters 
    for (int j = 0; j < paramNum; j++) 
    {
      
      //if (par[numTotParam+1+i]>0)  
      //  newParam[j] = tempFunc->GetParameter(j) - ( par[numTotParam+1+i] * (*diffVec)[0][i][j]) ;              
      //else  
      //  newParam[j] = tempFunc->GetParameter(j) + ( par[numTotParam+1+i] * (*diffVec)[1][i][j]) ; 
   
      //==// symmetrize variations
      double up = tempFunc->GetParameter(j) - ( par[numTotParam+1+i] * (*diffVec)[0][i][j]) ; 
      double dn = tempFunc->GetParameter(j) + ( par[numTotParam+1+i] * (*diffVec)[1][i][j]); 
      newParam[j] = (up + dn) / 2 ;      

      tempFunc -> FixParameter(j, newParam[j]);

    } 

  }
  //cout << endl; 

  //==// apply new signal shape to input TF1 
  for (int i = 1; i < paramNum; i++) // notice that we start with i = 0, since we only change the signal shape and not the normalization (i = 0)
  { 
    theFunc -> FixParameter(i, newParam[i]); 
    //par[i]  = newParam[i]; 
  } 
 
  delete tempFunc;  
}





























