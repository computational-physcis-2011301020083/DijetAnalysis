
////////////////////////////////////////////////////////////////////////////////                                                                             
//
// All the fitting functions used by SWiFt
//   
////////////////////////////////////////////////////////////////////////////////

using namespace std; 
using namespace TMath;

double COMenergy = 13000.0; 


TString getFunction(TString label)
{
 
  TString callFunc = ""; 
  //==// bkg only functions 
  if      (label == "dijet3Param")  
    callFunc = "dijet3Param(x,[0],[1],[2])";
  else if (label == "dijet4Param")  
    callFunc = "dijet4Param(x,[0],[1],[2],[3])";
  else if (label == "dijet5Param")  
    callFunc = "dijet5Param(x,[0],[1],[2],[3],[4])";
  else if (label == "dijet6Param")  
    callFunc = "dijet6Param(x,[0],[1],[2],[3],[4],[5])";
  else if (label == "dijet7Param")  
    callFunc = "dijet7Param(x,[0],[1],[2],[3],[4],[5],[6])";
  //==// signal only functions 
  else if (label == "CBFunc")       
    callFunc = "CBFunc(x, [0],[1],[2],[3],[4])"; 
  else if (label == "GFunc")       
    callFunc = "GFunc(x, [0],[1],[2])"; 
  else if (label == "GrLFunc")  
    callFunc = "GrLFunc(x,[0],[1],[2],[3],[4],[5])"; 
  //==// S+B functions 
  else if (label == "G_dijet3Param") 
    callFunc = "G_dijet3Param(x,[0],[1],[2],[3],[4],[5])"; 
  else if (label == "G_dijet4Param") 
    callFunc = "G_dijet4Param(x,[0],[1],[2],[3],[4],[5],[6])"; 
  else if (label == "G_dijet5Param") 
    callFunc = "G_dijet5Param(x,[0],[1],[2],[3],[4],[5],[6],[7])"; 
  else if (label == "GrL_dijet3Param") 
    callFunc = "GrL_dijet3Param(x,[0],[1],[2],[3],[4],[5],[6],[7],[8])"; 
  else if (label == "GrL_dijet4Param") 
    callFunc = "GrL_dijet4Param(x,[0],[1],[2],[3],[4],[5],[6],[7],[8],[9])"; 
  else if (label == "GrL_dijet5Param") 
    callFunc = "GrL_dijet5Param(x,[0],[1],[2],[3],[4],[5],[6],[7],[8],[9],[10])"; 
  else if (label == "CB_dijet3Param") 
    callFunc = "CB_dijet3Param(x,[0],[1],[2],[3],[4],[5],[6],[7],[8])"; 
  else 
  { 
    cout << "!!! I do not recognize the function !!!" << endl; 
    exit(1); 
  }

  return callFunc; 

}

TString getSBFunction(TString bkgLabel, TString sigLabel)
{
 
  TString callFunc = ""; 
  if      (bkgLabel == "dijet3Param" && sigLabel == "GFunc")   
    callFunc = "G_dijet3Param";  
  else if (bkgLabel == "dijet4Param" && sigLabel == "GFunc")   
    callFunc = "G_dijet4Param";  
  else if (bkgLabel == "dijet5Param" && sigLabel == "GFunc")   
    callFunc = "G_dijet5Param";  
  else if (bkgLabel == "dijet3Param" && sigLabel == "GrLFunc")   
    callFunc = "GrL_dijet3Param";  
  else if (bkgLabel == "dijet4Param" && sigLabel == "GrLFunc")   
    callFunc = "GrL_dijet4Param";  
  else if (bkgLabel == "dijet5Param" && sigLabel == "GrLFunc")   
    callFunc = "GrL_dijet5Param";  
  else 
  { 
    cout << "Bkg or signal or S+B function definition not found. Check defineFunctions.C." << endl; 
    exit(0); 
  }    
 
  return callFunc;   

} 

////////////////////////////////////////////////
//
// Background only functions 
//
////////////////////////////////////////////////

double dijet3Param(double x=9., double par0=9., double par1=9., double par2=9.) 
{ 
  double mX   = x/COMenergy; 
  double func = par0 * pow((1-mX), par1) * pow(mX, par2);
  return func;  
}

double dijet4Param(double x=9., double par0=9., double par1=9., double par2=9., double par3=9.) 
{ 
  double mX   = x/COMenergy;
  double func = par0 * pow((1-mX), par1) * pow( mX, ( par2+ Log(pow(mX, par3)) ) );
  //double func = par0 * pow((1-mX), par1) * pow( mX, ( par2+ par3*mX ) ); //Func1 
  //double func = par0 * pow((1-mX), par1+ par3*mX) * pow( mX, ( par2) ) ; //Func2
  //double func = par0 * pow((1-mX), par1) * pow( mX,  par2 )+par3 ; //Func3
  //double func=par0*pow(10,par1*mX+par2*mX*mX)+par3; //Func4
  //double func=par0*pow(10,par1*mX+par2*mX*mX+par3*mX*mX*mX); //Func5
  //double func=pow(10,par0+par1*mX)+pow(10,par2+par3*mX); //Func6
  //double func=par0*Landau(-x,par1,par2,1)+par3; //Func7
  return func;  
}

double dijet5Param(double x=9., double par0=9., double par1=9., double par2=9., double par3=9., double par4=9.) 
{ 
  double mX   = x/COMenergy; 
  double func = par0 * pow((1-mX), par1) * pow( mX, ( par2 + par3*Log(mX) + par4*(Log(mX)*Log(mX))  ) );
  return func;
}

double dijet6Param(double x=9., double par0=9., double par1=9., double par2=9., double par3=9., double par4=9., double par5=9.) 
{ 
  double mX   = x/COMenergy; 
  double func = par0 * pow((1-mX), par1) * pow( mX, ( par2 + par3*Log(mX) + par4*(Log(mX)*Log(mX)) + par5*(Log(mX)*Log(mX)*Log(mX)) ) );
  return func;
}


double dijet7Param(double x=9., double par0=9., double par1=9., double par2=9., double par3=9., double par4=9., double par5=9., double par6=9.) 
{ 
  double mX   = x/COMenergy; 
  double func = par0 * pow((1-mX), par1) * pow( mX, ( par2 + par3*Log(mX) + par4*(Log(mX)*Log(mX)) + par5*(Log(mX)*Log(mX)*Log(mX)) + par6*(Log(mX)*Log(mX)*Log(mX)*Log(mX))) );
  return func;
}

////////////////////////////////////////////////
//
// Signal only functions 
//
////////////////////////////////////////////////

//==// gaussian 
double GFunc(double x=9., double par0=9., double par1=9., double par2=9.) 
{
  double func = par0*Gaus(x,par1,par2,1); 
  return func;
} 

//==// gaussian + reverse landau 
double GrLFunc(double x=9., double par0=9., double par1=9., double par2=9., double par3=9., double par4=9., double par5=9.)
{
  double func = par0*((par3*Gaus(x,par1,par2,1)+(1-par3)*Landau(-x,par4,par5,1))); 
  return func; 
}

//==// double sided crystal ball 
double CBFunc(double x=9., double par0=9., double par1=9., double par2=9., double par3=9., double par4=9., double par5=9.) 
{ 
  // constants 
  double dM_a = -0.0151; 
  double dM_b = -0.048; 
  double dM_c = 2.94e-4;   

  double sig_a = 1.699; 
  double sig_b = 0.644; 

  double alpLow_a = 1.477; 
  double alpLow_b = -7.15e-3; 
  double alpLow_c = -8.8e-5; 

  double alpHigh_a = 1.903; 
  double alpHigh_b = -2.33e-3; 
  double alpHigh_c = 9.23e-4; 

  double nLow  = 12.1; 
  double nHigh = 11.6; 
 
  // compute variables
  double mXn     = (x-100)/100; 
  double dM      = dM_a + (dM_b*mXn) + (dM_c*mXn*mXn);
  double sig     = sig_a + (sig_b*mXn); 
  double alpLow  = alpLow_a + (alpLow_b*mXn) + (alpLow_c*mXn*mXn);
  double alpHigh = alpHigh_a + (alpHigh_b*mXn) + (alpHigh_c*mXn*mXn);        
  double t       = (x-par1)/sig;  

  // compute CB shape
  double CB = 0.0; 
  if (t<-alpLow)
  { 
    double num  = exp(-0.5*alpLow*alpLow); 
    double r    = alpLow/nLow;   
    double den  = pow( (r*( (1/r) - alpLow - t )), nLow); 
    CB = num/den; 
  }
  if (t>alpHigh)
  { 
    double num  = exp(-0.5*alpHigh*alpHigh); 
    double r    = alpHigh/nHigh;   
    double den  =pow( (r*( (1/r) - alpHigh + t )), nHigh); 
    CB = num/den; 
  }  
  if (t >= -alpLow && t <= alpHigh)
    CB = exp(-0.5*t*t);

  return CB;  
  
}

////////////////////////////////////////////////
//
// S+B fit functions 
//
////////////////////////////////////////////////

//==// gaussian + 3 parameter dijet
double G_dijet3Param(double x=9., double par0=9., double par1=9., double par2=9., double par3=9., double par4=9., double par5=9. ) 
{
  double mX   = x/COMenergy;  
  double func = par0 * Gaus(x,par1,par2,1) + par3 * pow( (1-mX), par4 ) * pow( mX, par5 ); 
  return func;
}

//==// gaussian + 4 parameter dijet 
double G_dijet4Param(double x=9., double par0=9., double par1=9., double par2=9., double par3=9., double par4=9., double par5=9., double par6=9.) 
{
  double mX   = x/COMenergy;  
  double func = par0 * Gaus(x,par1,par2,1) + par3 * pow( (1-mX), par4 ) * pow( mX, ( par5 + par6*Log(mX) ) ); 
  return func;
}

//==// gaussian + 5 parameter dijet 
double G_dijet5Param(double x=9., double par0=9., double par1=9., double par2=9., double par3=9., double par4=9., double par5=9., double par6=9., double par7=9.) 
{
  double mX   = x/COMenergy;  
  double func = par0 * Gaus(x,par1,par2,1) + par3 * pow((1-mX), par4) * pow( mX, ( par5 + par6*Log(mX) + par7*(Log(mX)*Log(mX)) ) ); 
  return func;
}



//==// gaussian + reverse landau + 3 parameter dijet 
double GrL_dijet3Param(double x=9., double par0=9., double par1=9., double par2=9., double par3=9., double par4=9., double par5=9., 
                      double par6=9., double par7=9., double par8=9.)
{
  double mX   = x/COMenergy;
  double func = par0*((par3*Gaus(x,par1,par2,1)+(1-par3)*Landau(-x,par4,par5,1)))+par6 * pow((1-mX),par7) * pow(mX, par8);
  return func; 
}

//==// gaussian + reverse landau + 4 parameter dijet 
double GrL_dijet4Param(double x=9., double par0=9., double par1=9., double par2=9., double par3=9., double par4=9., double par5=9., 
                      double par6=9., double par7=9., double par8=9., double par9=9.)
{
  double mX   = x/COMenergy;
  double func = par0*((par3*Gaus(x,par1,par2,1)+(1-par3)*Landau(-x,par4,par5,1)))+ par6 * pow((1-mX), par7) * pow( mX, ( par8+ par9*Log(mX) ) ); 
  return func; 
}

//==// gaussian + reverse landau + 5 parameter dijet 
double GrL_dijet5Param(double x=9., double par0=9., double par1=9., double par2=9., double par3=9., double par4=9., double par5=9., 
                      double par6=9., double par7=9., double par8=9., double par9=9., double par10=9.)
{
  double mX   = x/COMenergy;
  double func = par0*((par3*Gaus(x,par1,par2,1)+(1-par3)*Landau(-x,par4,par5,1)))+ par6 * pow((1-mX), par7) * pow( mX, ( par8 + par9*Log(mX) + par10*(Log(mX)*Log(mX)))); 
  return func; 
}

//==// double sided crystal ball + 3 parameter dijet 
double CB_dijet3Param(double x=9., double par0=9., double par1=9., double par2=9., double par3=9., double par4=9., double par5=9., 
                      double par6=9., double par7=9., double par8=9.)
{ 
  double mX   = x/COMenergy;
  double func = par0*CBFunc(x,par1,par2,par3,par4,par5)+par6*pow((1-mX),par7)*pow(mX, par8);
  return func; 

}













