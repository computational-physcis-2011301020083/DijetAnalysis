
////////////////////////////////////////////////////////////////////////////////                                                                             
//
// All the fitting functions used by SWiFt
//   
////////////////////////////////////////////////////////////////////////////////

using namespace std; 
using namespace TMath;


TString getFunction(TString label)
{
 
  TString callFunc = ""; 
  if (label == "CBFunc")       
    callFunc = "CBFunc(x, [0],[1],[2],[3],[4])"; 
  else if (label == "GFunc")       
    callFunc = "GFunc(x, [0],[1],[2])"; 
  else if (label == "GrLFunc")  
    callFunc = "GrLFunc(x,[0],[1],[2],[3],[4],[5])"; 
  else if (label == "GrLEFunc")  
    callFunc = "GrLEFunc(x,[0],[1],[2],[3],[4],[5],[6],[7],[8])"; 
  else 
  { 
    cout << "!!! I do not recognize the function !!!" << endl; 
    exit(1); 
  }
  return callFunc; 

}


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

//==// gaussian + reverse landau 
double GrLEFunc(double x=9., double par0=9., double par1=9., double par2=9., double par3=9., double par4=9., double par5=9., double par6=9., double par7=9., double par8=9.)
{
  double func = par0*((par3*Gaus(x,par1,par2,1)+(1-par3)*Landau(-x,par4,par5,1))) + par6*(pow(x/13000,par7))*(pow((1-(x/13000)), par8)); 
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









