////////////////////////////////////////////////////////////////////////////////                                                                      
// 
// Function to pick the window size 
// passes result to lowEdge and highEdge
//   
////////////////////////////////////////////////////////////////////////////////

void setWindowSize(double &lowEdge, double &highEdge, double windowCenter, vector<double> v_binLowEdge, double winSize, double winIdx)
{ 
  /*
  if (v_binLowEdge[winIdx]<=2215)
  {lowEdge=1200;
    highEdge=v_binLowEdge[winIdx+winSize];}
  else if (v_binLowEdge[winIdx]>=3827) 
  {lowEdge=v_binLowEdge[winIdx-winSize];
  highEdge=6165;}
 else
*/
  lowEdge=v_binLowEdge[winIdx-winSize];
  //highEdge=v_binLowEdge[winIdx+20];
  highEdge=v_binLowEdge[winIdx+winSize];
  /*
  double left  = v_binLowEdge[winIdx-winSize];
  double right=v_binLowEdge[winIdx+winSize];

  
  //==// calculate window sizes  
  double left  = windowCenter - (windowCenter * winSize/100);  
  double right = windowCenter + (windowCenter * winSize/100); 
  
  //==// find bin edge closest to window size calulated above 
  double diffLow  = 9999999.9;
  double diffHigh = 9999999.9; 
  for (int i = 0; i < v_binLowEdge.size(); i++)
  {  
    //==// low edge 
    double diffTempLow = abs(left - v_binLowEdge[i]); 
    if (diffTempLow < diffLow) 
    { 
      lowEdge = v_binLowEdge[i]; 
      diffLow = diffTempLow; 
    } 
      
    //==// high edge 
    double diffTempHigh = abs(right - v_binLowEdge[i]); 
    if (diffTempHigh < diffHigh) 
    {
      highEdge = v_binLowEdge[i]; 
      diffHigh = diffTempHigh; 
    }
  
  }
  */
  // if left edge is less than config.xminfit, set edge equal to config.xminfit 
  if (lowEdge < config.xminfit || lowEdge >windowCenter ) 
    lowEdge = config.xminfit;
      
  //==// if right edge is more than config.xmaxfit, set edge equal to config.xmaxfit 
  if (highEdge > config.xmaxfit || highEdge <windowCenter)
    highEdge = config.xmaxfit; 

/*

   if (lowEdge < 1200.0 || lowEdge >windowCenter )
    lowEdge = 1200.0;

  //==// if right edge is more than config.xmaxfit, set edge equal to config.xmaxfit 
  if (highEdge > 6165.0 || highEdge <windowCenter)
    highEdge =6165.0;
*/   

  //======================================================                                                                                        
  //  
  // Bunch of error checks
  //  
  //======================================================
  
  /*
  if (winSize == 0.0 )
  { 
    cout << endl;  
    printStuff("../ascii/kangaroo.txt");
    cout << "A window as a function of the signal width is asked for, however the width is set to zero in the config. Please pick a non-zero value. " << endl;
    cout << endl; 
    exit(0);
  }
  */


}

