////////////////////////////////////////////////////////////////////////////////                                                                      
// 
// Function to make SWiFt bkg                                                                                                  
//   
////////////////////////////////////////////////////////////////////////////////

#include "histPlotFormatFunctions.C"

void makeBkg(int winIdx, TH1D* swiftBkgHist, TH1D* bkgHistLocal )
{ 

  int start_bin    = swiftBkgHist -> FindBin(config.start);   // SWiFt starting point                                                                                        
  int end_bin      = swiftBkgHist -> FindBin(config.end);     // SWiFt ending point 
  int xminfit_bin  = swiftBkgHist -> FindBin(config.xminfit); // data hist min mass 
  int xmaxfit_bin  = swiftBkgHist -> FindBin(config.xmaxfit); // data hist max mass

  //==// set bin content as zero below xminfit and above xmaxfit 
  for (int i = 0; i < swiftBkgHist->GetNbinsX(); i++)
  {
    if (i < xminfit_bin || i > xmaxfit_bin) 
    { 
      swiftBkgHist -> SetBinContent(i, 0.0); 
      swiftBkgHist -> SetBinError(i, 0.0);  
    }
  }  

  //==// populate swiftBkgHist with swift bkg
  if ( winIdx == start_bin ) 
  { 
    for (int i = xminfit_bin; i <= start_bin ; i++) 
    {
      swiftBkgHist -> SetBinContent( i, bkgHistLocal -> GetBinContent(i) ); 
      swiftBkgHist -> SetBinError(i, sqrt(swiftBkgHist -> GetBinContent(i))); 
    }
  }  
  
  else if ( winIdx == end_bin )
  {
    for (int i = end_bin; i <= xmaxfit_bin ; i++) 
    {
      swiftBkgHist -> SetBinContent( i, bkgHistLocal -> GetBinContent(i) );          
      swiftBkgHist -> SetBinError(i, sqrt(swiftBkgHist -> GetBinContent(i)));
    }
  } 
 
  else
  {
    swiftBkgHist -> SetBinContent( winIdx, bkgHistLocal -> GetBinContent(winIdx) ); 
    swiftBkgHist -> SetBinError(winIdx, sqrt(swiftBkgHist -> GetBinContent(winIdx)));
  }
}





