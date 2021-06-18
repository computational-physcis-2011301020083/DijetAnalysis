////////////////////////////////////////////////////////////////////////////////                                                                       
// 
// Collection of statistical, goodness of fit, etc functions used by SWiFt
//   
////////////////////////////////////////////////////////////////////////////////

#ifndef __histPlotFormat
#define __histPlotFormat

using namespace std;
using namespace TMath;

#include "../config.C"

double lineY; 

//* 
//* Print some ascii animals
//*
void printStuff(TString file)
{ 
  string l; 
  ifstream f(file); 
  while ( getline (f,l) ) 
    cout << l << endl;  
  f.close();
}
 
//*
//* Function to write text onto canvas                                                             
//*
void fatal(TString msg) { printf("\nFATAL\n  %s\n\n",msg.Data()); abort(); }

//*
//* Open file 
//*
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

//*
//* Get histo from TFile
//*
TH1D* getHisto(TFile *f, TString hName)                                                            
{
  cout << "Hist Name: " << hName << endl;
  TH1D *h = (TH1D*)f->Get(hName);
  if (h==NULL) 
  {
    cout << endl; 
    printStuff("../ascii/cat.txt"); 
    fatal("Cannot access TH1D from file. Does it exist?");
  }
  return h; 
}

//*
//* Get histo from TDir 
//*
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

//*
//* Draw text 
//*
void drawText(double x, double y, TString txt, int col=kBlack, double size = 0.032, bool NDC = 1)
{
  static TLatex *tex = new TLatex(); 
  if (NDC) tex->SetNDC();
  tex->SetTextFont(42); 
  tex->SetTextSize(size); // 0.032 
  tex->SetTextColor(col);
  tex->DrawLatex(x,y,txt);
}

//*
//* Format TF1 
//*
void formatFunc(TF1 *func, int col = kBlack, int size = 1, int style = 1)
{ 
  func -> SetLineColor(col); 
  func -> SetLineWidth(size);
  func -> SetLineStyle(style);
}

//*
//* Format histogram     
//*                    
void formatHist(TH1D *h, TString ylbl = "", TString xlbl = "" )
{
  h -> SetLineWidth(1);
  h -> SetTitle(""); 
  h -> SetMarkerSize(0.5);

  // x-axis 
  h -> GetXaxis() -> SetTitle(xlbl);
  h -> GetXaxis() -> SetTitleFont(43);
  h -> GetXaxis() -> SetTitleSize(14);
  h -> GetXaxis() -> SetLabelFont(43);
  h -> GetXaxis() -> SetLabelSize(14); 
  h -> GetXaxis() -> SetMoreLogLabels(); 
  //h -> GetXaxis() -> SetTitleOffset(0.7);
 
  // y-axis
  h -> GetYaxis() -> SetTitle(ylbl);
  h -> GetYaxis() -> SetTitleFont(43);
  h -> GetYaxis() -> SetTitleSize(14);
  h -> GetYaxis() -> SetLabelFont(43);
  h -> GetYaxis() -> SetLabelSize(14); 
  //h -> GetYaxis() -> SetTitleOffset(0.7);

}

//*
//* Format diff plot
//*
void formatDiff(TH1D *diff, TString ylbl, TString xlbl = "m_{jj} [GeV]", bool setFill = 0, int fillCol = kBlue, bool xOff = 0) 
{ 
  
  if (setFill) diff  -> SetFillColor(fillCol); 

  diff -> GetYaxis() -> SetTitle(ylbl);
  diff -> GetYaxis() -> CenterTitle();
  diff -> GetYaxis() -> SetTitleFont(43);                                                                                                                      
  diff -> GetYaxis() -> SetTitleSize(10);
  diff -> GetYaxis() -> SetLabelFont(43);
  diff -> GetYaxis() -> SetLabelSize(10); 
  //diff -> GetYaxis() -> SetTitleOffset(1.8);  
  diff -> GetYaxis() -> SetNdivisions(407); //507

  diff -> GetXaxis() -> SetTitle(xlbl); 
  diff -> GetXaxis() -> SetTitleFont(43);
  diff -> GetXaxis() -> SetTitleSize(14);
  diff -> GetXaxis() -> SetLabelFont(43);
  diff -> GetXaxis() -> SetLabelSize(14);
  diff -> GetXaxis() -> SetTitle(xlbl); 
  diff -> GetXaxis() -> SetTitleOffset(6);

  diff -> SetStats(0);     
  
  if (xOff)
  { 
    diff->GetXaxis()->SetLabelOffset(999);
    diff->GetXaxis()->SetLabelSize(0);
    diff->GetXaxis()->SetTitle(" "); 
  }  

}

//*
//* Labels for histogram
//*
void label(TString label)
{ 
  drawText(0.70, 0.85,"#bf{#it{ATLAS}} Internal");
  drawText(0.70, 0.80, label);
} 

//*
//* Get histogram from file
//*
vector<TH1*> getHist(TFile *file, vector<TH1*> histVec, TString plotName)
{ 
  TIter next(file->GetListOfKeys());
  TH1 *histpart;
  TKey *key;
  while ((key=(TKey*)next()))
  {
    TString name = key->GetName();
    histpart     = (TH1*) file -> Get(name);
    if(name.Contains(plotName)) histVec.push_back(histpart);
  }

  return histVec;
}

//*
//* Make TGraph
//*
TGraph *makeTGraph(vector<double> VecX, vector<double> VecY, TString ylabel, int col=kBlack, TString xlabel = "m_{jj} [GeV]")
{ 
  TGraph * graph  = new TGraph(VecX.size(),&VecX[0],&VecY[0]);
  graph->SetMarkerStyle(20);
  graph->SetLineColor(col);
  graph -> SetMarkerSize(0.50); //0.25

  // x-axis 
  graph -> GetXaxis() -> SetTitle(xlabel);
  graph -> GetXaxis() -> SetTitleFont(43);
  graph -> GetXaxis() -> SetTitleSize(14);
  graph -> GetXaxis() -> SetLabelFont(43);
  graph -> GetXaxis() -> SetLabelSize(14); 
  graph -> GetXaxis() -> SetMoreLogLabels(); 
  //h -> GetXaxis() -> SetTitleOffset(0.7);
 
  // y-axis
  graph -> GetYaxis() -> SetTitle(ylabel);
  graph -> GetYaxis() -> SetTitleFont(43);
  graph -> GetYaxis() -> SetTitleSize(14);
  graph -> GetYaxis() -> SetLabelFont(43);
  graph -> GetYaxis() -> SetLabelSize(14); 
  //h -> GetYaxis() -> SetTitleOffset(0.7);

  return graph; 
}


//*
//* Make TGraph with bands 
//*
TGraph *makeBandTGraph(vector<double> VecX, vector<double> sigUp, vector<double> sigDown, int col = kGreen+1, 
                       TString ylabel = "JES", TString xlabel = "m_{jj} [GeV]")                         
{ 

  vector<double> Vec2X; 
  vector<double> Vec2Xsig; 
  for (int i = 0; i < VecX.size(); i++)
  {
    Vec2X.push_back(VecX[i]);
    Vec2Xsig.push_back(sigUp[i]);
  }
  for (int i = VecX.size()-1; i >= 0; i--)
  {   
    Vec2X.push_back(VecX[i]);
    Vec2Xsig.push_back(sigDown[i]); 
  }  

  TGraph* sigBand = makeTGraph(Vec2X, Vec2Xsig, ylabel, col, xlabel );
  sigBand -> SetFillColor(col);

  return sigBand; 

}

//*
//* Make TGraphErros
//*
TGraphErrors *makeTGraphErrors(vector<double> VecX, vector<double> VecY, vector<double> VecErrX, vector<double> VecErrY, TString ylabel, 
                         int col=kBlack, TString xlabel = "m_{jj} [GeV]")
{ 

  TGraphErrors* graph = new TGraphErrors(VecX.size(), &VecX[0], &VecY[0] , &VecErrX[0] , &VecErrY[0]);
  graph->SetMarkerStyle(20);
  graph->SetLineColor(col);
  graph -> SetMarkerSize(0.70); //0.25

  // x-axis
  graph -> GetXaxis() -> SetTitle(xlabel);
  graph -> GetXaxis() -> SetTitleFont(43);
  graph -> GetXaxis() -> SetTitleSize(14);
  graph -> GetXaxis() -> SetLabelFont(43);
  graph -> GetXaxis() -> SetLabelSize(14); 
  graph -> GetXaxis() -> SetMoreLogLabels(); 
  //h -> GetXaxis() -> SetTitleOffset(0.7);
 
  // y-axis
  graph -> GetYaxis() -> SetRangeUser(-2.5, 2.5); 
  graph -> GetYaxis() -> SetTitle(ylabel);
  graph -> GetYaxis() -> SetTitleFont(43);
  graph -> GetYaxis() -> SetTitleSize(14);
  graph -> GetYaxis() -> SetLabelFont(43);
  graph -> GetYaxis() -> SetLabelSize(14); 
  //h -> GetYaxis() -> SetTitleOffset(0.7);

  return graph; 
}

//*
//* Format TGraph
//*
void formatTGraph(TGraph* graph, TString xlabel, TString ylabel, int col=kBlack)
{
  graph -> GetXaxis() -> SetTitle(xlabel);
  graph -> GetXaxis() -> SetTitleFont(43);
  graph -> GetXaxis() -> SetTitleSize(14);
  graph -> GetXaxis() -> SetLabelFont(43);
  graph -> GetXaxis() -> SetLabelSize(14); 
  graph -> GetXaxis() -> SetMoreLogLabels(); 
  //h -> GetXaxis() -> SetTitleOffset(0.7);
 
  // y-axis
  graph -> GetYaxis() -> SetTitle(ylabel);
  graph -> GetYaxis() -> SetTitleFont(43);
  graph -> GetYaxis() -> SetTitleSize(14);
  graph -> GetYaxis() -> SetLabelFont(43);
  graph -> GetYaxis() -> SetLabelSize(14); 
  //h -> GetYaxis() -> SetTitleOffset(0.7);

  graph->SetLineColor(col);
  graph -> SetMarkerStyle(20);
  graph -> SetMarkerSize(0.5);

}

//*
//* Set best TGraph range 
//*
void setBestRangeTGraph (TGraph* graph1, TGraph* graph2)  
{ 

  double graph1Min = TMath::MinElement(graph1->GetN(), graph1->GetY()); 
  double graph1Max = TMath::MaxElement(graph1->GetN(), graph1->GetY()); 
  double graph2Min = TMath::MinElement(graph2->GetN(), graph2->GetY()); 
  double graph2Max = TMath::MaxElement(graph2->GetN(), graph2->GetY()); 

  //==// set low range  
  if (graph1Min < graph2Min)
  {
    graph1 -> SetMinimum(graph1Min);   
    graph2 -> SetMinimum(graph1Min);   
  }
  else
  {
    graph1 -> SetMinimum(graph2Min);   
    graph2 -> SetMinimum(graph2Min);   
  }
    
  //==// set high range  
  if (graph1Max > graph2Max)
  {
    graph1 -> SetMaximum(graph1Max);   
    graph2 -> SetMaximum(graph1Max);   
  }
  else
  {
    graph1 -> SetMaximum(graph2Max);   
    graph2 -> SetMaximum(graph2Max);   
  }
  
}

//*
//* Set Range for TH1
//*
void setBestRangeTH1 (TH1D* hist, double threshold)
{ 
  double max = hist -> GetMaximum(); 
  double min = hist -> GetMinimum(); 

  if (max < threshold || min > -threshold)
  { 
    hist -> SetMaximum( threshold ); 
    hist -> SetMinimum( -threshold ); 
  }

} 

//* 
//* Divide histogram by bin width 
//*
void fixBinWidth(TH1D *hist)
{ 
  for (int i = 1; i <= hist->GetNbinsX(); i++)
  { 
    hist->SetBinContent(i, hist->GetBinContent(i)/hist->GetBinWidth(i) );
    hist->SetBinError(i, hist->GetBinError(i)/hist->GetBinWidth(i) );
  } 
  
} 

//*
//* Convert fit to histogram 
//*
void fit_to_hist(TH1D* localHist, TF1* BfuncBkg, double min, double max)
{
  int bMin = localHist -> GetXaxis()->FindBin(min);
  int bMax = localHist -> GetXaxis()->FindBin(max);
  localHist -> SetLineColor(BfuncBkg->GetLineColor());
 
  for (int i = 1; i <= localHist->GetNbinsX(); i++)
  {
    if (i < bMin) 
    { 
      localHist -> SetBinContent(i, 0); // set bin content to 0 if bin < min
      localHist -> SetBinError(i, 0); 
    }
    else if (i > bMax) 
    {
      localHist -> SetBinContent(i, 0); // set bin content to 0 if bin > max 
      localHist -> SetBinError(i, 0);
    }
    else 
    {
      localHist -> SetBinContent(i, (BfuncBkg->Integral(localHist->GetBinLowEdge(i), localHist->GetBinLowEdge(i+1))) );
      localHist -> SetBinError(i, 0);
    }
  }        
}

//*
//* Initialize TF1 for bkg-only fits 
//*
void armFunc(TF1 *func, vector<double> paramVec)
{ 
  // initialize parameters 
  int numParam = paramVec.size(); 
  for (int i = 0; i < numParam; i++)  
    func -> SetParameter(i, paramVec[i]);
}

//*
//* Initialize TF1 for S+B fits 
//*
void armFunc(TF1 *func, TF1 *funcSig, TF1 *funcBkg )
{ 

  int numSigParam = funcSig -> GetNpar(); 
  int numBkgParam = funcBkg -> GetNpar(); 

  // initialize parameters for signal part 
  for (int i = 0; i < numSigParam; i++)
    if (i == 0)
      func -> SetParameter(i, config.norm); 
    else 
      func -> FixParameter(i, funcSig->GetParameter(i)); 

  // initialize parameters for bkg part 
  for (int i = 0; i < numBkgParam; i++)
    func -> SetParameter(numSigParam+i, funcBkg->GetParameter(i));
}

//*
//* Initialize TF1 for bkg component of S+B fits 
//*
void armFunc(TF1 *func, TF1* inFunc, int offset)
{ 
  int numParam = func->GetNpar(); 

  for (int i = 0; i < numParam; i++)                                                                                                                  
    func -> FixParameter(i, inFunc->GetParameter(offset+i));  
}

//*
//* Divide by uncertainty 
//*
void divideUncert(TH1D* diff, TH1D* data ) 
{ 
  for (int i = 1; i <= diff->GetNbinsX(); i++) 
  { 
    if (data->GetBinContent(i) > 0) diff -> SetBinContent(i, diff->GetBinContent(i)/sqrt(data->GetBinContent(i)) );     
    else diff -> SetBinContent(i, 0 );   
   }

} 

//*
//* Create a TLine 
//*
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

//*
//* Create a TText
//*
void makeTText(double x, double y, TString label = "", int col = kBlack, int trans = 0,  double size = 0.032)                                  
{
  TText *t = new TText(x,y,label);
  t -> SetTextColor(col);
  if (trans) t -> SetTextColorAlpha(col, 0.3);
  t -> SetTextAlign(12);
  t -> SetTextSize(size);
  t -> SetTextFont(122); 
  t -> Draw(); 

}

//*
//* Round off doubles 
//*
TString roundUp(double n, int precision = 1, bool fix =1)
{
  stringstream stream;
  if (fix) stream << fixed << setprecision(precision) << n;
  else 
  { 
    stream << scientific << setprecision(precision) << n;
  }
  TString s = stream.str();
  return s; 
  
} 

//*
//* Save info of fits that fail 
//*
void checkForFails( vector<double> tempfitQuality, vector<vector<double>> *fitQuality, int currentBin, int startBin, double windowCenter  )
{ 

  if(tempfitQuality[0] != 0 || tempfitQuality[1] != 0 || tempfitQuality[2] != 0) // if SIMPLEX or MIGRAD or HESS fail 
  { 
    // push back window center
    (*fitQuality)[currentBin-startBin].push_back(windowCenter); 

    // push back SIMPLEX and MIGRAD and HESS fit status 
    for (int j = 0; j < tempfitQuality.size(); j++)
      (*fitQuality)[currentBin-startBin].push_back(tempfitQuality[j]);
  }

}

//*
//* Make data & full background plots 
//*
void plotBkg (TCanvas* canvas, TString pdf, TH1D* histRaw, TH1D* bkgNomHist, TString legLabel1 = "Data", TString legLabel2 = "SWiFt Bkg", double chisqpVal = 999.9 ) 
{ 

  //==// make clones so that the originals are not modified 
  TH1D* hist    = (TH1D*) histRaw    -> Clone("hist"); 
  TH1D* bkgHist = (TH1D*) bkgNomHist -> Clone("bkgHist");

  //==// some formatting
  hist    -> SetMarkerSize(0.5);
  bkgHist -> SetLineWidth(2); 

  //==// set up pads 
  TPad *pad1 = new TPad("pad1","pad1 ", 0, 0.30, 1, 1);
  pad1 -> SetBottomMargin(0);
  pad1 -> SetLeftMargin(0.12); 
  pad1 -> SetLogy(1);
  //pad1 -> SetLogx(1); //
  pad1 -> Draw();

  TPad *pad2 = new TPad("pad2","pad2", 0, 0.01, 1, 0.30);
  pad2->SetBottomMargin(.4);
  pad2 -> SetLeftMargin(0.12); 
  pad2 -> SetLogy(0);
  //pad2 -> SetLogx(1); //
  pad2->Draw();
  
  //==// draw on pad1 
  pad1    -> cd(); 
  formatHist(hist, "Events");
  bkgHist -> Draw("][ hist");
  hist    -> Draw("e1 same");
  drawText( 0.75, 0.85, "#bf{#it{ATLAS}} Internal", kBlack, 0.050);
  drawText( 0.75, 0.80, legLabel1, hist->GetLineColor(), 0.040); 
  drawText( 0.75, 0.75, legLabel2, bkgHist->GetLineColor(), 0.040);
  drawText( 0.75, 0.70, "#chi^{2} pvalue: "+to_string(chisqpVal), kBlack, 0.040); 

  //==// draw on pad2 
  pad2 -> cd(); 
  TH1D *histDiff = (TH1D*)hist -> Clone("histDiff"); 
  histDiff -> SetLineWidth(2);
  histDiff -> Add(bkgHist, -1); 
  divideUncert(histDiff, hist ); 
  formatDiff(histDiff, "Significance", "m_{jj} [GeV]", 1, bkgHist -> GetLineColor());  
  setBestRangeTH1(histDiff, 1.0); 
  histDiff -> GetYaxis() -> SetTitleSize(12);
  histDiff -> GetYaxis() -> SetLabelSize(12);
  histDiff -> GetXaxis() -> SetMoreLogLabels();  
  histDiff -> Draw("][ hist"); 
  makeTLine(config.xminfit, 0, config.xmaxfit, 0, kBlack, "hist same", 1);  

  //==// print to canvas 
  canvas -> Print(pdf); 
  canvas -> Clear(); 

  //==// delete clones 
  delete hist; 
  delete bkgHist; 
  delete histDiff; 

}

void plotWindowFits(double windowCenter, double lowEdge, double highEdge, TH1D* histRaw, TH1D* BHistBkg, TCanvas *canvas, TString pdf, 
                    double bestpVal, double bestChisq, int windowNDF)
{
  //==// clear canvas 
  canvas -> Clear(); 

  //==// make pads for plots per window                                                                                    
  TPad *pad1 = new TPad("pad1","pad1", 0 , 0.40 , 1 , 1);
  pad1 -> SetBottomMargin(0);
  pad1 -> SetLeftMargin(0.12);  
  pad1 -> SetLogy(1);
  pad1 -> Draw(); 

  TPad *pad2 = new TPad("pad2","pad2", 0 , 0.16 , 1 , 0.40);
  pad2 -> SetBottomMargin(.4);
  pad2 -> SetLeftMargin(0.12);
  pad2 -> SetGridx(1);
  pad2 -> SetGridy(1);
  pad2 -> Draw();
  
  TPad *pad3 = new TPad("pad3","pad3", 0 , 0.03 , 1 , 0.25);
  pad3 -> SetBottomMargin(.4);
  pad3 -> SetLeftMargin(0.12);
  pad3 -> SetGridx(1);
  pad3 -> SetGridy(1);
  pad3 -> Draw();

  //==// set range on clone of input histogram 
  TH1D* histClone = (TH1D*) histRaw -> Clone("histClone");
  histClone -> GetXaxis() -> SetRangeUser(lowEdge, highEdge);

  //==// pad1 
  pad1 -> cd(); 
  canvas    -> SetLogy(); 
  histClone -> Draw("e1");                                                                                                         
  BHistBkg  -> Draw("][ hist same"); 
  histClone -> SetMinimum(BHistBkg->GetBinContent(BHistBkg->FindBin(highEdge)));
  double ymax = BHistBkg -> GetBinContent(BHistBkg->FindBin(windowCenter)); 
  makeTLine(windowCenter, 0, windowCenter, ymax, kGreen-10);
  drawText(0.68, 0.85,"#bf{#it{ATLAS}} Internal", kBlack, 0.07);
  drawText(0.68, 0.78, "Mass: "+roundUp(windowCenter,1)+" GeV", kBlack, 0.05 );                                                             
  drawText(0.68, 0.72, "Window: ["+roundUp(lowEdge,0)+", "+roundUp(highEdge,0)+"] GeV", kBlack, 0.05 );
  drawText(0.68, 0.66, "#chi^{2} p-value: "+roundUp(bestpVal, 5, 0), kBlack, 0.05 ); 
  drawText(0.68, 0.60, "#chi^{2} / NDF: "+roundUp(bestChisq)+" / "+to_string(windowNDF), kBlack, 0.05 ); 
  canvas    -> SetLogy(0); 

  //==// pad2
  pad2 -> cd();
  TH1D *diffBkgOnly = (TH1D*) histClone -> Clone("diffBkgOnly");
  diffBkgOnly -> Add(BHistBkg , -1);
  diffBkgOnly -> Draw("hist E1 P");
  diffBkgOnly -> SetMarkerSize(0.8);
  formatDiff(diffBkgOnly, "Data-Bkg_{only}", "m_{jj} [GeV]", 0, 1, 1);
  makeTLine(lowEdge, 0, highEdge, 0, kBlack, "same"); 

  //==// pad3
  pad3 -> cd();
  TH1D *diffBkgOnly2 = (TH1D*) histClone -> Clone("diffBfromSBUn");
  diffBkgOnly2 -> Add(BHistBkg , -1);
  divideUncert(diffBkgOnly2, histClone); 
  diffBkgOnly2 -> Draw("hist P");
  diffBkgOnly2 -> SetMarkerSize(0.8);
  formatDiff(diffBkgOnly2, "#frac{Data-Bkg_{only}}{#sqrt{Data}} ", "m_{jj} [GeV]", 0, 1, 0);
  makeTLine(lowEdge, 0, highEdge, 0, kBlack, "same"); 
  
  //==// print full canvas 
  canvas -> Print(pdf);

  //==// clear tmp allocated memory 
  delete histClone, diffBkgOnly, diffBkgOnly2; 
  delete pad1, pad2, pad3; 

}
















#endif



