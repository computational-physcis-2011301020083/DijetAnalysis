#!/usr/bin/env python
import os,math
from array import array
from ROOT import *
import argparse
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--path", dest='path',  default="", help="path")
parser.add_argument("--title", dest='title',  default="", help="title")
parser.add_argument("--outfile", dest='outfile',  default="", help="outfile")
parser.add_argument("--ftype", dest='ftype',  default="", help="ftype")
args = parser.parse_args()

def FitUn4p():
    #fbkg=TFile("/afs/cern.ch/user/d/ding/SWiFt/output/root/Btagged/Fulldata_mbj.root")
    fbkg=TFile("/afs/cern.ch/user/d/ding/SWiFt/output/root/Btagged/mbj_fulldataNew.root")
    hbkg=fbkg.Get("swiftBkgNominal")
    N=hbkg.GetNbinsX()
    #print N,hbkg.GetBinContent(1),hbkg.GetBinLowEdge(91)
        

    PE4dir="/afs/cern.ch/work/d/ding/RUN2/SWiFtFit/condor/"+args.path+"/root"
    PE4list=[]
    for f in os.listdir(PE4dir):     
      if os.path.getsize(PE4dir+"/"+f)>=1000:
        PE4list.append(PE4dir+"/"+f)
    print len(PE4list)
    
    PE5dir=PE4dir.replace("4p",args.ftype)
    PE5list=[]
    for f in os.listdir(PE5dir):
      if os.path.getsize(PE5dir+"/"+f)>=1000:
        PE5list.append(PE5dir+"/"+f)
    print len(PE5list)
   
   
    
    mean=[]    
    Uncern=[]
    for i in range(76):
      Uncern.append(0.0)
      mean.append(0.0)

    #print abs(-1)
    j=0.0
    for i in range(len(PE4list)):
     if PE4list[i].replace("4p",args.ftype) in PE5list: 
      fPE4=TFile(PE4list[i])
      hPE4=fPE4.Get("swiftBkgNominal")
      PE5=PE4list[i].replace("4p",args.ftype)
      fPE5=TFile(PE5)
      hPE5=fPE5.Get("swiftBkgNominal")
      j=j+1
      for i in range(1,77):
        Uncern[i-1]=Uncern[i-1]+hPE4.GetBinContent(i)-hPE5.GetBinContent(i)

      fPE4.Close()
      fPE5.Close()
    print j
    for i in range(76):
      Uncern[i]=abs(Uncern[i]/(j*hbkg.GetBinContent(i+1)))
    print Uncern
    fbkg.Close()
    
    bins=[1133, 1166, 1200, 1234, 1269, 1305, 1341, 1378, 1416, 1454, 1493, 1533, 1573, 1614, 1656, 1698, 1741, 1785, 1830, 1875, 1921, 1968, 2016, 2065, 2114, 2164, 2215, 2267, 2320, 2374, 2429, 2485, 2542, 2600, 2659, 2719, 2780, 2842, 2905, 2969, 3034, 3100, 3167, 3235, 3305, 3376, 3448, 3521, 3596, 3672, 3749, 3827, 3907, 3988, 4070, 4154, 4239, 4326, 4414, 4504, 4595, 4688, 4782, 4878, 4975, 5074, 5175, 5277, 5381, 5487, 5595, 5705, 5817, 5931, 6047, 6165]
    h2 = TH1D( args.title, args.title, len(bins)-1, array('f', bins) )
    for i in range(0,76):
      h2.SetBinContent(i+1,Uncern[i])
    outFile = TFile.Open(args.outfile,'UPDATE')
    outFile.cd()
    h2.Write(h2.GetName(),TObject.kOverwrite)
    outFile.Close()
    

   

   




if __name__=="__main__":
    FitUn4p()







