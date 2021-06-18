#!/usr/bin/env python
import os,math
from array import array
from ROOT import *
import argparse
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--path", dest='path',  default="", help="path")
parser.add_argument("--title", dest='title',  default="", help="title")
parser.add_argument("--outfile", dest='outfile',  default="", help="outfile")
args = parser.parse_args()
def FitUn4p():
    fbkg=TFile("/afs/cern.ch/user/d/ding/SWiFt/output/root/Btagged/mbj_fulldataNew.root")
    hbkg=fbkg.Get("swiftBkgNominal")
    N=hbkg.GetNbinsX()
    #for i in range(N):
      #print i,hbkg.GetBinContent(i),hbkg.GetBinLowEdge(i)
        
    print args.path
    PEdir="/afs/cern.ch/work/d/ding/RUN2/SWiFtFit/condor/"+args.path+"/root"
    PElist=[]
    Unlist=[]
    print PEdir
    for f in os.listdir(PEdir):     
      if os.path.getsize(PEdir+"/"+f)>=5000:
        PElist.append(PEdir+"/"+f)
    print len(PElist)
    
    mean=[]
    Uncern=[]
    for i in range(76):
      Uncern.append(0.0)
      mean.append(0.0)

    for f in PElist:
      fPE=TFile(f)
      hPE=fPE.Get("swiftBkgNominal")
      for i in range(1,77):
        mean[i-1]=mean[i-1]+hPE.GetBinContent(i)
      fPE.Close()
     
    for i in range(76):
      mean[i]=mean[i]/(len(PElist)*1.00)
     

    
    for f in PElist:
      fPE=TFile(f)
      hPE=fPE.Get("swiftBkgNominal")
      for i in range(1,77):
        Uncern[i-1]=Uncern[i-1]+pow((hPE.GetBinContent(i)-mean[i-1]),2)
      fPE.Close()
   
    for i in range(76):
      Uncern[i]=math.sqrt(Uncern[i]/(1.0*len(PElist)))/(hbkg.GetBinContent(i+1))
    for i in range(76):
      print Uncern[i],hbkg.GetBinLowEdge(i+1)
    #print Uncern
    fbkg.Close()
    
    
    bins= [1133, 1166, 1200, 1234, 1269, 1305, 1341, 1378, 1416, 1454, 1493, 1533, 1573, 1614, 1656, 1698, 1741, 1785, 1830, 1875, 1921, 1968, 2016, 2065, 2114, 2164, 2215, 2267, 2320, 2374, 2429, 2485, 2542, 2600, 2659, 2719, 2780, 2842, 2905, 2969, 3034, 3100, 3167, 3235, 3305, 3376, 3448, 3521, 3596, 3672, 3749, 3827, 3907, 3988, 4070, 4154, 4239, 4326, 4414, 4504, 4595, 4688, 4782, 4878, 4975, 5074, 5175, 5277, 5381, 5487, 5595, 5705, 5817, 5931, 6047, 6165]
    h1 = TH1D( args.title,args.title,len(bins)-1,array('f', bins) )
    for i in range(1,77):
    #  print i,h1.GetBinLowEdge(i)
      h1.SetBinContent(i,Uncern[i-1])
    outFile = TFile.Open(args.outfile,'UPDATE')
    outFile.cd()
    h1.Write(h1.GetName(),TObject.kOverwrite)
    outFile.Close()
    





   







if __name__=="__main__":
    FitUn4p()







