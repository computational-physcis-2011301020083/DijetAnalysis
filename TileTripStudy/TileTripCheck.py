from ROOT import *
import os,shutil,math,glob
from array import array


#Check tile trip
dictLow={"LBA":0.0,"LBC":-0.9,"EBA":0.8,"EBC":-1.7}
dictHigh={"LBA":0.9,"LBC":0.0,"EBA":1.7,"EBC":-0.8}
width=2.0*math.pi/64.0
def check(lumiblock,eta,phi,deadtitle):
    LumiblockStr=deadtitle.split(" ")[0]
    LumiblockLow=int(LumiblockStr.split(",")[0][1:])
    LumiblockHigh=int(LumiblockStr.split(",")[1][:-1])
    test=False
    for i in range(1,len(deadtitle.split(" "))):
      Mod=deadtitle.split(" ")[i][0:3]
      ModNum=int(deadtitle.split(" ")[i][3:5])
      EtaLow=dictLow[Mod]
      EtaHigh=dictHigh[Mod]
      if ModNum<=32:
        PhiLow=float(ModNum)*width
        PhiHigh=PhiLow+width
      elif ModNum<=64:
        PhiLow=float(ModNum)*width-2*math.pi
        PhiHigh=PhiLow+width
      else:
        print "Module Number bigger than 64!"
      DelPhi=(PhiHigh-PhiLow)/2.0
      PhiCenter=PhiLow+DelPhi if abs(DelPhi+PhiLow)<=math.pi else PhiHigh-DelPhi
      Delta=abs(phi-PhiCenter) if abs(phi-PhiCenter)<=math.pi else 2.0*math.pi-abs(phi-PhiCenter)
      if EtaLow<=eta<=EtaHigh and Delta<DelPhi:
          test=True
      else:
          test=False
    if LumiblockLow==1 and  LumiblockHigh==9999:
      return test
    elif LumiblockLow==2 and  LumiblockHigh==9999:
      return test
    else :
      return test and LumiblockLow<=lumiblock<=LumiblockHigh 

#Enter your analysis code below
#You need to confirm if the RunNumber of that event you search is in the bad tile info txt BAD_MODULES_RUN2.txt before you use the function check(lumiblock,eta,phi,deadtitle),  if it is not in BAD_MODULES_RUN2.txt, we can pass the event, it is good, if it is in BAD_MODULES_RUN2.txt, you can perform the check.

for File in datapaths:
  count=count+1
  print "Processing file count ",count
  f=TFile(File)
  t=f.Get("outTree")
  N=t.GetEntries()
  for j in range(0,N):
    t.GetEntry(j)
    if  "HLT_j420" in t.passedTriggers and t.jet_pt[0]>=420 and t.jet_pt[1]>=150 and (abs(t.jet_eta[0]-t.jet_eta[1])/2)<=0.6 and t.njet>=2 and t.mjj>=1100:
      count1=count1+1
      boolist=[]
      with open("BAD_MODULES_RUN2.txt","r") as fline:
        for line in fline:
          RDeadTile=''.join(line).split("  ")
          if int(RDeadTile[0])==t.runNumber:
            for k in range(1,len(RDeadTile)):
              boolist.append(check(t.lumiBlock,t.jet_eta[0],t.jet_phi[0],RDeadTile[k]))
              boolist.append(check(t.lumiBlock,t.jet_eta[1],t.jet_phi[1],RDeadTile[k]))
            if not True in boolist:
              count2=count2+1
              hist.Fill(t.mjj)
      fline.close()







