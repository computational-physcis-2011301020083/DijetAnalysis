from ROOT import *
import numpy as np
import argparse,math,os,glob,h5py
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--path", dest='path',  default="", help="path")
args = parser.parse_args()

datapaths=[]
datadir="/eos/atlas/atlascerngroupdisk/phys-exotics/jdm/dibjet/FullRUN2/LatestData/Data15"
paths=[]
for f in os.listdir(datadir):
  paths=paths+glob.glob(datadir+"/"+f+"/*.root")
datapaths=[]
for f in paths:
  F1=TFile(f)
  if F1.GetListOfKeys().Contains("outTree"):
    datapaths.append(f)

new_hdf5 = h5py.File("Info.h5", 'w')
j=0
for filepath in datapaths:
  print j
  j=j+1
  f=TFile(filepath)
  t=f.Get("outTree")
  a,b=t.AsMatrix(return_labels=True,columns=["runNumber","eventNumber"])
  if j==1:
    c=a
  else:
    c=np.vstack((a,c))
new_hdf5.create_dataset("data",data=c)
print c.shape


