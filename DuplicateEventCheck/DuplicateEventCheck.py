from ROOT import *
import numpy as np
import argparse,math,os,glob,h5py
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--path", dest='path',  default="", help="path")
args = parser.parse_args()

filepath=args.path
f=h5py.File(filepath)
info=f.get("data")

run_info=info[:,1]
print np.where(np.bincount(run_info) > 1)[0]



