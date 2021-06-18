import argparse,math,os
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--path", dest='path',  default="", help="path")
args = parser.parse_args()

f=open(args.path,"r")
outxt=list(f.readlines())
print outxt

name=[]
chi3=[]
chi4=[]

j=0
for i in outxt:
  j=j+1
  if j%4==1:
    #print i
    histname=i.split(" ")[-1].replace("\n","")
    name.append(histname)
  if j%4==3:
    chi3value=i.split(" ")[-1].replace("\n","")
    chi3.append(chi3value)
  if j%4==2:
    chi4value=i.split(" ")[-1].replace("\n","")
    chi4.append(chi4value)
    
for i in chi4:
  print i
'''
for i in outxt:
  if ".root" in i:
    histname=i.split(" ")[-1].replace("\n","")
    name.append(histname)
  elif "nominal" in i:
    chi3value=i.split(" ")[-1].replace("\n","")
    chi3.append(chi3value)
  elif "alternate" in i:
    chi4value=i.split(" ")[-1].replace("\n","")
    chi4.append(chi4value)

print "name"
for i in name:
  print i
print "chi3"
for i in chi3:
  print i
print "chi4"
for i in chi4:
  print i
'''

