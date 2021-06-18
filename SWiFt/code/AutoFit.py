from ROOT import *
import argparse,math,os
parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--path", dest='path',  default="", help="path")
parser.add_argument("--save", dest='save',  default="", help="save")
parser.add_argument("--txt", dest='txt',  default="", help="txt")
args = parser.parse_args()

histlist=[]
filepath=args.path
filename=filepath.split("/")[-1].split(".")[0]
f=TFile(filepath,"r")
histname="HT"
for key in f.GetListOfKeys():
  #if histname in key.GetName():
    histlist.append(key.GetName())


confdir="/afs/cern.ch/work/d/ding/ZplusX/SWiFt/configFiles"
rundir="/afs/cern.ch/work/d/ding/ZplusX/SWiFt/code"
conf="config_mass.txt"

savedirname=args.save
savedir1="/afs/cern.ch/work/d/ding/ZplusX/SWiFt/output/pdf"
savedir2="/afs/cern.ch/work/d/ding/ZplusX/SWiFt/output/root"
if not os.path.exists(savedir1+"/"+savedirname):
  os.mkdir(savedir1+"/"+savedirname)
  os.mkdir(savedir2+"/"+savedirname)

for hist in histlist:
  #outname=filename+"_"+hist
  outname=hist
  print outname
  keyline1="inputFile       : "+filepath+"\n"
  keyline2="inputHist       : "+hist+"\n"
  keyline3="pdfName         : pdf/"+savedirname+"/"+outname+"\n"
  keyline4="oRootName       : root/"+savedirname+"/"+outname+"\n"

  os.chdir(confdir)
  confnew=open("configTemp1.txt","w")
  with open(conf, "r") as fin:
    for line in fin:
      if line.startswith("inputFile"):
        line=keyline1
        confnew.write(line)
      elif line.startswith("inputHist "):
        line=keyline2
        confnew.write(line)
      elif line.startswith("pdfName"):
        line=keyline3
        confnew.write(line)
      elif line.startswith("oRootName"):
        line=keyline4
        confnew.write(line)
      else:
        confnew.write(line)
  fin.close()
  confnew.close()

  os.chdir(rundir)
  #os.system("./runSWiFt.sh ../configFiles/configTemp1.txt")
  output = os.popen("./runSWiFt.sh ../configFiles/configTemp1.txt")
  
  txtfile=args.txt+".txt"
  Record=open("./Txtfile/"+txtfile,"a")
  outxt=list(output.readlines())
  #print outxt[-2]  
  Line1=filepath+": "+hist+"\n"
  Line2=outxt[-1]
  Line3=outxt[-2]
  Line4=outxt[-3]
  Record.write(Line1)
  #Record.write(Line2)
  Record.write(Line3)
  Record.write(Line4)  
  Record.write("\n")  
  Record.close()

  os.remove(confdir+"/configTemp1.txt")












