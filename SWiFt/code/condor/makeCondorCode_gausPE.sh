#!/bin/bash

N_JOBS=100
ExeFile="job.sh"
JobFile="job"

#=====================
#
# Set paths 
#
#=====================

# main SWiFt code folder and code path 
#BasePath="/atlas/data19/ksekhon/SWiFt"
BasePath="$PWD/../.."
CodePath=${BasePath}"/code"

# template config 
Config=${BasePath}"/configFiles/condor/config_dijets_PE_condor.txt" 

# SWiFt background
swiftBkg=$BasePath"/output/root/data/dijet_gausRes_LLHScan.root"
swiftBkgName="swiftBkgBest"

# Signal spline path 
# ***** NOTE: PAY ATTENTION TO SIGNAL PICKED *****
SigPath=${BasePath}"/morphing/output/gaus_dijet"
SigFile="gaus_res.root"
SigJES1_1up=""
SigJES2_1up=""
SigJES3_1up=""
SigJES1_1down=""
SigJES2_1down=""
SigJES3_1down=""
SigType="GFunc" #  ******NOTE: This should be set to "GrLFunc" for q*, W', etc and "GFunc" for gaussian functions!!!!

# outputs
OutputPath="$PWD/../../condor"
OutputFolder="gausRes"
JobDir="${OutputPath}/$OutputFolder"

# Pseudo-experiment location
PEPath=${JobDir}"/PE"

#======================================================
#
# Print summary of paths and ask if they are correct 
#
#=====================================================
echo 
echo "Please check paths, signal files, etc:  "
echo "======================================= "
echo 
echo "    Base folder     = ${BasePath}"
echo
echo "    Template config = ${Config}"
echo 
echo "    SigType         = ${SigType}"
echo "    Sig spline file = ${SigPath}/${SigFile}" 
echo "           JES1_1up = ${SigPath}/${SigJES1_1up}"
echo "           JES2_1up = ${SigPath}/${SigJES2_1up}"
echo "           JES3_1up = ${SigPath}/${SigJES3_1up}"
echo "         JES1_1down = ${SigPath}/${SigJES1_1down}"
echo "         JES2_1down = ${SigPath}/${SigJES2_1down}"
echo "         JES3_1down = ${SigPath}/${SigJES3_1down}"
echo 
echo "    Output folder   = ${JobDir}"
echo 
echo -e "\e[31mPlease check and re-check the paths above before the $N_JOBS jobs are submitted.\e[0m "
echo -e "\e[31mIf everything looks good, please type 'yes'. If not, then type anything else to quit.\e[0m"
echo -n "[yes/no] "

read askPaths
if [ $askPaths != yes ]
then 
  echo "Exiting ... "
  exit
fi
  
#================================
#
# Make all the needed folders 
#
#================================
echo ""
echo "Making folders ..."
echo "=================="
mkdir -p $JobDir/condor
mkdir -p $JobDir/log 
mkdir -p $JobDir/pdf 
mkdir -p $JobDir/root 
mkdir -p $JobDir/configs

#=====================================
#
# Generating PEs from SWiFt background
#
#====================================
echo 
echo "Generating pseudo-experiments from SWiFt background ..."
echo "======================================================="
echo 
echo "    SWiFt Bkg       = $swiftBkg " 
echo "    Name            = $swiftBkgName " 
echo 
echo "    PE output path  = $PEPath "     
echo 
echo -e "\e[31mIs the SWiFt bkg information above correct?\e[0m"
echo -e "\e[31mIf yes, please type 'yes'. If no, then type anything else to quit.\e[0m"
echo -n "[yes/no] "

read askPE
if [ $askPE != yes ]
then 
  echo "Exiting ... "
  exit
fi 

mkdir -p $PEPath
root -l -b -q "doPE.C($N_JOBS, \"$swiftBkg\", \"$swiftBkgName\", \"$PEPath\")"
 
#================================
#
# Prepare config files  
#
#================================
echo 
echo "Making ${N_JOBS} config files in $JobDir/configs ..."
echo "====================================================" 

for i in `seq 1 $N_JOBS`; do
  PENUM=$i PEPATH=$PEPath SIGPATH=$SigPath SIGTYPE=$SigType SIGFILE=$SigFile SIGJES1_1up=$SigJES1_1up SIGJES2_1up=$SigJES2_1up SIGJES3_1up=$SigJES3_1up SIGJES1_1down=$SigJES1_1down SIGJES2_1down=$SigJES2_1down SIGJES3_1down=$SigJES3_1down OUTPATH=$JobDir envsubst < $Config > $JobDir/configs/config$i 
done

#================================
#
# Create Condor job  
#
#================================
echo 
echo "Creating condor job $JobDir/$JobFile ... "
echo "========================================="
  
# the following gets output-ed in the condor run file called "job" 
echo "Executable    = $JobDir/$ExeFile"                       > $JobDir/$JobFile
echo "Error         = $JobDir/condor/serr.\$(Process)"       >> $JobDir/$JobFile
echo "Output        = $JobDir/condor/sout.\$(Process)"       >> $JobDir/$JobFile
echo "Log           = $JobDir/condor/slog.\$(Process)"       >> $JobDir/$JobFile
echo                                                         >> $JobDir/$JobFile

for i in `seq 1 $N_JOBS` ; do
  echo "Arguments = $JobDir/configs/config$i $i"             >> $JobDir/$JobFile 
  echo '+JobFlavour = "workday"'                             >> $JobDir/$JobFile
  echo 'Queue'                                               >> $JobDir/$JobFile
  echo                                                       >> $JobDir/$JobFile 
done 
 

#================================
#
# Create runscript file  
#
#================================

# this is what the condor job will run
echo 
echo "Creating runscript file $JobDir/$ExeFile ..."
echo "============================================"

CMD="root -l -b -q \"SWiFt.C('-')\" < \$1" 

echo "#!/bin/bash"                                            > $JobDir/$ExeFile
echo ''                                                      >> $JobDir/$ExeFile
#echo "export HOME=$CodePath"                                 >> $JobDir/$ExeFile
#echo "source /usr/local/bin/setup/cvmfs_atlas.sh"            >> $JobDir/$ExeFile
echo "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase"  >> $JobDir/$ExeFile
echo "source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh" >> $JobDir/$ExeFile
echo "cd $CodePath"                                          >> $JobDir/$ExeFile
echo "localSetupROOT"                                        >> $JobDir/$ExeFile 
echo                                                         >> $JobDir/$ExeFile
echo "$CMD &> ${JobDir}/log/log\${2}.txt"                    >> $JobDir/$ExeFile
chmod a+x $JobDir/$ExeFile

 
#================================
#
# Launch onto Condor!   
#
#================================
echo 
echo "Launching jobs onto Condor!"
echo "==========================="
condor_submit $JobDir/$JobFile

echo 
echo "Done"
echo 

