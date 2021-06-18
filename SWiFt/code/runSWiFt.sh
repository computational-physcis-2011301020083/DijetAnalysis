#!/bin/bash

#================================
#
# Make output folders 
#
#================================

mkdir -p ../output/pdf/data
mkdir -p ../output/pdf/mc
mkdir -p ../output/root/data
mkdir -p ../output/root/mc
mkdir -p ../output/log

#=================================================================================
#
# Run SWiFt
# Eg: ./runSWiFt.sh ../configFiles/config_dijets_data.txt
# Eg: ./runSWiFt.sh ../configFiles/config_dijets_data.txt &> ../output/log/log.txt
#
#=================================================================================

# using config from commandline 
CONFIGFILE=$1

root -l -b -q "SWiFt.C(\"$CONFIGFILE\")"


