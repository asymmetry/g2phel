#!/bin/bash
EXP=g2p
TRACK=analysis
#TRACK=debug
OS=linux64

while read i
do
  if [[ -f "/mss/halla/g2p/raw/g2p_$i.dat.0" ]];then
    : > job.$i
    JOB=hel_$i
    echo "PROJECT:     $EXP" >> job.$i
    echo "JOBNAME:     $JOB" >> job.$i
    echo "TRACK:       $TRACK" >> job.$i
    echo "MEMORY:      2000 MB" >> job.$i
    echo "OS:          $OS" >> job.$i
    echo "COMMAND:     /w/halla-sfs62/g2p/chao/helicity/do_helicity" >> job.$i
    echo "OTHER_FILES: " >> job.$i
    echo "/w/halla-sfs62/g2p/chao/helicity/helicity" >> job.$i
    echo "/w/halla-sfs62/g2p/chao/helicity/align" >> job.$i
    echo "/w/halla-sfs62/g2p/chao/helicity/extract" >> job.$i
    echo "/w/halla-sfs62/g2p/chao/helicity/ring" >> job.$i
    echo "/w/halla-sfs62/g2p/chao/helicity/tir" >> job.$i
    echo "SINGLE_JOB:  true" >> job.$i
    echo "INPUT_FILES: " >> job.$i
    for file in $(ls /mss/halla/g2p/raw/g2p_$i.dat.*)
    do
      echo "$file" >> job.$i
    done
    echo "Job file for run $i has been created!"
    jsub job.$i
    mv job.$i jobfiles
  fi
done < list.dat