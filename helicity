#!/bin/bash
if [[ $# -eq 2 ]];then
  ./extract "$1" "$2"
  ./ring $1 RIN
  if [[ $1 -lt 40000 ]];then
    ./ring $1 HAP
  fi
  ./tir $1 1
  rm -f helRIN_$1.tmp
  rm -f helHAP_$1.tmp
  rm -f helTIR_$1.tmp
  ./align $1
  if [[ $1 -lt 40000 ]];then
      rm -f hel_$1.tmp
  fi
  tar cvzf "hel_$1.tar.gz" "hel_$1.dat" "helTIR_$1.dat" "helRIN_$1.dat" "helHAP_$1.dat"
fi