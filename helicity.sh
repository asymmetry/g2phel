#!/bin/bash
DIRPREFIX=${0%/*}

ARGS=$(getopt -o e:hr -l event:,help,root -n "helicity.sh" -- "$@")
ORI_ARGS="$@"
eval set -- "$ARGS"

function echo_help() {
  echo -e "usage: helicity.sh [options] RUNNUMBER"
  echo -e "  -e, --event=EVENT_AMOUNT   Set event limit"
  echo -e "  -h, --help                 This small usage guide"
  echo -e "  -r, --root                 Set to insert rootfile"
}

EVENT_AMOUNT=-1
DO_INSERT=0

while true; do
  case "$1" in       
    -e|--evt)         EVENT_AMOUNT="$2"; shift 2;;
    -h|--help)        echo_help; exit 0;;
    -r|--root)        DO_INSERT=1; shift;;
    --)               shift; break;;
    *)                echo_help; exit -1;;
  esac
done

if [[ ! -n "$1" ]]; then
  echo_help; exit -1
fi

if [[ "$1" -lt 20000 ]]; then
  CFGFILE=LHRS.cfg
elif [[ "$1" -lt 40000 ]]; then
  CFGFILE=RHRS.cfg
fi

echo -e "$DIRPREFIX/hel_decode -c $CFGFILE -e $EVENT_AMOUNT $1"
$DIRPREFIX/hel_decode -c $CFGFILE -e $EVENT_AMOUNT "$1"
echo
echo -e "$DIRPREFIX/hel_ring -c $CFGFILE $1"
$DIRPREFIX/hel_ring -c $CFGFILE "$1"
echo
echo -e "$DIRPREFIX/hel_happex -c $CFGFILE $1"
$DIRPREFIX/hel_happex -c $CFGFILE "$1"
echo
echo -e "$DIRPREFIX/hel_tir -r -c $CFGFILE $1"
$DIRPREFIX/hel_tir -r -c $CFGFILE "$1"
echo
echo -e "$DIRPREFIX/hel_alignr -c $CFGFILE $1"
$DIRPREFIX/hel_alignr -c $CFGFILE "$1"
echo
echo -e "$DIRPREFIX/hel_alignh -c $CFGFILE $1"
$DIRPREFIX/hel_alignh -c $CFGFILE "$1"
if [[ $DO_INSERT == 1 ]]; then
  echo
  echo -e "$DIRPREFIX/hel_insert -c $CFGFILE $1"
  $DIRPREFIX/hel_insert -c $CFGFILE "$1"
else
  echo
  echo -e "tar cvzf hel_$1.tar.gz hel_$1.dat helRIN_$1.dat helHAP_$1.dat"
  tar cvzf "hel_$1.tar.gz" "hel_$1.dat" "helRIN_$1.dat" "helHAP_$1.dat"
fi
