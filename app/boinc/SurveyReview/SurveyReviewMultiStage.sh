#!/bin/bash

# cool countdown function 
# found at http://www.unix.com/shell-programming-scripting/98889-display-runnning-countdown-bash-script.html
countdown()
(
  IFS=:
  set -- $*
  secs=$(( ${1#0} * 3600 + ${2#0} * 60 + ${3#0} ))
  while [ $secs -gt 0 ]
  do
    sleep 1 &
    printf "\rwill continue in %02d:%02d:%02d" $((secs/3600)) $(( (secs/60)%60)) $((secs%60))
    secs=$(( $secs - 1 ))
    wait
  done
  echo
)

# expect 1 argument
if [ $# -eq 0 ]
then
    echo "Usage: `basename $0` <list-of-field-files>"
    exit 1
fi

# first list all arguments, to double-check
MSG="calling SurveyReviewStage.sh on $# argument(s):"
for file in "$@"
  do
# first remove path
  file=`basename $file`
# then remove extension
  base=${file%%.*}
  MSG="$MSG $base"
done
echo $MSG
# wait 10 seconds
countdown "00:00:10"

for file in "$@"
  do
# first remove path
  file=`basename $file`
# then remove extension
  base=${file%%.*}
  STAGE_CALL="./SurveyReviewStage.sh $base"
  echo calling $STAGE_CALL
  $STAGE_CALL 
  if [ $? -ne 0 ] ; then
      echo error: failed to call $STAGE_CALL
      exit 1
  fi
done
