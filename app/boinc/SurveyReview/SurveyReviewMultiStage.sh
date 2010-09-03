#!/bin/bash

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
  base=${file%\.*}
  MSG="$MSG $base"
done
echo $MSG
echo will continue in 10 seconds
sleep 10

for file in "$@"
  do
# first remove path
  file=`basename $file`
# then remove extension
  base=${file%\.*}
  STAGE_CALL="./SurveyReviewStage.sh $base"
  echo calling $STAGE_CALL
  $STAGE_CALL 
  if [ $? -ne 0 ] ; then
      echo error: failed to call $STAGE_CALL
      exit 1
  fi
done
