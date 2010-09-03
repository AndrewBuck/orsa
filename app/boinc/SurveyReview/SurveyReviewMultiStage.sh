#!/bin/bash

# expect 1 argument
if [ $# -eq 0 ]
then
    echo "Usage: `basename $0` <list-of-field-files>"
    exit 1
fi

for file in "$@"
  do
# first remove path
  file=`basename $file`
# then remove extension
  base=${file%\.*}
  STAGE_CALL="./SurveyReviewStage.sh $base"
  echo calling $STAGE_CALL
  $STAGE_CALL
done
