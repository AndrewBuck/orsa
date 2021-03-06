#!/bin/bash

LOGFILE=SurveyReviewStage.log

msg() {
    echo $*
    echo `date` -- $* >> $LOGFILE
}

# expect 1 argument
if [ $# -ne 1 ]
then
    echo "Usage: `basename $0` <baseName>"
    exit 1
fi

# configure these parameters
BATCH_ID=R3A_1.9-2.3 # name of the batch
STAGE_DIR=./stage # directory containing the input files
# parameters for SurveyReviewMultipleJobsSubmission
MJS_baseName=$1_$BATCH_ID
MJS_samplesPerBin=100
MJS_a_AU_min=1.9
MJS_a_AU_max=2.3
MJS_e_min=0.0
MJS_e_max=1.0
MJS_i_DEG_min=00.0 
MJS_i_DEG_max=90.0
MJS_H_min=16
MJS_H_max=22

# lock, only one can run at any time
if [ -f .stage.lock ]
then
    msg "Another session of `basename $0` is running (or remove stale file .stage.lock)"
else
# lock
    touch .stage.lock 
# copy files from STAGE_DIR
    cp $STAGE_DIR/$1.txt field.dat
    if [ $? -ne 0 ] ; then
	msg error: cannot copy file $STAGE_DIR/$1.txt
	rm .stage.lock 
	exit 1
    else
	msg copied file $STAGE_DIR/$1.txt
    fi
    cp $STAGE_DIR/$1.fieldTime.dat fieldTime.dat
    if [ $? -ne 0 ] ; then
	msg error: cannot copy file $STAGE_DIR/$1.fieldTime.dat
	rm .stage.lock 
	exit 1
    else 
	msg copied file $STAGE_DIR/$1.fieldTime.dat
    fi
    cp $STAGE_DIR/$1.fit.dat fit.dat 
    if [ $? -ne 0 ] ; then
	msg error: cannot copy file $STAGE_DIR/$1.fit.dat
	rm .stage.lock 
	exit 1
    else
	msg copied file $STAGE_DIR/$1.fit.dat
    fi 
# call MJS
    MJS_CALL="./SurveyReviewMultipleJobsSubmission $MJS_baseName $MJS_samplesPerBin $MJS_a_AU_min $MJS_a_AU_max $MJS_e_min $MJS_e_max $MJS_i_DEG_min $MJS_i_DEG_max $MJS_H_min $MJS_H_max"
    msg calling $MJS_CALL
    $MJS_CALL
    if [ $? -ne 0 ] ; then
	msg error: failed to run SurveyReviewMultipleJobsSubmission
	rm .stage.lock 
	exit 1
    fi
#unlock
    rm .stage.lock 
    if [ $? -ne 0 ] ; then
	msg error: failed to delete lock file .stage.lock 
	exit 1
    fi
fi
