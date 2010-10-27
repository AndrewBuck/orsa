include(../../../orsa.pri)

TEMPLATE = subdirs

SUBDIRS += SurveyReview.app.pro 

unix:!macx {
	SUBDIRS += ExtractObservations.app.pro
	SUBDIRS += GenMake.app.pro
	SUBDIRS += GenSkipFile.app.pro
	SUBDIRS += DetectionEfficiency.app.pro
	SUBDIRS += DetectionEfficiencyFit.app.pro
	SUBDIRS += DetectionEfficiencyPlot.app.pro
	SUBDIRS += SurveyReviewWorkGenerator.app.pro
	SUBDIRS += SurveyReviewMultipleJobsSubmission.app.pro
	SUBDIRS += SurveyReviewValidator.app.pro
	SUBDIRS += SurveyReviewAssimilator.app.pro
	SUBDIRS += SurveyReviewMerge.app.pro
	SUBDIRS += SurveyReviewMergeNEOs.app.pro
	SUBDIRS += InspectMerge.app.pro
	SUBDIRS += InspectMergePlot.app.pro
	SUBDIRS += Inspect.app.pro
	SUBDIRS += InspectPlot.app.pro
	SUBDIRS += LunationCopy.app.pro
}
