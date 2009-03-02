include(../../../orsa.pri)

TEMPLATE = subdirs

SUBDIRS += SurveySimulator.app.pro

unix:!macx {
	SUBDIRS += SurveySimulatorWorkGenerator.app.pro multipleJobSubmission.app.pro SurveySimulatorAssimilator.app.pro
}

CONFIG  += ordered

