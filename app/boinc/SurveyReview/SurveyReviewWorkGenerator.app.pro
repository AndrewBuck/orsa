TEMPLATE = app

TARGET   = SurveyReviewWorkGenerator

CONFIG += qt

QT -= gui

CONFIG += boinc_include boinc_lib gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib spice_include spice_lib

include(../../../orsa.pri)

INCLUDEPATH += ../../../src/
DEPENDPATH  += ../../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../../bin/$${PLATFORM_NAME}

unix:!macx {
	INCLUDEPATH += /home/tricaric/boinc/ /home/tricaric/boinc/lib/ /home/tricaric/boinc/tools/ /home/tricaric/boinc/sched/ /home/tricaric/boinc/db/ /usr/include/mysql/
	LIBS += -L../../../lib/$${PLATFORM_NAME} -lorsa -lorsaInputOutput -lorsaSolarSystem -lorsaSPICE -lorsaEssentialOSG -lOpenThreads /home/tricaric/boinc/api/libboinc_api.a /home/tricaric/boinc/lib/libboinc.a -lmysqlclient -lssl -lcrypto
}

# BOINC sources, needed
SOURCES += /home/tricaric/boinc/db/boinc_db.cpp /home/tricaric/boinc/db/db_base.cpp /home/tricaric/boinc/sched/sched_msgs.cpp /home/tricaric/boinc/sched/sched_util.cpp /home/tricaric/boinc/sched/sched_config.cpp /home/tricaric/boinc/tools/backend_lib.cpp /home/tricaric/boinc/tools/process_result_template.cpp /home/tricaric/boinc/lib/crypt.cpp

HEADERS += SurveyReviewWorkGenerator.h
SOURCES += SurveyReviewWorkGenerator.cpp
