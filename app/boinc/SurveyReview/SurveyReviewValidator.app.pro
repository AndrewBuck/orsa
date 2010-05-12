TEMPLATE = app

TARGET   = SurveyReviewValidator

CONFIG -= qt

CONFIG += boinc_include boinc_lib 

win32 {
	CONFIG += osg_lib
}

include(../../../orsa.pri)

INCLUDEPATH += ../../../src/
DEPENDPATH  += ../../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../../bin/$${PLATFORM_NAME}

unix:!macx {
	INCLUDEPATH += /home/tricaric/boinc/ /home/tricaric/boinc/lib/ /home/tricaric/boinc/tools/ /home/tricaric/boinc/sched/ /home/tricaric/boinc/db/ /usr/include/mysql/
	LIBS += -lmysqlclient
}

# BOINC sources needed
#
SOURCES += /home/tricaric/boinc/sched/validator.cpp /home/tricaric/boinc/db/boinc_db.cpp /home/tricaric/boinc/db/db_base.cpp /home/tricaric/boinc/sched/sched_msgs.cpp /home/tricaric/boinc/sched/sched_config.cpp /home/tricaric/boinc/sched/credit.cpp /home/tricaric/boinc/sched/sched_util.cpp /home/tricaric/boinc/sched/validate_util2.cpp /home/tricaric/boinc/sched/validate_util.cpp 

# local
#
HEADERS +=   
SOURCES += SurveyReviewValidator.cpp 
