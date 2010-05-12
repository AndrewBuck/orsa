TEMPLATE = app

TARGET   = SurveyReviewValidator

CONFIG -= qt

CONFIG += boinc_include boinc_lib 

include(../../../orsa.pri)

INCLUDEPATH += ../../../src/
DEPENDPATH  += ../../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../../bin/$${PLATFORM_NAME}

unix:!macx {
	CONFIG += gmp_include gmp_lib
	INCLUDEPATH += /home/tricaric/boinc/ /home/tricaric/boinc/lib/ /home/tricaric/boinc/tools/ /home/tricaric/boinc/sched/ /home/tricaric/boinc/db/ /usr/include/mysql/
	LIBS += -lmysqlclient
	LIBS += -L/home/tricaric/sqlite -lsqlite3 -ldl -lgmp
}

# BOINC sources needed
#
SOURCES += /home/tricaric/boinc/sched/validator.cpp /home/tricaric/boinc/db/boinc_db.cpp /home/tricaric/boinc/db/db_base.cpp /home/tricaric/boinc/sched/sched_msgs.cpp /home/tricaric/boinc/sched/sched_config.cpp /home/tricaric/boinc/sched/credit.cpp /home/tricaric/boinc/sched/sched_util.cpp /home/tricaric/boinc/sched/validate_util2.cpp /home/tricaric/boinc/sched/validate_util.cpp 

# local
#
HEADERS += sqlite3.h sqlite3ext.h
SOURCES += SurveyReviewValidator.cpp 
