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
	INCLUDEPATH += /home/buck/Code/orsa/boinc/ /home/buck/Code/orsa/boinc/lib/ /home/buck/Code/orsa/boinc/tools/ /home/buck/Code/orsa/boinc/sched/ /home/buck/Code/orsa/boinc/db/ /usr/include/mysql/
	LIBS += -lmysqlclient
	LIBS += -L/home/buck/Code/orsa/sqlite -lsqlite3 -ldl -lgmp
}

# BOINC sources needed
#
SOURCES += /home/buck/Code/orsa/boinc/sched/validator.cpp /home/buck/Code/orsa/boinc/db/boinc_db.cpp /home/buck/Code/orsa/boinc/db/db_base.cpp /home/buck/Code/orsa/boinc/sched/sched_msgs.cpp /home/buck/Code/orsa/boinc/sched/sched_config.cpp /home/buck/Code/orsa/boinc/sched/credit.cpp /home/buck/Code/orsa/boinc/sched/sched_util.cpp /home/buck/Code/orsa/boinc/sched/validate_util2.cpp /home/buck/Code/orsa/boinc/sched/validate_util.cpp /home/buck/Code/orsa/boinc/sched/sched_shmem.cpp

# local
#
HEADERS += sqlite3.h sqlite3ext.h
SOURCES += SurveyReviewValidator.cpp 
