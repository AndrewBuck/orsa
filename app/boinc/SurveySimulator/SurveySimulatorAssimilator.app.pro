TEMPLATE = app

TARGET   = SurveySimulatorAssimilator

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
	INCLUDEPATH += /home/buck/Code/orsa/boinc/ /home/buck/Code/orsa/boinc/lib/ /home/buck/Code/orsa/boinc/tools/ /home/buck/Code/orsa/boinc/sched/ /home/buck/Code/orsa/boinc/db/ /usr/include/mysql/
	LIBS += -lmysqlclient
}

macx {
        LIBS += 
}

win32 {
        LIBS += -L../../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem -lorsaSPICE
}

# BOINC sources needed
#
SOURCES += /home/buck/Code/orsa/boinc/sched/assimilator.cpp /home/buck/Code/orsa/boinc/sched/validate_util.cpp /home/buck/Code/orsa/boinc/db/boinc_db.cpp /home/buck/Code/orsa/boinc/db/db_base.cpp /home/buck/Code/orsa/boinc/sched/sched_msgs.cpp /home/buck/Code/orsa/boinc/sched/sched_util.cpp /home/buck/Code/orsa/boinc/sched/sched_config.cpp

# local
#
HEADERS +=   
SOURCES += SurveySimulatorAssimilator.cpp 
