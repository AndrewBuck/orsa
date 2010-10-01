TEMPLATE = app

TARGET   = SurveySimulatorWorkGenerator

CONFIG += qt

QT -= gui

CONFIG += boinc_include boinc_lib gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib spice_include spice_lib

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
	LIBS += -L../../../lib/$${PLATFORM_NAME} -lorsa -lorsaInputOutput -lorsaSolarSystem -lorsaSPICE -lorsaEssentialOSG -lOpenThreads /home/buck/Code/orsa/boinc/api/libboinc_api.a /home/buck/Code/orsa/boinc/lib/libboinc.a -lmysqlclient -lssl -lcrypto
}

macx {
        LIBS += 
}

win32 {
        LIBS += -L../../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem -lorsaSPICE
}

# BOINC sources needed

SOURCES += /home/buck/Code/orsa/boinc/db/boinc_db.cpp /home/buck/Code/orsa/boinc/db/db_base.cpp /home/buck/Code/orsa/boinc/sched/sched_msgs.cpp /home/buck/Code/orsa/boinc/sched/sched_util.cpp /home/buck/Code/orsa/boinc/sched/sched_config.cpp /home/buck/Code/orsa/boinc/tools/backend_lib.cpp /home/buck/Code/orsa/boinc/tools/process_result_template.cpp /home/buck/Code/orsa/boinc/lib/crypt.cpp

HEADERS += parameter.h   telescope.h   SurveySimulator.h   SurveySimulatorWorkGenerator.h
SOURCES += parameter.cpp telescope.cpp SurveySimulator.cpp SurveySimulatorWorkGenerator.cpp
