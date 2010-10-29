TEMPLATE = app

TARGET = SurveyReview

CONFIG += qt

QT -= gui

CONFIG += boinc_include boinc_lib gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib spice_include spice_lib zlib_include

include(../../../orsa.pri)

INCLUDEPATH += ../../../src/
DEPENDPATH  += ../../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../../bin/$${PLATFORM_NAME}

unix:!macx {
	LIBS += -L../../../lib/$${PLATFORM_NAME} -lorsa -lorsaInputOutput -lorsaSolarSystem -lorsaSPICE -lorsaEssentialOSG -lOpenThreads -lorsaUtil /home/buck/Code/orsa/boinc/api/libboinc_api.a /home/buck/Code/orsa/boinc/lib/libboinc.a -L/home/buck/Code/orsa/sqlite -lsqlite3
}

macx {
	QMAKE_LFLAGS += -static-libgcc
	LIBS += -L../../../lib/$${PLATFORM_NAME} -lorsa -lorsaInputOutput -lorsaSolarSystem -lorsaSPICE -lorsaEssentialOSG -lOpenThreads -lorsaUtil
	SOURCES += /Users/buck/Code/orsa/sqlite/sqlite-3.6.23.1/sqlite3.c
}

win32 {
	QMAKE_LFLAGS += -static -static-libgcc
	LIBS += -L../../../lib/$${PLATFORM_NAME} -lorsa -lorsaInputOutput -lorsaSolarSystem -lorsaSPICE -lorsaEssentialOSG -lorsaUtil
	LIBS += C:\OpenSceneGraph\lib\libOpenThreads.a
	LIBS += C:\orsa\support\win32\cspice\lib\cspice.a C:\orsa\support\win32\cspice\lib\csupport.a
	LIBS += C:\orsa\support\win32\gmp\lib\libgmp.a C:\orsa\support\win32\gmp\lib\libgmpxx.a
	LIBS += C:\orsa\support\win32\gsl\lib\libgsl.a C:\orsa\support\win32\gsl\lib\libgslcblas.a
	LIBS += C:\orsa\lib\win32\liborsaSolarSystem.a
	LIBS += C:\orsa\lib\win32\liborsa.a C:\orsa\lib\win32\liborsaSPICE.a C:\orsa\lib\win32\liborsaInputOutput.a
	LIBS += C:\orsa\lib\win32\liborsaUtil.a
	LIBS += C:\orsa\support\win32\gmp\lib\libgmp.a C:\orsa\support\win32\gmp\lib\libgmpxx.a
	LIBS += C:\orsa\lib\win32\liborsaEssentialOSG.a
	LIBS += C:\orsa\support\win32\gsl\lib\libgsl.a C:\orsa\support\win32\gsl\lib\libgslcblas.a
	LIBS += C:\orsa\support\win32\cspice\lib\cspice.a C:\orsa\support\win32\cspice\lib\csupport.a
	LIBS += C:\orsa\lib\win32\liborsaSolarSystem.a
	LIBS += C:\Qt\2010.02.1\qt\lib\libQtCore.a
        HEADERS += C:\sqlite-amalgamation-3.6.23\sqlite3.h C:\sqlite-amalgamation-3.6.23\sqlite3ext.h
	SOURCES += C:\sqlite-amalgamation-3.6.23\sqlite3.c
}

HEADERS += skycoverage.h   SurveyReview.h   grain.h
SOURCES += skycoverage.cpp SurveyReview.cpp grain.cpp SurveyReviewMain.cpp

HEADERS += sqlite3.h sqlite3ext.h
