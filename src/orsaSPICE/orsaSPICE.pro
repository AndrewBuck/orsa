TEMPLATE = lib

# this staticlib option is needed expecially on x86-64 boxes
#CONFIG += staticlib
#
CONFIG += qt
QT -= gui

CONFIG += spice_include spice_lib gmp_include gmp_lib gsl_include gsl_lib osg_include tbb_include

include(../../orsa.pri)

INCLUDEPATH += ../
DEPENDPATH  += ../

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../lib/$${PLATFORM_NAME}

unix:!macx {
        LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem
}

macx {
        LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem
}

win32 {
        LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem
}

HEADERS = spice.h \
          spiceBodyAttitudeCallback.h \
          spiceBodyInterpolatedPosVelCallback.h \
          spiceBodyPosVelCallback.h

SOURCES = spice.cpp \
          spiceBodyAttitudeCallback.cpp \
          spiceBodyInterpolatedPosVelCallback.cpp \
          spiceBodyPosVelCallback.cpp
