TEMPLATE = lib

CONFIG += qt
QT -= gui

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib qt_include zlib_include zlib_lib

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

HEADERS = file.h \
          JPL_asteroid.h \
          JPL_comet.h \
          MPC.h \
          MPC_asteroid.h \
          MPC_comet.h \
          MPC_obscode.h \
          MPC_observations.h \
          RWO.h

SOURCES = file.cpp \
          JPL_asteroid.cpp \
          JPL_comet.cpp \
          MPC.cpp \
          MPC_asteroid.cpp \
          MPC_comet.cpp \
          MPC_obscode.cpp \
          MPC_observations.cpp \
          RWO.cpp
