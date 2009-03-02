TEMPLATE = lib

#CONFIG += staticlib
#CONFIG -= qt
#
CONFIG += qt
QT -= gui

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib tbb_include

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

HEADERS = bioCallback.h \
          FindNamedNodeVisitor.h \
          LightSourceUpdateCallback.h \
          OrbitNodeCallback.h \
          viz.h \
          AnimationTime.h \
          BodyTranslationCallback.h \
          CenterOfMassPositionCallback.h \
          BodyAttitudeCallback.h \
          Track.h \
          HUD.h \
          HUDCallback.h

SOURCES = bioCallback.cpp \
          OrbitNodeCallback.cpp \
          viz.cpp \
          AnimationTime.cpp \
          BodyTranslationCallback.cpp \
          Track.cpp \
          HUD.cpp \
          HUDCallback.cpp
