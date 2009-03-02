TEMPLATE = lib

#CONFIG += staticlib

CONFIG += qt
QT     += gui opengl

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib tbb_include

include(../../orsa.pri)

INCLUDEPATH += ../ /Users/tricaric/Universal/src/gmp/include
DEPENDPATH  += ../ /Users/tricaric/Universal/src/gmp/include

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

HEADERS = debug.h \
          event.h

SOURCES = debug.cpp
