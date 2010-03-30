TEMPLATE = lib

CONFIG += qt
QT -= gui

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib zlib_include zlib_lib

include(../../orsa.pri)

INCLUDEPATH += ../
DEPENDPATH  += ../

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../lib/$${PLATFORM_NAME}

unix:!macx {
        LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa
}

macx {
        LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa
}

win32 {
        LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa
}

HEADERS = attitude.h \
          data.h \
          datetime.h \
          gauss.h \
          gmst.h \
          obleq.h \
          observation.h \
          observatory.h \
          orbit.h \
		  print.h

SOURCES = attitude.cpp \
          datetime.cpp \
          gauss.cpp \
          gmst.cpp \
          obleq.cpp \
          observatory.cpp \
		  orbit.cpp \
          sdncal.h \
          dow.c \
          gregor.c \
          julian.c
