TEMPLATE = app

CONFIG += qt
QT     -= gui 
QT     += opengl

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib qt_include qwt_include qwt_lib spice_include spice_lib

include(../../orsa.pri)

TARGET   = DAWN

INCLUDEPATH += ../../src/
DEPENDPATH  += ../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../bin/$${PLATFORM_NAME}

unix:!macx {
	LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaOSG -lorsaQt -lorsaSolarSystem -lorsaSPICE -losgViewer -losgText -losgUtil -losgGA -losgdb_freetype -losgDB -losg -lOpenThreads -lgmp -lgmpxx
}

macx {
	LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaOSG -lorsaQt -lorsaSolarSystem -lorsaSPICE -losg -losgText -losgGA -losgViewer -lOpenThreads -losg -losgDB -losgUtil
}

HEADERS += DAWN.h   vesta.h  multiminPhase.h
SOURCES += DAWN.cpp main.cpp multiminPhase.cpp


