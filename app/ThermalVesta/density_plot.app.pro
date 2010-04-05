TEMPLATE = app

CONFIG += qt
QT     -= gui opengl

CONFIG += gmp_include gmp_lib gsl_include osg_include osg_lib gsl_lib spice_include spice_lib

include(../../orsa.pri)

TARGET   = density_plot

INCLUDEPATH += ../../src/
DEPENDPATH  += ../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = .

unix:!macx {
	INCLUDEPATH += /home/tricaric/dislin/
	LIBS += -L/home/tricaric/dislin/ -ldislin
	LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaEssentialOSG -lorsaSolarSystem -lorsaSPICE -lOpenThreads
}

HEADERS += 
SOURCES += density_plot.cpp
