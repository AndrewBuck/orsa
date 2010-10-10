TEMPLATE = app

TARGET   = qorsa

CONFIG += qt

CONFIG += qt_include     gmp_include gmp_lib      osg_include osg_lib         gsl_include gsl_lib 

include(../../orsa.pri)

INCLUDEPATH += ../../src/
DEPENDPATH  += ../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../bin/$${PLATFORM_NAME}

unix:!macx {
	LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem    -lgmp -lgmpxx          -losgViewer -losgText -losgUtil -losgGA -losgdb_freetype -losgDB -losg
}

macx {
        LIBS += 
}

win32 {
        LIBS += 
}

HEADERS = MainWindow.h \
	  AddObjectsWindow.h \
	  NewIntegrationWindow.h \
	  BodyTableModel.h

SOURCES = Qorsa.cpp \
	  MainWindow.cpp \
	  AddObjectsWindow.cpp \
	  NewIntegrationWindow.cpp \
	  BodyTableModel.cpp

