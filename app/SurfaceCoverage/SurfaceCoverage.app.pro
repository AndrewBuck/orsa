TEMPLATE = app

CONFIG += qt
QT     += gui
#QT += opengl

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include osg_lib qwt_include qwt_lib spice_include spice_lib

win32 {
	CONFIG += osg_lib
}

include(../../orsa.pri)

TARGET = SurfaceCoverage

INCLUDEPATH += ../../src/
DEPENDPATH  += ../../src/

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../bin/$${PLATFORM_NAME}

unix:!macx {
	LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaOSG -lorsaQt -lorsaSolarSystem -lorsaSPICE -losg -losgText -losgGA -losgViewer
}

macx {
        QMAKE_LFLAGS += -static-libgcc
        LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaOSG -lorsaQt -lorsaSolarSystem -lorsaSPICE -losg -losgText -losgGA -losgViewer
}

win32 {
        LIBS += -L../../lib/$${PLATFORM_NAME} -lorsa -lorsaSolarSystem -lorsaSPICE
}

HEADERS += SurfaceCoverage.h   SurfaceCoverageVersion.h   mainThread.h   multiminPhase.h   plot2D.h   plotTimeline.h   FaceColor.h   vesta.h
SOURCES += SurfaceCoverage.cpp SurfaceCoverageVersion.cpp mainThread.cpp multiminPhase.cpp plot2D.cpp plotTimeline.cpp FaceColor.cpp main.cpp

unix:!macx { BUILD_COUNTER_BIN = ./BuildCounter.$${PLATFORM_NAME} }
macx:      { BUILD_COUNTER_BIN = ./BuildCounter.$${PLATFORM_NAME}.app/Contents/MacOS/BuildCounter.$${PLATFORM_NAME} }

counter.depends  = $${SOURCES} $${HEADERS}
counter.target   = BuildCounter.h
counter.commands = $${BUILD_COUNTER_BIN}

QMAKE_EXTRA_TARGETS += counter
