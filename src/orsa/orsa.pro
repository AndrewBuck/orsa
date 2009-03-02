TEMPLATE = lib

#CONFIG += staticlib
#CONFIG -= qt
#
CONFIG += qt

QT -= gui

CONFIG += gmp_include gmp_lib gsl_include gsl_lib osg_include tbb_include

include(../../orsa.pri)

INCLUDEPATH += ../
DEPENDPATH  += ../

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../lib/$${PLATFORM_NAME}

unix:!macx {
        LIBS += -lrt
}

HEADERS = angle.h \
#          attitude.h \
          body.h \
          bodygroup.h \
#          borderies.h \
          box.h \
          cache.h \
          covariance.h \
          datetime.h \
          debug.h \
          double.h \
          euler.h \
          integrator.h \
          integrator_leapfrog.h \
          integrator_radau.h \
          interaction.h \
          interval.h \
          legendre.h \
          matrix.h \
          multifit.h \
          multimin.h \
          multipole.h \
# moved to orsaSolarSystem observation.h \
          orbit.h \
          paul.h \
          paulMoment.h \
          print.h \
          propulsion.h \
          quaternion.h \
          reference_frame.h \
          shape.h \
          slerp.h \
          spline.h \
          statistic.h \
          unit.h \
          util.h \
          vector.h

SOURCES = angle.cpp \
#          attitude.cpp \
          body.cpp \
          bodygroup.cpp \
#          borderies.cpp \
          box.cpp \ 
          covariance.cpp \
#          datetime.cpp \
          debug.cpp \
          double.cpp \
          euler.cpp \
          integrator.cpp \
          integrator_leapfrog.cpp \
          integrator_radau.cpp \
          legendre.cpp \
          matrix.cpp \
          multifit.cpp \
          multimin.cpp \
          multipole.cpp \
          orbit.cpp \
          paul.cpp \
          paulMoment.cpp \
          propulsion.cpp \
          quaternion.cpp \
          reference_frame.cpp \
          shape.cpp \
          unit.cpp \
          util.cpp \
          vector.cpp


isEmpty(USE_TBB_ALGO) {
	SOURCES += interaction.cpp
} else {
	SOURCES += interaction.TBB.cpp
}
