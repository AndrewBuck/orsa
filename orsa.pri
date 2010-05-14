# version, for the libraries
# when using VERSION, libs on win32 are called liborsa1 instead of liborsa
# VERSION = 1.0.0

#CONFIG += thread debug warn_on
CONFIG += thread release warn_on

# enable this to build static libs, useful to create binaries for BOINC
#CONFIG += staticlib

macx {
	CONFIG += staticlib
}

# not here... qt and gui are kept local on single directory's .pro files
#CONFIG += qt
#QT -= gui


# define ORSA_PRI as the absolute path to this file
unix:!macx {
	ORSA_PRI = /home/tricaric/orsa/orsa.pri
}
macx {
	ORSA_PRI = /Users/tricaric/orsa/orsa.pri
}
win32 {
	ORSA_PRI = C:\orsa\orsa.pri
}


# important tests, don't change
isEmpty(ORSA_PRI) {
	error(You need to set the ORSA_PRI variable correctly and to review all the other relevant settings in orsa.pri)
}
#
!exists($$ORSA_PRI) {
	error(ORSA_PRI is not set correctly [$$ORSA_PRI])
}


# platform name and other tweaking (unix,macx,win32)
unix:!macx {
	KERNEL_NAME   = $$system(uname -s)
	HARDWARE_NAME = $$system(uname -m)
	PLATFORM_NAME = $${KERNEL_NAME}_$${HARDWARE_NAME}

	DIR_SEP = "/"

	QMAKE_CXXFLAGS_RELEASE +=
	QMAKE_LFLAGS_RELEASE   +=

	QMAKE_CXXFLAGS_DEBUG += -pg -ggdb
	QMAKE_LFLAGS_DEBUG   += -pg -ggdb
}
macx {
	CONFIG += x86 # x86 ppc

	QMAKE_MAC_SDK = /Developer/SDKs/MacOSX10.5.sdk

        KERNEL_NAME   = $$system(uname -s)
        HARDWARE_NAME = $$system(uname -p)
        PLATFORM_NAME = $${KERNEL_NAME}_$${HARDWARE_NAME}

	DIR_SEP = "/"

        QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.5
#       QMAKE_LFLAGS += -undefined dynamic_lookup 

	QMAKE_CXXFLAGS += -mmacosx-version-min=10.5
	QMAKE_LFLAGS   += -mmacosx-version-min=10.5 

#	QMAKE_CXXFLAGS_DEBUG += -pg
#	QMAKE_LFLAGS_DEBUG   += -pg
}
win32 {
	PLATFORM_NAME = "win32"
	DIR_SEP = "\\"
}
#message([$$PLATFORM_NAME])


ORSA_BASE = $$dirname(ORSA_PRI)
#message(ORSA_BASE = $$ORSA_BASE)

SUPPORT_DIR = $$DIR_SEP"support"$$DIR_SEP$$PLATFORM_NAME

ORSA_SUPPORT = $$ORSA_BASE$$SUPPORT_DIR
#message(ORSA_SUPPORT = $$ORSA_SUPPORT)

# 3rd party library configuration, platform by platform (unix,macx,win32)

unix:!macx {
	boinc_include {
		INCLUDEPATH += /home/tricaric/boinc/ /home/tricaric/boinc/api/ /home/tricaric/boinc/lib
	}
	boinc_lib {
		LIBS += /home/tricaric/boinc/api/libboinc_api.a /home/tricaric/boinc/lib/libboinc.a
	}
	gmp_include {
		INCLUDEPATH +=
	}
	gmp_lib {
		LIBS += -lgmp -lgmpxx
	}
	gsl_include {
		INCLUDEPATH +=
	}
	gsl_lib {
		LIBS += -lgsl -lgslcblas
	}
	osg_include {
        	INCLUDEPATH += /home/tricaric/OpenSceneGraph/include		
	}
	osg_lib {
		LIBS += -L/home/tricaric/OpenSceneGraph/lib
	}
	osg_src {
		OSG_SRC = /home/tricaric/OpenSceneGraph/src/
	}
	qwt_include {
		INCLUDEPATH += /home/tricaric/qwt/src
	}
	qwt_lib {
		LIBS += -L/home/tricaric/qwt/lib/ -lqwt
	}
	spice_include {
		INCLUDEPATH += /home/tricaric/cspice/include
	}
	spice_lib {
		LIBS += /home/tricaric/cspice/lib/cspice.a /home/tricaric/cspice/lib/csupport.a
	}
	zlib_include {
        	INCLUDEPATH +=
	}
	zlib_lib {
		LIBS +=
	}
}

macx {
	boinc_include {
		INCLUDEPATH +=  /Users/tricaric/boinc/ /Users/tricaric/boinc/api/ /Users/tricaric/boinc/lib
	}
	boinc_lib {
		LIBS += /Users/tricaric/boinc/api/libboinc_api.a /Users/tricaric/boinc/lib/libboinc.a
	}
	gmp_include {
	     	INCLUDEPATH += $$ORSA_SUPPORT/gmp/include
	}
	gmp_lib {
		LIBS += $$ORSA_SUPPORT/gmp/lib/libgmpxx.a $$ORSA_SUPPORT/gmp/lib/libgmp.a
	}
	gsl_include {
		INCLUDEPATH += $$ORSA_SUPPORT/gsl/include
	}
	gsl_lib {
		LIBS += -L$$ORSA_SUPPORT/gsl/lib -lgsl -lgslcblas
	}
	osg_include {
        	INCLUDEPATH += /Users/tricaric/OpenSceneGraph/include
	}
	osg_lib {
		LIBS += -L/Users/tricaric/OpenSceneGraph/lib/
	}
	osg_src {
		OSG_SRC = /Users/tricaric/OpenSceneGraph/src/
	}
	qwt_include {
		INCLUDEPATH += /Users/tricaric/qwt/src
	}
	qwt_lib {
		LIBS += -L/Users/tricaric/qwt/lib/ -lqwt
	}
	spice_include {
		INCLUDEPATH += $$ORSA_SUPPORT/cspice/include
	}
	spice_lib {
		LIBS += $$ORSA_SUPPORT/cspice/lib/cspice.a $$ORSA_SUPPORT/cspice/lib/csupport.a
	}
	zlib_include {
        	INCLUDEPATH +=
	}
	zlib_lib {
		LIBS +=
	}
}

win32 {
	boinc_include {
		INCLUDEPATH += C:\boinc\boinc C:\boinc\boinc\api C:\boinc\boinc\lib	
	}
	boinc_lib {
		LIBS += C:\boinc\boinc\api\libboinc.a
	}
	gmp_include {
		INCLUDEPATH += $$ORSA_SUPPORT\gmp\include
	}
	gmp_lib {
		LIBS += $$ORSA_SUPPORT\gmp\lib\libgmpxx.a $$ORSA_SUPPORT\gmp\lib\libgmp.a 
	}
	gsl_include {
		INCLUDEPATH += $$ORSA_SUPPORT\gsl\include
	}
	gsl_lib {
		LIBS += -L$$ORSA_SUPPORT\gsl\lib -lgsl -lgslcblas
	}
	osg_include {
        	INCLUDEPATH += C:\OpenSceneGraph\include
#		QMAKE_CXXFLAGS += -DOT_LIBRARY_STATIC -DOSG_LIBRARY_STATIC
	}
	osg_lib {
		LIBS += -LC:\OpenSceneGraph\lib -lOpenThreads -losg
	}
	osg_src {
		OSG_SRC = 
	}
	qwt_include {
		INCLUDEPATH += C:\qwt\src
	}
	qwt_lib {
		LIBS += -LC:\qwt\lib\ -lqwt5
	}
	spice_include {
		INCLUDEPATH += $$ORSA_SUPPORT\cspice\include
	}
	spice_lib {
		LIBS += $$ORSA_SUPPORT\cspice\lib\cspice.a $$ORSA_SUPPORT\cspice\lib\csupport.a
	}
	zlib_include {
        	INCLUDEPATH += C:\Qt\2010.02.1\qt\src\3rdparty\zlib
	}
	zlib_lib {
		LIBS +=
	}
}

