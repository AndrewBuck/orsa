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

exists(local_qmake_include.pri) {
	message("local_qmake_include.pri file exists, reading it in...")
	include(local_qmake_include.pri)
}
else {
	message("local_qmake_include.pri does not exist.  This file is needed to create the Makefiles that will allow the project to be built.")
	message("")
	message("You can create this file by executing the following command:")
	message("")
	message("echo -en \"ORSA_PRI = $$PWD/orsa.pri\" > local_qmake_include.pri")
	message("")
	error("Exiting the build process because the local_qmake_include.pri file could not be found, please see the above instructions to remedy this.")
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

	QMAKE_CXXFLAGS_DEBUG += -g -pg -ggdb
	QMAKE_LFLAGS_DEBUG   += -g -pg -ggdb
}
macx {
	CONFIG += x86 # x86 ppc

	QMAKE_MAC_SDK = /Developer/SDKs/MacOSX10.5.sdk

        KERNEL_NAME   = $$system(uname -s)
        HARDWARE_NAME = $$system(uname -p)
        PLATFORM_NAME = $${KERNEL_NAME}_$${HARDWARE_NAME}

	DIR_SEP = "/"

        QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.5
	QMAKE_LFLAGS += -undefined dynamic_lookup 

	QMAKE_CXXFLAGS += -mmacosx-version-min=10.5 -read_only_relocs suppress
	QMAKE_LFLAGS   += -mmacosx-version-min=10.5 -read_only_relocs suppress

#	QMAKE_CXXFLAGS_DEBUG += -g -pg
#	QMAKE_LFLAGS_DEBUG   += -g -pg
}
win32 {
	PLATFORM_NAME = "win32"
	DIR_SEP = "\\"

#	QMAKE_CXXFLAGS_RELEASE += -g
#	QMAKE_LFLAGS_RELEASE   += -g

#	QMAKE_CXXFLAGS_DEBUG += -g -pg -ggdb
#	QMAKE_LFLAGS_DEBUG   += -g -pg -ggdb

}
#message([$$PLATFORM_NAME])


ORSA_BASE = $$dirname(ORSA_PRI)
message(ORSA_BASE = $$ORSA_BASE)

SUPPORT_DIR = $$DIR_SEP"support"$$DIR_SEP$$PLATFORM_NAME

ORSA_SUPPORT = $$ORSA_BASE$$SUPPORT_DIR
#message(ORSA_SUPPORT = $$ORSA_SUPPORT)

# 3rd party library configuration, platform by platform (unix,macx,win32)

unix:!macx {
	boinc_include {
		INCLUDEPATH += $$ORSA_BASE/../boinc/ $$ORSA_BASE/../boinc/api/ $$ORSA_BASE/../boinc/lib
	}
	boinc_lib {
		LIBS += $$ORSA_BASE/../boinc/api/libboinc_api.a $$ORSA_BASE/../boinc/lib/libboinc.a
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
        	INCLUDEPATH += $$ORSA_BASE/../OpenSceneGraph/include		
	}
	osg_lib {
		LIBS += -L$$ORSA_BASE/../OpenSceneGraph/lib -L$$ORSA_BASE/../OpenSceneGraph/lib/osgPlugins-2.8.3
	}
	osg_src {
		OSG_SRC = $$ORSA_BASE/../OpenSceneGraph/src/
	}
	qt_include {
		INCLUDEPATH += /usr/include/qt3 /usr/include/qt4
	}
	qt_lib {
		LIBS += -L$$ORSA_BASE/../qwt/lib
	}
	qwt_include {
		INCLUDEPATH += /usr/include/qwt-qt4
	}
	qwt_lib {
		LIBS += -L$$ORSA_BASE/../qwt/lib/ -lqwt
	}
	spice_include {
		INCLUDEPATH += $$ORSA_BASE/../cspice/include
	}
	spice_lib {
		LIBS += $$ORSA_BASE/../cspice/lib/cspice.a $$ORSA_BASE/../cspice/lib/csupport.a
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
		INCLUDEPATH +=  $$ORSA_BASE/../boinc/ $$ORSA_BASE/../boinc/api/ $$ORSA_BASE/../boinc/lib
	}
	boinc_lib {
		LIBS += $$ORSA_BASE/../boinc/mac_build/build/boinc.build/Deployment/libboinc.build/Objects-normal/i386/libboinc.a
		LIBS += $$ORSA_BASE/../boinc/mac_build/build/boinc.build/Deployment/api_libboinc.build/Objects-normal/i386/libboinc_api.a
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
        	INCLUDEPATH += $$ORSA_BASE/../OpenSceneGraph/include
	}
	osg_lib {
		LIBS += -L$$ORSA_BASE/../OpenSceneGraph/lib/
	}
	osg_src {
		OSG_SRC = $$ORSA_BASE/../OpenSceneGraph/src/
	}
	qwt_include {
		INCLUDEPATH += $$ORSA_BASE/../qwt/src
	}
	qwt_lib {
		LIBS += -L$$ORSA_BASE/../qwt/lib/ -lqwt
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
        QMAKE_CXXFLAGS += -DOSG_LIBRARY
	}
	osg_lib {
		LIBS += C:\OpenSceneGraph\lib\libOpenThreads.a C:\OpenSceneGraph\lib\libosg.a
	}
	osg_src {
		OSG_SRC = C:\OpenSceneGraph\src\
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

