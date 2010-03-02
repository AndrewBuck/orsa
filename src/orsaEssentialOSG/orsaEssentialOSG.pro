TEMPLATE = lib

CONFIG -= qt

CONFIG += osg_include osg_src

include(../../orsa.pri)

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../lib/$${PLATFORM_NAME}

SOURCES += $${OSG_SRC}/osg/ApplicationUsage.cpp \
           $${OSG_SRC}/osg/DeleteHandler.cpp \
           $${OSG_SRC}/osg/Matrixd.cpp \ 
           $${OSG_SRC}/osg/Matrixf.cpp \
           $${OSG_SRC}/osg/Notify.cpp \
           $${OSG_SRC}/osg/Observer.cpp \
           $${OSG_SRC}/osg/Quat.cpp \
           $${OSG_SRC}/osg/Referenced.cpp
