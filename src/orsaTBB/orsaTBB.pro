TEMPLATE = lib

CONFIG += qt

QT -= gui

CONFIG += tbb_include

include(../../orsa.pri)

INCLUDEPATH += ../
DEPENDPATH  += ../

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = ../../lib/$${PLATFORM_NAME}

HEADERS += malloc.h

SOURCES += malloc.cpp
