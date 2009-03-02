TEMPLATE = app

CONFIG -= qt

include(../../orsa.pri)

TARGET   = BuildCounter.$${PLATFORM_NAME}

UI_DIR      =  .ui/$${PLATFORM_NAME}
MOC_DIR     = .moc/$${PLATFORM_NAME}
OBJECTS_DIR = .obj/$${PLATFORM_NAME}
DESTDIR     = .

HEADERS = 
SOURCES = BuildCounter.cpp
