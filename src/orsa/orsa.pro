TEMPLATE = subdirs

!win32:!macx SUBDIRS += orsa.dynamiclib.pro
SUBDIRS += orsa.staticlib.pro
