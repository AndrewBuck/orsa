TEMPLATE = subdirs

!win32:!macx SUBDIRS += orsaQt.dynamiclib.pro
SUBDIRS += orsaQt.staticlib.pro
