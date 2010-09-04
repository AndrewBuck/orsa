TEMPLATE = subdirs

!win32:!macx SUBDIRS += orsaOSG.dynamiclib.pro
SUBDIRS += orsaOSG.staticlib.pro
