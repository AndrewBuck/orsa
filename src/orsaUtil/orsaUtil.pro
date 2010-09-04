TEMPLATE = subdirs

!win32:!macx SUBDIRS += orsaUtil.dynamiclib.pro
SUBDIRS += orsaUtil.staticlib.pro
