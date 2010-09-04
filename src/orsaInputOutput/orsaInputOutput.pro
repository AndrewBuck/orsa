TEMPLATE = subdirs

!win32:!macx SUBDIRS += orsaInputOutput.dynamiclib.pro
SUBDIRS += orsaInputOutput.staticlib.pro
