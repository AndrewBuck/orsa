TEMPLATE = subdirs

!win32: {
SUBDIRS += orsaInputOutput.dynamiclib.pro
}
SUBDIRS += orsaInputOutput.staticlib.pro
