TEMPLATE = subdirs

unix:!macx { SUBDIRS += orsaInputOutput.dynamiclib.pro }
SUBDIRS += orsaInputOutput.staticlib.pro
