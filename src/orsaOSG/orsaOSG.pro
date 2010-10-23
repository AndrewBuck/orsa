TEMPLATE = subdirs

unix:!macx { SUBDIRS += orsaOSG.dynamiclib.pro }
SUBDIRS += orsaOSG.staticlib.pro
