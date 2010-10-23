TEMPLATE = subdirs

unix:!macx { SUBDIRS += orsaSPICE.dynamiclib.pro }
SUBDIRS += orsaSPICE.staticlib.pro
