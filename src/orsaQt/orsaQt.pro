TEMPLATE = subdirs

unix:!macx { SUBDIRS += orsaQt.dynamiclib.pro }
SUBDIRS += orsaQt.staticlib.pro
