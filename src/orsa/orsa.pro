TEMPLATE = subdirs

unix:!macx { SUBDIRS += orsa.dynamiclib.pro }
SUBDIRS += orsa.staticlib.pro
