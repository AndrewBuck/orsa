TEMPLATE = subdirs

unix:!macx { SUBDIRS += orsaSolarSystem.dynamiclib.pro }
SUBDIRS += orsaSolarSystem.staticlib.pro
