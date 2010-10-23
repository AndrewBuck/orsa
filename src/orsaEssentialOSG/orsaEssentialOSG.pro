TEMPLATE = subdirs

unix:!macx { SUBDIRS += orsaEssentialOSG.dynamiclib.pro }
SUBDIRS += orsaEssentialOSG.staticlib.pro
