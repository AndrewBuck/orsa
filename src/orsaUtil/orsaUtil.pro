TEMPLATE = subdirs

unix:!macx { SUBDIRS += orsaUtil.dynamiclib.pro }
SUBDIRS += orsaUtil.staticlib.pro
