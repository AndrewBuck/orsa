CONFIG += ordered osg_src

include(../orsa.pri)

TEMPLATE = subdirs

!isEmpty(OSG_SRC) {
	SUBDIRS += orsaEssentialOSG
} 

SUBDIRS += orsa orsaSolarSystem orsaSPICE orsaInputOutput orsaOSG orsaQt orsaUtil
