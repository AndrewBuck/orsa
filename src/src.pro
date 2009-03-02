CONFIG += ordered osg_src

include(../orsa.pri)

TEMPLATE = subdirs

!isEmpty(OSG_SRC) {
	SUBDIRS += orsaEssentialOSG
} 

SUBDIRS += orsaTBB orsa orsaSolarSystem orsaInputOutput orsaOSG orsaQt orsaSPICE
