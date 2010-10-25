#ifndef UNITCOMBOBOX_H
#define UNITCOMBOBOX_H

#include <string>
using namespace std;

#include <QtGui/QComboBox>

class UnitComboBox : public QComboBox
{
	Q_OBJECT

	public:
		enum UnitType{TimeIntervalUnit=1, LengthUnit, MassUnit};

		UnitComboBox(UnitType nUnitType, string defaultUnitString = "", QWidget *nParent = NULL);
		double convertToCurrentUnit(double num);

	private:
		UnitType unitType;
};

#endif

