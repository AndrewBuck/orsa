#include <orsa/unit.h>

#include "UnitComboBox.h"

UnitComboBox::UnitComboBox(UnitType nUnitType, string defaultUnitString, QWidget *nParent)
	: QComboBox(nParent)
{
	unitType = nUnitType;

	switch(unitType)
	{
		case TimeIntervalUnit:
			addItem("Second", orsa::Unit::SECOND);
			addItem("Minute", orsa::Unit::MINUTE);
			addItem("Hour", orsa::Unit::HOUR);
			addItem("Day", orsa::Unit::DAY);
			addItem("Year", orsa::Unit::YEAR);

			if(defaultUnitString.length() > 0)
			{
				if(defaultUnitString.compare("second") == 0)
					setCurrentIndex(0);
				else if(defaultUnitString.compare("minute") == 0)
					setCurrentIndex(1);
				else if(defaultUnitString.compare("hour") == 0)
					setCurrentIndex(2);
				else if(defaultUnitString.compare("day") == 0)
					setCurrentIndex(3);
				else if(defaultUnitString.compare("year") == 0)
					setCurrentIndex(4);
			}
			break;
	}
}

double UnitComboBox::convertToCurrentUnit(double num)
{
	orsa::Unit::TimeUnit timeUnit;

	switch(unitType)
	{
		case TimeIntervalUnit:
			timeUnit = (orsa::Unit::TimeUnit)(itemData(currentIndex()).toInt());
			return orsa::Unit::FromUnits(num, timeUnit);
			break;
	}
}

