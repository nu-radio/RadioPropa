#ifndef RADIOPROPA_ICEMODEL_EXP3_H
#define RADIOPROPA_ICEMODEL_EXP3_H

#include "radiopropa/Vector3.h"
#include "radiopropa/ScalarField.h"
#include <cmath>          // For std::exp



namespace radiopropa {

class IceModel_Exp3 : public ScalarField 
{
	protected:
		double _cvac, _z1, _z2;
		double gradient_snow(double z) const;
        	double gradient_firn(double z) const;
        	double gradient_bubbly(double z) const;
	public:
		IceModel_Exp3();
		virtual ~IceModel_Exp3();
		virtual double getValue(const Vector3d &position) const;
		virtual double getAverageValue(const Vector3d &position1, const Vector3d &position2) const;
		virtual Vector3d getGradient(const Vector3d &position) const;
};
}
#endif //RADIOPROPA_ICEMODEL_EXP3_H
