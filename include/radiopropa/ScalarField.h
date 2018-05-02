#ifndef CRPROPA_SCALARFIELD_H
#define CRPROPA_SCALARFIELD_H

#include "radiopropa/Units.h"
#include "radiopropa/Vector3.h"
#include "radiopropa/Referenced.h"

#ifdef CRPROPA_HAVE_MUPARSER
#include "muParser.h"
#endif

namespace radiopropa {

/**
 @class ScalarField
 @brief Abstract base class fora scalar field
 */
class ScalarField: public Referenced {
public:
	virtual ~ScalarField() {
	}
	virtual double getValue(const Vector3d &position) const {
		return 1;
	};
	virtual Vector3d getGradient(const Vector3d &position) const {
		return Vector3d(0,0,0);
	};
};


class LinearIncrease: public ScalarField{
	private:
		Vector3d g0;
		double v0;

public:
	LinearIncrease(double _v0, const Vector3d &_g0) : v0(_v0), g0(_g0)
	{

	}

	virtual ~LinearIncrease() {
	}
	virtual double getValue(const Vector3d &position) const {
		double z = position.dot(g0);
		if ( position.z <= 0)
			return 1.;
		else
			return v0 * position.z;
	};
	virtual Vector3d getGradient(const Vector3d &position) const {
		return g0;
	};
};


class Gorham: public ScalarField
{
    
	/*
    Gorham Ice Model from https://icecube.wisc.edu/~mnewcomb/radio/
    n = 1.325 + 0.463 * (1.0 - math.exp(-0.0140*depth) )
   */
	private:
		double a,b,c;
	public:
		Gorham() : a(1.325), b(0.463), c(-0.0140)
		{
	
		}
		virtual ~Gorham()
		{
		}

	virtual double getValue(const Vector3d &position) const {
        return a + b * (1.0 - exp(-1.*c*position.z));
	};

	virtual Vector3d getGradient(const Vector3d &position) const {
        Vector3d v(0,0,0);
        v.z = -1.0 * b * c * exp(-1.*c*position.z);
        return v;
	}
};





} // namespace radiopropa

#endif // CRPROPA_MAGNETICFIELD_H
