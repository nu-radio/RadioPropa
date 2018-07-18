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
	virtual double getValue(const Vector3d &position) const {return 1.;};
	virtual Vector3d getGradient(const Vector3d &position) const {return Vector3d(1,0,0);};
};


class LinearIncrease: public ScalarField {
	private:
		Vector3d g0;
		double v0;
public:
	LinearIncrease(double _v0, const Vector3d &_g0);
	virtual ~LinearIncrease();
	virtual double getValue(const Vector3d &position) const;
	virtual Vector3d getGradient(const Vector3d &position) const;
};

/**
    Gorham Ice Model from https://icecube.wisc.edu/~mnewcomb/radio/
    n = 1.325 + 0.463 * (1.0 - math.exp(-0.0140*depth) )
**/
class GorhamIceModel: public ScalarField
{
		private:
		double a,b,c;
		double z0;
	public:
		GorhamIceModel(double z0 = 0, double _a = 1.325, double _b = 0.463, double _c =-0.0140);
		virtual ~GorhamIceModel();
		virtual double getValue(const Vector3d &position) const;
		virtual Vector3d getGradient(const Vector3d &position) const;
};


class n2linear : public ScalarField
{
	/* Refractivity resulting in linear decrease of the velocity in z starting
	 * from z=0
	 *
Tracing Analytic Ray Curves for Light and
Sound Propagation in Non-linear Media
Qi Mo* , Hengchin Yeh* , and Dinesh Manocha*
	 * Analytic solution from
	 *
	 */
	private:
		double n0, a;

public:
	n2linear(double _n0, double _a);

	virtual ~n2linear();
	virtual double getValue(const Vector3d &position) const;
	virtual Vector3d getGradient(const Vector3d &position) const;
};



} // namespace radiopropa

#endif // CRPROPA_MAGNETICFIELD_H
