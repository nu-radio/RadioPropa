#ifndef CRPROPA_SCALARFIELD_H
#define CRPROPA_SCALARFIELD_H

#include "radiopropa/Units.h"
#include "radiopropa/Vector3.h"
#include "radiopropa/Referenced.h"

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
 * @class N_constant
 * @brief constant refractive index
 */
class N_constant : public ScalarField {
	private:
		double n;
		double z0;
	public:
		N_constant(double _z0 = 0, double _n = 1.5);
		virtual double getValue(const Vector3d &position) const;
		virtual Vector3d getGradient(const Vector3d &position) const;
};

/**
 * @class Lin_grad
 * @brief linear refractive index change with offset
 */
class Lin_grad : public ScalarField {
	private:
		double step_n;
		double z0;
	public:
		Lin_grad(double _z0 = 2000, double _step_n = 0.1);
		virtual double getValue(const Vector3d &position) const;
		virtual Vector3d getGradient(const Vector3d &position) const;
};

/**
 * @class CloudModel_atm
 * @brief atmospheric refractive index model depending on temperature, pressure, and humidity with layers of 100%
 * humidity as cloud layers. From https://www.itu.int/rec/R-REC-P.453/en
 */
//atmospheric refractive index model, depending on temperature, pressure and humidity, saturated humidity as clouds
class CloudModel_atm : public ScalarField {
	private:
		double z_bottom;
		double z_top;
		double T0;
		double p0;
		double e;
		static double L;
		static double a;
		static double b;
		static double c;
		static double D;
		static double M;
		static double R;
		static double g;
	public:
		CloudModel_atm(double _z_bottom = 2000, double _z_top = 2500, double _T0 = 283, double _p0 = 870,
		               double _e = 0.78);
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

/**
 * @class surfaceDuct
 * @brief profile of a surface duct. The surface duct consists of a small region with a steep N decrease topped by a
 * regular atmospheric N-decrease. Specific values for Bishop used, from https://www.itu.int/rec/R-REC-P.453/en
 */
class surfaceDuct : public ScalarField {
	public:
		surfaceDuct();
		virtual double getValue(const Vector3d &position) const;
		virtual Vector3d getGradient(const Vector3d &position) const;
};

/**
 * @class elevatedDuct
 * @brief profile of an elevated duct. The elevated duct consists of a region with a steep N decrease,
 * a small region of shallow N decrease, followed  by normal atmospheric N decrease. Specific values for Bishop used,
 * from https://www.itu.int/rec/R-REC-P.453/en
 */
class elevatedDuct : public ScalarField {
	public:
		elevatedDuct();
		virtual double getValue(const Vector3d &position) const;
		virtual Vector3d getGradient(const Vector3d &position) const;
};

} // namespace radiopropa

#endif
