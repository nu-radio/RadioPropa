#ifndef CRPROPA_ICEMODEL_H
#define CRPROPA_ICEMODEL_H

#include "radiopropa/Units.h"
#include "radiopropa/Vector3.h"
#include "radiopropa/ScalarField.h"

#endif

namespace radiopropa {

class IceModel_Exponential: public ScalarField
{
	protected:
		double _n_ice, _delta_n, _z_0;
		double _z_surface;
	public:
		IceModel_Exponential(double z_surface = 0, double n_ice = 1 , double delta_n = 1, double z_0 = 1);
		virtual ~IceModel_Exponential();
		virtual double getValue(const Vector3d &position) const; 
		virtual Vector3d getGradient(const Vector3d &position) const;
};


class greenland_simple: public IceModel_Exponential
{
	public:
		greenland_simple(double z_surface = 0, double n_ice = 1.78 , double delta_n = 0.51, double z_0 = 37.25*meter);
		virtual ~greenland_simple();
};
class southpole_simple: public IceModel_Exponential
{
	public:
		southpole_simple(double z_surface = 0, double n_ice = 1.78 , double delta_n = 0.426, double z_0 = 71*meter);
		virtual ~southpole_simple();
};
class southpole_2015: public IceModel_Exponential
{
	public:
		southpole_2015(double z_surface = 0, double n_ice = 1.78 , double delta_n = 0.423, double z_0 = 77*meter);
		virtual ~southpole_2015();
};
class ARAsim_southpole: public IceModel_Exponential
{
	public:
		ARAsim_southpole(double z_surface = 0, double n_ice = 1.78 , double delta_n = 0.43, double z_0 = 75.75757575757576*meter);
		virtual ~ARAsim_southpole();
};
class mooresbay_simple: public IceModel_Exponential
{
	public:
		mooresbay_simple(double z_surface = 0, double n_ice = 1.78 , double delta_n = 0.46, double z_0 = 37.25*meter);
		virtual ~mooresbay_simple();
};


/**
    Gorham Ice Model from https://icecube.wisc.edu/~mnewcomb/radio/
    n = 1.325 + 0.463 * (1.0 - math.exp(-0.0140*depth) )
**/
class GorhamIceModel: public ScalarField
{
	private:
		double _a,_b,_c;
		double _z_surface;
	public:
		GorhamIceModel(double z_surface = 0, double a = 1.325, double b = 0.463, double c =-0.0140);
		virtual ~GorhamIceModel();
		virtual double getValue(const Vector3d &position) const;
		virtual Vector3d getGradient(const Vector3d &position) const;
};

}