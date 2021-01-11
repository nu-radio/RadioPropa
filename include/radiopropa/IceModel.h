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

/*
class IceModel_GreenlandSimple: public IceModel_Exponential
{
	public:
		IceModel_GreenlandSimple(double z_surface = 0, double n_ice = 1.78 , double delta_n = 0.51, double z_0 = 37.25*meter);
		virtual ~IceModel_GreenlandSimple();
};
class IceModel_SouthPoleSimple: public IceModel_Exponential
{
	public:
		IceModel_SouthPoleSimple(double z_surface = 0, double n_ice = 1.78 , double delta_n = 0.426, double z_0 = 71*meter);
		virtual ~IceModel_SouthPoleSimple();
};
class IceModel_SouthPole2015: public IceModel_Exponential
{
	public:
		IceModel_SouthPole2015(double z_surface = 0, double n_ice = 1.78 , double delta_n = 0.423, double z_0 = 77*meter);
		virtual ~IceModel_SouthPole2015();
};
class IceModel_SouthPoleARAsim: public IceModel_Exponential
{
	public:
		IceModel_SouthPoleARAsim(double z_surface = 0, double n_ice = 1.78 , double delta_n = 0.43, double z_0 = 75.75757575757576*meter);
		virtual ~IceModel_SouthPoleARAsim();
};
class IceModel_MooresbaySimple: public IceModel_Exponential
{
	public:
		IceModel_MooresbaySimple(double z_surface = 0, double n_ice = 1.78 , double delta_n = 0.46, double z_0 = 37.25*meter);
		virtual ~IceModel_MooresbaySimple();
};
*/

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