#ifndef RADIOPROPA_ICEMODEL_H
#define RADIOPROPA_ICEMODEL_H

#include "radiopropa/Units.h"
#include "radiopropa/Vector3.h"
#include "radiopropa/ScalarField.h"

#include <iostream>
#include <fstream>


namespace radiopropa {

class ExponentialIndex: public ScalarField
{
	protected:
		double _n_ice, _delta_n, _z_shift, _z_0;
	public:
		ExponentialIndex(double n_ice, double delta_n, double z_0, double z_shift = 0);
		virtual ~ExponentialIndex();
		virtual double getValue(const Vector3d &position) const;
		virtual double getAverageValue(const Vector3d &position1, const Vector3d &position2) const; 
		virtual Vector3d getGradient(const Vector3d &position) const;
};

class IceModel_Simple: public ScalarField
{
	protected:
		ExponentialIndex _ice;
		double _z_surface;
	public:
		IceModel_Simple(double n_ice, double delta_n, double z_0, double z_shift = 0, double z_surface = 0);
		virtual ~IceModel_Simple();
		virtual double getValue(const Vector3d &position) const;
		virtual double getAverageValue(const Vector3d &position1, const Vector3d &position2) const; 
		virtual Vector3d getGradient(const Vector3d &position) const;
};

class IceModel_Firn: public ScalarField
{
	protected:
		double _z_surface, _z_firn;
		ExponentialIndex _firn;
		ExponentialIndex _ice;
	public:
		IceModel_Firn(  
			double n_ice_firn, 
			double delta_n_firn, 
			double z_shift_firn, 
			double z_0_firn,
			double z_firn,
			double n_ice,  
			double delta_n, 
			double z_0,
			double z_shift = 0,
			double z_surface = 0);
		virtual ~IceModel_Firn();
		virtual double getValue(const Vector3d &position) const; 
		virtual double getAverageValue(const Vector3d &position1, const Vector3d &position2) const; 
		virtual Vector3d getGradient(const Vector3d &position) const;
};

/**
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
**/

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





/*class BirefringenceIceModel: public VectorField */
class BirefringenceIceModel: public Referenced
{
	private:

		double _n_ice,_delta_n,_z_0;

	public:
		BirefringenceIceModel(double n_ice = 1.78, double delta_n = 0.423, double z_0 =77);
		virtual ~BirefringenceIceModel();
		virtual Vector3d getValue(Vector3d &position);

        virtual double BSpline(double x, std::vector<double> t, const std::vector<double> c, const int k);
        virtual double Xindex(double x);
        virtual double Yindex(double x);
        virtual double Zindex(double x);
};



}
















#endif //RADIOPROPA_ICEMODEL_H
