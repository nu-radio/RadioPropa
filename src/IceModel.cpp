#include <radiopropa/IceModel.h>

namespace radiopropa {

IceModel_Exponential::IceModel_Exponential(double z_surface, double n_ice, double delta_n, double z_0): _z_surface(z_surface), _n_ice(n_ice), _delta_n(delta_n), _z_0(z_0)
{}

IceModel_Exponential::~IceModel_Exponential()
{}

double IceModel_Exponential::getValue(const Vector3d &position) const
{
	if (position.z - _z_surface <0)
	return _n_ice  - _delta_n  * exp(position.z / _z_0);
	else
	return 1.;
}

Vector3d IceModel_Exponential::getGradient(const Vector3d &position) const
{
	Vector3d v(0,0,0);
	if (position.z <0)
	{
	v.z = - _delta_n / _z_0 * exp(position.z/ _z_0);
	}
	return v;
}

/*
IceModel_GreenlandSimple::IceModel_GreenlandSimple(double z_surface, double n_ice, double delta_n, double z_0): IceModel_Exponential(z_surface,n_ice,delta_n,z_0)
{}
IceModel_GreenlandSimple::~IceModel_GreenlandSimple()
{}

IceModel_SouthPoleSimple::IceModel_SouthPoleSimple(double z_surface, double n_ice, double delta_n, double z_0): IceModel_Exponential(z_surface,n_ice,delta_n,z_0)
{}
IceModel_SouthPoleSimple::~IceModel_SouthPoleSimple()
{}

IceModel_SouthPole2015::IceModel_SouthPole2015(double z_surface, double n_ice, double delta_n, double z_0): IceModel_Exponential(z_surface,n_ice,delta_n,z_0)
{}
IceModel_SouthPole2015::~IceModel_SouthPole2015()
{}


IceModel_SouthPoleARAsim::IceModel_SouthPoleARAsim(double z_surface, double n_ice, double delta_n, double z_0): IceModel_Exponential(z_surface,n_ice,delta_n,z_0)
{}
IceModel_SouthPoleARAsim::~IceModel_SouthPoleARAsim()
{}


IceModel_MooresbaySimple::IceModel_MooresbaySimple(double z_surface, double n_ice, double delta_n, double z_0): IceModel_Exponential(z_surface,n_ice,delta_n,z_0)
{}
IceModel_MooresbaySimple::~IceModel_MooresbaySimple()
{}

*/


GorhamIceModel::GorhamIceModel(double z_surface, double a, double b, double c) : _z_surface(z_surface), _a(a), _b(b), _c(c)
{

}

GorhamIceModel::~GorhamIceModel()
{
}

double GorhamIceModel::getValue(const Vector3d &position) const
{
	if (position.z-_z_surface < 0)
      return _a + _b * (1.0 - exp(-1.*_c*(position.z-_z_surface)));
	else
      return 1.;

}

Vector3d GorhamIceModel::getGradient(const Vector3d &position) const
{
	Vector3d v(0,0,0);
	if (position.z < 0)
	{
      v.z = 1.0 * _b * _c * exp(-1.*_c*position.z);
	}
	// The gradient on discontinuities has to be 0 as these are handeled by a
	// different module!
	return v;
}
}