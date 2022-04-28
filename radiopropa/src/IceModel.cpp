#include <radiopropa/IceModel.h>

namespace radiopropa {

ExponentialIndex::ExponentialIndex(double n_ice, double delta_n, double z_0, double z_shift): 
	_n_ice(n_ice), 
	_delta_n(delta_n), 
	_z_shift(z_shift), 
	_z_0(z_0)
{}
ExponentialIndex::~ExponentialIndex()
{}
double ExponentialIndex::getValue(const Vector3d &position) const
{
	return _n_ice  - _delta_n  * exp((position.z-_z_shift) / _z_0);
}
double ExponentialIndex::getAverageValue(const Vector3d &position1, const Vector3d &position2) const
{
	return _n_ice - _delta_n * _z_0 / (position2.z - position1.z) * (exp((position2.z-_z_shift) / _z_0) - exp((position1.z-_z_shift) / _z_0));
}
Vector3d ExponentialIndex::getGradient(const Vector3d &position) const
{
	Vector3d v(0,0,0);
	v.z = - _delta_n / _z_0 * exp((position.z-_z_shift)/ _z_0);
	return v;
}


IceModel_Simple::IceModel_Simple(double n_ice, double delta_n, double z_0, double z_shift, double z_surface): 
	_z_surface(z_surface),
	_ice(n_ice, delta_n, z_0, z_shift)
{}
IceModel_Simple::~IceModel_Simple()
{}
double IceModel_Simple::getValue(const Vector3d &position) const
{
	if (position.z <= _z_surface) {
		return _ice.getValue(position);
	} else {
		return 1.;
	}
}
double IceModel_Simple::getAverageValue(const Vector3d &position1, const Vector3d &position2) const
{
	if ((position1.z <= _z_surface) and (position2.z <= _z_surface)) {
		return _ice.getAverageValue(position1,position2);
	} else {
		return 1.;
	}
}
Vector3d IceModel_Simple::getGradient(const Vector3d &position) const
{
	if (position.z <= _z_surface){
		return _ice.getGradient(position);
	} else {
		return Vector3d(0,0,0);
	}
}



IceModel_Firn::IceModel_Firn(double n_ice_firn, double delta_n_firn, double z_shift_firn, double z_0_firn, 
	double z_firn, double n_ice, double delta_n, double z_0, double z_shift, double z_surface) :
	_firn(n_ice_firn, delta_n_firn, z_0_firn, z_shift_firn),
	_ice(n_ice, delta_n, z_0, z_shift), 
	_z_surface(z_surface),
	_z_firn(z_firn)
{}
IceModel_Firn::~IceModel_Firn()
{}
double IceModel_Firn::getValue(const Vector3d &position) const
{
	if (position.z < _z_firn){
		return _ice.getValue(position);
	} else if (position.z <= _z_surface){
		return _firn.getValue(position);
	} else {
		return 1.;
	}
}
double IceModel_Firn::getAverageValue(const Vector3d &position1, const Vector3d &position2) const
{
	Vector3d p1 = position1;
	Vector3d p2 = position2;
	if (position1.z > position2.z){
		p1 = position2;
		p2 = position1;
	}

	if (p2.z <= _z_surface) {
		if (p2.z < _z_firn) {
			return _ice.getAverageValue(position1,position2);
		} else if (p1.z >= _z_firn) {
			return _firn.getAverageValue(position1, position2);
		} else {
			double n1 = _ice.getAverageValue(p1, Vector3d(0,0,_z_firn));
			double n2 = _firn.getAverageValue(Vector3d(0,0,_z_firn), p2);
			return (n1*(_z_firn - p1.z) + n2*(p2.z - _z_firn)) / 2;
		}
	} else if (p1.z <= _z_surface) {
		double n1 = _ice.getAverageValue(p1, Vector3d(0,0,_z_surface));
		double n2 = _firn.getAverageValue(Vector3d(0,0,_z_surface), p2);
		return (n1*(_z_surface - p1.z) + n2*(p2.z - _z_surface)) / 2;
	} else {
		return 1.;
	}
}
Vector3d IceModel_Firn::getGradient(const Vector3d &position) const
{
	if (position.z < _z_firn){
		return _ice.getGradient(position);
	} else if (position.z <= _z_surface) {
		return _firn.getGradient(position);
	} else {
		return Vector3d(0,0,0);
	}
}

/**
greenland_simple::greenland_simple(double z_surface, double n_ice, double delta_n, double z_0): IceModel_Exponential(z_surface,n_ice,delta_n,z_0)
{}
greenland_simple::~greenland_simple()
{}

southpole_simple::southpole_simple(double z_surface, double n_ice, double delta_n, double z_0): IceModel_Exponential(z_surface,n_ice,delta_n,z_0)
{}
southpole_simple::~southpole_simple()
{}

southpole_2015::southpole_2015(double z_surface, double n_ice, double delta_n, double z_0): IceModel_Exponential(z_surface,n_ice,delta_n,z_0)
{}
southpole_2015::~southpole_2015()
{}


ARAsim_southpole::ARAsim_southpole(double z_surface, double n_ice, double delta_n, double z_0): IceModel_Exponential(z_surface,n_ice,delta_n,z_0)
{}
ARAsim_southpole::~ARAsim_southpole()
{}


mooresbay_simple::mooresbay_simple(double z_surface, double n_ice, double delta_n, double z_0): IceModel_Exponential(z_surface,n_ice,delta_n,z_0)
{}
mooresbay_simple::~mooresbay_simple()
{}
**/



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