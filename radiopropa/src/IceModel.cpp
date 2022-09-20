#include <radiopropa/IceModel.h>

namespace radiopropa {

ExponentialIndex::ExponentialIndex()
{}
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
	Vector3d p1 = position1;
	Vector3d p2 = position2;
	if (position1.z > position2.z){
		p1 = position2;
		p2 = position1;
	}

	if (p2.z <= _z_surface) {
		return _ice.getAverageValue(position1,position2);
	} else if (p1.z <= _z_surface) {
		double n1 = _ice.getAverageValue(p1, Vector3d(0,0,_z_surface));
		double n2 = 1.;
		return (n1*(_z_surface - p1.z) + n2*(p2.z - _z_surface)) / (p2.z - p1.z);
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
		} else if (p1.z < _z_firn) {
			double n1 = _ice.getAverageValue(p1, Vector3d(0,0,_z_firn));
			double n2 = _firn.getAverageValue(Vector3d(0,0,_z_firn), p2);
			return (n1*(_z_firn - p1.z) + n2*(p2.z - _z_firn)) / (p2.z - p1.z);
		} else {
			return _firn.getAverageValue(position1, position2);
		}
	} else if (p1.z < _z_firn){
		double n1 = _ice.getAverageValue(p1, Vector3d(0,0,_z_firn));
		double n2 = _firn.getAverageValue(Vector3d(0,0,_z_firn), Vector3d(0,0,_z_surface));
		double n3 = 1;
		return (n1*(_z_firn - p1.z) + n2*(_z_surface - _z_firn) + n3*(p2.z - _z_surface)) / (p2.z - p1.z);
	} else if (p1.z <= _z_surface) {
		double n1 = _ice.getAverageValue(p1, Vector3d(0,0,_z_surface));
		double n2 = _firn.getAverageValue(Vector3d(0,0,_z_surface), p2);
		return (n1*(_z_surface - p1.z) + n2*(p2.z - _z_surface)) / (p2.z - p1.z);
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

IceModel_DoubleFirn::IceModel_DoubleFirn(ExponentialIndex firn_top, ExponentialIndex firn_bottom,
	ExponentialIndex ice, double z_firn_top, double z_firn_bottom, double z_surface)
{
	_z_surface = z_surface;
	_z_firn_top = z_firn_top;
	_z_firn_bottom = z_firn_bottom;
	_firn_top = firn_top;
	_firn_bottom = firn_bottom;
	_ice = ice;
}
IceModel_DoubleFirn::~IceModel_DoubleFirn()
{}
double IceModel_DoubleFirn::getValue(const Vector3d &position) const
{
	if (position.z < _z_firn_bottom){
		return _ice.getValue(position);
	} else if (position.z < _z_firn_top){
		return _firn_bottom.getValue(position);
	} else if (position.z <= _z_surface){
		return _firn_top.getValue(position);
	} else {
		return 1.;
	}
}
double IceModel_DoubleFirn::getAverageValue(const Vector3d &position1, const Vector3d &position2) const
{
	Vector3d p1 = position1;
	Vector3d p2 = position2;
	if (position1.z > position2.z){
		p1 = position2;
		p2 = position1;
	}

	if (p2.z < _z_firn_bottom) {
		return _ice.getAverageValue(position1,position2);
	} else if (p1.z > _z_surface) {
		return 1;
	} else if (p2.z < _z_firn_top) {
		if (p1.z >= _z_firn_bottom) {
			return _firn_bottom.getAverageValue(position1,position2);
		} else {
			double n1 = _ice.getAverageValue(p1, Vector3d(0,0,_z_firn_bottom));
			double n2 = _firn_bottom.getAverageValue(Vector3d(0,0,_z_firn_bottom), p2);
			return (n1*(_z_firn_bottom - p1.z) + n2*(p2.z - _z_firn_bottom)) / (p2.z - p1.z);	
		}
	} else if (p2.z <= _z_surface) {
		if (p1.z >= _z_firn_top) {
			return _firn_top.getAverageValue(position1,position2);
		} else if (p1.z >= _z_firn_bottom) {
			double n1 = _firn_bottom.getAverageValue(p1, Vector3d(0,0,_z_firn_top));
			double n2 = _firn_top.getAverageValue(Vector3d(0,0,_z_firn_top), p2);
			return (n1*(_z_firn_top - p1.z) + n2*(p2.z - _z_firn_top)) / (p2.z - p1.z);	
		} else {
			double n1 = _ice.getAverageValue(p1, Vector3d(0,0,_z_firn_bottom));
			double n2 = _firn_bottom.getAverageValue(Vector3d(0,0,_z_firn_bottom), Vector3d(0,0,_z_firn_top));
			double n3 = _firn_top.getAverageValue(Vector3d(0,0,_z_firn_top), p2);
			return (n1*(_z_firn_bottom - p1.z) + n2*(_z_firn_top - _z_firn_bottom) + n3*(p2.z - _z_firn_top)) / (p2.z - p1.z);
		}
	} else {
		if (p1.z >= _z_firn_top) {
			double n1 = _firn_top.getAverageValue(p1, Vector3d(0,0,_z_surface));
			double n2 = 1;
			return (n1*(_z_surface - p1.z) + n2*(p2.z - _z_surface)) / (p2.z - p1.z);
		} else if (p1.z >= _z_firn_bottom) {
			double n1 = _firn_bottom.getAverageValue(p1, Vector3d(0,0,_z_firn_top));
			double n2 = _firn_top.getAverageValue(Vector3d(0,0,_z_firn_top), Vector3d(0,0,_z_surface));
			double n3 = 1;
			return (n1*(_z_firn_top - p1.z) + n2*(_z_surface - _z_firn_top) + n3*(p2.z - _z_surface)) / (p2.z - p1.z);
		} else {
			double n1 = _ice.getAverageValue(p1, Vector3d(0,0,_z_firn_bottom));
			double n2 = _firn_bottom.getAverageValue(Vector3d(0,0,_z_firn_bottom), Vector3d(0,0,_z_firn_top));
			double n3 = _firn_top.getAverageValue(Vector3d(0,0,_z_firn_top), Vector3d(0,0,_z_surface));
			double n4 = 1;
			return (n1*(_z_firn_bottom - p1.z) + n2*(_z_firn_top - _z_firn_bottom) + n3*(_z_surface - _z_firn_top) + n4*(p2.z - _z_surface)) / (p2.z - p1.z);
		}
	}
}
Vector3d IceModel_DoubleFirn::getGradient(const Vector3d &position) const
{
	if (position.z < _z_firn_bottom){
		return _ice.getGradient(position);
	} else if (position.z < _z_firn_top){
		return _firn_bottom.getGradient(position);
	} else if (position.z <= _z_surface) {
		return _firn_top.getGradient(position);
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