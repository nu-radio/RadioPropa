#include <radiopropa/ScalarField.h>

namespace radiopropa {

LinearIncrease::LinearIncrease(double _v0, const Vector3d &_g0) : v0(_v0), g0(_g0)
{

}

LinearIncrease::~LinearIncrease() {
}

double LinearIncrease::getValue(const Vector3d &position) const
{
    return v0 * position.z;
};

Vector3d LinearIncrease::getGradient(const Vector3d &position) const
{
  return Vector3d(0, 0, v0);
};



GorhamIceModel::GorhamIceModel(double z0, double _a, double _b, double _c) : z0(z0), a(_a), b(_b), c(_c)
{

}

GorhamIceModel::~GorhamIceModel()
{
}

double GorhamIceModel::getValue(const Vector3d &position) const
{
	if (position.z-z0 < 0)
      return a + b * (1.0 - exp(-1.*c*(position.z-z0)));
	else
      return 1.;

};

Vector3d GorhamIceModel::getGradient(const Vector3d &position) const
{
	Vector3d v(0,0,0);
	if (position.z < 0)
	{
      v.z = 1.0 * b * c * exp(-1.*c*position.z);
	}
	// The gradient on discontinuities has to be 0 as these are handeled by a
	// different module!
	return v;
}



n2linear::n2linear(double _n0, double _a) : n0(_n0), a(_a) { }

n2linear::~n2linear() { };

double n2linear::getValue(const Vector3d &position) const
{
    return sqrt(n0*n0 + a * position.z);
}

Vector3d n2linear::getGradient(const Vector3d &position) const
{
		return Vector3d(0,0,a/2. / getValue(position));
}


} // namespace
