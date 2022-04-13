#include <radiopropa/IceModel.h>

#include <vector>

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
	if ((position1.z <= _z_surface) and (position2.z <= _z_surface)) {
		if ((position1.z < _z_firn) and (position2.z < _z_firn)) {
			return _ice.getAverageValue(position1,position2);
		} else if ((position1.z >= _z_firn) and (position2.z >= _z_firn)) {
			return _firn.getAverageValue(position1,position2);
		} else {
			double n1 = _ice.getAverageValue(position1,position2);
			double n2 = _firn.getAverageValue(position1,position2);
			return (n1+n2)/2;
		}
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









BirefringenceIceModel::BirefringenceIceModel(double n_ice, double delta_n, double z_0) : _n_ice(n_ice), _delta_n(delta_n), _z_0(z_0)
{

}

BirefringenceIceModel::~BirefringenceIceModel()
{
}






  // returns the value of a BSpline function at position x for internal knots t, coefficients c and order k
  double BirefringenceIceModel::BSpline(double x, std::vector<double> t, const std::vector<double> c, const int k) 
  {
    for (int i = 0; i < k; ++i) {  // add outer knots
      t.insert(t.begin(), t[0]);
      t.push_back(t[t.size() - 1]);
    }
    const int m = t.size();
    double B[m - 1][k + 1];


    for (int i = 0; i < (m - 1); ++i) {
      if ((x >= t[i]) && (x < t[i + 1])) {
        B[i][0] = 1.;
      } else {
        B[i][0] = 0.;
      }
    }

    for (int j = 1; j < (k + 1); ++j) {
      for (int i = 0; i < (m - j - 1); ++i) {
        double first_term = 0;
        double second_term = 0;
        if (t[i + j] - t[i] == 0.0) {
          first_term = 0.0;
        } else {
          first_term = ((x - t[i]) / (t[i + j] - t[i])) * B[i][j - 1];
        }
        if (t[i + j + 1] - t[i + 1] == 0.0) {
          second_term = 0.0;
        } else {
          second_term = ((t[i + j + 1] - x) / (t[i + j + 1] - t[i + 1])) * B[i + 1][j - 1];
        }
        B[i][j] = first_term + second_term;
      }
    }
    double y = 0;
    for (int i = 0; i < (m - k - 1); ++i) {
      y += c[i] * B[i][k];
    }

    return y;
  }


  double BirefringenceIceModel::Xindex(double x) 
  {
    static const double t[] = {0.     ,      73.01301301 ,  91.23123123 , 109.44944945 , 114.07407407,
  118.55855856 , 127.66766767 , 136.77677678 , 138.03803804  ,139.15915916,
  160.      ,    320.   ,       480.     ,     636.   ,       797.,
  957.     ,    2135.63563564, 2500.      };
    std::vector<double> knots(t, t + sizeof(t) / sizeof(t[0]));
    static const double c[] = { 1.77905532, 1.77905889 ,1.7790641,  1.77906385 ,1.77906251, 1.77906099,
 1.7790598 , 1.77905798, 1.77905602, 1.77905472, 1.77905203, 1.7790262,
 1.77890482 ,1.77861256, 1.77828408, 1.77793062, 1.77676224, 1.77700331,
 1.77692707, 1.77694786   };
    std::vector<double> coeffs(c, c + sizeof(c) / sizeof(c[0]));
    return BSpline(x, knots, coeffs, 3);
  }


  double BirefringenceIceModel::Yindex(double x) 
  {
    static const double t[] = { 0.      ,     73.01301301 ,  91.23123123  ,109.44944945 , 118.55855856,
  127.66766767,  132.29229229 , 136.77677678  ,320.    ,      480.,
  636.      ,    797. ,         957. ,        1816.81681682, 1862.36236236,
 1953.45345345, 2135.63563564, 2500.       };
    std::vector<double> knots(t, t + sizeof(t) / sizeof(t[0]));
    static const double c[] = { 1.77976748, 1.77977201, 1.77976748, 1.77976479 ,1.77976486, 1.77976603,
 1.77976731, 1.7797687,  1.77978596, 1.77991725, 1.78037522, 1.7806736,
 1.78077385, 1.78119828, 1.78126133, 1.78125755, 1.78125539, 1.7812571,
 1.78125605, 1.78125644 };
    std::vector<double> coeffs(c, c + sizeof(c) / sizeof(c[0]));
    return BSpline(x, knots, coeffs, 3);
  }

  double BirefringenceIceModel::Zindex(double x) 
  {
    static const double t[] = {  0.    ,       73.01301301 , 109.44944945, 127.66766767 , 136.77677678,
  320.     ,     957.   ,      1953.45345345, 2135.63563564, 2500.         };
    std::vector<double> knots(t, t + sizeof(t) / sizeof(t[0]));
    static const double c[] = { 1.78111645, 1.78111532, 1.78111053, 1.78111663, 1.78111896 ,1.78111715,
 1.78075855, 1.78177168, 1.78174358, 1.78173987, 1.78174342, 1.78174101 };
    std::vector<double> coeffs(c, c + sizeof(c) / sizeof(c[0]));
    return BSpline(x, knots, coeffs, 3);
  }




Vector3d BirefringenceIceModel::getValue(Vector3d &position) 
{

      
       double n = _n_ice - _delta_n * (exp(position.z / _z_0));
       double nx = n + Xindex( - position.z)  - 1.78 ;
       double ny = n + Yindex( - position.z)  - 1.78;
       double nz = n + Zindex( - position.z)  - 1.78;

   
    Vector3d nv (nx, ny, nz);

       return nv;
        


}








}
