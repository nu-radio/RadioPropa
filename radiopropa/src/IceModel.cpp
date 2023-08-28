#include <radiopropa/IceModel.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

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




IceModel_Data1D::IceModel_Data1D(std::string interpolation, int axis):
interpolation(interpolation), axis(axis)
{}
IceModel_Data1D::~IceModel_Data1D()
{}
double IceModel_Data1D::getValue(const Vector3d &position) const
{
	std::vector<double> pos = {position.x, position.y, position.z};
	double n;

	if ((coordinates.front() > pos[axis]) or (coordinates.back() < pos[axis])) {
		throw std::domain_error("The position must lie inside the range of the data.");
	}
	if (interpolation == "linear"){
		for (int i=0; i < coordinates.size(); i++){
			double x0 = coordinates[i];
			double x1 = coordinates[i+1];
			double y0 = indices_of_refraction[i];
			double y1 = indices_of_refraction[i+1];

			if ((x0 <= pos[axis]) and (x1 >= pos[axis])) {
				n = y0 + ((y1-y0)/(x1-x0)) * (pos[axis] - x0);
				continue;
			}
		}
	}
	return n;
} 
double IceModel_Data1D::getAverageValue(const Vector3d &position1, const Vector3d &position2) const
{
	double n_avg;
	std::vector<double> pos1 = {position1.x, position1.y, position1.z};
	std::vector<double> pos2 = {position2.x, position2.y, position2.z};
	if (pos1[axis] > pos2[axis]){
		pos2 = {position1.x, position1.y, position1.z};
		pos1 = {position2.x, position2.y, position2.z};
	}

	if ((coordinates.front() > pos1[axis]) or (coordinates.back() < pos2[axis])) {
		throw std::domain_error("The positions must lie inside the range of the data.");
	}

	if (interpolation == "linear"){
		double xp1 = pos1[axis];
		double yp1 = 0;
		double xp2 = pos2[axis];
		double yp2 = 0;

		double integral = 0;

		for (int i=0; i < coordinates.size(); i++){
			double x0 = coordinates[i];
			double x1 = coordinates[i+1];

			if ((xp1 > x1) or (xp2 < x0)) {
				continue;
			}

			double y0 = indices_of_refraction[i];
			double y1 = indices_of_refraction[i+1];

			if ((x0 <= xp1) and (x1 >= xp1)) {
				yp1 = y0 + ((y1-y0)/(x1-x0)) * (xp1 - x0);
				integral += ((x1 - xp1) * (y1 + yp1) / 2);
			}
			if ((x0 > xp1) and (x0 < xp2)){
				integral += ((x1 - x0) * (y1 + y0) / 2);
			}
			if ((x0 <= xp2) and (x1 >= xp2)) {
				yp2 = y0 + ((y1-y0)/(x1-x0)) * (xp2 - x0);
				integral -= ((x1 - xp2) * (y1 + yp2) / 2);
			}
		}

		n_avg = (integral / (xp2 - xp1));
	}
	return n_avg;
}
Vector3d IceModel_Data1D::getGradient(const Vector3d &position) const
{
	Vector3d gradient = Vector3d(0,0,0);
	std::vector<double> pos = {position.x, position.y, position.z};

	if ((coordinates.front() >= pos[axis]) or (coordinates.back() <= pos[axis])) {
		throw std::domain_error("The position must lie inside the range of the data.");
	}

	if (interpolation == "linear"){
		double g = 0;
		for (int i=0; i < coordinates.size(); i++){
			if ((coordinates[i] < pos[axis]) and (coordinates[i+1] > pos[axis])){
				g = (indices_of_refraction[i+1]-indices_of_refraction[i])/(coordinates[i+1]-coordinates[i]);
				continue;
			} else if (coordinates[i] == pos[axis]) {
				g  = (indices_of_refraction[i+1]-indices_of_refraction[i])/(coordinates[i+1]-coordinates[i]);
				g += (indices_of_refraction[i]-indices_of_refraction[i-1])/(coordinates[i]-coordinates[i-1]);
				g /= 2;
				continue;
			}
		}

		if (axis == 0) {
			gradient.x = g;
		} else if (axis == 1) {
			gradient.y = g;
		} else if (axis == 2) {
			gradient.z = g;
		}
	}
	return gradient;
}

void IceModel_Data1D::loadDataFromCSV(std::string filename, char delimeter, int header_lineindex)
{
	struct data_point {
	  double coord;
	  double ior;
	  static bool compare(data_point dp1, data_point dp2)
	  {
			return (dp1.coord < dp2.coord);
		}
	};

	std::ifstream data_file;
	data_file.open(filename);
	if (!data_file.is_open())
	{
		throw std::ios_base::failure("File failed to open.");
	}

	std::vector<data_point> data;
	std::string line;
	std::string temp_coord;
	std::string temp_ior;
	int lineindex = 0;
	while (std::getline(data_file, line)) {
		if ((line.empty()) or (line.front() == '#')){
			lineindex ++;
			continue;
		}

		if (lineindex == header_lineindex){
			lineindex ++;
			continue;
		}

		std::stringstream ss(line);
		std::getline(ss,temp_coord, delimeter);
		std::getline(ss,temp_ior,delimeter);
		data.push_back({std::atof(temp_coord.c_str()),std::atof(temp_ior.c_str())});
	}

	std::sort(data.begin(), data.end(), data_point::compare);
	coordinates.clear();
	indices_of_refraction.clear();
	for (int i=0; i < data.size(); i++){
		coordinates.push_back(data[i].coord);
		indices_of_refraction.push_back(data[i].ior);
	}
}

void IceModel_Data1D::loadDataFromVector(std::vector<double> coord, std::vector<double> ior)
{
	struct data_point {
	  double coord;
	  double ior;
	  static bool compare(data_point dp1, data_point dp2)
	  {
			return (dp1.coord < dp2.coord);
		}
	};
		
	std::vector<data_point> data;

	for (int i=0; i < coord.size(); i++){
		data.push_back({coord[i],ior[i]});
	}

	std::sort(data.begin(), data.end(), data_point::compare);

	coordinates.clear();
	indices_of_refraction.clear();
	for (int i=0; i < data.size(); i++){
		coordinates.push_back(data[i].coord);
		indices_of_refraction.push_back(data[i].ior);
	}
}

void IceModel_Data1D::loadDataFromVector(std::vector<std::vector<double>> data)
{
	struct data_point {
	  static bool compare(std::vector<double> dp1, std::vector<double> dp2)
	  {
			return (dp1[0] < dp2[0]);
		}
	};
		
	std::sort(data.begin(), data.end(), data_point::compare);

	coordinates.clear();
	indices_of_refraction.clear();
	for (int i=0; i < data.size(); i++){
		coordinates.push_back(data[i][0]);
		indices_of_refraction.push_back(data[i][1]);
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


  double BirefringenceIceModel::A_Xindex(double x) 
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


  double BirefringenceIceModel::A_Yindex(double x) 
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

  double BirefringenceIceModel::A_Zindex(double x) 
  {
    static const double t[] = {  0.    ,       73.01301301 , 109.44944945, 127.66766767 , 136.77677678,
  320.     ,     957.   ,      1953.45345345, 2135.63563564, 2500.         };
    std::vector<double> knots(t, t + sizeof(t) / sizeof(t[0]));
    static const double c[] = { 1.78111645, 1.78111532, 1.78111053, 1.78111663, 1.78111896 ,1.78111715,
 1.78075855, 1.78177168, 1.78174358, 1.78173987, 1.78174342, 1.78174101 };
    std::vector<double> coeffs(c, c + sizeof(c) / sizeof(c[0]));
    return BSpline(x, knots, coeffs, 3);
  }




Vector3d BirefringenceIceModel::getValue(Vector3d &position, char bire_model) 
{

	double n = _n_ice - _delta_n * (exp(position.z / _z_0));

	double nx = n + A_Xindex( - position.z)  - 1.78;
	double ny = n + A_Yindex( - position.z)  - 1.78;
	double nz = n + A_Zindex( - position.z)  - 1.78;


	Vector3d nv (nx, ny, nz);

       return nv;
        


}








}
