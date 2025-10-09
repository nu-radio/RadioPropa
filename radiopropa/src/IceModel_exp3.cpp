#include <radiopropa/IceModel_exp3.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

namespace radiopropa {

IceModel_Exp3::IceModel_Exp3()
  	: _cvac(1),
	_z1(-14.9 / _cvac),
    	_z2(-80.5 / _cvac)
{
}
IceModel_Exp3::~IceModel_Exp3()
{}
double IceModel_Exp3::getValue(const Vector3d &position) const 
{
	if (position.z > 0) {
		return 1.0;
	}
    	else if (position.z > _z1) {
        	return 1.51188 - 0.271579 * std::exp(0.114553 * position.z * _cvac);
    	}
    	else if (position.z > _z2) {
        	return 1.89957 - 0.529715 * std::exp(0.0129175 * position.z * _cvac);
    	}
    	else {
        	return 1.77468 - 1.41573 * std::exp(0.0387882 * position.z * _cvac);
    	}
}
double IceModel_Exp3::getAverageValue(const Vector3d &position1, const Vector3d &position2) const
{
	Vector3d p1 = position1;
	Vector3d p2 = position2;
	if (position1.z > position2.z){
		p1 = position2;
		p2 = position1;
	}
	double n1 = getValue(position1);
	double n2 = getValue(position2);
	double n_z2 = getValue(Vector3d(0,0,_z2));
	double n_z1 = getValue(Vector3d(0,0,_z1));
	
	if (p2.z <= 0) {
		if (p2.z <= _z2) {
			return (n1 + n2) / 2.0;
		} else if (p1.z <= _z2) {
			if (p2.z > _z2 && p2.z <= _z1) {
				double avg1 = (n1 + n_z2)/2.0;	
				double avg2 = (n_z2 + n2)/2.0;
				return (avg1*(_z2 - p1.z) + avg2*(p2.z - _z2)) / (p2.z - p1.z);
			} else {
				double avg1 = (n1 + n_z2)/2.0;
			  	double avg2 = (n_z2 + n_z1)/2.0;
				double avg3 = (n_z1 + n2)/2.0;
				return (avg1*(_z2 - p1.z) + avg2*(_z1 - _z2) + avg3*(p2.z - _z1)) / (p2.z - p1.z);
			}
		} else if (p1.z > _z2 && p1.z <= _z1) {
			if (p2.z > _z2 && p2.z <= _z1) {
				return (n1 + n2) / 2.0;
			} else {
				double avg1 = (n1 + n_z1)/2.0;
				double avg2 = (n_z1 + n2)/2.0;
				return (avg1*(_z1 - p1.z) + avg2*(p2.z - _z1)) / (p2.z - p1.z);
			}
		} else {
		       return (n1 + n2) / 2.0;
		}	
	} else if (p1.z <= _z2) {
		double avg1 = (n1 + n_z2)/2.0;
		double avg2 = (n_z2 + n_z1)/2.0;
		double avg3 = (n_z1 + 1.0)/2.0;
		double avg4 = (1.0 + n2)/2.0;
		return (avg1*(_z2 - p1.z) + avg2*(_z1 - _z2) + avg3*(0 - _z1) + avg4*(p2.z))/(p2.z - p1.z);
	} else if (p1.z > _z2 && p1.z <= _z1) {
		double avg1 = (n1 + n_z1)/2.0;
		double avg2 = (n_z1 + 1.0)/2.0;
		double avg3 = (1.0 + n2)/2.0;
		return (avg1*(_z1 - p1.z) + avg2*(0 - _z1) + avg3*(p2.z))/(p2.z - p1.z);
	} else if (p1.z > _z1 && p1.z <= 0) {
		double avg1 = (n1 + 1.0)/2.0;
		double avg2 = (1.0 + n2)/2.0;
		return (avg1*(0 - p1.z) + avg2*(p2.z))/(p2.z - p1.z);
	} else {
		return 1.0;
	}
}

double IceModel_Exp3::gradient_snow(double z) const {
    return -0.271579 * 0.114553 * _cvac * std::exp(0.114553 * z * _cvac);
}

double IceModel_Exp3::gradient_firn(double z) const {
    return -0.529715 * 0.0129175 * _cvac * std::exp(0.0129175 * z * _cvac);
}

double IceModel_Exp3::gradient_bubbly(double z) const {
    return -1.41573 * 0.0387882 * _cvac * std::exp(0.0387882 * z * _cvac);
}

Vector3d IceModel_Exp3::getGradient(const Vector3d &position) const
{
	if (position.z <= _z2) {
		return Vector3d(0,0,gradient_bubbly(position.z));
	} else if (position.z > _z2 && position.z <= _z1) {
		return Vector3d(0,0,gradient_firn(position.z));
	} else if (position.z > _z1 && position.z <= 0) {
		return Vector3d(0,0,gradient_snow(position.z));
	} else {
		return Vector3d(0,0,0);
	}
}
}


