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



N_constant::N_constant(double _z0, double _n): z0(_z0), n(_n)
{

}

double N_constant::getValue(const Vector3d &position) const
{
	if (position.z - z0 > 0)
		return n;
	else
		return 1.0;
}

Vector3d N_constant::getGradient(const Vector3d &position) const
{
	return Vector3d(0,0,0);
}



Lin_grad::Lin_grad(double _z0, double _step_n): z0(_z0), step_n(_step_n)
{

}

double Lin_grad::getValue(const Vector3d &position) const
{
	double d = position.z - z0;
	if (d > 0)
                return 1.0 + d*step_n;
        else
                return 1.0;

}

Vector3d Lin_grad::getGradient(const Vector3d &position) const
{
       if (position.z - z0 > 0)
	       return Vector3d(0, 0, step_n);
       else
	       return Vector3d(0,0,0);
}



CloudModel_atm::CloudModel_atm(double _z0, double _T0): z0(_z0), T0(_T0)
{

}

double CloudModel_atm::L = 6.5e-3;
double CloudModel_atm::p0 = 870; //mbar Bishop is in the mountains!
double CloudModel_atm::M = 0.02896; //kg/mol
double CloudModel_atm::R = 8.314; //J/K/mol
double CloudModel_atm::g = 9.807; //m/s**2
double CloudModel_atm::D = 0.61121;
double CloudModel_atm::a = 18.678;
double CloudModel_atm::b = 234.5;
double CloudModel_atm::c = 257.14;
double CloudModel_atm::e = 0.78;



double CloudModel_atm::getValue(const Vector3d &position) const
{

	double h = position.z;
        double T = T0 - L*h;
	double T_C = T - 273.15;
        double P_air = p0*exp(-M*g*h/R/T);
        double P_sat = D*exp((a-T_C/b)*(T_C/(c+T_C)))*10; //Bucks equation gives kPa

	if (position.z - z0 > 0)
	{
		double N= 77.6*P_air/T + 3.73*1e5*P_sat/pow(T, 2);
		return N*1e-6+1.0;
	}
	else
	{
		double N= 77.6*P_air/T + 3.73*1e5*e*P_sat/pow(T, 2);
		return N*1e-6+1.0;
	}

}

Vector3d CloudModel_atm::getGradient(const Vector3d &position) const
{
	double h = position.z;
        double T = T0 - L*h;
	double T_C = T - 273.15;
        double P_air = p0*exp(-M*g*h/R/T);
        double P_sat = D*exp((a-T_C/b)*(T_C/(c+T_C)))*10;

	double dP_air = P_air *(-M*g/R)*T0/pow(T,2);
	double dP_sat = -L* P_sat*T_C/(T_C+c)*(-2/b+a/T_C-(a-T_C/b)/(T_C+c));
        double dT = -L;

	if (position.z - z0 > 0)
	{
		double dN= 77.6*(1/T*dP_air - P_air/pow(T,2)*dT) + 3.73*1e5*(1/pow(T,2)*dP_sat-2*P_sat/pow(T,3)*dT);
		double dn =  dN*1e-6;
		return Vector3d(0,0,dn);
	}
	else
	{
		double dN= 77.6*(1/T*dP_air - P_air/pow(T,2)*dT) + 3.73*1e5*(1/pow(T,2)*e*dP_sat-2*e*P_sat/pow(T,3)*dT);
                double dn =  dN*1e-6;
		return Vector3d(0,0,dn);
	}
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
