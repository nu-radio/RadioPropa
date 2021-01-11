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
}

Vector3d LinearIncrease::getGradient(const Vector3d &position) const
{
  return Vector3d(0, 0, v0);
}



N_constant::N_constant(double _z0, double _n) : z0(_z0), n(_n) {}

double N_constant::getValue(const Vector3d &position) const
{
	if (position.z - z0 > 0)
		return n;
	else
		return 1.0;
}

Vector3d N_constant::getGradient(const Vector3d &position) const
{
	return Vector3d(0, 0, 0);
}


Lin_grad::Lin_grad(double _z0, double _step_n) : z0(_z0), step_n(_step_n) {}

double Lin_grad::getValue(const Vector3d &position) const
{
	double d = position.z - z0;
	if (d > 0)
		return 1.0 + d * step_n;
	else
		return 1.0;
}

Vector3d Lin_grad::getGradient(const Vector3d &position) const
{
	if (position.z - z0 > 0)
		return Vector3d(0, 0, step_n);
	else
		return Vector3d(0, 0, 0);
}


CloudModel_atm::CloudModel_atm(double _z_bottom, double _z_top, double _T0, double _p0, double _e) :
		z_bottom(_z_bottom), z_top(_z_top), T0(_T0), p0(_p0), e(_e) {}

	double CloudModel_atm::L = 6.5e-3; //K/m, tropospheric lapse rate
	double CloudModel_atm::M = 0.02896; //kg/mol, molar mass of air
	double CloudModel_atm::R = 8.314; //J/K/mol, universal gas constant
	double CloudModel_atm::g = 9.807; //m/s**2, gravitational acceleration
	double CloudModel_atm::D = 6.1121; //constants for Bucks equation
	double CloudModel_atm::a = 18.678;
	double CloudModel_atm::b = 234.5;
	double CloudModel_atm::c = 257.14;

double CloudModel_atm::getValue(const Vector3d &position) const
{
	double h = position.z;
	double T = T0 - L * h; //tropospheric temperature decrease
	double T_C = T - 273.15;
	double P_air = p0 * pow(T / T0, g * M / L / R); //barometric height equation with linear temperature decrease
	double P_sat = D * exp((a - T_C / b) * (T_C / (c + T_C))); //Bucks equation

	if (z_bottom < position.z && position.z < z_top)
	{
		double N = 77.6 * P_air / T - 5.6 * P_sat / T + 3.75 * 1e5 * P_sat / pow(T, 2); //according to ITU Rec. 453
		return N * 1e-6 + 1.0;
	}
	else
	{
		double N = 77.6 * P_air / T - 5.6 * e * P_sat / T + 3.75 * 1e5 * e * P_sat / pow(T, 2);
		return N * 1e-6 + 1.0;
	}

}

Vector3d CloudModel_atm::getGradient(const Vector3d &position) const
{
	double h = position.z;
	double T = T0 - L * h;
	double T_C = T - 273.15;
	double P_air = p0 * pow(T/T0,g*M/L/R);
	double P_sat = D * exp((a - T_C / b) * (T_C / (c + T_C)));
	double dP_air = - P_air /T *M*g/R; //analytic derivations of formulas above
	double dP_sat = -L * P_sat * T_C / (T_C + c) * (-2 / b + a / T_C - (a - T_C / b) / (T_C + c));
	double dT = -L;

	if (z_bottom < position.z && position.z < z_top)
	{
		double dN = 77.6 * (1 / T * dP_air - P_air / pow(T, 2) * dT) - 5.6 *(1/T*dP_sat - 1/pow(T,2)*P_sat*dT) +
                        3.75 * 1e5 * (1 / pow(T, 2) * dP_sat - 2 * P_sat / pow(T, 3) * dT); //analytic derivation
		double dn = dN * 1e-6;
		return Vector3d(0, 0, dn);
	}
	else
	{
		double dN = 77.6 * (1 / T * dP_air - P_air / pow(T, 2) * dT) - 5.6 *(1/T*e*dP_sat - 1/pow(T,2)*e*P_sat*dT)+
                        3.75 * 1e5 * (1 / pow(T, 2) * e * dP_sat - 2 * e * P_sat / pow(T, 3) * dT);
		double dn = dN * 1e-6;
		return Vector3d(0, 0, dn);
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

surfaceDuct::surfaceDuct() {}

double surfaceDuct::getValue(const Vector3d &position) const
{
	if (position.z < 64.2) // steep N decrease
	{
		double M =  300 - 8.74/64.21*position.z;
		double N = M - 157.0/1000*position.z;
		return N*1e-6 + 1;
	}
	else // regular N decrease
	{
		double M = 300 - 8.74 + 117.0/1000*(position.z - 64.21);
		double N = M - 157.0/1000*position.z;
		return N*1e-6 + 1;
	}
}

Vector3d surfaceDuct::getGradient(const Vector3d &position) const
{
        if (position.z < 64.2) // steep N decrease
        {
            double dM = -8.74/64.21;
            double dN =  dM - 157.0/1000.0;
            return Vector3d(0,0,dN*1e-6);
        }
        else // regular N decrease
        {
            double dN = -40;
            return Vector3d(0,0,dN*1e-9);
        }
}

elevatedDuct::elevatedDuct(){}

double elevatedDuct::getValue(const Vector3d &position) const
{
	double M_e = 9.29; // constants defining slope and z value for slope change, average values for Bishop
	double u = 1166.48;
	double b = 1100.25;
	double th = 147.81;
	double dM = 117.0/1000;
	double M;
	double dm_1 = M_e/(u-b);
	double dm_2 = M_e/(b + th -u);
	if (position.z < 1166.48) // steep N decrease
	{
		M = 300 + dm_1*position.z;
	}
	else if(position.z < 1248.06) // shallow N decrease
	{
		M = 300 + dm_1*u - dm_2*(position.z - u);
	}
	else // regular N decrease
	{
		M = 300 + dm_1*u - dm_2*(b + th - u) + dM*(position.z - b - th);
	}
        double N = M - 157.0/1000*position.z;
        return N*1e-6 + 1;
}

Vector3d elevatedDuct::getGradient(const Vector3d &position) const
{
	double dm;
	if(position.z < 1166.48) // steep N decrease
	{
		dm = 9.29/(1166.48 - 1100.25);
	}
	else if(position.z < 1248.06) // shallow N decrease
	{
		dm = -9.29/(1248.06 - 1166.48);
	}
	else // regular N decrease
	{
		dm = 117.0/1000;
	}
	double dN = dm - 157.0/1000;
	return Vector3d(0,0,dN*1e-6);
}

} // namespace
