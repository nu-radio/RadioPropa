#include "radiopropa/module/PropagationCK.h"

#include <functional> 
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <radiopropa/Trace.h>


namespace radiopropa {

// Cash-Karp coefficients
const double cash_karp_a[] = { 0., 0., 0., 0., 0., 0., 1. / 5., 0., 0., 0., 0.,
		0., 3. / 40., 9. / 40., 0., 0., 0., 0., 3. / 10., -9. / 10., 6. / 5.,
		0., 0., 0., -11. / 54., 5. / 2., -70. / 27., 35. / 27., 0., 0., 1631.
				/ 55296., 175. / 512., 575. / 13824., 44275. / 110592., 253.
				/ 4096., 0. };

const double cash_karp_b[] = { 37. / 378., 0, 250. / 621., 125. / 594., 0., 512.
		/ 1771. };

const double cash_karp_bs[] = { 2825. / 27648., 0., 18575. / 48384., 13525.
		/ 55296., 277. / 14336., 1. / 4. };

void PropagationCK::tryStep(const Y &y, Y &out, Y &error, double h,
		ParticleState &particle, double z) const {
	std::vector<Y> k;
	k.reserve(6);

	out = y;
	error = Y(0);

	// calculate the sum of b_i * k_i
	for (size_t i = 0; i < 6; i++) {

		Y y_n = y;
		for (size_t j = 0; j < i; j++)
			y_n += k[j] * a[i * 6 + j] * h;

		// update k_i
		k[i] = dYdt(y_n, particle, z);

		out += k[i] * b[i] * h;
		error += k[i] * (b[i] - bs[i]) * h;
	}
}


PropagationCK::Y PropagationCK::dYdt(const Y &y, ParticleState &p, double z) const {
	// normalize direction vector to prevent numerical losses
	double n = field->getValue(y.x);
	Vector3d velocity = y.u.getUnitVector() * c_light / n;
	Vector3d dudt =  field->getGradient(y.x) / n/n * c_light; 
	return Y(velocity, dudt);
}

PropagationCK::PropagationCK(ref_ptr<ScalarField> field, double tolerance,
		double minStep, double maxStep) :
		minStep(0) {
	setField(field);
	setTolerance(tolerance);
	setMaximumStep(maxStep);
	setMinimumStep(minStep);

	// load Cash-Karp coefficients
	a.assign(cash_karp_a, cash_karp_a + 36);
	b.assign(cash_karp_b, cash_karp_b + 6);
	bs.assign(cash_karp_bs, cash_karp_bs + 6);
}

void PropagationCK::process(Candidate *candidate) const {
	// save the new previous particle state
	ParticleState &current = candidate->current;
	candidate->previous = current;

	double step = clip(candidate->getNextStep(), minStep, maxStep);

	Y yIn(current.getPosition(), current.getDirection());
	Y yOut, yErr;
	double newStep = step;
	double r = 42;  // arbitrary value > 1
	double z = 0; // RedShift to 0.

	// try performing step until the target error (tolerance) or the minimum step size has been reached
	while (r > 1) {
		step = newStep;
		tryStep(yIn, yOut, yErr, step / c_light, current, z);

		r = yErr.u.getR() / tolerance;  // ratio of absolute direction error and tolerance
		newStep = step * 0.95 * pow(r, -0.2);  // update step size to keep error close to tolerance
		newStep = clip(newStep, 0.1 * step, 5 * step);  // limit the step size change
		newStep = clip(newStep, minStep, maxStep);

		if (step == minStep)
			break;  // performed step already at the minimum
	}

	current.setPosition(yOut.x);
	current.setDirection(yOut.u.getUnitVector());
	//actual step can be smaller  du to non c velocity in medium
	Vector3d actual_step = yOut.x - yIn.x;
	candidate->setCurrentStep(actual_step.getR());
	candidate->setNextStep(newStep);
	double n = field->getValue(yOut.x);
	candidate->setPropagationTime(candidate->getPropagationTime() + step / c_light);
	candidate->appendPathPosition(yOut.x);
}

void PropagationCK::setField(ref_ptr<ScalarField> f) {
	field = f;
}

void PropagationCK::setTolerance(double tol) {
	if ((tol > 1) or (tol < 0))
		throw std::runtime_error(
				"PropagationCK: target error not in range 0-1");
	tolerance = tol;
}

void PropagationCK::setMinimumStep(double min) {
	if (min < 0)
		throw std::runtime_error("PropagationCK: minStep < 0 ");
	if (min > maxStep)
		throw std::runtime_error("PropagationCK: minStep > maxStep");
	minStep = min;
}

void PropagationCK::setMaximumStep(double max) {
	if (max < minStep)
		throw std::runtime_error("PropagationCK: maxStep < minStep");
	maxStep = max;
}

double PropagationCK::getTolerance() const {
	return tolerance;
}

double PropagationCK::getMinimumStep() const {
	return minStep;
}

double PropagationCK::getMaximumStep() const {
	return maxStep;
}


double PropagationCK::getdeterminant(double n, Vector3d dir, Vector3d n_vec) 
{

    double A = (pow(n_vec.x, 2) - pow(n, 2))*(pow(n_vec.y, 2) - pow(n, 2))*(pow(n_vec.z, 2) - pow(n, 2));
    double B = pow(n, 2) * ( pow(dir.x, 2) *(pow(n_vec.y, 2) - pow(n, 2))*(pow(n_vec.z, 2) - pow(n, 2)) + pow(dir.y, 2) *(pow(n_vec.x, 2) - pow(n, 2))*(pow(n_vec.z, 2) - pow(n, 2)) + pow(dir.z, 2) *(pow(n_vec.x, 2) - pow(n, 2))*(pow(n_vec.y, 2) - pow(n, 2)));

       return  A + B ;
}


double PropagationCK::bisection(double a, double b, Vector3d dir, Vector3d n_vec)
{

    double EPSILON = 0.0000000001;
    double c = a;
    while ((b-a) >= EPSILON)
    {
        c = (a+b)/2;

        if (getdeterminant(c, dir,  n_vec) == 0.0)
            break;

        else if (getdeterminant(c, dir, n_vec)*getdeterminant(a,  dir,  n_vec) < 0)
            b = c;
        else
            a = c;
    }
    return (a+b)/2;
}


double PropagationCK::minFinder(double a, double b, Vector3d dir, Vector3d n_vec)
{

    double EPSILON = 0.0000000001;
    double c = a;
    while ((b-a) >= EPSILON)
    {
        c = (a+b)/2;

        if (getdeterminant(a,  dir,  n_vec) > getdeterminant(b,  dir,  n_vec))
            a = c;

        else if (getdeterminant(a,  dir,  n_vec) < getdeterminant(b,  dir,  n_vec))
            b = c;
        else if (getdeterminant(a,  dir,  n_vec) == getdeterminant(b,  dir,  n_vec))
            break;
    }
    return (a+b)/2;
}




Vector3d PropagationCK::getEffectiveIndices(Vector3d dir, Vector3d n_vec) {

    /*
    dir.x = 1 / sqrt(3);    
    dir.y = 1 / sqrt(3); 
    dir.z = 1 / sqrt(3); 

    n_vec.x = 1.77;    
    n_vec.y = 1.78; 
    n_vec.z = 1.79; 
    */


    double minimum = minFinder(1, 2,  dir,  n_vec);
    double N1 = bisection(1, minimum,  dir,  n_vec);
    double N2 = bisection(minimum, 2,  dir,  n_vec);

    /*
    std::cout.precision (15);
    std::cout << N1<<std::endl;
    std::cout << N2<<std::endl;
    */

    Vector3d N (N1, N2, 1);

    return N;
}



double PropagationCK::getTimeDelay(double n1, double n2, double l)
{
    double c = 299792458;
    return l/c * (n1 - n2);
}

Vector3d PropagationCK::getPolarization(double n, Vector3d dir, Vector3d n_vec)
{
    Vector3d e_v;
    e_v.x = dir.x/(pow(n, 2) - pow(n_vec.x, 2));
    e_v.y = dir.y/(pow(n, 2) - pow(n_vec.y, 2));
    e_v.z = dir.z/(pow(n, 2) - pow(n_vec.z, 2));
    return e_v;
}


ElectricField PropagationCK::apply_birefringence(ElectricField Pulse, Vector3d dir, Vector3d n_vec)
//std::vector<std::vector<double>> PropagationCK::apply_birefringence(ElectricField Pulse, Vector3d dir, Vector3d n_vec)
{

    std::vector<std::vector<double>> Efield = Pulse.getFrequencySpectrum_real();

    double rate = Pulse.getSamplingRate();

    float l = 1;     /*has to be set but how? */

    Vector3d N_eff = getEffectiveIndices(dir, n_vec);
    double dt = getTimeDelay(N_eff.x, N_eff.y, l);
    Vector3d pol_1 = getPolarization(N_eff.x, dir, n_vec);
    Vector3d pol_2 = getPolarization(N_eff.y, dir, n_vec);

    double a = pol_1.y;
    double b = pol_1.z;
    double c = pol_2.y;
    double d = pol_2.z;

    //Trace N1 = Efield[1] * a + Efield[2] * b;
    //Trace N2 = Efield[1] * c + Efield[2] * d;

    //std::vector<double> X = {1.1, 2.2, 3.3, 4.4};
    //std::vector<double> Y = {10, 20, 30, 40};

    //std::vector<double> pass = X;
    //double scale = 2;

    std::vector<double> Th_a = Efield[1];
    std::vector<double> Ph_b = Efield[2];
    std::vector<double> Th_c = Efield[1];
    std::vector<double> Ph_d = Efield[2];

    std::vector<double> N1 = Efield[1];
    std::vector<double> N2 = Efield[1];


    std::transform(Efield[1].begin(), Efield[1].end(), Th_a.begin(), [&a](double element) { return element *= a; });
    std::transform(Efield[2].begin(), Efield[2].end(), Ph_b.begin(), [&b](double element) { return element *= b; });
    std::transform(Efield[1].begin(), Efield[1].end(), Th_c.begin(), [&c](double element) { return element *= c; });
    std::transform(Efield[2].begin(), Efield[2].end(), Ph_d.begin(), [&d](double element) { return element *= d; });

    std::transform(Th_a.begin(), Th_a.end(), Ph_b.begin(), N1.begin(), std::plus<double>());
    std::transform(Th_c.begin(), Th_c.end(), Ph_d.begin(), N2.begin(), std::plus<double>());

    Trace N_shift;
    N_shift.setFrequencySpectrum(N1, N1, rate);
    
    N_shift.applyTimeShift(dt);

    N1 = N_shift.getFrequencySpectrum_real();

    double det = a*d - b*c;

    double a_inv = 1 / det * d;
    double b_inv = -1 / det * b;
    double c_inv = -1 / det * c;
    double d_inv = 1 / det * a;


    std::vector<double> N1_a = N1;
    std::vector<double> N2_b = N2;
    std::vector<double> N1_c = N1;
    std::vector<double> N2_d = N2;

    //std::vector<double> N1 = Efield[1];
    //std::vector<double> N2 = Efield[1];

    std::transform(N1.begin(), N1.end(), N1_a.begin(), [&a_inv](double element) { return element *= a_inv; });
    std::transform(N2.begin(), N2.end(), N2_b.begin(), [&b_inv](double element) { return element *= b_inv; });
    std::transform(N1.begin(), N1.end(), N1_c.begin(), [&c_inv](double element) { return element *= c_inv; });
    std::transform(N2.begin(), N2.end(), N2_d.begin(), [&d_inv](double element) { return element *= d_inv; });

    std::transform(N1_a.begin(), N1_a.end(), N2_b.begin(), Efield[1].begin(), std::plus<double>());
    std::transform(N1_c.begin(), N1_c.end(), N2_d.begin(), Efield[2].begin(), std::plus<double>());


    //Trace Th = N1 * a_inv + N2 * b_inv;
    //Trace Ph = N1 * c_inv + N2 * d_inv;


    Pulse.setFrequencySpectrum(Efield[0], Efield[0], Efield[1], Efield[1], Efield[2], Efield[2], rate);

    return Pulse;
}


std::vector<double> PropagationCK::apply_birefringence_1(Vector3d dir, Vector3d n_vec)
{

    Trace Theta;
    Trace Phi;
    std::vector<double> real{0,0,0,0,0,1,0,0,0,0,0};
    std::vector<double> imag{0,0,0,0,0,2,0,0,0,0,0};
    std::vector<double> realr{0,0,0,0,0,1,0,0,0,0,0};
    std::vector<double> imagr{0,0,0,0,0,1,0,0,0,0,0};

    Theta.setFrequencySpectrum(real, real, 0.001);
    Phi.setFrequencySpectrum(imag, imag, 0.001);

    Theta.addTraces(Phi);

    Theta.multiplyConstant(5);

    return Theta.getFrequencySpectrum_real();
}





std::string PropagationCK::getDescription() const {
	std::stringstream s;
	s << "Propagation in magnetic fields using the Cash-Karp method.";
	s << " Target error: " << tolerance;
	s << ", Minimum Step: " << minStep / kpc << " kpc";
	s << ", Maximum Step: " << maxStep / kpc << " kpc";
	return s.str();
}

} // namespace radiopropa
