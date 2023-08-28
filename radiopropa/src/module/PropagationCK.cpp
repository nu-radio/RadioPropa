#include "radiopropa/module/PropagationCK.h"

#include <functional> 
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <radiopropa/Trace.h>
#include <radiopropa/IceModel.h>
#include <fstream>
#include <math.h>

#include <chrono>


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
		double minStep, double maxStep, char Birefringence):
		minStep(0) {
	setField(field);
	setTolerance(tolerance);
	setMaximumStep(maxStep);
	setMinimumStep(minStep);
    setBirefringenceState(Birefringence);
    



	// load Cash-Karp coefficients
	a.assign(cash_karp_a, cash_karp_a + 36);
	b.assign(cash_karp_b, cash_karp_b + 6);
	bs.assign(cash_karp_bs, cash_karp_bs + 6);
}

void PropagationCK::process(Candidate *candidate) const {

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
	//actual step can be smaller due to non c velocity in medium
	Vector3d actual_step = yOut.x - yIn.x;
	candidate->setCurrentStep(actual_step.getR());
	candidate->setNextStep(newStep);
	double n = field->getValue(yOut.x);
	candidate->setPropagationTime(candidate->getPropagationTime() + step / c_light);
	candidate->appendPathPosition(yOut.x);

    if (getBirefringenceState() != '0')
        {
        ElectricField birefringenceField = current.getElectricField();
        Vector3d birefringenceDirection = current.getDirection();
        Vector3d birefringencePosition = current.getPosition();
        Vector3d birefringencePreviousPosition = candidate->previous.getPosition();
        Vector3d birefringenceNindex = BirefringenceIceModel().getValue(birefringencePosition, getBirefringenceState());
        birefringenceField = apply_birefringence(birefringenceField, birefringenceDirection, birefringenceNindex, birefringencePosition, birefringencePreviousPosition);
        candidate->current.setElectricField(birefringenceField);
        }
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

void PropagationCK::setBirefringenceState(char bir) {
	Birefringence_ = bir;
}

char PropagationCK::getBirefringenceState() const {
	return Birefringence_;
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


std::vector<double> PropagationCK::getEffectiveIndices_analytical(Vector3d dir, Vector3d n_vec) const 
{
    // calculated the effective refractive indices from the propagation direction and the dielectric tensor
    double N1 = sqrt((-2 * pow(n_vec.x, 2) * pow(n_vec.y, 2) * pow(n_vec.z, 2)) / (pow(n_vec.y, 2) * pow(n_vec.z, 2) * (-1 + pow(dir.x, 2)) + pow(n_vec.x, 2) * (pow(n_vec.z, 2) * (-1 + pow(dir.y, 2)) + pow(n_vec.y, 2) * (-1 + pow(dir.z, 2))) - sqrt(4 * pow(n_vec.x, 2) * pow(n_vec.y, 2) * pow(n_vec.z, 2) * (pow(n_vec.z, 2) * (-1 + pow(dir.x, 2) + pow(dir.y, 2)) + pow(n_vec.y, 2) * (-1 + pow(dir.x, 2) + pow(dir.z, 2)) + pow(n_vec.x, 2) * (-1 + pow(dir.z, 2) + pow(dir.y, 2))) + pow(pow(n_vec.y, 2) * pow(n_vec.z, 2) * (-1 + pow(dir.x, 2)) + pow(n_vec.x, 2) * (pow(n_vec.z, 2) * (-1 + pow(dir.y, 2)) + pow(n_vec.y, 2) * (-1 + pow(dir.z, 2))),2 ))));
    double N2 = sqrt((-2 * pow(n_vec.x, 2) * pow(n_vec.y, 2) * pow(n_vec.z, 2)) / (pow(n_vec.y, 2) * pow(n_vec.z, 2) * (-1 + pow(dir.x, 2)) + pow(n_vec.x, 2) * (pow(n_vec.z, 2) * (-1 + pow(dir.y, 2)) + pow(n_vec.y, 2) * (-1 + pow(dir.z, 2))) + sqrt(4 * pow(n_vec.x, 2) * pow(n_vec.y, 2) * pow(n_vec.z, 2) * (pow(n_vec.z, 2) * (-1 + pow(dir.x, 2) + pow(dir.y, 2)) + pow(n_vec.y, 2) * (-1 + pow(dir.x, 2) + pow(dir.z, 2)) + pow(n_vec.x, 2) * (-1 + pow(dir.z, 2) + pow(dir.y, 2))) + pow(pow(n_vec.y, 2) * pow(n_vec.z, 2) * (-1 + pow(dir.x, 2)) + pow(n_vec.x, 2) * (pow(n_vec.z, 2) * (-1 + pow(dir.y, 2)) + pow(n_vec.y, 2) * (-1 + pow(dir.z, 2))),2 ))));
    std::vector<double> N_vector(2);
    N_vector[0] = N1;
    N_vector[1] = N2;
    return N_vector;
}



double PropagationCK::getTimeDelay(std::vector<double> N_vector, double l) const
{
    //calculates the time delay from the effective refractive indices and the propagation length
    double c = 0.299792458;
    double time_delay = l/c * (N_vector[0] - N_vector[1]);
    return time_delay;
}

Vector3d PropagationCK::getPolarization(double n, Vector3d dir, Vector3d n_vec, double prec) const
{
    //calculates the polarization vector in spherical coordinates from the effective refractive indices, the propagation direction and the dielectric tensor
    Vector3d e_v;
    Vector3d e_c;
    double p = 1 / prec;
    double n_round = ceil(n * p) / p;

    if (ceil(n_vec.x * p) / p == n_round)
        {
        e_v.x = 1.0;
        e_v.y = 0.0;
        e_v.z = 0.0;
        }

    else if (ceil(n_vec.y * p) / p == n_round)
        {
        e_v.x = 0.0;
        e_v.y = 1.0;
        e_v.z = 0.0;
        }

    else if (ceil(n_vec.z * p) / p == n_round)
        {
        e_v.x = 0.0;
        e_v.y = 0.0;
        e_v.z = 1.0;
        }

    else
        {
        e_v.x = dir.x/(pow(n, 2) - pow(n_vec.x, 2));
        e_v.y = dir.y/(pow(n, 2) - pow(n_vec.y, 2));
        e_v.z = dir.z/(pow(n, 2) - pow(n_vec.z, 2));

        double norm = sqrt(pow(e_v.x, 2)+pow(e_v.y, 2)+pow(e_v.z, 2));
        e_v.x = e_v.x / norm;
        e_v.y = e_v.y / norm;
        e_v.z = e_v.z / norm;
        }

    double theta = acos(dir.z);
    double phi =  atan2(dir.y, dir.x);

    e_c.x = sin(theta) * cos(phi) * e_v.x + sin(theta) * sin(phi) * e_v.y + cos(theta) * e_v.z;
    e_c.y = cos(theta) * cos(phi) * e_v.x + cos(theta) * sin(phi) * e_v.y - sin(theta) * e_v.z;
    e_c.z = -sin(phi) * e_v.x + cos(phi) * e_v.y;

    return e_c;
}


ElectricField PropagationCK::apply_birefringence(ElectricField Pulse, Vector3d dir, Vector3d n_vec, Vector3d cur_pos, Vector3d pre_pos) const
{
    // applies the birefringence effect to the given pulse using the propagation direction and the dielectric tensor

    std::vector<std::vector<std::complex<double>>> Efield = Pulse.getFrequencySpectrum();

    double rate = Pulse.getSamplingRate();
    double prec = 1e-08;
    double l = sqrt(pow(pre_pos.x - cur_pos.x, 2) + pow(pre_pos.y - cur_pos.y, 2) + pow(pre_pos.z - cur_pos.z, 2));    

    std::vector<double> N_eff = getEffectiveIndices_analytical(dir, n_vec);
    Vector3d pol_1 = getPolarization(N_eff[0], dir, n_vec, prec);
    Vector3d pol_2 = getPolarization(N_eff[1], dir, n_vec, prec);

    double dt;
  
    if (isnan(pol_1.y) or isnan(pol_1.z) or isnan(pol_2.y) or isnan(pol_2.z)) 
    {
    return Pulse;
    }

    double a;
    double b;
    double c;
    double d;

    double a_inv;
    double b_inv;
    double c_inv;
    double d_inv;

    double det = pol_1.y * pol_2.z - pol_1.z * pol_2.y;

    if (det == 0)

        {
        a = 1;
        b = 0;
        c = 0;
        d = 1;

        dt = 0;

        a_inv = 1;
        b_inv = 0;
        c_inv = 0;
        d_inv = 1;
        }

    else 
        {
        a = pol_1.y;
        b = pol_1.z;
        c = pol_2.y;
        d = pol_2.z;

        dt = getTimeDelay(N_eff, l);

        a_inv = 1 / det * d;
        b_inv = -1 / det * b;
        c_inv = -1 / det * c;
        d_inv = 1 / det * a;
        }

    std::vector<std::complex<double>> Th_a(Efield[1].begin(), Efield[1].end());
    std::vector<std::complex<double>> Ph_b(Efield[1].begin(), Efield[1].end());
    std::vector<std::complex<double>> Th_c(Efield[1].begin(), Efield[1].end());
    std::vector<std::complex<double>> Ph_d(Efield[1].begin(), Efield[1].end());

    std::transform(Efield[1].begin(), Efield[1].end(), Th_a.begin(), [&a](std::complex<double> element) { return element *= a; });
    std::transform(Efield[2].begin(), Efield[2].end(), Ph_b.begin(), [&b](std::complex<double> element) { return element *= b; });
    std::transform(Efield[1].begin(), Efield[1].end(), Th_c.begin(), [&c](std::complex<double> element) { return element *= c; });
    std::transform(Efield[2].begin(), Efield[2].end(), Ph_d.begin(), [&d](std::complex<double> element) { return element *= d; });

    std::transform(Th_a.begin(), Th_a.end(), Ph_b.begin(), Efield[1].begin(), std::plus<std::complex<double>>());
    std::transform(Th_c.begin(), Th_c.end(), Ph_d.begin(), Efield[2].begin(), std::plus<std::complex<double>>());

    Trace N_shift;
    N_shift.setFrequencySpectrum(Efield[2], rate);
    N_shift.applyTimeShift(dt);
    Efield[2] = N_shift.getFrequencySpectrum();

    std::transform(Efield[1].begin(), Efield[1].end(), Th_a.begin(), [&a_inv](std::complex<double> element) { return element *= a_inv; });
    std::transform(Efield[2].begin(), Efield[2].end(), Ph_b.begin(), [&b_inv](std::complex<double> element) { return element *= b_inv; });
    std::transform(Efield[1].begin(), Efield[1].end(), Th_c.begin(), [&c_inv](std::complex<double> element) { return element *= c_inv; });
    std::transform(Efield[2].begin(), Efield[2].end(), Ph_d.begin(), [&d_inv](std::complex<double> element) { return element *= d_inv; });

    std::transform(Th_a.begin(), Th_a.end(), Ph_b.begin(), Efield[1].begin(), std::plus<std::complex<double>>());
    std::transform(Th_c.begin(), Th_c.end(), Ph_d.begin(), Efield[2].begin(), std::plus<std::complex<double>>());

    ElectricField Pulse_new;
    Pulse_new.setFrequencySpectrum(Efield[0], Efield[1], Efield[2], rate);    

    return Pulse_new;
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
