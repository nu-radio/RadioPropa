#ifndef RADIOPROPA_PROPAGATIONCK_H
#define RADIOPROPA_PROPAGATIONCK_H

#include "radiopropa/Module.h"
#include "radiopropa/Units.h"
#include "radiopropa/ScalarField.h"


namespace radiopropa {

/**
 @class PropagationCK
 @brief Propagation through scalar using the Cash-Karp method.

 This module solves the Eikonal equation of motion of a ray propagating through a refractivity field.\n
 It uses the Runge-Kutta integration method with Cash-Karp coefficients.\n
 The step size control tries to keep the relative error close to, but smaller than the designated tolerance.
 Additionally a minimum and maximum size for the steps can be set.
 */
class PropagationCK: public Module {
public:
	class Y {
	public:
		Vector3d x, u; /*< phase-point: position and direction */

		Y() {
		}

		Y(const Vector3d &x, const Vector3d &u) :
				x(x), u(u) {
		}

		Y(double f) :
				x(Vector3d(f, f, f)), u(Vector3d(f, f, f)) {
		}

		Y operator *(double f) const {
			return Y(x * f, u * f);
		}

		Y &operator +=(const Y &y) {
			x += y.x;
			u += y.u;
			return *this;
		}
	};

private:
	std::vector<double> a, b, bs; /*< Cash-Karp coefficients */
	ref_ptr<ScalarField> field;
	double tolerance; /*< target relative error of the numerical integration */
	double minStep; /*< minimum step size of the propagation */
	double maxStep; /*< maximum step size of the propagation */




public:
	PropagationCK(ref_ptr<ScalarField> field = NULL, double tolerance = 1e-4,
			double minStep = (1E-3 * meter), double maxStep = (1 * meter));
	void process(Candidate *candidate) const;

	// derivative of phase point, dY/dt = d/dt(x, u) = (v, du/dt)
	// du/dt = q*c^2/E * (u x B)
	Y dYdt(const Y &y, ParticleState &p, double z) const;

	void tryStep(const Y &y, Y &out, Y &error, double t,
			ParticleState &p, double z) const;

	void setField(ref_ptr<ScalarField> field);
	void setTolerance(double tolerance);
	void setMinimumStep(double minStep);
	void setMaximumStep(double maxStep);

	double getTolerance() const;
	double getMinimumStep() const;
	double getMaximumStep() const;
    double getdeterminant(double n, Vector3d dir, Vector3d n_vec);
    double bisection(double a, double b, Vector3d dir, Vector3d n_vec);
    double minFinder(double a, double b, Vector3d dir, Vector3d n_vec);
    double getTimeDelay(double n1, double n2, double l);
    Vector3d getPolarization(double n, Vector3d dir, Vector3d n_vec);
	Vector3d getEffectiveIndices(Vector3d dir, Vector3d n_vec);
    ElectricField apply_birefringence_1(Vector3d dir, Vector3d n_vec);
    ElectricField apply_birefringence(ElectricField Pulse, Vector3d dir, Vector3d n_vec);
	std::string getDescription() const;
};

} // namespace radiopropa

#endif // RADIOPROPA_PROPAGATIONCK_H
