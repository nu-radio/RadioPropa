#include "radiopropa/module/Discontinuity.h"
#include <cmath>

namespace radiopropa{

Discontinuity::Discontinuity(Surface *_surface, double _n1, double _n2) : surface(_surface), n1(_n1), n2(_n2)
{

}

void Discontinuity::process(Candidate *candidate) const
{
		double cx = surface->distance(candidate->current.getPosition());
		double px = surface->distance(candidate->previous.getPosition());

		if (std::signbit(cx) == std::signbit(px))
		{
			candidate->limitNextStep(fabs(cx));
			return;
		}
		// Crossed the boundary, the secondary propagates further, while the
		// candidate is reflected.

		// calculate inersection point
		Vector3d dp = candidate->current.getPosition() - candidate->previous.getPosition();
		Vector3d intersectionPoint = candidate->previous.getPosition() + dp.getUnitVector() * px;

		// surface normal in intersection point 
		Vector3d localNormal= surface->normal(intersectionPoint);

		// correct n1, n2 ratio according to crossing direction
		double NR;
		if (px < 0)
				 NR = n1 / n2;
		else
				 NR = n2 / n1;

		// reflection according to Snell's law
		// angle to the surface normal alpha, beta and
		// sin/cos (alpha, beta) = salpha, sbeta, calpha, cbeta
		double calpha = fabs(candidate->current.getDirection().dot(localNormal));
		double salpha2 = 1 - calpha*calpha;
		double sbeta2 = NR * NR * salpha2;

		// Reflection coefficents
		double R_perp = 1.;
		double R_para = 1.;

		// vector in plane ip perpendicular to direction and normal
		Vector3d ip = localNormal.cross(candidate->current.getDirection());
		// Calculate amplitudes parallel and perpendicular for Fresnell
		const Vector3d &A = candidate->current.getAmplitude();
		//Vector3d Aperp = localNormal * A.dot(localNormal);
		const Vector3d Apara = ip * A.dot(ip);
		const Vector3d Aperp = A - Apara;

		if (fabs(sbeta2) < 1)
		{ // Partial reflection, calculate reflection and transmission
			// coefficient from Fresnell equations
			double cbeta = sqrt(1-sbeta2);

			// calculate reflection coefficients
			double T_perp = 2 * NR * calpha / (NR * calpha + cbeta);
			R_perp = (NR * calpha - cbeta) / (NR * calpha + cbeta);

			double T_para = 2 * NR * calpha / (calpha + NR * cbeta);
			R_para = (calpha - NR * cbeta) / (calpha + NR * cbeta);

			// add transmitted particle
			ref_ptr<Candidate> c2 = candidate->clone(false);

			Vector3d Vnew = -1. * NR * (localNormal.cross(ip)) + localNormal * (sqrt(1 - NR*NR * ip.dot(ip)));
			c2->current.setDirection(Vnew);

			// calculate + set amplitude
			double alpha = acos(calpha);
			double beta = acos(cbeta);
			const Vector3d Aperp_p = Aperp.getRotated(ip, beta - alpha);
			Vector3d TransmittedAmplitude = Apara * T_para + Aperp_p * T_perp;
			// correction factor to account for the increased beamwidth
			double c = sqrt(1. / NR * cbeta / calpha);
			TransmittedAmplitude *= c;

			c2->current.setAmplitude(TransmittedAmplitude);
			candidate->addSecondary(c2);
		}

		const Vector3d &V = candidate->current.getDirection();
		const Vector3d u = localNormal * (V.dot(localNormal));
		const Vector3d new_direction = V - u*2;
		candidate->current.setDirection(new_direction);

		// Reflected Amplitude
		const Vector3d Aperp_p = Aperp - localNormal * (Aperp.dot(localNormal)) * 2;
		const Vector3d ReflectedAmplitude = Apara * R_para + Aperp_p * R_perp;

		candidate->current.setAmplitude(ReflectedAmplitude);

		// update position slightly to move on correct side of plane
		Vector3d X = candidate->current.getPosition();
		candidate->current.setPosition(X + new_direction * candidate->getCurrentStep());

}

}
