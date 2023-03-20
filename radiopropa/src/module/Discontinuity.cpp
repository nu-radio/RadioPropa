#include "radiopropa/module/Discontinuity.h"
#include <cmath>
#include <complex>

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
			// check direction of ray
			if (candidate -> current.getDirection().z < 0)
			{
				localNormal.z = -1.*localNormal.z;
			}

			// reflection according to Snell's law
			// angle to the surface normal alpha, beta and
			// sin/cos (alpha, beta) = salpha, sbeta, calpha, cbeta
			double calpha = fabs(candidate->current.getDirection().dot(localNormal));
			double salpha2 = 1 - calpha*calpha;
			double sbeta2 = NR * NR * salpha2;

			candidate->appendReflectionAngle(acos(calpha));

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
	std::string Discontinuity::getDescription() const {
		std::stringstream ss;
		ss << "Discontinuity";
		ss << "\n    " << surface->getDescription() << "\n";
		ss << "    n1: " << n1 << "\n";
		ss << "    n2: " << n2;

		return ss.str();
	}


	TransmissiveLayer::TransmissiveLayer(Surface *surface, double transmission) : 
		surface(surface), transmission(transmission)
	{
	}
	void TransmissiveLayer::process(Candidate *candidate) const
	{
		double cx = surface->distance(candidate->current.getPosition());
		double px = surface->distance(candidate->previous.getPosition());

		if (std::signbit(cx) == std::signbit(px)){
			candidate->limitNextStep(fabs(cx));
			return;
		} else {
			candidate->current.setAmplitude(candidate->current.getAmplitude()*transmission);
		}
	}
	std::string TransmissiveLayer::getDescription() const {
		std::stringstream ss;
		ss << "TransmissiveLayer";
		ss << "\n    " << surface->getDescription() << "\n";
		ss << "    transmission coefficent: " << transmission << "\n";

		return ss.str();
	}


	ReflectiveLayer::ReflectiveLayer(Surface *surface, double reflection) : 
		surface(surface), reflection(reflection), times_reflectedoff()
	{
	}
	void ReflectiveLayer::process(Candidate *candidate) const
	{
		double cx = surface->distance(candidate->current.getPosition());
		double px = surface->distance(candidate->previous.getPosition());

		if (std::signbit(cx) == std::signbit(px)){
			candidate->limitNextStep(fabs(cx));
			return;
		} else {
			candidate->current.setAmplitude(candidate->current.getAmplitude()*reflection);

			Vector3d normal = surface->normal(candidate->current.getPosition());
            Vector3d v = candidate->current.getDirection();
            double cos_theta = v.dot(normal);
            Vector3d u = normal * (cos_theta);
            Vector3d new_direction = v - u*2; //new direction due to reflection of surface
            candidate->current.setDirection(new_direction);

            // update position slightly to move on correct side of plane
            Vector3d x = candidate->current.getPosition();
            candidate->current.setPosition(x + new_direction * candidate->getCurrentStep());

            //keeping track of the amount of reflections is currently not working
            /*int amount = 0;
            if (times_reflectedoff.find(candidate) != times_reflectedoff.end()) {
               amount = 1;
            } else {
                amount = this->getTimesReflectedOff(candidate) + 1;
            }
            this->setTimesReflectedOff(candidate, amount);*/
		}
	}
	std::string ReflectiveLayer::getDescription() const {
		std::stringstream ss;
		ss << "ReflectiveLayer";
		ss << "\n    " << surface->getDescription() << "\n";
		ss << "    reflection coefficent: " << reflection << "\n";

		return ss.str();
	}
	
	//keeping track of the amount of reflections is currently not working
	/*int ReflectiveLayer::getTimesReflectedOff(Candidate *candidate) const{
		return times_reflectedoff.at(candidate);
	}
	void ReflectiveLayer::setTimesReflectedOff(Candidate *candidate, int amount){
		times_reflectedoff[candidate] = amount;	
	}*/
}
