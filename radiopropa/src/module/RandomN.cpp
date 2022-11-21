#include "radiopropa/module/RandomN.h"
#include "radiopropa/Random.h"
#include <cmath>
#include <complex>

namespace radiopropa{

	RandomN::RandomN(Surface *_surface, double _dn, bool _surfacemode, double _fraction) : 
		surface(_surface), 
		dn(_dn), 
		surfacemode(_surfacemode),
		fraction(_fraction)

	{
	}

	void RandomN::process(Candidate *candidate) const
	{
		if (candidate -> current.getDirection().z < 0)	{
			return;
		}
		else {
			double cx = surface->distance(candidate->current.getPosition());
			double px = surface->distance(candidate->previous.getPosition());

			const Vector3d normal = surface->normal(candidate->current.getPosition());
			const Vector3d v = candidate->current.getDirection();
			const double cos_theta = v.dot(normal);
			const Vector3d u = normal * (cos_theta*v.getR());
			const Vector3d surface_direction = (v - u) / (v-u).getR();
			
			Random &random = Random::instance();
			double zCurr = candidate->current.getPosition().z;
			//double n2 = random.randNorm(n1, n1 / 2); //+-0.05
			double n1 = 1.78-0.481*exp(zCurr/37);		//calculate n for current z pos (ROSS ICE SHELF)
			/*if (zCurr < 14.9) {					//calculate n for current z pos (GREENLAND FIRN)
				double n1 = 1.78-0.31*exp(z/40.9)
			}
			else {
				double n1 = 1.78-0.502*exp(z/30.8)
			}
			*/

        	double n2 = random.randNorm(n1, dn); //+-0.05
			if (n2 < n1-2*dn or n2 > n1+2*dn){
				while (n2 < n1-2*dn or n2 > n1+2*dn){
					n2 = random.randNorm(n1, dn);
				}
			}
			
			//double n2 = random.randUniform(n1-0.05,n1+0.05);
			if (this->createdAtSurface(candidate)){
				//new direction parallel to layer
				candidate->current.setDirection(surface_direction);

				/*Propagation module bends the ray slightly downwards,
				resulting in a straigt line with a small negative slope
				with respect to the layer. Adjusting for the position 
				overcomes this*/
				this->positionCorrection(candidate, surface_direction);
				candidate->limitNextStep(tolerance);

			} 
			else {
				if (std::signbit(cx) == std::signbit(px))
				{
					candidate->limitNextStep(fabs(cx));
					return;
				} //else {
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

				//if (fabs(sbeta2) < 1)
				//if (1>0)
				{ // Partial reflection, calculate reflection and transmission
					// coefficient from Fresnell equations
				Random &random = Random::instance();
				double n2 = random.randNorm(n1, dn); //+-0.05
				if (n2 < n1-2*dn or n2 > n1+2*dn){
					while (n2 < n1-2*dn or n2 > n1+2*dn){
						n2 = random.randNorm(n1, dn);
					}
				}
				double NR;
				if (px < 0)
					NR = n1 / n2;
				else
					NR = n2 / n1;

				double cbeta = sqrt(1-sbeta2);

				// calculate reflection coefficients
				double T_perp = 2 * NR * calpha / (NR * calpha + cbeta);
				R_perp = (NR * calpha - cbeta) / (NR * calpha + cbeta);

				double T_para = 2 * NR * calpha / (calpha + NR * cbeta);
				R_para = (calpha - NR * cbeta) / (calpha + NR * cbeta);

				// add transmitted particle

				//ref_ptr<Candidate> c2 = candidate->clone(false);

				Vector3d transmitted_direction = -1. * NR * (localNormal.cross(ip)) + localNormal * (sqrt(1 - NR*NR * ip.dot(ip)));
				candidate->current.setDirection(transmitted_direction);

				// calculate + set amplitude
				double alpha = acos(calpha);
				double beta = acos(cbeta);
				const Vector3d Aperp_p = Aperp.getRotated(ip, beta - alpha);
				Vector3d TransmittedAmplitude = Apara * T_para + Aperp_p * T_perp;
				// correction factor to account for the increased beamwidth
				double c = sqrt(1. / NR * cbeta / calpha);
				TransmittedAmplitude *= c;

				candidate->current.setAmplitude(TransmittedAmplitude);
				} 
			} 
		}
	}
	std::string RandomN::getDescription() const {
		std::stringstream ss;
		ss << "RandomN";
		ss << "\n    " << surface->getDescription() << "\n";
		ss << "    dn: " << dn << "\n";

		return ss.str();
	}
	double RandomN::getFraction() const{
		return fraction;
	}
	void RandomN::setFraction(double new_fraction){
		fraction = new_fraction;
	}

	bool RandomN::parallelToSurface(Vector3d position, Vector3d direction) const{
			Vector3d normal = surface->normal(position);
	        double cos_theta = direction.dot(normal);
	        return (abs(cos_theta) < 0.001);
	}
	bool RandomN::atSurface(Vector3d position) const{
			double distance = surface->distance(position);
	        return (abs(distance) <= tolerance);
	}
	bool RandomN::createdAtSurface(Candidate *candidate) const{
		Vector3d position = candidate->created.getPosition();
	    return this->atSurface(position);
	}
	void RandomN::positionCorrection(Candidate* candidate, Vector3d new_direction) const{
		Vector3d c = candidate->current.getPosition();
		Vector3d p = candidate->previous.getPosition();
		double step_size = (c-p).getR();
		Vector3d new_position = p + new_direction*step_size;
		candidate->current.setPosition(new_position);
	}
	void RandomN::setSurfacemode(bool mode){
		surfacemode = mode;
	}


	TransmissiveRandomNLayer::TransmissiveRandomNLayer(Surface *surface, double transmission) : 
		surface(surface), transmission(transmission)
	{
	}
	void TransmissiveRandomNLayer::process(Candidate *candidate) const
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
	std::string TransmissiveRandomNLayer::getDescription() const {
		std::stringstream ss;
		ss << "TransmissiveLayer";
		ss << "\n    " << surface->getDescription() << "\n";
		ss << "    transmission coefficent: " << transmission << "\n";

		return ss.str();
	}


	ReflectiveRandomNLayer::ReflectiveRandomNLayer(Surface *surface, double reflection) : 
		surface(surface), reflection(reflection)
	{
	}
	void ReflectiveRandomNLayer::process(Candidate *candidate)
	{
		double cx = surface->distance(candidate->current.getPosition());
		double px = surface->distance(candidate->previous.getPosition());

		if (std::signbit(cx) == std::signbit(px)){
			candidate->limitNextStep(fabs(cx));
			return;
		} 
	}
	std::string ReflectiveRandomNLayer::getDescription() const {
		std::stringstream ss;
		ss << "ReflectiveLayer";
		ss << "\n    " << surface->getDescription() << "\n";
		ss << "    reflection coefficent: " << reflection << "\n";

		return ss.str();
	}
	int ReflectiveRandomNLayer::getTimesReflectedoff(Candidate *candidate){
		return times_reflectedoff[candidate];
	}
}
