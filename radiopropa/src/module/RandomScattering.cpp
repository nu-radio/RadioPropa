#include "radiopropa/module/RandomScattering.h"
#include "radiopropa/Random.h"
#include <cmath>
#include <complex>

namespace radiopropa{

	RandomScattering::RandomScattering(Surface *_surface, double _delta) : 
		surface(_surface),
		delta(_delta)
	{
	}

	void RandomScattering::process(Candidate *candidate) const
	{
		if (candidate -> current.getDirection().z <= 0)	{
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

			
			if (this->createdAtSurface(candidate)){
				//new direction parallel to layer
				//candidate->current.setDirection(surface_direction);

				/*Propagation module bends the ray slightly downwards,
				resulting in a straigt line with a small negative slope
				with respect to the layer. Adjusting for the position 
				overcomes this*/
				//this->positionCorrection(candidate, surface_direction);
				candidate->limitNextStep(tolerance);

			} else {
				if (std::signbit(cx) == std::signbit(px))
				{
					candidate->limitNextStep(fabs(cx));
					return;
				} else {

					//create random Vector and set direction of candidate to direction of random Vector
					Random &random = Random::instance();
        			double rand = random.randNorm(0, delta); //+-0.05
					if (rand < -2*delta or rand > 2*delta){
						while (rand < -2*delta or rand > 2*delta){
							rand = random.randNorm(0, delta);
						}
					}
					
					Vector3d new_direction = candidate->current.getDirection()+Vector3d(0,0,rand);
					candidate->current.setDirection(new_direction);

					// update position slightly to move on correct side of plane
					Vector3d X = candidate->current.getPosition();
					candidate->current.setPosition(X + new_direction * candidate->getCurrentStep());
				}
			}
		}
	}
    
	std::string RandomScattering::getDescription() const {
		std::stringstream ss;
		ss << "RandomScattering";
		ss << "\n    " << surface->getDescription() << "\n";
		return ss.str();
	}
    
	bool RandomScattering::parallelToSurface(Vector3d position, Vector3d direction) const{
			Vector3d normal = surface->normal(position);
	        double cos_theta = direction.dot(normal);
	        return (abs(cos_theta) < 0.001);
	}
	bool RandomScattering::atSurface(Vector3d position) const{
			double distance = surface->distance(position);
	        return (abs(distance) <= tolerance);
	}
	bool RandomScattering::createdAtSurface(Candidate *candidate) const{
		Vector3d position = candidate->created.getPosition();
	    return this->atSurface(position);
	}
	void RandomScattering::positionCorrection(Candidate* candidate, Vector3d new_direction) const{
		Vector3d c = candidate->current.getPosition();
		Vector3d p = candidate->previous.getPosition();
		double step_size = (c-p).getR();
		Vector3d new_position = p + new_direction*step_size;
		candidate->current.setPosition(new_position);
	}

	RandomScatteringLayer::RandomScatteringLayer(Surface *surface, double transmission) : 
		surface(surface), transmission(transmission)
	{
	}
	void RandomScatteringLayer::process(Candidate *candidate) const
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
    	std::string RandomScatteringLayer::getDescription() const {
		std::stringstream ss;
		ss << "RandomScatteringLayer";
		return ss.str();
	    }
}
