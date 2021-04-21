#include "radiopropa/module/DensityPerturbation.h"
#include <cmath>
#include <stdlib.h>

namespace radiopropa{

	PerturbationLayer::PerturbationLayer(Surface *_surface, double _thickness, double _threshold, double _fraction): 
		surface(_surface),
		thickness(_thickness),
		threshold(_threshold),
		fraction(_fraction)
	{}
	PerturbationLayer::~PerturbationLayer(){}

	void PerturbationLayer::process(Candidate* candidate) const{
		candidate->limitNextStep(abs(abs(surface->distance(candidate->current.getPosition()))-thickness/2));

        if (this->inLayer(candidate->current.getPosition())){
            const Vector3d normal = surface->normal(candidate->current.getPosition());
            const Vector3d v = candidate->current.getDirection();
            const double cos_theta = v.dot(normal);

            const Vector3d normal_previous = surface->normal(candidate->previous.getPosition());
            const double cos_theta_previous = candidate->previous.getDirection().dot(normal_previous);
            
            if (this->createdInLayer(candidate)){
                const Vector3d u = normal * (cos_theta*v.getR());
                const Vector3d new_direction = (v - u) / (v-u).getR(); //new direction parallel to layer
                candidate->current.setDirection(new_direction);

                /*Propagation module bends the ray slightly downwards,
                resulting in a straigt line with a small negative slope
                with respect to the layer. Adjusting for the position 
                overcomes this*/
                this->positionCorrection(candidate, new_direction);
                candidate->limitNextStep(thickness/4);

            } else if (this->inLayer(candidate->previous.getPosition()) and (abs(cos_theta) < threshold) and (std::signbit(cos_theta) != std::signbit(cos_theta_previous))) {
                bool has_daugther_in_layer = false;
                for (auto& secondary: candidate->secondaries){
                    has_daugther_in_layer = (has_daugther_in_layer or this->createdInLayer(secondary));
                }

                if (not has_daugther_in_layer) {
                    //The secondary propagates further in layer because of very small fraction
                    ref_ptr<Candidate> secondary = candidate->clone(false);
                    secondary->created = candidate->previous;
                    Vector3d E = candidate->current.getAmplitude();
                    secondary->current.setAmplitude(E*fraction);
                    
                    const Vector3d u = normal * (cos_theta*v.getR());
                    const Vector3d new_direction = (v - u) / (v-u).getR(); //new direction parallel to layer
                    secondary->current.setDirection(new_direction);

                    /*Propagation module bends the ray slightly downwards,
                    resulting in a straigt line with a small negative slope
                    with respect to the layer. Adjusting for the position 
                    overcomes this.*/
                    positionCorrection(secondary, new_direction);

                    secondary->limitNextStep(thickness/4);
                    candidate->addSecondary(secondary);
                }
            }
        }
	}

	std::string PerturbationLayer::getDescription() const {
		std::stringstream ss;
		ss << "Density Perturbation";
		ss << "\n    " << surface->getDescription() << "\n";
		ss << "    thickness: " << thickness << "\n";
		ss << "    threshold: " << threshold << "\n";
		ss << "    fraction : " << fraction;

		return ss.str();
	}


	ref_ptr<Surface> PerturbationLayer::getSurface() const{
		return surface;
	}
	double PerturbationLayer::getThickness() const{
		return thickness;
	}
	void PerturbationLayer::setThickness(double new_thickness){
		thickness = new_thickness;
	}
	double PerturbationLayer::getThreshold() const{
		return threshold;
	}
	void PerturbationLayer::setThreshold(double new_threshold){
		threshold = new_threshold;
	}
	double PerturbationLayer::getFraction() const{
		return fraction;
	}
	void PerturbationLayer::setFraction(double new_fraction){
		fraction = new_fraction;
	}


	void PerturbationLayer::positionCorrection(Candidate* candidate, Vector3d new_direction) const{
		Vector3d c = candidate->current.getPosition();
		Vector3d p = candidate->previous.getPosition();
		double step_size = (c-p).getR();
		Vector3d new_position = p + new_direction*step_size;
		candidate->current.setPosition(new_position);
	}
	bool PerturbationLayer::parallelToLayer(Vector3d position, Vector3d direction) const{
		Vector3d normal = surface->normal(position);
        double cos_theta = direction.dot(normal);
        return (abs(cos_theta) < threshold);
	}
	bool PerturbationLayer::inLayer(Vector3d position) const{
		double distance = surface->distance(position);
        return (abs(distance) <= thickness/2);
	}
    bool PerturbationLayer::createdInLayer(Candidate* candidate) const{
        return inLayer(candidate->created.getPosition());
    }
    PerturbationLayer* PerturbationLayer::clone() const{
    	PerturbationLayer* cloned = new PerturbationLayer;
    	cloned->surface = surface;
    	cloned->thickness = thickness;
    	cloned->threshold = threshold;
    	cloned->fraction = fraction;
        return cloned;
    }


    
    PerturbationHorizontal::PerturbationHorizontal(double _z, double _thickness, double _threshold, double _fraction):
    	PerturbationLayer(new Plane(Vector3d(0,0,_z), Vector3d(0,0,1)), _thickness, _threshold, _fraction), z(_z) {
	}
    PerturbationHorizontal::~PerturbationHorizontal(){}

    void PerturbationHorizontal::positionCorrection(Candidate *candidate, Vector3d new_direction) const{
    	Vector3d p = candidate->previous.getPosition();
        Vector3d c = candidate->current.getPosition();
        Vector3d new_position = Vector3d(c.x,c.y,p.z);
        candidate->current.setPosition(new_position);
    }

    PerturbationLayer* PerturbationHorizontal::clone() const{
    	PerturbationHorizontal* cloned = new PerturbationHorizontal;
    	cloned->surface = surface;
    	cloned->thickness = thickness;
    	cloned->threshold = threshold;
    	cloned->fraction = fraction;
        return cloned;
    }
	
}