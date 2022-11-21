#ifndef RandomN_H
#define RandomN_H

#include "radiopropa/Candidate.h"
#include "radiopropa/Module.h"
#include "radiopropa/Geometry.h"
#include "radiopropa/Referenced.h"
#include "radiopropa/Units.h"
#include "radiopropa/Random.h"

namespace radiopropa
{
/**
 @class RandomN
 @brief Dcsiontinuity in the refractive index at which secondary rays are created according to   Fresnell equations.
 */
class RandomN: public Module {
private:
	ref_ptr<Surface> surface;
	double dn, fraction;
    bool surfacemode;
    double tolerance = 0.01*meter;


public:
    RandomN(Surface *_surface, double _dn, bool surfacemode=false, double fraction=0.023); 
    void process(Candidate *candidate) const;
    std::string getDescription() const;
    
    double getFraction() const;
    void setFraction(double new_fraction);
    
    bool parallelToSurface(Vector3d position, Vector3d direction) const;
    bool atSurface(Vector3d position) const;
    bool createdAtSurface(Candidate *candidate) const;
    void positionCorrection(Candidate* candidate, Vector3d new_direction) const;
    void setSurfacemode(bool mode);
};

class TransmissiveRandomNLayer: public Module {
private:
    ref_ptr<Surface> surface;
    double transmission;
public:
    TransmissiveRandomNLayer(Surface *surface, double transmission=1); 
    void process(Candidate *candidate) const;
    std::string getDescription() const;
};

class ReflectiveRandomNLayer: public Module {
private:
    ref_ptr<Surface> surface;
    double reflection;
    std::map<Candidate*, int> times_reflectedoff;
public:
    ReflectiveRandomNLayer(Surface *surface, double reflection=1); 
    void process(Candidate *candidate);
    std::string getDescription() const;

    int getTimesReflectedoff(Candidate *candidate);
};


}

#endif // RandomN_H
