#ifndef CRPROPA_OBSERVER_H
#define CRPROPA_OBSERVER_H

#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "../Candidate.h"
#include "../Module.h"
#include "../Referenced.h"
#include "../Vector3.h"
#include "../Geometry.h"

namespace radiopropa {

enum DetectionState {
	DETECTED, VETO, NOTHING
};

/**
 @class ObserverFeature
 @brief Abstract base class for features of cosmic ray observers
 */
class ObserverFeature: public Referenced {
protected:
	std::string description;
public:
	virtual DetectionState checkDetection(Candidate *candidate) const;
	virtual void onDetection(Candidate *candidate) const;
	virtual std::string getDescription() const;
};

/**
 @class Observer
 @brief General cosmic ray observer
 */
class Observer: public Module {
	std::string flagKey;
	std::string flagValue;
private:
	std::vector<ref_ptr<ObserverFeature> > features;
	ref_ptr<Module> detectionAction;
	bool clone;
	bool makeInactive;
public:
	Observer();
	void add(ObserverFeature *feature);
	void onDetection(Module *action, bool clone = false);
	void process(Candidate *candidate) const;
	std::string getDescription() const;
	void setFlag(std::string key, std::string value);
	void setDeactivateOnDetection(bool deactivate);
};

/**
 @class ObserverDetectAll
 @brief Detects all particles
 */
class ObserverDetectAll: public ObserverFeature {
public:
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};



/**
 @class ObserverSurface
 @brief Detects particles crossing the durface
 */
class ObserverSurface: public ObserverFeature {
	private:
		ref_ptr<Surface> surface;

	public:
		ObserverSurface(Surface* _surface);
		DetectionState checkDetection(Candidate *candidate) const;
		std::string getDescription() const;
};


/**
 @class ObserverTracking
 @brief Tracks particles inside a sphere
 */
class ObserverTracking: public ObserverFeature {
private:
	Vector3d center;
	double radius;
    double stepSize;
public:
	ObserverTracking(Vector3d center, double radius, double stepSize = 0);
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};


/**
 @class ObserverPoint
 @brief Detects particles when reaching x = 0

 This module limits the next step size to prevent candidates from overshooting.
 Should be renamed to Observer1D, once old observer-scheme is removed.
 */
class ObserverPoint: public ObserverFeature {
public:
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};


/**
 @class ObserverInactiveVeto
 @brief Veto for inactive candidates
 */
class ObserverInactiveVeto: public ObserverFeature {
public:
	DetectionState checkDetection(Candidate *candidate) const;
	std::string getDescription() const;
};


/**
 @class ObserverTimeEvolution
 @brief Observes the time evolution of the candidates (phase-space elements)
 This observer is very useful if the time evolution of the particle density is needed. It detects all candidates in regular timeintervals and limits the nextStep of candidates to prevent overshooting of detection intervals.
 */
class ObserverTimeEvolution: public ObserverFeature {
private:
  std::vector<double> detList;
public:
  ObserverTimeEvolution();
  ObserverTimeEvolution(double min, double dist, double numb);
  void addTime(const double &position);
  const std::vector<double>& getTimes() const;
  DetectionState checkDetection(Candidate *candidate) const;
  std::string getDescription() const;
};

}

#endif // CRPROPA_OBSERVER_H
