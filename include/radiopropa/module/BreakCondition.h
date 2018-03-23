#ifndef CRPROPA_BREAKCONDITION_H
#define CRPROPA_BREAKCONDITION_H

#include "radiopropa/Module.h"

namespace radiopropa {

/**
 @class MaximumTrajectoryLength
 @brief Deactivates the candidate beyond a maximum trajectory length

 This modules deactivates the candidate at a given maximum trajectory length.
 In that case the property ("Deactivated", module::description) is set.
 It also limits the candidates next step size to ensure the maximum trajectory length is no exceeded.
 */
class MaximumTrajectoryLength: public AbstractCondition {
	double maxLength;
	std::vector<Vector3d> observerPositions;
public:
	MaximumTrajectoryLength(double length = 0);
	void setMaximumTrajectoryLength(double length);
	double getMaximumTrajectoryLength() const;
	void addObserverPosition(const Vector3d &position);
	const std::vector<Vector3d>& getObserverPositions() const;
	std::string getDescription() const;
	void process(Candidate *candidate) const;
};

/**
 @class MinimumFrequency
 @brief Deactivates the candidate below a minimum frequency

 This modules deactivates the candidate below a given minimum frequency.
 In that case the property ("Deactivated", module::description) is set.
 */
class MinimumFrequency: public AbstractCondition {
	double minFrequency;
public:
	MinimumFrequency(double minFrequency = 0);
	void setMinimumFrequency(double frequency);
	double getMinimumFrequency() const;
	std::string getDescription() const;
	void process(Candidate *candidate) const;
};


/**
 @class MinimumAmplitude
 @brief Deactivates the candidate below a minimum rigidity

 This modules deactivates the candidate below a given minimum rigidity (E/Z in EeV).
 In that case the property ("Deactivated", module::description) is set.
 */
class MinimumAmplitude: public AbstractCondition {
	double minAmplitude;
public:
	MinimumAmplitude(double minAmplitude = 0);
	void setMinimumAmplitude(double minAmplitude);
	double getMinimumAmplitude() const;
	std::string getDescription() const;
	void process(Candidate *candidate) const;
};

/**
 @class DetectionLength
 @brief Detects the candidate at a given trajectoryLength
 
 This break condition can be used for non-regular time observation of the particle density. See also TimeEvolutionObserver.
 */
class DetectionLength: public AbstractCondition {
	double detLength;
public:
	DetectionLength(double length = 0);
	void setDetectionLength(double length);
	double getDetectionLength() const;
	std::string getDescription() const;
	void process(Candidate *candidate) const;
};

} // namespace radiopropa

#endif // CRPROPA_BREAKCONDITION_H
