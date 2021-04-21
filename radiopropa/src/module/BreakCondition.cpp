#include "radiopropa/module/BreakCondition.h"
#include "radiopropa/Units.h"

#include <sstream>

namespace radiopropa {

MaximumTrajectoryLength::MaximumTrajectoryLength(double maxLength) :
		maxLength(maxLength) {
}

void MaximumTrajectoryLength::setMaximumTrajectoryLength(double length) {
	maxLength = length;
}

double MaximumTrajectoryLength::getMaximumTrajectoryLength() const {
	return maxLength;
}

void MaximumTrajectoryLength::addObserverPosition(const Vector3d& position) {
	observerPositions.push_back(position);
}

const std::vector<Vector3d>& MaximumTrajectoryLength::getObserverPositions() const {
	return observerPositions;
}

std::string MaximumTrajectoryLength::getDescription() const {
	std::stringstream s;
	s << "Maximum trajectory length: " << maxLength / meter  << " m, ";
	s << "Flag: '" << rejectFlagKey << "' -> '" << rejectFlagValue << "', ";
	s << "MakeInactive: " << (makeRejectedInactive ? "yes" : "no");
	if (rejectAction.valid())
		s << ", Action: " << rejectAction->getDescription();
	s << "\n  Observer positions: \n";
	for (size_t i = 0; i < observerPositions.size(); i++)
		s << "    - " << observerPositions[i] / Mpc << " Mpc\n";
	return s.str();
}

void MaximumTrajectoryLength::process(Candidate *c) const {
	double length = c->getTrajectoryLength();
	Vector3d position = c->current.getPosition();

	if(observerPositions.size()) {
		bool inRange = false;
		for (size_t i = 0; i < observerPositions.size(); i++) {
			double distance = position.getDistanceTo(observerPositions[i]);
			if (distance + length < maxLength)
				inRange = true;
		}
		if (!inRange) {
			reject(c);
			return;
		}
	}

	if (length >= maxLength) {
		reject(c);
	} else {
		c->limitNextStep(maxLength - length);
	}
}

//*****************************************************************************
MinimumFrequency::MinimumFrequency(double minFrequency) :
		minFrequency(minFrequency) {
}

void MinimumFrequency::setMinimumFrequency(double frequency) {
	minFrequency = frequency;
}

double MinimumFrequency::getMinimumFrequency() const {
	return minFrequency;
}

void MinimumFrequency::process(Candidate *c) const {
	if (c->current.getFrequency() > minFrequency)
		return;
	else
		reject(c);
}

std::string MinimumFrequency::getDescription() const {
	std::stringstream s;
	s << "Minimum frequency: " << minFrequency << ", ";
	s << "Flag: '" << rejectFlagKey << "' -> '" << rejectFlagValue << "', ";
	s << "MakeInactive: " << (makeRejectedInactive ? "yes" : "no");
	if (rejectAction.valid())
		s << ", Action: " << rejectAction->getDescription();
	return s.str();
}

//*****************************************************************************
MinimumAmplitude::MinimumAmplitude(double minAmplitude) :
		minAmplitude(minAmplitude) {
}

void MinimumAmplitude::setMinimumAmplitude(double minAmplitude) {
	this->minAmplitude = minAmplitude;
}

double MinimumAmplitude::getMinimumAmplitude() const {
	return minAmplitude;
}

void MinimumAmplitude::process(Candidate *c) const {
	if (c->current.getAmplitude().getR() < minAmplitude)
		reject(c);
}

std::string MinimumAmplitude::getDescription() const {
	std::stringstream s;
	s << "Minimum amplitude: " << minAmplitude << " [a.u.], ";
	s << "Flag: '" << rejectFlagKey << "' -> '" << rejectFlagValue << "', ";
	s << "MakeInactive: " << (makeRejectedInactive ? "yes" : "no");
	if (rejectAction.valid())
		s << ", Action: " << rejectAction->getDescription();
	return s.str();
}

//*****************************************************************************


//*****************************************************************************
DetectionLength::DetectionLength(double detLength) :
		detLength(detLength) {
}

void DetectionLength::setDetectionLength(double length) {
	detLength = length;
}

double DetectionLength::getDetectionLength() const {
	return detLength;
}


std::string DetectionLength::getDescription() const {
	std::stringstream s;
	s << "Detection length: " << detLength / meter<< " m, ";
	s << "Flag: '" << rejectFlagKey << "' -> '" << rejectFlagValue << "', ";
	s << "MakeInactive: " << (makeRejectedInactive ? "yes" : "no");
	if (rejectAction.valid())
		s << ", Action: " << rejectAction->getDescription();
	return s.str();
}

void DetectionLength::process(Candidate *c) const {
	double length = c->getTrajectoryLength();
	double step = c->getCurrentStep();

	if (length >= detLength && length - step < detLength) {
		reject(c);
	} else {
		c->limitNextStep(detLength - length);
	}
}


} // namespace radiopropa
