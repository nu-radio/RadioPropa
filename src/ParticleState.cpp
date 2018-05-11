#include "radiopropa/ParticleState.h"
#include "radiopropa/Units.h"
#include "radiopropa/Common.h"

#include <stdlib.h>
#include <sstream>

namespace radiopropa {

ParticleState::ParticleState(int id, double E, Vector3d pos, Vector3d dir) {
	setId(id);
	setFrequency(E);
	setPosition(pos);
	setDirection(dir);
}

void ParticleState::setPosition(const Vector3d &pos) {
	position = pos;
}

const Vector3d &ParticleState::getPosition() const {
	return position;
}

void ParticleState::setDirection(const Vector3d &dir) {
	direction = dir / dir.getR();
}

const Vector3d &ParticleState::getDirection() const {
	return direction;
}

void ParticleState::setAmplitude(double newAmplitude) {
	amplitude = newAmplitude; 
}

double ParticleState::getAmplitude() const {
	return amplitude;
}


void ParticleState::setFrequency(double newFrequency) {
	frequency = std::max(0., newFrequency); // prevent negative energies
}

double ParticleState::getFrequency() const {
	return frequency;
}

double ParticleState::getRigidity() const {
	return fabs(frequency / charge);
}

void ParticleState::setId(int newId) {

}

int ParticleState::getId() const {
	return id;
}

double ParticleState::getMass() const {
	return pmass;
}

double ParticleState::getCharge() const {
	return charge;
}

double ParticleState::getLorentzFactor() const {
	return frequency / (pmass * c_squared);
}

void ParticleState::setLorentzFactor(double lf) {
	lf = std::max(0., lf); // prevent negative Lorentz factors
	frequency = lf * pmass * c_squared;
}

Vector3d ParticleState::getVelocity() const {
	return direction * c_light;
}

Vector3d ParticleState::getMomentum() const {
	return direction * (frequency / c_light);
}

std::string ParticleState::getDescription() const {
	std::stringstream ss;
	ss << "Particle " << id << ", ";
	ss << "E = " << frequency / EeV << " EeV, ";
	ss << "x = " << position / Mpc << " Mpc, ";
	ss << "p = " << direction;
	return ss.str();
}

} // namespace radiopropa
