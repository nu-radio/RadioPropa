#include "radiopropa/Source.h"
#include "radiopropa/Random.h"
#include "radiopropa/Cosmology.h"
#include "radiopropa/Common.h"
#include "radiopropa/Units.h"

#ifdef CRPROPA_HAVE_MUPARSER
#include "muParser.h"
#endif

#include <sstream>
#include <stdexcept>

namespace radiopropa {

// Source ---------------------------------------------------------------------
void Source::add(SourceFeature* property) {
	features.push_back(property);
}

ref_ptr<Candidate> Source::getCandidate() const {
	ref_ptr<Candidate> candidate = new Candidate();
	for (int i = 0; i < features.size(); i++)
		(*features[i]).prepareCandidate(*candidate);
	return candidate;
}

std::string Source::getDescription() const {
	std::stringstream ss;
	ss << "Cosmic ray source\n";
	for (int i = 0; i < features.size(); i++)
		ss << "    " << features[i]->getDescription();
	return ss.str();
}

// SourceList------------------------------------------------------------------
void SourceList::add(Source* source, double weight) {
	sources.push_back(source);
	if (cdf.size() > 0)
		weight += cdf.back();
	cdf.push_back(weight);
}

ref_ptr<Candidate> SourceList::getCandidate() const {
	if (sources.size() == 0)
		throw std::runtime_error("SourceList: no sources set");
	size_t i = Random::instance().randBin(cdf);
	return (sources[i])->getCandidate();
}

std::string SourceList::getDescription() const {
	std::stringstream ss;
	ss << "List of cosmic ray sources\n";
	for (int i = 0; i < sources.size(); i++)
		ss << "  " << sources[i]->getDescription();
	return ss.str();
}

// SourceFeature---------------------------------------------------------------
void SourceFeature::prepareCandidate(Candidate& candidate) const {
	ParticleState &source = candidate.source;
	prepareParticle(source);
	candidate.created = source;
	candidate.current = source;
	candidate.previous = source;
}

std::string SourceFeature::getDescription() const {
	return description;
}


// ----------------------------------------------------------------------------
SourceFrequency::SourceFrequency(double frequency) :
		E(frequency) {
	setDescription();
}

void SourceFrequency::prepareParticle(ParticleState& p) const {
	p.setFrequency(E);
}

void SourceFrequency::setDescription() {
	std::stringstream ss;
	ss << "SourceFrequency: " << E / EeV << " EeV\n";
	description = ss.str();
}



// ----------------------------------------------------------------------------
SourcePosition::SourcePosition(Vector3d position) :
		position(position) {
	setDescription();
}

SourcePosition::SourcePosition(double d) :
		position(Vector3d(d, 0, 0)) {
	setDescription();
}

void SourcePosition::prepareParticle(ParticleState& particle) const {
	particle.setPosition(position);
}

void SourcePosition::setDescription() {
	std::stringstream ss;
	ss << "SourcePosition: " << position / Mpc << " Mpc\n";
	description = ss.str();
}

// ----------------------------------------------------------------------------
SourceMultiplePositions::SourceMultiplePositions() {
	setDescription();
}

void SourceMultiplePositions::add(Vector3d pos, double weight) {
	positions.push_back(pos);
	if (cdf.size() > 0)
		weight += cdf.back();
	cdf.push_back(weight);
}

void SourceMultiplePositions::prepareParticle(ParticleState& particle) const {
	if (positions.size() == 0)
		throw std::runtime_error("SourceMultiplePositions: no position set");
	size_t i = Random::instance().randBin(cdf);
	particle.setPosition(positions[i]);
}

void SourceMultiplePositions::setDescription() {
	std::stringstream ss;
	ss << "SourceMultiplePositions: Random position from list\n";
	for (int i = 0; i < positions.size(); i++)
		ss << "  " << positions[i] / Mpc << " Mpc\n";
	description = ss.str();
}

// ----------------------------------------------------------------------------
SourceUniformSphere::SourceUniformSphere(Vector3d center, double radius) :
		center(center), radius(radius) {
	setDescription();
}

void SourceUniformSphere::prepareParticle(ParticleState& particle) const {
	Random &random = Random::instance();
	double r = pow(random.rand(), 1. / 3.) * radius;
	particle.setPosition(center + random.randVector() * r);
}

void SourceUniformSphere::setDescription() {
	std::stringstream ss;
	ss << "SourceUniformSphere: Random position within a sphere at ";
	ss << center / Mpc << " Mpc with";
	ss  << radius / Mpc << " Mpc radius\n";
	description = ss.str();
}

// ----------------------------------------------------------------------------
SourceUniformShell::SourceUniformShell(Vector3d center, double radius) :
		center(center), radius(radius) {
	setDescription();
}

void SourceUniformShell::prepareParticle(ParticleState& particle) const {
	Random &random = Random::instance();
	particle.setPosition(center + random.randVector() * radius);
}

void SourceUniformShell::setDescription() {
	std::stringstream ss;
	ss << "SourceUniformShell: Random position on a spherical shell at ";
	ss << center / Mpc << " Mpc with ";
	ss << radius / Mpc << " Mpc radius\n";
	description = ss.str();
}

// ----------------------------------------------------------------------------
SourceUniformBox::SourceUniformBox(Vector3d origin, Vector3d size) :
		origin(origin), size(size) {
	setDescription();
}

void SourceUniformBox::prepareParticle(ParticleState& particle) const {
	Random &random = Random::instance();
	Vector3d pos(random.rand(), random.rand(), random.rand());
	particle.setPosition(pos * size + origin);
}

void SourceUniformBox::setDescription() {
	std::stringstream ss;
	ss << "SourceUniformBox: Random uniform position in box with ";
	ss << "origin = " << origin / Mpc << " Mpc and ";
	ss << "size = " << size / Mpc << " Mpc\n";
	description = ss.str();
}

// ---------------------------------------------------------------------------
SourceUniformCylinder::SourceUniformCylinder(Vector3d origin, double height, double radius) :
    origin(origin), height(height), radius(radius) {
}

void SourceUniformCylinder::prepareParticle(ParticleState& particle) const {
  Random &random = Random::instance();
  double phi = 2*M_PI*random.rand();
  double RandRadius = radius*pow(random.rand(), 1. / 2.);
  Vector3d pos(cos(phi)*RandRadius, sin(phi)*RandRadius, (-0.5+random.rand())*height);
  particle.setPosition(pos + origin);
  }

void SourceUniformCylinder::setDescription() {
	std::stringstream ss;
	ss << "SourceUniformCylinder: Random uniform position in cylinder with ";
	ss << "origin = " << origin / Mpc << " Mpc and ";
	ss << "radius = " << radius / Mpc << " Mpc and";
	ss << "height = " << height / Mpc << " Mpc\n";
	description = ss.str();
}


// ----------------------------------------------------------------------------
SourceUniform1D::SourceUniform1D(double minD, double maxD, bool withCosmology) {
	this->withCosmology = withCosmology;
	if (withCosmology) {
		this->minD = comoving2LightTravelDistance(minD);
		this->maxD = comoving2LightTravelDistance(maxD);
	} else {
		this->minD = minD;
		this->maxD = maxD;
	}
	setDescription();
}

void SourceUniform1D::prepareParticle(ParticleState& particle) const {
	Random& random = Random::instance();
	double d = random.rand() * (maxD - minD) + minD;
	if (withCosmology)
		d = lightTravel2ComovingDistance(d);
	particle.setPosition(Vector3d(d, 0, 0));
}

void SourceUniform1D::setDescription() {
	std::stringstream ss;
	ss << "SourceUniform1D: Random uniform position in D = ";
	ss << minD / Mpc << " - " << maxD / Mpc << " Mpc";
	if (withCosmology)
		ss << " (including cosmology)";
	ss << "\n";
	description = ss.str();
}

// ----------------------------------------------------------------------------
SourceDensityGrid::SourceDensityGrid(ref_ptr<ScalarGrid> grid) :
		grid(grid) {
	float sum = 0;
	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				sum += grid->get(ix, iy, iz);
				grid->get(ix, iy, iz) = sum;
			}
		}
	}
	setDescription();
}

void SourceDensityGrid::prepareParticle(ParticleState& particle) const {
	Random &random = Random::instance();

	// draw random bin
	size_t i = random.randBin(grid->getGrid());
	Vector3d pos = grid->positionFromIndex(i);

	// draw uniform position within bin
	double dx = random.rand() - 0.5;
	double dy = random.rand() - 0.5;
	double dz = random.rand() - 0.5;
	pos += Vector3d(dx, dy, dz) * grid->getSpacing();

	particle.setPosition(pos);
}

void SourceDensityGrid::setDescription() {
	description = "SourceDensityGrid: 3D source distribution according to density grid\n";
}

// ----------------------------------------------------------------------------
SourceDensityGrid1D::SourceDensityGrid1D(ref_ptr<ScalarGrid> grid) :
		grid(grid) {
	if (grid->getNy() != 1)
		throw std::runtime_error("SourceDensityGrid1D: Ny != 1");
	if (grid->getNz() != 1)
		throw std::runtime_error("SourceDensityGrid1D: Nz != 1");

	float sum = 0;
	for (int ix = 0; ix < grid->getNx(); ix++) {
		sum += grid->get(ix, 0, 0);
		grid->get(ix, 0, 0) = sum;
	}
	setDescription();
}

void SourceDensityGrid1D::prepareParticle(ParticleState& particle) const {
	Random &random = Random::instance();

	// draw random bin
	size_t i = random.randBin(grid->getGrid());
	Vector3d pos = grid->positionFromIndex(i);

	// draw uniform position within bin
	double dx = random.rand() - 0.5;
	pos.x += dx * grid->getSpacing();

	particle.setPosition(pos);
}

void SourceDensityGrid1D::setDescription() {
	description = "SourceDensityGrid1D: 1D source distribution according to density grid\n";
}

// ----------------------------------------------------------------------------
SourceIsotropicEmission::SourceIsotropicEmission() {
	setDescription();
}

void SourceIsotropicEmission::prepareParticle(ParticleState& particle) const {
	Random &random = Random::instance();
	particle.setDirection(random.randVector());
}

void SourceIsotropicEmission::setDescription() {
	description = "SourceIsotropicEmission: Random isotropic direction\n";
}

// ----------------------------------------------------------------------------
SourceDirection::SourceDirection(Vector3d direction) :
		direction(direction) {
	setDescription();
}

void SourceDirection::prepareParticle(ParticleState& particle) const {
	particle.setDirection(direction);
}

void SourceDirection::setDescription() {
	std::stringstream ss;
	ss <<  "SourceDirection: Emission direction = " << direction << "\n";
	description = ss.str();
}

// ----------------------------------------------------------------------------
SourceEmissionMap::SourceEmissionMap(EmissionMap *emissionMap) : emissionMap(emissionMap) {
	setDescription();
}

void SourceEmissionMap::prepareCandidate(Candidate &candidate) const {
	if (emissionMap) {
		bool accept = emissionMap->checkDirection(candidate.source);
		candidate.setActive(accept);
	}
}

void SourceEmissionMap::setDescription() {
	description = "SourceEmissionMap: accept only directions from emission map\n";
}

void SourceEmissionMap::setEmissionMap(EmissionMap *emissionMap) {
	this->emissionMap = emissionMap;
}

// ----------------------------------------------------------------------------
SourceEmissionCone::SourceEmissionCone(Vector3d direction, double aperture) :
		direction(direction), aperture(aperture) {
	setDescription();
}

void SourceEmissionCone::prepareParticle(ParticleState& particle) const {
	Random &random = Random::instance();
	particle.setDirection(random.randConeVector(direction, aperture));
}

void SourceEmissionCone::setDescription() {
	std::stringstream ss;
	ss << "SourceEmissionCone: Jetted emission in ";
	ss << "direction = " << direction << " with ";
	ss << "half-opening angle = " << aperture << " rad\n";
	description = ss.str();
}

// ----------------------------------------------------------------------------
SourceAmplitude::SourceAmplitude(double z) :
		A(z) {
	setDescription();
}

void SourceAmplitude::prepareParticle(ParticleState& p) const {
	Vector3d v =p.getAmplitude();

	p.setAmplitude(v / v.getR() * A);
}

void SourceAmplitude::setDescription() {
	std::stringstream ss;
	ss << "SourceAmplitude: Amplitude A = " << A << "\n";
	description = ss.str();
}




} // namespace radiopropa
