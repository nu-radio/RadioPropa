#include "HepPID/ParticleIDMethods.hh"
#include "radiopropa/Random.h"
#include "radiopropa/magneticLens/ParticleMapsContainer.h"
#include "radiopropa/Units.h"

#include <iostream>
#include <fstream>
namespace radiopropa 
{

ParticleMapsContainer::~ParticleMapsContainer()
{
	for(std::map<int, std::map<int, double*> >::iterator pid_iter = _data.begin(); 
			pid_iter != _data.end(); ++pid_iter) {
		for(std::map<int, double*>::iterator frequency_iter = pid_iter->second.begin();
			frequency_iter != pid_iter->second.end(); ++frequency_iter) {
				delete[] (frequency_iter->second);
		}
	}
}

int ParticleMapsContainer::frequency2Idx(double frequency) const
{
	double lE = log10(frequency / eV);
	return int((lE - _bin0lowerEdge) / _deltaLogE);
}

double ParticleMapsContainer::idx2Frequency(int idx) const
{
	return pow(10, idx * _deltaLogE + _bin0lowerEdge + _deltaLogE / 2) * eV;
}

		
double* ParticleMapsContainer::getMap(const int particleId, double frequency)
{
	_weightsUpToDate = false;
	if (_data.find(particleId) == _data.end())
	{
		std::cerr << "No map for ParticleID " << particleId << std::endl;
		return NULL;
	}
	int frequencyIdx	= frequency2Idx(frequency);
	if (_data[particleId].find(frequencyIdx) == _data[particleId].end())
	{
		std::cerr << "No map for ParticleID and frequency" << frequency / eV << " eV" << std::endl;
		return NULL;
	}
	return _data[particleId][frequency2Idx(frequency)];
}
			
			
void ParticleMapsContainer::addParticle(const int particleId, double frequency, double galacticLongitude, double galacticLatitude, double weight)
{
	_weightsUpToDate = false;
	if (_data.find(particleId) == _data.end())
	{
		map<int, double*> M;
		_data[particleId] = M;
	}

	int frequencyIdx	= frequency2Idx(frequency);
	if (_data[particleId].find(frequencyIdx) == _data[particleId].end())
	{
		_data[particleId][frequencyIdx] = new double[_pixelization.getNumberOfPixels()];
		std::fill(_data[particleId][frequencyIdx], _data[particleId][frequencyIdx] + _pixelization.getNumberOfPixels(), 0);
	}

	uint32_t pixel = _pixelization.direction2Pix(galacticLongitude, galacticLatitude);
	_data[particleId][frequencyIdx][pixel] += weight;
}


void ParticleMapsContainer::addParticle(const int particleId, double frequency, const Vector3d &p, double weight)
{
	double galacticLongitude = atan2(-p.y, -p.x);
	double galacticLatitude =	M_PI / 2 - acos(-p.z / p.getR());
	addParticle(particleId, frequency, galacticLongitude, galacticLatitude, weight);
}


// returns a vector of all particle ids in th
std::vector<int> ParticleMapsContainer::getParticleIds()
{
	std::vector<int> ids;
	for(std::map<int, std::map<int, double*> >::iterator pid_iter = _data.begin(); 
			pid_iter != _data.end(); ++pid_iter) 
	{
		ids.push_back(pid_iter->first);
	}
	return ids;
}


std::vector<double> ParticleMapsContainer::getEnergies(int pid)
{
	std::vector<double> energies;
	if (_data.find(pid) != _data.end())
	{
		for(std::map<int, double*>::iterator iter = _data[pid].begin(); 
			iter != _data[pid].end(); ++iter) 
		{
			energies.push_back( idx2Frequency(iter->first) / eV );
		}
	}
	return energies;
}


void ParticleMapsContainer::applyLens(MagneticLens &lens)
{
	// if lens is normalized, this should not be necessary.
	_weightsUpToDate = false;

	for(std::map<int, std::map<int, double*> >::iterator pid_iter = _data.begin(); 
			pid_iter != _data.end(); ++pid_iter) {
		for(std::map<int, double*>::iterator frequency_iter = pid_iter->second.begin();
			frequency_iter != pid_iter->second.end(); ++frequency_iter) {
		//	// transform only nuclei
			double frequency = idx2Frequency(frequency_iter->first);
			int chargeNumber = HepPID::Z(pid_iter->first);
			if (chargeNumber != 0 && lens.rigidityCovered(frequency / chargeNumber))
			{
				lens.transformModelVector(frequency_iter->second, frequency / chargeNumber);
			}
			else
			{ // still normalize the vectors 
				for(size_t j=0; j< _pixelization.getNumberOfPixels() ; j++)
				{
					frequency_iter->second[j]/=lens.getNorm();
				}
			}
		}
	}
}


void ParticleMapsContainer::_updateWeights()
{
	if (_weightsUpToDate)
		return;

	for(std::map<int, std::map<int, double*> >::iterator pid_iter = _data.begin(); 
			pid_iter != _data.end(); ++pid_iter) 
	{
		_weightsPID[pid_iter->first] = 0;

		for(std::map<int, double*>::iterator frequency_iter = pid_iter->second.begin();
			frequency_iter != pid_iter->second.end(); ++frequency_iter) 
		{

			_weights_pidFrequency[pid_iter->first][frequency_iter->first] = 0;
			for(size_t j=0; j< _pixelization.getNumberOfPixels() ; j++)
			{
				_weights_pidFrequency[pid_iter->first][frequency_iter->first] +=frequency_iter->second[j];
					
				_weightsPID[pid_iter->first]+=frequency_iter->second[j];
			}
		_sumOfWeights+=_weights_pidFrequency[pid_iter->first][frequency_iter->first];
		}
	}
	_weightsUpToDate = true;
}


void ParticleMapsContainer::getRandomParticles(size_t N, vector<int> &particleId, 
	vector<double> &frequency, vector<double> &galacticLongitudes,
	vector<double> &galacticLatitudes)
{
	_updateWeights();

	particleId.resize(N);
	frequency.resize(N);
	galacticLongitudes.resize(N);
	galacticLatitudes.resize(N);

	for(size_t i=0; i< N; i++)
	{
		//get particle
		double r = Random::instance().rand() * _sumOfWeights;
		std::map<int, double>::iterator iter = _weightsPID.begin();
		while ((r-= iter->second) > 0)
		{
		 ++iter; 
		}
		particleId[i] = iter->first;
	
		//get frequency
		r = Random::instance().rand() * iter->second;
		iter = _weights_pidFrequency[particleId[i]].begin();
		while ((r-= iter->second) > 0)
		{
		 ++iter; 
		}
		frequency[i] = idx2Frequency(iter->first) / eV;

		placeOnMap(particleId[i], frequency[i] * eV, galacticLongitudes[i], galacticLatitudes[i]);
	}
}


bool ParticleMapsContainer::placeOnMap(int pid, double frequency, double &galacticLongitude, double &galacticLatitude)
{
	_updateWeights();

	if (_data.find(pid) == _data.end())
	{
		return false;
	}
	int frequencyIdx	= frequency2Idx(frequency);
	if (_data[pid].find(frequencyIdx) == _data[pid].end())
	{
		return false;
	}

	double r = Random::instance().rand() * _weights_pidFrequency[pid][frequencyIdx];

	for(size_t j=0; j< _pixelization.getNumberOfPixels(); j++)
	{
		r-= _data[pid][frequencyIdx][j];
		if (r <=0)
		{
			_pixelization.getRandomDirectionInPixel(j, galacticLongitude, galacticLatitude );
			return true;
		}
	}
	return false;
}


void ParticleMapsContainer::forceWeightUpdate()
{
	_weightsUpToDate = false;
}

} // namespace parsec
