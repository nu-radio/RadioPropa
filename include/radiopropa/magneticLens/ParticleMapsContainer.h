#ifndef PARTICLEMAPSCONTAINER_HH
#define PARTICLEMAPSCONTAINER_HH 

#include <map>
#include <vector>
#include "radiopropa/magneticLens/Pixelization.h"
#include "radiopropa/magneticLens/MagneticLens.h"

#include "radiopropa/Vector3.h"

namespace radiopropa{

/// Container for particlemaps
/// The maps are stored with discrete energies on a logarithmic scale. The
/// default frequency width is 0.02 with an frequency bin from 10**17.99 - 10**18.01
/// eV
class ParticleMapsContainer 
{
	private:
    std::map< int , std::map <int , double*> > _data;				
		Pixelization _pixelization;
    double _deltaLogE;
    double _bin0lowerEdge;

		// get the bin number of the frequency
		int frequency2Idx(double frequency) const;
		double idx2Frequency(int idx) const;


		// weights of the particles
		double _sumOfWeights;
		std::map< int , double > _weightsPID;				
		std::map< int , map<int, double> > _weights_pidFrequency;				

		// lazy update of weights 
		bool _weightsUpToDate;
		void _updateWeights();
  public:

		ParticleMapsContainer(double deltaLogE = 0.02, double bin0lowerEdge = 17.99) : _deltaLogE(deltaLogE), _bin0lowerEdge(bin0lowerEdge), _pixelization(6), _weightsUpToDate(false), _sumOfWeights(0)
		{
		}

    ~ParticleMapsContainer();

    size_t getNumberOfPixels()
    {
      return _pixelization.getNumberOfPixels();
    }

		/// returns the map for the particleId with the given frequency,. frequency in
		/// Joule
		double *getMap(const int particleId, double frequency);
			
		/// adds a particle to the map container
    /// particleId is HEP particleId, frequency [Joule], galacticLongitude and
    /// galacticLatitude in [rad]
		void addParticle(const int particleId, double frequency, double galacticLongitude, double galacticLatitude, double weight = 1);
			
		void addParticle(const int particleId, double frequency, const Vector3d &v, double weight = 1);
		
		// returns a vector of all particle ids in th
		std::vector<int> getParticleIds();

		// the energies are in eV
		std::vector<double> getEnergies(int pid);

		void applyLens(MagneticLens &lens);;

		// frequency in eV , galacticLongitude in rad [-pi ... pi], galacticLatitudes in rad [-pi/2 ... pi/2]
		void getRandomParticles(size_t N, vector<int> &particleId, 
			vector<double> &frequency, vector<double> &galacticLongitudes,
			vector<double> &galacticLatitudes);

		// places a cosmic ray with given PID and frequency according to the
		// probability maps. Returns false if not possible.
		bool placeOnMap(int pid, double frequency, double &galacticLongitude, double &galacticLatitude);

		// force weight update prior to get random particles. Only necessary when
		// reusing pointer to maps after calculating weights
		void forceWeightUpdate();

		double getSumOfWeights()
		{
			if (!_weightsUpToDate)
				_updateWeights();
			return _sumOfWeights;
		}

		double getWeight(int pid, double frequency)
		{
			if (!_weightsUpToDate)
				_updateWeights();
			return _weights_pidFrequency[pid][frequency2Idx(frequency)];				
		}
};


/*
Python code that we have to write







*/




} // namespace parsec

#endif // PARTICLEMAPSCONTAINER_HH
