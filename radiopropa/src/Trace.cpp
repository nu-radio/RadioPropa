#include "radiopropa/Trace.h"
#include "radiopropa/Units.h"
#include <math.h>

namespace radiopropa {

Trace::Trace():
	samplingRate(NULL),
	traceStartTime(0.0),
	frequencySpectrum(NULL) {
}

std::vector<std::complex<double>> Trace::getFrequencySpectrum() const {
	return frequencySpectrum;
}


std::vector<double> Trace::getFrequencySpectrum_real() const {
	std::vector<double> spectrum_real;
	for (int i=0; i < frequencySpectrum.size(); i++) {
		spectrum_real.push_back(std::real(frequencySpectrum[i]));	}
	return spectrum_real;
}

std::vector<double> Trace::getFrequencySpectrum_imag() const {
	std::vector<double> spectrum_imag;
	for (int i=0; i < frequencySpectrum.size(); i++) {
		spectrum_imag.push_back(std::imag(frequencySpectrum[i]));
	}
	return spectrum_imag;
}

std::vector<double> Trace::getFrequencies() const {
	int length_t = this->getNumberOfSamples();
	int length_f = (length_t / 2) + 1;
	std::vector<double> frequencies;
	double f0 = 1.0 / samplingRate;
	for (int i=0; i < length_f; i++) {
		double f = i * f0;
		frequencies.push_back(f);
	}
	return frequencies;
}

double Trace::getSamplingRate() const {
	return samplingRate;
}

int Trace::getNumberOfSamples() const {
	int length = (frequencySpectrum.size() - 1) * 2;    
	return length;
}

double Trace::getTraceStartTime() const {
	return traceStartTime;
}

void Trace::setFrequencySpectrum(std::vector<double> spectrum_real, std::vector<double> spectrum_imag, double rate) {
	if (spectrum_real.size() != spectrum_imag.size()) {
		std::cout << "Real and imaginary part of spectrum should have the same length";
		throw std::exception();
	}
	for (int i=0; i < spectrum_real.size(); i++) {
		std::complex<double> A(spectrum_real[i], spectrum_imag[i]);
		frequencySpectrum.push_back(A);
	}
	samplingRate = rate;
}

void Trace::setFrequencySpectrum(std::vector<std::complex<double>> spectrum, double rate) {
	frequencySpectrum = spectrum;
    samplingRate = rate;
}

void Trace::setTraceStartTime(double start_time) {
	traceStartTime = start_time;
}

void Trace::addTraceStartTime(double start_time) {
	traceStartTime += start_time;
}

void Trace::applyTimeShift(double delta_t, bool silent) {
	if ((delta_t > 0.1 * this->getNumberOfSamples() / this->getSamplingRate()) and not silent) {
         std::cout << "Trace is shifted by more than 10% of its length";
	}
    std::vector<std::complex<double>> new_spectrum;
    for (int i=0; i < frequencySpectrum.size(); i++) {
    	std::complex<double> A = frequencySpectrum[i];
    	double freq = this->getFrequencies()[i];
    	double phase = 2.0 * M_PI * delta_t * freq;
    	std::complex<double> phase_vector = std::polar(1.0, phase);
    	new_spectrum.push_back(A * phase_vector);
    }
    this->setFrequencySpectrum(new_spectrum, samplingRate);
}




ElectricField::ElectricField():
	samplingRate(NULL),
	traceStartTime(0.0),
	r(),
	theta(),
	phi() {
}

std::vector<std::vector<std::complex<double>>> ElectricField::getFrequencySpectrum() const {
	std::vector<std::vector<std::complex<double>>> spectrum;
	spectrum.push_back(r.getFrequencySpectrum());
	spectrum.push_back(theta.getFrequencySpectrum());
	spectrum.push_back(phi.getFrequencySpectrum());
	return spectrum;
}

std::vector<std::vector<double>> ElectricField::getFrequencySpectrum_real() const {
	std::vector<std::vector<double>> spectrum_real;
	spectrum_real.push_back(r.getFrequencySpectrum_real());
	spectrum_real.push_back(theta.getFrequencySpectrum_real());
	spectrum_real.push_back(phi.getFrequencySpectrum_real());
	return spectrum_real;
}

std::vector<std::vector<double>> ElectricField::getFrequencySpectrum_imag() const {
	std::vector<std::vector<double>> spectrum_imag;
	spectrum_imag.push_back(r.getFrequencySpectrum_imag());
	spectrum_imag.push_back(theta.getFrequencySpectrum_imag());
	spectrum_imag.push_back(phi.getFrequencySpectrum_imag());
	return spectrum_imag;
}

std::vector<double> ElectricField::getFrequencies() const {
	if ((not (r.getSamplingRate() == theta.getSamplingRate())) or (not (r.getSamplingRate() == phi.getSamplingRate()))) {
		std::cout << "Sampling rates of r, theta and phi component should be the same";
		throw std::exception();
	}
	if ((not (r.getFrequencies().size() == theta.getFrequencies().size())) or (not (r.getFrequencies().size() == phi.getFrequencies().size()))) {
		std::cout << "The r, theta and phi component should have the same amount of frequencies be the same";
		throw std::exception();
	}

	return r.getFrequencies();
}

double ElectricField::getSamplingRate() const {
	return samplingRate;
}

int ElectricField::getNumberOfSamples() const {
	int length_r = (r.getFrequencySpectrum().size() - 1) * 2;
	int length_theta = (theta.getFrequencySpectrum().size() - 1) * 2;
	int length_phi = (phi.getFrequencySpectrum().size() - 1) * 2;   
	if ((length_r != length_theta) and (length_r != length_phi)) {
		throw std::exception();
	}
	return length_r;
}

double ElectricField::getTraceStartTime() const {
	return traceStartTime;
}

void ElectricField::setFrequencySpectrumComponent(std::string component, std::vector<double> spectrum_real, std::vector<double> spectrum_imag, double rate) {
	if (spectrum_real.size() != spectrum_imag.size()) {
		std::cout << "Real and imaginary part of spectrum should have the same length";
		throw std::exception();
	} 

	if (component == "r") {
		r.setFrequencySpectrum(spectrum_real, spectrum_imag, rate);
	} else if (component == "theta") {
		theta.setFrequencySpectrum(spectrum_real, spectrum_imag, rate);
	} else if (component == "phi") {
		phi.setFrequencySpectrum(spectrum_real, spectrum_imag, rate);
	} else {
		std::cout << "Component should be r, theta or phi";
		throw std::exception();
	}
	samplingRate = rate;
}

void ElectricField::setFrequencySpectrumComponent(std::string component, std::vector<std::complex<double>> spectrum, double rate) {
	if (component == "r") {
		r.setFrequencySpectrum(spectrum, rate);
	} else if (component == "theta") {
		theta.setFrequencySpectrum(spectrum, rate);
	} else if (component == "phi") {
		phi.setFrequencySpectrum(spectrum, rate);
	} else {
		std::cout << "Component should be r, theta or phi";
		throw std::exception();
	}
	samplingRate = rate;
}

void ElectricField::setFrequencySpectrum(std::vector<double> spectrum_r_real, std::vector<double> spectrum_r_imag, 
							 	  std::vector<double> spectrum_theta_real, std::vector<double> spectrum_theta_imag, 
								  std::vector<double> spectrum_phi_real, std::vector<double> spectrum_phi_imag, 
								  double rate) {
	if ((not (spectrum_r_real.size() == spectrum_r_imag.size())) or (not (spectrum_theta_real.size() == spectrum_theta_imag.size())) or (not (spectrum_phi_real.size() == spectrum_phi_imag.size()))) {
		std::cout << "Real and imaginary part of spectrum should have the same length";
		throw std::exception();
	} 

	if ((not (spectrum_r_real.size() == spectrum_theta_real.size())) or (not (spectrum_r_real.size() == spectrum_phi_real.size()))) {
		std::cout << "The r, theta and phi component of spectrum should have the same length";
		throw std::exception();
	}

	r.setFrequencySpectrum(spectrum_r_real, spectrum_r_imag, rate);
	theta.setFrequencySpectrum(spectrum_theta_real, spectrum_theta_imag, rate);
	phi.setFrequencySpectrum(spectrum_phi_real, spectrum_phi_imag, rate);
	samplingRate = rate;
}

void ElectricField::setFrequencySpectrum(std::vector<std::complex<double>> spectrum_r, 
							  	  std::vector<std::complex<double>> spectrum_theta,
								  std::vector<std::complex<double>> spectrum_phi,
								  double rate) {
	if ((not (spectrum_r.size() == spectrum_theta.size())) or (not (spectrum_r.size() == spectrum_phi.size()))) {
		std::cout << "The r, theta and phi component should have the same length";
		throw std::exception();
	}
	r.setFrequencySpectrum(spectrum_r, rate);
	theta.setFrequencySpectrum(spectrum_theta, rate);
	phi.setFrequencySpectrum(spectrum_phi, rate);
	samplingRate = rate;
}

void ElectricField::setTraceStartTime(double start_time) {
	traceStartTime = start_time;
	r.setTraceStartTime(start_time);
	theta.setTraceStartTime(start_time);
	phi.setTraceStartTime(start_time);
}

void ElectricField::addTraceStartTime(double start_time) {
	traceStartTime += start_time;
	r.addTraceStartTime(start_time);
	theta.addTraceStartTime(start_time);
	phi.addTraceStartTime(start_time);
}

void ElectricField::applyTimeShift(double delta_t, bool silent) {
    r.applyTimeShift(delta_t, silent);
    theta.applyTimeShift(delta_t, silent);
    phi.applyTimeShift(delta_t, silent);
}

} // namespace radiopropa