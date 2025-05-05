#ifndef RADIOPROPA_TRACE_H
#define RADIOPROPA_TRACE_H

#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include <complex>
#include <tuple>
#include "radiopropa/Referenced.h"

namespace radiopropa {

class Trace: public Referenced{
protected:
	double samplingRate;
	double traceStartTime;
	std::vector<std::complex<double>> frequencySpectrum;

public:
	Trace();

	std::vector<std::complex<double>> getFrequencySpectrum() const;
	std::vector<double> getFrequencySpectrum_real() const;
	std::vector<double> getFrequencySpectrum_imag() const;
	std::vector<double> getFrequencies() const;
	double getSamplingRate() const;
	int getNumberOfSamples() const;
    double getTraceStartTime() const;

	void setFrequencySpectrum(std::vector<double> spectrum_real, std::vector<double> spectrum_imag, double rate);
	void setFrequencySpectrum(std::vector<std::complex<double>> spectrum, double rate);
	void setTraceStartTime(double start_time);

    void addTraceStartTime(double start_time);
    void applyTimeShift(double delta_t, bool silent = false);
	
};

class ElectricField: public Referenced {
protected:
	double samplingRate;
	double traceStartTime;

	Trace r;
	Trace theta;
	Trace phi;

public:
	ElectricField();
	//ElectricField(const ElectricField &efield);

	std::vector<Trace> getTraces() const;
	std::vector<std::vector<std::complex<double>>> getFrequencySpectrum() const;
	std::vector<std::vector<double>> getFrequencySpectrum_real() const;
	std::vector<std::vector<double>> getFrequencySpectrum_imag() const;
	std::vector<double> getFrequencies() const;
	double getSamplingRate() const;
	int getNumberOfSamples() const;
    double getTraceStartTime() const;

	void setFrequencySpectrumComponent(std::string component, std::vector<double> spectrum_real, std::vector<double> spectrum_imag, double rate);
	void setFrequencySpectrumComponent(std::string component, std::vector<std::complex<double>> spectrum, double rate);
	void setFrequencySpectrum(std::vector<double> spectrum_r_real, std::vector<double> spectrum_r_imag, 
							  std::vector<double> spectrum_theta_real, std::vector<double> spectrum_theta_imag, 
							  std::vector<double> spectrum_phi_real, std::vector<double> spectrum_phi_imag, 
							  double rate);
	void setFrequencySpectrum(std::vector<std::complex<double>> spectrum_r, 
							  std::vector<std::complex<double>> spectrum_theta,
							  std::vector<std::complex<double>> spectrum_phi,
							  double rate);
	void setTraceStartTime(double start_time);

    void addTraceStartTime(double start_time);
    void applyTimeShift(double delta_t, bool silent = false);

};

} // namespace radiopropa

#endif // RADIOPROPA_TRACE_H
