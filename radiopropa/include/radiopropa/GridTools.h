#ifndef CRPROPA_GRIDTOOLS_H
#define CRPROPA_GRIDTOOLS_H

#include "radiopropa/Grid.h"
#include <string>

/**
 @file
 @brief Grid related functions: load, dump, save, create turbulent field ...

 This file contains a number of functions related to scalar and vector grids (Grid.h).

 Dump/load functions are available for saving/loading grids to/from and binary and plain text files.
 In the files the grid points are stored from (0, 0, 0) to (Nx, Ny, Nz) with the z-index changing the fastest.
 Vector components are stored per grid point in xyz-order.
 In case of plain-text files the vector components are separated by a blank or tab and grid points are stored one per line.
 All functions offer a conversion factor that is multiplied to all values.
 */

namespace radiopropa {

/** Evaluate the mean vector of all grid points */
Vector3f meanFieldVector(ref_ptr<VectorGrid> grid);

/** Evaluate the mean of all grid points */
double meanFieldStrength(ref_ptr<ScalarGrid> grid);
/** Evaluate the mean of all grid points */
double meanFieldStrength(ref_ptr<VectorGrid> grid);

/** Evaluate the RMS of all grid points */
double rmsFieldStrength(ref_ptr<ScalarGrid> grid);
/** Evaluate the RMS of all grid points */
double rmsFieldStrength(ref_ptr<VectorGrid> grid);

/** Multiply all grid values by a given factor */
void scaleGrid(ref_ptr<ScalarGrid> grid, double a);
/** Multiply all grid values by a given factor */
void scaleGrid(ref_ptr<VectorGrid> grid, double a);

/** Analytically calculate the correlation length of a turbulent field */
double turbulentCorrelationLength(double lMin, double lMax,
		double alpha = (-11./3.));

/** Load a VectorGrid from a binary file with single precision */
void loadGrid(ref_ptr<VectorGrid> grid, std::string filename,
		double conversion = 1);

/** Load a ScalarGrid from a binary file with single precision */
void loadGrid(ref_ptr<ScalarGrid> grid, std::string filename,
		double conversion = 1);

/** Dump a VectorGrid to a binary file */
void dumpGrid(ref_ptr<VectorGrid> grid, std::string filename,
		double conversion = 1);

/** Dump a ScalarGrid to a binary file with single precision */
void dumpGrid(ref_ptr<ScalarGrid> grid, std::string filename,
		double conversion = 1);

/** Load a VectorGrid grid from a plain text file */
void loadGridFromTxt(ref_ptr<VectorGrid> grid, std::string filename,
		double conversion = 1);

/** Load a ScalarGrid from a plain text file */
void loadGridFromTxt(ref_ptr<ScalarGrid> grid, std::string filename,
		double conversion = 1);

/** Dump a VectorGrid to a plain text file */
void dumpGridToTxt(ref_ptr<VectorGrid> grid, std::string filename,
		double conversion = 1);

/** Dump a ScalarGrid to a plain text file */
void dumpGridToTxt(ref_ptr<ScalarGrid> grid, std::string filename,
		double conversion = 1);

} // namespace radiopropa

#endif // CRPROPA_GRIDTOOLS_H
