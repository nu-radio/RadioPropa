#include "radiopropa/GridTools.h"
#include "radiopropa/Random.h"
#include "radiopropa/Units.h"

#include <fstream>
#include <sstream>

namespace radiopropa {

void scaleGrid(ref_ptr<ScalarGrid> grid, double a) {
	for (int ix = 0; ix < grid->getNx(); ix++)
		for (int iy = 0; iy < grid->getNy(); iy++)
			for (int iz = 0; iz < grid->getNz(); iz++)
				grid->get(ix, iy, iz) *= a;
}

void scaleGrid(ref_ptr<VectorGrid> grid, double a) {
	for (int ix = 0; ix < grid->getNx(); ix++)
		for (int iy = 0; iy < grid->getNy(); iy++)
			for (int iz = 0; iz < grid->getNz(); iz++)
				grid->get(ix, iy, iz) *= a;
}

Vector3f meanFieldVector(ref_ptr<VectorGrid> grid) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	Vector3f mean(0.);
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				mean += grid->get(ix, iy, iz);
	return mean / Nx / Ny / Nz;
}

double meanFieldStrength(ref_ptr<VectorGrid> grid) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	double mean = 0;
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				mean += grid->get(ix, iy, iz).getR();
	return mean / Nx / Ny / Nz;
}

double meanFieldStrength(ref_ptr<ScalarGrid> grid) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	double mean = 0;
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				mean += grid->get(ix, iy, iz);
	return mean / Nx / Ny / Nz;
}

double rmsFieldStrength(ref_ptr<VectorGrid> grid) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	double sumV2 = 0;
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				sumV2 += grid->get(ix, iy, iz).getR2();
	return std::sqrt(sumV2 / Nx / Ny / Nz);
}

double rmsFieldStrength(ref_ptr<ScalarGrid> grid) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	double sumV2 = 0;
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				sumV2 += pow(grid->get(ix, iy, iz), 2);
	return std::sqrt(sumV2 / Nx / Ny / Nz);
}

double turbulentCorrelationLength(double lMin, double lMax, double alpha) {
	double r = lMin / lMax;
	double a = -alpha - 2;
	return lMax / 2 * (a - 1) / a * (1 - pow(r, a)) / (1 - pow(r, a - 1));
}

void loadGrid(ref_ptr<VectorGrid> grid, std::string filename, double c) {
	std::ifstream fin(filename.c_str(), std::ios::binary);
	if (!fin) {
		std::stringstream ss;
		ss << "load VectorGrid: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}

	// get length of file and compare to size of grid
	fin.seekg(0, fin.end);
	size_t length = fin.tellg() / sizeof(float);
	fin.seekg (0, fin.beg);

	size_t nx = grid->getNx();
	size_t ny = grid->getNy();
	size_t nz = grid->getNz();

	if (length != (3 * nx * ny * nz))
		throw std::runtime_error("loadGrid: file and grid size do not match");

	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				Vector3f &b = grid->get(ix, iy, iz);
				fin.read((char*) &(b.x), sizeof(float));
				fin.read((char*) &(b.y), sizeof(float));
				fin.read((char*) &(b.z), sizeof(float));
				b *= c;
			}
		}
	}
	fin.close();
}

void loadGrid(ref_ptr<ScalarGrid> grid, std::string filename, double c) {
	std::ifstream fin(filename.c_str(), std::ios::binary);
	if (!fin) {
		std::stringstream ss;
		ss << "load ScalarGrid: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}

	// get length of file and compare to size of grid
	fin.seekg(0, fin.end);
	size_t length = fin.tellg() / sizeof(float);
	fin.seekg (0, fin.beg);

	size_t nx = grid->getNx();
	size_t ny = grid->getNy();
	size_t nz = grid->getNz();

	if (length != (nx * ny * nz))
		throw std::runtime_error("loadGrid: file and grid size do not match");

	for (int ix = 0; ix < nx; ix++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int iz = 0; iz < nz; iz++) {
				float &b = grid->get(ix, iy, iz);
				fin.read((char*) &b, sizeof(float));
				b *= c;
			}
		}
	}
	fin.close();
}

void dumpGrid(ref_ptr<VectorGrid> grid, std::string filename, double c) {
	std::ofstream fout(filename.c_str(), std::ios::binary);
	if (!fout) {
		std::stringstream ss;
		ss << "dump VectorGrid: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}
	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				Vector3f b = grid->get(ix, iy, iz) * c;
				fout.write((char*) &(b.x), sizeof(float));
				fout.write((char*) &(b.y), sizeof(float));
				fout.write((char*) &(b.z), sizeof(float));
			}
		}
	}
	fout.close();
}

void dumpGrid(ref_ptr<ScalarGrid> grid, std::string filename, double c) {
	std::ofstream fout(filename.c_str(), std::ios::binary);
	if (!fout) {
		std::stringstream ss;
		ss << "dump ScalarGrid: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}
	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				float b = grid->get(ix, iy, iz) * c;
				fout.write((char*) &b, sizeof(float));
			}
		}
	}
	fout.close();
}

void loadGridFromTxt(ref_ptr<VectorGrid> grid, std::string filename, double c) {
	std::ifstream fin(filename.c_str());
	if (!fin) {
		std::stringstream ss;
		ss << "load VectorGrid: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}
	// skip header lines
	while (fin.peek() == '#')
		fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				Vector3f &b = grid->get(ix, iy, iz);
				fin >> b.x >> b.y >> b.z;
				b *= c;
				if (fin.eof())
					throw std::runtime_error("load VectorGrid: file too short");
			}
		}
	}
	fin.close();
}

void loadGridFromTxt(ref_ptr<ScalarGrid> grid, std::string filename, double c) {
	std::ifstream fin(filename.c_str());
	if (!fin) {
		std::stringstream ss;
		ss << "load ScalarGrid: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}
	// skip header lines
	while (fin.peek() == '#')
		fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				float &b = grid->get(ix, iy, iz);
				fin >> b;
				b *= c;
				if (fin.eof())
					throw std::runtime_error("load ScalarGrid: file too short");
			}
		}
	}
	fin.close();
}

void dumpGridToTxt(ref_ptr<VectorGrid> grid, std::string filename, double c) {
	std::ofstream fout(filename.c_str());
	if (!fout) {
		std::stringstream ss;
		ss << "dump VectorGrid: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}
	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				Vector3f b = grid->get(ix, iy, iz) * c;
				fout << b << "\n";
			}
		}
	}
	fout.close();
}

void dumpGridToTxt(ref_ptr<ScalarGrid> grid, std::string filename, double c) {
	std::ofstream fout(filename.c_str());
	if (!fout) {
		std::stringstream ss;
		ss << "dump ScalarGrid: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}
	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				float b = grid->get(ix, iy, iz) * c;
				fout << b << "\n";
			}
		}
	}
	fout.close();
}

} // namespace radiopropa
