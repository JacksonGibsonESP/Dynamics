#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <random>
#include <chrono>
#include "annealing.h"
#include "GRS.h"
#include "NelderMead.h"
#include "fitting.h"
using namespace std;

bool toofar(Atom part, Atom part2, double cutoff) {
	if (part.x == part2.x && part.y == part2.y && part.z == part2.z)
		return true;
	return sqrt((part2.x - part.x) * (part2.x - part.x)
		+ (part2.y - part.y) * (part2.y - part.y)
		+ (part2.z - part.z) * (part2.z - part.z)) > cutoff ? true : false;
}

double full_energy(vector <vector<Atom>> &contains, int n_x, int n_y, int n_z, double cutoff, RGL potentials[3], double translation[9], bool surface) {

	double E_full = 0;
	double bulk_size_x = translation[0] + translation[3] + translation[6];
	double bulk_size_y = translation[1] + translation[4] + translation[7];
	double bulk_size_z = translation[2] + translation[5] + translation[8];

	//#pragma omp parallel for reduction(+:E_full) num_threads(3) 
	for (int i = 0; i < n_x * n_y * n_z; ++i) {
		for (Atom& particle : contains[i]) {
			int atomcell_x = i % n_x;
			int atomcell_y = i / n_x % n_y;
			int atomcell_z = i / (n_x * n_y);
			double E_b_temp = 0;
			double E_r_temp = 0;
			for (int dx = -1; dx <= 1; ++dx) {
				for (int dy = -1; dy <= 1; ++dy) {
					for (int dz = -1; dz <= 1; ++dz) {
						for (Atom particle2 : contains[
							(atomcell_x + n_x + dx) % n_x +
								(atomcell_y + n_y + dy) % n_y * n_x +
								(atomcell_z + n_z + dz) % n_z * n_y * n_x]) {
							//check for translation
							if (abs(particle2.x - particle.x) > bulk_size_x / 2) { //more than one cell far
								if (particle2.x > particle.x) { //move to neighbor cell
									particle2.x -= translation[0];
									particle2.y -= translation[1];
									particle2.z -= translation[2];
								}
								else {
									particle2.x += translation[0];
									particle2.y += translation[1];
									particle2.z += translation[2];
								}
							}
							if (abs(particle2.y - particle.y) > bulk_size_y / 2) { //more than one cell far
								if (particle2.y > particle.y) { //move to neighbor cell
									particle2.x -= translation[3];
									particle2.y -= translation[4];
									particle2.z -= translation[5];
								}
								else {
									particle2.x += translation[3];
									particle2.y += translation[4];
									particle2.z += translation[5];
								}
							}
							if (!surface) {
								if (abs(particle2.z - particle.z) > bulk_size_z / 2) { //more than one cell far
									if (particle2.z > particle.z) { //move to neighbor cell
										particle2.x -= translation[6];
										particle2.y -= translation[7];
										particle2.z -= translation[8];
									}
									else {
										particle2.x += translation[6];
										particle2.y += translation[7];
										particle2.z += translation[8];
									}
								}
							}
							if (toofar(particle, particle2, cutoff)) {
								continue;
							}

							int choice = 0;
							if (particle.impurity && particle2.impurity) {
								choice = 2; //A-A
							} else if (!particle.impurity && !particle2.impurity) {
								choice = 0; //B-B
							} else {
								choice = 1; //A-B
							}

							//distance
							double r = sqrt((particle2.x - particle.x) * (particle2.x - particle.x)
								+ (particle2.y - particle.y) * (particle2.y - particle.y)
								+ (particle2.z - particle.z) * (particle2.z - particle.z));

							//repulsive and bond energies
							E_r_temp += (potentials[choice].A_1 * (r - potentials[choice].r0) + potentials[choice].A_0)
								* exp(-potentials[choice].p * (r / potentials[choice].r0 - 1));
							E_b_temp += potentials[choice].s * potentials[choice].s
								* exp(-2 * potentials[choice].q * (r / potentials[choice].r0 - 1));
						}
					}
				}
			}
			E_full += E_r_temp - sqrt(E_b_temp);
		}
	}
	return E_full;
}

void stretch(vector <vector<Atom>> &contains, int n_x, int n_y, int n_z, double alpha) {
	for (int i = 0; i < n_x * n_y * n_z; ++i) {
		for (Atom& particle : contains[i]) {
			particle.x = particle.x + alpha * particle.x;
			particle.y = particle.y + alpha * particle.y;
			particle.z = particle.z + alpha * particle.z;
		}
	}
}

void elastic_C11(vector <vector<Atom>> &contains, int n_x, int n_y, int n_z, double alpha) {
	for (int i = 0; i < n_x * n_y * n_z; ++i) {
		for (Atom& particle : contains[i]) {
			particle.x = particle.x + alpha * particle.x;
			particle.y = particle.y + alpha * particle.y;
		}
	}
}

void elastic_C12(vector <vector<Atom>> &contains, int n_x, int n_y, int n_z, double alpha) {
	for (int i = 0; i < n_x * n_y * n_z; ++i) {
		for (Atom& particle : contains[i]) {
			particle.x = particle.x + alpha * particle.x;
			particle.y = particle.y - alpha * particle.y;
		}
	}
}

void elastic_C44(vector <vector<Atom>> &contains, int n_x, int n_y, int n_z, double alpha) {
	for (int i = 0; i < n_x * n_y * n_z; ++i) {
		for (Atom& particle : contains[i]) {
			double tmp_x = particle.x + alpha * particle.y;
			double tmp_y = alpha * particle.x + particle.y;
			particle.x = tmp_x;
			particle.y = tmp_y;
			particle.z = particle.z / (1 - alpha * alpha);
		}
	}
}

void dump(string filename, int size, vector <vector<Atom>> &contains) {
	ofstream out(filename);
	out << size << '\n';
	out << "one frame\n";
	for (unsigned int i = 0; i < contains.size(); ++i) {
		for (unsigned int j = 0; j < contains[i].size(); ++j) {
			out << contains[i][j].x << ' ' << contains[i][j].y << ' ' << contains[i][j].z << '\n';
		}
	}
	out.close();
}

int convert(string filename, int n_x, int n_y, int n_z, double lattice_constant, vector <vector<Atom>> &out, int bulk_atom_count) {

	out.resize(n_x * n_y * n_z);
	vector<Atom> init;

	ifstream fin(filename, ifstream::in);

	if (!fin.is_open()) {
		std::cout << "File can not be opened\n";
		exit(0);
	}

	int size = 0;
	fin >> size;
	
	if (bulk_atom_count != 0 && size > bulk_atom_count) {
		out.resize(n_x * n_y * (n_z + 1));
	}

	while (fin.good()) {
		Atom tmp;
		fin >> tmp.impurity;
		fin >> tmp.x;
		fin >> tmp.y;
		fin >> tmp.z;
		if (tmp.impurity) {
			tmp.mass = 58.6934;
		}
		init.push_back(tmp);
	}
	fin.close();

	for (unsigned int i = 0; i < init.size(); ++i) {
		Atom tmp = init[i];
		int x = int(tmp.x / lattice_constant); // step by lattice constant
		int y = int(tmp.y / lattice_constant);
		int z = int(tmp.z / lattice_constant);
		out[z * n_y * n_x + y * n_x + x].push_back(tmp);
	}
	return size;
}

void calc_E_coh(CalcData &calcData) {
	calcData.tr.make_it_pure();
	calcData.E_full = full_energy(calcData.bulk, calcData.n_x, calcData.n_y, calcData.n_z, calcData.cutoff, calcData.potentials, calcData.tr.translation, false);
	calcData.E_coh = calcData.E_full / calcData.bulk_atom_count;
}

void calc_B(CalcData &calcData) {
	vector <vector<Atom>> curr = calcData.bulk;
	stretch(curr, calcData.n_x, calcData.n_y, calcData.n_z, calcData.alpha);
	calcData.tr.alpha = calcData.alpha;
	calcData.tr.make_it_B();
	double E_full_stretched = full_energy(curr, calcData.n_x, calcData.n_y, calcData.n_z, calcData.cutoff, calcData.potentials, calcData.tr.translation, false);
	curr = calcData.bulk;

	stretch(curr, calcData.n_x, calcData.n_y, calcData.n_z, -calcData.alpha);
	calcData.tr.alpha = -calcData.alpha;
	calcData.tr.make_it_B();
	double E_full_compressed = full_energy(curr, calcData.n_x, calcData.n_y, calcData.n_z, calcData.cutoff, calcData.potentials, calcData.tr.translation, false);
	curr = calcData.bulk;
	calcData.B = (2 * (E_full_stretched - 2 * calcData.E_full + E_full_compressed) * calcData.qsi) / (9 * calcData.V_0 * calcData.alpha * calcData.alpha);
}

void calc_C11_C12(CalcData &calcData) {
	vector <vector<Atom>> curr = calcData.bulk;
	elastic_C11(curr, calcData.n_x, calcData.n_y, calcData.n_z, calcData.alpha);
	calcData.tr.alpha = calcData.alpha;
	calcData.tr.make_it_C11();
	double E_full_stretched_C11 = full_energy(curr, calcData.n_x, calcData.n_y, calcData.n_z, calcData.cutoff, calcData.potentials, calcData.tr.translation, false);
	curr = calcData.bulk;

	elastic_C11(curr, calcData.n_x, calcData.n_y, calcData.n_z, -calcData.alpha);
	calcData.tr.alpha = -calcData.alpha;
	calcData.tr.make_it_C11();
	double E_full_compressed_C11 = full_energy(curr, calcData.n_x, calcData.n_y, calcData.n_z, calcData.cutoff, calcData.potentials, calcData.tr.translation, false);
	curr = calcData.bulk;

	elastic_C12(curr, calcData.n_x, calcData.n_y, calcData.n_z, calcData.alpha);
	calcData.tr.alpha = calcData.alpha;
	calcData.tr.make_it_C12();
	double E_full_stretched_C12 = full_energy(curr, calcData.n_x, calcData.n_y, calcData.n_z, calcData.cutoff, calcData.potentials, calcData.tr.translation, false);
	curr = calcData.bulk;

	elastic_C12(curr, calcData.n_x, calcData.n_y, calcData.n_z, -calcData.alpha);
	calcData.tr.alpha = -calcData.alpha;
	calcData.tr.make_it_C12();
	double E_full_compressed_C12 = full_energy(curr, calcData.n_x, calcData.n_y, calcData.n_z, calcData.cutoff, calcData.potentials, calcData.tr.translation, false);
	curr = calcData.bulk;

	calcData.C11 = ((E_full_stretched_C11 - 2 * calcData.E_full + E_full_compressed_C11) / (calcData.alpha * calcData.alpha) + (E_full_stretched_C12 - 2 * calcData.E_full + E_full_compressed_C12) / (calcData.alpha * calcData.alpha)) *calcData.qsi / (2 * calcData.V_0);

	calcData.C12 = ((E_full_stretched_C11 - 2 * calcData.E_full + E_full_compressed_C11) / (calcData.alpha * calcData.alpha) - (E_full_stretched_C12 - 2 * calcData.E_full + E_full_compressed_C12) / (calcData.alpha * calcData.alpha)) * calcData.qsi / (2 * calcData.V_0);
}

void calc_C44(CalcData &calcData) {
	vector <vector<Atom>> curr = calcData.bulk;
	elastic_C44(curr, calcData.n_x, calcData.n_y, calcData.n_z, calcData.alpha);
	calcData.tr.alpha = calcData.alpha;
	calcData.tr.make_it_C44();
	double E_full_stretched_C44 = full_energy(curr, calcData.n_x, calcData.n_y, calcData.n_z, calcData.cutoff, calcData.potentials, calcData.tr.translation, false);
	curr = calcData.bulk;

	elastic_C44(curr, calcData.n_x, calcData.n_y, calcData.n_z, -calcData.alpha);
	calcData.tr.alpha = -calcData.alpha;
	calcData.tr.make_it_C44();
	double E_full_compressed_C44 = full_energy(curr, calcData.n_x, calcData.n_y, calcData.n_z, calcData.cutoff, calcData.potentials, calcData.tr.translation, false);
	curr = calcData.bulk;

	calcData.C44 = ((E_full_stretched_C44 - 2 * calcData.E_full + E_full_compressed_C44)  * calcData.qsi) / (2 * calcData.V_0 * calcData.alpha * calcData.alpha);
}

void calc_E_sol(CalcData &calcData) {
	calcData.tr.make_it_pure();
	double E_AB = full_energy(calcData.imp, calcData.n_x, calcData.n_y, calcData.n_z, calcData.cutoff, calcData.potentials, calcData.tr.translation, false);
	double E_B = full_energy(calcData.bulk, calcData.n_x, calcData.n_y, calcData.n_z, calcData.cutoff, calcData.potentials, calcData.tr.translation, false);
	calcData.E_sol = E_AB - E_B - calcData.E_coh_imp + calcData.E_coh;
}

void calc_E_in_dim(CalcData &calcData) {
	calcData.tr.make_it_pure();
	double E_dim_in_surf = full_energy(calcData.dim_in_surf, calcData.n_x, calcData.n_y, calcData.n_z, calcData.cutoff, calcData.potentials, calcData.tr.translation, true);
	double E_surf = full_energy(calcData.bulk, calcData.n_x, calcData.n_y, calcData.n_z, calcData.cutoff, calcData.potentials, calcData.tr.translation, true);
	double E_adatom_in_surf = full_energy(calcData.adatom_in_surf, calcData.n_x, calcData.n_y, calcData.n_z, calcData.cutoff, calcData.potentials, calcData.tr.translation, true);
	calcData.E_in_dim = (E_dim_in_surf - E_surf) - 2 * (E_adatom_in_surf - E_surf);
}

void calc_E_on_dim(CalcData &calcData) {
	calcData.tr.make_it_pure();
	double E_dim_on_surf = full_energy(calcData.dim_on_surf, calcData.n_x, calcData.n_y, calcData.n_z + 1, calcData.cutoff, calcData.potentials, calcData.tr.translation, true);
	double E_surf = full_energy(calcData.bulk, calcData.n_x, calcData.n_y, calcData.n_z, calcData.cutoff, calcData.potentials, calcData.tr.translation, true);
	double E_adatom_on_surf = full_energy(calcData.adatom_on_surf, calcData.n_x, calcData.n_y, calcData.n_z + 1, calcData.cutoff, calcData.potentials, calcData.tr.translation, true);
	calcData.E_on_dim = (E_dim_on_surf - E_surf) - 2 * (E_adatom_on_surf - E_surf);
}

double fitting_BB(double x[6], CalcData &calcData) {
	calcData.potentials[0].A_1 = x[0];
	calcData.potentials[0].A_0 = x[1];
	calcData.potentials[0].s = x[2];
	calcData.potentials[0].p = x[3];
	calcData.potentials[0].q = x[4];
	calcData.potentials[0].r0 = x[5];
	calc_E_coh(calcData);
	calc_B(calcData);
	calc_C11_C12(calcData);
	calc_C44(calcData);
	std::cout << calcData.E_coh << " " << calcData.B << " " << calcData.C11 << " " << calcData.C12 << " " << calcData.C44 << '\n';
	return sqrt(((calcData.E_coh_r - calcData.E_coh) * (calcData.E_coh_r - calcData.E_coh) / calcData.E_coh_r * calcData.E_coh_r +
		(calcData.B_r - calcData.B) * (calcData.B_r - calcData.B) / calcData.B_r * calcData.B_r +
		(calcData.C11_r - calcData.C11) * (calcData.C11_r - calcData.C11) / calcData.C11_r * calcData.C11_r +
		(calcData.C12_r - calcData.C12) * (calcData.C12_r - calcData.C12) / calcData.C12_r * calcData.C12_r +
		(calcData.C44_r - calcData.C44) * (calcData.C44_r - calcData.C44) / calcData.C44_r * calcData.C44_r) / 5.0);
}

double fitting_AB(double x[6], CalcData &calcData) {
	calcData.potentials[1].A_1 = x[0];
	calcData.potentials[1].A_0 = x[1];
	calcData.potentials[1].s = x[2];
	calcData.potentials[1].p = x[3];
	calcData.potentials[1].q = x[4];
	calcData.potentials[1].r0 = x[5];
	calc_E_sol(calcData);
	std::cout << calcData.E_sol << "\n";
	return sqrt((calcData.E_sol_r - calcData.E_sol) * (calcData.E_sol_r - calcData.E_sol) / calcData.E_sol_r * calcData.E_sol_r);
}

double fitting_AA(double x[6], CalcData &calcData) {
	calcData.potentials[2].A_1 = x[0];
	calcData.potentials[2].A_0 = x[1];
	calcData.potentials[2].s = x[2];
	calcData.potentials[2].p = x[3];
	calcData.potentials[2].q = x[4];
	calcData.potentials[2].r0 = x[5];
	calc_E_in_dim(calcData);
	calc_E_on_dim(calcData);
	std::cout << calcData.E_in_dim << " " << calcData.E_on_dim << '\n';
	return sqrt(((calcData.E_in_dim_r - calcData.E_in_dim) * (calcData.E_in_dim_r - calcData.E_in_dim) / calcData.E_in_dim_r * calcData.E_in_dim_r +
		(calcData.E_on_dim_r - calcData.E_on_dim) * (calcData.E_on_dim_r - calcData.E_on_dim) / calcData.E_on_dim_r * calcData.E_on_dim_r) / 2.0);
}

void print_potentials(RGL potentials[3], CalcData &calcData) {
	ofstream fout("BB.txt", ofstream::out);

	if (!fout.is_open()) {
		std::cout << "BB.txt can not be opened\n";
		exit(0);
	}

	const int size = 100;
	double x[size];
	double y[size];

	double step = calcData.lattice_constant / 100.0;
	double curr_x = 1.6;

	for (int i = 0; i < size; ++i) {
		x[i] = curr_x;
		curr_x += step;
		double E_r = (potentials[0].A_1 * (x[i] - potentials[0].r0) + potentials[0].A_0)
			* exp(-potentials[0].p * (x[i] / potentials[0].r0 - 1));
		double E_b = - sqrt(potentials[0].s * potentials[0].s
			* exp(-2 * potentials[0].q * (x[i] / potentials[0].r0 - 1)));
		y[i] = E_r + E_b;
	}
	fout << size << '\n';
	fout << "BB potential\n";
	for (int i = 0; i < size; ++i) {
		fout << x[i] << " " << y[i] << '\n';
	}
	fout.close();

	fout.open("AB.txt", ofstream::out);

	if (!fout.is_open()) {
		std::cout << "AB.txt can not be opened\n";
		exit(0);
	}

	curr_x = 1.6;

	for (int i = 0; i < size; ++i) {
		x[i] = curr_x;
		curr_x += step;
		double E_r = (potentials[1].A_1 * (x[i] - potentials[1].r0) + potentials[1].A_0)
			* exp(-potentials[1].p * (x[i] / potentials[1].r0 - 1));
		double E_b = -sqrt(potentials[1].s * potentials[1].s
			* exp(-2 * potentials[1].q * (x[i] / potentials[1].r0 - 1)));
		y[i] = E_r + E_b;
	}
	fout << size << '\n';
	fout << "AB potential\n";
	for (int i = 0; i < size; ++i) {
		fout << x[i] << " " << y[i] << '\n';
	}
	fout.close();

	fout.open("AA.txt", ofstream::out);

	if (!fout.is_open()) {
		std::cout << "AA.txt can not be opened\n";
		exit(0);
	}

	curr_x = 2.25;

	for (int i = 0; i < size; ++i) {
		x[i] = curr_x;
		curr_x += step * 2.0;
		double E_r = (potentials[2].A_1 * (x[i] - potentials[2].r0) + potentials[2].A_0)
			* exp(-potentials[2].p * (x[i] / potentials[2].r0 - 1));
		double E_b = -sqrt(potentials[2].s * potentials[2].s
			* exp(-2 * potentials[2].q * (x[i] / potentials[2].r0 - 1)));
		y[i] = E_r + E_b;
	}
	fout << size << '\n';
	fout << "AA potential\n";
	for (int i = 0; i < size; ++i) {
		fout << x[i] << " " << y[i] << '\n';
	}
	fout.close();
}

/*const double D = 933;//eV
const double a = 4.085 / sqrt(2);*/

void calc(vector <vector<Atom>> &contains, int n_x, int n_y, int n_z, double cutoff, RGL potentials[3], double translation[9], bool surface) {

	double bulk_size_x = translation[0] + translation[3] + translation[6];
	double bulk_size_y = translation[1] + translation[4] + translation[7];
	double bulk_size_z = translation[2] + translation[5] + translation[8];

	//#pragma omp parallel for reduction(+:E_full) num_threads(3) 
	for (int i = 0; i < n_x * n_y * n_z; ++i) {
		for (Atom& particle : contains[i]) {
			int atomcell_x = i % n_x;
			int atomcell_y = i / n_x % n_y;
			int atomcell_z = i / (n_x * n_y);
			double Eb = 0;
			double dEb_dx = 0;
			double dEb_dy = 0;
			double dEb_dz = 0;
			double dEr_dx = 0;
			double dEr_dy = 0;
			double dEr_dz = 0;
			for (int dx = -1; dx <= 1; ++dx) {
				for (int dy = -1; dy <= 1; ++dy) {
					for (int dz = -1; dz <= 1; ++dz) {
						for (Atom particle2 : contains[
							(atomcell_x + n_x + dx) % n_x +
								(atomcell_y + n_y + dy) % n_y * n_x +
								(atomcell_z + n_z + dz) % n_z * n_y * n_x]) {
							//check for translation
							if (abs(particle2.x - particle.x) > bulk_size_x / 2) { //more than one cell far
								if (particle2.x > particle.x) { //move to neighbor cell
									particle2.x -= translation[0];
									particle2.y -= translation[1];
									particle2.z -= translation[2];
								}
								else {
									particle2.x += translation[0];
									particle2.y += translation[1];
									particle2.z += translation[2];
								}
							}
							if (abs(particle2.y - particle.y) > bulk_size_y / 2) { //more than one cell far
								if (particle2.y > particle.y) { //move to neighbor cell
									particle2.x -= translation[3];
									particle2.y -= translation[4];
									particle2.z -= translation[5];
								}
								else {
									particle2.x += translation[3];
									particle2.y += translation[4];
									particle2.z += translation[5];
								}
							}
							if (!surface) {
								if (abs(particle2.z - particle.z) > bulk_size_z / 2) { //more than one cell far
									if (particle2.z > particle.z) { //move to neighbor cell
										particle2.x -= translation[6];
										particle2.y -= translation[7];
										particle2.z -= translation[8];
									}
									else {
										particle2.x += translation[6];
										particle2.y += translation[7];
										particle2.z += translation[8];
									}
								}
							}
							if (toofar(particle, particle2, cutoff)) {
								continue;
							}

							int choice = 0;
							if (particle.impurity && particle2.impurity) {
								choice = 2; //A-A
							}
							else if (!particle.impurity && !particle2.impurity) {
								choice = 0; //B-B
							}
							else {
								choice = 1; //A-B
							}

							//distance
							double r = sqrt((particle2.x - particle.x) * (particle2.x - particle.x)
								+ (particle2.y - particle.y) * (particle2.y - particle.y)
								+ (particle2.z - particle.z) * (particle2.z - particle.z));

							//repulsive and bond energies
							double tmp1 = potentials[choice].s * potentials[choice].s * exp(-2. * potentials[choice].q * (r / potentials[choice].r0 - 1));
							double tmp2 = potentials[choice].q / (potentials[choice].r0 * r);
							dEb_dx += tmp1 * tmp2 * (particle.x - particle2.x);
							dEb_dy += tmp1 * tmp2 * (particle.y - particle2.y);
							dEb_dz += tmp1 * tmp2 * (particle.z - particle2.z);

							Eb += tmp1;

							tmp1 = exp(-potentials[choice].p * (r / potentials[choice].r0 - 1)) * (potentials[choice].A_1
								* (1. / r + (-potentials[choice].r0 + r) * (-potentials[choice].p / (potentials[choice].r0 * r)))
								+ potentials[choice].A_0 * (-potentials[choice].p / (potentials[choice].r0 * r)));

							dEr_dx += tmp1 * (particle.x - particle2.x);

							dEr_dy += tmp1 * (particle.y - particle2.y);

							dEr_dz += tmp1 * (particle.z - particle2.z);
						}
					}
				}
			}
			Eb = sqrt(Eb);
			dEb_dx /= Eb;
			dEb_dy /= Eb;
			dEb_dz /= Eb;
			particle.f_x = -dEr_dx - dEb_dx;
			particle.f_y = -dEr_dy - dEb_dy;
			particle.f_z = -dEr_dz - dEb_dz;
		}
	}
}

void dump(string filename, vector <vector<Atom>> &contains) {
	ofstream out(filename);
	unsigned int size = 0;
	for (unsigned int i = 0; i < contains.size(); ++i) {
		size += contains[i].size();
	}
	out << size << '\n';
	out << "one frame of relaxation\n";
	for (unsigned int i = 0; i < contains.size(); ++i) {
		for (unsigned int j = 0; j < contains[i].size(); ++j) {
			//out << setw(15) << contains[i][j].x << setw(15) << contains[i][j].y << setw(15) << contains[i][j].z << '\n';
			out << contains[i][j].x << ' ' << contains[i][j].y << ' ' << contains[i][j].z << ' ' << contains[i][j].impurity << '\n';
		}
	}
	out.close();
}

int main(int argc, char* argv[]) {

	if (argc != 8) {
		std::cout << "Insufficient number of arguments\n";
		return 0;
	}

	CalcData calcData;

	//cell size
	calcData.n_x = 3;
	calcData.n_y = 3;
	calcData.n_z = 3;

	ifstream fin(argv[7], ifstream::in);

	if (!fin.is_open()) {
		std::cout << "File can not be opened\n";
		return 0;
	}

	string dumb;
	getline(fin, dumb);
	fin >> calcData.E_coh_r;
	fin >> calcData.B_r;
	fin >> calcData.C11_r;
	fin >> calcData.C12_r;
	fin >> calcData.C44_r;
	fin >> calcData.E_sol_r;
	fin >> calcData.E_in_dim_r;
	fin >> calcData.E_on_dim_r;
	fin >> calcData.E_coh_imp;
	fin >> calcData.lattice_constant;

	calcData.cutoff = 1.7 * calcData.lattice_constant;

	calcData.bulk_atom_count = convert(argv[1], calcData.n_x, calcData.n_y, calcData.n_z, calcData.lattice_constant, calcData.bulk, 0);
	convert(argv[2], calcData.n_x, calcData.n_y, calcData.n_z, calcData.lattice_constant, calcData.imp, calcData.bulk_atom_count);
	convert(argv[3], calcData.n_x, calcData.n_y, calcData.n_z, calcData.lattice_constant, calcData.dim_in_surf, calcData.bulk_atom_count);
	convert(argv[4], calcData.n_x, calcData.n_y, calcData.n_z, calcData.lattice_constant, calcData.adatom_in_surf, calcData.bulk_atom_count);
	convert(argv[5], calcData.n_x, calcData.n_y, calcData.n_z, calcData.lattice_constant, calcData.dim_on_surf, calcData.bulk_atom_count);
	convert(argv[6], calcData.n_x, calcData.n_y, calcData.n_z, calcData.lattice_constant, calcData.adatom_on_surf, calcData.bulk_atom_count);

	calcData.tr.lattice_constant = calcData.lattice_constant;
	calcData.tr.alpha = calcData.alpha;
	calcData.V_0 = calcData.lattice_constant * calcData.lattice_constant * calcData.lattice_constant * calcData.n_x * calcData.n_y * calcData.n_z;

	/*auto start_time = chrono::high_resolution_clock::now();

	int n = 6;
	double *start;
	double *xmin;
	double ynewlo;

	start = new double[n];
	xmin = new double[n];
	///////////////////////////////////////////////////////
	std::cout << "\n";
	std::cout << "BB fitting:\n";
	std::cout << "\n";

	start[0] = 0;
	start[1] = 0.1028;
	start[2] = 1.178;
	start[3] = 10.928;
	start[4] = 3.139;
	start[5] = calcData.lattice_constant / sqrt(2);

	std::cout << "\n";
	std::cout << "  Starting point X:\n";
	std::cout << "\n";

	for (int i = 0; i < n; i++)
	{
		std::cout << "  " << setw(14) << start[i] << "\n";
	}

	ynewlo = fitting_BB(start, calcData);

	std::cout << "\n";
	std::cout << "  F(X) = " << ynewlo << "\n";
	std::cout << "\n";

	/*	annealing(fitting_BB, n, start, xmin);
	ynewlo = fitting_BB(xmin);*/

	//*
	/*GRS(fitting_BB, n, start, xmin, calcData);
	ynewlo = fitting_BB(xmin, calcData);
	/*/
	/*NelderMead(fitting_BB, n, start, xmin);
	ynewlo = fitting_BB(xmin);
	//*/

	/*std::cout << "\n";
	std::cout << "  Estimate of minimizing value X*:\n";
	std::cout << "\n";

	for (int i = 0; i < n; i++)
	{
		std::cout << "  " << setw(14) << xmin[i] << "\n";
	}

	std::cout << "\n";
	std::cout << "  F(X*) = " << ynewlo << "\n";
	std::cout << "\n";

	calcData.potentials[0].A_1 = xmin[0];
	calcData.potentials[0].A_0 = xmin[1];
	calcData.potentials[0].s = xmin[2];
	calcData.potentials[0].p = xmin[3];
	calcData.potentials[0].q = xmin[4];
	calcData.potentials[0].r0 = xmin[5];
	///////////////////////////////////////////////////////
	std::cout << "\n";
	std::cout << "AB fitting:\n";
	std::cout << "\n";

	start[0] = 0;
	start[1] = 0.0376;
	start[2] = 1.070;
	start[3] = 16.999;
	start[4] = 1.189;
	start[5] = 3.523 / sqrt(2);

	std::cout << "\n";
	std::cout << "  Starting point X:\n";
	std::cout << "\n";
	for (int i = 0; i < n; i++)
	{
		std::cout << "  " << setw(14) << start[i] << "\n";
	}

	ynewlo = fitting_AB(start, calcData);

	std::cout << "\n";
	std::cout << "  F(X) = " << ynewlo << "\n";
	std::cout << "\n";

	//GRS(fitting_AB, n, start, xmin, calcData);
		
	std::cout << "\n";
	std::cout << "  Estimate of minimizing value X*:\n";
	std::cout << "\n";

	for (int i = 0; i < n; i++)
	{
		std::cout << "  " << setw(14) << xmin[i] << "\n";
	}

	std::cout << "\n";
	std::cout << "  F(X*) = " << ynewlo << "\n";
	std::cout << "\n";

	calcData.potentials[1].A_1 = xmin[0];
	calcData.potentials[1].A_0 = xmin[1];
	calcData.potentials[1].s = xmin[2];
	calcData.potentials[1].p = xmin[3];
	calcData.potentials[1].q = xmin[4];
	calcData.potentials[1].r0 = xmin[5];
	///////////////////////////////////////////////////////
	std::cout << "\n";
	std::cout << "AA fitting:\n";
	std::cout << "\n";

	start[0] = 0;
	start[1] = 0.0376;
	start[2] = 1.070;
	start[3] = 16.999;
	start[4] = 1.189;
	start[5] = 3.523 / sqrt(2);

	std::cout << "\n";
	std::cout << "  Starting point X:\n";
	std::cout << "\n";

	for (int i = 0; i < n; i++)
	{
		std::cout << "  " << setw(14) << start[i] << "\n";
	}

	ynewlo = fitting_AA(start, calcData);

	std::cout << "\n";
	std::cout << "  F(X) = " << ynewlo << "\n";

	//GRS(fitting_AA, n, start, xmin, calcData);

	std::cout << "\n";
	std::cout << "  Estimate of minimizing value X*:\n";
	std::cout << "\n";
	for (int i = 0; i < n; i++)
	{
		std::cout << "  " << setw(14) << xmin[i] << "\n";
	}

	std::cout << "\n";
	std::cout << "  F(X*) = " << ynewlo << "\n";

	std::cout << "\n";

	calcData.potentials[2].A_1 = xmin[0];
	calcData.potentials[2].A_0 = xmin[1];
	calcData.potentials[2].s = xmin[2];
	calcData.potentials[2].p = xmin[3];
	calcData.potentials[2].q = xmin[4];
	calcData.potentials[2].r0 = xmin[5];

	delete[] start;
	delete[] xmin;
	
	std::cout << "\nB-B potential:\n" <<
		"A1 = " << calcData.potentials[0].A_1 << '\n' <<
		"A0 = " << calcData.potentials[0].A_0 << '\n' <<
		"s = " << calcData.potentials[0].s << '\n' <<
		"p = " << calcData.potentials[0].p << '\n' <<
		"q = " << calcData.potentials[0].q << '\n' <<
		"r0 = " << calcData.potentials[0].r0 << '\n';

	std::cout << "\nA-B potential:\n" <<
		"A1 = " << calcData.potentials[1].A_1 << '\n' <<
		"A0 = " << calcData.potentials[1].A_0 << '\n' <<
		"s = " << calcData.potentials[1].s << '\n' <<
		"p = " << calcData.potentials[1].p << '\n' <<
		"q = " << calcData.potentials[1].q << '\n' <<
		"r0 = " << calcData.potentials[1].r0 << '\n';

	std::cout << "\nA-A potential:\n" <<
		"A1 = " << calcData.potentials[2].A_1 << '\n' <<
		"A0 = " << calcData.potentials[2].A_0 << '\n' <<
		"s = " << calcData.potentials[2].s << '\n' <<
		"p = " << calcData.potentials[2].p << '\n' <<
		"q = " << calcData.potentials[2].q << '\n' <<
		"r0 = " << calcData.potentials[2].r0 << '\n';

	print_potentials(calcData.potentials, calcData);

	auto end_time = chrono::high_resolution_clock::now();
	

	std::cout << '\n' << chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() << " ms\n";
	*/
	//B-B potential:
	calcData.potentials[0].A_1 = 0;
	calcData.potentials[0].A_0 = 0.1028;
	calcData.potentials[0].s = 1.178;
	calcData.potentials[0].p = 10.928;
	calcData.potentials[0].q = 3.14603;
	calcData.potentials[0].r0 = 2.88853;

	//A-B potential:
	calcData.potentials[1].A_1 = 0;
	calcData.potentials[1].A_0 = 0.042866;
	calcData.potentials[1].s = 0.990775;
	calcData.potentials[1].p = 16.3642;
	calcData.potentials[1].q = 1.189;
	calcData.potentials[1].r0 = 2.45985;

	//A-A potential:
	calcData.potentials[2].A_1 = 0;
	calcData.potentials[2].A_0 = 0.0377356;
	calcData.potentials[2].s = 1.07419;
	calcData.potentials[2].p = 17.214;
	calcData.potentials[2].q = 1.42031;
	calcData.potentials[2].r0 = 3.04121;

	calcData.tr.make_it_pure();
	//calc(calcData.bulk, calcData.n_x, calcData.n_y, calcData.n_z, calcData.cutoff, calcData.potentials, calcData.tr.translation, false);

	vector <vector<Atom>> prev = calcData.dim_on_surf;
	vector <vector<Atom>> curr = calcData.dim_on_surf;
	//calculate prev forces
	calc(prev, calcData.n_x, calcData.n_y, calcData.n_z + 1, calcData.cutoff, calcData.potentials, calcData.tr.translation, true);
	//dump iteration
	cout << "frame_0.txt\n";
	dump("frame_0.txt", curr);
	
	int iterations_count;
	cout << "iterations count: ";
	cin >> iterations_count;
	double dt;
	cout << "dt: ";
	cin >> dt;
	for (int it_cnt = 0; it_cnt < iterations_count; ++it_cnt) {
		//count positions
		for (unsigned int i = 0; i < prev.size(); ++i) {
			for (unsigned int j = 0; j < prev[i].size(); ++j) {
				curr[i][j].x = prev[i][j].x + prev[i][j].v_x * dt + prev[i][j].f_x * dt * dt / (2 * curr[i][j].mass);
				curr[i][j].y = prev[i][j].y + prev[i][j].v_y * dt + prev[i][j].f_y * dt * dt / (2 * curr[i][j].mass);
				curr[i][j].z = prev[i][j].z + prev[i][j].v_z * dt + prev[i][j].f_z * dt * dt / (2 * curr[i][j].mass);
			}
		}
		//calculate
		calc(curr, calcData.n_x, calcData.n_y, calcData.n_z + 1, calcData.cutoff, calcData.potentials, calcData.tr.translation, true);
		//count velocities
		for (unsigned int i = 0; i < curr.size(); ++i) {
			for (unsigned int j = 0; j < curr[i].size(); ++j) {
				curr[i][j].v_x = prev[i][j].v_x + (prev[i][j].f_x + curr[i][j].f_x) / (2 * curr[i][j].mass) * dt;
				curr[i][j].v_y = prev[i][j].v_y + (prev[i][j].f_y + curr[i][j].f_y) / (2 * curr[i][j].mass) * dt;
				curr[i][j].v_z = prev[i][j].v_z + (prev[i][j].f_z + curr[i][j].f_z) / (2 * curr[i][j].mass) * dt;
			}
		}
		//dump iteration
		cout << "frame_" << it_cnt + 1 << ".txt\n";
		dump("frame_" + to_string(it_cnt + 1) + ".txt", curr);
		//swap
		prev.swap(curr);
	}
	
	system("pause");
	return 0;
}
