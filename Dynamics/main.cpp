#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <chrono>
#include "asa047.hpp"
using namespace std;

class Atom {
public:
	double x = 0;
	double y = 0;
	double z = 0;
	bool impurity = false;
};

class RGL {
public:
	double A_1 = 0;
	double A_0 = 0;
	double s = 0;
	double p = 0;
	double q = 0;
	double r0 = 0;
};

class Translation {
public:
	double translation[9];
	double lattice_constant = 0;
	double alpha = 0;
	void make_it_pure() {
		//x
		translation[0] = 3 * lattice_constant;
		translation[1] = 0;
		translation[2] = 0;
		//y
		translation[3] = 0;
		translation[4] = 3 * lattice_constant;
		translation[5] = 0;
		//z
		translation[6] = 0;
		translation[7] = 0;
		translation[8] = 3 * lattice_constant;
	}
	void make_it_B() {
		//x
		translation[0] = 3 * lattice_constant * (1 + alpha);
		translation[1] = 0;
		translation[2] = 0;
		//y
		translation[3] = 0;
		translation[4] = 3 * lattice_constant * (1 + alpha);
		translation[5] = 0;
		//z
		translation[6] = 0;
		translation[7] = 0;
		translation[8] = 3 * lattice_constant * (1 + alpha);
	}
	void make_it_C11() {
		//x
		translation[0] = 3 * lattice_constant * (1 + alpha);
		translation[1] = 0;
		translation[2] = 0;
		//y
		translation[3] = 0;
		translation[4] = 3 * lattice_constant * (1 + alpha);
		translation[5] = 0;
		//z
		translation[6] = 0;
		translation[7] = 0;
		translation[8] = 3 * lattice_constant;
	}
	void make_it_C12() {
		//x
		translation[0] = 3 * lattice_constant * (1 + alpha);
		translation[1] = 0;
		translation[2] = 0;
		//y
		translation[3] = 0;
		translation[4] = 3 * lattice_constant * (1 - alpha);
		translation[5] = 0;
		//z
		translation[6] = 0;
		translation[7] = 0;
		translation[8] = 3 * lattice_constant;
	}
	void make_it_C44() {
		//x
		translation[0] = 3 * lattice_constant;
		translation[1] = alpha * 3 * lattice_constant;
		translation[2] = 0;
		//y
		translation[3] = alpha * 3 * lattice_constant;
		translation[4] = 3 * lattice_constant;
		translation[5] = 0;
		//z
		translation[6] = 0;
		translation[7] = 0;
		translation[8] = 3 * lattice_constant / (1 - alpha * alpha);
	}
};

//without impurities
vector <vector<Atom>> bulk;
//with one impurity for E_sol
vector <vector<Atom>> imp;
//with dimer in surf
vector <vector<Atom>> dim;
//with adatom in surf
vector <vector<Atom>> adatom;
//with dimer on surf
vector <vector<Atom>> dim_on_surf;
//with adatom on surf
vector <vector<Atom>> adatom_on_surf;
//surf energy
vector <vector<Atom>> surf;

Translation tr;
int n_x, n_y, n_z;
double alpha = 0.01;
double lattice_constant;
double cutoff;
RGL potentials[3];
double V_0;
double E_full;
double qsi = 0.8018993929636421;
int atom_count;
double E_coh;
double B;
double C11;
double C12;
double C44;
double E_sol;
double E_in_dim;
double E_on_dim;

//reference values
double E_coh_r;
double B_r;
double C11_r;
double C12_r;
double C44_r;
double E_sol_r;
double E_in_dim_r;
double E_on_dim_r;

double E_coh_imp;

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

	#pragma omp parallel for num_threads(3) 
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
								choice = 2; //B-B
							} else if (!particle.impurity && !particle2.impurity) {
								choice = 0; //A-A
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
			#pragma omp critical
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

int convert(string filename, int n_x, int n_y, int n_z, double lattice_constant, vector <vector<Atom>> &out) {

	out.resize(n_x * n_y * n_z);
	vector<Atom> init;

	ifstream fin(filename, ifstream::in);

	if (!fin.is_open()) {
		cout << "File can not be opened\n";
		exit(0);
	}

	int size = 0;
	fin >> size;

	while (fin.good()) {
		Atom tmp;
		fin >> tmp.impurity;
		fin >> tmp.x;
		fin >> tmp.y;
		fin >> tmp.z;
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

void calc_E_coh() {
	tr.make_it_pure();
	E_full = full_energy(bulk, n_x, n_y, n_z, cutoff, potentials, tr.translation, false);
	E_coh = E_full / atom_count;
}

void calc_B() {
	vector <vector<Atom>> curr = bulk;
	stretch(curr, n_x, n_y, n_z, alpha);
	tr.alpha = alpha;
	tr.make_it_B();
	double E_full_stretched = full_energy(curr, n_x, n_y, n_z, cutoff, potentials, tr.translation, false);
	curr = bulk;

	stretch(curr, n_x, n_y, n_z, -alpha);
	tr.alpha = -alpha;
	tr.make_it_B();
	double E_full_compressed = full_energy(curr, n_x, n_y, n_z, cutoff, potentials, tr.translation, false);
	curr = bulk;
	B = (2 * (E_full_stretched - 2 * E_full + E_full_compressed) * qsi) / (9 * V_0 * alpha * alpha);
}

void calc_C11_C12() {
	vector <vector<Atom>> curr = bulk;
	elastic_C11(curr, n_x, n_y, n_z, alpha);
	tr.alpha = alpha;
	tr.make_it_C11();
	double E_full_stretched_C11 = full_energy(curr, n_x, n_y, n_z, cutoff, potentials, tr.translation, false);
	curr = bulk;

	elastic_C11(curr, n_x, n_y, n_z, -alpha);
	tr.alpha = -alpha;
	tr.make_it_C11();
	double E_full_compressed_C11 = full_energy(curr, n_x, n_y, n_z, cutoff, potentials, tr.translation, false);
	curr = bulk;

	elastic_C12(curr, n_x, n_y, n_z, alpha);
	tr.alpha = alpha;
	tr.make_it_C12();
	double E_full_stretched_C12 = full_energy(curr, n_x, n_y, n_z, cutoff, potentials, tr.translation, false);
	curr = bulk;

	elastic_C12(curr, n_x, n_y, n_z, -alpha);
	tr.alpha = -alpha;
	tr.make_it_C12();
	double E_full_compressed_C12 = full_energy(curr, n_x, n_y, n_z, cutoff, potentials, tr.translation, false);
	curr = bulk;

	C11 = ((E_full_stretched_C11 - 2 * E_full + E_full_compressed_C11) / (alpha * alpha) + (E_full_stretched_C12 - 2 * E_full + E_full_compressed_C12) / (alpha * alpha)) * qsi / (2 * V_0);

	C12 = ((E_full_stretched_C11 - 2 * E_full + E_full_compressed_C11) / (alpha * alpha) - (E_full_stretched_C12 - 2 * E_full + E_full_compressed_C12) / (alpha * alpha)) * qsi / (2 * V_0);
}

void calc_C44() {
	vector <vector<Atom>> curr = bulk;
	elastic_C44(curr, n_x, n_y, n_z, alpha);
	tr.alpha = alpha;
	tr.make_it_C44();
	double E_full_stretched_C44 = full_energy(curr, n_x, n_y, n_z, cutoff, potentials, tr.translation, false);
	curr = bulk;

	elastic_C44(curr, n_x, n_y, n_z, -alpha);
	tr.alpha = -alpha;
	tr.make_it_C44();
	double E_full_compressed_C44 = full_energy(curr, n_x, n_y, n_z, cutoff, potentials, tr.translation, false);
	curr = bulk;

	C44 = ((E_full_stretched_C44 - 2 * E_full + E_full_compressed_C44)  * qsi) / (2 * V_0 * alpha * alpha);
}

void calc_E_sol() {
	tr.make_it_pure();
	double E_AB = full_energy(imp, n_x, n_y, n_z, cutoff, potentials, tr.translation, false);
	double E_B = full_energy(bulk, n_x, n_y, n_z, cutoff, potentials, tr.translation, false);
	E_sol = E_AB - E_B - E_coh_imp + E_coh;
}

void calc_E_in_dim() {
	double E_dim_surf = full_energy(dim, n_x, n_y, n_z, cutoff, potentials, tr.translation, true);
	double E_surf = full_energy(bulk, n_x, n_y, n_z, cutoff, potentials, tr.translation, true);
	double E_adatom_surf = full_energy(adatom, n_x, n_y, n_z, cutoff, potentials, tr.translation, true);
	E_in_dim = (E_dim_surf - E_surf) - 2 * (E_adatom_surf - E_surf);
}

void calc_E_on_dim() {
	double E_dim_surf2 = full_energy(dim_on_surf, n_x, n_y, n_z, cutoff, potentials, tr.translation, true);
	double E_surf2 = full_energy(surf, n_x, n_y, n_z, cutoff, potentials, tr.translation, true);
	double E_adatom_surf2 = full_energy(adatom_on_surf, n_x, n_y, n_z, cutoff, potentials, tr.translation, true);
	E_on_dim = (E_dim_surf2 - E_surf2) - 2 * (E_adatom_surf2 - E_surf2);
}

double fitting_BB(double x[6]) {
	potentials[0].A_1 = x[0];
	potentials[0].A_0 = x[1];
	potentials[0].s = x[2];
	potentials[0].p = x[3];
	potentials[0].q = x[4];
	potentials[0].r0 = x[5];
	calc_E_coh();
	calc_B();
	calc_C11_C12();
	calc_C44();
	cout << E_coh << " " << B << " " << C11 << " " << C12 << " " << C44 << '\n';
	return sqrt(((E_coh_r - E_coh) * (E_coh_r - E_coh) / E_coh_r * E_coh_r +
		(B_r - B) * (B_r - B) / B_r * B_r +
		(C11_r - C11) * (C11_r - C11) / C11_r * C11_r +
		(C12_r - C12) * (C12_r - C12) / C12_r * C12_r +
		(C44_r - C44) * (C44_r - C44) / C44_r * C44_r) / 5.0);
}

double fitting_AB(double x[6]) {
	potentials[1].A_1 = x[0];
	potentials[1].A_0 = x[1];
	potentials[1].s = x[2];
	potentials[1].p = x[3];
	potentials[1].q = x[4];
	potentials[1].r0 = x[5];
	calc_E_sol();
	cout << E_sol << "\n";
	return sqrt((E_sol_r - E_sol) * (E_sol_r - E_sol) / E_sol_r * E_sol_r);
}

double fitting_AA(double x[6]) {
		potentials[2].A_1 = x[0];
		potentials[2].A_0 = x[1];
		potentials[2].s = x[2];
		potentials[2].p = x[3];
		potentials[2].q = x[4];
		potentials[2].r0 = x[5];
		calc_E_in_dim();
		calc_E_on_dim();
		cout << E_in_dim << " " << E_on_dim << '\n';
		return sqrt(((E_in_dim_r - E_in_dim) * (E_in_dim_r - E_in_dim) / E_in_dim_r * E_in_dim_r +
			(E_on_dim_r - E_on_dim) * (E_on_dim_r - E_on_dim) / E_on_dim_r * E_on_dim_r) / 2.0);
	}

void print_potentials(RGL potentials[3]) {
	ofstream fout("BB.txt", ofstream::out);

	if (!fout.is_open()) {
		cout << "BB.txt can not be opened\n";
		exit(0);
	}

	const int size = 100;
	double x[size];
	double y[size];

	double step = lattice_constant / 50.0;
	double curr_x = step;

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
		cout << "AB.txt can not be opened\n";
		exit(0);
	}

	curr_x = step;

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
		cout << "AA.txt can not be opened\n";
		exit(0);
	}

	curr_x = step;

	for (int i = 0; i < size; ++i) {
		x[i] = curr_x;
		curr_x += step;
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

int main(int argc, char* argv[]) {

	if (argc != 8) {
		cout << "Insufficient number of arguments\n";
		return 0;
	}

	//cell size
	n_x = 3;
	n_y = 3;
	n_z = 3;

	lattice_constant = 4.085;
	cutoff = 1.7 * lattice_constant;


	atom_count = convert(argv[1], n_x, n_y, n_z, lattice_constant, bulk);
	convert(argv[2], n_x, n_y, n_z, lattice_constant, imp);
	convert(argv[3], n_x, n_y, n_z, lattice_constant, dim);
	convert(argv[4], n_x, n_y, n_z, lattice_constant, adatom);
	convert(argv[5], n_x, n_y, n_z, lattice_constant, dim_on_surf);
	convert(argv[6], n_x, n_y, n_z, lattice_constant, adatom_on_surf);
	convert(argv[7], n_x, n_y, n_z, lattice_constant, surf);

	//Ag-Ni
	E_coh_r = -2.960;
	B_r = 1.08;
	C11_r = 1.32;
	C12_r = 0.97;
	C44_r = 0.51;
	E_sol_r = 0.539;
	E_in_dim_r = 0.06;
	E_on_dim_r = -0.56;
	E_coh_imp = -4.435;

	tr.lattice_constant = lattice_constant;
	tr.alpha = alpha;
	V_0 = lattice_constant * lattice_constant * lattice_constant * n_x * n_y * n_z;

	auto start_time = chrono::high_resolution_clock::now();

	int icount;
	int ifault;
	int kcount;
	int konvge;
	int n;
	int numres;
	double reqmin;
	double *start;
	double *step;
	double *xmin;
	double ynewlo;

	n = 6;

	start = new double[n];
	step = new double[n];
	xmin = new double[n];
	///////////////////////////////////////////////////////
	cout << "\n";
	cout << "BB fitting:\n";
	cout << "\n";

	start[0] = 0;
	start[1] = 0.1028;
	start[2] = 1.178;
	start[3] = 10.928;
	start[4] = 3.139;
	start[5] = lattice_constant / sqrt(2);

	reqmin = 1.0E-08;

	step[0] = 1.0;
	step[1] = 1.0;
	step[2] = 1.0;
	step[3] = 1.0;
	step[4] = 1.0;
	step[5] = 1.0;

	konvge = 10;
	kcount = 500;

	cout << "\n";
	cout << "  Starting point X:\n";
	cout << "\n";
	for (int i = 0; i < n; i++)
	{
		cout << "  " << setw(14) << start[i] << "\n";
	}

	ynewlo = fitting_BB(start);

	cout << "\n";
	cout << "  F(X) = " << ynewlo << "\n";

	nelmin(fitting_BB, n, start, xmin, &ynewlo, reqmin, step,
		konvge, kcount, &icount, &numres, &ifault);

	cout << "\n";
	cout << "  Return code IFAULT = " << ifault << "\n";
	cout << "\n";
	cout << "  Estimate of minimizing value X*:\n";
	cout << "\n";
	for (int i = 0; i < n; i++)
	{
		cout << "  " << setw(14) << xmin[i] << "\n";
	}

	cout << "\n";
	cout << "  F(X*) = " << ynewlo << "\n";

	cout << "\n";
	cout << "  Number of iterations = " << icount << "\n";
	cout << "  Number of restarts =   " << numres << "\n";
	
	if (ifault != 0) {
		return 0;
	}

	potentials[0].A_1 = xmin[0];
	potentials[0].A_0 = xmin[1];
	potentials[0].s = xmin[2];
	potentials[0].p = xmin[3];
	potentials[0].q = xmin[4];
	potentials[0].r0 = xmin[5];
	///////////////////////////////////////////////////////
	cout << "\n";
	cout << "AB fitting:\n";
	cout << "\n";

	start[0] = 0;
	start[1] = 0.1028;
	start[2] = 1.178;
	start[3] = 10.928;
	start[4] = 3.139;
	start[5] = lattice_constant / sqrt(2);

	reqmin = 1.0E-08;

	step[0] = 1.0;
	step[1] = 1.0;
	step[2] = 1.0;
	step[3] = 1.0;
	step[4] = 1.0;
	step[5] = 1.0;

	konvge = 10;
	kcount = 500;

	cout << "\n";
	cout << "  Starting point X:\n";
	cout << "\n";
	for (int i = 0; i < n; i++)
	{
		cout << "  " << setw(14) << start[i] << "\n";
	}

	ynewlo = fitting_AB(start);

	cout << "\n";
	cout << "  F(X) = " << ynewlo << "\n";

	nelmin(fitting_AB, n, start, xmin, &ynewlo, reqmin, step,
		konvge, kcount, &icount, &numres, &ifault);

	cout << "\n";
	cout << "  Return code IFAULT = " << ifault << "\n";
	cout << "\n";
	cout << "  Estimate of minimizing value X*:\n";
	cout << "\n";
	for (int i = 0; i < n; i++)
	{
		cout << "  " << setw(14) << xmin[i] << "\n";
	}

	cout << "\n";
	cout << "  F(X*) = " << ynewlo << "\n";

	cout << "\n";
	cout << "  Number of iterations = " << icount << "\n";
	cout << "  Number of restarts =   " << numres << "\n";

	if (ifault != 0) {
		return 0;
	}

	potentials[1].A_1 = xmin[0];
	potentials[1].A_0 = xmin[1];
	potentials[1].s = xmin[2];
	potentials[1].p = xmin[3];
	potentials[1].q = xmin[4];
	potentials[1].r0 = xmin[5];
	///////////////////////////////////////////////////////
	cout << "\n";
	cout << "AA fitting:\n";
	cout << "\n";

	start[0] = 0;
	start[1] = 0.1028;
	start[2] = 1.178;
	start[3] = 10.928;
	start[4] = 3.139;
	start[5] = lattice_constant / sqrt(2);

	reqmin = 1.0E-08;

	step[0] = 1.0;
	step[1] = 1.0;
	step[2] = 1.0;
	step[3] = 1.0;
	step[4] = 1.0;
	step[5] = 1.0;

	konvge = 10;
	kcount = 500;

	cout << "\n";
	cout << "  Starting point X:\n";
	cout << "\n";
	for (int i = 0; i < n; i++)
	{
		cout << "  " << setw(14) << start[i] << "\n";
	}

	ynewlo = fitting_AA(start);

	cout << "\n";
	cout << "  F(X) = " << ynewlo << "\n";

	nelmin(fitting_AA, n, start, xmin, &ynewlo, reqmin, step,
		konvge, kcount, &icount, &numres, &ifault);

	cout << "\n";
	cout << "  Return code IFAULT = " << ifault << "\n";
	cout << "\n";
	cout << "  Estimate of minimizing value X*:\n";
	cout << "\n";
	for (int i = 0; i < n; i++)
	{
		cout << "  " << setw(14) << xmin[i] << "\n";
	}

	cout << "\n";
	cout << "  F(X*) = " << ynewlo << "\n";

	cout << "\n";
	cout << "  Number of iterations = " << icount << "\n";
	cout << "  Number of restarts =   " << numres << "\n";

	if (ifault != 0) {
		return 0;
	}

	potentials[2].A_1 = xmin[0];
	potentials[2].A_0 = xmin[1];
	potentials[2].s = xmin[2];
	potentials[2].p = xmin[3];
	potentials[2].q = xmin[4];
	potentials[2].r0 = xmin[5];

	delete[] start;
	delete[] step;
	delete[] xmin;

	cout << "\nB-B potential:\n" <<
		"A1 = " << potentials[0].A_1 << '\n' <<
		"A0 = " << potentials[0].A_0 << '\n' <<
		"s = " << potentials[0].s << '\n' <<
		"p = " << potentials[0].p << '\n' <<
		"q = " << potentials[0].q << '\n' <<
		"r0 = " << potentials[0].r0 << '\n';

	cout << "\nA-B potential:\n" <<
		"A1 = " << potentials[1].A_1 << '\n' <<
		"A0 = " << potentials[1].A_0 << '\n' <<
		"s = " << potentials[1].s << '\n' <<
		"p = " << potentials[1].p << '\n' <<
		"q = " << potentials[1].q << '\n' <<
		"r0 = " << potentials[1].r0 << '\n';

	cout << "\nA-A potential:\n" <<
		"A1 = " << potentials[2].A_1 << '\n' <<
		"A0 = " << potentials[2].A_0 << '\n' <<
		"s = " << potentials[2].s << '\n' <<
		"p = " << potentials[2].p << '\n' <<
		"q = " << potentials[2].q << '\n' <<
		"r0 = " << potentials[2].r0 << '\n';

	print_potentials(potentials);

	auto end_time = chrono::high_resolution_clock::now();
	cout << '\n' << chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() << " ms\n";
	return 0;
}
