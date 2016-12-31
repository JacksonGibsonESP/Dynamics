#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

const double lattice_constant = 3.615;
const double cutoff = 1.7 * lattice_constant;

/*
//Cu
const double A_1 = 0;
const double A_0 = 0.0855;
const double s = 1.2240;
const double p = 10.96;
const double q = 2.278;
const double r0 = lattice_constant / sqrt(2);
*/
/*
//Ag
const double lattice_constant = 4.085;
const double cutoff = 1.7 * lattice_constant;
const double A_1 = 0;
const double A_0 = 0.1028;
const double s = 1.178;
const double p = 10.928;
const double q = 3.139;
const double r0 = lattice_constant / sqrt(2);*/

const double qsi = 0.8018993929636421;

class Atom {
public:
	double x = 0;
	double y = 0;
	double z = 0;
	double E_r = 0;
	double E_b = 0;
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

bool toofar(Atom part, Atom part2, double cutoff) {
	if (part.x == part2.x && part.y == part2.y && part.z == part2.z)
		return true;
	return sqrt((part2.x - part.x) * (part2.x - part.x)
		+ (part2.y - part.y) * (part2.y - part.y)
		+ (part2.z - part.z) * (part2.z - part.z)) > cutoff ? true : false;
}

double full_energy(vector <vector<Atom>> &contains, int n_x, int n_y, int n_z, RGL potentials[3], bool impurity, double translation[9]) {

	double E_full = 0;
	double bulk_size_x = translation[0] + translation[3] + translation[6];
	double bulk_size_y = translation[1] + translation[4] + translation[7];
	double bulk_size_z = translation[2] + translation[5] + translation[8];

	for (int i = 0; i < n_x * n_y * n_z; ++i) {
		for (Atom& particle : contains[i]) {
			int atomcell_x = i % n_x;
			int atomcell_y = i / n_x % n_y;
			int atomcell_z = i / (n_x * n_y);
			double E_b_temp = 0;
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
							if (toofar(particle, particle2, cutoff)) {
								continue;
							}

							int choice = 0;
							if (!impurity) { //turn on switcher
								choice = 0;
							} else if (particle.impurity && particle2.impurity) {
								choice = 2; //B-B
							} else if (!particle.impurity && !particle2.impurity) {
								choice = 0; //A-A
							} else if (!particle.impurity || !particle2.impurity) {
								choice = 1; //A-B
							}

							//distance
							double r = sqrt((particle2.x - particle.x) * (particle2.x - particle.x)
								+ (particle2.y - particle.y) * (particle2.y - particle.y)
								+ (particle2.z - particle.z) * (particle2.z - particle.z));

							//repulsive and bond energies
							particle.E_r += (potentials[choice].A_1 * (r - potentials[choice].r0) + potentials[choice].A_0)
								* exp(-potentials[choice].p * (r / potentials[choice].r0 - 1));
							E_b_temp += potentials[choice].s * potentials[choice].s
								* exp(-2 * potentials[choice].q * (r / potentials[choice].r0 - 1));
						}
					}
				}
			}
			particle.E_b -= sqrt(E_b_temp);
			E_full += particle.E_r + particle.E_b;
			particle.E_r = 0;
			particle.E_b = 0;
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

int main(int argc, char* argv[]) {

	vector<Atom> init;

	ifstream fin(argv[1], ifstream::in);

	if (!fin.is_open()) {
		cout << "File can not be opened\n";
		return 0;
	}

	string dumb;
	int size = 0;
	fin >> size;
	
	double bulk_size_x = 3 * lattice_constant;
	double bulk_size_y = 3 * lattice_constant;
	double bulk_size_z = 3 * lattice_constant;

	while (fin.good()) {
		Atom tmp;
		fin >> tmp.impurity;
		fin >> tmp.x;
		fin >> tmp.y;
		fin >> tmp.z;
		init.push_back(tmp);
	}
	fin.close();

	int n_x = 3;
	int n_y = 3;
	int n_z = 3;
	double step = lattice_constant;

	vector <vector<Atom>> prev(n_x * n_y * n_z);
	vector <vector<Atom>> curr(n_x * n_y * n_z);
	for (unsigned int i = 0; i < init.size(); ++i) {
		Atom tmp = init[i];
		int x = int(tmp.x / step);
		int y = int(tmp.y / step);
		int z = int(tmp.z / step);
		prev[z * n_y * n_x + y * n_x + x].push_back(tmp);
		curr[z * n_y * n_x + y * n_x + x].push_back(tmp);
	}

	RGL potentials[3];

	//A-A Cu
	potentials[0].A_1 = 0;
	potentials[0].A_0 = 0.0854;
	potentials[0].s = 1.2243;
	potentials[0].p = 10.939;
	potentials[0].q = 2.2799;
	potentials[0].r0 = 2.5563;

	//impurity
	//A-B Co-Cu
	potentials[1].A_1 = -0.7922;
	potentials[1].A_0 = -0.0487;
	potentials[1].s = 0.7356;
	potentials[1].p = 8.1825;
	potentials[1].q = 3.344;
	potentials[1].r0 = 2.4049;
	
	//B-B Co
	potentials[2].A_1 = -0.3583;
	potentials[2].A_0 = 0.1385;
	potentials[2].s = 1.5247;
	potentials[2].p = 7.6788;
	potentials[2].q = 2.139;
	potentials[2].r0 = 2.378;

	double translation[9];
	//x
	translation[0] = bulk_size_x;
	translation[1] = 0;
	translation[2] = 0;
	//y
	translation[3] = 0;
	translation[4] = bulk_size_y;
	translation[5] = 0;
	//z
	translation[6] = 0;
	translation[7] = 0;
	translation[8] = bulk_size_z;

	double E_full = full_energy(prev, n_x, n_y, n_z, potentials, true, translation);
	double E_full_without_impurity = full_energy(prev, n_x, n_y, n_z, potentials, false, translation); //for E_sol
	cout << "E_full = " << E_full << endl;
	double E_coh = E_full / size;
	cout << "E_coh = " << E_coh << endl;
	cout << "alpha = ";
	double alpha;
	cin >> alpha;

	double V_0 = lattice_constant * lattice_constant * lattice_constant * n_x * n_y * n_z;

	//B////////////////////////////////////////////////////////////////////////////////////
	stretch(curr, n_x, n_y, n_z, alpha);
	//x
	translation[0] = bulk_size_x * (1 + alpha);
	translation[1] = 0;
	translation[2] = 0;
	//y
	translation[3] = 0;
	translation[4] = bulk_size_y * (1 + alpha);
	translation[5] = 0;
	//z
	translation[6] = 0;
	translation[7] = 0;
	translation[8] = bulk_size_z * (1 + alpha);
	double E_full_stretched = full_energy(curr, n_x, n_y, n_z, potentials, true, translation);
	curr = prev;

	stretch(curr, n_x, n_y, n_z, -alpha);
	//x
	translation[0] = bulk_size_x * (1 - alpha);
	translation[1] = 0;
	translation[2] = 0;
	//y
	translation[3] = 0;
	translation[4] = bulk_size_y * (1 - alpha);
	translation[5] = 0;
	//z
	translation[6] = 0;
	translation[7] = 0;
	translation[8] = bulk_size_z * (1 - alpha);
	double E_full_compressed = full_energy(curr, n_x, n_y, n_z, potentials, true, translation);
	curr = prev;

	double B = (2 * (E_full_stretched - 2 * E_full + E_full_compressed) * qsi) / (9 * V_0 * alpha * alpha);
	cout << "B = " << B << endl;

	//C11 C12//////////////////////////////////////////////////////////////////////////////
	elastic_C11(curr, n_x, n_y, n_z, alpha);
	//x
	translation[0] = bulk_size_x * (1 + alpha);
	translation[1] = 0;
	translation[2] = 0;
	//y
	translation[3] = 0;
	translation[4] = bulk_size_y * (1 + alpha);
	translation[5] = 0;
	//z
	translation[6] = 0;
	translation[7] = 0;
	translation[8] = bulk_size_z;
	double E_full_stretched_C11 = full_energy(curr, n_x, n_y, n_z, potentials, true, translation);
	curr = prev;

	elastic_C11(curr, n_x, n_y, n_z, -alpha);
	//x
	translation[0] = bulk_size_x * (1 - alpha);
	translation[1] = 0;
	translation[2] = 0;
	//y
	translation[3] = 0;
	translation[4] = bulk_size_y * (1 - alpha);
	translation[5] = 0;
	//z
	translation[6] = 0;
	translation[7] = 0;
	translation[8] = bulk_size_z;
	double E_full_compressed_C11 = full_energy(curr, n_x, n_y, n_z, potentials, true, translation);
	curr = prev;

	elastic_C12(curr, n_x, n_y, n_z, alpha);
	//x
	translation[0] = bulk_size_x * (1 + alpha);
	translation[1] = 0;
	translation[2] = 0;
	//y
	translation[3] = 0;
	translation[4] = bulk_size_y * (1 - alpha);
	translation[5] = 0;
	//z
	translation[6] = 0;
	translation[7] = 0;
	translation[8] = bulk_size_z;
	double E_full_stretched_C12 = full_energy(curr, n_x, n_y, n_z, potentials, true, translation);
	curr = prev;

	elastic_C12(curr, n_x, n_y, n_z, -alpha);
	//x
	translation[0] = bulk_size_x * (1 - alpha);
	translation[1] = 0;
	translation[2] = 0;
	//y
	translation[3] = 0;
	translation[4] = bulk_size_y * (1 + alpha);
	translation[5] = 0;
	//z
	translation[6] = 0;
	translation[7] = 0;
	translation[8] = bulk_size_z;
	double E_full_compressed_C12 = full_energy(curr, n_x, n_y, n_z, potentials, true, translation);
	curr = prev;

	double C11 = ((E_full_stretched_C11 - 2 * E_full + E_full_compressed_C11) / (alpha * alpha) + (E_full_stretched_C12 - 2 * E_full + E_full_compressed_C12) / (alpha * alpha)) * qsi / (2 * V_0);
	cout << "C11 = " << C11 << endl;

	double C12 = ((E_full_stretched_C11 - 2 * E_full + E_full_compressed_C11) / (alpha * alpha) - (E_full_stretched_C12 - 2 * E_full + E_full_compressed_C12) / (alpha * alpha)) * qsi / (2 * V_0);
	cout << "C12 = " << C12 << endl;
	
	//C44//////////////////////////////////////////////////////////////////////////////////
	elastic_C44(curr, n_x, n_y, n_z, alpha);
	//x
	translation[0] = bulk_size_x;
	translation[1] = alpha * bulk_size_x;
	translation[2] = 0;
	//y
	translation[3] = alpha * bulk_size_y;
	translation[4] = bulk_size_y;
	translation[5] = 0;
	//z
	translation[6] = 0;
	translation[7] = 0;
	translation[8] = bulk_size_z / (1 - alpha * alpha);

	double E_full_stretched_C44 = full_energy(curr, n_x, n_y, n_z, potentials, true, translation);
	curr = prev;

	elastic_C44(curr, n_x, n_y, n_z, -alpha);
	//x
	translation[0] = bulk_size_x;
	translation[1] = -alpha * bulk_size_x;
	translation[2] = 0;
	//y
	translation[3] = -alpha * bulk_size_y;
	translation[4] = bulk_size_y;
	translation[5] = 0;
	//z
	translation[6] = 0;
	translation[7] = 0;
	translation[8] = bulk_size_z / (1 - alpha * alpha);

	double E_full_compressed_C44 = full_energy(curr, n_x, n_y, n_z, potentials, true, translation);
	curr = prev;

	double C44 = ((E_full_stretched_C44 - 2 * E_full + E_full_compressed_C44)  * qsi) / (2 * V_0 * alpha * alpha);
	cout << "C44 = " << C44 << endl;

	//E_sol for Co
	double E_sol = E_full - E_full_without_impurity + 4.386 + E_coh;
	cout << "E_sol = " << E_sol << endl;

	return 0;
}
