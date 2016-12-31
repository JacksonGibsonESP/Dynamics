#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
using namespace std;

/*
//test commit
//Cu
const double lattice_constant = 3.615;
const double cutoff = 1.7 * lattice_constant;
const double atom_mass = 63.54;
const double A = 0.0855;
const double s = 1.2240;
const double p = 10.96;
const double q = 2.278;
const double r0 = lattice_constant / sqrt(2);
*/

//Ag
const double lattice_constant = 4.085;
const double cutoff = 1.7 * lattice_constant;
const double A = 0.1028;
const double s = 1.178;
const double p = 10.928;
const double q = 3.139;
const double r0 = lattice_constant / sqrt(2);

const double qsi = 0.8018993929636421;

class Atom {
public:
	double x = 0;
	double y = 0;
	double z = 0;
	double E_r = 0;
	double E_b = 0;
};

bool toofar(Atom part, Atom part2, double cutoff) {
	if (part.x == part2.x && part.y == part2.y && part.z == part2.z)
		return true;
	return sqrt((part2.x - part.x) * (part2.x - part.x)
		+ (part2.y - part.y) * (part2.y - part.y)
		+ (part2.z - part.z) * (part2.z - part.z)) > cutoff ? true : false;
}

double full_energy(vector <vector<Atom>> &contains, int n_x, int n_y, int n_z, double bulk_size_x, double bulk_size_y, double bulk_size_z) {

	double E_coh = 0;

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
							if (abs(particle2.x - particle.x) > bulk_size_x / 2) { //дальше одной клетки
								if (particle2.x > particle.x) {
									particle2.x -= bulk_size_x; //в соседнюю клетку справа перейдём
								}
								else {
									particle2.x += bulk_size_x; //в соседнюю клетку слева перейдём
								}
							}
							if (abs(particle2.y - particle.y) > bulk_size_y / 2) { //дальше одной клетки
								if (particle2.y > particle.y) {
									particle2.y -= bulk_size_y; //в соседнюю клетку справа перейдём
								}
								else {
									particle2.y += bulk_size_y; //в соседнюю клетку слева перейдём
								}
							}
							if (abs(particle2.z - particle.z) > bulk_size_z / 2) { //дальше одной клетки
								if (particle2.z > particle.z) {
									particle2.z -= bulk_size_z; //в соседнюю клетку справа перейдём
								}
								else {
									particle2.z += bulk_size_z; //в соседнюю клетку слева перейдём
								}
							}
							if (toofar(particle, particle2, cutoff)) {
								continue;
							}

							//distance
							double r = sqrt((particle2.x - particle.x) * (particle2.x - particle.x)
								+ (particle2.y - particle.y) * (particle2.y - particle.y)
								+ (particle2.z - particle.z) * (particle2.z - particle.z));

							//Coh
							particle.E_r += A * exp(-p * (r / r0 - 1));
							E_b_temp += s * s * exp(-2 * q * (r / r0 - 1));
						}
					}
				}
			}
			particle.E_b -= sqrt(E_b_temp);
			E_coh += particle.E_r + particle.E_b;
			particle.E_r = 0;
			particle.E_b = 0;
		}
	}
	return E_coh;
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
	getline(fin, dumb);

	double max_x = 0;
	double max_y = 0;
	double max_z = 0;

	while (fin.good()) {
		Atom tmp;
		fin >> tmp.x;
		fin >> tmp.x;
		if (tmp.x > max_x) {
			max_x = tmp.x;
		}
		fin >> tmp.y;
		if (tmp.y > max_y) {
			max_y = tmp.y;
		}
		fin >> tmp.z;
		if (tmp.z > max_z) {
			max_z = tmp.z;
		}
		init.push_back(tmp);
	}

	fin.close();

	int n_x = 3;
	int n_y = 3;
	int n_z = 3;
	double bulk_size_x = 3 * lattice_constant;
	double bulk_size_y = 3 * lattice_constant;
	double bulk_size_z = 3 * lattice_constant;
	/*double bulk_size_x = max_x;
	double bulk_size_y = max_y;
	double bulk_size_z = max_z;
	double step = lattice_constant;
	int n_x = int(bulk_size_x / step) + 1;
	int n_y = int(bulk_size_y / step) + 1;
	int n_z = int(bulk_size_z / step) + 1;*/

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
	
	double E_full = full_energy(prev, n_x, n_y, n_z, bulk_size_x, bulk_size_y, bulk_size_z);
	cout << "E_full = " << E_full << endl;
	double E_coh = E_full / size;
	cout << "E_coh = " << E_coh << endl;
	cout << "alpha = ";
	double alpha;
	cin >> alpha;
	
	//B
	stretch(curr, n_x, n_y, n_z, alpha);
	double E_full_stretched = full_energy(curr, n_x, n_y, n_z, bulk_size_x * (1 + alpha), bulk_size_y * (1 + alpha), bulk_size_z * (1 + alpha));
	curr = prev;

	stretch(curr, n_x, n_y, n_z, -alpha);
	double E_full_compressed = full_energy(curr, n_x, n_y, n_z, bulk_size_x * (1 - alpha), bulk_size_y * (1 - alpha), bulk_size_z * (1 - alpha));
	curr = prev;

	double V_0 = lattice_constant * lattice_constant * lattice_constant * n_x * n_y * n_z;
	double B = (2 * (E_full_stretched - 2 * E_full + E_full_compressed) * qsi) / (9 * V_0 * alpha * alpha);
	cout << "B = " << B << endl;

	//C11 C12
	elastic_C11(curr, n_x, n_y, n_z, alpha);
	double E_full_stretched_C11 = full_energy(curr, n_x, n_y, n_z, bulk_size_x * (1 + alpha), bulk_size_y * (1 + alpha), bulk_size_z);
	curr = prev;

	elastic_C11(curr, n_x, n_y, n_z, -alpha);
	double E_full_compressed_C11 = full_energy(curr, n_x, n_y, n_z, bulk_size_x * (1 - alpha), bulk_size_y * (1 - alpha), bulk_size_z);
	curr = prev;

	elastic_C12(curr, n_x, n_y, n_z, alpha);
	double E_full_stretched_C12 = full_energy(curr, n_x, n_y, n_z, bulk_size_x * (1 + alpha), bulk_size_y * (1 - alpha), bulk_size_z);
	curr = prev;

	elastic_C12(curr, n_x, n_y, n_z, -alpha);
	double E_full_compressed_C12 = full_energy(curr, n_x, n_y, n_z, bulk_size_x * (1 - alpha), bulk_size_y * (1 + alpha), bulk_size_z);
	curr = prev;

	double C11 = ((E_full_stretched_C11 - 2 * E_full + E_full_compressed_C11) / (alpha * alpha) + (E_full_stretched_C12 - 2 * E_full + E_full_compressed_C12) / (alpha * alpha)) * qsi / (2 * V_0);
	cout << "C11 = " << C11 << endl;

	double C12 = ((E_full_stretched_C11 - 2 * E_full + E_full_compressed_C11) / (alpha * alpha) - (E_full_stretched_C12 - 2 * E_full + E_full_compressed_C12) / (alpha * alpha)) * qsi / (2 * V_0);
	cout << "C12 = " << C12 << endl;

	//C44
	elastic_C44(curr, n_x, n_y, n_z, alpha);
	double E_full_stretched_C44 = full_energy(curr, n_x, n_y, n_z,
		bulk_size_x + alpha * bulk_size_y,
		alpha * bulk_size_x + bulk_size_y,
		bulk_size_z / (1 - alpha * alpha));
	curr = prev;

	elastic_C44(curr, n_x, n_y, n_z, -alpha);
	double E_full_compressed_C44 = full_energy(curr, n_x, n_y, n_z,
		bulk_size_x - alpha * bulk_size_y,
		-alpha * bulk_size_x + bulk_size_y,
		bulk_size_z / (1 - alpha * alpha));
	curr = prev;

	double C44 = ((E_full_stretched_C44 - 2 * E_full + E_full_compressed_C44)  * qsi) / (2 * V_0 * alpha * alpha);
	cout << "C44 = " << C44 << endl;

	return 0;
}
