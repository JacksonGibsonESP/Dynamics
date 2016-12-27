#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
using namespace std;

//const int iterations_count = 10;
//const double dt = 1e-3;
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

// double D = 933;//eV
//const double a = 3.615 / 1.44;

const double qsi = 0.8018993929636421;
class Atom {
public:
	double x = 0;
	double y = 0;
	double z = 0;
	double E_r = 0;
	double E_b = 0;
	double f_x = 0;
	double f_y = 0;
	double f_z = 0;
	double v_x = 0;
	double v_y = 0;
	double v_z = 0;
	double mass = 0;
	/*Atom() {
	}
	Atom(double mass) {
	Atom();
	this->mass = mass;
	}*/
};

bool toofar(Atom part, Atom part2, double cutoff) {
	if (part.x == part2.x && part.y == part2.y && part.z == part2.z)
		return true;
	return sqrt((part2.x - part.x) * (part2.x - part.x)
		+ (part2.y - part.y) * (part2.y - part.y)
		+ (part2.z - part.z) * (part2.z - part.z)) > cutoff ? true : false;
}

/*void calc(vector <vector<Atom>> &contains, int n_x, int n_y, int n_z, double cell_size_x, double cell_size_y, double cell_size_z) {
	for (int i = 0; i < n_x * n_y * n_z; ++i) {
		for (Atom& particle : contains[i]) {
			particle.f_x = 0;
			particle.f_y = 0;
			particle.f_z = 0;
			int atomcell_x = i % n_x;
			int atomcell_y = i / n_x % n_y;
			int atomcell_z = i / (n_x * n_y);
			double F_r = 0;
			double F_b_1 = 0;
			double F_b_2 = 0;
			double F_b_res = 0;
			for (int dx = -1; dx <= 1; ++dx) {
				for (int dy = -1; dy <= 1; ++dy) {
					for (int dz = -1; dz <= 1; ++dz) {
						for (Atom particle2 : contains[
							(atomcell_x + n_x + dx) % n_x +
								(atomcell_y + n_y + dy) % n_y * n_x +
								(atomcell_z + n_z + dz) % n_z * n_y * n_x]) {
							//check for translation
							if (abs(particle2.x - particle.x) > cell_size_x / 2) { //äàëüøå îäíîé êëåòêè
								if (particle2.x > particle.x) {
									particle2.x -= cell_size_x; //â ñîñåäíþþ êëåòêó ñïðàâà ïåðåéä¸ì
								}
								else {
									particle2.x += cell_size_x; //â ñîñåäíþþ êëåòêó ñëåâà ïåðåéä¸ì
								}
							}
							if (abs(particle2.y - particle.y) > cell_size_y / 2) { //äàëüøå îäíîé êëåòêè
								if (particle2.y > particle.y) {
									particle2.y -= cell_size_y; //â ñîñåäíþþ êëåòêó ñïðàâà ïåðåéä¸ì
								}
								else {
									particle2.y += cell_size_y; //â ñîñåäíþþ êëåòêó ñëåâà ïåðåéä¸ì
								}
							}
							if (abs(particle2.z - particle.z) > cell_size_z / 2) { //äàëüøå îäíîé êëåòêè
								if (particle2.z > particle.z) {
									particle2.z -= cell_size_z; //â ñîñåäíþþ êëåòêó ñïðàâà ïåðåéä¸ì
								}
								else {
									particle2.z += cell_size_z; //â ñîñåäíþþ êëåòêó ñëåâà ïåðåéä¸ì
								}
							}
							if (toofar(particle, particle2, cutoff)) {
								continue;
							}

							//distance
							double r = sqrt((particle2.x - particle.x) * (particle2.x - particle.x)
								+ (particle2.y - particle.y) * (particle2.y - particle.y)
								+ (particle2.z - particle.z) * (particle2.z - particle.z));
							/*
							//RGL potential
							double F = A * p * exp(p * (1 - r / r0)) / r0
							- 0.5 / sqrt(s * s * exp(-2 * q * (r / r0 - 1))) * 2 * q * s * s * exp(2 * q * (1 - r / r0)) / r0;
							double dir_x = particle2.x - particle.x;
							double dir_y = particle2.y - particle.y;
							double dir_z = particle2.z - particle.z;
							double len = sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
							dir_x /= len;
							dir_y /= len;
							dir_z /= len;
							dir_x *= F;
							dir_y *= F;
							dir_z *= F;
							particle.f_x -= dir_x;
							particle.f_y -= dir_y;
							particle.f_z -= dir_z;
							*/
							/*//Ïîòåíöèàë Ëåííàðä-Äæîíñà. Ðàáîòàåò
							double F = 12 * D / a * (pow(a / r, 13) - pow(a / r, 7));

							double dir_x = particle2.x - particle.x;
							double dir_y = particle2.y - particle.y;
							double dir_z = particle2.z - particle.z;
							double len = sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);
							dir_x /= len;
							dir_y /= len;
							dir_z /= len;
							dir_x *= F;
							dir_y *= F;
							dir_z *= F;

							particle.f_x -= dir_x;
							particle.f_y -= dir_y;
							particle.f_z -= dir_z;
						}
					}
				}
			}
		}
	}
}*/

void dump(string filename, unsigned int size, vector <vector<Atom>> &contains) {
	ofstream out(filename);
	out << size << '\n';
	out << "one frame of relaxation\n";
	for (unsigned int i = 0; i < contains.size(); ++i) {
		for (unsigned int j = 0; j < contains[i].size(); ++j) {
			//out << setw(15) << contains[i][j].x << setw(15) << contains[i][j].y << setw(15) << contains[i][j].z << '\n';
			out << contains[i][j].x << ' ' << contains[i][j].y << ' ' << contains[i][j].z << '\n';
		}
	}
	out.close();
}

double full_energy(vector <vector<Atom>> &contains, int n_x, int n_y, int n_z, double cell_size_x, double cell_size_y, double cell_size_z) {

	double E_coh = 0;

	for (int i = 0; i < n_x * n_y * n_z; ++i) {
		for (Atom& particle : contains[i]) {
			particle.f_x = 0;
			particle.f_y = 0;
			particle.f_z = 0;
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
							if (abs(particle2.x - particle.x) > cell_size_x / 2) { //äàëüøå îäíîé êëåòêè
								if (particle2.x > particle.x) {
									particle2.x -= cell_size_x; //â ñîñåäíþþ êëåòêó ñïðàâà ïåðåéä¸ì
								}
								else {
									particle2.x += cell_size_x; //â ñîñåäíþþ êëåòêó ñëåâà ïåðåéä¸ì
								}
							}
							if (abs(particle2.y - particle.y) > cell_size_y / 2) { //äàëüøå îäíîé êëåòêè
								if (particle2.y > particle.y) {
									particle2.y -= cell_size_y; //â ñîñåäíþþ êëåòêó ñïðàâà ïåðåéä¸ì
								}
								else {
									particle2.y += cell_size_y; //â ñîñåäíþþ êëåòêó ñëåâà ïåðåéä¸ì
								}
							}
							if (abs(particle2.z - particle.z) > cell_size_z / 2) { //äàëüøå îäíîé êëåòêè
								if (particle2.z > particle.z) {
									particle2.z -= cell_size_z; //â ñîñåäíþþ êëåòêó ñïðàâà ïåðåéä¸ì
								}
								else {
									particle2.z += cell_size_z; //â ñîñåäíþþ êëåòêó ñëåâà ïåðåéä¸ì
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

void stretch(vector <vector<Atom>> &contains, int n_x, int n_y, int n_z, double cell_size_x, double cell_size_y, double cell_size_z, double alpha) {

	double E_coh = 0;

	for (int i = 0; i < n_x * n_y * n_z; ++i) {
		for (Atom& particle : contains[i]) {
			particle.x = particle.x + alpha * particle.x;
			particle.y = particle.y + alpha * particle.y;
			particle.z = particle.z + alpha * particle.z;
		}
	}
}

void elastic(vector <vector<Atom>> &contains, int n_x, int n_y, int n_z, double cell_size_x, double cell_size_y, double cell_size_z, double alpha) {

	double E_coh = 0;

	for (int i = 0; i < n_x * n_y * n_z; ++i) {
		for (Atom& particle : contains[i]) {
			particle.x = particle.x + alpha * particle.x;
			particle.y = particle.y + alpha * particle.y;
		}
	}
}

int main(int argc, char* argv[]) {

	vector<Atom> init;

	/*if (argc == 1 || argc >= 4) {
		cout << "Incorrect number of arguments\n";
		return 0;
	}

	if (strcmp(argv[2], "-l") != 0) {
		cout << "Incorrect key\n";
	}*/

	ifstream fin(argv[1], ifstream::in);

	if (!fin.is_open()) {
		cout << "File can not be opened\n";
		return 0;
	}

	string dumb;
	//double cell_size_x = 0;
	//double cell_size_y = 0;
	//double cell_size_z = 0;
	//int xy = 0;
	//int z = 0;
	//if (argc == 3) {
		//getline(fin, dumb);
		//fin >> cell_size_x;
		//fin >> cell_size_y;
		//fin >> cell_size_z;
		//cell_size_z /= 2;
		//getline(fin, dumb);
		//fin >> xy >> z;
	//}

	//int size = xy * z;

	int size = 0;
	fin >> size;
	getline(fin, dumb);


	/*double cell_x = 0;
	double max_y = 0;
	double max_z = 0;*/

	while (fin.good()) {
		Atom tmp;
		fin >> tmp.x;
		fin >> tmp.x;
		/*if (tmp.x > max_x) {
			max_x = tmp.x;
		}*/
		fin >> tmp.y;
		/*if (tmp.y > max_y) {
			max_y = tmp.y;
		}*/
		fin >> tmp.z;
		/*if (tmp.z > max_z) {
			max_z = tmp.z;
		}*/
		/*if (tmp.z > cell_size_z - a / 10) {
		getline(fin, dumb);
		continue;
		}*/
		init.push_back(tmp);
		//getline(fin, dumb);
	}

	fin.close();

	/*int n_x = (int)(max_x / step + 1);
	int n_y = (int)(max_y / step + 1);
	int n_z = (int)(max_z / step + 1);
	double cell_x = max_x / n_x;
	double cell_y = max_y / n_y;
	double cell_z = max_z / n_z;*/

	int n_x = 3;
	int n_y = 3;
	int n_z = 3;
	double cell_x = 3 * lattice_constant;
	double cell_y = 3 * lattice_constant;
	double cell_z = 3 * lattice_constant;

	int step = lattice_constant;

	vector <vector<Atom>> prev(n_x * n_y * n_z);
	vector <vector<Atom>> curr(n_x * n_y * n_z);
	for (unsigned int i = 0; i < init.size(); ++i) {
		Atom tmp = init[i];
		int x = int(tmp.x / step);
		int y = int(tmp.y / step);
		int z = int(tmp.z / step);
//		tmp.mass = atom_mass;
		prev[z * n_y * n_x + y * n_x + x].push_back(tmp);
		curr[z * n_y * n_x + y * n_x + x].push_back(tmp);
	}
	/*
	//calculate prev forces
	calc(prev, n_x, n_y, n_z, cell_size_x, cell_size_y, cell_size_z);
	//dump iteration
	cout << "frame_0.txt\n";
	dump("frame_0.txt", init.size(), curr);
	cout << "E_coh = " << coh_energy(curr, n_x, n_y, n_z, cell_size_x, cell_size_y, cell_size_z) << endl;
	int iterations_count;
	cout << "iterations count: ";
	cin >> iterations_count;
	double dt;
	cout << "dt: ";
	cin >> dt;
	for (int it_cnt = 0; it_cnt < iterations_count; ++it_cnt) {
	//count positions
	for (int i = 0; i < n_x * n_y * n_z; ++i) {
	for (unsigned int j = 0; j < prev[i].size(); ++j) {
	curr[i][j].x = prev[i][j].x + prev[i][j].v_x * dt + prev[i][j].f_x * dt * dt / (2 * curr[i][j].mass);
	curr[i][j].y = prev[i][j].y + prev[i][j].v_y * dt + prev[i][j].f_y * dt * dt / (2 * curr[i][j].mass);
	curr[i][j].z = prev[i][j].z + prev[i][j].v_z * dt + prev[i][j].f_z * dt * dt / (2 * curr[i][j].mass);
	}
	}
	//calculate
	calc(curr, n_x, n_y, n_z, cell_size_x, cell_size_y, cell_size_z);
	//count velocities
	for (int i = 0; i < n_x * n_y * n_z; ++i) {
	for (unsigned int j = 0; j < curr[i].size(); ++j) {
	curr[i][j].v_x = prev[i][j].v_x + (prev[i][j].f_x + curr[i][j].f_x) / (2 * curr[i][j].mass) * dt;
	curr[i][j].v_y = prev[i][j].v_y + (prev[i][j].f_y + curr[i][j].f_y) / (2 * curr[i][j].mass) * dt;
	curr[i][j].v_z = prev[i][j].v_z + (prev[i][j].f_z + curr[i][j].f_z) / (2 * curr[i][j].mass) * dt;
	}
	}
	//dump iteration
	cout << "frame_" << it_cnt + 1 << ".txt\n";
	dump("frame_" + to_string(it_cnt + 1) + ".txt", init.size(), curr);
	cout << "E_coh = " << coh_energy(curr, n_x, n_y, n_z, cell_size_x, cell_size_y, cell_size_z) << endl;
	//swap
	prev.swap(curr);
	//count full energy and break
	}
	*/
	double E_full = full_energy(prev, n_x, n_y, n_z, cell_x, cell_y, cell_z);
	cout << "E_full = " << E_full << endl;
	double E_coh = E_full / size;
	cout << "E_coh = " << E_coh << endl;
	double alpha = 0.01;
	cout << "alpha = " << alpha << endl;
	//B
	stretch(curr, n_x, n_y, n_z, cell_x, cell_y, cell_z, alpha);
	double E_full_stretched = full_energy(curr, n_x, n_y, n_z, cell_x * (1 + alpha), cell_y * (1 + alpha), cell_z * (1 + alpha));
	cout << "E_full_stretched = " << E_full_stretched << endl;
	curr = prev;

	stretch(curr, n_x, n_y, n_z, cell_x, cell_y, cell_z, -alpha);
	double E_full_compressed = full_energy(curr, n_x, n_y, n_z, cell_x * (1 - alpha), cell_y * (1 - alpha), cell_z * (1 - alpha));
	cout << "E_full_compressed = " << E_full_compressed << endl;
	curr = prev;

	double V_0 = lattice_constant * lattice_constant * lattice_constant * n_x * n_y * n_z;
	double B = (2 * (E_full_stretched - 2 * E_full + E_full_compressed) * qsi) / (9 * V_0 * alpha * alpha);
	cout << "B = " << B << endl;
	curr = prev;
	//Нелдера-мида деформируемого многогранника
	//C11
	elastic(curr, n_x, n_y, n_z, cell_x, cell_y, cell_z, alpha);
	E_full_stretched = full_energy(curr, n_x, n_y, n_z, cell_x * (1 + alpha), cell_y * (1 + alpha), cell_z);
	cout << "E_full_stretched = " << E_full_stretched << endl;
	curr = prev;

	elastic(curr, n_x, n_y, n_z, cell_x, cell_y, cell_z, -alpha);
	E_full_compressed = full_energy(curr, n_x, n_y, n_z, cell_x * (1 - alpha), cell_y * (1 - alpha), cell_z);
	cout << "E_full_compressed = " << E_full_compressed << endl;
	curr = prev;



	return 0;
}
