#pragma once
#include<vector>

class Atom {
public:
	double x = 0;
	double y = 0;
	double z = 0;
	double f_x = 0;
	double f_y = 0;
	double f_z = 0;
	double v_x = 0;
	double v_y = 0;
	double v_z = 0;
	double mass = 107.8682;
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

class CalcData {
public:
	//without impurities
	std::vector <std::vector<Atom>> bulk;
	int bulk_atom_count;
	//with one impurity for E_sol
	std::vector <std::vector<Atom>> imp;
	//with dimer in surf
	std::vector <std::vector<Atom>> dim_in_surf;
	//with adatom in surf
	std::vector <std::vector<Atom>> adatom_in_surf;
	//with dimer on surf
	std::vector <std::vector<Atom>> dim_on_surf;
	//with adatom on surf
	std::vector <std::vector<Atom>> adatom_on_surf;

	Translation tr;
	int n_x, n_y, n_z;
	double alpha = 0.01;
	double lattice_constant;
	double cutoff;
	RGL potentials[3];
	double V_0;
	double E_full;
	double qsi = 0.8018993929636421;
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
};
