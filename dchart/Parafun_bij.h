#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include<math.h>
#include<vector>
#include<set>
#include<iostream>
#include<fstream>
#include"Eigen/Dense"
#include"Eigen/Sparse"
#include "dchart/ScafData.h"
#include "dchart/PardisoSolver.h"
#include<time.h>

using namespace Eigen;
using namespace std;

double get_smallest_pos_quad_zero_bij(double a, double b, double c);
class Parafun_bij
{

public:
	Parafun_bij(ScafData & data) :d_(data) {
		pardiso = NULL;
	};
	~Parafun_bij();


	void after_mesh_improve();

	//void BPE();
	//void calc_gradient_norm(const VectorXd &x);

	void init();
	void init_area();

	void run_bpe();
	double perform_iteration_cm(bool use_CM);
	void id_vs_index();

	//void Update_source_same_t();

	//void Tutte();
	void Pre_calculate();

	void CM();
	void SLIM();

	void Energy(const VectorXd &x, double &energy);
	void Energysource();
	void MyEnergySource(std::vector<double>& face_enery);
	double compute_energy(const Eigen::MatrixXd & x, bool whole=false);
	void adjust_scaf_weight(double new_weight);
	void handle_mintri();


	double newton_equation(const double & a, const double & b, const double & K);

	void backtracking_line_search(const VectorXd &x, const VectorXd &d, const VectorXd &negetive_grad, double &alpha);

	void local_coordinate_inverse(int i, double &p00, double &p01, double &p10, double &p11);
	void local_coordinate_inverse_scaf(int i, double &p00, double &p01, double &p10, double &p11);



	void max_step(const VectorXd &xx, const VectorXd &dd, double &step);


	ScafData &d_;
	MatrixXd fix_pos;
	int total_num;
	double area_threshold;
	vector<int> var_ids;
	vector<int> id2index;

	double time_consumption;

	int F_N;
	int V_N;
	int dim = 2;


	vector<double> area;
	vector<double> area_scaf;
	vector<double> area_src;


	vector<double> source_p00;
	vector<double> source_p01;
	vector<double> source_p10;
	vector<double> source_p11;

	vector<int> F0;
	vector<int> F1;
	vector<int> F2;

	PardisoSolver* pardiso;
	vector<int> pardiso_ia;
	vector<int> pardiso_ja;
	vector<double> pardiso_a;
	vector<double> pardiso_b;

	double energy_uniform;
	double energy_area;

	VectorXd position_of_mesh;

	vector<int> id_h00; vector<int> id_h01; vector<int> id_h02; vector<int> id_h03; vector<int> id_h04; vector<int> id_h05;
	vector<int> id_h11; vector<int> id_h12; vector<int> id_h13; vector<int> id_h14; vector<int> id_h15;
	vector<int> id_h22; vector<int> id_h23; vector<int> id_h24; vector<int> id_h25;
	vector<int> id_h33; vector<int> id_h34; vector<int> id_h35;
	vector<int> id_h44; vector<int> id_h45;
	vector<int> id_h55;
};

