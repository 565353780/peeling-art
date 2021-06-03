#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include<math.h>
#include<vector>
#include<iostream>
#include<fstream>
#include"Eigen\src\Eigenvalues\EigenSolver.h"
#include"Eigen/Dense"
#include"Eigen/Sparse"
#include"Eigen/SparseLU"
#include "Eigen/SVD"
#include "MeshViewer/MeshDefinition.h"
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "dchart/PardisoSolver.h"

using namespace Eigen;
using namespace std;
double get_smallest_pos_quad_zero(double a, double b, double c);

class Parafun
{
public:
	Parafun();
	Parafun(Mesh& mesh_, int p_type, vector<double>& s_p00, vector<double>& s_p01, vector<double>& s_p10, vector<double>& s_p11);
	Parafun(Mesh& mesh_, int p_type);
	~Parafun();
	bool ChangeVertPoint(Mesh& out_mesh);
	bool ChangeVertRef(Mesh& out_mesh, int zero_idx);

	void calc_gradient_norm(const VectorXd &x);

	void init_Deformation(std::set<int>& my_fix_id, std::vector<int>& my_handle_id, std::vector<int>& my_face_in_out);
	void init_handle(std::vector<double>& handle_new_position_x, std::vector<double>& handle_new_position_y);

	void init_Parametric(std::set<int>& my_fix_id = std::set<int>() );

	void recover_to_src();
	void init_area();

	bool run_cm();

	void Update_source_same_t();

	void Tutte();
	void Pre_calculate(Eigen::MatrixX3i& s_T);

	void CM();
	bool SLIM();

	void Energy(const VectorXd &x, double &energy);
	bool Energysource();

	double newton_equation(const double & a, const double & b, const double & K);

	void backtracking_line_search(const VectorXd &x, const VectorXd &d, const VectorXd &negetive_grad, double &alpha);

	void local_coordinate_inverse(int i, double &p00, double &p01, double &p10, double &p11);

	void max_step(const VectorXd &xx, const VectorXd &dd, double &step);

	void id_vs_index();


	Mesh mesh;
	double Intp_T_Min;
	double changetocm_flag;
	double convgence_con_rate;
	double time_consumption;
	int MAX_ITER_NUM;

	//要进行的操作标识：1.Parametric; 2.Deformation with handel and fix; 3.Paranetrize with fix;
	int process_type;
	double originmesh_area_sqrt;

	int total_num;

	vector<int> var_ids;
	vector<int> id2index;

	int			fix_number;
	VectorXi	fixed_ids;
	int dim = 2;

	vector<int> handle_ids;
	MatrixXd handle_target_pos;  //2D

	double handle_weight=150;
	double scaf_weight = 0.001;
	
	vector<int> face_in_out;
	int F_N;
	int V_N;

	VectorXd negative_grad_norm;
	double g_norm;

	vector<double> area;
	vector<double> area_uniform;
	vector<double> area_src;

	vector<double> source_p00;
	vector<double> source_p01;
	vector<double> source_p10;
	vector<double> source_p11;

	vector<double> update_p00;
	vector<double> update_p01;
	vector<double> update_p10;
	vector<double> update_p11;

	vector<vector<double>> source_fixed_p;

	vector<int> F0;
	vector<int> F1;
	vector<int> F2;

	PardisoSolver* pardiso;
	vector<int> pardiso_ia;
	//vector<int> pardiso_i;
	vector<int> pardiso_ja;
	vector<double> pardiso_a;
	vector<double> pardiso_b;

	double energy_uniform;
	double energy_area;

	double bound_distortion_K;
	VectorXd position_of_mesh;

	vector<int> id_h00; vector<int> id_h01; vector<int> id_h02; vector<int> id_h03; vector<int> id_h04; vector<int> id_h05;
	vector<int> id_h11; vector<int> id_h12; vector<int> id_h13; vector<int> id_h14; vector<int> id_h15;
	vector<int> id_h22; vector<int> id_h23; vector<int> id_h24; vector<int> id_h25;
	vector<int> id_h33; vector<int> id_h34; vector<int> id_h35;
	vector<int> id_h44; vector<int> id_h45;
	vector<int> id_h55;
};

