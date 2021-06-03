#pragma once
#include "dchart/PardisoSolver.h"
#include <Eigen/Dense>
#include <set>
#include <list>

struct MeshStructure
{
	MeshStructure(std::vector<std::vector<int>> &adj_,
		std::vector<std::vector<int>> &FVs_,
		std::vector<std::vector<int>> &VVs_,
		std::vector<std::vector<int>> &VFs_,
		Eigen::MatrixXd &Vs_,
		std::vector<Eigen::Vector3d> &normals_,
		std::vector<double> &areas_) :adj(adj_), FVs(FVs_), VVs(VVs_), VFs(VFs_), Vs(Vs_),
		normals(normals_), areas(areas_) {};

	std::vector<std::vector<int>> &adj;
	std::vector<std::vector<int>> &FVs;
	std::vector<std::vector<int>> &VVs;
	std::vector<std::vector<int>> &VFs;
	Eigen::MatrixXd &Vs;
	std::vector<Eigen::Vector3d> &normals;
	std::vector<double> &areas;

};

struct Cut
{
	Cut(const int& c1, const int& c2, const list<int>& br_, const int& f, const int& b) :c1_id(c1), c2_id(c2), s_es(br_), front_border(f), back_border(b) {};

	int c1_id;
	int c2_id;
	list<int> s_es;
	int front_border;
	int back_border;
};

using namespace std;
void orderTheBoundary(vector<list<int>>& order_boundary, const vector<int>& nextlocation);
bool growFromP(int p, set<int>& isused, list<int>& order_boundary, const vector<int>& nextlocation);

void map_vertices_to_circle(const Eigen::MatrixXd& V, const Eigen::VectorXi& bnd, Eigen::MatrixXd& UV);

void Tutte(int V_N, const Eigen::MatrixXi& F, const Eigen::VectorXi& bnd, const Eigen::MatrixXd& bnd_uv, Eigen::MatrixXd & uv_init);
void preCalc_pardiso(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, PardisoSolver & pardiso);

void boundary_loop(const Eigen::MatrixXi &F_ref, std::vector<std::vector<int>>& boundaryEdges);