#pragma once

#include "MeshViewer/MeshDefinition.h"

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/SparseCore>
#include <Eigen/PardisoSupport>


class ARAPPara_3D
{
public:
	ARAPPara_3D();
	~ARAPPara_3D();

	inline float cotan(const Eigen::Vector3f &v1, const Eigen::Vector3f &v2)
	{
		return v1.dot(v2) / (v1.cross(v2)).norm();
	}

public:
	void evaluateCoff();

	void evaluateL();
	double solveMatrix_xy();
	double solveMatrixz();

	void CreatAllMatrix();

	void optimizeMesh();

	void Position_Up();

public:
	Mesh mesh_plane;
	Mesh mesh_surface;

private:

	/*mesh_÷–∞Î±ﬂ∂‘Ω«µƒ”‡«–*/
	std::vector<double> cot_of_edge_;

	std::vector<Eigen::Matrix3d> matrix_L_;

	int center_id;
	std::set<int> fix_id;

	Eigen::SparseMatrix<double>	matrix_Axy;
	Eigen::SparseMatrix<double>	matrix_Az;
	Eigen::PardisoLU<Eigen::SparseMatrix<double>>	solver_xy;
	Eigen::PardisoLU<Eigen::SparseMatrix<double>>	solver_z;
};

