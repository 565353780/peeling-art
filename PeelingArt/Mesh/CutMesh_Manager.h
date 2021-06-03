#pragma once
#include "MeshViewer/MeshDefinition.h"
#include <Eigen/Dense>
#include <vector>

class BPolyline;

class CutMesh_Manager
{
public:
	CutMesh_Manager();
	~CutMesh_Manager();

public:
	void InlineBoundary();

	void getColor(int ci, double& R, double& G, double& B);

public:
	Mesh cut_mesh_a; //part surface mesh
	Mesh cut_mesh_b; //part parametric mesh, no boundary condition;
	std::vector<Mesh> cut_mesh_b_list;
	std::vector<bool> is_show_list;

	std::vector<BPolyline*> cut_b_boundary;

	std::vector<std::pair<int, int> >cut_b_endpoints_id;
	std::vector<std::pair<Eigen::Vector3f, Eigen::Vector3f> >cut_b_endpoints;
	
	int is_over_load;
	int user_select;
};

