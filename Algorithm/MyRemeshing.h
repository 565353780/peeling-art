#pragma once
#include "MeshViewer/MeshDefinition.h"
#include <Eigen/Dense>

class Surface_Manager;

class MyRemeshing
{
public:
	MyRemeshing();
	~MyRemeshing();

public:
	std::vector<int> my_edge_cost;

	void My_Remeshing(double mu, Mesh& mesh);

	void My_Remeshing_Out_Part(double mu, Mesh& mesh);
	bool My_edge_Collapse_Out(Mesh::EdgeHandle& eh, Mesh& mesh);
	bool My_edge_Split_Out(Mesh::EdgeHandle& eh, Mesh& mesh);
	bool My_edge_Filp(Mesh::EdgeHandle& eh, Mesh& mesh);

	void getVerticesArea(std::vector< float > &vertice_area, std::vector< float >& arealist, Mesh& mesh);

public:
	Surface_Manager* mesh_suface_;
};

