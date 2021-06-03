#pragma once
#include <Eigen/Dense>
#include <vector>
#include "MeshViewer/MeshDefinition.h"

class BPolyline;
class Surface_Manager;

class MeshGeneration
{
public:
	MeshGeneration();
	~MeshGeneration();

public:
    void meshFromUserPolyline(std::vector<Eigen::Vector3f>& use_point, Mesh& fix_mesh);
	
	void meshFromChangePolyline(BPolyline* bp, Mesh& ref_mesh);

	void PlaneMesh_To_Surface(int id);
	void ReSetEye(Mesh& mymesh);

	void mesh_scaffold_Inline(BPolyline* bp, Mesh& fix_mesh, Mesh& fin_mesh);

	
	void mesh_scaffold_Inline_new(BPolyline* bp, Mesh& fix_mesh, Mesh& fin_mesh, std::vector<std::vector<int>>& insert_list);

	Eigen::Vector3f plane_eye_position;
	BPolyline* user_pc;

	Mesh mesh;
    bool property_added = false;

	Surface_Manager* surface_manager;

	void meshFromPolyline(BPolyline* bp, Mesh& mesh_);
	void meshFromPolylineList(std::vector<BPolyline* >& bp_list);

private:
	void divide_line(Eigen::Vector3f& line_start, Eigen::Vector3f& line_end, int n, std::vector<Eigen::Vector3f>& dive_point);
	void divide_line_nobe(Eigen::Vector3f& line_start, Eigen::Vector3f& line_end, int n, std::vector<Eigen::Vector3f>& dive_point);

private:
	void meshInOut();
	void InMeshInline(Mesh& in_mesh);

	void SharpMeshInit();

	void plane_to_Octahedron();
	void WeldingOctahedron_to_Sharp();
	bool find_first(int &id, std::vector<int>& first_list, std::vector<int> second_list);

	void Project_Surface();

	void DelateWrongPoint();

	void FigureMesh();
	

	Eigen::Vector3d m_center;
	double half_length;

	double use_curve_area;

	int n_boundary;

	BPolyline* rec_boundary;

public:
	void meshFromUserPolyline_picture(std::vector<Eigen::Vector3f>& use_point, Mesh& fix_mesh);
};

