#pragma once

#include "MeshViewer/MeshDefinition.h"

#include <Eigen/Dense>

#include "Mesh/surface_manager.h"


class ARAPParametricManger
{
public:
	ARAPParametricManger();

	~ARAPParametricManger();

public:
	//----------------------------------------------------------------------------------------------
	//�����и�����߽�
	//----------------------------------------------------------------------------------------------
	int  cutmeshBoundary_Inline();

	void cutmeshBoundary_fix(Mesh& ref_mesh, int times, std::pair < int, int >& before_end_id);

	void compute_cut_Boundary_id(); //����߽类�и�ɼ���

	void change_cut_ref(Mesh& ref_mesh);

public:
	Mesh mesh_;
	Mesh cut_mesh_ref;

	std::set< int > fixed_id;
	
	double area_lambda;
	int ref_fix_id;

private:	
	std::vector< std::vector< int > > all_boundary_id; //all boundary part id

	
};
