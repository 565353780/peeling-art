#pragma once
#include <Eigen/Dense>

#include "MeshViewer/MeshDefinition.h"

#include "VT_vert.h"

class ZipUp
{
public:
	ZipUp();

	~ZipUp();

	void getBoundaryIndex();

	/*useless*/
	void getBoundaryAngle();

	void getVertex_list()
	{
		vertex_list.clear();
		vertex_list.resize(mesh_.n_vertices());
		for (Mesh::VertexIter v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); ++v_it)
		{
			int id = (*v_it).idx();
			vertex_list[id] = mesh_.point(*v_it);
		}
	}

	/*计算树结构*/
	void Prim(bool plan_B = false, bool is_uniform = false);

	/*计算树分支*/
	void evaluateComponent();

	/*清理细小分支*/
	void cleanComponent(float tol);

	/*输出树分支*/
	void writeComponent(std::string a);

	/*输出边界点下标*/
	void writeBoundaryIndex(std::string a);

	/*输出边界点坐标标*/
	void writeBoundaryPosition(std::string a);

	/*输出树结构*/
	void writeBoundaryStructure(std::string a);

	/*输出树的所有边*/
	void writeTreeEdge(std::string a);

//	Mesh3D *mesh_;
	Mesh mesh_;
	std::vector<OpenMesh::Vec3d> vertex_list;

	std::vector<int> boundary_vertex_index_;

	std::vector<float> boundary_angle_;

	/*依次记录树的每条边，第一行为子节点*/
	Eigen::Matrix2Xf tree_edges_;

	/*依次记录每条边界边上，落的其他边界点*/
	std::vector<std::vector<int>> bedge_structure_;

	// 以下三者元素对应

	/*依次记录树的折线分支，记录的是局部下标*/
	std::vector<std::vector<int>> tree_component_;

	/*折线分支长度*/
	std::vector<float> component_length_;

	/*树的顶点列表*/
	std::vector<VT_vert> tvert_list_;

	/*树的根节点位置*/
	VT_vert *tree_root_;
};