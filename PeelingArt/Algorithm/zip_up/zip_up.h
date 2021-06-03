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

	/*�������ṹ*/
	void Prim(bool plan_B = false, bool is_uniform = false);

	/*��������֧*/
	void evaluateComponent();

	/*����ϸС��֧*/
	void cleanComponent(float tol);

	/*�������֧*/
	void writeComponent(std::string a);

	/*����߽���±�*/
	void writeBoundaryIndex(std::string a);

	/*����߽�������*/
	void writeBoundaryPosition(std::string a);

	/*������ṹ*/
	void writeBoundaryStructure(std::string a);

	/*����������б�*/
	void writeTreeEdge(std::string a);

//	Mesh3D *mesh_;
	Mesh mesh_;
	std::vector<OpenMesh::Vec3d> vertex_list;

	std::vector<int> boundary_vertex_index_;

	std::vector<float> boundary_angle_;

	/*���μ�¼����ÿ���ߣ���һ��Ϊ�ӽڵ�*/
	Eigen::Matrix2Xf tree_edges_;

	/*���μ�¼ÿ���߽���ϣ���������߽��*/
	std::vector<std::vector<int>> bedge_structure_;

	// ��������Ԫ�ض�Ӧ

	/*���μ�¼�������߷�֧����¼���Ǿֲ��±�*/
	std::vector<std::vector<int>> tree_component_;

	/*���߷�֧����*/
	std::vector<float> component_length_;

	/*���Ķ����б�*/
	std::vector<VT_vert> tvert_list_;

	/*���ĸ��ڵ�λ��*/
	VT_vert *tree_root_;
};