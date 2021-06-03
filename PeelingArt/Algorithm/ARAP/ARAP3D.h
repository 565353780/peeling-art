#pragma once
//#define EIGEN_USE_MKL_ALL

#include "MeshViewer/MeshDefinition.h"

#include <Eigen/Dense>
#include <Eigen/StdVector>

class Surface_Manager;
struct VFrame;

class ARAP3D
{
public:
	ARAP3D();

	~ARAP3D();

	inline float cotan(const Eigen::Vector3f &v1, const Eigen::Vector3f &v2 )
	{
		return v1.dot(v2) / (v1.cross(v2)).norm();
	}

	inline OpenMesh::Vec4d Matrix_product(const OpenMesh::Vec4d& M1, const OpenMesh::Vec4d& M2)
	{
		double a = M1[0] * M2[0] + M1[1] * M2[2];
		double b = M1[0] * M2[1] + M1[1] * M2[3];
		double c = M1[2] * M2[0] + M1[3] * M2[2];
		double d = M1[2] * M2[1] + M1[3] * M2[3];

		return OpenMesh::Vec4d(a, b, c, d);
	}
	inline OpenMesh::Vec2d Matrix_product_vector(const OpenMesh::Vec4d& M1, const OpenMesh::Vec2d& V2)
	{
		double a = M1[0] * V2[0] + M1[1] * V2[1];
		double b = M1[2] * V2[0] + M1[3] * V2[1];

		return OpenMesh::Vec2d(a, b);
	}
	inline OpenMesh::Vec4d Matrix_transpose(const OpenMesh::Vec4d& M1)
	{
		return OpenMesh::Vec4d(M1[0], M1[2], M1[1], M1[3]);
	}

public:
	double myfun(double s)
	{
		return (s - 1)*(s - 1);
	}

public:

	/*����ֻ����������ڵ�λ�õ�ϵ��*/
	void evaluateCoff();
	
	/*����frame, Lt, Rt*/
	void evaluateFLR();

	/*���ɾ���*/
	bool createMatrix();

	/*�����Ż�*/
	void optimizeMesh();

//	void computeEnergy(double area_lambda);

public:
	/*�ض���Ĳο�λ��*/
	std::map<int, Eigen::Vector3f> ref_map_;
	/*�����Ƿ�Ӧ�ڵ�������������*/
	std::vector<bool> is_singular_face_;

	Mesh mesh_;
	Surface_Manager* mesh_surface_;

	std::vector<double> energy;

	double shrink_weight;
	double special_weight;

private:
	
	/*mesh_�а�߶Խǵ�����*/
	std::vector<double> cot_of_hedge_;

	/*mesh_�а�����Խڵ��Ŀ����ƽ������*/
	std::vector<OpenMesh::Vec2d> x_of_hedge_new;

	/*mesh_�а�����Խڵ��ƽ������*/
	std::vector<OpenMesh::Vec2d> u_of_hedge_new;

	/*��ǰ��������Lt���������˳������*/
	std::vector<OpenMesh::Vec4d> matrix_L_new;

	/*��ǰ��������Rt_ij���󣬣��԰�ߵ�˳������*/
	std::vector<OpenMesh::Vec4d> matrix_R_new;

	/*mesh_�нڵ㴦��һ���������*/
	std::vector<VFrame> vert_frame;

	/*�ݶȷ���*/
	Eigen::VectorXd	vector_du_;
};