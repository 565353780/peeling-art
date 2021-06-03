#pragma once
#include <Eigen/Dense>
#include <list>
#include "MeshViewer/MeshDefinition.h"

struct VFrame
{
	OpenMesh::Vec3d e0;
	OpenMesh::Vec3d e1;
	OpenMesh::Vec3d n;
};

class Surface_Manager
{
public:
	Surface_Manager();

	~Surface_Manager();

	enum {
		MESH, SURFACE
	};

	/*����һ��pt���Ļ��ܣ�������������r_u,r_v,n*/
	void frameAtPoint(Mesh::Point& pt, VFrame& frm);

	/*pt���ֲ�����ϵ�У�����x��ʾ�ĵ�*/
	Mesh::Point PointWithCoordinate(Mesh::Point &pt, VFrame &frm, const Eigen::Vector2d &x);
	void PointProjection(Mesh::Point& pt, Mesh::Point& proj_p);

	int surface_type;

private:
	void CompressionProjiect(Mesh::Point& p, Mesh::Point& new_p);

	void NormAt(Mesh::Point& p, OpenMesh::Vec3d& n);

};