#include "surface_manager.h"

#include <Eigen/Geometry>
#include <iostream>


Surface_Manager::Surface_Manager()
{
	surface_type = SURFACE;
}

Surface_Manager::~Surface_Manager()
{

}

void Surface_Manager::frameAtPoint(Mesh::Point &pt, VFrame& frm)
{
	OpenMesh::Vec3d n;
	assert(surface_type == SURFACE);
	NormAt(pt, n);
	frm.e0 = pt;
	frm.e1 = OpenMesh::cross(n,pt);
	frm.n = n;
}

Mesh::Point Surface_Manager::PointWithCoordinate(Mesh::Point &pt, VFrame &frm, const Eigen::Vector2d &x)
{
	Mesh::Point fin = pt + x(0) * frm.e0 + x(1) * frm.e1;
	Mesh::Point close_p;

	assert(surface_type == SURFACE);
	CompressionProjiect(fin, close_p);

	return close_p;
}

void Surface_Manager::PointProjection(Mesh::Point& pt, Mesh::Point& proj_p)
{
	assert(surface_type == SURFACE);
	CompressionProjiect(pt, proj_p);
}

void Surface_Manager::CompressionProjiect(Mesh::Point& p, Mesh::Point& new_p)
{
	double lambda = 0.76;
	double a = 1.618;

	double pr = p.norm();
	double prxy = sqrt(p[0] * p[0] + p[1] * p[1]);
	double sin_phi = p[2] / pr;
	double cos_phi = prxy / pr;
	double cos_theta = p[0] / prxy;
	double sin_theta = p[1] / prxy;
	if (prxy < 1e-8)
	{
		cos_phi = 0;
		cos_theta = 0;
		sin_theta = 0;
	}

	double r = (1 - lambda)* pow(cos_phi, a) + lambda;
	double new_x = r*cos_theta*cos_phi;
	double new_y = r*sin_theta*cos_phi;
	double new_z = r*sin_phi;

	new_p = Mesh::Point(new_x, new_y, new_z);
}

void Surface_Manager::NormAt(Mesh::Point& p, OpenMesh::Vec3d& n)
{
	double lambda = 0.76;
	double a = 1.618;

	double pr = p.norm();
	double prxy = sqrt(p[0] * p[0] + p[1] * p[1]);
	double sin_phi = p[2] / pr;
	double cos_phi = prxy / pr;
	double cos_theta = p[0] / prxy;
	double sin_theta = p[1] / prxy;
	if (prxy < 1e-8)
	{
		n = p;
		n.normalize();
		p = Mesh::Point(1, 0, 0);
		return;
	}

	double dr = (1 - lambda)*a* pow(cos_phi, a - 1)*sin_phi;
	n = p + dr * OpenMesh::Vec3d(cos_theta * sin_phi, sin_phi * sin_theta, -cos_phi);
	n.normalize();
	p = OpenMesh::Vec3d(-n[1], n[0], 0);
	p.normalize();

	return;
}



