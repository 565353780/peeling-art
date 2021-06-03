#include "ARAPPara_3D.h"

using namespace Eigen;
using namespace std;

ARAPPara_3D::ARAPPara_3D()
{
}


ARAPPara_3D::~ARAPPara_3D()
{
	matrix_L_.clear();
}

void ARAPPara_3D::evaluateCoff()
{
	int edge_num = mesh_surface.n_edges();
	int face_num = mesh_surface.n_faces();

	cot_of_edge_.resize(edge_num, 0);

	for (Mesh::FaceIter f_it = mesh_surface.faces_begin(); f_it != mesh_surface.faces_end(); f_it++)
	{
		int f_id = f_it.handle().idx();
		Mesh::FaceHandle face = f_it.handle();
		Mesh::HalfedgeHandle h_edge = mesh_surface.halfedge_handle(face);
		vector<int>			ind(3);
		vector<int>			edge_ind(3);
		vector<Vector3f>	vert(3);


		/*找到三角形的三个半边与顶点*/
		for (int i(0); i < 3; i++)
		{
			Mesh::HalfedgeHandle eh = mesh_surface.next_halfedge_handle(h_edge);
			Mesh::VertexHandle vb = mesh_surface.to_vertex_handle(eh);
			ind[i] = h_edge.idx();
			edge_ind[i] = mesh_surface.edge_handle(h_edge).idx();
			Mesh::Point p = mesh_surface.point(vb);
			vert[i] = Eigen::Vector3f(p[0], p[1], p[2]);
			h_edge = mesh_surface.next_halfedge_handle(h_edge);
		}


		int i1, i2;
		/*计算每条边对角的余切值*/
		for (int i(0); i < 3; i++)
		{
			i1 = (i + 1) % 3;
			i2 = (i + 2) % 3;
			double coti = cotan(vert[i1] - vert[i], vert[i2] - vert[i]);
			cot_of_edge_[edge_ind[i]] += coti/2;
		}
	}

	matrix_L_.resize(mesh_surface.n_vertices());
	Mesh::Point m_center(0,0,0);
	for (int i = 0; i < mesh_surface.n_vertices(); i++)
	{
		matrix_L_[i].setIdentity();
	}
	
	center_id = 0;
	std::vector<int> fix_id_list(4,0);
	std::vector<double> fix_dis(4, 0);
	double center_dis = (mesh_plane.point(mesh_plane.vertex_handle(0)) - m_center).norm();
	for (int i = 1; i < mesh_surface.n_vertices(); i++)
	{
		Mesh::Point p = mesh_plane.point(mesh_plane.vertex_handle(i));
		double dis = (p - m_center).norm();
		if (dis < center_dis)
		{
			center_dis = dis;
			center_id = i;
		}

		if (p[0] < fix_dis[0])
		{
			fix_dis[0] = p[0];
			fix_id_list[0] = i;
		}
		if (p[0] > fix_dis[1])
		{
			fix_dis[1] = p[0];
			fix_id_list[1] = i;
		}
		if (p[0] < 0 && p[1] < fix_dis[2])
		{
			fix_dis[2] = p[1];
			fix_id_list[2] = i;
		}
		if (p[0] > 0 && p[1] < fix_dis[3])
		{
			fix_dis[3] = p[1];
			fix_id_list[3] = i;
		}
	}

	fix_id.clear();
	for (int i = 0; i < fix_id_list.size(); i++)
	{
		fix_id.insert(fix_id_list[i]);
		std::cout << fix_id_list[i] << std::endl;
	}

	for (Mesh::EdgeIter e_it = mesh_surface.edges_begin(); e_it != mesh_surface.edges_end(); e_it++)
	{
		if (mesh_surface.is_boundary(*e_it))
		{
			cot_of_edge_[e_it->idx()] *= 1;
		}
	}


}

void ARAPPara_3D::evaluateL()
{
	Matrix3f			mat_temp;
	vector<int>			he_ind(3);
	vector<int>			vert_ind(3);
	vector<double>      begin_edge;
	vector<double>      end_edge;

	for (Mesh::VertexIter v_it = mesh_surface.vertices_begin(); v_it != mesh_surface.vertices_end(); v_it++)
	{
		int vi_id = v_it->idx();
		Matrix3d S = Matrix3d::Zero();

		double error1 = 0;
		begin_edge.clear();
		for (Mesh::VertexOHalfedgeIter voh_it = mesh_surface.voh_begin(*v_it); voh_it != mesh_surface.voh_end(*v_it); voh_it++)
		{
			Mesh::VertexHandle vj = mesh_surface.to_vertex_handle(voh_it);
			int vj_id = vj.idx();
			int eij_id = mesh_surface.edge_handle(voh_it).idx();

			Mesh::Point pi1 = mesh_surface.point(mesh_surface.vertex_handle(vi_id));
			Mesh::Point pj1 = mesh_surface.point(mesh_surface.vertex_handle(vj_id));
			Mesh::Point pi2 = mesh_plane.point(mesh_plane.vertex_handle(vi_id));
			Mesh::Point pj2 = mesh_plane.point(mesh_plane.vertex_handle(vj_id));

			Vector3d e_ij_S = Vector3d(pi1[0] - pj1[0], pi1[1] - pj1[1], pi1[2] - pj1[2]);
			Vector3d e_ij_P = Vector3d(pi2[0] - pj2[0], pi2[1] - pj2[1], pi2[2] - pj2[2]);
			double a = cot_of_edge_[eij_id];
			S += cot_of_edge_[eij_id]*e_ij_S*(e_ij_P.transpose());

			Vector3d er_ij = e_ij_P - matrix_L_[vi_id] * e_ij_S;
			begin_edge.push_back(er_ij.norm());
			error1 += cot_of_edge_[eij_id] * er_ij.norm();
		}

		/*计算Lt*/
		JacobiSVD<Eigen::Matrix3d> svd(S, Eigen::ComputeFullU | Eigen::ComputeFullV);
		matrix_L_[vi_id] = svd.matrixV() * (svd.matrixU().transpose());
	}
}

double ARAPPara_3D::solveMatrix_xy()
{
	double				wgt_j, wgt_i;

	int					vert_num = mesh_surface.n_vertices();
	VectorXd vector_x(vert_num);
	VectorXd vector_y(vert_num);
	vector_x.fill(0.0f);
	vector_y.fill(0.0f);

	for (Mesh::VertexIter v_it = mesh_surface.vertices_begin(); v_it != mesh_surface.vertices_end(); v_it++)
	{
		int vi_id = v_it.handle().idx();
		if (vi_id == center_id)
		{
			vector_x(vi_id) = mesh_plane.point(mesh_plane.vertex_handle(vi_id))[0];
			vector_y(vi_id) = mesh_plane.point(mesh_plane.vertex_handle(vi_id))[1];
			continue;
		}

		wgt_i = 0;
		double lambda = 0;
		for (Mesh::VertexOHalfedgeIter voh_it = mesh_surface.voh_begin(*v_it); voh_it != mesh_surface.voh_end(*v_it); voh_it++)
		{
			Mesh::VertexHandle vj = mesh_surface.to_vertex_handle(voh_it);
			int vj_id = vj.idx();
			int eij_id = mesh_surface.edge_handle(voh_it).idx();

			Mesh::Point pi1 = mesh_surface.point(mesh_surface.vertex_handle(vi_id));
			Mesh::Point pj1 = mesh_surface.point(mesh_surface.vertex_handle(vj_id));

			wgt_i += cot_of_edge_[eij_id];

			Vector3d e_ij_S = Vector3d(pi1[0] - pj1[0], pi1[1] - pj1[1], pi1[2] - pj1[2]);
			vector_x(vi_id) += 0.5*cot_of_edge_[eij_id] * ((matrix_L_[vi_id] + matrix_L_[vj_id])*e_ij_S).x();
			vector_y(vi_id) += 0.5*cot_of_edge_[eij_id] * ((matrix_L_[vi_id] + matrix_L_[vj_id])*e_ij_S).y();
		}
		vector_x(vi_id) += lambda * mesh_plane.point(mesh_plane.vertex_handle(vi_id))[0];
		vector_y(vi_id) += lambda * mesh_plane.point(mesh_plane.vertex_handle(vi_id))[1];
	}

	VectorXd result_x(vert_num);
	VectorXd result_y(vert_num);

	result_x = solver_xy.solve(vector_x);
	result_y = solver_xy.solve(vector_y);

	for (Mesh::VertexIter v_it = mesh_plane.vertices_sbegin(); v_it != mesh_plane.vertices_end(); ++v_it)
	{
		int id = v_it->idx();
		Mesh::Point p = mesh_plane.point(*v_it);
		p[0] = result_x(id);
		p[1] = result_y(id);
		mesh_plane.set_point(*v_it, p);
	}

	double error = 0;
/*
	for (Mesh::VertexIter v_it = mesh_surface.vertices_begin(); v_it != mesh_surface.vertices_end(); v_it++)
	{
		int vi_id = v_it->idx();
		Matrix3d S;
		for (Mesh::VertexOHalfedgeIter voh_it = mesh_surface.voh_begin(*v_it); voh_it != mesh_surface.voh_end(*v_it); voh_it++)
		{
			Mesh::VertexHandle vj = mesh_surface.to_vertex_handle(voh_it);
			int vj_id = vj.idx();
			int eij_id = mesh_surface.edge_handle(voh_it).idx();

			Mesh::Point pi1 = mesh_surface.point(mesh_surface.vertex_handle(vi_id));
			Mesh::Point pj1 = mesh_surface.point(mesh_surface.vertex_handle(vj_id));
			Mesh::Point pi2 = mesh_plane.point(mesh_plane.vertex_handle(vi_id));
			Mesh::Point pj2 = mesh_plane.point(mesh_plane.vertex_handle(vj_id));

			Vector3d e_ij_S = Vector3d(pi1[0] - pj1[0], pi1[1] - pj1[1], pi1[2] - pj1[2]);
			Vector3d e_ij_P = Vector3d(pi2[0] - pj2[0], pi2[1] - pj2[1], pi2[2] - pj2[2]);

			Vector3d er_ij = e_ij_P - matrix_L_[vi_id] * e_ij_S;
			error += cot_of_edge_[eij_id] * er_ij.squaredNorm();
		}
	}*/

	//std::cout << "the error after globalxy is =======" << error << "========="  << std::endl;

	return error;
}

double ARAPPara_3D::solveMatrixz()
{
	double				wgt_j, wgt_i;

	int					vert_num = mesh_surface.n_vertices();

	VectorXd vector_z(vert_num);
	vector_z.fill(0.0f);

	for (Mesh::VertexIter v_it = mesh_surface.vertices_begin(); v_it != mesh_surface.vertices_end(); v_it++)
	{
		int vi_id = v_it.handle().idx();
		if (fix_id.find(vi_id) != fix_id.end())
		{
			vector_z(vi_id) = 0;
		}

		wgt_i = 0;
		for (Mesh::VertexOHalfedgeIter voh_it = mesh_surface.voh_begin(*v_it); voh_it != mesh_surface.voh_end(*v_it); voh_it++)
		{
			Mesh::VertexHandle vj = mesh_surface.to_vertex_handle(voh_it);
			int vj_id = vj.idx();
			int eij_id = mesh_surface.edge_handle(voh_it).idx();

			Mesh::Point pi1 = mesh_surface.point(mesh_surface.vertex_handle(vi_id));
			Mesh::Point pj1 = mesh_surface.point(mesh_surface.vertex_handle(vj_id));

			wgt_i += cot_of_edge_[eij_id];

			Vector3d e_ij_S = Vector3d(pi1[0] - pj1[0], pi1[1] - pj1[1], pi1[2] - pj1[2]);
			vector_z(vi_id) += 0.5*cot_of_edge_[eij_id]*((matrix_L_[vi_id] + matrix_L_[vj_id])*e_ij_S).z();
		}
	}
	VectorXd result_z(vert_num);
	result_z = solver_z.solve(vector_z);
	
	for (Mesh::VertexIter v_it = mesh_plane.vertices_sbegin(); v_it != mesh_plane.vertices_end(); ++v_it)
	{
		int id = v_it->idx();
		Mesh::Point p = mesh_plane.point(*v_it);

		p[2] = result_z(id);
		mesh_plane.set_point(*v_it, p);
	}

	double error = 0;
	for (Mesh::VertexIter v_it = mesh_surface.vertices_begin(); v_it != mesh_surface.vertices_end(); v_it++)
	{
		int vi_id = v_it->idx();
		Matrix3d S;
		for (Mesh::VertexOHalfedgeIter voh_it = mesh_surface.voh_begin(*v_it); voh_it != mesh_surface.voh_end(*v_it); voh_it++)
		{
			Mesh::VertexHandle vj = mesh_surface.to_vertex_handle(voh_it);
			int vj_id = vj.idx();
			int eij_id = mesh_surface.edge_handle(voh_it).idx();

			Mesh::Point pi1 = mesh_surface.point(mesh_surface.vertex_handle(vi_id));
			Mesh::Point pj1 = mesh_surface.point(mesh_surface.vertex_handle(vj_id));
			Mesh::Point pi2 = mesh_plane.point(mesh_plane.vertex_handle(vi_id));
			Mesh::Point pj2 = mesh_plane.point(mesh_plane.vertex_handle(vj_id));

			Vector3d e_ij_S = Vector3d(pi1[0] - pj1[0], pi1[1] - pj1[1], pi1[2] - pj1[2]);
			Vector3d e_ij_P = Vector3d(pi2[0] - pj2[0], pi2[1] - pj2[1], pi2[2] - pj2[2]);

			Vector3d er_ij = e_ij_P - matrix_L_[vi_id] * e_ij_S;
			error += cot_of_edge_[eij_id] * er_ij.squaredNorm();
		}
	}

	std::cout << "the error after globalz is =======" << error << "=========" << std::endl;

	return error;
}

void ARAPPara_3D::CreatAllMatrix()
{
	double				wgt_j, wgt_i;
	int					vert_num = mesh_surface.n_vertices();
	matrix_Axy = Eigen::SparseMatrix<double>(vert_num, vert_num);
	matrix_Az = Eigen::SparseMatrix<double>(vert_num, vert_num);

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList_xy;
	tripletList_xy.clear();
	for (Mesh::VertexIter v_it = mesh_surface.vertices_begin(); v_it != mesh_surface.vertices_end(); v_it++)
	{
		int vi_id = v_it.handle().idx();
		if (vi_id == center_id)
		{
			tripletList_xy.push_back(T(vi_id, vi_id, 1.0));
			continue;
		}
		wgt_i = 0;
		double lambda = 0;
		for (Mesh::VertexOHalfedgeIter voh_it = mesh_surface.voh_begin(*v_it); voh_it != mesh_surface.voh_end(*v_it); voh_it++)
		{
			Mesh::VertexHandle vj = mesh_surface.to_vertex_handle(voh_it);
			int vj_id = vj.idx();
			int eij_id = mesh_surface.edge_handle(voh_it).idx();

			Mesh::Point pi1 = mesh_surface.point(mesh_surface.vertex_handle(vi_id));
			Mesh::Point pj1 = mesh_surface.point(mesh_surface.vertex_handle(vj_id));

			wgt_i += cot_of_edge_[eij_id];
			tripletList_xy.push_back(T(vi_id, vj_id, -cot_of_edge_[eij_id]));

		}
		tripletList_xy.push_back(T(vi_id, vi_id, wgt_i + lambda));
	}

	matrix_Axy.setFromTriplets(tripletList_xy.begin(), tripletList_xy.end());
	PardisoLU<SparseMatrix<double>>	solver;
	solver_xy.compute(matrix_Axy);
	if (solver_xy.info() != Success)
	{
		cout << "PardisoLU: Decomposition failed." << endl;
	}
	tripletList_xy.clear();

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList_z;
	tripletList_z.clear();
	for (Mesh::VertexIter v_it = mesh_surface.vertices_begin(); v_it != mesh_surface.vertices_end(); v_it++)
	{
		int vi_id = v_it.handle().idx();
		if (fix_id.find(vi_id) != fix_id.end())
		{
			tripletList_z.push_back(T(vi_id, vi_id, 1.0));
			continue;
		}

		wgt_i = 0;
		for (Mesh::VertexOHalfedgeIter voh_it = mesh_surface.voh_begin(*v_it); voh_it != mesh_surface.voh_end(*v_it); voh_it++)
		{
			Mesh::VertexHandle vj = mesh_surface.to_vertex_handle(voh_it);
			int vj_id = vj.idx();
			int eij_id = mesh_surface.edge_handle(voh_it).idx();

			Mesh::Point pi1 = mesh_surface.point(mesh_surface.vertex_handle(vi_id));
			Mesh::Point pj1 = mesh_surface.point(mesh_surface.vertex_handle(vj_id));

			wgt_i += cot_of_edge_[eij_id];
			tripletList_z.push_back(T(vi_id, vj_id, -cot_of_edge_[eij_id]));
		}
		tripletList_z.push_back(T(vi_id, vi_id, wgt_i));
	}

	matrix_Az.setFromTriplets(tripletList_z.begin(), tripletList_z.end());
	solver_z.compute(matrix_Az);
	if (solver_z.info() != Success)
	{
		cout << "BiCGSTAB: Decomposition failed." << endl;
	}
	tripletList_z.clear();
}

void ARAPPara_3D::optimizeMesh()
{
	evaluateCoff();
	CreatAllMatrix();

	double old_e = 1000;

	int i = 0;
	for (; i < 100; i++)
	{
		Position_Up();
		evaluateL();
		double error = solveMatrixz();
		if (abs(error - old_e) < 0.0001 || error<0) break;
		evaluateL();
		solveMatrix_xy();

		old_e = error;
	}
	if (i == 100)
	{
		evaluateL();
		solveMatrixz();
	}

	double maxh = 0;

	OpenMesh::IO::write_mesh(mesh_plane, "Interaction\\CurvedMesh.obj");

	std::cout << "max iter is " << i << std::endl;
	std::cout << "max high is " << maxh << std::endl;
}

void ARAPPara_3D::Position_Up()
{
	for (Mesh::VertexIter v_it = mesh_plane.vertices_begin(); v_it != mesh_plane.vertices_end(); v_it++)
	{
		Mesh::Point p = mesh_plane.point(*v_it);
		p[2] = abs(p[2]);
		mesh_plane.set_point(*v_it, p);
	}
	mesh_plane.update_normals();
}