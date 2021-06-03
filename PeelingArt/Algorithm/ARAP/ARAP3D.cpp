#define EIGEN_VECTORIZE_SSE4_2
#define NO_USE_MESH

#include "ARAP3D.h"

#include "Mesh/surface_manager.h"
#include <fstream>
#include <Eigen/SparseCore>
#include <Eigen/PardisoSupport>
#include <float.h>


using namespace std;
using namespace Eigen;

extern OpenMesh::FPropHandleT<int> face_type;
extern OpenMesh::VPropHandleT<Eigen::Vector3f> user_vert_ref;
extern OpenMesh::VPropHandleT<int> vert_type;

ARAP3D::ARAP3D()
{
	shrink_weight = 1;
	special_weight = 1;
}

ARAP3D::~ARAP3D()
{
    delete mesh_surface_;
}

void ARAP3D::evaluateCoff()
{
	int he_num = mesh_.n_halfedges();
	int face_num = mesh_.n_faces();
	const double cot60deg = sqrt(3.0) / 3;
	const double sin60deg = sqrt(3.0) / 2;

	cot_of_hedge_.resize(he_num);
	x_of_hedge_new.resize(he_num);

	vector<int>			ind(3);
	vector<Vector3f>	vert(3);
	for (Mesh::FaceIter f_it = mesh_.faces_begin(); f_it != mesh_.faces_end(); f_it++)
	{
		int f_id = f_it.handle().idx();

		Mesh::FaceHalfedgeIter fh_it = mesh_.fh_begin(*f_it); ind[0] = (*fh_it).idx();
		++fh_it; ind[1] = (*fh_it).idx();
		++fh_it; ind[2] = (*fh_it).idx();

		if (is_singular_face_[f_id])
		{
			cot_of_hedge_[ind[0]] = cot60deg; cot_of_hedge_[ind[1]] = cot60deg; cot_of_hedge_[ind[2]] = cot60deg;

			x_of_hedge_new[ind[0]] = OpenMesh::Vec2d(0, 0);
			x_of_hedge_new[ind[1]] = OpenMesh::Vec2d(1.0, 0);
			x_of_hedge_new[ind[2]] = OpenMesh::Vec2d(0.5, sin60deg);
		}
		else
		{
			Mesh::HalfedgeHandle edge = mesh_.halfedge_handle(*f_it);
			for (int i(0); i < 3; i++)
			{
				Mesh::HalfedgeHandle eh = mesh_.next_halfedge_handle(edge);
				Mesh::VertexHandle vb = mesh_.to_vertex_handle(eh);
				vert[i] = ref_map_[vb.idx()];
				edge = mesh_.next_halfedge_handle(edge);
			}

			cot_of_hedge_[ind[0]] = cotan(vert[1] - vert[0], vert[2] - vert[0]);
			cot_of_hedge_[ind[1]] = cotan(vert[2] - vert[1], vert[0] - vert[1]);
			cot_of_hedge_[ind[2]] = cotan(vert[0] - vert[2], vert[1] - vert[2]);

			Vector3f v10 = vert[1] - vert[0];
			Vector3f v20 = vert[2] - vert[0];

			x_of_hedge_new[ind[0]] = OpenMesh::Vec2d(0, 0);
			x_of_hedge_new[ind[1]] = OpenMesh::Vec2d((vert[1] - vert[0]).norm(), 0);
			x_of_hedge_new[ind[2]] = OpenMesh::Vec2d(v10.dot(v20), (v10.cross(v20)).norm()) / v10.norm();
		}
	}
}

void ARAP3D::evaluateFLR()
{
	vector<int>			he_ind(3);
	vector<int>			vert_ind(3);
	vector<OpenMesh::Vec3d>	u_of_hedge_3d(3);

	int face_num = mesh_.n_faces();
	int vert_num = mesh_.n_vertices();
	int he_num = mesh_.n_halfedges();

	matrix_R_new.resize(he_num);
	matrix_L_new.resize(face_num);

	u_of_hedge_new.resize(he_num);
	vert_frame.resize(vert_num);

	int n_v = mesh_.n_vertices();
	for (int v_id = 0; v_id < n_v; v_id++)
	{
		Mesh::Point p = mesh_.point(mesh_.vertex_handle(v_id));
		mesh_surface_->frameAtPoint(p, vert_frame[v_id]);
	}

	for (int f_id = 0; f_id < face_num; f_id++)
	{
		Mesh::HalfedgeHandle edge = mesh_.halfedge_handle(mesh_.face_handle(f_id));

		for (int i(0); i < 3; i++)
		{
			Mesh::HalfedgeHandle eh = mesh_.next_halfedge_handle(edge);
			Mesh::VertexHandle vb = mesh_.to_vertex_handle(eh);

			he_ind[i] = edge.idx();
			vert_ind[i] = vb.idx();
			u_of_hedge_3d[i] = mesh_.point(vb);

			edge = mesh_.next_halfedge_handle(edge);
		}

		OpenMesh::Vec3d n = (u_of_hedge_3d[0] + u_of_hedge_3d[1] + u_of_hedge_3d[2]) / 3;
		n.normalize();
		OpenMesh::Vec3d e0(-n[1],n[0],0);
		if (n[0] * n[0] + n[1] * n[1] < 1e-6)
		{
			e0 = OpenMesh::Vec3d(1, 0, 0);
		}
		e0.normalize();
		OpenMesh::Vec3d e1 = OpenMesh::cross(n, e0);

//		VFrame face_frame;
//		face_frame.e0 = e0; face_frame.e1 = e1; face_frame.n = n;

		u_of_hedge_new[he_ind[0]] = OpenMesh::Vec2d(0, 0);
		u_of_hedge_new[he_ind[1]] = OpenMesh::Vec2d(OpenMesh::dot(u_of_hedge_3d[1] - u_of_hedge_3d[0], e0), OpenMesh::dot(u_of_hedge_3d[1] - u_of_hedge_3d[0], e1));
		u_of_hedge_new[he_ind[2]] = OpenMesh::Vec2d(OpenMesh::dot(u_of_hedge_3d[2] - u_of_hedge_3d[0], e0), OpenMesh::dot(u_of_hedge_3d[2] - u_of_hedge_3d[0], e1));


		double a = 0, b = 0, c = 0, d = 0;
		for (int i(0); i < 3; i++)
		{
			VFrame& vf = vert_frame[vert_ind[i]];
			a = OpenMesh::dot(e0, vf.e0); b = OpenMesh::dot(e0, vf.e1);
			c = OpenMesh::dot(e1, vf.e0); d = OpenMesh::dot(e1, vf.e1);

			matrix_R_new[he_ind[i]] = OpenMesh::Vec4d(a, b, c, d);
		}

		int i1, i2;
		a = 0, b = 0, c = 0, d = 0;
		for (int i(0); i < 3; i++)
		{
			i1 = (i + 1) % 3;
			i2 = (i + 2) % 3;
			OpenMesh::Vec2d v1 = u_of_hedge_new[he_ind[i1]] - u_of_hedge_new[he_ind[i2]];
			OpenMesh::Vec2d v2 = x_of_hedge_new[he_ind[i1]] - x_of_hedge_new[he_ind[i2]];
			a += cot_of_hedge_[he_ind[i]] * v1[0]*v2[0]; b += cot_of_hedge_[he_ind[i]] * v1[0]*v2[1]; 
			c += cot_of_hedge_[he_ind[i]] * v1[1]*v2[0]; d += cot_of_hedge_[he_ind[i]] * v1[1]*v2[1];
		}
		Matrix2d mat_temp;
		mat_temp << a, b, c, d;

		JacobiSVD<Eigen::Matrix2d> svd(mat_temp, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Matrix2d U = svd.matrixU();
		Matrix2d VT = svd.matrixV().transpose();

		if (is_singular_face_[f_id])
		{
			double s = svd.singularValues()(0);
			a = s*U(0, 0) * VT(0, 0); b = s*U(0, 0)*VT(0, 1);
			c = s*U(1, 0) * VT(0, 0); d = s*U(1, 0)*VT(0, 1);
			matrix_L_new[f_id] = OpenMesh::Vec4d(a, b, c, d);
		}
		else
		{
			Matrix2d L = U * VT;
			matrix_L_new[f_id] = OpenMesh::Vec4d(L(0,0), L(0,1), L(1,0), L(1,1));
		}
	}

}

/*
void ARAP3D::computeEnergy(double area_lambda)
{
	Matrix3d			face_frame;
	vector<int>			he_ind(3);
	vector<int>			vert_ind(3);
	vector<Vector3d>	u_of_hedge_3d(3);
	vector<Vector3f>	vert(3);

	int face_num = mesh_.n_faces();
	int he_num = mesh_.n_halfedges();
	u_of_hedge_.resize(he_num);
	energy.resize(face_num);

	double all_area = 0;

	for (Mesh::FaceIter f_it = mesh_.faces_begin(); f_it != mesh_.faces_end(); f_it++)
	{
		int f_id = f_it.handle().idx();
		Mesh::FaceHandle face = f_it.handle();
		Mesh::HalfedgeHandle edge = mesh_.halfedge_handle(face);

		for (int i(0); i < 3; i++)
		{
			Mesh::HalfedgeHandle eh = mesh_.next_halfedge_handle(edge);
			Mesh::VertexHandle vb = mesh_.to_vertex_handle(eh);

			he_ind[i] = edge.idx();

			Mesh::Point p = mesh_.point(vb);
			u_of_hedge_3d[i] = Vector3d(p[0], p[1], p[2]);

			vert[i] = area_lambda*mesh_.property(user_vert_ref, vb);
			vert_ind[i] = vb.idx();

			edge = mesh_.next_halfedge_handle(edge);
		}

		face_frame.col(0) = (u_of_hedge_3d[1] - u_of_hedge_3d[0]);
		float base_length = face_frame.col(0).norm();
		u_of_hedge_[he_ind[0]] = Vector2d(0, 0);
		u_of_hedge_[he_ind[1]] = Vector2d(base_length, 0);
		face_frame.col(1) = (u_of_hedge_3d[2] - u_of_hedge_3d[0]);
		face_frame.col(2) = face_frame.col(0).cross(face_frame.col(1));
		u_of_hedge_[he_ind[2]] = Vector2d(face_frame.col(0).dot(face_frame.col(1)), face_frame.col(2).norm()) / base_length;

		x_of_hedge_[he_ind[0]] = Vector2d(0, 0);
		x_of_hedge_[he_ind[1]] = Vector2d((vert[1] - vert[0]).norm(), 0);
		Vector3f v10 = vert[1] - vert[0];
		Vector3f v20 = vert[2] - vert[0];
		x_of_hedge_[he_ind[2]] = Vector2d(v10.dot(v20), (v10.cross(v20)).norm()) / v10.norm();
		double area = (v10.cross(v20).norm() / 2);

		Matrix2d A,R;
		A << (x_of_hedge_[he_ind[1]] - x_of_hedge_[he_ind[0]])(0), (x_of_hedge_[he_ind[1]] - x_of_hedge_[he_ind[0]])(1),
			(x_of_hedge_[he_ind[2]] - x_of_hedge_[he_ind[0]])(0), (x_of_hedge_[he_ind[2]] - x_of_hedge_[he_ind[0]])(1);
		Vector2d b1, b2, r1, r2;
		b1 << (u_of_hedge_[he_ind[1]] - u_of_hedge_[he_ind[0]])(0), (u_of_hedge_[he_ind[2]] - u_of_hedge_[he_ind[0]])(0);
		b2 << (u_of_hedge_[he_ind[1]] - u_of_hedge_[he_ind[0]])(1), (u_of_hedge_[he_ind[2]] - u_of_hedge_[he_ind[0]])(1);
		r1 = A.lu().solve(b1);
		r2 = A.lu().solve(b2);
		R << r1(0), r1(1), r2(0), r2(1);

		JacobiSVD<Eigen::Matrix2d> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);

		if (is_singular_face_[f_id])
		{
			energy[f_id] = 0;
		}
		else
		{
			Vector2d sigv = svd.singularValues();
			energy[f_id] = myfun(sigv(0)) + myfun(sigv(1));
			energy[f_id] *= area;
			all_area += area;
		}
	}

	double a = 0;
	for (int i = 0; i < energy.size(); i++)
	{
		a = a + energy[i];
	}
	a = a * 11 / all_area;

	std::cout << "energy is =======================================================================" << std::endl;
	std::cout << a << "face is " << energy.size()  << std::endl;
	std::cout << "=================================================================================";
}*/

bool ARAP3D::createMatrix()
{
	evaluateFLR();

	int					j, t;

	Eigen::SparseMatrix<double>	matrix_A(2 * mesh_.n_vertices(), 2 * mesh_.n_vertices());
	VectorXd vector_b(2 * mesh_.n_vertices());
	vector_b.setZero();
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;

	for (Mesh::VertexIter v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); v_it++)
	{
		int i = v_it.handle().idx();
		Mesh::VertexHandle vert = v_it.handle();
		Mesh::HalfedgeHandle edge = mesh_.halfedge_handle(vert);
		Mesh::HalfedgeHandle p_edge;

		/*处理边界点的起始边*/
		bool on_boundary = mesh_.is_boundary(vert);
		if (on_boundary)
		{
			while (mesh_.face_handle(edge).is_valid())
			{
				edge = mesh_.prev_halfedge_handle(edge);
				edge = mesh_.opposite_halfedge_handle(edge);
			}
		}

		OpenMesh::Vec4d wgt_i(0,0,0,0);
		do
		{
			OpenMesh::Vec4d wgt_j(0, 0, 0, 0);
			Mesh::VertexHandle vb = mesh_.to_vertex_handle(edge);			
			j = vb.idx();

			int ind_i, ind_j, ind_ij;
			if (mesh_.face_handle(edge).is_valid())
			{
				Mesh::FaceHandle ef = mesh_.face_handle(edge);
				t = ef.idx();
				ind_j = mesh_.prev_halfedge_handle(edge).idx();
				ind_i = mesh_.next_halfedge_handle(edge).idx();
				ind_ij = edge.idx();

				double sigma = 1;
				if (is_singular_face_[t])
				{
					sigma = shrink_weight;
				}
				else
				{
					sigma = 1;
					if (mesh_.property(face_type, mesh_.face_handle(edge)) == 2)
					{
						sigma = special_weight;
					}
				}

				OpenMesh::Vec4d matrix_RT_new = Matrix_transpose(matrix_R_new[ind_i]);
				wgt_j -= sigma*cot_of_hedge_[ind_ij] * Matrix_product(matrix_RT_new, matrix_R_new[ind_j]);
				wgt_i += sigma*cot_of_hedge_[ind_ij] * Matrix_product(matrix_RT_new, matrix_R_new[ind_i]);
				
				OpenMesh::Vec2d nx = Matrix_product_vector(matrix_L_new[t], x_of_hedge_new[ind_i] - x_of_hedge_new[ind_j]);
				OpenMesh::Vec2d nu = u_of_hedge_new[ind_i] - u_of_hedge_new[ind_j];
				OpenMesh::Vec2d nb = sigma*cot_of_hedge_[ind_ij] * Matrix_product_vector(matrix_RT_new, nx - nu);
				vector_b.segment<2>(2 * i) += Eigen::Vector2d(nb[0], nb[1]);
			}

			p_edge = mesh_.opposite_halfedge_handle(edge);
			if (mesh_.face_handle(p_edge).is_valid())
			{
				Mesh::FaceHandle pef = mesh_.face_handle(p_edge);
				t = pef.idx();
				ind_j = mesh_.next_halfedge_handle(p_edge).idx();
				ind_i = mesh_.prev_halfedge_handle(p_edge).idx();
				ind_ij = p_edge.idx();

				double sigma = 1;
				if (is_singular_face_[t])
				{
					sigma = shrink_weight;
				}
				else
				{
					sigma = 1;
					if (mesh_.property(face_type, mesh_.face_handle(p_edge)) == 2)
					{
						sigma = special_weight;
					}
				}

				OpenMesh::Vec4d matrix_RT_new = Matrix_transpose(matrix_R_new[ind_i]);
				wgt_j -= sigma*cot_of_hedge_[ind_ij] * Matrix_product(matrix_RT_new, matrix_R_new[ind_j]);
				wgt_i += sigma*cot_of_hedge_[ind_ij] * Matrix_product(matrix_RT_new, matrix_R_new[ind_i]);

				OpenMesh::Vec2d nx = Matrix_product_vector(matrix_L_new[t], x_of_hedge_new[ind_i] - x_of_hedge_new[ind_j]);
				OpenMesh::Vec2d nu = u_of_hedge_new[ind_i] - u_of_hedge_new[ind_j];
				OpenMesh::Vec2d nb = sigma*cot_of_hedge_[ind_ij] * Matrix_product_vector(matrix_RT_new, nx - nu);
				vector_b.segment<2>(2 * i) += Eigen::Vector2d(nb[0], nb[1]);
			}
			tripletList.push_back(T(2 * i, 2 * j, wgt_j[0]));
			tripletList.push_back(T(2 * i, 2 * j + 1, wgt_j[1]));
			tripletList.push_back(T(2 * i + 1, 2 * j, wgt_j[2]));
			tripletList.push_back(T(2 * i + 1, 2 * j + 1, wgt_j[3]));

			edge = mesh_.opposite_halfedge_handle(edge);
			edge = mesh_.next_halfedge_handle(edge);

		} while (edge.is_valid() &&  edge != mesh_.halfedge_handle(vert));

		tripletList.push_back(T(2 * i, 2 * i, wgt_i[0]));
		tripletList.push_back(T(2 * i, 2 * i + 1, wgt_i[1]));
		tripletList.push_back(T(2 * i + 1, 2 * i, wgt_i[2]));
		tripletList.push_back(T(2 * i + 1, 2 * i + 1, wgt_i[3]));
	}

    matrix_A.setFromTriplets(tripletList.begin(), tripletList.end());

	PardisoLU<SparseMatrix<double>>	solver;

    solver.compute(matrix_A);

	if (solver.info() != Success)
	{
		cout << "BiCGSTAB: Decomposition failed." << endl;
		return false;
	}

	vector_du_.resize(vector_b.size());
	vector_du_ = solver.solve(vector_b);

	return true;
}

void ARAP3D::optimizeMesh()
{
	clock_t t1;
	t1 = clock();
    //std::cout << "begin_sharp_optimize " << t1 << std::endl;
	evaluateCoff();
	createMatrix();
	
	float infnorm = max(-(vector_du_.minCoeff()), vector_du_.maxCoeff());

	vector_du_ *= 0.1f / max(infnorm, 0.1f);

	int vert_num = mesh_.n_vertices();

    for (Mesh::VertexIter v_it = mesh_.vertices_sbegin(); v_it != mesh_.vertices_end(); ++v_it)
	{
		int id = v_it->idx();
		Mesh::Point p = mesh_.point(v_it);
		Mesh::Point new_p = mesh_surface_->PointWithCoordinate(p, vert_frame[id], vector_du_.segment<2>(2 * id));
		mesh_.set_point(v_it.handle(), new_p);
	}

	energy.resize(mesh_.n_faces());

	t1 = clock();
    //std::cout << "end_sharp_optimize1 " << t1 << std::endl;
}
