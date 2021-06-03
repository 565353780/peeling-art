#include "ARAPParametricManger.h"

#include <fstream>
#include <Eigen/SparseCore>
#include <Eigen/PardisoSupport>

using namespace std;
using namespace Eigen;

extern OpenMesh::FPropHandleT<bool> is_singular;
extern OpenMesh::FPropHandleT<int> face_type;
extern OpenMesh::VPropHandleT<Eigen::Vector3f> vert_ref;
extern OpenMesh::EPropHandleT<int> edge_type;
extern OpenMesh::VPropHandleT<int> vert_type;
extern OpenMesh::VPropHandleT<int> ori_vert_idx;

ARAPParametricManger::ARAPParametricManger()
{
	ref_fix_id = 1;
}

ARAPParametricManger::~ARAPParametricManger()
{

}

int  ARAPParametricManger::cutmeshBoundary_Inline()
{
	mesh_ = cut_mesh_ref;

	compute_cut_Boundary_id();	
	return all_boundary_id.size();
}

void ARAPParametricManger::cutmeshBoundary_fix(Mesh& ref_mesh, int times, std::pair < int, int >& before_end_id)
{
	mesh_ = cut_mesh_ref;
	fixed_id.clear();
	for (int i = 0; i < all_boundary_id[times].size(); i++)
	{
		fixed_id.insert(all_boundary_id[times][i]);
	}
	before_end_id.first = all_boundary_id[times][0];
	before_end_id.second = all_boundary_id[times][fixed_id.size() - 1];
	change_cut_ref(ref_mesh);
}


void ARAPParametricManger::compute_cut_Boundary_id()
{
	Mesh::HalfedgeHandle hfe;
	Mesh::HalfedgeHandle hfe_begin;
	for (Mesh::HalfedgeIter he_it = mesh_.halfedges_begin(); he_it != mesh_.halfedges_end(); he_it++)
	{
		if (mesh_.is_boundary(*he_it) && mesh_.property(edge_type, mesh_.edge_handle(*he_it)) == 1)
		{
			hfe = he_it.handle();
			hfe_begin = hfe;
			break;
		}
	}
	all_boundary_id.clear();

	do 
	{
		std::vector<int> one_part_bdy;
		one_part_bdy.clear();
		do
		{
			Mesh::VertexHandle vh = mesh_.to_vertex_handle(hfe);
			if (mesh_.property(vert_type, vh) == 0)
			{
				one_part_bdy.push_back(vh.idx());
			}			
			hfe = mesh_.next_halfedge_handle(hfe);
		} while (mesh_.property(edge_type, mesh_.edge_handle(hfe)) != 1 && hfe != hfe_begin);

		if (one_part_bdy.size() >= 3)
		{
			std::cout << one_part_bdy.size() << std::endl;

			//if one part is too long, may be it need be split
			if (one_part_bdy.size() > 170)
			{
				std::vector<int> part_b1, part_b2, part_b3;
				part_b1.clear();
				for (int i = 0; i < 40; i++)
				{
					part_b1.push_back(one_part_bdy[i]);
				}
				part_b2.clear();
				for (int i = part_b1.size(); i < one_part_bdy.size(); i++)
				{
					part_b2.push_back(one_part_bdy[i]);
				}
				all_boundary_id.push_back(part_b1);
				all_boundary_id.push_back(part_b2);
			}
			else
			{
				all_boundary_id.push_back(one_part_bdy);
			}
			
		}	
	} while (hfe != hfe_begin);

}

void ARAPParametricManger::change_cut_ref(Mesh& ref_mesh)
{
	std::vector<Vector3f> p_list;
	std::vector<Vector3f> q_list;
	p_list.clear();
	q_list.clear();
	for (set<int>::iterator it = fixed_id.begin(); it != fixed_id.end(); it++)
	{	
		Mesh::VertexHandle vh = cut_mesh_ref.vertex_handle(*it);
		for (Mesh::VertexIter vi_it = ref_mesh.vertices_begin(); vi_it != ref_mesh.vertices_end(); vi_it++)
		{
			if (cut_mesh_ref.property(ori_vert_idx, vh) == ref_mesh.property(ori_vert_idx, *vi_it))
			{
				p_list.push_back(cut_mesh_ref.property(vert_ref, vh));
				q_list.push_back(ref_mesh.property(vert_ref, *vi_it));
			}
		}
	}
	assert(p_list.size() == fixed_id.size());
	
	Vector3f mu_p = Vector3f::Zero();
	Vector3f mu_q = Vector3f::Zero();
	for (int i = 0; i < p_list.size(); i++)
	{
		mu_p += p_list[i];
		mu_q += q_list[i];
	}
	mu_p /= p_list.size();
	mu_q /= q_list.size();

	Matrix3f H = Matrix3f::Zero();
	for (int i = 0; i < p_list.size(); i++)
	{
		H += (p_list[i] - mu_p) * (q_list[i] - mu_q).transpose();
	}
	Matrix2f HH = H.topLeftCorner(2,2);
	JacobiSVD<Eigen::Matrix2f> svd(HH, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Matrix2f R =  svd.matrixV() * (svd.matrixU().transpose());

	Matrix3f Rr = Matrix3f::Identity();
	Rr.topLeftCorner(2, 2) = R;

	Vector3f t = -Rr*mu_p + mu_q;

	bool is_over = false;
	for (Mesh::VertexIter v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); v_it++)
	{
		int v_id = v_it.handle().idx();
		Vector3f old_p = cut_mesh_ref.property(vert_ref, cut_mesh_ref.vertex_handle(v_id));
		Vector3f new_p = Rr*old_p + t;
		mesh_.property(vert_ref, *v_it) = new_p;
	}

	int i = 0;
	for (set<int>::iterator it = fixed_id.begin(); it != fixed_id.end(); it++)
	{
		Mesh::VertexHandle vh = mesh_.vertex_handle(*it);
		mesh_.property(vert_ref, vh) = q_list[i];
		mesh_.property(vert_type, vh) = 1;
		i++;
	}
}