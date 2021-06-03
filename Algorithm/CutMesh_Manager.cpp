#include "CutMesh_Manager.h"
#include "Algorithm/BPoly/BPolyline.h"
#include "Algorithm/MeshGeneration.h"
#include <Eigen/Dense>

extern OpenMesh::VPropHandleT<Eigen::Vector3f> vert_ref;

CutMesh_Manager::CutMesh_Manager()
{
}


CutMesh_Manager::~CutMesh_Manager()
{
	cut_mesh_a.clear();
	cut_mesh_b.clear();
	for (int i = 0; i < cut_mesh_b_list.size(); i++)
	{
		cut_mesh_b_list[i].clear();
	}
	cut_mesh_b_list.clear();
	is_show_list.clear();
	for (int i = 0; i < cut_b_boundary.size(); i++)
	{
		delete cut_b_boundary[i];
	}
	cut_b_boundary.clear();
}

void CutMesh_Manager::InlineBoundary()
{
	std::vector<int> cut_polyline_id;
	cut_b_boundary.resize(cut_mesh_b_list.size());
	cut_b_endpoints.resize(cut_mesh_b_list.size());

	cut_polyline_id.clear();
	Mesh& my_mesh = cut_mesh_b;
	auto heit = my_mesh.halfedges_begin();
	while (!my_mesh.is_boundary(*heit))
		heit++;
	auto he_start = *heit;
	auto he_it = he_start;
	do
	{
		he_it = my_mesh.next_halfedge_handle(he_it);
		cut_polyline_id.push_back(my_mesh.to_vertex_handle(he_it).idx());
	} while (he_it != he_start);

	for (int i = 0; i < cut_mesh_b_list.size(); i++)
	{
		my_mesh = cut_mesh_b_list[i];
		cut_b_boundary[i] = new BPolyline();
		if (is_show_list[i])
		{
			std::vector<Eigen::Vector3f> boundary(cut_polyline_id.size());
			for (int j = 0; j < cut_polyline_id.size(); j++)
			{
				boundary[j] = my_mesh.property(vert_ref, my_mesh.vertex_handle(cut_polyline_id[j]));
			}

			cut_b_boundary[i]->IntPolyline(boundary, true);

			cut_b_endpoints[i].first = my_mesh.property(vert_ref, my_mesh.vertex_handle(cut_b_endpoints_id[i].first));
			cut_b_endpoints[i].second = my_mesh.property(vert_ref, my_mesh.vertex_handle(cut_b_endpoints_id[i].second));

			for (int j = 0; j < cut_polyline_id.size(); j++)
			{
				if (cut_b_endpoints_id[i].first == cut_polyline_id[j]) cut_b_endpoints_id[i].first = j;
				if (cut_b_endpoints_id[i].second == cut_polyline_id[j]) cut_b_endpoints_id[i].second = j;
			}

			for (Mesh::VertexIter v_it = cut_mesh_b_list[i].vertices_begin(); v_it != cut_mesh_b_list[i].vertices_end(); v_it++)
			{
				Eigen::Vector3f ref = cut_mesh_b_list[i].property(vert_ref, *v_it);
				Mesh::Point p(ref[0], ref[1], ref[2]);
				cut_mesh_b_list[i].set_point(v_it.handle(), p);
			}
		}
	}

	
}

void CutMesh_Manager::getColor(int ci, double& R, double& G, double& B)
{
	double nn = 255;

	switch (ci)
	{
	case 1:
		R = 253 / nn; G = 96 / nn;  B = 95 / nn; break;
	case 2:
		R = 232 / nn; G = 132 / nn; B = 6 / nn; break;
	case 3:
		R = 177 / nn; G = 85 / nn; B = 240 / nn; break;
	case 4:
		R = 44 / nn; G = 195 / nn; B = 175 / nn; break;
	case 5:
		R = 247 / nn; G = 217 / nn; B = 3 / nn; break;
	case 6:
		R = 253 / nn; G = 171 / nn; B = 1 / nn; break;
	default:
		R = 0; G = 0; B = 0; break;
	}
}
