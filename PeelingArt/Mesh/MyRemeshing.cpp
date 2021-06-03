#include "MyRemeshing.h"
#include "Mesh/surface_manager.h"

extern OpenMesh::FPropHandleT<bool> is_singular;
extern OpenMesh::VPropHandleT<Eigen::Vector3f> vert_ref;
extern OpenMesh::EPropHandleT<int> edge_type;
extern OpenMesh::VPropHandleT<int> vert_type;


MyRemeshing::MyRemeshing()
{
}


MyRemeshing::~MyRemeshing()
{
}

void MyRemeshing::My_Remeshing(double mu, Mesh& mesh)
{
	My_Remeshing_Out_Part(mu, mesh);
	mesh.update_normals();
}

void  MyRemeshing::My_Remeshing_Out_Part(double mu, Mesh& mesh)
{
	float target_length_out = 0;
	float all_in_edge_length = 0;
	int in_edge_number = 0;
	for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
	{
		if (mesh.property(edge_type, *e_it) < 1)
		{
			float len = mesh.calc_edge_length(e_it.handle());
			all_in_edge_length = all_in_edge_length + len;
			in_edge_number++;
		}
	}
	target_length_out = mu*all_in_edge_length / in_edge_number;

	for (int split_time = 0; split_time < 3; split_time++)
	{
		int i = 0;
		int ne = mesh.n_edges();
		for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end() && i<ne; e_it++, i++)
		{
			float e_length = mesh.calc_edge_length(e_it.handle());
			if (e_length > target_length_out *4/3 && mesh.property(edge_type, *e_it) == 1)
			{
				My_edge_Split_Out(e_it.handle(), mesh);
			}
		}
	}

	my_edge_cost.resize(mesh.n_edges(), 0);
	for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
	{
		int i = e_it.handle().idx();
		if (mesh.property(edge_type, *e_it) <= 0)
		{
			my_edge_cost[i] = 1;
		}
	}
	for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
	{
		float e_length = mesh.calc_edge_length(e_it.handle());
		int e_id = e_it.handle().idx();
		if (e_length < target_length_out * 3 / 4 && my_edge_cost[e_id] < 1 && mesh.property(edge_type, *e_it) == 1)
		{
			My_edge_Collapse_Out(e_it.handle(), mesh);
		}
	}
	mesh.garbage_collection();

	mesh.update_normals();

	std::vector< float > vertice_area;
	std::vector< float > arealist;
	getVerticesArea(vertice_area, arealist, mesh);

	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		if (mesh.property(vert_type, *v_it) == 1)
		{
			Mesh::Point new_p = Mesh::Point(0, 0, 0);
			float weight = 0;
			for (Mesh::VertexVertexIter vv_it = mesh.vv_begin(v_it.handle()); vv_it != mesh.vv_end(v_it.handle()); vv_it++)
			{
				new_p = new_p + mesh.point(vv_it.handle())*vertice_area[vv_it.handle().idx()];
				weight = weight + vertice_area[vv_it.handle().idx()];
			}
			new_p = new_p / weight;

			Mesh::Point project_p;
			mesh_suface_->PointProjection(new_p, project_p);
			mesh.set_point(v_it.handle(), project_p);
		}
	}
}

bool  MyRemeshing::My_edge_Collapse_Out(Mesh::EdgeHandle& eh, Mesh& mesh)
{
	bool is_collapse = false;
	Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh, 0);
	Mesh::HalfedgeHandle heh1 = mesh.halfedge_handle(eh, 1);
	Mesh::VertexHandle ve = mesh.to_vertex_handle(heh);
	Mesh::VertexHandle vb = mesh.from_vertex_handle(heh);

	Mesh::Point new_p = (mesh.point(ve) + mesh.point(vb)) / 2;
	Mesh::Point project_p;
	mesh_suface_->PointProjection(new_p, project_p);

	if (mesh.property(vert_type, vb) == 1 && mesh.property(vert_type, ve) == 1)
	{
		if (mesh.is_collapse_ok(heh))
		{
			mesh.collapse(heh);
			mesh.set_point(ve, project_p);
			is_collapse = true;
		}
		else if (mesh.is_collapse_ok(heh1))
		{
			mesh.collapse(heh1);
			mesh.set_point(vb, project_p);
			is_collapse = true;
		}
	}
	else if (mesh.property(vert_type, vb) == 1 && mesh.property(vert_type, ve) == 0)
	{
		if (mesh.is_collapse_ok(heh))
		{
			mesh.collapse(heh);
			is_collapse = true;
		}
	}
	else if (mesh.property(vert_type, vb) == 0 && mesh.property(vert_type, ve) == 1)
	{
		if (mesh.is_collapse_ok(heh1))
		{
			mesh.collapse(heh1);
			is_collapse = true;
		}
	}

	if (is_collapse == true)
	{
		//for (Mesh::VertexEdgeIter ve_it = mesh.ve_begin(ve); ve_it != mesh.ve_end(ve); ve_it++)
		//{
		//	int veid = ve_it.handle().idx();
		//	my_edge_cost[veid] += 1;
		//}
	}
	return is_collapse;
}

bool  MyRemeshing::My_edge_Split_Out(Mesh::EdgeHandle& eh, Mesh& mesh)
{
	Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh, 0);
	Mesh::VertexHandle ve = mesh.to_vertex_handle(heh);
	Mesh::VertexHandle vb = mesh.from_vertex_handle(heh);

	Mesh::Point new_p = (mesh.point(ve) + mesh.point(vb)) / 2;
	Mesh::Point project_p;
	mesh_suface_->PointProjection(new_p, project_p);

	if (mesh.property(edge_type, eh) = 1)
	{
		Mesh::VertexHandle new_evh = mesh.split(eh, project_p);
		mesh.property(vert_type, new_evh) = 1;
		for (Mesh::VertexEdgeIter ve_it = mesh.ve_begin(new_evh); ve_it != mesh.ve_end(new_evh); ve_it++)
		{
			mesh.property(edge_type, *ve_it) = 1;
		}
		for (Mesh::VertexFaceIter vf_it = mesh.vf_begin(new_evh); vf_it != mesh.vf_end(new_evh); vf_it++)
		{
			mesh.property(is_singular, *vf_it) = true;
		}
	}
	return true;
}

void MyRemeshing::getVerticesArea(std::vector<float> &vertice_area, std::vector<float>& arealist, Mesh& mesh)
{
	arealist.resize(mesh.n_faces());
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		int fid = f_it.handle().idx();
		Mesh::FaceVertexIter fv_it = mesh.fv_begin(f_it); Mesh::Point p0 = mesh.point(*fv_it);
		++fv_it; Mesh::Point p1 = mesh.point(*fv_it);
		++fv_it; Mesh::Point p2 = mesh.point(*fv_it);
		arealist[fid] = (cross(p1 - p0, p2 - p0)).norm() / 2;
	}
	vertice_area.resize(mesh.n_vertices());
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		int vid = v_it.handle().idx();
		float varea = 0;
		for (Mesh::VertexFaceIter vf_it = mesh.vf_begin(v_it.handle()); vf_it.is_valid(); ++vf_it)
		{
			int vfid = vf_it.handle().idx();
			varea = varea + arealist[vfid];
		}
		vertice_area[vid] = varea;
	}
}

bool  MyRemeshing::My_edge_Filp(Mesh::EdgeHandle& eh, Mesh& mesh)
{
	int a = eh.idx();

	Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh, 0);
	Mesh::HalfedgeHandle oheh = mesh.opposite_halfedge_handle(heh);
	Mesh::HalfedgeHandle h1 = mesh.next_halfedge_handle(heh);
	Mesh::HalfedgeHandle h2 = mesh.next_halfedge_handle(oheh);

	Mesh::VertexHandle v1 = mesh.to_vertex_handle(h1);
	Mesh::VertexHandle v2 = mesh.from_vertex_handle(h1);
	Mesh::VertexHandle v3 = mesh.to_vertex_handle(h2);
	Mesh::VertexHandle v4 = mesh.from_vertex_handle(h2);

	std::vector<Mesh::VertexHandle> vlist;
	std::vector<int> vldeg;
	vlist.resize(4);
	vldeg.resize(4);
	vlist[0] = mesh.from_vertex_handle(h1);
	vlist[1] = mesh.to_vertex_handle(h1);
	vlist[2] = mesh.from_vertex_handle(h2);
	vlist[3] = mesh.to_vertex_handle(h2);

	for (int i = 0; i < 4; i++)
	{
		vldeg[i] = mesh.valence(vlist[i]);
	}

	float nv1 = abs(vldeg[0] - 6) + abs(vldeg[1] - 6) + abs(vldeg[2] - 6) + abs(vldeg[3] - 6);
	float nv2 = abs(vldeg[0] - 7) + abs(vldeg[1] - 5) + abs(vldeg[2] - 7) + abs(vldeg[3] - 5);
	if (nv2 < nv1 && vldeg[0]>4 && vldeg[2]>4)
	{
		if (mesh.is_flip_ok(eh))
		{
			mesh.flip(eh);
		}
	}

	return true;
}
