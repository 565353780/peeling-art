#include "MeshGeneration.h"
#include "Mesh/TriangleInterface.h"
#include "Mesh/surface_manager.h"
#include "Algorithm/BPoly/BPolyline.h"
#include <iostream>
#include <fstream>

extern OpenMesh::FPropHandleT<bool> is_singular;
extern OpenMesh::FPropHandleT<int> face_type;
extern OpenMesh::VPropHandleT<Eigen::Vector3f> vert_ref;
extern OpenMesh::EPropHandleT<int> edge_type;
extern OpenMesh::VPropHandleT<int> vert_type;
extern OpenMesh::VPropHandleT<Eigen::Vector3f> user_vert_ref;
extern OpenMesh::VPropHandleT<int> ori_vert_idx;

MeshGeneration::MeshGeneration()
{
	use_curve_area = 4 * 3.1415926;
	n_boundary = 30;
	rec_boundary = NULL;
	user_pc = NULL;
}


MeshGeneration::~MeshGeneration()
{
	if (user_pc != NULL) delete user_pc;
	if (rec_boundary != NULL) delete rec_boundary;
}

void MeshGeneration::meshFromPolyline(BPolyline* bp, Mesh& mesh_)
{
	user_pc = new BPolyline();
	rec_boundary = new BPolyline();

	Eigen::MatrixXd vert_list;
	Eigen::MatrixXi face_list;

	double area_max = 0.01;

	Eigen::MatrixXd bnd_pts;
	Eigen::MatrixXi edges;

	if (!bp->is_close)
	{
		std::cout << "the curve must be close" << std::endl;
		return;
	}

	int n_all_point = bp->polyline_x.size() - 1;
	int n_all_edge = bp->polyline_x.size() - 1;
	bnd_pts.resize(n_all_point, 2);
	edges.resize(n_all_edge, 2);

	int np = bp->polyline_x.size() - 1;
	for (int j = 0; j < np; j++)
	{
		bnd_pts(j, 0) = bp->polyline_x[j];
		bnd_pts(j, 1) = bp->polyline_y[j];
		edges.row(j) << j, j + 1;
	}
	edges.row(np - 1) = Eigen::RowVector2i(np - 1, 0);

	triangulate(bnd_pts, edges, area_max, vert_list, face_list);

    std::ofstream file("output/obj/loadmesh_p.obj");
	file << "g_object use_mesh" << std::endl;
	for (int i = 0; i < vert_list.rows(); i++)
	{
		file << "v " << vert_list(i, 0) << " " << vert_list(i, 1) << " " << 0 << std::endl;
	}
	for (int i = 0; i < face_list.rows(); i++)
	{
		file << "f " << face_list(i, 0) + 1 << " " << face_list(i, 1) + 1 << " " << face_list(i, 2) + 1 << std::endl;
	}
    file.close();
	
    bool read_OK = OpenMesh::IO::read_mesh(mesh_, "output/obj/loadmesh_p.obj");
	if (!read_OK)
	{
		std::cout << "get mesh fill" << std::endl;
		return;
	}
}

void MeshGeneration::meshFromPolylineList(std::vector<BPolyline* >& bp_list)
{
	Eigen::MatrixXd vert_list;
	Eigen::MatrixXi face_list;

	double area_max = 0.01;

	Eigen::MatrixXd bnd_pts;
	Eigen::MatrixXi edges;

	int n_all_point = 0;
	int n_all_edge = 0;
	for (int i = 0; i < bp_list.size(); i++)
	{
		if (bp_list[i]->is_close)
		{
			int n_p = bp_list[i]->polyline_x.size() - 1;
			n_all_point += n_p;
			n_all_edge += n_p;
		}
		else
		{
			int n_p = bp_list[i]->polyline_x.size();
			n_all_point += n_p;
			n_all_edge += n_p - 1;
		}
	}
	bnd_pts.resize(n_all_point, 2);
	edges.resize(n_all_edge, 2);

	int nn_p = 0, nn_e = 0;
	for (int i = 0; i < bp_list.size(); i++)
	{
		if (bp_list[i]->is_close)
		{
			int np = bp_list[i]->polyline_x.size() - 1;
			for (int j = 0; j < np; j++)
			{
				bnd_pts(j+nn_p, 0) = bp_list[i]->polyline_x[j];
				bnd_pts(j + nn_p, 1) = bp_list[i]->polyline_y[j];
				edges.row(j + nn_e) << j + nn_p, j + nn_p + 1;
			}
			edges.row(nn_e+np-1) = Eigen::RowVector2i(nn_p + np - 1, nn_p);
			nn_p += np;
			nn_e += np;
		}
		else
		{
			int np = bp_list[i]->polyline_x.size();
			for (int j = 0; j < np-1; j++)
			{
				bnd_pts(j + nn_p, 0) = bp_list[i]->polyline_x[j];
				bnd_pts(j + nn_p, 1) = bp_list[i]->polyline_y[j];
				edges.row(j + nn_e) << j + nn_p, j + nn_p + 1;
			}
			bnd_pts(np + nn_p - 1, 0) = bp_list[i]->polyline_x[np - 1];
			bnd_pts(np + nn_p - 1, 1) = bp_list[i]->polyline_y[np - 1];
			nn_p += np;
			nn_e += np - 1;
		}
	}

	triangulate(bnd_pts,edges, area_max, vert_list, face_list);

    std::ofstream file("output/obj/Load_Mesh.obj");
	file << "g_object use_mesh" << std::endl;
	for (int i = 0; i < vert_list.rows(); i++)
	{
		file << "v " << vert_list(i, 0) << " " << vert_list(i, 1) << " " << 0 << std::endl;
	}
	for (int i = 0; i < face_list.rows(); i++)
	{
		file << "f " << face_list(i, 0)+1 << " " << face_list(i, 1)+1 << " " << face_list(i, 2)+1 << std::endl;
	}
    file.close();
}

void MeshGeneration::meshFromChangePolyline(BPolyline* bp, Mesh& ref_mesh)
{
	std::vector<BPolyline* > bp_list(1);
	bp_list[0] = bp;
	meshFromPolylineList(bp_list);
    bool read_OK = OpenMesh::IO::read_mesh(ref_mesh, "output/obj/Load_Mesh.obj");

	for (Mesh::VertexIter v_it = ref_mesh.vertices_begin(); v_it != ref_mesh.vertices_end(); v_it++)
	{
		if (ref_mesh.is_boundary(v_it.handle()))
		{
			ref_mesh.property(vert_type, *v_it) = 0;
		}
		else
		{
			ref_mesh.property(vert_type, *v_it) = -1;
		}
		Mesh::Point mp = ref_mesh.point(v_it.handle());
		Eigen::Vector3f p = Eigen::Vector3f(mp[0], mp[1], mp[2]);
		ref_mesh.property(vert_ref, *v_it) = p;
		ref_mesh.property(user_vert_ref, *v_it) = p;
	}
	for (Mesh::FaceIter f_it = ref_mesh.faces_begin(); f_it != ref_mesh.faces_end(); f_it++)
	{
		ref_mesh.property(is_singular, *f_it) = false;
		ref_mesh.property(face_type, *f_it) = 0;
	}
}
void MeshGeneration::meshFromUserPolyline(std::vector<Eigen::Vector3f>& use_point, Mesh& fix_mesh)
{
	n_boundary = 30;

	double area = 4 * 3.1415926;
	double factor = 1;
	std::vector<BPolyline* > bp_list;
	bp_list.clear();
    if(user_pc != NULL)
    {
        delete user_pc;
    }
	user_pc = new BPolyline();
	user_pc->IntPolyline(use_point);
	user_pc->Centralization(factor);
	user_pc->ChangePolylineAndReduce();

	m_center = (user_pc->rectangle_max + user_pc->rectangle_min) / 2;
	half_length = (user_pc->rectangle_max.x() - user_pc->rectangle_min.x()) / 2;

	Eigen::Vector3f line_start;
	Eigen::Vector3f line_endl;
	std::vector<Eigen::Vector3f> boundary_p;
	boundary_p.clear();
	{
		line_start = Eigen::Vector3f(user_pc->rectangle_min.x(), user_pc->rectangle_min.y(), 0);
		line_endl = Eigen::Vector3f(user_pc->rectangle_min.x(), user_pc->rectangle_max.y(), 0);
		divide_line(line_start, line_endl, n_boundary, boundary_p);
		boundary_p.pop_back();

		line_start = line_endl;
		line_endl = Eigen::Vector3f(user_pc->rectangle_max.x(), user_pc->rectangle_max.y(), 0);
		divide_line(line_start, line_endl, n_boundary, boundary_p);
		boundary_p.pop_back();

		line_start = line_endl;
		line_endl = Eigen::Vector3f(user_pc->rectangle_max.x(), user_pc->rectangle_min.y(), 0);
		divide_line(line_start, line_endl, n_boundary, boundary_p);
		boundary_p.pop_back();

		line_start = line_endl;
		line_endl = Eigen::Vector3f(user_pc->rectangle_min.x(), user_pc->rectangle_min.y(), 0);
		divide_line(line_start, line_endl, n_boundary, boundary_p);
		boundary_p.pop_back();
	}

	bp_list.push_back(user_pc);
	
	meshFromPolylineList(bp_list);

    bool read_OK = OpenMesh::IO::read_mesh(fix_mesh, "output/obj/Load_Mesh.obj");

	assert(read_OK);
	InMeshInline(fix_mesh);

    mesh.add_property(vert_type);
    mesh.add_property(vert_ref);
    mesh.add_property(user_vert_ref);
    mesh.add_property(is_singular);
    mesh.add_property(face_type);

	mesh_scaffold_Inline(user_pc, fix_mesh, mesh);

	mesh.request_vertex_normals();
	mesh.request_face_normals();
	mesh.update_normals();
}



void MeshGeneration::mesh_scaffold_Inline_new(BPolyline* bp, Mesh& fix_mesh, Mesh& fin_mesh, std::vector<std::vector<int>>& insert_list)
{
	n_boundary = 30;

	bp->Centralization(1);

	if (rec_boundary == NULL)
	{
		m_center = (bp->rectangle_max + bp->rectangle_min) / 2;
		half_length = (bp->rectangle_max.x() - bp->rectangle_min.x()) / 2;

		std::vector<Eigen::Vector3f> regect_b;
		regect_b.clear();
		Eigen::Vector3f line_start;
		Eigen::Vector3f line_endl;
		std::vector<Eigen::Vector3f> boundary_p;
		boundary_p.clear();
		{
			line_start = Eigen::Vector3f(bp->rectangle_min.x(), bp->rectangle_min.y(), 0);
			line_endl = Eigen::Vector3f(bp->rectangle_min.x(), bp->rectangle_max.y(), 0);
			divide_line(line_start, line_endl, n_boundary, boundary_p);
			boundary_p.pop_back();

			line_start = line_endl;
			line_endl = Eigen::Vector3f(bp->rectangle_max.x(), bp->rectangle_max.y(), 0);
			divide_line(line_start, line_endl, n_boundary, boundary_p);
			boundary_p.pop_back();

			line_start = line_endl;
			line_endl = Eigen::Vector3f(bp->rectangle_max.x(), bp->rectangle_min.y(), 0);
			divide_line(line_start, line_endl, n_boundary, boundary_p);
			boundary_p.pop_back();

			line_start = line_endl;
			line_endl = Eigen::Vector3f(bp->rectangle_min.x(), bp->rectangle_min.y(), 0);
			divide_line(line_start, line_endl, n_boundary, boundary_p);
			boundary_p.pop_back();
		}

		rec_boundary = new BPolyline();
		rec_boundary->IntPolyline(boundary_p, true);
	}

	std::vector<BPolyline* > bp_list(2);
	bp_list[0] = rec_boundary;
	bp_list[1] = bp;

	Eigen::MatrixXd vert_list;
	Eigen::MatrixXi face_list;

	Eigen::MatrixXd bnd_pts;
	Eigen::MatrixXi edges;


	int n_all_point = 0;
	int n_all_edge = 0;
	for (int i = 0; i < bp_list.size(); i++)
	{
		int n_p = bp_list[i]->polyline_x.size() - 1;
		n_all_point += n_p;
		n_all_edge += n_p;
	}

	int n_all_point_boundy = n_all_point;

	for (int i = 0; i < insert_list.size(); i++)
	{
		insert_list[i].push_back(-1);
	}
	for (int i = 4; i < insert_list.size(); i++)
	{
		std::vector<int> n_id;
		n_id.push_back(-1);
		for (int j = 0; j < insert_list[i].size(); j++)
		{
			n_id.push_back(insert_list[i][j]);
		}
		insert_list[i] = n_id;
	}


	std::vector<Eigen::Vector2i> skeleton_edg;
	std::vector<Eigen::Vector2d> skleton_pts;
	for (int i = 0; i < insert_list.size(); i++)
	{
		for (int j = insert_list[i].size()-1; j > 0; j = j - 2)
		{
			double start_x, start_y;
			double end_x, end_y;
			if (insert_list[i][j] == -1)
			{			
				if (i == 0)
				{
					start_x = 0; start_y = half_length; skeleton_edg.push_back(Eigen::Vector2i(45, n_all_point));
				}
				if (i == 1)
				{
					start_x = half_length; start_y = 0; skeleton_edg.push_back(Eigen::Vector2i(75, n_all_point));
				}
				if (i == 2)
				{
					start_x = 0; start_y = -half_length; skeleton_edg.push_back(Eigen::Vector2i(105, n_all_point));
				}
				if (i == 3)
				{
					start_x = -half_length; start_y = 0; skeleton_edg.push_back(Eigen::Vector2i(15, n_all_point));
				}
				if (i == 4)
				{
					start_x = half_length; start_y = 0; skeleton_edg.push_back(Eigen::Vector2i(75, n_all_point));
				}
				if (i == 5)
				{
					start_x = 0; start_y = -half_length; skeleton_edg.push_back(Eigen::Vector2i(105, n_all_point));
				}
				if (i == 6)
				{
					start_x = -half_length; start_y = 0; skeleton_edg.push_back(Eigen::Vector2i(15, n_all_point));
				}
				if (i == 7)
				{
					start_x = 0; start_y = half_length; skeleton_edg.push_back(Eigen::Vector2i(45, n_all_point));
				}
			}
			else
			{
				start_x = user_pc->polyline_x[insert_list[i][j]];
				start_y = user_pc->polyline_y[insert_list[i][j]];
				skeleton_edg.push_back(Eigen::Vector2i(insert_list[i][j]+120, n_all_point));
			}
			if (insert_list[i][j - 1] == -1)
			{
				if (i == 4)
				{
					end_x = 0; end_y = half_length;
				}
				if (i == 5)
				{
					end_x = half_length; end_y = 0;
				}
				if (i == 6)
				{
					end_x = 0; end_y = -half_length;
				}
				if (i == 7)
				{
					end_x = -half_length; end_y = 0;
				}
			}
			else
			{
				end_x = user_pc->polyline_x[insert_list[i][j - 1]];
				end_y = user_pc->polyline_y[insert_list[i][j - 1]];
			}

			double dis = sqrt((end_x - start_x)*(end_x - start_x) + (end_y - start_y)*(end_y - start_y));
			int d_time = ceil(dis / 0.05);


			for (int k = 1; k < d_time; k++)
			{
				skleton_pts.push_back(Eigen::Vector2d(start_x + (end_x - start_x)*double(k) / d_time, start_y + (end_y - start_y)*double(k) / d_time));
				skeleton_edg.push_back(Eigen::Vector2i(k + n_all_point - 1, k + n_all_point));
			}
			if (insert_list[i][j - 1] == -1)
			{
				if (i == 4)
				{
					skeleton_edg[skeleton_edg.size() - 1](1) = 45;
				}
				if (i == 5)
				{
					skeleton_edg[skeleton_edg.size() - 1](1) = 75;
				}
				if (i == 6)
				{
					skeleton_edg[skeleton_edg.size() - 1](1) = 105;
				}
				if (i == 7)
				{
					skeleton_edg[skeleton_edg.size() - 1](1) = 15;
				}
			}
			else
			{
				skeleton_edg[skeleton_edg.size() - 1](1) = insert_list[i][j - 1] + 120;
			}			
			n_all_point += d_time - 1;
			n_all_edge += d_time;
		}
	}

	bnd_pts.resize(n_all_point, 2);
	edges.resize(n_all_edge, 2);

	int nn_p = 0, nn_e = 0;
	for (int i = 0; i < bp_list.size(); i++)
	{
		int np = bp_list[i]->polyline_x.size() - 1;
		for (int j = 0; j < np; j++)
		{
			bnd_pts(j + nn_p, 0) = bp_list[i]->polyline_x[j];
			bnd_pts(j + nn_p, 1) = bp_list[i]->polyline_y[j];
			edges.row(j + nn_e) << j + nn_p, j + nn_p + 1;
		}
		edges.row(nn_e + np - 1) = Eigen::RowVector2i(nn_p + np - 1, nn_p);
		nn_p += np;
		nn_e += np;
	}

	for (int i = 0; i < skleton_pts.size(); i++)
	{
		bnd_pts(i + nn_p, 0) = skleton_pts[i](0);
		bnd_pts(i + nn_p, 1) = skleton_pts[i](1);

	}
	for (int i = 0; i < skeleton_edg.size(); i++)
	{
		edges(i + nn_p, 0) = skeleton_edg[i](0);
		edges(i + nn_p, 1) = skeleton_edg[i](1);
	}

	int out_bud = rec_boundary->polyline_x.size() - 1;
	int in_bud = bp->polyline_x.size() - 1;
	Mesh::Point p = fix_mesh.point(fix_mesh.vertex_handle(in_bud + 1));

	Eigen::MatrixXd hole;
	hole.resize(1, 2);
//	hole(0, 0) = p[0];
//	hole(0, 1) = p[1];

	hole(0, 0) = 0;
	hole(0, 1) = 0;

	triangluate(bnd_pts, edges, hole, vert_list, face_list);



    std::ofstream file1("output/obj/Load_Mesh.obj");
	file1 << "g_object use_mesh" << std::endl;
	for (int i = 0; i < vert_list.rows(); i++)
	{
		file1 << "v " << vert_list(i, 0) << " " << vert_list(i, 1) << " " << 0 << std::endl;
	}
	for (int i = 0; i < face_list.rows(); i++)
	{
		file1 << "f " << face_list(i, 0) + 1 << " " << face_list(i, 1) + 1 << " " << face_list(i, 2) + 1 << std::endl;
	}
    file1.close();

	int ori_vert_num = fix_mesh.n_vertices();
	int ori_face_num = fix_mesh.n_faces();


    std::ofstream file("output/obj/Load_Mesh.obj");
	file << "g_object use_mesh" << std::endl;
	for (int i = 0; i < out_bud; i++)
	{
		file << "v " << vert_list(i, 0) << " " << vert_list(i, 1) << " " << 0 << std::endl;
	}
	for (int i = 0; i < ori_vert_num; i++)
	{
		Mesh::Point p = fix_mesh.point(fix_mesh.vertex_handle(i));
		file << "v " << p[0] << " " << p[1] << " " << 0 << std::endl;
	}
	for (int i = n_all_point_boundy; i < vert_list.rows(); i++)
	{
		file << "v " << vert_list(i, 0) << " " << vert_list(i, 1) << " " << 0 << std::endl;
	}

	for (Mesh::FaceIter f_it = fix_mesh.faces_begin(); f_it != fix_mesh.faces_end(); f_it++)
	{
		std::vector<int> face_vert_id;
		face_vert_id.clear();
		for (Mesh::FaceVertexIter fv_it = fix_mesh.fv_begin(*f_it); fv_it != fix_mesh.fv_end(*f_it); fv_it++)
		{
			face_vert_id.push_back(fv_it.handle().idx() + out_bud + 1);
		}
		file << "f " << face_vert_id[0] << " " << face_vert_id[1] << " " << face_vert_id[2] << std::endl;
	}

	for (int i = 0; i < face_list.rows(); i++)
	{
		std::vector<int> face_vert_id(3);
		for (int j = 0; j < 3; j++)
		{
			if (face_list(i, j) < n_all_point_boundy)
			{
				face_vert_id[j] = face_list(i, j) + 1;
			}
			else
			{
				face_vert_id[j] = face_list(i, j) - in_bud + ori_vert_num + 1;
			}
		}
		file << "f " << face_vert_id[0] << " " << face_vert_id[1] << " " << face_vert_id[2] << std::endl;
	}
    file.close();


    bool read_OK = OpenMesh::IO::read_mesh(fin_mesh, "output/obj/Load_Mesh.obj");
	assert(read_OK);

	for (Mesh::VertexIter v_it = fin_mesh.vertices_begin(); v_it != fin_mesh.vertices_end(); v_it++)
	{
		int id = v_it.handle().idx();
		if (id >= out_bud && id < out_bud + in_bud)
		{
			fin_mesh.property(vert_type, *v_it) = 0;
			Mesh::VertexHandle fix_vh = fix_mesh.vertex_handle(id - out_bud);
			fin_mesh.property(vert_ref, *v_it) = fix_mesh.property(vert_ref, fix_vh);
			fin_mesh.property(user_vert_ref, *v_it) = fix_mesh.property(user_vert_ref, fix_vh);
		}
		else if (id >= n_all_point && id < ori_vert_num + out_bud)
		{
			fin_mesh.property(vert_type, *v_it) = -1;
			Mesh::VertexHandle fix_vh = fix_mesh.vertex_handle(id - out_bud);
			fin_mesh.property(vert_ref, *v_it) = fix_mesh.property(vert_ref, fix_vh);
			fin_mesh.property(user_vert_ref, *v_it) = fix_mesh.property(user_vert_ref, fix_vh);
		}
		else
		{
			fin_mesh.property(vert_type, *v_it) = 1;
			Mesh::Point mp = fin_mesh.point(v_it.handle());
			Eigen::Vector3f p = Eigen::Vector3f(mp[0], mp[1], mp[2]);
			fin_mesh.property(vert_ref, *v_it) = p;
			fin_mesh.property(user_vert_ref, *v_it) = p;
		}

	}
	for (Mesh::FaceIter f_it = fin_mesh.faces_begin(); f_it != fin_mesh.faces_end(); f_it++)
	{
		fin_mesh.property(face_type, *f_it) = 0;
		int fid = f_it.handle().idx();
		if (fid < ori_face_num) fin_mesh.property(is_singular, *f_it) = false;
		else fin_mesh.property(is_singular, *f_it) = true;
	}


	//mesh.request_vertex_normals();
	//mesh.request_face_normals();
	//mesh.update_normals();
}

void MeshGeneration::mesh_scaffold_Inline(BPolyline* bp, Mesh& fix_mesh, Mesh& fin_mesh)
{
	n_boundary = 30;

	double area = 4 * 3.1415926;
	double factor = 1;
	bp->Centralization(factor);

	if (rec_boundary == NULL)
	{
		m_center = (bp->rectangle_max + bp->rectangle_min) / 2;
		half_length = (bp->rectangle_max.x() - bp->rectangle_min.x()) / 2;

		std::vector<Eigen::Vector3f> regect_b;
		regect_b.clear();
		Eigen::Vector3f line_start;
		Eigen::Vector3f line_endl;
		std::vector<Eigen::Vector3f> boundary_p;
		boundary_p.clear();
		{
			line_start = Eigen::Vector3f(bp->rectangle_min.x(), bp->rectangle_min.y(), 0);
			line_endl = Eigen::Vector3f(bp->rectangle_min.x(), bp->rectangle_max.y(), 0);
			divide_line(line_start, line_endl, n_boundary, boundary_p);
			boundary_p.pop_back();

			line_start = line_endl;
			line_endl = Eigen::Vector3f(bp->rectangle_max.x(), bp->rectangle_max.y(), 0);
			divide_line(line_start, line_endl, n_boundary, boundary_p);
			boundary_p.pop_back();

			line_start = line_endl;
			line_endl = Eigen::Vector3f(bp->rectangle_max.x(), bp->rectangle_min.y(), 0);
			divide_line(line_start, line_endl, n_boundary, boundary_p);
			boundary_p.pop_back();

			line_start = line_endl;
			line_endl = Eigen::Vector3f(bp->rectangle_min.x(), bp->rectangle_min.y(), 0);
			divide_line(line_start, line_endl, n_boundary, boundary_p);
			boundary_p.pop_back();
		}

		rec_boundary = new BPolyline();
		rec_boundary->IntPolyline(boundary_p, true);
	}

	std::vector<BPolyline* > bp_list(2);
	bp_list[0] = rec_boundary;
	bp_list[1] = bp;

	Eigen::MatrixXd vert_list;
	Eigen::MatrixXi face_list;

	double area_max = 0.01;

	Eigen::MatrixXd bnd_pts;
	Eigen::MatrixXi edges;

	int n_all_point = 0;
	int n_all_edge = 0;
	for (int i = 0; i < bp_list.size(); i++)
	{
		int n_p = bp_list[i]->polyline_x.size() - 1;
		n_all_point += n_p;
		n_all_edge += n_p;
	}
	bnd_pts.resize(n_all_point, 2);
	edges.resize(n_all_edge, 2);

	int nn_p = 0, nn_e = 0;
	for (int i = 0; i < bp_list.size(); i++)
	{
		int np = bp_list[i]->polyline_x.size() - 1;
		for (int j = 0; j < np; j++)
		{
			bnd_pts(j + nn_p, 0) = bp_list[i]->polyline_x[j];
			bnd_pts(j + nn_p, 1) = bp_list[i]->polyline_y[j];
			edges.row(j + nn_e) << j + nn_p, j + nn_p + 1;
		}
		edges.row(nn_e + np - 1) = Eigen::RowVector2i(nn_p + np - 1, nn_p);
		nn_p += np;
		nn_e += np;
	}

	int out_bud = rec_boundary->polyline_x.size() - 1;
	int in_bud = bp->polyline_x.size() - 1;
	Mesh::Point p = fix_mesh.point(fix_mesh.vertex_handle(in_bud + 1));

	Eigen::MatrixXd hole;
	hole.resize(1, 2);
	hole(0, 0) = p[0];
	hole(0, 1) = p[1];

	triangluate(bnd_pts, edges, hole, vert_list, face_list);

    std::ofstream file1("output/obj/Load_Mesh.obj");
	file1 << "g_object use_mesh" << std::endl;
	for (int i = 0; i < vert_list.rows(); i++)
	{
		file1 << "v " << vert_list(i, 0) << " " << vert_list(i, 1) << " " << 0 << std::endl;
	}
	for (int i = 0; i < face_list.rows(); i++)
	{
		file1 << "f " << face_list(i, 0) + 1 << " " << face_list(i, 1) + 1 << " " << face_list(i, 2) + 1 << std::endl;
	}
    file1.close();

	int ori_vert_num = fix_mesh.n_vertices();
	int ori_face_num = fix_mesh.n_faces();
	

    std::ofstream file("output/obj/Load_Mesh.obj");
	file << "g_object use_mesh" << std::endl;
	for (int i = 0; i < out_bud; i++)
	{
		file << "v " << vert_list(i, 0) << " " << vert_list(i, 1) << " " << 0 << std::endl;
	}
	for (int i = 0; i < ori_vert_num; i++)
	{
		Mesh::Point p = fix_mesh.point(fix_mesh.vertex_handle(i));
		file << "v " << p[0] << " " << p[1] << " " << 0 << std::endl;
	}
	for (int i = n_all_point; i < vert_list.rows(); i++)
	{
		file << "v " << vert_list(i, 0) << " " << vert_list(i, 1) << " " << 0 << std::endl;
	}

	for (Mesh::FaceIter f_it = fix_mesh.faces_begin(); f_it != fix_mesh.faces_end(); f_it++)
	{
		std::vector<int> face_vert_id;
		face_vert_id.clear();
		for (Mesh::FaceVertexIter fv_it = fix_mesh.fv_begin(*f_it); fv_it != fix_mesh.fv_end(*f_it); fv_it++)
		{
			face_vert_id.push_back(fv_it.handle().idx() + out_bud + 1);
		}
		file << "f " << face_vert_id[0] << " " << face_vert_id[1] << " " << face_vert_id[2] << std::endl;
	}

	for (int i = 0; i < face_list.rows(); i++)
	{
		std::vector<int> face_vert_id(3);
		for (int j = 0; j < 3; j++)
		{
			if (face_list(i, j) < n_all_point)
			{
				face_vert_id[j] = face_list(i, j) + 1;
			}
			else
			{
				face_vert_id[j] = face_list(i, j) - in_bud + ori_vert_num + 1;
			}
		}
		file << "f " << face_vert_id[0] << " " << face_vert_id[1] << " " << face_vert_id[2] << std::endl;
	}
    file.close();

	
    bool read_OK = OpenMesh::IO::read_mesh(fin_mesh, "output/obj/Load_Mesh.obj");
	assert(read_OK);

	for (Mesh::VertexIter v_it = fin_mesh.vertices_begin(); v_it != fin_mesh.vertices_end(); v_it++)
	{
		int id = v_it.handle().idx();
		if (id >= out_bud && id < out_bud + in_bud)
		{
			fin_mesh.property(vert_type, *v_it) = 0;
			Mesh::VertexHandle fix_vh = fix_mesh.vertex_handle(id - out_bud);
			fin_mesh.property(vert_ref, *v_it) = fix_mesh.property(vert_ref, fix_vh);
			fin_mesh.property(user_vert_ref, *v_it) = fix_mesh.property(user_vert_ref, fix_vh);
		}
		else if (id >= n_all_point && id < ori_vert_num + out_bud)
		{
			fin_mesh.property(vert_type, *v_it) = -1;
			Mesh::VertexHandle fix_vh = fix_mesh.vertex_handle(id - out_bud);
			fin_mesh.property(vert_ref, *v_it) = fix_mesh.property(vert_ref, fix_vh);
			fin_mesh.property(user_vert_ref, *v_it) = fix_mesh.property(user_vert_ref, fix_vh);
		}
		else
		{
			fin_mesh.property(vert_type, *v_it) = 1;
			Mesh::Point mp = fin_mesh.point(v_it.handle());
			Eigen::Vector3f p = Eigen::Vector3f(mp[0], mp[1], mp[2]);
			fin_mesh.property(vert_ref, *v_it) = p;
			fin_mesh.property(user_vert_ref, *v_it) = p;
		}
		
	}
	for (Mesh::FaceIter f_it = fin_mesh.faces_begin(); f_it != fin_mesh.faces_end(); f_it++)
	{
		int fid = f_it.handle().idx();
		if (fid < ori_face_num)
		{
			fin_mesh.property(face_type, *f_it) = fix_mesh.property(face_type, fix_mesh.face_handle(fid));
			fin_mesh.property(is_singular, *f_it) = false;
		}
		else
		{
			fin_mesh.property(face_type, *f_it) = 0;
			fin_mesh.property(is_singular, *f_it) = true;
		}
	}
}

void MeshGeneration::divide_line(Eigen::Vector3f& line_start, Eigen::Vector3f& line_end, int n, std::vector<Eigen::Vector3f>& dive_point)
{
	dive_point.reserve(dive_point.size() + n);
	for (int i = 0; i <= n; i++)
	{
		double lambda = double(n - i) / n;
		dive_point.push_back(lambda*line_start + (1 - lambda)*line_end);
	}
}

void MeshGeneration::divide_line_nobe(Eigen::Vector3f& line_start, Eigen::Vector3f& line_end, int n, std::vector<Eigen::Vector3f>& dive_point)
{
	dive_point.reserve(dive_point.size() + n);
	for (int i = 1; i <= n-1; i++)
	{
		double lambda = double(n - i) / n;
		dive_point.push_back(lambda*line_start + (1 - lambda)*line_end);
	}
}

void MeshGeneration::meshInOut()
{
	mesh.add_property(vert_type);
	mesh.add_property(vert_ref);
	mesh.add_property(user_vert_ref);
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		Mesh::Point mp = mesh.point(v_it.handle());
		Eigen::Vector3f p = Eigen::Vector3f(mp[0], mp[1], mp[2]);
		mesh.property(vert_ref, *v_it) = p;
		mesh.property(user_vert_ref, *v_it) = p;
		mesh.property(vert_type, *v_it) = -1*user_pc->PointInCurve(p);		
	}

	mesh.add_property(is_singular);
	mesh.add_property(face_type);
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		mesh.property(face_type, *f_it) = 0;
		Mesh::Point center = Mesh::Point(0, 0, 0);
		double v_in_out = 0;
		for (Mesh::FaceVertexIter fv_it = mesh.fv_begin(*f_it); fv_it != mesh.fv_end(*f_it); fv_it++)
		{
			center += mesh.point(fv_it.handle());
			v_in_out += mesh.property(vert_type, *fv_it);
		}
		center /= 3;
		if (v_in_out == -3) mesh.property(is_singular, *f_it) = false;
		else if (v_in_out == 3) mesh.property(is_singular, *f_it) = true;
		else
		{
			Eigen::Vector3f p = Eigen::Vector3f(center[0], center[1], center[2]);
			int a = user_pc->PointInCurve(p);
			if (a > 0) mesh.property(is_singular, *f_it) = false;
			if (a < 0) mesh.property(is_singular, *f_it) = true;
		}
	}
}

void MeshGeneration::InMeshInline(Mesh& in_mesh)
{
    if(property_added)
    {
        in_mesh.remove_property(vert_type);
        in_mesh.remove_property(vert_ref);
        in_mesh.remove_property(user_vert_ref);
        in_mesh.remove_property(is_singular);
        in_mesh.remove_property(face_type);
    }
	in_mesh.add_property(vert_type);
	in_mesh.add_property(vert_ref);
	in_mesh.add_property(user_vert_ref);

	for (Mesh::VertexIter v_it = in_mesh.vertices_begin(); v_it != in_mesh.vertices_end(); v_it++)
	{
		Mesh::Point mp = in_mesh.point(v_it.handle());
		Eigen::Vector3f p = Eigen::Vector3f(mp[0], mp[1], mp[2]);
		in_mesh.property(vert_ref, *v_it) = p;
		in_mesh.property(user_vert_ref, *v_it) = p;
		if (in_mesh.is_boundary(*v_it)) in_mesh.property(vert_type, *v_it) = 0;
		else in_mesh.property(vert_type, *v_it) = -1;
	}

	in_mesh.add_property(is_singular);
	in_mesh.add_property(face_type);
	for (Mesh::FaceIter f_it = in_mesh.faces_begin(); f_it != in_mesh.faces_end(); f_it++)
	{
		in_mesh.property(face_type, *f_it) = 0;
		in_mesh.property(is_singular, *f_it) = false;
	}
}

void MeshGeneration::plane_to_Octahedron()
{
	Mesh::Point y = Mesh::Point(0, 0, 0);
	double sqrt2 = sqrt(2);
	double length = half_length;
	double high_c = length;
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		Mesh::Point p = mesh.point(v_it.handle());
		if (p[0] + p[1] <= length && p[0] >= 0 && p[1] >= 0)
		{
			double pnorm = p[0] + p[1];
			Mesh::Point v = Mesh::Point(p[0], p[1], high_c - pnorm);
			mesh.set_point(v_it.handle(), v);
		}
		else if (p[1] <= length + p[0] && p[0] <= 0 && p[1] >= 0)
		{
			double pnorm = p[1] - p[0];
			Mesh::Point v = Mesh::Point(p[0], p[1], high_c - pnorm);
			mesh.set_point(v_it.handle(), v);
		}
		else if (p[0] <= length + p[1] && p[0] >= 0 && p[1] <= 0)
		{
			double pnorm = p[0] - p[1];
			Mesh::Point v = Mesh::Point(p[0], p[1], high_c - pnorm);
			mesh.set_point(v_it.handle(), v);
		}
		else if (p[0] + p[1] >= -length && p[0] <= 0 && p[1] <= 0)
		{
			double pnorm = -p[0] - p[1];
			Mesh::Point v = Mesh::Point(p[0], p[1], high_c - pnorm);
			mesh.set_point(v_it.handle(), v);
		}
		else if (p[0] + p[1] >= length && p[0] >= 0 && p[1] >= 0)
		{
			Mesh::Point q = Mesh::Point(y[0] + length, y[1] + length, 0);
			double pnorm = (q[0] - p[0]) + (q[1] - p[1]);
			Mesh::Point v = Mesh::Point(q[1] - p[1], q[0] - p[0], pnorm -high_c);
			mesh.set_point(v_it.handle(), v);
		}
		else if (p[0] + p[1] <= -length && p[0] <= 0 && p[1] <= 0)
		{
			Mesh::Point q = Mesh::Point(y[0] - length, y[1] - length, 0);
			double pnorm = -(q[0] - p[0]) - (q[1] - p[1]);
			Mesh::Point v = Mesh::Point(q[1] - p[1], q[0] - p[0], pnorm - high_c);
			mesh.set_point(v_it.handle(), v);
		}
		else if (p[0] >= p[1] + length && p[0] >= 0 && p[1] <= 0)
		{
			Mesh::Point q = Mesh::Point(y[0] + length, y[1] - length, 0);
			double pnorm = (q[0] - p[0]) - (q[1] - p[1]);
			Mesh::Point v = Mesh::Point(-(q[1] - p[1]), -(q[0] - p[0]), pnorm - high_c);
			mesh.set_point(v_it.handle(), v);
		}
		else if (p[1] >= p[0] + length && p[0] <= 0 && p[1] >= 0)
		{
			Mesh::Point q = Mesh::Point(y[0] - length, y[1] + length, 0);
			double pnorm = -(q[0] - p[0]) + (q[1] - p[1]);
			Mesh::Point v = Mesh::Point(-(q[1] - p[1]), -(q[0] - p[0]), pnorm - high_c);
			mesh.set_point(v_it.handle(), v);
		}
		else
		{
			std::cout << "liuhao" << std::endl;
			Mesh::Point v = Mesh::Point(0, 0, 100);
			mesh.set_point(v_it.handle(), v);
		}
	}
}

void MeshGeneration::WeldingOctahedron_to_Sharp()
{
	std::vector<int> vert_first;
	std::vector<int> vert_second;
	std::vector<std::pair<int, int> > vert_pari;
	vert_pari.clear();
	std::pair<int, int> pa;
	for (int i = 0; i < n_boundary / 2; i++)
	{
		vert_first.push_back(i);
		vert_second.push_back(n_boundary - i);
		vert_first.push_back(i + n_boundary);
		vert_second.push_back(2 * n_boundary - i);
		vert_first.push_back(i + 2 * n_boundary);
		vert_second.push_back(3 * n_boundary - i);
		vert_first.push_back(i + 3 * n_boundary);
		vert_second.push_back(4 * n_boundary - i);
	}
	vert_first[1] = 0;
	vert_first[2] = 0;
	vert_second[3] = vert_first[3];
	vert_first[3] = 0;

	mesh.add_property(ori_vert_idx);
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		mesh.property(ori_vert_idx, *v_it) = v_it.handle().idx();
	}

	std::vector<int> new_idx(mesh.n_vertices(), -1);
	std::vector< std::vector<int> > new_face_ori_id;
	new_face_ori_id.clear();
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		if (mesh.is_boundary(f_it.handle(), true))
		{
			std::vector<int> face_vert(0);
			std::vector<bool> vert_bool(3);
			for (Mesh::FaceVertexIter fv_it = mesh.fv_begin(f_it); fv_it != mesh.fv_end(f_it); fv_it++)
			{
				face_vert.push_back(fv_it.handle().idx());
			}
			for (int i = 0; i < 3; i++)
			{
				vert_bool[i] = find_first(face_vert[i], vert_first, vert_second);
			}
			if (vert_bool[0] || vert_bool[1] || vert_bool[2])
			{
				mesh.delete_face(f_it.handle());
				new_face_ori_id.push_back(face_vert);
			}
		}
	}
	mesh.garbage_collection();

	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		new_idx[mesh.property(ori_vert_idx, *v_it)] = v_it.handle().idx();
	}

	std::cout << mesh.n_vertices() << std::endl;
	for (int i = 0; i < new_face_ori_id.size(); i++)
	{
		int a = new_idx[new_face_ori_id[i][0]];
		int b = new_idx[new_face_ori_id[i][1]];
		int c = new_idx[new_face_ori_id[i][2]];
		Mesh::FaceHandle fh = mesh.add_face(mesh.vertex_handle(a), mesh.vertex_handle(b), mesh.vertex_handle(c));
		mesh.property(is_singular, fh) = true;
		mesh.property(face_type, fh) = 0;
	}
	mesh.remove_property(ori_vert_idx);

	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		Mesh::Point p = mesh.point(v_it.handle());
		p.normalize();
		mesh.set_point(v_it.handle(), p);
	}
}

void MeshGeneration::SharpMeshInit()
{
	mesh.add_property(edge_type);
	for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
	{
		Mesh::EdgeHandle eh = e_it.handle();
		Mesh::HalfedgeHandle heh = mesh.halfedge_handle(eh, 0);
		Mesh::HalfedgeHandle heh1 = mesh.halfedge_handle(eh, 1);
		Mesh::FaceHandle face1 = mesh.face_handle(heh);
		Mesh::FaceHandle face2 = mesh.face_handle(heh1);
		if (mesh.property(is_singular, face1) && mesh.property(is_singular, face2))
		{
			mesh.property(edge_type, *e_it) = 1;
		}
		else if (!mesh.property(is_singular, face1) && !mesh.property(is_singular, face2))
		{
			mesh.property(edge_type, *e_it) = -1;
		}
		else
		{
			Mesh::VertexHandle vh1 = mesh.from_vertex_handle(heh);
			Mesh::VertexHandle vh2 = mesh.to_vertex_handle(heh);
			if (mesh.property(vert_type, vh1) == 0 && mesh.property(vert_type, vh2) == 0)
			{
				mesh.property(edge_type, *e_it) = 0;
			}
		}
	}

	mesh.update_normals();
}

bool MeshGeneration::find_first(int &id, std::vector<int>& first_list, std::vector<int> second_list)
{
	for (int i = 0; i < second_list.size(); i++)
	{
		if (id == second_list[i])
		{
			id = first_list[i];
			return true;
		}
	}
	return false;
}

void MeshGeneration::Project_Surface()
{
	for (Mesh::VertexIter v_it = mesh.vertices_sbegin(); v_it != mesh.vertices_end(); ++v_it)
	{
		Mesh::Point p = mesh.point(v_it.handle());
		Mesh::Point new_p;
		surface_manager->PointProjection(p, new_p);
		mesh.set_point(v_it.handle(), new_p);
	}
}

void MeshGeneration::PlaneMesh_To_Surface(int id)
{
	DelateWrongPoint();
	plane_to_Octahedron();
	WeldingOctahedron_to_Sharp();
	SharpMeshInit();

	mesh.update_normals();
	if (id >= 0 && id < mesh.n_vertices())
	{
		std::cout << "== root ==" << id << "==" << std::endl;
		plane_eye_position = mesh.property(vert_ref, mesh.vertex_handle(id));

		Mesh::Point p0 = mesh.point(mesh.vertex_handle(id));
		Mesh::Point z0(0, 0, 1);

		Mesh::Point az = cross(p0, z0);

		if (az.norm() > 1e-6)
		{
			az.normalize();
			Mesh::Point ax = p0;
			Mesh::Point ay = cross(az, ax);

			Mesh::Point axt = z0;
			Mesh::Point ayt = cross(az, axt);

			for (int i(0); i < mesh.n_vertices(); i++)
			{
				Mesh::VertexHandle vhi = mesh.vertex_handle(i);
				Mesh::Point xi = mesh.point(vhi);

				xi = dot(xi, az) * az + dot(xi, ax) * axt + dot(xi, ay) * ayt;
				mesh.set_point(vhi, xi);
			}
		}
	}
	else
	{
        //std::cout << "No Fix id" << std::endl;
		plane_eye_position = Eigen::Vector3f(0, 0, 100);
	}

	Project_Surface();
	mesh.update_normals();
}

void MeshGeneration::ReSetEye(Mesh& mymesh)
{
	mesh = mymesh;
	int id = -1;
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		if ((mesh.property(vert_ref, *v_it) - plane_eye_position).norm() < 1e-4) 
			id = v_it.handle().idx();
		Mesh::Point p = mesh.point(v_it.handle());
		p.normalize();
		mesh.set_point(v_it.handle(), p);
	}

	mesh.update_normals();
	if (id >= 0 && id < mesh.n_vertices())
	{
		std::cout << "==" << id << "==" << std::endl;
		plane_eye_position = mesh.property(vert_ref, mesh.vertex_handle(id));

		Mesh::VertexHandle vh = mesh.vertex_handle(id);
		Mesh::Point p0 = mesh.point(vh);

		Mesh::Point z0(0, 0, 1);

		Mesh::Point az = cross(p0, z0);

		if (az.norm() > 1e-6)
		{
			az.normalize();
			Mesh::Point ax = p0;
			Mesh::Point ay = cross(az, ax);

			Mesh::Point axt = z0;
			Mesh::Point ayt = cross(az, axt);

			for (int i(0); i < mesh.n_vertices(); i++)
			{
				Mesh::VertexHandle vhi = mesh.vertex_handle(i);
				Mesh::Point xi = mesh.point(vhi);

				xi = dot(xi, az) * az + dot(xi, ax) * axt + dot(xi, ay) * ayt;
				mesh.set_point(vhi, xi);
			}
		}
	}
	else
	{
        //std::cout << "No Fix id" << std::endl;
	}

	Project_Surface();
	mesh.update_normals();
}

void MeshGeneration::DelateWrongPoint()
{
	mesh.add_property(ori_vert_idx);
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		mesh.property(ori_vert_idx, *v_it) = v_it.handle().idx();
	}

	std::vector< std::vector<int> > new_face_ori_id;
	std::vector<int> new_idx(mesh.n_vertices(), -1);
	new_face_ori_id.clear();
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		std::vector<int> v_id;
		v_id.clear();
		for (Mesh::VertexVertexIter vv_it = mesh.vv_begin(v_it); vv_it != mesh.vv_end(v_it); vv_it++)
		{
			v_id.push_back(vv_it.handle().idx());
		}
		bool isonb = v_id[0] < 120 || v_id[1] < 120 || v_id[2] < 120;
		if (v_id.size() == 3 && isonb && !mesh.is_boundary(v_it))
		{
			new_face_ori_id.push_back(v_id);
			mesh.delete_vertex(v_it.handle());
		}
	}
	mesh.garbage_collection();

	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		new_idx[mesh.property(ori_vert_idx, *v_it)] = v_it.handle().idx();
	}
	for (int i = 0; i < new_face_ori_id.size(); i++)
	{
		int a = new_idx[new_face_ori_id[i][0]];
		int b = new_idx[new_face_ori_id[i][1]];
		int c = new_idx[new_face_ori_id[i][2]];
		Mesh::FaceHandle fh = mesh.add_face(mesh.vertex_handle(c), mesh.vertex_handle(b), mesh.vertex_handle(a));
		mesh.property(is_singular, fh) = true;
		mesh.property(face_type, fh) = 0;
	}
	mesh.remove_property(ori_vert_idx);
}

void MeshGeneration::meshFromUserPolyline_picture(std::vector<Eigen::Vector3f>& use_point, Mesh& fix_mesh)
{
	n_boundary = 30;

	double area = 4 * 3.1415926;
	double factor = 1;
	std::vector<BPolyline* > bp_list;
	bp_list.clear();
	user_pc = new BPolyline();
	user_pc->IntPolyline(use_point);
	user_pc->Centralization(factor);
	user_pc->ChangePolylineAndReduce();

	std::vector<std::vector<int>> insert_list;
	user_pc->Add_insert_point(insert_list);

	m_center = (user_pc->rectangle_max + user_pc->rectangle_min) / 2;
	half_length = (user_pc->rectangle_max.x() - user_pc->rectangle_min.x()) / 2;

	Eigen::MatrixXd vert_list;
	Eigen::MatrixXi face_list;
	Eigen::MatrixXd bnd_pts;
	Eigen::MatrixXi edges;

	int n_all_point = user_pc->polyline_x.size() - 1;
	int n_all_edge = user_pc->polyline_x.size() - 1;

	std::vector<Eigen::Vector2i> skeleton_edg;
	std::vector<Eigen::Vector2d> skleton_pts;
	skleton_pts.push_back(Eigen::Vector2d(0, 0));
	n_all_point++;
	int n_spline = n_all_point - 1;

	for (int i = 0; i < insert_list.size(); i++)
	{
		for (int j = 0; j < insert_list[i].size(); j = j + 2)
		{
			double start_x, start_y;
			double end_x, end_y;
			if (insert_list[i][j] == -1)
			{
				start_x = 0;
				start_y = 0;
				skeleton_edg.push_back(Eigen::Vector2i(n_spline, n_all_point));
			}
			else
			{
				start_x = user_pc->polyline_x[insert_list[i][j]];
				start_y = user_pc->polyline_y[insert_list[i][j]];
				skeleton_edg.push_back(Eigen::Vector2i(insert_list[i][j], n_all_point));
			}
			end_x = user_pc->polyline_x[insert_list[i][j + 1]];
			end_y = user_pc->polyline_y[insert_list[i][j + 1]];

			double dis = sqrt((end_x - start_x)*(end_x - start_x) + (end_y - start_y)*(end_y - start_y));
			int d_time = ceil(dis / 0.04);


			for (int k = 1; k < d_time; k++)
			{
				skleton_pts.push_back(Eigen::Vector2d(start_x + (end_x - start_x)*double(k) / d_time, start_y + (end_y - start_y)*double(k) / d_time));
				skeleton_edg.push_back(Eigen::Vector2i(k + n_all_point - 1, k + n_all_point));
			}
			skeleton_edg[skeleton_edg.size() - 1](1) = insert_list[i][j + 1];
			n_all_point += d_time - 1;
			n_all_edge += d_time;
		}
	}

	bnd_pts.resize(n_all_point, 2);
	edges.resize(n_all_edge, 2);

	int np = user_pc->polyline_x.size() - 1;
	for (int j = 0; j < np; j++)
	{
		bnd_pts(j, 0) = user_pc->polyline_x[j];
		bnd_pts(j, 1) = user_pc->polyline_y[j];
		edges.row(j) << j, j + 1;
	}
	edges.row(np - 1) = Eigen::RowVector2i(np - 1, 0);

	for (int i = 0; i < skleton_pts.size(); i++)
	{
		bnd_pts(i + np, 0) = skleton_pts[i](0);
		bnd_pts(i + np, 1) = skleton_pts[i](1);

	}
	for (int i = 0; i < skeleton_edg.size(); i++)
	{
		edges(i + np, 0) = skeleton_edg[i](0);
		edges(i + np, 1) = skeleton_edg[i](1);
	}

	std::cout << bnd_pts.row(bnd_pts.rows() - 1);
	std::cout << edges.row(edges.rows() - 1);

	triangulate(bnd_pts, edges, 0.01, vert_list, face_list);

    std::ofstream file("output/obj/Load_Mesh.obj");
	file << "g_object use_mesh" << std::endl;
	for (int i = 0; i < vert_list.rows(); i++)
	{
		file << "v " << vert_list(i, 0) << " " << vert_list(i, 1) << " " << 0 << std::endl;
	}
	for (int i = 0; i < face_list.rows(); i++)
	{
		file << "f " << face_list(i, 0) + 1 << " " << face_list(i, 1) + 1 << " " << face_list(i, 2) + 1 << std::endl;
	}
    file.close();

    bool read_OK = OpenMesh::IO::read_mesh(fix_mesh, "output/obj/Load_Mesh.obj");

	assert(read_OK);
	InMeshInline(fix_mesh);

	mesh.add_property(vert_type);
	mesh.add_property(vert_ref);
	mesh.add_property(user_vert_ref);
	mesh.add_property(is_singular);
	mesh.add_property(face_type);

	mesh_scaffold_Inline_new(user_pc, fix_mesh, mesh, insert_list);


	mesh.request_vertex_normals();
	mesh.request_face_normals();
	mesh.update_normals();
}
