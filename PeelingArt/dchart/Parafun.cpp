#include "Parafun.h"
#include "Mesh/TriangleInterface.h"

extern OpenMesh::VPropHandleT<Eigen::Vector3f> vert_ref;

Parafun::Parafun()
{

}

Parafun::Parafun(Mesh& mesh_, int p_type, vector<double>& s_p00, vector<double>& s_p01, vector<double>& s_p10, vector<double>& s_p11)
{
	mesh = mesh_;
	process_type = p_type;

	F_N = mesh.n_faces();
	total_num = mesh.n_vertices();

	area.resize(F_N);
	F0.resize(F_N);
	F1.resize(F_N);
	F2.resize(F_N);

	negative_grad_norm.resize(2 * total_num);

	mesh.request_face_normals();
	mesh.update_normals();
	pardiso = NULL;
	convgence_con_rate = 1e-6;
	MAX_ITER_NUM = 5000;
	bound_distortion_K = 250;


	source_p00 = s_p00;
	source_p01 = s_p01;
	source_p10 = s_p10;
	source_p11 = s_p11;
}

Parafun::Parafun(Mesh& mesh_, int p_type)
{
	mesh = mesh_;
	process_type = p_type;

	F_N = mesh.n_faces();
	total_num = mesh.n_vertices();

	area.resize(F_N);
	F0.resize(F_N);
	F1.resize(F_N);
	F2.resize(F_N);

	negative_grad_norm.resize(2 * total_num);

	mesh.request_face_normals();
	mesh.update_normals();
	pardiso = NULL;
	convgence_con_rate = 1e-6;
	MAX_ITER_NUM = 5000;
	bound_distortion_K = 250;
}

Parafun::~Parafun()
{
}

bool Parafun::ChangeVertPoint(Mesh& out_mesh)
{
	for (auto it1 = out_mesh.vertices_begin(); it1 != out_mesh.vertices_end(); ++it1)
	{
		int id = it1->idx();
		OpenMesh::Vec3d pos(originmesh_area_sqrt*position_of_mesh(id), originmesh_area_sqrt*position_of_mesh(id +total_num), 0);
		out_mesh.set_point(*it1, pos);
	}
	out_mesh.update_normals();
	return true;
}
bool Parafun::ChangeVertRef(Mesh& out_mesh, int zero_idx)
{
	Vector3f move;
	if (zero_idx == -1)
	{
		move = Vector3f(0, 0, 0);
	}
	else
	{
		move = Vector3f(originmesh_area_sqrt*position_of_mesh(zero_idx), originmesh_area_sqrt*position_of_mesh(zero_idx + total_num), 0);
	}
	for (auto it1 = out_mesh.vertices_begin(); it1 != out_mesh.vertices_end(); ++it1)
	{
		int id = it1->idx();
		Vector3f pos(originmesh_area_sqrt*position_of_mesh(id), originmesh_area_sqrt*position_of_mesh(id + total_num), 0);
		out_mesh.property(vert_ref, *it1) = pos - move;
	}
	return true;
}

void Parafun::init_Deformation(std::set<int>& my_fix_id, std::vector<int>& my_handle_id, std::vector<int>& my_face_in_out)
{
	originmesh_area_sqrt = 1.0;
	Eigen::MatrixX3i s_T;
	s_T.resize(mesh.n_faces(), 3);
	for (int i = 0; i < F_N; ++i)
	{
		auto face = mesh.face_handle(i);
		int dd = 0;
		for (Mesh::FaceVertexIter it2 = mesh.fv_begin(face); it2 != mesh.fv_end(face); ++it2)
		{
			auto vertex_ = it2.handle();
			switch (dd)
			{
			case 0: F0[i] = vertex_.idx(); break;
			case 1: F1[i] = vertex_.idx(); break;
			case 2: F2[i] = vertex_.idx(); break;
			default:
				break;
			}
			dd++;
		}
		s_T.row(i) = RowVector3i(F0[i], F1[i], F2[i]);
	}

	position_of_mesh.resize(2 * total_num);
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		int i = v_it.handle().idx();
		position_of_mesh(i) = mesh.point(v_it.handle())[0];
		position_of_mesh(i + total_num) = mesh.point(v_it.handle())[1];
	}

	total_num = mesh.n_vertices();
	//��ʼ����̶�������� and fix_id
	fix_number = my_fix_id.size();
	if (fix_number != 0)
	{
		fixed_ids.resize(my_fix_id.size());
		int i = 0;
		for (set<int>::iterator iter = my_fix_id.begin(); iter != my_fix_id.end(); iter++)
		{
			fixed_ids(i) = *iter;
			i++;
		}
	}
	V_N = total_num - my_fix_id.size();


	assert(my_face_in_out.size() == F_N);
	face_in_out = my_face_in_out;

	id_vs_index();

	Pre_calculate(s_T);
	//Tutte();
	init_area();

	handle_ids.clear();
	for (int i = 0; i < my_handle_id.size(); i++)
	{
		handle_ids.push_back(my_handle_id[i]);
	}
}
void Parafun::init_handle(std::vector<double>& handle_new_position_x, std::vector<double>& handle_new_position_y)
{
	int n = handle_new_position_x.size();
	handle_target_pos.resize(n, 2);
	for (int i = 0; i < n; i++)
	{
		handle_target_pos(i, 0) = handle_new_position_x[i];
		handle_target_pos(i, 1) = handle_new_position_y[i];
	}
}

void Parafun::init_Parametric(std::set<int>& my_fix_id)
{
	originmesh_area_sqrt = 1.0;
	Eigen::MatrixX3i s_T;
	s_T.resize(mesh.n_faces(), 3);
	for (int i = 0; i < F_N; ++i)
	{
		auto face = mesh.face_handle(i);
		Mesh::FaceVertexIter fv_it = mesh.fv_begin(face); F0[i] = (*fv_it).idx();
		++fv_it; F1[i] = (*fv_it).idx();
		++fv_it; F2[i] = (*fv_it).idx();

		s_T.row(i) = RowVector3i(F0[i], F1[i], F2[i]);
	}

	position_of_mesh.resize(2 * total_num);
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		int i = v_it.handle().idx();
		Vector3f ref = mesh.property(vert_ref, *v_it);
		position_of_mesh(i) = ref.x();
		position_of_mesh(i + total_num) = ref.y();
	}

	//��ʼ����̶���� fix_id
	fix_number = my_fix_id.size();
	if (fix_number != 0)
	{
		fixed_ids.resize(my_fix_id.size());
		int i = 0;
		for (set<int>::iterator iter = my_fix_id.begin(); iter != my_fix_id.end(); iter++)
		{
			fixed_ids(i) = *iter;
			i++;
		}	
	}
	V_N = total_num - my_fix_id.size();

	face_in_out.resize(F_N, 1);

	id_vs_index();

	Pre_calculate(s_T);
	if (process_type <= 2)
	{
		Tutte();
	}
	
	init_area();
}

void Parafun::recover_to_src()
{
	area = area_src;
	update_p00 = source_p00;
	update_p01 = source_p01;
	update_p10 = source_p10;
	update_p11 = source_p11;
}

void Parafun::init_area()
{
	area.resize(F_N);
	area_uniform.resize(F_N);
	area_src.resize(F_N);
	for (int i = 0; i < F_N; ++i)
	{
		auto face = mesh.face_handle(i);
		if (face_in_out[i] > 0)
		{
			area_src[i] = mesh.calc_sector_area(mesh.halfedge_handle(face));
		}
		else
		{
			area_src[i] = 0.0001;
		}	
	}
	fill(area_uniform.begin(), area_uniform.end(), 1.0 / F_N);
	//fill(area.begin(), area.end(), 1.0 / F_N);
	area = area_src;
}

bool Parafun::run_cm()
{
	if (pardiso != NULL)
	{
		delete pardiso;
		pardiso = NULL;
	}
	pardiso = new PardisoSolver();
	pardiso->ia = pardiso_ia;
	pardiso->ja = pardiso_ja;
	pardiso->a.resize(pardiso_ja.size());
	pardiso->nnz = pardiso_ja.size();
	pardiso->num = 2 * V_N;

	pardiso->pardiso_init();

	vector<double> energy_area_process;
	energy_area_process.reserve(MAX_ITER_NUM);

	Energysource();

	energy_area_process.push_back(energy_area);
	double energy_pre = 0;
	double energy_cur = energy_uniform;

	int iter_num_cur = 0;

	int sum_iter_num = 0;

	double conv_percent = 1;

	g_norm = 1.0;

	long time_beg, time_end;
	time_beg = clock();

	bool is_cm_ok = true;

	recover_to_src();
	energy_cur = energy_area;
	while (iter_num_cur < MAX_ITER_NUM)
	{
		iter_num_cur++;
		energy_pre = energy_cur;

		if (process_type == 1 || process_type == 3)
		{
			CM();
		}
		else if (process_type == 2 || process_type == 4)
		{
//			CM();
			is_cm_ok = SLIM();
			if (!is_cm_ok)
			{
				break;
			}
		}
		else
		{
			std::cout << "liuhao wrong type" << std::endl;
		}
		energy_area_process.push_back(energy_area);
		sum_iter_num++;

		energy_cur = energy_area;
		conv_percent = abs(energy_cur - energy_pre) / energy_pre;
		//calc_gradient_norm(position_of_mesh);
		if (conv_percent <= convgence_con_rate /*|| g_norm <= convgence_con_rate8*/)
		{
			break;
		}
	}

	time_end = clock();
	time_consumption = (time_end - time_beg) / 1000.0;

	cout << "PP ====== time_consumption: " << time_consumption << " s; sum_iter: " << sum_iter_num << endl;

	delete pardiso;
	pardiso = NULL;	

	return is_cm_ok;
}

void Parafun::Tutte()
{
	int boundary_num = 0;
	auto it1 = mesh.halfedges_begin();
	while (!mesh.is_boundary(*it1))
		it1++;
	auto he_start = *it1;
	auto he_it = he_start;
	do
	{
		he_it = mesh.next_halfedge_handle(he_it);
		boundary_num++;
	} while (he_it != he_start);

	double delta_angle = 2 * M_PI / boundary_num;
	double area_1_factor = sqrt(1.0 / M_PI);

	position_of_mesh.resize(2 * V_N);
	for (int i = 0; i < boundary_num; ++i)
	{
		auto v_h = mesh.to_vertex_handle(he_start);
		position_of_mesh(v_h.idx()) = area_1_factor*cos(i * delta_angle);
		position_of_mesh(v_h.idx() + V_N) = area_1_factor*sin(-i * delta_angle);		
		he_start = mesh.next_halfedge_handle(he_start);
	}

	vector<int> pardiso_it;
	vector<int> pardiso_jt;
	vector<double> pardiso_t;
	vector<double> pardiso_tu;
	vector<double> pardiso_tv;

	pardiso_it.reserve(V_N+1);
	pardiso_jt.reserve(6*V_N);
	pardiso_t.reserve(6*V_N);
	pardiso_tu.resize(V_N, 0.0);
	pardiso_tv.resize(V_N, 0.0);
	for (size_t i = 0; i < V_N; i++)
	{
		pardiso_it.push_back(pardiso_jt.size());

		auto v_h = mesh.vertex_handle(i);
		if (mesh.is_boundary(v_h))
		{
			pardiso_jt.push_back(i);
			pardiso_t.push_back(1);

			pardiso_tu[i]=position_of_mesh(i);
			pardiso_tv[i]=position_of_mesh(i + V_N);

		}
		else
		{
			pardiso_jt.push_back(i);
			pardiso_t.push_back(mesh.valence(v_h));
			vector<int> row_id;
			row_id.reserve(mesh.valence(v_h));
			double bu = 0.0; double bv = 0.0;
			for (auto it2 = mesh.vv_begin(v_h); it2 != mesh.vv_end(v_h); ++it2)
			{
				int vv_id = it2->idx();
				if (mesh.is_boundary(*it2))
				{
					bu += position_of_mesh(vv_id);
					bv += position_of_mesh(vv_id+V_N);
				}
				else
				{
					if (vv_id>i)
					{
						row_id.push_back(vv_id);
					}
				}
			}
			sort(row_id.begin(), row_id.end(), less<int>());
			for (size_t j = 0; j < row_id.size(); j++)
			{
				pardiso_jt.push_back(row_id[j]);
				pardiso_t.push_back(-1);
			}
			pardiso_tu[i] = bu;
			pardiso_tv[i] = bv;
		}
	}
	pardiso_it.push_back(pardiso_jt.size());

	if (pardiso != NULL)
	{
		delete pardiso;
		pardiso = NULL;
	}
	pardiso = new PardisoSolver();
	pardiso->ia = pardiso_it;
	pardiso->ja = pardiso_jt;
	pardiso->nnz = pardiso_jt.size();
	pardiso->num = V_N;

	pardiso->pardiso_init();

	pardiso->a = pardiso_t;

	pardiso->rhs = pardiso_tu;
	pardiso->factorize();


	pardiso->pardiso_solver();

	for (size_t i = 0; i < V_N; i++)
	{
		position_of_mesh(i) = (pardiso->result)[i];
	}

	pardiso->rhs = pardiso_tv;
	pardiso->pardiso_solver();
	for (size_t i = 0; i < V_N; i++)
	{
		position_of_mesh(i+V_N) = (pardiso->result)[i];
	}

	delete pardiso;
	pardiso = NULL;
}

void Parafun::Pre_calculate(Eigen::MatrixX3i& s_T)
{
	if (process_type != 2)
	{
		source_p00.resize(F_N);
		source_p01.resize(F_N);
		source_p10.resize(F_N);
		source_p11.resize(F_N);

		for (int i = 0; i < F_N; ++i)
		{
			double p00, p01, p10, p11;
			local_coordinate_inverse(i, p00, p01, p10, p11);

			source_p00[i] = p00;
			source_p01[i] = p01;
			source_p10[i] = p10;
			source_p11[i] = p11;
		}
	}


	pardiso_ia.clear(); pardiso_ia.reserve(2 * V_N + 1);
	pardiso_ja.clear(); pardiso_ja.reserve(8 * V_N);

	typedef Triplet<int> T;
	std::vector<T> tripletlist;

	std::vector<std::set<int>> VV_tmp;
	VV_tmp.resize(V_N);
	//for (size_t i = 0; i < d_.m_T.rows(); i++)
	//{
	//	int vid[3];

	//	for (size_t j = 0; j < d_.m_T.cols(); j++)
	//	{
	//		vid[j] = id2index[d_.m_T(i, j)];
	//	}
	//	VV_tmp[vid[0]].insert(vid[1]);
	//	VV_tmp[vid[0]].insert(vid[2]);

	//	VV_tmp[vid[1]].insert(vid[0]);
	//	VV_tmp[vid[1]].insert(vid[2]);

	//	VV_tmp[vid[2]].insert(vid[0]);
	//	VV_tmp[vid[2]].insert(vid[1]);
	//}

	vector<int> s_vid;
	for (size_t i = 0; i < F_N; i++)
	{
		s_vid.clear();
		for (size_t j = 0; j < 3; j++)
		{
			int s_id = id2index[s_T(i, j)];
			if (s_id != -1)
				s_vid.push_back(s_id);
		}
		if (s_vid.size() <= 1)
			continue;
		if (s_vid.size() == 2)
		{
			VV_tmp[s_vid[0]].insert(s_vid[1]);
			VV_tmp[s_vid[1]].insert(s_vid[0]);
		}
		else
		{
			VV_tmp[s_vid[0]].insert(s_vid[1]);
			VV_tmp[s_vid[0]].insert(s_vid[2]);

			VV_tmp[s_vid[1]].insert(s_vid[0]);
			VV_tmp[s_vid[1]].insert(s_vid[2]);

			VV_tmp[s_vid[2]].insert(s_vid[0]);
			VV_tmp[s_vid[2]].insert(s_vid[1]);
		}
	}

	for (int i = 0; i < V_N; i++)
	{
		pardiso_ia.push_back(pardiso_ja.size());
		VV_tmp[i].insert(i);
		vector<int> row_id;
		for (auto&var : VV_tmp[i])
		{
			row_id.push_back(var);
		}

		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);

		int dd = 0;
		for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
		{
			pardiso_ja.push_back(row_id[k]);
			tripletlist.push_back(T(i, row_id[k], dd));
			++dd;
		}
		for (int k = 0; k < row_id.size(); k++)
		{
			pardiso_ja.push_back(row_id[k] + V_N);
			tripletlist.push_back(T(i, row_id[k] + V_N, dd));
			++dd;
		}
	}
	for (int i = V_N; i < 2 * V_N; i++)
	{
		pardiso_ia.push_back(pardiso_ja.size());
		vector<int> row_id;
		for (auto&var : VV_tmp[i - V_N])
		{
			row_id.push_back(var);
		}
		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i - V_N);

		int dd = 0;
		for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
		{
			pardiso_ja.push_back(row_id[k] + V_N);
			tripletlist.push_back(T(i, row_id[k] + V_N, dd));
			++dd;
		}
	}

	SparseMatrix<int> find_id_in_rows;
	find_id_in_rows.resize(2 * V_N, 2 * V_N);
	find_id_in_rows.setFromTriplets(tripletlist.begin(), tripletlist.end());

	pardiso_ia.push_back(pardiso_ja.size());

	id_h00.resize(F_N, -1); id_h01.resize(F_N, -1); id_h02.resize(F_N, -1); id_h03.resize(F_N, -1); id_h04.resize(F_N, -1); id_h05.resize(F_N, -1);
	id_h11.resize(F_N, -1); id_h12.resize(F_N, -1); id_h13.resize(F_N, -1); id_h14.resize(F_N, -1); id_h15.resize(F_N, -1);
	id_h22.resize(F_N, -1); id_h23.resize(F_N, -1); id_h24.resize(F_N, -1); id_h25.resize(F_N, -1);
	id_h33.resize(F_N, -1); id_h34.resize(F_N, -1); id_h35.resize(F_N, -1);
	id_h44.resize(F_N, -1); id_h45.resize(F_N, -1);
	id_h55.resize(F_N, -1);

	//for (int i = 0; i < d_.m_T.rows(); i++)
	//{
	//	int f0 = id2index[F0[i]]; int f1 = id2index[F1[i]]; int f2 = id2index[F2[i]]; int f3 = f0 + V_N; int f4 = f1 + V_N; int f5 = f2 + V_N;
	//	if (f0 == -1 || f1 == -1 || f2 == -1)
	//	{
	//		cout << "fuckfuckfuck!!!!" << endl;
	//	}
	//	int min01 = min(f0, f1); int max01 = f0 + f1 - min01;
	//	int min02 = min(f0, f2); int max02 = f0 + f2 - min02;
	//	int min12 = min(f1, f2); int max12 = f1 + f2 - min12;
	//	id_h00[i] = pardiso_ia[f0]; id_h01[i] = pardiso_ia[min01] + find_id_in_rows.coeff(min01, max01); id_h02[i] = pardiso_ia[min02] + find_id_in_rows.coeff(min02, max02);
	//	id_h03[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f3); id_h04[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f4); id_h05[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f5);
	//	id_h11[i] = pardiso_ia[f1]; id_h12[i] = pardiso_ia[min12] + find_id_in_rows.coeff(min12, max12);
	//	id_h13[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f3); id_h14[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f4); id_h15[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f5);
	//	id_h22[i] = pardiso_ia[f2];
	//	id_h23[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f3); id_h24[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f4); id_h25[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f5);
	//	id_h33[i] = pardiso_ia[f3]; id_h34[i] = pardiso_ia[min01 + V_N] + find_id_in_rows.coeff(min01 + V_N, max01 + V_N); id_h35[i] = pardiso_ia[min02 + V_N] + find_id_in_rows.coeff(min02 + V_N, max02 + V_N);
	//	id_h44[i] = pardiso_ia[f4]; id_h45[i] = pardiso_ia[min12 + V_N] + find_id_in_rows.coeff(min12 + V_N, max12 + V_N);
	//	id_h55[i] = pardiso_ia[f5];
	//}

	for (int i = 0; i < F_N; i++)
	{
		int f0 = id2index[F0[i]]; int f1 = id2index[F1[i]]; int f2 = id2index[F2[i]];
		int f3 = f0 + V_N; int f4 = f1 + V_N; int f5 = f2 + V_N;
		if (f0 != -1)
		{
			id_h00[i] = pardiso_ia[f0];
			id_h33[i] = pardiso_ia[f3];
			id_h03[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f3);
		}
		if (f1 != -1)
		{
			id_h11[i] = pardiso_ia[f1];
			id_h44[i] = pardiso_ia[f4];
			id_h14[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f4);
		}
		if (f2 != -1)
		{
			id_h22[i] = pardiso_ia[f2];
			id_h55[i] = pardiso_ia[f5];
			id_h25[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f5);
		}

		if (f1 != -1 && f2 != -1)
		{
			int min12 = min(f1, f2); int max12 = f1 + f2 - min12;
			id_h12[i] = pardiso_ia[min12] + find_id_in_rows.coeff(min12, max12);
			id_h15[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f5);
			id_h24[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f4);
			id_h45[i] = pardiso_ia[min12 + V_N] + find_id_in_rows.coeff(min12 + V_N, max12 + V_N);
		}
		if (f0 != -1 && f2 != -1)
		{
			int min02 = min(f0, f2); int max02 = f0 + f2 - min02;
			id_h02[i] = pardiso_ia[min02] + find_id_in_rows.coeff(min02, max02);
			id_h05[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f5);
			id_h23[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f3);
			id_h35[i] = pardiso_ia[min02 + V_N] + find_id_in_rows.coeff(min02 + V_N, max02 + V_N);
		}
		if (f1 != -1 && f0 != -1)
		{
			int min01 = min(f0, f1); int max01 = f0 + f1 - min01;
			id_h01[i] = pardiso_ia[min01] + find_id_in_rows.coeff(min01, max01);
			id_h04[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f4);
			id_h13[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f3);
			id_h34[i] = pardiso_ia[min01 + V_N] + find_id_in_rows.coeff(min01 + V_N, max01 + V_N);
		}

	}


}

//void Parafun::Update_source_same_t()
//{
//	double t_min = 1;
//	int geqK = -1;
//
//	vector<double> all_s0; all_s0.resize(F_N);
//	vector<double> all_s1; all_s1.resize(F_N);
//
//	vector<double> all_w00; all_w00.resize(F_N);
//	vector<double> all_w01; all_w01.resize(F_N);
//	vector<double> all_w10; all_w10.resize(F_N);
//	vector<double> all_w11; all_w11.resize(F_N);
//
//
//	int f0, f1, f2;
//	double x0, y0, x1, y1, x2, y2;
//	double det;
//	double E_d;
//	double tt;
//	double new_sig0, new_sig1;
//	double j00, j01, j10, j11;
//	double p00, p01, p10, p11;
//	double q00, q01, q10, q11;
//
//	double *position = position_of_mesh.data();
//
//	for (int i = 0; i < F_N; ++i)
//	{
//		f0 = F0[i];
//		f1 = F1[i];
//		f2 = F2[i];
//
//		x0 = position[f0];
//		y0 = position[f0 + V_N];
//
//		x1 = position[f1];
//		y1 = position[f1 + V_N];
//
//		x2 = position[f2];
//		y2 = position[f2 + V_N];
//
//		q00 = x1 - x0; q01 = x2 - x0;
//		q10 = y1 - y0; q11 = y2 - y0;
//
//		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];
//
//		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;
//
//
//		det = j00*j11 - j01*j10;
//		E_d = (1 + 1 / (det*det)) * (j00*j00 + j01*j01 + j10*j10 + j11*j11);
//
//		double alpha_0 = j00 + j11; double alpha_1 = j10 - j01;
//		double beta_0 = j00 - j11; double beta_1 = j10 + j01;
//
//		double alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1*alpha_1);
//		double beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1*beta_1);
//
//		double sig0 = alpha_norm + beta_norm;
//		double sig1 = alpha_norm - beta_norm;
//		all_s0[i] = sig0;
//		all_s1[i] = sig1;
//
//		double temp = 1 / (sig1*sig1 - sig0*sig0);
//
//		all_w00[i] = temp*(j00*j00 + j10*j10 - 0.5*(sig0*sig0 + sig1*sig1));
//		all_w01[i] = temp*(j00*j01 + j10*j11);
//		all_w10[i] = temp*(j01*j00 + j11*j10);
//		all_w11[i] = temp*(j01*j01 + j11*j11 - 0.5*(sig0*sig0 + sig1*sig1));
//
//		if (E_d<=bound_distortion_K)
//		{
//			geqK++;
//		}
//		else
//		{
//			tt = newton_equation(sig0, sig1, bound_distortion_K);
//			if (tt < t_min)
//			{
//				t_min = tt;
//			}
//		}
//	}
//
//	changetocm_flag = (geqK + 1.0) / F_N;
//
//	for (int i = 0; i < F_N; ++i)
//	{
//		double sig0 = all_s0[i];
//		double sig1 = all_s1[i];
//
//		new_sig0 = pow(sig0, t_min - 1);
//		new_sig1 = pow(sig1, t_min - 1);
//
//		double delta_new = new_sig1 - new_sig0;
//		double plus_new = 0.5*(new_sig1 + new_sig0);
//
//		double w00 = delta_new*all_w00[i] + plus_new;
//		double w01 = delta_new*all_w01[i];
//		double w10 = delta_new*all_w10[i];
//		double w11 = delta_new*all_w11[i] + plus_new;
//
//		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];
//
//		update_p00[i] = p00*w00 + p01*w10;
//		update_p01[i] = p00*w01 + p01*w11;
//		update_p10[i] = p10*w00 + p11*w10;
//		update_p11[i] = p10*w01 + p11*w11;
//	}
//
//	Intp_T_Min = t_min;
//}

bool Parafun::SLIM()
{
	cout << "perform_iteration_slim begin-----" << endl;

	double area_now;
	int f0, f1, f2;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	double x0, y0, x1, y1, x2, y2;

	double alpha_norm, beta_norm;

	double alpha_0, alpha_1, beta_0, beta_1;

	double sig0, sig1;

	double det, tr;
	double r0, r1, r2, r3;
	double d00, d01, d02,
		d10, d11, d12;

	double new_sig0, new_sig1;
	double temp;
	double w00, w01, w10, w11;
	double p1, p2, p3, w1, w2, w3;

	double h00, h01, h02, h03, h04, h05,
		h11, h12, h13, h14, h15,
		h22, h23, h24, h25,
		h33, h34, h35,
		h44, h45,
		h55;
	double *position = position_of_mesh.data();

	int nnz = pardiso_ja.size();
	pardiso_a.clear(); pardiso_b.clear();
	pardiso_a.resize(nnz, 0.0);
	pardiso_b.resize(2 * V_N, 0.0);

	for (int i = 0; i < F_N; ++i)
	{
		if (face_in_out[i] > 0)
		{
			area_now = area[i];
		}
		else
		{
			area_now = scaf_weight*area[i];
		}
//		area_now = area[i];
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];


		x0 = position[f0];
		y0 = position[f0 + total_num];

		x1 = position[f1];
		y1 = position[f1 + total_num];

		x2 = position[f2];
		y2 = position[f2 + total_num];

		f0 = id2index[f0];
		f1 = id2index[f1];
		f2 = id2index[f2];
		if (f0 == -1 && f1 == -1 && f2 == -1)
			continue;

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
		beta_0 = j00 - j11; beta_1 = j10 + j01;

		alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1 * alpha_1);
		beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1 * beta_1);

		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;

		new_sig0 = sqrt(1 + 1 / sig0 + 1 / (sig0*sig0) + 1 / (sig0*sig0*sig0)); new_sig1 = sqrt(1 + 1 / sig1 + 1 / (sig1*sig1) + 1 / (sig1*sig1*sig1));

		if (abs(sig1 - sig0) < 1e-8)
		{
			temp = 0;
		}
		else
		{
			temp = (new_sig1 - new_sig0) / (sig1*sig1 - sig0 * sig0);
		}

		w00 = temp * (j00*j00 + j01 * j01 - 0.5*(sig0*sig0 + sig1 * sig1)) + 0.5*(new_sig0 + new_sig1);
		w01 = temp * (j00*j10 + j01 * j11);
		w10 = temp * (j10*j00 + j11 * j01);
		w11 = temp * (j10*j10 + j11 * j11 - 0.5*(sig0*sig0 + sig1 * sig1)) + 0.5*(new_sig0 + new_sig1);

		p1 = p00 * p00 + p01 * p01; p2 = p00 * p10 + p01 * p11; p3 = p10 * p10 + p11 * p11;
		w1 = w00 * w00 + w10 * w10; w2 = w00 * w01 + w10 * w11; w3 = w01 * w01 + w11 * w11;

		//area_now *= 2;

		h00 = area_now * (p1 + p2 + p2 + p3)*w1; h01 = -area_now * (p1 + p2)*w1; h02 = -area_now * (p2 + p3)*w1; h03 = area_now * (p1 + p2 + p2 + p3)*w2; h04 = -area_now * (p1 + p2)*w2; h05 = -area_now * (p2 + p3)*w2;
		h11 = area_now * p1*w1;                  h12 = area_now * p2*w1;    	 h13 = -area_now * (p1 + p2)*w2; h14 = area_now * p1*w2;                  h15 = area_now * p2*w2;
		h22 = area_now * p3*w1;                  h23 = -area_now * (p2 + p3)*w2; h24 = area_now * p2*w2;         h25 = area_now * p3*w2;
		h33 = area_now * (p1 + p2 + p2 + p3)*w3; h34 = -area_now * (p1 + p2)*w3; h35 = -area_now * (p2 + p3)*w3;
		h44 = area_now * p1*w3;                  h45 = area_now * p2*w3;
		h55 = area_now * p3*w3;


		det = j00 * j11 - j01 * j10;
		tr = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);

		d00 = -p00 - p10; d01 = p00; d02 = p10;
		d10 = -p01 - p11; d11 = p01; d12 = p11;

		r0 = area_now * ((1 + 1 / (det*det))*j00 - tr * j11 / (det*det*det));
		r1 = area_now * ((1 + 1 / (det*det))*j01 + tr * j10 / (det*det*det));
		r2 = area_now * ((1 + 1 / (det*det))*j10 + tr * j01 / (det*det*det));
		r3 = area_now * ((1 + 1 / (det*det))*j11 - tr * j00 / (det*det*det));

		if (f0 != -1)
		{
			pardiso_b[f0] -= r0 * d00 + r1 * d10;
			pardiso_b[f0 + V_N] -= r2 * d00 + r3 * d10;
			pardiso_a[id_h00[i]] += h00;
			pardiso_a[id_h33[i]] += h33;
			pardiso_a[id_h03[i]] += h03;
		}
		if (f1 != -1)
		{
			pardiso_b[f1] -= r0 * d01 + r1 * d11;
			pardiso_b[f1 + V_N] -= r2 * d01 + r3 * d11;
			pardiso_a[id_h11[i]] += h11;
			pardiso_a[id_h44[i]] += h44;
			pardiso_a[id_h14[i]] += h14;
		}
		if (f2 != -1)
		{
			pardiso_b[f2] -= r0 * d02 + r1 * d12;
			pardiso_b[f2 + V_N] -= r2 * d02 + r3 * d12;
			pardiso_a[id_h22[i]] += h22;
			pardiso_a[id_h55[i]] += h55;
			pardiso_a[id_h25[i]] += h25;
		}

		if (f0 != -1 && f1 != -1)
		{
			pardiso_a[id_h01[i]] += h01;
			pardiso_a[id_h04[i]] += h04;
			pardiso_a[id_h13[i]] += h13;
			pardiso_a[id_h34[i]] += h34;
		}
		if (f0 != -1 && f2 != -1)
		{
			pardiso_a[id_h02[i]] += h02;
			pardiso_a[id_h05[i]] += h05;
			pardiso_a[id_h23[i]] += h23;
			pardiso_a[id_h35[i]] += h35;
		}
		if (f2 != -1 && f1 != -1)
		{
			pardiso_a[id_h12[i]] += h12;
			pardiso_a[id_h15[i]] += h15;
			pardiso_a[id_h24[i]] += h24;
			pardiso_a[id_h45[i]] += h45;
		}

	}

	//handle term

	for (size_t i = 0; i < handle_ids.size(); i++)
	{
		f0 = handle_ids[i];
		x0 = position[f0];
		y0 = position[f0 + total_num];

		f0 = id2index[f0];

		pardiso_b[f0] -= 2 * handle_weight*(x0 - handle_target_pos(i, 0));
		pardiso_b[f0 + V_N] -= 2 * handle_weight*(y0 - handle_target_pos(i, 1));

		pardiso_a[pardiso_ia[f0]] += 2 * handle_weight;
		pardiso_a[pardiso_ia[f0 + V_N]] += 2 * handle_weight;

	}
	pardiso->a = pardiso_a;
	pardiso->rhs = pardiso_b;

	pardiso->factorize();
	pardiso->pardiso_solver();

	vector<double> result_d = pardiso->result;

	VectorXd negative_grad(2 * total_num), d(2 * total_num);
	negative_grad.setZero();
	d.setZero();

	for (int i = 0; i < V_N; i++)
	{
		negative_grad(var_ids[i]) = pardiso_b[i];
		negative_grad(var_ids[i] + total_num) = pardiso_b[i + V_N];
		d(var_ids[i]) = result_d[i];
		d(var_ids[i] + total_num) = result_d[i + V_N];
	}


	double temp_t;
	max_step(position_of_mesh, d, temp_t);

	double alpha = min(1.0, 0.8 * temp_t);
	backtracking_line_search(position_of_mesh, d, negative_grad, alpha);

	position_of_mesh += alpha * d;
	cout << "slim step length : " << alpha << endl;

	return Energysource();
}

void Parafun::CM()
{
	double area_now;
	int f0, f1, f2;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	double x0, y0, x1, y1, x2, y2;

	double hi_0, hi_1;

	double alpha_0, alpha_1, beta_0, beta_1;

	double s1, s2, sig0, sig1;

	double alpha_norm, beta_norm;
	double h_u, h_v, walpha, wbeta;

	double a1x0, a1x1, a1x2, a1x3, a1x4, a1x5,
		a2x0, a2x1, a2x2, a2x3, a2x4, a2x5;

	double aa, bb;
	double uu, vv, uv;
	double u, v;

	double h00, h01, h02, h03, h04, h05,
		h11, h12, h13, h14, h15,
		h22, h23, h24, h25,
		h33, h34, h35,
		h44, h45,
		h55;
	double *position = position_of_mesh.data();
	int nnz = pardiso_ja.size();
	pardiso_a.clear(); pardiso_b.clear();
	pardiso_a.resize(nnz, 0.0);
	pardiso_b.resize(2 * V_N, 0.0);


	//scaffold computation
	for (int i = 0; i < F_N; i++)
	{
		if (face_in_out[i] > 0)
		{
			area_now = area[i];
		}
		else
		{
			area_now = scaf_weight*area[i];
		}

		area_now = area[i];
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];
		x0 = position[f0];
		y0 = position[f0 + total_num];

		x1 = position[f1];
		y1 = position[f1 + total_num];

		x2 = position[f2];
		y2 = position[f2 + total_num];

		f0 = id2index[f0];
		f1 = id2index[f1];
		f2 = id2index[f2];
		int flag_ = 0;
		if (f0 != -1)flag_++;
		if (f1 != -1)flag_++;
		if (f2 != -1)flag_++;

		if (flag_ == 0)continue;

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
		beta_0 = j00 - j11;  beta_1 = j10 + j01;

		alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1 * alpha_1);
		beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1 * beta_1);

		if (beta_norm<1e-30)
		{
			beta_norm = 1e-15;
		}

		s1 = (p00)*(p00 + p10) + (p01)*(p01 + p11);
		s2 = (p10)*(p00 + p10) + (p11)*(p01 + p11);

		double h1 = p00 * p00 + p01 * p01;
		double h2 = p00 * p10 + p01 * p11;
		double h3 = p10 * p10 + p11 * p11;
		double h4 = p00 * p11 - p01 * p10;

		a1x0 = alpha_0 * (-p00 - p10) + alpha_1 * (p01 + p11);  a1x1 = alpha_0 * p00 - alpha_1 * p01; a1x2 = alpha_0 * p10 - alpha_1 * p11;
		a1x3 = alpha_0 * (-p01 - p11) + alpha_1 * (-p00 - p10); a1x4 = alpha_0 * p01 + alpha_1 * p00; a1x5 = alpha_0 * p11 + alpha_1 * p10;

		a2x0 = beta_0 * (-p00 - p10) + beta_1 * (-p01 - p11);   a2x1 = beta_0 * p00 + beta_1 * p01;   a2x2 = beta_0 * p10 + beta_1 * p11;
		a2x3 = beta_0 * (p01 + p11) + beta_1 * (-p00 - p10);    a2x4 = -beta_0 * p01 + beta_1 * p00;  a2x5 = -beta_0 * p11 + beta_1 * p10;

		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;

		hi_0 = 2 + 6 * 1 / (sig0*sig0*sig0*sig0); hi_1 = 2 + 6 * 1 / (sig1*sig1*sig1*sig1);

		aa = 0.25 / alpha_norm; bb = 0.25 / beta_norm;

		uu = aa * aa*(area_now*hi_0 + area_now * hi_1);
		vv = bb * bb*(area_now*hi_0 + area_now * hi_1);
		uv = aa * bb*(area_now*hi_0 - area_now * hi_1);

		h_u = area_now * (2 * sig0 - 2 * 1 / (sig0*sig0*sig0));
		h_v = area_now * (2 * sig1 - 2 * 1 / (sig1*sig1*sig1));

		walpha = h_u + h_v;
		wbeta = h_u - h_v;

		double hwa1 = (walpha * 0.25 / alpha_norm); double hwa2 = -(walpha * 0.25*0.25 / (alpha_norm*alpha_norm*alpha_norm));
		double hwb1 = (wbeta * 0.25 / beta_norm); double hwb2 = -(wbeta *0.25*0.25 / (beta_norm*beta_norm*beta_norm));


		h00 = uu * a1x0*a1x0 + vv * a2x0*a2x0 + uv * a1x0*a2x0 + uv * a2x0*a1x0; h01 = uu * a1x0*a1x1 + vv * a2x0*a2x1 + uv * a1x0*a2x1 + uv * a2x0*a1x1; h02 = uu * a1x0*a1x2 + vv * a2x0*a2x2 + uv * a1x0*a2x2 + uv * a2x0*a1x2; h03 = uu * a1x0*a1x3 + vv * a2x0*a2x3 + uv * a1x0*a2x3 + uv * a2x0*a1x3; h04 = uu * a1x0*a1x4 + vv * a2x0*a2x4 + uv * a1x0*a2x4 + uv * a2x0*a1x4; h05 = uu * a1x0*a1x5 + vv * a2x0*a2x5 + uv * a1x0*a2x5 + uv * a2x0*a1x5;

		h11 = uu * a1x1*a1x1 + vv * a2x1*a2x1 + uv * a1x1*a2x1 + uv * a2x1*a1x1; h12 = uu * a1x1*a1x2 + vv * a2x1*a2x2 + uv * a1x1*a2x2 + uv * a2x1*a1x2; h13 = uu * a1x1*a1x3 + vv * a2x1*a2x3 + uv * a1x1*a2x3 + uv * a2x1*a1x3; h14 = uu * a1x1*a1x4 + vv * a2x1*a2x4 + uv * a1x1*a2x4 + uv * a2x1*a1x4; h15 = uu * a1x1*a1x5 + vv * a2x1*a2x5 + uv * a1x1*a2x5 + uv * a2x1*a1x5;

		h22 = uu * a1x2*a1x2 + vv * a2x2*a2x2 + uv * a1x2*a2x2 + uv * a2x2*a1x2; h23 = uu * a1x2*a1x3 + vv * a2x2*a2x3 + uv * a1x2*a2x3 + uv * a2x2*a1x3; h24 = uu * a1x2*a1x4 + vv * a2x2*a2x4 + uv * a1x2*a2x4 + uv * a2x2*a1x4; h25 = uu * a1x2*a1x5 + vv * a2x2*a2x5 + uv * a1x2*a2x5 + uv * a2x2*a1x5;

		h33 = uu * a1x3*a1x3 + vv * a2x3*a2x3 + uv * a1x3*a2x3 + uv * a2x3*a1x3; h34 = uu * a1x3*a1x4 + vv * a2x3*a2x4 + uv * a1x3*a2x4 + uv * a2x3*a1x4; h35 = uu * a1x3*a1x5 + vv * a2x3*a2x5 + uv * a1x3*a2x5 + uv * a2x3*a1x5;

		h44 = uu * a1x4*a1x4 + vv * a2x4*a2x4 + uv * a1x4*a2x4 + uv * a2x4*a1x4; h45 = uu * a1x4*a1x5 + vv * a2x4*a2x5 + uv * a1x4*a2x5 + uv * a2x4*a1x5;

		h55 = uu * a1x5*a1x5 + vv * a2x5*a2x5 + uv * a1x5*a2x5 + uv * a2x5*a1x5;

		if (walpha >= 0)
		{
			h00 += hwa1 * (s1 + s2) + hwa2 * a1x0*a1x0; h01 += hwa1 * (-s1) + hwa2 * a1x0*a1x1; h02 += hwa1 * (-s2) + hwa2 * a1x0*a1x2; h03 += hwa2 * a1x0*a1x3; h04 += hwa1 * (h4)+hwa2 * a1x0*a1x4; h05 += hwa1 * (-h4) + hwa2 * a1x0*a1x5;
			h11 += hwa1 * (h1)+hwa2 * a1x1*a1x1;        h12 += hwa1 * (h2)+hwa2 * a1x1*a1x2;    h13 += hwa1 * (-h4) + hwa2 * a1x1*a1x3; h14 += hwa2 * a1x1*a1x4; h15 += hwa1 * (h4)+hwa2 * a1x1*a1x5;
			h22 += hwa1 * (h3)+hwa2 * a1x2*a1x2;        h23 += hwa1 * (h4)+hwa2 * a1x2*a1x3;    h24 += hwa1 * (-h4) + hwa2 * a1x2*a1x4; h25 += hwa2 * a1x2*a1x5;
			h33 += hwa1 * (s1 + s2) + hwa2 * a1x3*a1x3; h34 += hwa1 * (-s1) + hwa2 * a1x3*a1x4; h35 += hwa1 * (-s2) + hwa2 * a1x3*a1x5;
			h44 += hwa1 * (h1)+hwa2 * a1x4*a1x4;        h45 += hwa1 * (h2)+hwa2 * a1x4*a1x5;
			h55 += hwa1 * (h3)+hwa2 * a1x5*a1x5;
		}
		h00 += hwb1 * (s1 + s2) + hwb2 * a2x0*a2x0; h01 += hwb1 * (-s1) + hwb2 * a2x0*a2x1; h02 += hwb1 * (-s2) + hwb2 * a2x0*a2x2; h03 += hwb2 * a2x0*a2x3; h04 += hwb1 * (-h4) + hwb2 * a2x0*a2x4; h05 += hwb1 * (h4)+hwb2 * a2x0*a2x5;
		h11 += hwb1 * (h1)+hwb2 * a2x1*a2x1;        h12 += hwb1 * (h2)+hwb2 * a2x1*a2x2;    h13 += hwb1 * (h4)+hwb2 * a2x1*a2x3;    h14 += hwb2 * a2x1*a2x4; h15 += hwb1 * (-h4) + hwb2 * a2x1*a2x5;
		h22 += hwb1 * (h3)+hwb2 * a2x2*a2x2;        h23 += hwb1 * (-h4) + hwb2 * a2x2*a2x3; h24 += hwb1 * (h4)+hwb2 * a2x2*a2x4;    h25 += hwb2 * a2x2*a2x5;
		h33 += hwb1 * (s1 + s2) + hwb2 * a2x3*a2x3; h34 += hwb1 * (-s1) + hwb2 * a2x3*a2x4; h35 += hwb1 * (-s2) + hwb2 * a2x3*a2x5;
		h44 += hwb1 * (h1)+hwb2 * a2x4*a2x4;        h45 += hwb1 * (h2)+hwb2 * a2x4*a2x5;
		h55 += hwb1 * (h3)+hwb2 * a2x5*a2x5;

		u = aa * walpha; v = bb * wbeta;

		if (f0 != -1)
		{
			pardiso_b[f0] -= (u*a1x0 + v * a2x0);
			pardiso_b[f0 + V_N] -= (u*a1x3 + v * a2x3);
			pardiso_a[id_h00[i]] += h00;
			pardiso_a[id_h33[i]] += h33;
			pardiso_a[id_h03[i]] += h03;
		}
		if (f1 != -1)
		{
			pardiso_b[f1] -= (u*a1x1 + v * a2x1);
			pardiso_b[f1 + V_N] -= (u*a1x4 + v * a2x4);
			pardiso_a[id_h11[i]] += h11;
			pardiso_a[id_h44[i]] += h44;
			pardiso_a[id_h14[i]] += h14;
		}
		if (f2 != -1)
		{
			pardiso_b[f2] -= (u*a1x2 + v * a2x2);
			pardiso_b[f2 + V_N] -= (u*a1x5 + v * a2x5);
			pardiso_a[id_h22[i]] += h22;
			pardiso_a[id_h55[i]] += h55;
			pardiso_a[id_h25[i]] += h25;
		}

		if (flag_ == 3)
		{
			pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02;  pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
			pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h15[i]] += h15;
			pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24;
			pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
			pardiso_a[id_h45[i]] += h45;
		}
		else if (flag_ == 2)
		{
			if (f0 == -1)
			{
				pardiso_a[id_h12[i]] += h12;
				pardiso_a[id_h15[i]] += h15;
				pardiso_a[id_h24[i]] += h24;
				pardiso_a[id_h45[i]] += h45;
			}
			else if (f1 == -1)
			{
				pardiso_a[id_h02[i]] += h02;
				pardiso_a[id_h05[i]] += h05;
				pardiso_a[id_h23[i]] += h23;
				pardiso_a[id_h35[i]] += h35;
			}
			else
			{
				pardiso_a[id_h01[i]] += h01;
				pardiso_a[id_h04[i]] += h04;
				pardiso_a[id_h13[i]] += h13;
				pardiso_a[id_h34[i]] += h34;
			}
		}
	}

	//handle term

	for (size_t i = 0; i < handle_ids.size(); i++)
	{
		f0 = handle_ids[i];
		x0 = position[f0];
		y0 = position[f0 + total_num];

		f0 = id2index[f0];

		pardiso_b[f0] -= 2 * handle_weight*(x0 - handle_target_pos(i, 0));
		pardiso_b[f0 + V_N] -= 2 * handle_weight*(y0 - handle_target_pos(i, 1));

		pardiso_a[pardiso_ia[f0]] += 2 * handle_weight;
		pardiso_a[pardiso_ia[f0+V_N]] += 2 * handle_weight;

	}


	pardiso->a = pardiso_a;
	pardiso->rhs = pardiso_b;

	pardiso->factorize();
	pardiso->pardiso_solver();

	vector<double> result_d = pardiso->result;

	VectorXd negative_grad(2 * total_num), d(2 * total_num);
	negative_grad.setZero();
	d.setZero();

	for (int i = 0; i < V_N; i++)
	{
		negative_grad(var_ids[i]) = pardiso_b[i];
		negative_grad(var_ids[i] + total_num) = pardiso_b[i + V_N];
		d(var_ids[i]) = result_d[i];
		d(var_ids[i] + total_num) = result_d[i + V_N];
	}

	pardiso->free_numerical_factorization_memory();

	double temp_t;
	max_step(position_of_mesh, d, temp_t);

	double alpha = 0.95 * temp_t;

	backtracking_line_search(position_of_mesh, d, negative_grad, alpha);

	position_of_mesh += alpha * d;
//	cout << "cm step length: " << alpha << endl;
	Energysource();
}

void Parafun::id_vs_index()
{
	auto& bnd_ids = fixed_ids;
	const auto bnd_n = bnd_ids.size();

	if (bnd_n == 0)
	{
		var_ids.clear();
		var_ids.reserve(total_num);
		for (int i = 0; i < total_num; i++)
		{
			var_ids.push_back(i);
		}
	}
	else
	{
		var_ids.clear();
		var_ids.reserve(total_num);
		int assign = 0, i = 0;
		for (int get = 0; i < total_num && get < bnd_ids.size(); i++)
		{
			if (bnd_ids(get) == i)
				get++;
			else
				var_ids.push_back(i);
		}
		while (i < total_num)
		{
			var_ids.push_back(i);
			i++;
		}
	}

	id2index.resize(total_num, -1);
	for (size_t i = 0; i < var_ids.size(); i++)
	{
		id2index[var_ids[i]] = i;
	}

}

void Parafun::max_step(const VectorXd &xx, const VectorXd &dd, double &step)
{
	double temp_t = numeric_limits<double>::infinity();
	int f0, f1, f2;
	double a, b, c, b1, b2, tt, tt1, tt2;
	double x0, x1, x2, x3, x4, x5, d0, d1, d2, d3, d4, d5;
	const double *x = xx.data();
	const double *d = dd.data();
	for (int i = 0; i < F_N; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = x[f0]; x1 = x[f1]; x2 = x[f2]; x3 = x[f0 + total_num]; x4 = x[f1 + total_num]; x5 = x[f2 + total_num];
		d0 = d[f0]; d1 = d[f1]; d2 = d[f2]; d3 = d[f0 + total_num]; d4 = d[f1 + total_num]; d5 = d[f2 + total_num];

		a = (d1 - d0) * (d5 - d3) - (d4 - d3) * (d2 - d0);
		b1 = (d1 - d0) * (x5 - x3) + (x1 - x0) * (d5 - d3);
		b2 = (x4 - x3) * (d2 - d0) + (x2 - x0) * (d4 - d3);
		b = b1 - b2;
		c = (x1 - x0) * (x5 - x3) - (x4 - x3) * (x2 - x0);
		tt = get_smallest_pos_quad_zero(a, b, c);

		if (temp_t > tt)
		{
			temp_t = tt;
		}

	}
	if (temp_t == INFINITY)
	{
		temp_t = 100;
	}
	step = temp_t;
}

double get_smallest_pos_quad_zero(double a, double b, double c)
{
	using namespace std;
	double t1, t2;
	if (std::abs(a) > 1.0e-10)
	{
		double delta_in = pow(b, 2) - 4 * a * c;
		if (delta_in <= 0)
		{
			return INFINITY;
		}

		double delta = sqrt(delta_in); // delta >= 0
		if (b >= 0) // avoid subtracting two similar numbers
		{
			double bd = -b - delta;
			t1 = 2 * c / bd;
			t2 = bd / (2 * a);
		}
		else
		{
			double bd = -b + delta;
			t1 = bd / (2 * a);
			t2 = (2 * c) / bd;
		}

		assert(std::isfinite(t1));
		assert(std::isfinite(t2));

		if (a < 0) std::swap(t1, t2); // make t1 > t2
		// return the smaller positive root if it exists, otherwise return infinity
		if (t1 > 0)
		{
			return t2 > 0 ? t2 : t1;
		}
		else
		{
			return INFINITY;
		}
	}
	else
	{
		if (b == 0) return INFINITY; // just to avoid divide-by-zero
		t1 = -c / b;
		return t1 > 0 ? t1 : INFINITY;
	}
}

//void Parafun::calc_gradient_norm(const VectorXd &x)
//{
//	double area_now;
//	int f0, f1, f2;
//	double j00, j01, j10, j11;
//	double p00, p01, p10, p11;
//	double q00, q01, q10, q11;
//
//	double x0, y0, x1, y1, x2, y2;
//
//	double det, tr;
//	double r0, r1, r2, r3;
//	double d00, d01, d02,
//		d10, d11, d12;
//
//	negative_grad_norm.setZero();
//
//	const double *position = x.data();
//
//	for (int i = 0; i < F_N; ++i)
//	{
//		area_now = area[i];
//		f0 = F0[i];
//		f1 = F1[i];
//		f2 = F2[i];
//
//		x0 = position[f0];
//		y0 = position[f0 + V_N];
//
//		x1 = position[f1];
//		y1 = position[f1 + V_N];
//
//		x2 = position[f2];
//		y2 = position[f2 + V_N];
//
//		q00 = x1 - x0; q01 = x2 - x0;
//		q10 = y1 - y0; q11 = y2 - y0;
//
//		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];
//
//		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;
//
//		det = j00*j11 - j01*j10;
//		tr = (j00*j00 + j01*j01 + j10*j10 + j11*j11);
//
//		d00 = -p00 - p10; d01 = p00; d02 = p10;
//		d10 = -p01 - p11; d11 = p01; d12 = p11;
//
//		r0 = 2 * area_now * ((1 + 1 / (det*det))*j00 - tr*j11 / (det*det*det));
//		r1 = 2 * area_now * ((1 + 1 / (det*det))*j01 + tr*j10 / (det*det*det));
//		r2 = 2 * area_now * ((1 + 1 / (det*det))*j10 + tr*j01 / (det*det*det));
//		r3 = 2 * area_now * ((1 + 1 / (det*det))*j11 - tr*j00 / (det*det*det));
//
//		negative_grad_norm(f0) -= r0*d00 + r1*d10;
//		negative_grad_norm(f1) -= r0*d01 + r1*d11;
//		negative_grad_norm(f2) -= r0*d02 + r1*d12;
//		negative_grad_norm(f0 + V_N) -= r2*d00 + r3*d10;
//		negative_grad_norm(f1 + V_N) -= r2*d01 + r3*d11;
//		negative_grad_norm(f2 + V_N) -= r2*d02 + r3*d12;
//	}
//
//	g_norm = negative_grad_norm.norm();
//}

void Parafun::backtracking_line_search(const VectorXd &x, const VectorXd &d, const VectorXd &negetive_grad, double &alpha)
{
	double h = 0.5;
	double tt = -(negetive_grad.transpose()*d)(0, 0);
	double c = 0.2; 
	double ex;
	Energy(x, ex);
	double e;
	VectorXd x_new = x + alpha * d;
	Energy(x_new, e);
	while (e > ex + alpha * c * tt)
	{
		alpha = h*alpha;
		x_new = x + alpha * d;
		Energy(x_new, e);
	}
}

void Parafun::Energy(const VectorXd &position, double &energyupdate)
{
	double energy = 0;


	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_d;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;
	const double *pos = position.data();
	for (int i = 0; i < F_N; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = pos[f0];
		y0 = pos[f0 + total_num];

		x1 = pos[f1];
		y1 = pos[f1 + total_num];

		x2 = pos[f2];
		y2 = pos[f2 + total_num];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = update_p00[i]; p01 = update_p01[i]; p10 = update_p10[i]; p11 = update_p11[i];

		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;


		det = j00*j11 - j01*j10;
		E_d = (1 + 1 / (det*det)) * (j00*j00 + j01*j01 + j10*j10 + j11*j11);

		energy += area[i] * E_d;
	}
	for (size_t i = 0; i < handle_ids.size(); i++)
	{
		f0 = handle_ids[i];
		x0 = pos[f0];
		y0 = pos[f0 + total_num];
		energy += handle_weight*((x0 - handle_target_pos(i, 0))*(x0 - handle_target_pos(i, 0)) + (y0 - handle_target_pos(i, 1))*(y0 - handle_target_pos(i, 1)));
	}
	energyupdate = energy;
}

bool Parafun::Energysource()
{
	double end_e_one_temp = 0, end_e_area = 0;

	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_1, E_2;

	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	const double *pos = position_of_mesh.data();
	for (int i = 0; i < F_N; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = pos[f0];
		y0 = pos[f0 + total_num];

		x1 = pos[f1];
		y1 = pos[f1 + total_num];

		x2 = pos[f2];
		y2 = pos[f2 + total_num];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;

		det = j00*j11 - j01*j10;

		E_1 = (j00*j00 + j01*j01 + j10*j10 + j11*j11);
		E_2 = 1.0 / (det*det)* E_1;

		end_e_one_temp += E_1;
		end_e_one_temp += E_2;
		end_e_area += ((E_1 + E_2)*area_src[i]);
	}
	energy_uniform = end_e_one_temp /F_N;
	for (size_t i = 0; i < handle_ids.size(); i++)
	{
		f0 = handle_ids[i];
		x0 = pos[f0];
		y0 = pos[f0 + total_num];
		double handle_energy = handle_weight*((x0 - handle_target_pos(i, 0))*(x0 - handle_target_pos(i, 0)) + (y0 - handle_target_pos(i, 1))*(y0 - handle_target_pos(i, 1)));
		cout << "handle_e: " << handle_energy << endl;
	}
	energy_area = end_e_area;
//	cout << "meshdis: " << energy_area / originmesh_area_sqrt / originmesh_area_sqrt << endl;

	if (!isfinite(energy_area))
	{
		return false;
	}
	return true;


}

void Parafun::local_coordinate_inverse(int i, double &p00, double &p01, double &p10, double &p11)
{
	int f0 = F0[i];
	int f1 = F1[i];
	int f2 = F2[i];

	OpenMesh::Vec3d x_ = (mesh.point(mesh.vertex_handle(f1)) - mesh.point(mesh.vertex_handle(f0)));
	double x1_0 = x_.length();
	OpenMesh::Vec3d l_ = mesh.point(mesh.vertex_handle(f2)) - mesh.point(mesh.vertex_handle(f0));
	OpenMesh::Vec3d y_ = mesh.normal(mesh.face_handle(i)) % (1 / x1_0 * x_);
	double x2_0 = 1 / x1_0*l_ | x_;
	double y2_0 = l_ | y_;

	p00 = 1 / x1_0;
	p01 = -x2_0 / (x1_0*y2_0);
	p10 = 0;
	p11 = 1 / y2_0;
}

double Parafun::newton_equation(const double & a, const double & b, const double & K)
{
	double tt = 1;
	double E_d = pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K;
	while (abs(E_d) > 1e-5)
	{
		tt = tt - 1 / (2 * log(a)*pow(a, 2 * tt) + 2 * log(b)* pow(b, 2 * tt) + 2 * log(1 / a)* pow(1 / a, 2 * tt) + 2 * log(1 / b)*pow(1 / b, 2 * tt))*(pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K);
		E_d = pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K;
	}
	return tt;
}
