#include "dchart/Parafun_bij.h"

void Parafun_bij::after_mesh_improve()
{

	int mv_n = d_.mv_num;
	int mf_n = d_.mf_num;
	int sv_n = d_.sv_num;
	int sf_n = d_.sf_num;

	int v_n = mv_n + sv_n;
	int f_n = mf_n + sf_n;

	var_ids.clear();


	total_num = d_.v_num;
	F_N = d_.f_num;
	V_N = d_.v_num-d_.frame_ids.size();
	F0.resize(F_N);
	F1.resize(F_N);
	F2.resize(F_N);

	position_of_mesh.resize(2 * total_num);

	for (size_t i = 0; i < dim; i++)
	{
		position_of_mesh.block(i*total_num, 0, total_num, 1) = d_.w_uv.col(i);
	}
	//{
	//	ofstream of_ycy;
	//	of_ycy.open("uv_now.txt", ios::trunc);
	//	for (size_t i = 0; i < 2 * total_num; i++)
	//	{
	//		of_ycy <<position_of_mesh(i) << endl;
	//	}
	//	of_ycy.close();
	//}

	if (pardiso)
	{
		delete pardiso;
		pardiso = NULL;
	}
	init();
}

Parafun_bij::~Parafun_bij()
{
}

//void Parafun::BPE()
//{
//	if (pardiso != NULL)
//	{
//		delete pardiso;
//		pardiso = NULL;
//	}
//	pardiso = new PardisoSolver();
//	pardiso->ia = pardiso_ia;
//	pardiso->ja = pardiso_ja;
//	pardiso->a.resize(pardiso_ja.size());
//	pardiso->nnz = pardiso_ja.size();
//	pardiso->num = 2 * V_N;
//
//	pardiso->pardiso_init();
//
//	vector<double> energy_area_process;
//	energy_area_process.reserve(MAX_ITER_NUM);
//	Energysource();
//
//	energy_area_process.push_back(energy_area);
//	double energy_pre = 0;
//	double energy_cur = energy_uniform;
//	
//	int iter_num_cur = 0;
//	Intp_T_Min = 0;
//	changetocm_flag = 0;
//
//	int slim_iter_num = 0;
//	int sum_iter_num = 0;
//
//	double conv_percent = 1;
//
//	g_norm=1.0;
//
//	long time_beg, time_end;
//	time_beg = clock();
//
//	while (iter_num_cur < MAX_ITER_NUM)
//	{
//		iter_num_cur++;
//		energy_pre = energy_cur;
//		Update_source_same_t();
//		if (changetocm_flag < 0.99&&conv_percent>0.1&&Intp_T_Min<0.999)
//		{
//			SLIM();
//			energy_area_process.push_back(energy_area);
//			slim_iter_num++;
//			sum_iter_num++;
//
//			energy_cur = energy_uniform;
//			conv_percent = abs(energy_cur - energy_pre) / energy_pre;
//			calc_gradient_norm(position_of_mesh);
//			if (conv_percent <= convgence_con_rate||g_norm<= convgence_con_rate)
//			{
//				break;
//			}
//		}
//		else
//		{
//			break;
//		}
//	}
//
//	int cm_iter_num = 0;
//	while (iter_num_cur < MAX_ITER_NUM)
//	{
//		iter_num_cur++;
//		energy_pre = energy_cur;
//		Update_source_same_t();
//		if (conv_percent > 0.01&&Intp_T_Min<0.999)
//		{
//			CM();
//			energy_area_process.push_back(energy_area);
//			cm_iter_num++;
//			sum_iter_num++;
//
//			energy_cur = energy_uniform;
//			conv_percent = abs(energy_cur - energy_pre) / energy_pre;
//			calc_gradient_norm(position_of_mesh);
//			if (conv_percent <= convgence_con_rate || g_norm <= convgence_con_rate)
//			{
//				break;
//			}
//		}
//		else
//		{
//			recover_to_src();
//			energy_cur = energy_area;
//			while (iter_num_cur < MAX_ITER_NUM)
//			{
//				iter_num_cur++;
//				energy_pre = energy_cur;
//				
//				CM();
//				energy_area_process.push_back(energy_area);
//				sum_iter_num++;
//
//				energy_cur = energy_area;
//				conv_percent = abs(energy_cur - energy_pre) / energy_pre;
//				calc_gradient_norm(position_of_mesh);
//				if (conv_percent <= convgence_con_rate || g_norm <= convgence_con_rate)
//				{
//					break;
//				}
//			}
//			break;
//		}
//
//	}
//
//	time_end = clock();
//	time_consumption = (time_end - time_beg) / 1000.0;
//	
//	ofstream of_energy;
//	string energy_str = modelname + "_comp_energy_result.txt";
//	cout << "COMP ====== time_consumption: " << time_consumption << " s;slim_iter: " << slim_iter_num << "; cm_iter: " << cm_iter_num << "; sum_iter: " << sum_iter_num << endl;
//
//	of_energy.open(energy_str, ios::trunc);
//	for (size_t i = 0; i < energy_area_process.size(); i++)
//	{
//		of_energy <<fixed<<setprecision(8)<<energy_area_process[i]<< endl;
//	}
//
//	of_energy.close();
//
//	delete pardiso;
//	pardiso = NULL;
//}

void Parafun_bij::init()
{
	for (int i = 0; i < F_N; ++i)
	{
		F0[i] = d_.surface_F(i, 0);
		F1[i] = d_.surface_F(i, 1);
		F2[i] = d_.surface_F(i, 2);
	}
	//var_ids ??????????????????index???????? ?????????????????? id??
	// id2index ????????????????id ??????????????????index???????????????????????? -1??
	id_vs_index();
	handle_mintri();
	Pre_calculate();
	init_area();
}

void Parafun_bij::init_area()
{
	area.resize(F_N);
	int src_t_num = d_.m_T.rows();
	area_scaf.resize(F_N- src_t_num);
	area_src.resize(src_t_num);
	for (int i = 0; i < src_t_num; ++i)
	{
		area_src[i] = d_.m_M(i);
		area[i] = d_.m_M(i);
	}
	for (int i = src_t_num; i < F_N; ++i)
	{
		area_scaf[i- src_t_num] = d_.s_M(i- src_t_num);
		area[i] = d_.s_M(i - src_t_num);
	}
}

void Parafun_bij::run_bpe()
{
	init();
	//BPE();
}

double Parafun_bij::perform_iteration_cm(bool use_CM)
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
	long time_beg, time_end;
	time_beg = clock();
	if(use_CM)
		CM();
	else
		SLIM();

	//for (size_t i = 0; i < dim; i++)
	//{
	//	d_.w_uv.col(i) = position_of_mesh.block(i*total_num, 0, total_num, 1);
	//}
	d_.w_uv = Map<Matrix<double, -1, -1, Eigen::ColMajor>>(position_of_mesh.data(), total_num, dim);

	time_end = clock();
	time_consumption = (time_end - time_beg) / 1000.0;

	//cout << "PP ====== time_consumption: " << time_consumption << " s; energy: " << energy_area << endl;

	delete pardiso;
	pardiso = NULL;
	return energy_area/d_.mesh_measure;
}


void Parafun_bij::id_vs_index()
{
	auto& bnd_ids = d_.frame_ids;
	const auto bnd_n = bnd_ids.size();
	fix_pos.resize(bnd_n, dim);
	for (size_t i = 0; i < bnd_n; i++)
	{
		fix_pos.row(i) = d_.w_uv.row(bnd_ids(i));
	}

	{ 
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

//void Parafun::Tutte()
//{
//	int boundary_num = 0;
//	auto it1 = mesh.halfedges_begin();
//	while (!mesh.is_boundary(*it1))
//		it1++;
//	auto he_start = *it1;
//	auto he_it = he_start;
//	do
//	{
//		he_it = mesh.next_halfedge_handle(he_it);
//		boundary_num++;
//	} while (he_it != he_start);
//
//	double delta_angle = 2 * M_PI / boundary_num;
//	double area_1_factor = sqrt(1.0 / M_PI);
//
//	position_of_mesh.resize(2 * V_N);
//	for (int i = 0; i < boundary_num; ++i)
//	{
//		auto v_h = mesh.to_vertex_handle(he_start);
//		position_of_mesh(v_h.idx()) = area_1_factor*cos(i * delta_angle);
//		position_of_mesh(v_h.idx() + V_N) = area_1_factor*sin(-i * delta_angle);		
//		he_start = mesh.next_halfedge_handle(he_start);
//	}
//
//	vector<int> pardiso_it;
//	vector<int> pardiso_jt;
//	vector<double> pardiso_t;
//	vector<double> pardiso_tu;
//	vector<double> pardiso_tv;
//
//	pardiso_it.reserve(V_N+1);
//	pardiso_jt.reserve(6*V_N);
//	pardiso_t.reserve(6*V_N);
//	pardiso_tu.resize(V_N, 0.0);
//	pardiso_tv.resize(V_N, 0.0);
//	for (size_t i = 0; i < V_N; i++)
//	{
//		pardiso_it.push_back(pardiso_jt.size());
//
//		auto v_h = mesh.vertex_handle(i);
//		if (mesh.is_boundary(v_h))
//		{
//			pardiso_jt.push_back(i);
//			pardiso_t.push_back(1);
//
//			pardiso_tu[i]=position_of_mesh(i);
//			pardiso_tv[i]=position_of_mesh(i + V_N);
//
//		}
//		else
//		{
//			pardiso_jt.push_back(i);
//			pardiso_t.push_back(mesh.valence(v_h));
//			vector<int> row_id;
//			row_id.reserve(mesh.valence(v_h));
//			double bu = 0.0; double bv = 0.0;
//			for (auto it2 = mesh.vv_begin(v_h); it2 != mesh.vv_end(v_h); ++it2)
//			{
//				int vv_id = it2->idx();
//				if (mesh.is_boundary(*it2))
//				{
//					bu += position_of_mesh(vv_id);
//					bv += position_of_mesh(vv_id+V_N);
//				}
//				else
//				{
//					if (vv_id>i)
//					{
//						row_id.push_back(vv_id);
//					}
//				}
//			}
//			sort(row_id.begin(), row_id.end(), less<int>());
//			for (size_t j = 0; j < row_id.size(); j++)
//			{
//				pardiso_jt.push_back(row_id[j]);
//				pardiso_t.push_back(-1);
//			}
//			pardiso_tu[i] = bu;
//			pardiso_tv[i] = bv;
//		}
//	}
//	pardiso_it.push_back(pardiso_jt.size());
//
//	if (pardiso != NULL)
//	{
//		delete pardiso;
//		pardiso = NULL;
//	}
//	pardiso = new PardisoSolver();
//	pardiso->ia = pardiso_it;
//	pardiso->ja = pardiso_jt;
//	pardiso->nnz = pardiso_jt.size();
//	pardiso->num = V_N;
//
//	pardiso->pardiso_init();
//
//	pardiso->a = pardiso_t;
//
//	pardiso->rhs = pardiso_tu;
//	pardiso->factorize();
//
//
//	pardiso->pardiso_solver();
//
//	for (size_t i = 0; i < V_N; i++)
//	{
//		position_of_mesh(i) = (pardiso->result)[i];
//	}
//
//	pardiso->rhs = pardiso_tv;
//	pardiso->pardiso_solver();
//	for (size_t i = 0; i < V_N; i++)
//	{
//		position_of_mesh(i+V_N) = (pardiso->result)[i];
//	}
//
//	delete pardiso;
//	pardiso = NULL;
//}

void Parafun_bij::Pre_calculate()
{
	source_p00.resize(F_N);
	source_p01.resize(F_N);
	source_p10.resize(F_N);
	source_p11.resize(F_N);

	for (int i = 0; i < d_.m_T.rows(); ++i)
	{
		double p00, p01, p10, p11;
		local_coordinate_inverse(i, p00, p01, p10, p11);

		source_p00[i] = p00;
		source_p01[i] = p01;
		source_p10[i] = p10;
		source_p11[i] = p11;
	}
	for (int i = d_.m_T.rows(); i < F_N; ++i)
	{
		double p00, p01, p10, p11;
		local_coordinate_inverse_scaf(i, p00, p01, p10, p11);

		source_p00[i] = p00;
		source_p01[i] = p01;
		source_p10[i] = p10;
		source_p11[i] = p11;
	}

	pardiso_ia.clear(); pardiso_ia.reserve(2 * V_N + 1);
	pardiso_ja.clear(); pardiso_ja.reserve(8 * V_N);

	typedef Triplet<int> T;
	std::vector<T> tripletlist;

	std::vector<std::set<int>> VV_tmp;
	VV_tmp.resize(V_N);
	for (size_t i = 0; i < d_.m_T.rows(); i++)
	{
		int vid[3];

		for (size_t j = 0; j < d_.m_T.cols(); j++)
		{
			vid[j] = id2index[d_.m_T(i, j)];
		}
		VV_tmp[vid[0]].insert(vid[1]);
		VV_tmp[vid[0]].insert(vid[2]);

		VV_tmp[vid[1]].insert(vid[0]);
		VV_tmp[vid[1]].insert(vid[2]);

		VV_tmp[vid[2]].insert(vid[0]);
		VV_tmp[vid[2]].insert(vid[1]);
	}


	vector<int> s_vid;
	for (size_t i = 0; i < d_.s_T.rows(); i++)
	{
		s_vid.clear();
		for (size_t j = 0; j < d_.s_T.cols(); j++)
		{
			int s_id = id2index[d_.s_T(i, j)];
			if (s_id !=-1)
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
			pardiso_ja.push_back(row_id[k]+V_N);
			tripletlist.push_back(T(i, row_id[k]+V_N, dd));
			++dd;
		}
	}
	for (int i = V_N; i < 2*V_N; i++)
	{
		pardiso_ia.push_back(pardiso_ja.size());
		vector<int> row_id;
		for (auto&var : VV_tmp[i-V_N])
		{
			row_id.push_back(var);
		}
		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i-V_N);

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
	
	for (int i = 0; i < d_.m_T.rows(); i++)
	{
		int f0 = id2index[F0[i]]; int f1 = id2index[F1[i]]; int f2 = id2index[F2[i]]; int f3 = f0 + V_N; int f4 = f1 + V_N; int f5 = f2 + V_N;
		if (f0 == -1 || f1 == -1 || f2 == -1)
		{
			cout << "fuckfuckfuck!!!!" << endl;
		}

		int min01 = min(f0, f1); int max01 = f0 + f1 - min01;
		int min02 = min(f0, f2); int max02 = f0 + f2 - min02;
		int min12 = min(f1, f2); int max12 = f1 + f2 - min12;

		id_h00[i] = pardiso_ia[f0]; id_h01[i] = pardiso_ia[min01] + find_id_in_rows.coeff(min01, max01); id_h02[i] = pardiso_ia[min02] + find_id_in_rows.coeff(min02, max02);
		id_h03[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f3); id_h04[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f4); id_h05[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f5);

		id_h11[i] = pardiso_ia[f1]; id_h12[i] = pardiso_ia[min12] + find_id_in_rows.coeff(min12, max12);
		id_h13[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f3); id_h14[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f4); id_h15[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f5);

		id_h22[i] = pardiso_ia[f2];
		id_h23[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f3); id_h24[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f4); id_h25[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f5);

		id_h33[i] = pardiso_ia[f3]; id_h34[i] = pardiso_ia[min01 + V_N] + find_id_in_rows.coeff(min01 + V_N, max01 + V_N); id_h35[i] = pardiso_ia[min02 + V_N] + find_id_in_rows.coeff(min02 + V_N, max02 + V_N);

		id_h44[i] = pardiso_ia[f4]; id_h45[i] = pardiso_ia[min12 + V_N] + find_id_in_rows.coeff(min12 + V_N, max12 + V_N);

		id_h55[i] = pardiso_ia[f5];

	}

	for (int i = d_.m_T.rows(); i < F_N; i++)
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

void Parafun_bij::SLIM()
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
	int src_t_num = d_.m_T.rows();
	for (int i = 0; i < src_t_num; ++i)
	{
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

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;

		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
		beta_0 = j00 - j11; beta_1 = j10 + j01;

		alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1*alpha_1);
		beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1*beta_1);

		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;

		new_sig0 = sqrt(1 + 1 / sig0 + 1 / (sig0*sig0) + 1 / (sig0*sig0*sig0)); new_sig1 = sqrt(1 + 1 / sig1 + 1 / (sig1*sig1) + 1 / (sig1*sig1*sig1));

		if (beta_norm < 1e-6)
		{
			temp = 0;
		}
		else
		{
			temp = (new_sig1 - new_sig0) / (sig1*sig1 - sig0 * sig0);
		}

		w00 = temp*(j00*j00 + j01*j01 - 0.5*(sig0*sig0 + sig1*sig1)) + 0.5*(new_sig0 + new_sig1);
		w01 = temp*(j00*j10 + j01*j11);
		w10 = temp*(j10*j00 + j11*j01);
		w11 = temp*(j10*j10 + j11*j11 - 0.5*(sig0*sig0 + sig1*sig1)) + 0.5*(new_sig0 + new_sig1);

		p1 = p00*p00 + p01*p01; p2 = p00*p10 + p01*p11; p3 = p10*p10 + p11*p11;
		w1 = w00*w00 + w10*w10; w2 = w00*w01 + w10*w11; w3 = w01*w01 + w11*w11;

		//area_now *= 2;

		h00 = area_now *(p1 + p2 + p2 + p3)*w1; h01 = -area_now *(p1 + p2)*w1; h02 = -area_now *(p2 + p3)*w1; h03 = area_now *(p1 + p2 + p2 + p3)*w2; h04 = -area_now *(p1 + p2)*w2; h05 = -area_now *(p2 + p3)*w2;
		h11 = area_now *p1*w1;                  h12 = area_now *p2*w1;    	 h13 = -area_now *(p1 + p2)*w2; h14 = area_now *p1*w2;                  h15 = area_now *p2*w2;
		h22 = area_now *p3*w1;                  h23 = -area_now *(p2 + p3)*w2; h24 = area_now *p2*w2;         h25 = area_now *p3*w2;
		h33 = area_now *(p1 + p2 + p2 + p3)*w3; h34 = -area_now *(p1 + p2)*w3; h35 = -area_now *(p2 + p3)*w3;
		h44 = area_now *p1*w3;                  h45 = area_now *p2*w3;
		h55 = area_now *p3*w3;


		det = j00*j11 - j01*j10;
		tr = (j00*j00 + j01*j01 + j10*j10 + j11*j11);

		d00 = -p00 - p10; d01 = p00; d02 = p10;
		d10 = -p01 - p11; d11 = p01; d12 = p11;

		r0 = area_now * ((1 + 1 / (det*det))*j00 - tr*j11 / (det*det*det));
		r1 = area_now * ((1 + 1 / (det*det))*j01 + tr*j10 / (det*det*det));
		r2 = area_now * ((1 + 1 / (det*det))*j10 + tr*j01 / (det*det*det));
		r3 = area_now * ((1 + 1 / (det*det))*j11 - tr*j00 / (det*det*det));


		pardiso_b[f0] -= r0*d00 + r1*d10;
		pardiso_b[f1] -= r0*d01 + r1*d11;
		pardiso_b[f2] -= r0*d02 + r1*d12;
		pardiso_b[f0 + V_N] -= r2*d00 + r3*d10;
		pardiso_b[f1 + V_N] -= r2*d01 + r3*d11;
		pardiso_b[f2 + V_N] -= r2*d02 + r3*d12;

		pardiso_a[id_h00[i]] += h00; pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02; pardiso_a[id_h03[i]] += h03; pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
		pardiso_a[id_h11[i]] += h11; pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h14[i]] += h14; pardiso_a[id_h15[i]] += h15;
		pardiso_a[id_h22[i]] += h22; pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24; pardiso_a[id_h25[i]] += h25;
		pardiso_a[id_h33[i]] += h33; pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
		pardiso_a[id_h44[i]] += h44; pardiso_a[id_h45[i]] += h45;
		pardiso_a[id_h55[i]] += h55;

	}

	for (int i = src_t_num; i < F_N; ++i)
	{
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
		if (f0 == -1&& f1 == -1&& f2 == -1)
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

	Energysource();
}

void Parafun_bij::CM()
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

	for (int i = 0; i < d_.m_T.rows(); i++)
	{
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


		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;

		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
		beta_0 = j00 - j11;  beta_1 = j10 + j01;

		alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1*alpha_1);
		beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1*beta_1);
		s1 = (p00)*(p00 + p10) + (p01)*(p01 + p11);
		s2 = (p10)*(p00 + p10) + (p11)*(p01 + p11);

		double h1 = p00*p00 + p01*p01;
		double h2 = p00*p10 + p01*p11;
		double h3 = p10*p10 + p11*p11;
		double h4 = p00*p11 - p01*p10;

		a1x0 = alpha_0*(-p00 - p10) + alpha_1*(p01 + p11);  a1x1 = alpha_0*p00 - alpha_1*p01; a1x2 = alpha_0*p10 - alpha_1*p11;
		a1x3 = alpha_0*(-p01 - p11) + alpha_1*(-p00 - p10); a1x4 = alpha_0*p01 + alpha_1*p00; a1x5 = alpha_0*p11 + alpha_1*p10;

		a2x0 = beta_0*(-p00 - p10) + beta_1*(-p01 - p11);   a2x1 = beta_0*p00 + beta_1*p01;   a2x2 = beta_0*p10 + beta_1*p11;
		a2x3 = beta_0*(p01 + p11) + beta_1*(-p00 - p10);    a2x4 = -beta_0*p01 + beta_1*p00;  a2x5 = -beta_0*p11 + beta_1*p10;

		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;

		hi_0 = 2 + 6 * 1 / (sig0*sig0*sig0*sig0); hi_1 = 2 + 6 * 1 / (sig1*sig1*sig1*sig1);

		aa = 0.25 / alpha_norm; bb = 0.25 / beta_norm;

		uu = aa*aa*(area_now*hi_0 + area_now*hi_1);
		vv = bb*bb*(area_now*hi_0 + area_now*hi_1);
		uv = aa*bb*(area_now*hi_0 - area_now*hi_1);

		h_u = area_now * (2 * sig0 - 2 * 1 / (sig0*sig0*sig0));
		h_v = area_now * (2 * sig1 - 2 * 1 / (sig1*sig1*sig1));

		walpha = h_u + h_v;
		wbeta = h_u - h_v;

		double hwa1 = (walpha * 0.25 / alpha_norm); double hwa2 = -(walpha * 0.25*0.25 / (alpha_norm*alpha_norm*alpha_norm));
		double hwb1 = (wbeta * 0.25 / beta_norm); double hwb2 = -(wbeta *0.25*0.25 / (beta_norm*beta_norm*beta_norm));


		h00 = uu*a1x0*a1x0 + vv*a2x0*a2x0 + uv*a1x0*a2x0 + uv*a2x0*a1x0; h01 = uu*a1x0*a1x1 + vv*a2x0*a2x1 + uv*a1x0*a2x1 + uv*a2x0*a1x1; h02 = uu*a1x0*a1x2 + vv*a2x0*a2x2 + uv*a1x0*a2x2 + uv*a2x0*a1x2; h03 = uu*a1x0*a1x3 + vv*a2x0*a2x3 + uv*a1x0*a2x3 + uv*a2x0*a1x3; h04 = uu*a1x0*a1x4 + vv*a2x0*a2x4 + uv*a1x0*a2x4 + uv*a2x0*a1x4; h05 = uu*a1x0*a1x5 + vv*a2x0*a2x5 + uv*a1x0*a2x5 + uv*a2x0*a1x5;

		h11 = uu*a1x1*a1x1 + vv*a2x1*a2x1 + uv*a1x1*a2x1 + uv*a2x1*a1x1; h12 = uu*a1x1*a1x2 + vv*a2x1*a2x2 + uv*a1x1*a2x2 + uv*a2x1*a1x2; h13 = uu*a1x1*a1x3 + vv*a2x1*a2x3 + uv*a1x1*a2x3 + uv*a2x1*a1x3; h14 = uu*a1x1*a1x4 + vv*a2x1*a2x4 + uv*a1x1*a2x4 + uv*a2x1*a1x4; h15 = uu*a1x1*a1x5 + vv*a2x1*a2x5 + uv*a1x1*a2x5 + uv*a2x1*a1x5;

		h22 = uu*a1x2*a1x2 + vv*a2x2*a2x2 + uv*a1x2*a2x2 + uv*a2x2*a1x2; h23 = uu*a1x2*a1x3 + vv*a2x2*a2x3 + uv*a1x2*a2x3 + uv*a2x2*a1x3; h24 = uu*a1x2*a1x4 + vv*a2x2*a2x4 + uv*a1x2*a2x4 + uv*a2x2*a1x4; h25 = uu*a1x2*a1x5 + vv*a2x2*a2x5 + uv*a1x2*a2x5 + uv*a2x2*a1x5;

		h33 = uu*a1x3*a1x3 + vv*a2x3*a2x3 + uv*a1x3*a2x3 + uv*a2x3*a1x3; h34 = uu*a1x3*a1x4 + vv*a2x3*a2x4 + uv*a1x3*a2x4 + uv*a2x3*a1x4; h35 = uu*a1x3*a1x5 + vv*a2x3*a2x5 + uv*a1x3*a2x5 + uv*a2x3*a1x5;

		h44 = uu*a1x4*a1x4 + vv*a2x4*a2x4 + uv*a1x4*a2x4 + uv*a2x4*a1x4; h45 = uu*a1x4*a1x5 + vv*a2x4*a2x5 + uv*a1x4*a2x5 + uv*a2x4*a1x5;

		h55 = uu*a1x5*a1x5 + vv*a2x5*a2x5 + uv*a1x5*a2x5 + uv*a2x5*a1x5;

		if (walpha >= 0)
		{
			h00 += hwa1*(s1 + s2) + hwa2*a1x0*a1x0; h01 += hwa1*(-s1) + hwa2*a1x0*a1x1; h02 += hwa1*(-s2) + hwa2*a1x0*a1x2; h03 += hwa2*a1x0*a1x3; h04 += hwa1*(h4)+hwa2*a1x0*a1x4; h05 += hwa1*(-h4) + hwa2*a1x0*a1x5;
			h11 += hwa1*(h1)+hwa2*a1x1*a1x1;        h12 += hwa1*(h2)+hwa2*a1x1*a1x2;    h13 += hwa1*(-h4) + hwa2*a1x1*a1x3; h14 += hwa2*a1x1*a1x4; h15 += hwa1*(h4)+hwa2*a1x1*a1x5;
			h22 += hwa1*(h3)+hwa2*a1x2*a1x2;        h23 += hwa1*(h4)+hwa2*a1x2*a1x3;    h24 += hwa1*(-h4) + hwa2*a1x2*a1x4; h25 += hwa2*a1x2*a1x5;
			h33 += hwa1*(s1 + s2) + hwa2*a1x3*a1x3; h34 += hwa1*(-s1) + hwa2*a1x3*a1x4; h35 += hwa1*(-s2) + hwa2*a1x3*a1x5;
			h44 += hwa1*(h1)+hwa2*a1x4*a1x4;        h45 += hwa1*(h2)+hwa2*a1x4*a1x5;
			h55 += hwa1*(h3)+hwa2*a1x5*a1x5;

		}
		h00 += hwb1*(s1 + s2) + hwb2*a2x0*a2x0; h01 += hwb1*(-s1) + hwb2*a2x0*a2x1; h02 += hwb1*(-s2) + hwb2*a2x0*a2x2; h03 += hwb2*a2x0*a2x3; h04 += hwb1*(-h4) + hwb2*a2x0*a2x4; h05 += hwb1*(h4)+hwb2*a2x0*a2x5;
		h11 += hwb1*(h1)+hwb2*a2x1*a2x1;        h12 += hwb1*(h2)+hwb2*a2x1*a2x2;    h13 += hwb1*(h4)+hwb2*a2x1*a2x3;    h14 += hwb2*a2x1*a2x4; h15 += hwb1*(-h4) + hwb2*a2x1*a2x5;
		h22 += hwb1*(h3)+hwb2*a2x2*a2x2;        h23 += hwb1*(-h4) + hwb2*a2x2*a2x3; h24 += hwb1*(h4)+hwb2*a2x2*a2x4;    h25 += hwb2*a2x2*a2x5;
		h33 += hwb1*(s1 + s2) + hwb2*a2x3*a2x3; h34 += hwb1*(-s1) + hwb2*a2x3*a2x4; h35 += hwb1*(-s2) + hwb2*a2x3*a2x5;
		h44 += hwb1*(h1)+hwb2*a2x4*a2x4;        h45 += hwb1*(h2)+hwb2*a2x4*a2x5;
		h55 += hwb1*(h3)+hwb2*a2x5*a2x5;

		u = aa*walpha; v = bb*wbeta;

		pardiso_b[f0] -= (u*a1x0 + v*a2x0);
		pardiso_b[f1] -= (u*a1x1 + v*a2x1);
		pardiso_b[f2] -= (u*a1x2 + v*a2x2);
		pardiso_b[f0 + V_N] -= (u*a1x3 + v*a2x3);
		pardiso_b[f1 + V_N] -= (u*a1x4 + v*a2x4);
		pardiso_b[f2 + V_N] -= (u*a1x5 + v*a2x5);

		pardiso_a[id_h00[i]] += h00; pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02; pardiso_a[id_h03[i]] += h03; pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
		pardiso_a[id_h11[i]] += h11; pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h14[i]] += h14; pardiso_a[id_h15[i]] += h15;
		pardiso_a[id_h22[i]] += h22; pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24; pardiso_a[id_h25[i]] += h25;
		pardiso_a[id_h33[i]] += h33; pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
		pardiso_a[id_h44[i]] += h44; pardiso_a[id_h45[i]] += h45;
		pardiso_a[id_h55[i]] += h55;
	}

	//scaffold computation
	for (int i = d_.m_T.rows(); i < F_N; i++)
	{
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
		beta_norm = 1e-10;

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
		else if(flag_==2)
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

	pardiso->a = pardiso_a;
	pardiso->rhs = pardiso_b;

	pardiso->factorize();
	pardiso->pardiso_solver();

	vector<double> result_d = pardiso->result;

	VectorXd negative_grad(2 * total_num), d(2 * total_num);
	negative_grad.setZero();
	d.setZero();

	for (int i = 0; i <V_N; i++)
	{
		negative_grad(var_ids[i]) = pardiso_b[i];
		negative_grad(var_ids[i]+total_num) = pardiso_b[i+V_N];
		d(var_ids[i]) = result_d[i];
		d(var_ids[i]+total_num) = result_d[i+V_N];
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

void Parafun_bij::max_step(const VectorXd &xx, const VectorXd &dd, double &step)
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
		tt = get_smallest_pos_quad_zero_bij( a, b, c);

		if (temp_t > tt)
		{
			temp_t = tt;
		}

	}
	step = temp_t;
}

double get_smallest_pos_quad_zero_bij(double a, double b, double c)
{
	using namespace std;
	double t1, t2;
	if (std::abs(a) < 1.0e-10)
	{
		a *= 1e6;
		b *= 1e6;
		c *= 1e6;
	}
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

void Parafun_bij::backtracking_line_search(const VectorXd &x, const VectorXd &d, const VectorXd &negetive_grad, double &alpha)
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

void Parafun_bij::Energy(const VectorXd &position, double &energyupdate)
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

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;


		det = j00*j11 - j01*j10;
		E_d = (1 + 1 / (det*det)) * (j00*j00 + j01*j01 + j10*j10 + j11*j11);

		energy += area[i] * E_d;
	}
	energyupdate = energy;
}

void Parafun_bij::Energysource()
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
		end_e_area += ((E_1 + E_2)*area[i]);
	}
	energy_uniform = end_e_one_temp /F_N;

	energy_area = end_e_area;
}
void Parafun_bij::MyEnergySource(std::vector<double>& face_enery)
{
	double end_e_one_temp = 0, end_e_area = 0;

	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_1, E_2;

	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	const double *pos = position_of_mesh.data();
	for (int i = 0; i < face_enery.size(); ++i)
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
		end_e_area += ((E_1 + E_2)*area[i]);
		face_enery[i] = ((E_1 + E_2)*area[i]);
	}
	energy_uniform = end_e_one_temp / F_N;

	energy_area = end_e_area;
}

double Parafun_bij::compute_energy(const Eigen::MatrixXd & x, bool whole)
{
	double end_e_one_temp = 0, end_e_area = 0;

	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_1, E_2;

	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	const double *pos = x.data();
	int src_t_num = d_.m_T.rows();

	for (int i = 0; i < src_t_num; ++i)
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

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

		det = j00 * j11 - j01 * j10;

		E_1 = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);
		E_2 = 1.0 / (det*det)* E_1;

		end_e_one_temp += E_1;
		end_e_one_temp += E_2;
		end_e_area += ((E_1 + E_2)*area_src[i]);
	}
	//cout << "compute energy " << end_e_area << endl;
	if (whole)
	{
		for (int i = src_t_num; i < F_N; ++i)
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

			j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

			det = j00 * j11 - j01 * j10;

			E_1 = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);
			E_2 = 1.0 / (det*det)* E_1;

			end_e_one_temp += E_1;
			end_e_one_temp += E_2;
			end_e_area += ((E_1 + E_2)*area_scaf[i- src_t_num]);
		}
	}
	//cout << "compute energy with scaf "<<whole<<" "<< end_e_area << endl;
	//cout << "compute energy with scaf/area " << whole << " " << end_e_area/d_.mesh_measure << endl;

	return end_e_area;
}

void Parafun_bij::local_coordinate_inverse(int i, double &p00, double &p01, double &p10, double &p11)
{
	int f0 = F0[i];
	int f1 = F1[i];
	int f2 = F2[i];
	
	Vector3d x_(d_.m_V(f1,0)- d_.m_V(f0,0), d_.m_V(f1, 1) - d_.m_V(f0, 1), d_.m_V(f1, 2) - d_.m_V(f0, 2));
	double x1_0 = x_.norm();
	x_/=x1_0;
	Vector3d l_(d_.m_V(f2, 0) - d_.m_V(f0, 0), d_.m_V(f2, 1) - d_.m_V(f0, 1), d_.m_V(f2, 2) - d_.m_V(f0, 2));

	Vector3d n_ = x_.cross(l_);
	n_.normalize();
	Vector3d y_ = n_.cross(x_);
	double x2_0 = l_.dot(x_);
	double y2_0 = l_.dot(y_);

	p00 = 1 / x1_0;
	p01 = -x2_0 / (x1_0*y2_0);
	p10 = 0;
	p11 = 1 / y2_0;
}

void Parafun_bij::local_coordinate_inverse_scaf(int i, double & p00, double & p01, double & p10, double & p11)
{
	int f0 = F0[i];
	int f1 = F1[i];
	int f2 = F2[i];

	Vector2d x_(d_.w_uv(f1, 0) - d_.w_uv(f0, 0), d_.w_uv(f1, 1) - d_.w_uv(f0, 1));
	Vector2d l_(d_.w_uv(f2, 0) - d_.w_uv(f0, 0), d_.w_uv(f2, 1) - d_.w_uv(f0, 1));

	double area_tri = abs(x_(0)*l_(1) - x_(1)*l_(0));
	double x1_0, x2_0, y2_0;
	if (area_tri > area_threshold)
	{
		x1_0 = x_.norm();
		x_ /= x1_0;
		Vector2d y_(-x_(1), x_(0));
		x2_0 = l_.dot(x_);
		y2_0 = l_.dot(y_);
	}
	else
	{
		//cout << "area too small!!!!!!!!!!!!! " << endl;
		double h = sqrt((2*area_threshold) / sqrt(3.0));
		x1_0 = h;
		x2_0 = h / 2.0;
		y2_0 = sqrt(3.0)*h / 2.0;
	}
	p00 = 1 / x1_0;
	p01 = -x2_0 / (x1_0*y2_0);
	p10 = 0;
	p11 = 1 / y2_0;
}

double Parafun_bij::newton_equation(const double & a, const double & b, const double & K)
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

void Parafun_bij::adjust_scaf_weight(double new_weight)
{
	d_.scaffold_factor = new_weight;
	d_.update_scaffold();
	init_area();
}

void Parafun_bij::handle_mintri()
{
	double min_bnd_edge_len = numeric_limits<double>::infinity();
	int acc_bnd = 0;
	for (int i = 0; i < d_.bnd_sizes.size(); i++)
	{
		int current_size = d_.bnd_sizes[i];

		for (int e = acc_bnd; e < acc_bnd + current_size - 1; e++)
		{
			min_bnd_edge_len = (std::min)(min_bnd_edge_len,
				(d_.w_uv.row(d_.internal_bnd(e)) -
					d_.w_uv.row(d_.internal_bnd(e + 1)))
				.squaredNorm());
		}
		min_bnd_edge_len = (std::min)(min_bnd_edge_len,
			(d_.w_uv.row(d_.internal_bnd(acc_bnd)) -
				d_.w_uv.row(d_.internal_bnd(acc_bnd + current_size -
					1)))
			.squaredNorm());
		acc_bnd += current_size;
	}

	area_threshold = min_bnd_edge_len / 4.0;
}
