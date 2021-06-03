//
// Created by Zhongshi Jiang on 2/12/17.
//

#include "dchart/ScafData.h"
#include <iostream>
#include "Algorithm/TriangleInterface.h"
//#include <igl/boundary_loop.h>


using namespace std;
using namespace Eigen;

void ScafData::update_scaffold()
{
  mv_num = m_V.rows();
  mf_num = m_T.rows();

  v_num = w_uv.rows();
  sf_num = s_T.rows();

  sv_num = v_num - mv_num;
  f_num = sf_num + mf_num;

  s_M = Eigen::VectorXd::Constant(sf_num, scaffold_factor);
}

void ScafData::mesh_improve(bool in_packing = false) {

	MatrixXd m_uv = w_uv.topRows(mv_num);
	MatrixXd V_bnd;
	V_bnd.resize(internal_bnd.size(), 2);
	for (int i = 0; i < internal_bnd.size(); i++) // redoing step 1.
	{
		V_bnd.row(i) = m_uv.row(internal_bnd(i));
	}

	if (rect_frame_V.size() == 0) {
		Matrix2d ob;// = rect_corners;
		{
			VectorXd uv_max = m_uv.colwise().maxCoeff();
			VectorXd uv_min = m_uv.colwise().minCoeff();
			VectorXd uv_mid = (uv_max + uv_min) / 2.;

			//double scaf_range = 3;
			Eigen::Array2d scaf_range(3, 3);
			ob.row(0) = uv_mid.array() + scaf_range * ((uv_min - uv_mid).array());
			ob.row(1) = uv_mid.array() + scaf_range * ((uv_max - uv_mid).array());
		}
		Vector2d rect_len;
		rect_len << ob(1, 0) - ob(0, 0), ob(1, 1) - ob(0, 1);
		int frame_points = 20;

		rect_frame_V.resize(4 * frame_points, 2);
		for (int i = 0; i < frame_points; i++) {
			// 0,0;0,1
			rect_frame_V.row(i) << ob(0, 0), ob(0, 1) + i * rect_len(1) / frame_points;
			// 0,0;1,1
			rect_frame_V.row(i + frame_points)
				<< ob(0, 0) + i * rect_len(0) / frame_points, ob(1, 1);
			// 1,0;1,1
			rect_frame_V.row(i + 2 * frame_points) << ob(1, 0), ob(1, 1)
				- i * rect_len(1) /
				frame_points;
			// 1,0;0,1
			rect_frame_V.row(i + 3 * frame_points)
				<< ob(1, 0) - i * rect_len(0) / frame_points, ob(0, 1);
			// 0,0;0,1
		}
		frame_ids = Eigen::VectorXi::LinSpaced(rect_frame_V.rows(), mv_num,
			mv_num + rect_frame_V.rows());
	}

	// Concatenate Vert and Edge
	MatrixXd V;
	MatrixXi E;

	{
		V.resize(V_bnd.rows() + rect_frame_V.rows(), V_bnd.cols());
		V << V_bnd, rect_frame_V;
	}
	E.resize(V.rows(), 2);
	for (int i = 0; i < E.rows(); i++)
		E.row(i) << i, i + 1;
	int acc_bs = 0;
	for (auto bs : bnd_sizes) {
		E(acc_bs + bs - 1, 1) = acc_bs;
		acc_bs += bs;
	}
	E(V.rows() - 1, 1) = acc_bs;
	assert(acc_bs == internal_bnd.size());


	//ycy H 存储的是每个chart的第一个三角形的重心
	MatrixXd H = MatrixXd::Zero(component_sizes.size(), 2);
	{
		int hole_f = 0;
		int hole_i = 0;
		for (auto cs : component_sizes) {
			for (int i = 0; i < 3; i++)
				H.row(hole_i) += m_uv.row(m_T(hole_f, i)); // redoing step 2
			hole_f += cs;
			hole_i++;
		}
	}
	H /= 3.;

	MatrixXd uv2;
	triangluate_bij(V, E, H, uv2, s_T);
	auto bnd_n = internal_bnd.size();

	for (auto i = 0; i < s_T.rows(); i++)
		for (auto j = 0; j < s_T.cols(); j++) {
			auto &x = s_T(i, j);
			if (x < bnd_n) x = internal_bnd(x);
			else x += m_uv.rows() - bnd_n;
		}
	
	
	{
		surface_F.resize(m_T.rows() + s_T.rows(), 3);
		surface_F << m_T, s_T;
	}


	w_uv.conservativeResize(m_uv.rows() - bnd_n + uv2.rows(), 2);
	w_uv.bottomRows(uv2.rows() - bnd_n) = uv2.bottomRows(-bnd_n + uv2.rows());
	update_scaffold();
}

void ScafData::add_new_patch(const Eigen::MatrixXd &V_in,
                             const Eigen::MatrixXi &F_ref,
                             const Eigen::RowVectorXd &center) {
  using namespace std;
  using namespace Eigen;

  VectorXd M;
  M.resize(F_ref.rows());
  {
	  Eigen::Vector3d v0, v1, v2, e01, e02, normal_f;
	  double area_f;
	  for (size_t i = 0; i < F_ref.rows(); i++)
	  {
		  v0 = V_in.row(F_ref(i,0));
		  v1 = V_in.row(F_ref(i,1));
		  v2 = V_in.row(F_ref(i,2));
		  e01 = v1 - v0;
		  e02 = v2 - v0;
		  normal_f = e01.cross(e02);
		  area_f = normal_f.norm() / 2.0;
		  M(i) = area_f;
	  }
  }


  Eigen::MatrixXd V_ref = V_in;
  Eigen::MatrixXd uv_init;
  Eigen::VectorXi bnd;
  Eigen::MatrixXd bnd_uv;
  std::vector<std::vector<int>> all_bnds;
  boundary_loop(F_ref, all_bnds);
  int num_holes = all_bnds.size() - 1;

  std::sort(all_bnds.begin(), all_bnds.end(), [](std::vector<int>& a, std::vector<int>&b){
    return a.size() > b.size();
  });

  bnd =  Map<Eigen::VectorXi>(all_bnds[0].data(),
                              all_bnds[0].size());

  map_vertices_to_circle(V_ref, bnd, bnd_uv);
  bnd_uv *= sqrt(M.sum()/M_PI);
  bnd_uv.rowwise() += center;
  mesh_measure += M.sum()/2;
  std::cout<<"Mesh Measure"<< M.sum()/2<<"number holes "<< num_holes <<std::endl;

  if(num_holes == 0) {
    if (bnd.rows() == V_ref.rows()) {
      std::cout << "All vert on boundary" << std::endl;
      uv_init.resize(V_ref.rows(), 2);
      for (int i = 0; i < bnd.rows(); i++) {
        uv_init.row(bnd(i)) = bnd_uv.row(i);
      }
    } else {
      Tutte(V_ref.rows(), F_ref, bnd, bnd_uv, uv_init);
    }
  } else {
    auto &F = F_ref;
    auto &V = V_in;
    auto &primary_bnd = bnd;
    // fill holes
    int n_filled_faces = 0;
    int real_F_num = F.rows();
    for (int i = 0; i < num_holes; i++) {
      n_filled_faces += all_bnds[i + 1].size();
    }
    MatrixXi F_filled(n_filled_faces + real_F_num, 3);
    F_filled.topRows(real_F_num) = F;

    int new_vert_id = V.rows();
    int new_face_id = real_F_num;

    for (int i = 0; i < num_holes; i++) {
      int cur_bnd_size = all_bnds[i + 1].size();
      auto it = all_bnds[i + 1].begin();
      auto back = all_bnds[i + 1].end() - 1;
      F_filled.row(new_face_id++) << *it, *back, new_vert_id;
      while (it != back) {
        F_filled.row(new_face_id++)
            << *(it + 1), *(it), new_vert_id;
        it++;
      }
      new_vert_id++;
    }
    assert(new_face_id == F_filled.rows());
    assert(new_vert_id == V.rows() + num_holes);

    Tutte(V_ref.rows()+num_holes,F_filled, primary_bnd, bnd_uv, uv_init);
    uv_init.conservativeResize(V.rows(), 2);
  }

  component_sizes.push_back(F_ref.rows());

  MatrixXd m_uv = w_uv.topRows(mv_num);

  w_uv.resize(m_uv.rows() + uv_init.rows(), 2);
  if (mv_num == 0)
  {
	  w_uv = uv_init;
  }
  else
  {
	  w_uv << m_uv, uv_init;
  } 
  m_M.conservativeResize(mf_num + M.size());
  m_M.bottomRows(M.size()) = M;


  for(auto cur_bnd : all_bnds) {
    internal_bnd.conservativeResize(internal_bnd.size()+ cur_bnd.size());
    internal_bnd.bottomRows(cur_bnd.size()) =
        Map<ArrayXi>(cur_bnd.data(),cur_bnd.size()) + mv_num;
    bnd_sizes.push_back(cur_bnd.size());
  }

  m_T.conservativeResize(mf_num + F_ref.rows(), 3);
  m_T.bottomRows(F_ref.rows()) = F_ref.array() + mv_num;
  mf_num += F_ref.rows();

  m_V.conservativeResize(mv_num + V_ref.rows(), 3);
  m_V.bottomRows(V_ref.rows()) = V_ref;
  mv_num += V_ref.rows();

  rect_frame_V.resize(0, 0);

  mesh_improve(false);
}
//
//void boundary_loop(const Eigen::MatrixXi &F_ref, std::vector<std::vector<int>>& boundaryloop)
//{
//	std::vector<std::vector<int>> boundaryEdges;
//	std::vector<std::vector<int>> edges;
//	int n_fvs = F_ref.cols();
//
//	for (int it = 0; it < F_ref.rows(); it++)
//	{
//		for (int i = 0; i < n_fvs; i++)
//		{
//			int var = F_ref(it, i);
//			int var_n = F_ref(it, (i + 1) % n_fvs);
//			if (var > var_n) std::swap(var, var_n);
//			std::vector<int> edge(4);
//			edge[0] = var;
//			edge[1] = var_n;
//			edge[2] = it;
//			edge[3] = i;
//			edges.emplace_back(edge);
//		}
//	}
//	std::sort(edges.begin(), edges.end());
//	int i = 1;
//	for (; i < edges.size();)
//	{
//		auto& r1 = edges[i - 1];
//		auto& r2 = edges[i];
//		if ((r1[0] == r2[0]) && (r1[1] == r2[1]))
//		{
//			i += 2;
//		}
//		else
//		{
//			boundaryEdges.emplace_back(edges[i - 1]);
//			i++;
//		}
//	}
//	if (i == edges.size())
//		boundaryEdges.emplace_back(edges.back());
//
//	for (auto&var : boundaryEdges)
//	{
//		var[0] = F_ref(var[2], var[3]);
//		var[1] = F_ref(var[2], (var[3] + 1) % n_fvs);
//	}
//	int ev0 = boundaryEdges.front()[0];
//	int ev1 = boundaryEdges.front()[1];
//
//	vector<int> visited;
//	visited.resize(boundaryEdges.size(), 0);
//	visited[0] = 1;
//	vector<int> loop0;
//	loop0.push_back(ev1);
//	while (ev1 != ev0)
//	{
//		for (int i = 1; i < boundaryEdges.size(); i++)
//		{
//			if (visited[i] == 1)
//				continue;
//			if (boundaryEdges[i][0] == ev1)
//			{
//				visited[i] = 1;
//				ev1 = boundaryEdges[i][1];
//				loop0.push_back(ev1);
//				break;
//			}
//		}
//	}
//	boundaryloop.emplace_back(loop0);
//}
//
//void map_vertices_to_circle(const Eigen::MatrixXd& V, const Eigen::VectorXi& bnd, Eigen::MatrixXd& UV)
//{
//	// Get sorted list of boundary vertices
//	std::vector<int> interior, map_ij;
//	map_ij.resize(V.rows());
//
//	std::vector<bool> isOnBnd(V.rows(), false);
//	for (int i = 0; i < bnd.size(); i++)
//	{
//		isOnBnd[bnd[i]] = true;
//		map_ij[bnd[i]] = i;
//	}
//
//	for (int i = 0; i < (int)isOnBnd.size(); i++)
//	{
//		if (!isOnBnd[i])
//		{
//			map_ij[i] = interior.size();
//			interior.push_back(i);
//		}
//	}
//
//	// Map boundary to unit circle
//	std::vector<double> len(bnd.size());
//	len[0] = 0.;
//
//	for (int i = 1; i < bnd.size(); i++)
//	{
//		len[i] = len[i - 1] + (V.row(bnd[i - 1]) - V.row(bnd[i])).norm();
//	}
//	double total_len = len[len.size() - 1] + (V.row(bnd[0]) - V.row(bnd[bnd.size() - 1])).norm();
//
//	UV.resize(bnd.size(), 2);
//	for (int i = 0; i < bnd.size(); i++)
//	{
//		double frac = len[i] * 2. * M_PI / total_len;
//		UV.row(map_ij[bnd[i]]) << cos(frac), sin(frac);
//	}
//
//}
