#include "BijectivePara.h"

extern OpenMesh::VPropHandleT<Eigen::Vector3f> vert_ref;
extern OpenMesh::VPropHandleT<Eigen::Vector3f> user_vert_ref;

BijectivePara::BijectivePara()
{
}

BijectivePara::~BijectivePara()
{
}

void BijectivePara::parameterization(Mesh& mesh_, std::vector<double>& face_energy)
{
	using namespace Eigen;
	using namespace std;

	bool auto_weight = true;
	bool use_CM = true;

	bool adjust_frame = true;
	int iteration_count =0;

	parafun_solver->after_mesh_improve();

	static double last_mesh_energy = parafun_solver->compute_energy(scaf_data.w_uv, false) / scaf_data.mesh_measure - 2 * parafun_solver->dim;

	for (int i = 0; i < MAX_ITER_NUM; i++) {
		std::cout << "=============" << std::endl;
		std::cout << "Iteration:" << iteration_count++ << '\t';

		//if(!s_.demotype)
			//d_.rect_frame_V.resize(0, 0);
		scaf_data.mesh_improve(true);
		parafun_solver->after_mesh_improve();

		if (auto_weight)
			parafun_solver->adjust_scaf_weight(
			(last_mesh_energy)*scaf_data.mesh_measure / (scaf_data.sf_num) / 100.0);
		scaf_data.energy = parafun_solver->perform_iteration_cm(use_CM);


		double current_mesh_energy =
			parafun_solver->compute_energy(scaf_data.w_uv, false) / scaf_data.mesh_measure - 2 * parafun_solver->dim;
		double mesh_energy_decrease = last_mesh_energy - current_mesh_energy;

		cout << "Energy After:"
			<< scaf_data.energy - 2 * scaf_data.dim
			<< "\tMesh Energy:"
			<< current_mesh_energy
			<< "\tEnergy Decrease"
			<< mesh_energy_decrease
			<< endl;
		cout << "V_num: " << scaf_data.v_num << " F_num: " << scaf_data.f_num << endl;

		double conv_rate = abs(mesh_energy_decrease) / last_mesh_energy;
		last_mesh_energy = current_mesh_energy;
		if (conv_rate < convgence_con_rate)
		{
            printf("Conv_rate get %.9f, so break out!\n", conv_rate);
			break;
		}
	}


	face_energy.clear();

	OpenMesh::Vec2f m_center = OpenMesh::Vec2f(0, 0);
	for (int i = 0; i < mesh_.n_vertices(); i++)
	{
		m_center += OpenMesh::Vec2f(scaf_data.w_uv(i, 0), scaf_data.w_uv(i, 1));
	}
	m_center = m_center / mesh_.n_vertices();

	int center_id = 0;
	double center_dis = (OpenMesh::Vec2f(scaf_data.w_uv(0, 0), scaf_data.w_uv(0, 1)) - m_center).norm();
	for (int i = 1; i < mesh_.n_vertices(); i++)
	{
		OpenMesh::Vec2f p(scaf_data.w_uv(i, 0), scaf_data.w_uv(i, 1));
		double dis = (p - m_center).norm();
		if (dis < center_dis)
		{
			center_dis = dis;
			center_id = i;
		}
	}

	Eigen::Vector3f t(scaf_data.w_uv(center_id, 0), scaf_data.w_uv(center_id, 1), 0);
	Eigen::Vector3f t1 = mesh_.property(vert_ref, mesh_.vertex_handle(center_id));

	Eigen::Vector3f r1(scaf_data.w_uv(0, 0), scaf_data.w_uv(0, 1), 0);
	Eigen::Vector3f r2 = mesh_.property(user_vert_ref, mesh_.vertex_handle(0));
	r1 = r1 - t;
	r2 = r2 - t1;
	double c1 = r1.cross(r2).z();
	double theta = acos(r1.dot(r2) / (r1.norm()*r2.norm()));
	Matrix3f R;
	R << cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1;
	if (c1 < 0)
	{
		Matrix3f R2 = R.inverse();
		R = R2;
	}

	for (int i = 0; i < scaf_data.mv_num; i++)
	{
		Eigen::Vector3f r(scaf_data.w_uv(i, 0), scaf_data.w_uv(i, 1), 0);
		mesh_.property(vert_ref, mesh_.vertex_handle(i)) = R*(r - t);
	}
}

void BijectivePara::load(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
	scaf_data = ScafData();
	scaf_data.add_new_patch(V, F, Eigen::RowVector2d(0, 0));
	parafun_solver.reset(new Parafun_bij(scaf_data));

}

void BijectivePara::load(Mesh& mesh_)
{
	Eigen::MatrixXi F(mesh_.n_faces(), 3);
	for (Mesh::FaceIter f_it = mesh_.faces_begin(); f_it != mesh_.faces_end(); f_it++)
	{
		int f_id = (*f_it).idx();
		Mesh::FaceVertexIter fv_it = mesh_.fv_begin(*f_it); F(f_id, 0) = (*fv_it).idx();
		++fv_it; F(f_id, 1) = (*fv_it).idx();
		++fv_it; F(f_id, 2) = (*fv_it).idx();
	}
	Eigen::MatrixXd V(mesh_.n_vertices(), 3);
	for (Mesh::VertexIter v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); v_it++)
	{
		int id = v_it.handle().idx();
		Mesh::Point p = mesh_.point(v_it.handle());
		V(id, 0) = p[0];
		V(id, 1) = p[1];
		V(id, 2) = p[2];
	}

	load(V, F);	
}

