#include "BSpline.h"
#include <fstream>

using namespace Eigen;

BSpline::BSpline()
{
}


BSpline::~BSpline()
{
}

double BSpline::B_base(std::vector<double>& tlist, double t, int k, int i)
{
	double result = 1;
	if (k == 1)
	{
		if (tlist[i] <= t && t < tlist[i + 1])
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	
	double alpha, beta;
	if (tlist[i + k - 1] - tlist[i] == 0)
	{
		alpha = 0;
	}
	else
	{
		alpha = (t - tlist[i]) / (tlist[i + k - 1] - tlist[i]);
	}

	if (tlist[i + k] - tlist[i + 1] == 0)
	{
		beta = 0;
	}
	else
	{
		beta = (tlist[i + k] - t) / (tlist[i + k] - tlist[i + 1]);
	}

	result = alpha*B_base(tlist, t, k - 1, i) + beta*B_base(tlist, t, k - 1, i + 1);
	return result;
}
Vector3f BSpline::DB_base(double t)
{
	int n = time_points.size();
	double alltime = time_points[time_points.size() - 1];

	std::vector<double> Dtime_point(time_points.size() - 2);
	for (int i = 1; i < time_points.size()-1; i++)
	{
		Dtime_point[i - 1] = time_points[i];
	}

	int s = 0;
	for (int u = 0; u < n; u++)
	{
		if (time_points[u] < t+1e-6 && time_points[u + 1] > t)
		{
			s = u;
			break;
		}
	}
	double b = 0;
	Vector3f dt = Vector3f::Zero();
	for (int i = s + 2 - 4; i <= s; i++)
	{
		b = B_base(time_points, t, 3, i);
		dt = dt + b*(Contral_Point[i] - Contral_Point[i - 1]);
	}
	return dt;
}

void BSpline::CalculatePointlist()
{
	draw_Contral_Point.clear();
	int t_num = time_points.size() / 5;
	int nn = t_num * 100 + 100;


	double t_fin = time_points[time_points.size() - 1] - 1e-4;
	double t_be = 0 + 1e-4;
	double dt = (t_fin - t_be) / nn;
	for (int i = 0; i <= nn; i++)
	{
		double ft = t_be + i*dt;
		Eigen::Vector3f r;
		CalculatePoint(Contral_Point_ori, ft, r);
		draw_Contral_Point.push_back(r);
	}
}
void BSpline::CalculatePointlist_finally()
{
	draw_first_Point.clear();
	draw_second_Point.clear();

	int t_num = time_points.size() / 5;
	int nn = t_num * 100 + 100;

	double t_fin = time_points[time_points.size() - 1] - 1e-4;
	double t_be = 0 + 1e-4;
	double dt = (t_fin - t_be) / nn;
	for (int i = 0; i <= nn; i++)
	{
		double ft = t_be + i*dt;
		Eigen::Vector3f r;
		CalculatePoint(finally_first_Point, ft, r);
		draw_first_Point.push_back(r);
		CalculatePoint(finally_second_Point, ft, r);
		draw_second_Point.push_back(r);
	}
}
void BSpline::CalculatePoint(std::vector<Eigen::Vector3f>& control, double t, Eigen::Vector3f& r)
{
	int n_size = control.size();
	r.setZero();
	for (int i = 0; i < n_size; i++)
	{
		double b = 0;
		if (t < time_points[i] || (i + 4 <= time_points.size() - 1 && t > time_points[i + 4]))
		{
			b = 0;
			continue;
		}
		else
		{
			b = B_base(time_points, t, 4, i);
		}
		r = r + b*control[i];
	}
}

void BSpline::get_sharp_spline(double theta, Vector3f& left_fix, Vector3f& right_fix)
{
	Matrix3f R1;
	R1 << cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1;
	Matrix3f R2 = R1.inverse();
	
	std::vector<Vector3f> left_control;
	std::vector<Vector3f> right_control;

	int n = Contral_Point.size();
	left_control.resize(n);
	right_control.resize(n);
	for (int i = 0; i < n; i++)
	{
		left_control[i] = Contral_Point[0] + R1*(Contral_Point[i] - Contral_Point[0]);
		right_control[i] = Contral_Point[0] + R2*(Contral_Point[i] - Contral_Point[0]);
	}
	left_control[n - 1] = left_fix;
	right_control[n - 1] = right_fix;

	std::vector<Vector3f> fix_infor_left(3);
	std::vector<Vector3f> fix_infor_right(3);

	fix_infor_left[0] = left_control[0];
	fix_infor_left[1] = left_control[n - 1];
	fix_infor_left[2] = R1 * DB_base(1e-5);

	fix_infor_right[0] = right_control[0];
	fix_infor_right[1] = right_control[n - 1];
	fix_infor_right[2] = R2*DB_base(1e-5);

	find_Colset(left_control, fix_infor_left);
	find_Colset(right_control, fix_infor_right);

	finally_first_Point = left_control;
	finally_contral_Point = Contral_Point;
	finally_second_Point = right_control;
}

void BSpline::find_Colset(std::vector<Vector3f>& new_control, std::vector<Vector3f>& fix_infor)
{
	double error = 1;
	while (error > 1e-2)
	{
		error = 0;

		std::vector<Eigen::Vector2f> p_list;
		std::vector<Eigen::Vector2f> q_list;
		p_list.clear();
		q_list.clear();
		for (int i = 1; i < Contral_Point.size(); i++)
		{
			p_list.push_back((new_control[i] - new_control[i - 1]).head(2));
			q_list.push_back((Contral_Point[i] - Contral_Point[i - 1]).head(2));
		}

		Vector2f mu_p = Vector2f::Zero();
		Vector2f mu_q = Vector2f::Zero();
		for (int i = 0; i < p_list.size(); i++)
		{
			mu_p += p_list[i];
			mu_q += q_list[i];
		}
		mu_p /= p_list.size();
		mu_q /= q_list.size();

		Matrix2f H = Matrix2f::Zero();
		for (int i = 0; i < p_list.size(); i++)
		{
			H += (p_list[i] - mu_p) * (q_list[i] - mu_q).transpose();
			//		H += p_list[i] * q_list[i].transpose();
		}
		JacobiSVD<Eigen::Matrix2f> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Matrix2f R = svd.matrixV() * (svd.matrixU().transpose());

		int n = Contral_Point.size();
		MatrixXf A = MatrixXf(2 * (n - 1), 2 * n);
		A.setZero();
		VectorXf b(2 * (n-1) );
		for (int i = 0; i < n - 1; i++)
		{
			A.block(2 * i, 2 * i, 2, 2) = -1 * R;
			A.block(2 * i, 2 * (i + 1), 2, 2) = R;
			b(2 * i) = q_list[i].x();
			b(2 * i + 1) = q_list[i].y();
		}
		MatrixXf M = MatrixXf(6, 2 * n);
		M.setZero();
		M.block(0, 0, 2, 2) = Matrix2f::Identity();
		M.block(2, 2 * (n - 1), 2, 2) = Matrix2f::Identity();

		M(4, 0) = -B_base(time_points, 1e-5, 3, 1);
		M(4, 2) = B_base(time_points, 1e-5, 3, 1) - B_base(time_points, 1e-5, 3, 2);
		M(4, 4) = B_base(time_points, 1e-5, 3, 2) - B_base(time_points, 1e-5, 3, 3);
		M(4, 6) = B_base(time_points, 1e-5, 3, 3);
		M(5, 1) = -B_base(time_points, 1e-5, 3, 1);
		M(5, 3) = B_base(time_points, 1e-5, 3, 1) - B_base(time_points, 1e-5, 3, 2);
		M(5, 5) = B_base(time_points, 1e-5, 3, 2) - B_base(time_points, 1e-5, 3, 3);
		M(5, 7) = B_base(time_points, 1e-5, 3, 3);

		VectorXf c(6);
		c.segment(0, 2) = fix_infor[0].head(2);
		c.segment(2, 2) = fix_infor[1].head(2);
		c.segment(4, 2) = fix_infor[2].head(2);

		int m = 2*n;
		MatrixXf K(m + 6, m + 6);
		K.setZero();
		K.block(0, 0, m, m) = A.transpose()*A;
		K.block(0, m, m, 6) = M.transpose();
		K.block(m, 0, 6, m) = M;
		VectorXf f(m + 6);
		f.segment(0, m) = A.transpose()*b;
		f.segment(m, 6) = c;

		VectorXf u;
		u = K.lu().solve(f);

		for (int i = 0; i < new_control.size(); i++)
		{
			error += (Vector3f(u(2 * i), u(2 * i + 1), 0) - new_control[i]).norm();
			new_control[i] = Vector3f(u(2 * i), u(2 * i + 1), 0);
		}
	}
}

void BSpline::GetSpline_from_data(std::vector<Vector3f>& control_data, std::vector<double>& knock)
{
	Contral_Point.resize(control_data.size());
	for (int i = 0; i < control_data.size(); i++)
	{
		Contral_Point[i] = control_data[i];
	}
	time_points.resize(knock.size());
	for (int i = 0; i < knock.size(); i++)
	{
		time_points[i] = knock[i];
	}
}
void BSpline::GetSpline_from_data(std::vector<Vector3f>& control_p)
{
	if (control_p.size() < 2) return;
	if (control_p.size() == 2)
	{
		Contral_Point_ori.resize(4);
		Contral_Point_ori[0] = control_p[0];
		Contral_Point_ori[1] = 0.66*control_p[0] + 0.33*control_p[1];
		Contral_Point_ori[2] = 0.33*control_p[0] + 0.66*control_p[1];
		Contral_Point_ori[3] = control_p[1];
	}else if (control_p.size() == 3)
	{
		Contral_Point_ori.resize(4);
		Contral_Point_ori[0] = control_p[0];
		Contral_Point_ori[1] = 0.5*control_p[0] + 0.5*control_p[1];
		Contral_Point_ori[2] = 0.5*control_p[1] + 0.5*control_p[2];
		Contral_Point_ori[3] = control_p[2];
	}
	else
	{
		Contral_Point_ori.resize(control_p.size());
		for (int i = 0; i < control_p.size(); i++)
		{
			Contral_Point_ori[i] = control_p[i];
		}
	}
	time_points.clear();
	time_points.push_back(0);
	time_points.push_back(0);
	time_points.push_back(0);
	for (int i = 0; i <= Contral_Point_ori.size() - 3; i++)
	{
		time_points.push_back(i);
	}
	double fin = Contral_Point_ori.size() - 3;
	time_points.push_back(fin);
	time_points.push_back(fin);
	time_points.push_back(fin);

}