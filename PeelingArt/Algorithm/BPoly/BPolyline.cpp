#include "BPolyline.h"
#include <fstream>
#include <iostream>
#include <set>

const int scrford = 3;
//const int scrford = 1;


BPolyline::BPolyline()
{
	polyline_x.clear();
	polyline_y.clear();
}


BPolyline::~BPolyline()
{
	polyline_x.clear();
	polyline_y.clear();
}

bool BPolyline::IntPolyline(std::vector<Eigen::Vector3f>& polint_list)
{
	if (polint_list.size() < 0)
	{
		return false;
	}
	polyline_x.clear();
	polyline_y.clear();

	double x_max = -1, y_max = -1;
	double x_min = 1, y_min = 1;

	for (int i = 0; i < polint_list.size(); i++)
	{
		polyline_x.push_back(polint_list[i].x());
		polyline_y.push_back(polint_list[i].y());

		if (x_max < polyline_x[i]) x_max = polyline_x[i];
		if (x_min > polyline_x[i]) x_min = polyline_x[i];
		if (y_max < polyline_y[i]) y_max = polyline_y[i];
		if (y_min > polyline_y[i]) y_min = polyline_y[i];
	}
	if ((polint_list[polint_list.size() - 1] - polint_list[0]).norm() < 1e-6)
	{
		polyline_x[polint_list.size() - 1] = polyline_x[0];
		polyline_y[polint_list.size() - 1] = polyline_y[0];
		is_close = true;
	}
	else
	{
		is_close = false;
	}

	double p_max = std::max(x_max, y_max);
	double p_min = std::min(x_min, y_min);
	double r_max = ceil(p_max * 10 + (p_max - p_min) / scrford);
	double r_min = floor(p_min * 10 - (p_max - p_min) / scrford);
	double lenth = (r_max - r_min)/10;

	double center_x = (x_max + x_min) / 2;
	double center_y = (y_max + y_min) / 2;
	Centralization(center_x, center_y);
	rectangle_max = Eigen::Vector3d(lenth / 2, lenth / 2, 0);
	rectangle_min = Eigen::Vector3d(-lenth / 2, -lenth / 2, 0);

	return is_close;
}
bool BPolyline::IntPolyline(std::vector<Eigen::Vector3f>& polint_list, bool is_close_)
{
	if (polint_list.size() < 0)
	{
		return false;
	}

	double x_max = -1, y_max = -1;
	double x_min = 1, y_min = 1;

	polyline_x.clear();
	polyline_y.clear();
	for (int i = 0; i < polint_list.size(); i++)
	{
		polyline_x.push_back(polint_list[i].x());
		polyline_y.push_back(polint_list[i].y());
		if (x_max < polyline_x[i]) x_max = polyline_x[i];
		if (x_min > polyline_x[i]) x_min = polyline_x[i];
		if (y_max < polyline_y[i]) y_max = polyline_y[i];
		if (y_min > polyline_y[i]) y_min = polyline_y[i];
	}
	is_close = is_close_;
	if (is_close_)
	{
		polyline_x.push_back(polyline_x[0]);
		polyline_y.push_back(polyline_y[0]);
	}

	double p_max = std::max(x_max, y_max);
	double p_min = std::min(x_min, y_min);
	double r_max = ceil(p_max * 10 + (p_max - p_min) / scrford);
	double r_min = floor(p_min * 10 - (p_max - p_min) / scrford);
	double lenth = (r_max - r_min) / 10;

	rectangle_max = Eigen::Vector3d(r_max / 10, r_max / 10, 0);
	rectangle_min = Eigen::Vector3d(r_min / 10, r_min / 10, 0);

	return is_close;
}

int BPolyline::PointInCurve(Eigen::Vector3f& p)
{
	if (!is_close)
	{
		std::cout << "the curve is not close" << std::endl;
		return -3;
	}
	if (p(2) != 0)
	{
		std::cout << "only computer xy Plane" << std::endl;
	}

	double theta = 0;
	for (int i = 0; i < polyline_x.size() - 1; i++)
	{
		OpenMesh::Vec3d v1(polyline_x[i] - p.x(), polyline_y[i] - p.y(), 0);
		OpenMesh::Vec3d v2(polyline_x[i + 1] - p.x(), polyline_y[i + 1] - p.y(), 0);
		double th = OpenMesh::dot(v1,v2) / (v1.norm()*v2.norm());
		if (th == -1 || v1.norm() < 1e-5)
		{
			return 0;
		}
		
		double sig = (OpenMesh::cross(v1, v2)[2])>0 ? 1 : -1;
		th = sig*acos(th);
		theta += th;
	}
	theta = abs(theta);

	if (theta / (2 * 3.1415926) > 0.9)
	{
		return 1;
	}
	else
	{
		return -1;
	}
	return -2;
}

void BPolyline::Centralization(double factor)
{
	double center_x = (rectangle_min + rectangle_max).x() / 2;
	double center_y = (rectangle_min + rectangle_max).y() / 2;

	for (int i = 0; i < polyline_x.size(); i++)
	{
		polyline_x[i] = factor*(polyline_x[i] - center_x);
		polyline_y[i] = factor*(polyline_y[i] - center_y);
	}

	rectangle_max = factor*(rectangle_max - Eigen::Vector3d(center_x, center_y, 0));
	rectangle_min = factor*(rectangle_min - Eigen::Vector3d(center_x, center_y, 0));
}
void BPolyline::Centralization(double center_x, double center_y)
{
	for (int i = 0; i < polyline_x.size(); i++)
	{
		polyline_x[i] = (polyline_x[i] - center_x);
		polyline_y[i] = (polyline_y[i] - center_y);
	}
}
void BPolyline::Centralization(double center_x, double center_y, double factor)
{
	for (int i = 0; i < polyline_x.size(); i++)
	{
		polyline_x[i] = factor*(polyline_x[i] - center_x);
		polyline_y[i] = factor*(polyline_y[i] - center_y);
	}

	rectangle_max = factor*(rectangle_max - Eigen::Vector3d(center_x, center_y, 0));
	rectangle_min = factor*(rectangle_min - Eigen::Vector3d(center_x, center_y, 0));
}


void BPolyline::ChangePolylineAndReduce(std::vector<int>& index_ori, std::vector<Eigen::Vector3f>& first_part, std::vector<Eigen::Vector3f>& secont_part)
{
	std::vector<double> new_x;
	std::vector<double> new_y;
	new_x.clear();
	new_y.clear();
	if (index_ori[0] < index_ori[2])
	{
		for (int i = 0; i < index_ori[0]; i++)
		{
			new_x.push_back(polyline_x[i]);
			new_y.push_back(polyline_y[i]);
		}
		int n = first_part.size();
		for (int i = 1; i < n; i++)
		{
			new_x.push_back(first_part[n - i].x());
			new_y.push_back(first_part[n - i].y());
		}
		n = secont_part.size();
		for (int i = 0; i < n; i++)
		{
			new_x.push_back(secont_part[i].x());
			new_y.push_back(secont_part[i].y());
		}
		for (int i = index_ori[2]; i < polyline_x.size(); i++)
		{
			new_x.push_back(polyline_x[i]);
			new_y.push_back(polyline_y[i]);
		}
	}

	if (index_ori[2] < index_ori[0])
	{
		for (int i = index_ori[2]; i <= index_ori[0]; i++)
		{
			new_x.push_back(polyline_x[i]);
			new_y.push_back(polyline_y[i]);
		}
		int n = first_part.size();
		for (int j = 1; j < n; j++)
		{
			new_x.push_back(first_part[n - j].x());
			new_y.push_back(first_part[n - j].y());
		}
		n = secont_part.size();
		for (int j = 0; j < n; j++)
		{
			new_x.push_back(secont_part[j].x());
			new_y.push_back(secont_part[j].y());
		}
	}

	polyline_x.clear();
	polyline_y.clear();
	polyline_x.push_back(new_x[0]);
	polyline_y.push_back(new_y[0]);
	int n = 0;
	for (int i = 0; i < new_x.size()-1; i++)
	{
		if (abs(new_x[i] - polyline_x[n]) + abs(new_y[i] - polyline_y[n]) > 0.03)
		{
			polyline_x.push_back(new_x[i]);
			polyline_y.push_back(new_y[i]);
			n++;
		}
	}
	double dis = sqrt((polyline_x[0] - polyline_x[n])*(polyline_x[0] - polyline_x[n]) + (polyline_y[0] - polyline_y[n])*(polyline_y[0] - polyline_y[n]));
	if (is_close == true && dis > 1e-5)
	{
		polyline_x.push_back(polyline_x[0]);
		polyline_y.push_back(polyline_y[0]);
	}
}

void BPolyline::ChangePolylineAndReduce(std::vector<int>& index_ori, std::vector<Eigen::Vector3f>& one_part)
{
	std::vector<double> new_x;
	std::vector<double> new_y;
	new_x.clear();
	new_y.clear();
	if (index_ori[0] < index_ori[2])
	{
		for (int i = 0; i < index_ori[0]; i++)
		{
			new_x.push_back(polyline_x[i]);
			new_y.push_back(polyline_y[i]);
		}
		int n = one_part.size();
		for (int i = 0; i < n; i++)
		{
			new_x.push_back(one_part[i].x());
			new_y.push_back(one_part[i].y());
		}
		for (int i = index_ori[2]; i < polyline_x.size(); i++)
		{
			new_x.push_back(polyline_x[i]);
			new_y.push_back(polyline_y[i]);
		}
	}

	if (index_ori[2] < index_ori[0])
	{
		for (int i = index_ori[2]; i <= index_ori[0]; i++)
		{
			new_x.push_back(polyline_x[i]);
			new_y.push_back(polyline_y[i]);
		}
		int n = one_part.size();
		for (int i = 0; i < n; i++)
		{
			new_x.push_back(one_part[i].x());
			new_y.push_back(one_part[i].y());
		}
	}

	polyline_x.clear();
	polyline_y.clear();
	polyline_x.push_back(new_x[0]);
	polyline_y.push_back(new_y[0]);
	int n = 0;
	for (int i = 0; i < new_x.size() - 1; i++)
	{
		if (abs(new_x[i] - polyline_x[n]) + abs(new_y[i] - polyline_y[n]) > 0.03)
		{
			polyline_x.push_back(new_x[i]);
			polyline_y.push_back(new_y[i]);
			n++;
		}
	}
	double dis = sqrt((polyline_x[0] - polyline_x[n])*(polyline_x[0] - polyline_x[n]) + (polyline_y[0] - polyline_y[n])*(polyline_y[0] - polyline_y[n]));
	if (is_close == true && dis > 1e-5)
	{
		polyline_x.push_back(polyline_x[0]);
		polyline_y.push_back(polyline_y[0]);
	}
}

void BPolyline::ChangePolylineAndReduce()
{
	std::vector<double> new_x;
	std::vector<double> new_y;
	new_x.clear();
	new_y.clear();
	for (int i = 0; i < polyline_x.size(); i++)
	{
		new_x.push_back(polyline_x[i]);
		new_y.push_back(polyline_y[i]);
	}
	polyline_x.clear();
	polyline_y.clear();
	polyline_x.push_back(new_x[0]);
	polyline_y.push_back(new_y[0]);
	int n = 0;
	for (int i = 0; i < new_x.size() - 1; i++)
	{
		if (abs(new_x[i] - polyline_x[n]) + abs(new_y[i] - polyline_y[n]) > 0.03)
		{
			polyline_x.push_back(new_x[i]);
			polyline_y.push_back(new_y[i]);
			n++;
		}
	}
	double dis = sqrt((polyline_x[0] - polyline_x[n])*(polyline_x[0] - polyline_x[n]) + (polyline_y[0] - polyline_y[n])*(polyline_y[0] - polyline_y[n]));
	if (is_close == true && dis > 1e-5)
	{
		polyline_x.push_back(polyline_x[0]);
		polyline_y.push_back(polyline_y[0]);
	}
}

void BPolyline::FindCenter(double& center_x, double& center_y)
{
	if (polyline_x.size() == 0)
	{
		return;
	}
	double x_max = -1, y_max = -1;
	double x_min = 1, y_min = 1;

	for (int i = 0; i < polyline_x.size(); i++)
	{
		if (x_max < polyline_x[i]) x_max = polyline_x[i];
		if (x_min > polyline_x[i]) x_min = polyline_x[i];
		if (y_max < polyline_y[i]) y_max = polyline_y[i];
		if (y_min > polyline_y[i]) y_min = polyline_y[i];
	}

	double p_max = std::max(x_max, y_max);
	double p_min = std::min(x_min, y_min);
	double r_max = ceil(p_max * 10 + (p_max - p_min) / scrford);
	double r_min = floor(p_min * 10 - (p_max - p_min) / scrford);
	double lenth = (r_max - r_min) / 10;

	center_x = (x_max + x_min) / 2;
	center_y = (y_max + y_min) / 2;
	Centralization(center_x, center_y);
	rectangle_max = Eigen::Vector3d(lenth / 2, lenth / 2, 0);
	rectangle_min = Eigen::Vector3d(-lenth / 2, -lenth / 2, 0);
}

double BPolyline::Id_Geodesic_Distance(int fidx, int sidx)
{

	int nn;
	if (is_close)
	{
		nn = polyline_x.size() - 1;
	}
	else
	{
		nn = polyline_x.size();
	}

	if (fidx < 0 || fidx > nn || sidx < 0 || sidx > nn)
	{
		return 99999;
	}

	int d1 = ((sidx - fidx) + nn) % nn;
	int d2 = ((fidx - sidx) + nn) % nn;

	if (d1 > d2)
	{
		double a = fidx;
		fidx = sidx;
		sidx = a;
	}

	double dis = 0;
	for (int i = (fidx + 1) % nn; i != sidx % nn; i = (i+1)%nn)
	{
		int j = (i - 1) % nn;
		if (j == -1)
		{
			j = nn;
		}
		dis += sqrt((polyline_x[i] - polyline_x[j])*(polyline_x[i] - polyline_x[j]) + (polyline_y[i] - polyline_y[j])*(polyline_y[i] - polyline_y[j]));
	}

	return dis;
}

void BPolyline::Add_insert_point(std::vector< std::vector<int> >& insert_list)
{
	std::vector<double> new_x;
	std::vector<double> new_y;
	new_x.clear();
	new_y.clear();
	insert_list.resize(8);
	for (int i = 0; i < 4; i++)
	{
		insert_list[i].push_back(-1);
	}


	for (int i = 0; i < polyline_x.size() - 1; i++)
	{
		new_x.push_back(polyline_x[i]);
		new_y.push_back(polyline_y[i]);

		Eigen::Vector3f a(polyline_x[i], polyline_y[i], 0);
		Eigen::Vector3f b(polyline_x[i + 1], polyline_y[i + 1], 0);
		double lambda;

		if (compute_insert(a, b, Eigen::Vector3f(0, 0, 0), Eigen::Vector3f(0, 10, 0)) > 0)
		{
			lambda = compute_insert(a, b, Eigen::Vector3f(0, 0, 0), Eigen::Vector3f(0, 10, 0));
			new_x.push_back(lambda*polyline_x[i] + (1 - lambda)*polyline_x[i + 1]);
			new_y.push_back(lambda*polyline_y[i] + (1 - lambda)*polyline_y[i + 1]);
			insert_list[0].push_back(new_x.size()-1);
		}
		else if (compute_insert(a, b, Eigen::Vector3f(0, 0, 0), Eigen::Vector3f(10, 0, 0)) > 0)
		{
			lambda = compute_insert(a, b, Eigen::Vector3f(0, 0, 0), Eigen::Vector3f(10, 0, 0));
			new_x.push_back(lambda*polyline_x[i] + (1 - lambda)*polyline_x[i + 1]);
			new_y.push_back(lambda*polyline_y[i] + (1 - lambda)*polyline_y[i + 1]);
			insert_list[1].push_back(new_x.size() - 1);
		}
		else if (compute_insert(a, b, Eigen::Vector3f(0, 0, 0), Eigen::Vector3f(0, -10, 0)) > 0)
		{
			lambda = compute_insert(a, b, Eigen::Vector3f(0, 0, 0), Eigen::Vector3f(0, -10, 0));
			new_x.push_back(lambda*polyline_x[i] + (1 - lambda)*polyline_x[i + 1]);
			new_y.push_back(lambda*polyline_y[i] + (1 - lambda)*polyline_y[i + 1]);
			insert_list[2].push_back(new_x.size() - 1);
		}
		else if (compute_insert(a, b, Eigen::Vector3f(0, 0, 0), Eigen::Vector3f(-10, 0, 0)) > 0)
		{
			lambda = compute_insert(a, b, Eigen::Vector3f(0, 0, 0), Eigen::Vector3f(-10, 0, 0));
			new_x.push_back(lambda*polyline_x[i] + (1 - lambda)*polyline_x[i + 1]);
			new_y.push_back(lambda*polyline_y[i] + (1 - lambda)*polyline_y[i + 1]);
			insert_list[3].push_back(new_x.size() - 1);
		}
		else if (compute_insert(a, b, Eigen::Vector3f(0, rectangle_max.y(), 0), Eigen::Vector3f(rectangle_max.x(), 0, 0)) > 0)
		{
			lambda = compute_insert(a, b, Eigen::Vector3f(0, rectangle_max.y(), 0), Eigen::Vector3f(rectangle_max.x(), 0, 0));
			new_x.push_back(lambda*polyline_x[i] + (1 - lambda)*polyline_x[i + 1]);
			new_y.push_back(lambda*polyline_y[i] + (1 - lambda)*polyline_y[i + 1]);
			insert_list[4].push_back(new_x.size() - 1);
		}
		else if (compute_insert(a, b, Eigen::Vector3f(rectangle_max.x(), 0, 0), Eigen::Vector3f(0, rectangle_min.y(), 0)) > 0)
		{
			lambda = compute_insert(a, b, Eigen::Vector3f(rectangle_max.x(), 0, 0), Eigen::Vector3f(0, rectangle_min.y(), 0));
			new_x.push_back(lambda*polyline_x[i] + (1 - lambda)*polyline_x[i + 1]);
			new_y.push_back(lambda*polyline_y[i] + (1 - lambda)*polyline_y[i + 1]);
			insert_list[5].push_back(new_x.size() - 1);
		}
		else if (compute_insert(a, b, Eigen::Vector3f(0, rectangle_min.y(), 0), Eigen::Vector3f(rectangle_min.x(), 0, 0)) > 0)
		{
			lambda = compute_insert(a, b, Eigen::Vector3f(0, rectangle_min.y(), 0), Eigen::Vector3f(rectangle_min.x(), 0, 0));
			new_x.push_back(lambda*polyline_x[i] + (1 - lambda)*polyline_x[i + 1]);
			new_y.push_back(lambda*polyline_y[i] + (1 - lambda)*polyline_y[i + 1]);
			insert_list[6].push_back(new_x.size() - 1);
		}
		else if (compute_insert(a, b, Eigen::Vector3f(rectangle_min.x(), 0, 0), Eigen::Vector3f(0, rectangle_max.y(), 0)) > 0)
		{
			lambda = compute_insert(a, b, Eigen::Vector3f(rectangle_min.x(), 0, 0), Eigen::Vector3f(0, rectangle_max.y(), 0));
			new_x.push_back(lambda*polyline_x[i] + (1 - lambda)*polyline_x[i + 1]);
			new_y.push_back(lambda*polyline_y[i] + (1 - lambda)*polyline_y[i + 1]);
			insert_list[7].push_back(new_x.size() - 1);
		}
	}

	new_x.push_back(polyline_x[polyline_x.size()-1]);
	new_y.push_back(polyline_y[polyline_y.size()-1]);

	polyline_x = new_x;
	polyline_y = new_y;

	for (int i = 0; i < insert_list.size(); i++)
	{
		std::vector<double > dis;
		std::vector<int> new_insert = insert_list[i];

		if (i < 4)
		{
			for (int j = 0; j < insert_list[i].size(); j++)
			{
				if (insert_list[i][j] == -1)
				{
					dis.push_back(0);
				}
				else
				{
					dis.push_back(polyline_x[insert_list[i][j]] * polyline_x[insert_list[i][j]] + polyline_y[insert_list[i][j]] * polyline_y[insert_list[i][j]]);
				}
			}
		}
		if (i == 4)
		{
			for (int j = 0; j < insert_list[i].size(); j++)
			{
				double x = polyline_x[insert_list[i][j]];
				double y = polyline_y[insert_list[i][j]] - rectangle_max.y();
				dis.push_back(x*x + y*y);
			}
		}
		if (i == 5)
		{
			for (int j = 0; j < insert_list[i].size(); j++)
			{
				double x = polyline_x[insert_list[i][j]] - rectangle_max.x();
				double y = polyline_y[insert_list[i][j]];
				dis.push_back(x*x + y*y);
			}
		}
		if (i == 6)
		{
			for (int j = 0; j < insert_list[i].size(); j++)
			{
				double x = polyline_x[insert_list[i][j]];
				double y = polyline_y[insert_list[i][j]] - rectangle_min.y();
				dis.push_back(x*x + y*y);
			}
		}
		if (i == 7)
		{
			for (int j = 0; j < insert_list[i].size(); j++)
			{
				double x = polyline_x[insert_list[i][j]] - rectangle_min.x();
				double y = polyline_y[insert_list[i][j]];
				dis.push_back(x*x + y*y);
			}
		}

		for (int k = 0; k < insert_list[i].size(); k++)
		{
			for (int mk = 0; mk < insert_list[i].size() - 1; mk++)
			{
				if (dis[mk] > dis[mk + 1])
				{
					int a = insert_list[i][mk];
					double b = dis[mk];
					insert_list[i][mk] = insert_list[i][mk + 1];
					insert_list[i][mk + 1] = a;

					dis[mk] = dis[mk + 1];
					dis[mk + 1] = b;
				}
			}
		}

	}
}

double BPolyline::compute_insert(Eigen::Vector3f& line11, Eigen::Vector3f& line12, Eigen::Vector3f& line21, Eigen::Vector3f& line22)
{
	Eigen::Matrix2d A;
	A << (line11.x() - line12.x()), -(line21.x() - line22.x()), (line11.y() - line12.y()), -(line21.y() - line22.y());
	Eigen::Vector2d b;
	b << line22.x() - line12.x(), line22.y() - line12.y();
	Eigen::Vector2d lambda;
	lambda = A.lu().solve(b);
	if (lambda(0) >= 0 && lambda(0) <= 1 && lambda(1) >= 0 && lambda(1) <= 1)
	{
		return lambda(0);
	}
	else
	{
		return -1;
	}

}

