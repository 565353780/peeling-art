#pragma once
#include <vector>
#include "MeshViewer\MeshDefinition.h"
#include <Eigen/Dense>
#include <Eigen/LU>

class BSpline
{
public:
	BSpline();
	~BSpline();

public:
	std::vector<double> time_points;
	std::vector<Eigen::Vector3f> Contral_Point;
	std::vector<Eigen::Vector3f> Contral_Point_ori;

	std::vector<Eigen::Vector3f> finally_contral_Point;
	std::vector<Eigen::Vector3f> finally_first_Point;
	std::vector<Eigen::Vector3f> finally_second_Point;

	std::vector<Eigen::Vector3f> draw_Contral_Point;
	std::vector<Eigen::Vector3f> draw_first_Point;
	std::vector<Eigen::Vector3f> draw_second_Point;

public:
	double B_base(std::vector<double>& tlist, double t, int k, int i);
	Eigen::Vector3f DB_base(double t);

	void CalculatePointlist();
	void CalculatePointlist_finally();
	void CalculatePoint(std::vector<Eigen::Vector3f>& control, double t, Eigen::Vector3f& r);

	void get_sharp_spline(double theta, Eigen::Vector3f& left_fix, Eigen::Vector3f& right_fix);

	void find_Colset(std::vector<Eigen::Vector3f>& new_control, std::vector<Eigen::Vector3f>& fix_infor);

	void GetSpline_from_data(std::vector<Eigen::Vector3f>& control_data);
	void GetSpline_from_data(std::vector<Eigen::Vector3f>& control_data, std::vector<double>& knock);

};

