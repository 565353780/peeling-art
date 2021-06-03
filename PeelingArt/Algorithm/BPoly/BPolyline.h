#pragma once
#include <vector>
#include <Eigen/Dense>
#include <OpenMesh/Core/Geometry/VectorT.hh>

class BPolyline
{
public:
	BPolyline();
	~BPolyline();

	std::vector<double> polyline_x;
	std::vector<double> polyline_y;

	Eigen::Vector3d rectangle_max;
	Eigen::Vector3d rectangle_min;

	bool is_close;

public:
	bool IntPolyline(std::vector<Eigen::Vector3f>& polint_list);
	bool IntPolyline(std::vector<Eigen::Vector3f>& polint_list, bool is_close_);

	//get relationship between point and curve: in return 1£¬boundary return 0£¬out return -1.
	int PointInCurve(Eigen::Vector3f& p);

	//reduce some point in Polyline for easy comput;
	void ChangePolylineAndReduce();
	void ChangePolylineAndReduce(std::vector<int>& index_ori, std::vector<Eigen::Vector3f>& one_part);
	void ChangePolylineAndReduce(std::vector<int>& index_ori, std::vector<Eigen::Vector3f>& first_part, std::vector<Eigen::Vector3f>& secont_part);

	void Centralization(double factor);
	void Centralization(double center_x, double center_y);
	void Centralization(double center_x, double center_y, double factor);

	void FindCenter(double& center_x, double& center_y);

	double Id_Geodesic_Distance(int fidx, int sidx);

public:
	void Add_insert_point(std::vector< std::vector<int> >& insert_list);
	double compute_insert(Eigen::Vector3f& line11, Eigen::Vector3f& line12, Eigen::Vector3f& line21, Eigen::Vector3f& line22);
};

