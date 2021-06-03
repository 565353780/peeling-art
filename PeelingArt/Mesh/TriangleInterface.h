#pragma once
#include<Eigen/Core>

#include "triangle.h" 

void triangulate(const Eigen::MatrixXd &bnd_pts, const Eigen::MatrixXi &E, const double& area_threshold,Eigen::MatrixXd & pts,Eigen::MatrixXi & FV);
void triangluate(const Eigen::MatrixXd &bnd_pts, const Eigen::MatrixXi &E, const Eigen::MatrixXd &hole, Eigen::MatrixXd & pts, Eigen::MatrixXi & FV);
void triangluate_bij(const Eigen::MatrixXd &bnd_pts, const Eigen::MatrixXi &E, const Eigen::MatrixXd &hole, Eigen::MatrixXd & pts, Eigen::MatrixXi & FV);


