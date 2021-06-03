#pragma once
#include"dchart/Parafun_bij.h"
#include"dchart/ScafData.h"
#include "MeshViewer/MeshDefinition.h"
#include <memory>

class BijectivePara
{
public:
	BijectivePara();
	~BijectivePara();

	void parameterization(Mesh& mesh_, std::vector<double>& face_energy);
	void load(const Eigen::MatrixXd&, const Eigen::MatrixXi&);
	void load(Mesh& mesh_);

public:
	ScafData scaf_data;
	std::shared_ptr<Parafun_bij> parafun_solver = nullptr;
	double convgence_con_rate = 1e-5;
	int MAX_ITER_NUM = 500;
};

