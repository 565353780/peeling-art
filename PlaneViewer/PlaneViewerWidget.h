#pragma once
#include "PlaneViewer\childwidget.h"
#include "MeshViewer\MeshDefinition.h"
#include <Eigen/Dense>

#include "dchart/Parafun.h"

class BPolyline;
class QPainter;
class CutMesh_Manager;
class BSpline;
class QPushButton;

class PlaneViewerWidget :
	public Childwidget
{
	Q_OBJECT

public:
	PlaneViewerWidget(QWidget *parent = 0);
	~PlaneViewerWidget();

public:
	void setBoundary(BPolyline* bp);
	void setCutManager(std::vector<CutMesh_Manager*>& cut_list);
	//mouse_mode_
	enum {
		TRANS, CATSPLINE, FSPLINE, CHSPLINE, PICK, DEFHANDLE,  DEFMOVE, WEIGHT, T2_MODE, N_MODE
	};
	//spline_type_
	enum {
		CUT_SPLINE, FREE_SPLINE, SPLIT_SPLINE
	};
	//draw_type_
	enum {
		NOIC_PLANE, PLANE, CURVED
	};
	void setMouseMode(int mm);
	void changeRefMesh();

private:
	int mouse_mode_;

protected:
	virtual void mousePressEvent(QMouseEvent *_event);
	virtual void mouseReleaseEvent(QMouseEvent *_event);
	virtual void mouseMoveEvent(QMouseEvent *_event);
	virtual void wheelEvent(QWheelEvent* _event);

	virtual void draw_scene(int drawmode);

protected:
	void pick_point(int x, int y, Eigen::Vector3f& spoint);

private:
	bool nodo_spline_;
	std::vector<int> index_in_ori; //the index in original boundary index

private:
	void draw_Refmesh_And_Boundary();
	void draw_CurvedMesh_And_Boundary();
	void draw_All_Cut();

	void draw_Rectangel();
	void draw_Cut_Mesh();
	void draw_Spline();
	void draw_Handel();
	void draw_weight_bound();

	std::vector<std::vector<Vector3f> > all_curve_reduct_cut;

public:
	void select_interation_Part(Eigen::Vector3f& np);
	void select_interation_Part(int i, int j);
	void ClearAllPart();

	void ChangeCutPart();

public:
	void ComputerCutSpline(Eigen::Vector3f& bound_b, Eigen::Vector3f& bound_e);
	void ComputerSplitSpline(Eigen::Vector3f& bound_b, Eigen::Vector3f& bound_e);

	void ChangeShowSpline();

	void ComputeMaxH();

private:
	double compute_insert(Eigen::Vector3f& line1_b, Eigen::Vector3f& line1_e, Eigen::Vector3f& line2_b, Eigen::Vector3f& line2_e);
	//当点击样条第一个点时，初始化区域的两个端点和样条起点，以及样条的类型。
	void Inline_Cutspline_First(Eigen::Vector3f& np);
	void Inline_Freespline_First(Eigen::Vector3f& np);
	void Inline_Splitspline_First(Eigen::Vector3f& np);

	void Change_Spline_FirstPoint(Eigen::Vector3f& np, int s_type);
	void Change_Spline_FinallyPoint(Eigen::Vector3f& np);

	void swapSequence(std::vector<Eigen::Vector3f>& Pointlist); //调换样条的顺序

	void get_show_spline_choosed_control(Eigen::Vector3f& np);
	void get_boundary_closed_id(Eigen::Vector3f& np, int& id);

	
private:
	void SetMeshWeight();
	void SetWeightPolyline();

public:
	void ChangeRefOK();
	double computer_mesh_area(Mesh& my_mesh);

public:
	BPolyline* boundary_poly;
	BPolyline* curved_bopundary_poly; //save the curved mesh boundary, not change.
	
private:
	BSpline* show_spline;	
	std::vector<Eigen::Vector3f> Spline_ControlPoint;
	int spline_type;
	bool have_spline;
	int choose_spline_control_id;

	std::vector< BPolyline* > Weight_poly_list;
	std::vector< Eigen::Vector3f > Weight_bond_Point;
	
public:
	std::vector<CutMesh_Manager* > cut_Manager;
	CutMesh_Manager* user_select_part;
	bool is_user_select_part;

public:
	bool is_curved_inline;
	int planeview_draw_model;
	double max_curved_h;

public:
	//这里的网格是平面参考网格。
	Mesh mesh;
	Mesh curved_mesh;
	
	bool All_Change_OK;

private:
	int  FindMeshClosetid(Eigen::Vector3f np);
	void My_Deformation_prepare(int handle_id);
	void My_Deformotion(Mesh::Point P);
	void My_Deformotion_finally();

	std::set<int> deformation_fix_id;
	int handel_id;
	Parafun	parafun_m;  //for deformation
};

