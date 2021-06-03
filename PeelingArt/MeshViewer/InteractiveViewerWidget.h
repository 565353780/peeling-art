#ifndef INTERACTIVE_VIEWER_WIDGET
#define INTERACTIVE_VIEWER_WIDGET

#include "MeshViewerWidget.h"
#include "ANN\ANN.h"
#include <queue>

#include <Eigen/Dense>

#include "Algorithm/ARAP/ARAPParametricManger.h"
#include "dchart/Parafun.h"
#include "dchart/BijectivePara.h"
#include "Mesh/surface_manager.h"

#include "EasyToShow/EdgeToGraph.h"

class MyRemeshing;
class ARAP3D;
class BPolyline;
class CutMesh_Manager;
class MeshGeneration;

using namespace Eigen;

class InteractiveViewerWidget : public MeshViewerWidget
{
	Q_OBJECT
public:
	InteractiveViewerWidget(QWidget* parent = 0);
	InteractiveViewerWidget(QGLFormat& _fmt, QWidget* _parent);
	~InteractiveViewerWidget();

	void clearSelectedData()
	{
		handle_point.clear();
	};

	virtual void clearAllMesh()
	{
		draw_new_mesh = false;
		clearSelectedData();
		MeshViewerWidget::clearAllMesh();
	}

signals:
	void mouse_press_signal(Mesh::Point P);
	void mouse_move_signal(OpenMesh::Vec3d xy);
	void mouse_release_signal(Mesh::Point  P);
	void draw_from_out_signal();
	void setMouseMode_signal(int);

public:
	enum { TRANS, MOVE, T2_MODE, N_MODE, HANDLE, EYE };
	void setMouseMode(int mm);
	int mouseMode() const { return mouse_mode_; }

protected:
	virtual void mousePressEvent(QMouseEvent *_event);
	virtual void mouseReleaseEvent(QMouseEvent *_event);
	virtual void mouseMoveEvent(QMouseEvent *_event);
	virtual void wheelEvent(QWheelEvent* _event);
	int mouse_mode_;

protected:
	int	 pick_vertex(int x, int y);
	void pick_point(int x,int y);
	int find_vertex_using_selected_point();
	void move_point_based_lastVertex(int x,int y);	
	void buildIndex();
	ANNkd_tree* kdTree;

protected:
	double selectedPoint[3];
	int lastestVertex;

protected:
	void draw_different_model(int drawmode);
	void draw_selected_vertex();
	void draw_selected_face();
	void draw_selected_edge();
	void draw_handle_point();
	virtual void draw_scene(int drawmode);
	bool draw_new_mesh;
	

public:
    void  LoadFace(bool property_add);
    void  RelizeARAP3d_Ref(double lambda);

	
	
	void MyARAP3d();
		
	void Parameterization();
	
	void add_mesh_ori_info();
	void Transition_Mesh_To_Ref();
	void Transition_Ref_To_Mesh(); //更改参考为参数化的参考，对结果影响很大,never usr。

	void My_Deformation_preapre(int handle_id_in_fix);
	void My_Deformation(Mesh::Point P);
	void My_Deformotion_finally(); //re scaffold

	void SetEye();
	void GetCut();

	
	//将没有覆盖的部分分成不相交的块。
	void Change_cutmesh();
	void FillArea(Mesh::FaceHandle& fh, std::vector<int>& cut_face_type_list, int type);

	Mesh getRefMesh()
	{
		return ref_mesh;
	}
	Mesh getCurvedMesh(){
		return curved_mesh;
	}
	void setChangeMesh(Mesh mesh_, bool All_change_ok);

	BPolyline*					ref_boundary_;

	ARAP3D*						arap3d_solver_;
	MyRemeshing*				remeshing_;
	ARAPParametricManger		para_manger;
	Parafun						parafun_m;
	Surface_Manager				surface_manager;
	MeshGeneration*				def_initial_mg;

	bool						have_first_load;
	//标识显示参数化还是显示球面网格。
	bool						is_show_papamaterization;
	//标识是否进行了参数化。
	bool						NO_Papamaterization;
	//标识是否改变了网格。
	bool						is_finally_change_shape;
	bool                        is_deformation_prepare;
	bool						draw_plane_mode;


	int							papamaterization_times;
	double                      area_lambda;

	int							plane_Eye_id;

	int							User_Interaction_Time; //标记用户进行了多少次交互

	std::vector< CutMesh_Manager* > cut_mesh_mg;


protected:
	void ARAPresult_to_OBJ();
	void Parameterzation_to_OBJ();	
	void boundary_poly_to_file();

	double  overlode();

	double comput_mesh_area_face(Mesh& my_mesh, Mesh::FaceHandle face);
	void comput_cover_area(Mesh& my_mesh, double &area_r, double& area_g);

	void get_ref_boundary_poly(Mesh& my_mesh);
	void get_curved_boundary_poly(Mesh& my_mesh);
	
	void comput_in_out_time(Mesh& my_mesh, int& in_time, int& out_time);

	void Parameterization_3D(Mesh& plane_mesh);

private:
	std::vector<int>  handle_point;

public:
    int png_idx = 0;
    EdgeToGraph* graph;
    bool graph_init = false;
	
};

#endif
