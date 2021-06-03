#include <QMouseEvent>
#include <QLineEdit>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QtCore>
#include <QUrl>
#include <fstream>
#include <math.h>

#include "Algorithm/BPoly/BSpline.h"
#include "Algorithm/BPoly/BPolyline.h"
#include "Algorithm/MeshGeneration.h"
#include "Algorithm/MyRemeshing.h"
#include "Algorithm/ARAP/ARAP3D.h"
#include "Algorithm/CutMesh_Manager.h"
#include "Algorithm/ARAP/ARAPPara_3D.h"

#include "InteractiveViewerWidget.h"

#include "Algorithm/zip_up/zip_up.h"
#include "dchart/BijectivePara.h"
#include <QDebug>
#include <QFileDialog>
#include <QFileInfo>
#include <QTextStream>

const double PI = 3.1415926;
//外部面标记为ture， 内部面标记为false
OpenMesh::FPropHandleT<bool> is_singular;
OpenMesh::FPropHandleT<int> face_type;
//球面点对应的当前的平面参考
OpenMesh::VPropHandleT<Eigen::Vector3f> vert_ref;
//球面点对应的用户的平面参考
OpenMesh::VPropHandleT<Eigen::Vector3f> user_vert_ref;
//边的标识 外部标记为1，内部标记为-1，边界标记为0
OpenMesh::EPropHandleT<int> edge_type;
//点的标识 外部标记为1，内部标记为-1，边界标记为0
OpenMesh::VPropHandleT<int> vert_type;
//新的网格对应原始网格的序号。
OpenMesh::VPropHandleT<int> ori_vert_idx;

InteractiveViewerWidget::InteractiveViewerWidget(QWidget* parent /* = 0 */)
	:MeshViewerWidget(parent)
{
	draw_new_mesh = false;
	clearSelectedData();
	kdTree = NULL;
	have_first_load = false;
	is_show_papamaterization = false;
	NO_Papamaterization = true;
	papamaterization_times = 0;
	is_finally_change_shape = false;
	is_deformation_prepare = false;
	plane_Eye_id = -2;
	draw_plane_mode = false;
}

InteractiveViewerWidget::InteractiveViewerWidget(QGLFormat& _fmt, QWidget* _parent)
:MeshViewerWidget(_fmt, _parent)
{
	draw_new_mesh = false;
	clearSelectedData();
	kdTree = NULL;
	have_first_load = false;
	is_show_papamaterization = false;
	NO_Papamaterization = true;
	papamaterization_times = 0;
	is_finally_change_shape = false;
	is_deformation_prepare = false;
	plane_Eye_id = -2;
	draw_plane_mode = false;
}

InteractiveViewerWidget::~InteractiveViewerWidget()
{
	if(kdTree) delete kdTree;
}

//=====================================================================================
void InteractiveViewerWidget::setMouseMode(int mm)
{
	if(mouse_mode_ != T2_MODE)
	{
		mouse_mode_ = mm;
		if( TRANS != mouse_mode_ )
		{ buildIndex(); }
		emit setMouseMode_signal(mm);
	}
}

void InteractiveViewerWidget::mousePressEvent(QMouseEvent *_event)
{
	if (mouse_mode_ == TRANS)
	{
		MeshViewerWidget::mousePressEvent(_event);
	}
	else
	{
		if(mouse_mode_ != T2_MODE)
		{
			pick_point( _event->x(), _event->y() );
			if( mouse_mode_ == MOVE )
			{
				if (_event->button() == Qt::RightButton)
				{
					setMouseMode(TRANS);
				}
				else
				{
					pick_vertex(_event->x(), _event->y());
					int id_in_handle = -1;
					for (int i = 0; i < handle_point.size(); i++)
					{
						if (lastestVertex == handle_point[i])
						{
							id_in_handle = i;
							My_Deformation_preapre(id_in_handle);
						}					
					}
					if (id_in_handle == -1)
					{
						is_deformation_prepare = false;
					}
				}
			}
			else if (mouse_mode_ == HANDLE)
			{
				if (_event->button() == Qt::RightButton)
				{
					mouse_mode_ = TRANS;
				}
				else
				{
					int v_index = pick_vertex(_event->x(), _event->y());
					if (v_index >= 0 && mesh.property(vert_type, mesh.vertex_handle(v_index)) <= 0)
					{
						std::vector<int>::iterator it;
						it = std::find(handle_point.begin(), handle_point.end(), v_index);
						if (it == handle_point.end())
						{
							handle_point.push_back(v_index);
						}
						else
						{
							handle_point.erase(it);
						}						
					}
				}
			}
			else if (mouse_mode_ == EYE)
			{
				int v_index = pick_vertex(_event->x(), _event->y());
				if (v_index >= 0 && mesh.property(vert_type, mesh.vertex_handle(v_index)) <= 0)
				{
					plane_Eye_id = v_index;
				}
				else
				{
					plane_Eye_id = -1;
				}
				SetEye();
				mouse_mode_ = TRANS;
			}
		}
	}
	updateGL();
}

void InteractiveViewerWidget::mouseMoveEvent(QMouseEvent *_event)
{
	if(mouse_mode_ == TRANS)
	{
		MeshViewerWidget::mouseMoveEvent(_event);
	}
	else
	{
		if( mouse_mode_ != T2_MODE)
		{
			if( mouse_mode_ == MOVE && is_deformation_prepare)
			{
					move_point_based_lastVertex(_event->x(), _event->y());
					Mesh::Point P(selectedPoint[0], selectedPoint[1], selectedPoint[2]);
					My_Deformation(P);
					updateGL();
			}
		}
		else
		{
			
		}
		
	}
}

void InteractiveViewerWidget::mouseReleaseEvent(QMouseEvent *_event)
{
	if(mouse_mode_ == TRANS)
	{
		MeshViewerWidget::mouseMoveEvent(_event);
	}
	else
	{
		if(mouse_mode_ != T2_MODE )
		{
			if( mouse_mode_ == MOVE && is_deformation_prepare)
			{
				My_Deformotion_finally();
				updateGL();
			}
		}
		else
		{
		}
	}
	
}

void InteractiveViewerWidget::wheelEvent(QWheelEvent* _event)
{
	if(mouse_mode_ != N_MODE && mouse_mode_ != T2_MODE)
	{
		MeshViewerWidget::wheelEvent(_event);
	}
}

void InteractiveViewerWidget::pick_point(int x,int y)
{
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLdouble winX = double(x);
	GLdouble winY = double( height() - y );
	GLfloat winZ = 0.0;
	glReadPixels((int)winX, (int)winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
	gluUnProject(winX, winY, (GLdouble)winZ, &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &selectedPoint[0], &selectedPoint[1], &selectedPoint[2]);
}
int InteractiveViewerWidget::pick_vertex(int x, int y)
{
	int r = find_vertex_using_selected_point();
	lastestVertex = r;
	printf("Select Vertex : %d\n", r);
	return r;
}

void InteractiveViewerWidget::move_point_based_lastVertex(int x,int y)
{
	if(lastestVertex<0 || lastestVertex>=mesh.n_vertices())
	{
		return;
	}
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLdouble winX = 0.0;
	GLdouble winY = 0.0;
	GLdouble winZ = 0.0;
	OpenMesh::Vec3d p = mesh.point(mesh.vertex_handle(lastestVertex));
	gluProject(p[0], p[1], p[2],  &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &winX, &winY, &winZ);
	
	gluUnProject((GLdouble)(x), (GLdouble)( height() - y ), winZ,  &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &selectedPoint[0], &selectedPoint[1], &selectedPoint[2]);
}

int InteractiveViewerWidget::find_vertex_using_selected_point()
{
	ANNpoint tp = annAllocPt(3); tp[0] = selectedPoint[0]; tp[1] = selectedPoint[1]; tp[2] = selectedPoint[2];
	ANNidxArray nnIdx = new ANNidx[1]; ANNdistArray dists = new ANNdist[1];
	kdTree->annkSearch(tp, 1, nnIdx, dists);
	return nnIdx[0];
}

void InteractiveViewerWidget::buildIndex()
{
	if(mesh.n_vertices() == 0)
		return;

	Mesh::VertexIter v_it(mesh.vertices_begin());
	Mesh::VertexIter v_end(mesh.vertices_end());
	Mesh::Point p;
	unsigned nv = mesh.n_vertices();
	ANNpointArray dataPts = annAllocPts(nv, 3);
	int count = 0;
	for(; v_it != v_end; ++v_it)
	{
		p = mesh.point(v_it);
		dataPts[count][0] = p[0]; dataPts[count][1] = p[1]; dataPts[count][2] = p[2];
		++count;
	}

	if(kdTree) delete kdTree;
	kdTree = new ANNkd_tree(dataPts, nv, 3);
}

// this is other function
void InteractiveViewerWidget::draw_different_model(int drawmode)
{
	glViewport ( 0,0, width(),height());
	glMatrixMode( GL_PROJECTION );
	glLoadMatrixd( &ProjectionMatrix[0] );
	glMatrixMode( GL_MODELVIEW );
	glLoadMatrixd( &ModelViewMatrix[0] );

	emit draw_from_out_signal();

	{
		//draw select vertex, face, edge.
		glDisable(GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);

		glPointSize(1);

		draw_selected_vertex();
		draw_selected_face();
		draw_selected_edge();

		draw_handle_point();
	}

	if(draw_new_mesh)
	{
		draw_scene_mesh(drawmode);
	}
}

void InteractiveViewerWidget::draw_selected_vertex()
{
	if( have_first_load && 0)
	{
		Mesh::Point p;
		glColor3f(1.0, 0.0, 0.0);
		glPointSize(8);
		glBegin(GL_POINTS);

		for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
		{
			if (mesh.property(vert_type, v_it) == 0)
			{
				p = mesh.point(v_it.handle());
				glVertex3dv(p.data());
			}				
		}
		glEnd();
	}
	if (is_show_papamaterization && 0)
	{
		for (int i = 0; i < cut_mesh_mg.size(); i++)
		{
			Mesh::Point p;
			glColor3f(1.0, 0, 0);
			glPointSize(8);

			Mesh::ConstFaceVertexIter fv_it;
			Mesh::FaceHandle f_handle;
			Mesh& mymesh = cut_mesh_mg[i]->cut_mesh_a;

			glBegin(GL_POINTS);
			for (Mesh::VertexIter v_it = mymesh.vertices_begin(); v_it != mymesh.vertices_end(); v_it++)
			{
				if (mymesh.is_boundary(v_it.handle()))
				{
					p = mymesh.point(v_it.handle());
					glVertex3dv(p.data());
				}
			}
			glEnd();
		}
	}
}

void InteractiveViewerWidget::draw_selected_face()
{
	double cm = 255;
	if (have_first_load && !is_show_papamaterization)
	{
		Mesh::Point p;
		Mesh::ConstFaceVertexIter fv_it;
		Mesh::FaceHandle f_handle;
		for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
		{
			if (!mesh.property(is_singular, *f_it) && mesh.property(face_type, *f_it) != 2)
			{
				glColor3f(183 / cm, 205 / cm, 84 / cm);
				f_handle = f_it.handle();
				fv_it = mesh.fv_iter(f_handle);
				glBegin(GL_POLYGON);
				for (fv_it; fv_it; ++fv_it)
				{
					glVertex3dv(&mesh.point(fv_it)[0]);
				}
				glEnd();
			}
			if (mesh.property(is_singular, *f_it))
			{
				glColor3f(255/cm, 167/cm, 48/cm);
				f_handle = f_it.handle();
				fv_it = mesh.fv_iter(f_handle);
				glBegin(GL_POLYGON);
				for (fv_it; fv_it; ++fv_it)
				{
					glVertex3dv(&mesh.point(fv_it)[0]);
				}
				glEnd();
			}

			if (mesh.property(face_type, *f_it) == 2)
			{
				glColor3f(184 / cm, 131 / cm, 186 / cm);
				f_handle = f_it.handle();
				fv_it = mesh.fv_iter(f_handle);
				glBegin(GL_POLYGON);
				for (fv_it; fv_it; ++fv_it)
				{
					glVertex3dv(&mesh.point(fv_it)[0]);
				}
				glEnd();
				glColor3f(1.0, 0.5, 1.0);
			}
		}
	}
	if (is_show_papamaterization)
	{
		glDepthFunc(GL_LEQUAL);
		glClear(GL_DEPTH_BUFFER_BIT);

		double cn = 255;
		glColor3f(171 / cn, 197 / cn, 38 / cn);
		Mesh::Point p;
		Mesh::ConstFaceVertexIter fv_it;
		Mesh::FaceHandle f_handle;
		for (Mesh::FaceIter f_it = cut_mesh.faces_begin(); f_it != cut_mesh.faces_end(); f_it++)
		{
			if (!cut_mesh.property(is_singular, *f_it))
			{
				f_handle = f_it.handle();
				fv_it = cut_mesh.fv_iter(f_handle);
				glBegin(GL_POLYGON);
				for (fv_it; fv_it; ++fv_it)
				{
					glVertex3dv(&cut_mesh.point(fv_it)[0]);
				}
				glEnd();
			}
		}

		for (int i = 0; i < cut_mesh_mg.size(); i++)
		{
			double r = 0, g = 0, b = 0;
			cut_mesh_mg[i]->getColor(i + 1, r, g, b);
			glColor3f(r, g, b);

			Mesh::ConstFaceVertexIter fv_it;
			Mesh::FaceHandle f_handle;
			Mesh& mymesh = cut_mesh_mg[i]->cut_mesh_a;
			for (Mesh::FaceIter f_it = mymesh.faces_begin(); f_it != mymesh.faces_end(); f_it++)
			{
				f_handle = f_it.handle();
				fv_it = mymesh.fv_iter(f_handle);
				glBegin(GL_POLYGON);
				for (fv_it; fv_it; ++fv_it)
				{
					glVertex3dv(&mymesh.point(fv_it)[0]);
				}
				glEnd();
			}
		}
	}
}

void InteractiveViewerWidget::draw_selected_edge()
{
	if (is_show_papamaterization)
	{
		glColor3f(0.0, 0.0, 0.0);
		Mesh::EdgeHandle e_handle;
		Mesh::HalfedgeHandle he_handle;

		for (Mesh::EdgeIter e_it = cut_mesh.edges_begin(); e_it != cut_mesh.edges_end(); e_it++)
		{
			if (cut_mesh.property(edge_type, e_it) <= 0)
			{
				e_handle = e_it.handle();
				he_handle = cut_mesh.halfedge_handle(e_handle, 0);
				glBegin(GL_LINES);
				glVertex3dv(&cut_mesh.point(cut_mesh.from_vertex_handle(he_handle))[0]);
				glVertex3dv(&cut_mesh.point(cut_mesh.to_vertex_handle(he_handle))[0]);
				glEnd();
			}
		}
	}
}

void InteractiveViewerWidget::draw_scene(int drawmode)
{
	if (!mesh.n_vertices()) { return; }
	draw_different_model(drawmode);

	if (is_show_papamaterization || draw_plane_mode)
	{
		MeshViewerWidget::draw_scene_mesh(-100);
	}
	else if (!draw_new_mesh && !is_show_papamaterization)
	{
		MeshViewerWidget::draw_scene(drawmode);
	}
}


void InteractiveViewerWidget::draw_handle_point()
{
	if (handle_point.size() > 0)
	{
		Mesh::Point p;
		glColor3f(0.0, 0.0, 0.0);
		glPointSize(10);
		glBegin(GL_POINTS);

		glColor3f(0.0, 0.0, 0.0);
		for (unsigned int i = 0; i < handle_point.size(); ++i)
		{
			p = mesh.point(mesh.vertex_handle(handle_point[i]));
			glVertex3dv(p.data());
		}

		glEnd();
		glPointSize(1);
	}
}
//main process function
void InteractiveViewerWidget::MyARAP3d()
{
	mouse_mode_ = TRANS;

	if (mesh.n_vertices() < 1)
	{
		std::cout << "Must A Well Define Mesh" << std::endl;
		return;
	}

	arap3d_solver_ = new ARAP3D();
	remeshing_ = new MyRemeshing();
	arap3d_solver_->mesh_surface_ = &surface_manager;
	remeshing_->mesh_suface_ = &surface_manager;
	arap3d_solver_->special_weight = 15;

	double area_g = 1000;
	double area_r = 20;

	if (papamaterization_times == 0)
	{
		double mu = 0.7;
		arap3d_solver_->shrink_weight = 0.01; 
		//arap3d_solver_->shrink_weight = 0.1;
		double over_area = 0;
		double v_lambda = area_lambda;

		int max_iter_time = 101;

		for (int k = 1; k < max_iter_time; k++)
		{
			RelizeARAP3d_Ref(v_lambda);
			arap3d_solver_->mesh_ = mesh;
			arap3d_solver_->optimizeMesh();
			mesh = arap3d_solver_->mesh_;

			//over
			over_area = overlode();
			if (over_area > 0.5 && User_Interaction_Time == 0)
			{
				v_lambda = v_lambda*0.95;
			}

			if (k % 3 == 0)
			{
				remeshing_->My_Remeshing(mu, mesh);
				def_initial_mg->ReSetEye(mesh);
				mesh = def_initial_mg->mesh;
			}

			if (area_g < 0.1)
			{
				NO_Papamaterization = true;
				ARAPresult_to_OBJ();
				updateGL();
				return;
			}

			double old_arear = area_r;
			comput_cover_area(mesh, area_r, area_g);
			if (abs(area_r - old_arear) < 0.00001) break;

			updateGL();
		}
		remeshing_->My_Remeshing(mu, mesh);
		def_initial_mg->ReSetEye(mesh);
		mesh = def_initial_mg->mesh;

		NO_Papamaterization = true;
		area_lambda = v_lambda;
		is_show_papamaterization = false;

		Parameterization();
	}

	if (papamaterization_times > 0 && is_finally_change_shape)
	{
		int max_step_iter_time = 51;
		int one_step_iter_time = 0;
		int all_step_iter_time = 0;
		int step_time = 1;
		arap3d_solver_->shrink_weight = 0.1;

		while (true)
		{
			RelizeARAP3d_Ref(1);
			arap3d_solver_->mesh_ = mesh;
			arap3d_solver_->optimizeMesh();
			mesh = arap3d_solver_->mesh_;

			double old_arear = area_r;
			comput_cover_area(mesh, area_r, area_g);

			if (all_step_iter_time % 3 == 0 && all_step_iter_time > 0)
			{
				remeshing_->My_Remeshing(1.0, mesh);
			}
			def_initial_mg->ReSetEye(mesh);
			mesh = def_initial_mg->mesh;

			double myover_area = overlode();
			double a = myover_area < 0.06 || arap3d_solver_->shrink_weight > 25;
			double b = area_g < 0.1 || all_step_iter_time > 1000 || arap3d_solver_->shrink_weight > 8;
			if (a && b)
			{
				ARAPresult_to_OBJ();
				NO_Papamaterization = true;
				cut_mesh_mg.clear();

				papamaterization_times = 0;
				Parameterization();

				GetCut();
				updateGL();
				return;
			}

			one_step_iter_time++;
			all_step_iter_time++;

			int step_it = (10 * one_step_iter_time / arap3d_solver_->shrink_weight) - 2;
			if ((abs(area_r - old_arear) < 0.0005 && step_it > 2) || one_step_iter_time >= max_step_iter_time)
			{
				step_time++;
				double bl = 1.0 + 0.15;
				arap3d_solver_->shrink_weight *= bl;
				one_step_iter_time = 0;

				NO_Papamaterization = true;
				Parameterization();
			}

			updateGL();
		}
	}
	delete arap3d_solver_;
	delete remeshing_;
	plane_Eye_id = -2;
	updateGL();
}

void InteractiveViewerWidget::Parameterization()
{
	if (NO_Papamaterization && papamaterization_times == 0)
	{
		Transition_Mesh_To_Ref();

		//compute the parameterize of the S^m------------------------
		std::vector<double> face_energy;
		BijectivePara* bpara = new BijectivePara();
		bpara->load(ref_mesh);
		bpara->parameterization(ref_mesh, face_energy);
		delete bpara;

		//compute curved mesh------------------------
		curved_mesh = ref_mesh;
		if (!is_finally_change_shape)
		{
			for (Mesh::VertexIter v_it = curved_mesh.vertices_begin(); v_it != curved_mesh.vertices_end(); v_it++)
			{
				Mesh::Point p = Mesh::Point(curved_mesh.property(vert_ref, *v_it).x()*0.75, curved_mesh.property(vert_ref, *v_it).y()*0.75, 0);
				curved_mesh.set_point(*v_it, p);
			}
			Parameterization_3D(curved_mesh);
		}
		else
		{
			cut_mesh_mg.clear();
		}

		//using curved mesh may be lead Self-intersection, it is not solved.
//		get_ref_boundary_poly(ref_mesh);		
		get_curved_boundary_poly(curved_mesh);
				
		for (int ai = 0; ai < cut_mesh_mg.size(); ai++)
		{
			//compute the local parametric
			parafun_m = Parafun(cut_mesh_mg[ai]->cut_mesh_a, 1);
			parafun_m.init_Parametric();
			parafun_m.run_cm();
			parafun_m.ChangeVertRef(parafun_m.mesh, 0);

			para_manger.cut_mesh_ref = parafun_m.mesh;
			int cupy_n = para_manger.cutmeshBoundary_Inline();

			cut_mesh_mg[ai]->cut_mesh_b = para_manger.mesh_;

			cut_mesh_mg[ai]->cut_mesh_b_list.clear();
			cut_mesh_mg[ai]->is_show_list.clear();
			cut_mesh_mg[ai]->cut_b_endpoints_id.clear();
			

			for (int i = 0; i < cupy_n; i++)
			{
				//计算同一个区域的不同拷贝的参数化
				std::pair < int, int > bound_end_id;
//				para_manger.cutmeshBoundary_fix(ref_mesh, i, bound_end_id);
				para_manger.cutmeshBoundary_fix(curved_mesh, i, bound_end_id);
				parafun_m = Parafun(para_manger.mesh_, 3);
				parafun_m.init_Parametric(para_manger.fixed_id);
				parafun_m.run_cm();
				parafun_m.ChangeVertRef(parafun_m.mesh, -1);

				int in_time, out_time;
				comput_in_out_time(parafun_m.mesh, in_time, out_time);
				if (cut_mesh_mg[ai]->is_over_load > 0)
				{
					if (in_time < 0.3*parafun_m.mesh.n_vertices())
					{
						cut_mesh_mg[ai]->cut_mesh_b_list.push_back(parafun_m.mesh);
						cut_mesh_mg[ai]->cut_b_endpoints_id.push_back(bound_end_id);
						cut_mesh_mg[ai]->is_show_list.push_back(true);
					}
				}
				else
				{
					cut_mesh_mg[ai]->cut_mesh_b_list.push_back(parafun_m.mesh);
					cut_mesh_mg[ai]->cut_b_endpoints_id.push_back(bound_end_id);
					cut_mesh_mg[ai]->is_show_list.push_back(true);
				}

			}
			cut_mesh_mg[ai]->InlineBoundary();
		}
		papamaterization_times = 1;
		NO_Papamaterization = false;
		is_show_papamaterization = true;

		Parameterzation_to_OBJ();
	}

	if (NO_Papamaterization && papamaterization_times > 0)
	{
		Transition_Mesh_To_Ref();

		parafun_m = Parafun(ref_mesh, 1);
		parafun_m.init_Parametric();
		parafun_m.run_cm();
		parafun_m.ChangeVertRef(ref_mesh, para_manger.ref_fix_id);

		get_ref_boundary_poly(ref_mesh);

		cut_mesh_mg.clear();

		papamaterization_times++;
		NO_Papamaterization = false;
		is_show_papamaterization = false;
	}

	updateGL();
}




double InteractiveViewerWidget::overlode()
{
	mesh.update_face_normals();
	double overload_area = 0;
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		Mesh::Normal f_normal = mesh.normal(f_it.handle());
		Mesh::FaceVertexIter fv_it = mesh.fv_begin(*f_it); Mesh::Point p0 = mesh.point(*fv_it); 
		++fv_it;	Mesh::Point p1 = mesh.point(*fv_it); 
		++fv_it;	Mesh::Point p2 = mesh.point(*fv_it);
		Mesh::Point f_center = (p0 + p1 + p2) / 3;
		if (dot(f_normal, f_center) < 0)
		{
			OpenMesh::Vec3d v1 = p1 - p0;
			OpenMesh::Vec3d v2 = p2 - p0;
			double area_f = (cross(v1, v2)).norm() / 2;
			overload_area += area_f;
		}
	}
	return overload_area;
}

void  InteractiveViewerWidget::RelizeARAP3d_Ref(double lambda)
{
	arap3d_solver_->is_singular_face_.resize(mesh.n_faces());
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		arap3d_solver_->is_singular_face_[f_it.handle().idx()] = mesh.property(is_singular, *f_it);
	}
	arap3d_solver_->ref_map_.clear();
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		int i = v_it.handle().idx();
		Eigen::Vector3f refx = mesh.property(vert_ref, *v_it) * lambda;
		arap3d_solver_->ref_map_.insert(std::pair<int, Eigen::Vector3f>(i, refx));
	}
}

void  InteractiveViewerWidget::LoadFace()
{	
    if (1)
	{
		surface_manager.surface_type = Surface_Manager::SURFACE;

        QString filename=QFileDialog::getOpenFileName(0,tr(""),"./Samples","txt (*.txt);;SVG (*.svg)");
        QFileInfo info(filename);
        QFile file(info.absoluteFilePath());
        std::vector< Vector3f > handle_list;
        int handlesize=0;
        if(file.open(QFile::ReadOnly)){
            int index=0;
            QString line;
            QStringList data_list;
            QTextStream text_stream(&file);
            while(!text_stream.atEnd()){
                line=text_stream.readLine();
                data_list=line.split(" ");
                if(index==0){
                    if(data_list.size()==1){
                        handlesize=data_list[0].toInt();
                        handle_list.reserve(handlesize);
                        continue;
                    }
                }
                if(data_list.size()==3){
                    handle_list.push_back(Vector3f(data_list[0].toFloat(),data_list[1].toFloat(),data_list[2].toFloat()));
                }
                index++;
            }
            file.close();
        }
        if(handlesize==0 && handle_list.size()>0){
            handlesize=handle_list.size();
        }
/*old load file*/
//		std::ifstream filehandle("Mesh\\AMesh.txt", std::ios::in);//the input is by rhino
//		int handlesize;
//		filehandle >> handlesize;
//		std::vector< Vector3f > handle_list(handlesize);
//		for (int i = 0; i < handlesize; ++i)
//		{
//			double a, b, c;
//			filehandle >> a;
//			filehandle >> b;
//			filehandle >> c;
//			handle_list[i] = Vector3f(a, b, c);
//		}
//		filehandle.close();

        qDebug()<<"handlesize:"<<handlesize;
		def_initial_mg = new MeshGeneration();
		def_initial_mg->surface_manager = &surface_manager;

		def_initial_mg->meshFromUserPolyline(handle_list, ref_mesh);
//		def_initial_mg->meshFromUserPolyline_picture(handle_list, ref_mesh);

		ref_boundary_ = new BPolyline();
		ref_boundary_->IntPolyline(handle_list);
		ref_boundary_->Centralization(1);
		ref_boundary_->ChangePolylineAndReduce();

		mesh = def_initial_mg->mesh;
		
		initMesh();
		mesh_vector.push_back(mesh); mesh_vector_index = 0;

		have_first_load = true;
		area_lambda = 1;
	}
	draw_plane_mode = true;

	papamaterization_times = 0;
	User_Interaction_Time = 0;

	//进行直接编辑时需要把这项加上
//	ref_mesh.request_face_normals();
//	ref_mesh.request_vertex_normals();
//	ref_mesh.request_halfedge_normals();

	//进行不编辑直接做时需要加上这项
	is_finally_change_shape = true;
	papamaterization_times = 1;

	updateGL();
}

void InteractiveViewerWidget::add_mesh_ori_info()
{
	ref_mesh = mesh;
	cut_mesh = mesh;
	ref_mesh.add_property(ori_vert_idx);
	cut_mesh.add_property(ori_vert_idx);
	for (Mesh::VertexIter v_it = ref_mesh.vertices_begin(); v_it != ref_mesh.vertices_end(); v_it++)
	{
		ref_mesh.property(ori_vert_idx, *v_it) = v_it.handle().idx();
	}
	for (Mesh::VertexIter v_it = cut_mesh.vertices_begin(); v_it != cut_mesh.vertices_end(); v_it++)
	{
		cut_mesh.property(ori_vert_idx, *v_it) = v_it.handle().idx();
	}

}
void InteractiveViewerWidget::Transition_Mesh_To_Ref()
{
	add_mesh_ori_info();

	if (NO_Papamaterization && papamaterization_times == 0)
	{
		std::vector<bool> cut_face_list(cut_mesh.n_faces(), false);
		std::vector<bool> cut_add_face_list(cut_mesh.n_faces(), true);
		for (Mesh::FaceIter f_it = cut_mesh.faces_begin(); f_it != cut_mesh.faces_end(); f_it++)
		{
			std::vector<int> is_delate;
			is_delate.clear();
			for (Mesh::FaceVertexIter fv_it = cut_mesh.fv_begin(*f_it); fv_it != cut_mesh.fv_end(*f_it); fv_it++)
			{
				is_delate.push_back(cut_mesh.property(vert_type, *fv_it));
			}
			if (is_delate[0] < 1 && is_delate[1] < 1 && is_delate[2] < 1)
			{
				cut_face_list[f_it.handle().idx()] = true;
			}
		}
		for (Mesh::FaceIter f_it = cut_mesh.faces_begin(); f_it != cut_mesh.faces_end(); f_it++)
		{
			if (cut_face_list[f_it.handle().idx()] && cut_mesh.property(is_singular, *f_it))
			{
				std::vector<bool> is_delate;
				int around_face_type = 0;
				is_delate.clear();

				//120是最初设计的边界点的数目；
				Mesh::FaceVertexIter fv_it = mesh.fv_begin(*f_it);
				int id1 = fv_it.handle().idx() - 120; Mesh::Point p1 = mesh.point(*fv_it);
				++fv_it;  int id2 = fv_it.handle().idx() - 120; Mesh::Point p2 = mesh.point(*fv_it);
				++fv_it;  int id3 = fv_it.handle().idx() - 120; Mesh::Point p3 = mesh.point(*fv_it);

				std::vector<double> dis(3);
				dis[0] = ref_boundary_->Id_Geodesic_Distance(id1, id2);
				dis[1] = ref_boundary_->Id_Geodesic_Distance(id2, id3);
				dis[2] = ref_boundary_->Id_Geodesic_Distance(id1, id3);
				std::vector<double> dis_m(3);
				dis_m[0] = (p2 - p1).norm();
				dis_m[1] = (p3 - p2).norm();
				dis_m[2] = (p3 - p1).norm();

				int max_id = 0;
				if (dis[1] > dis[0]) max_id = 1;
				if (dis[2] > dis[max_id]) max_id = 2;

				if (dis[max_id]/ dis_m[max_id] < 3)
				{
					cut_add_face_list[f_it.handle().idx()] = false;
				}
			}
		}
		for (Mesh::FaceIter f_it = cut_mesh.faces_begin(); f_it != cut_mesh.faces_end(); f_it++)
		{
			if (cut_face_list[f_it.handle().idx()] && cut_add_face_list[f_it.handle().idx()])
			{
				cut_mesh.delete_face(f_it.handle());
			}
		}

		cut_mesh.garbage_collection();
		Change_cutmesh();
	}

	if (NO_Papamaterization)
	{
		for (Mesh::FaceIter f_it = ref_mesh.faces_begin(); f_it != ref_mesh.faces_end(); f_it++)
		{
			if (ref_mesh.property(is_singular, *f_it) == true)
			{
				ref_mesh.delete_face(f_it.handle());
			}
		}
		ref_mesh.garbage_collection();
		double mymin = 100;
		int min_idx = -1;
		for (Mesh::VertexIter v_it = ref_mesh.vertices_begin(); v_it != ref_mesh.vertices_end(); v_it++)
		{
			int i = v_it.handle().idx();
			Eigen::Vector3f refx = ref_mesh.property(vert_ref, *v_it);
			if (mymin > refx.norm())
			{
				mymin = refx.norm();
				min_idx = i;
			}
		}
		para_manger.ref_fix_id = min_idx;
		para_manger.area_lambda = area_lambda;

		if (papamaterization_times == 0)
		{
			cut_mesh = ref_mesh;
		}
	}
}


void InteractiveViewerWidget::Change_cutmesh()
{
	std::vector<int> cut_vert_type_list(cut_mesh.n_vertices(), 1);
	std::vector<int> cut_face_type_list(cut_mesh.n_faces(), 1);
	std::vector<double> cut_face_area(cut_mesh.n_faces(), 0);

	int type = -1;
	for (Mesh::FaceIter f_it = cut_mesh.faces_begin(); f_it != cut_mesh.faces_end(); f_it++)
	{
		int id = f_it.handle().idx();
		if (cut_face_type_list[id] > 0)
		{
			FillArea(f_it.handle(), cut_face_type_list, type);
			type--;
		}
	}

	for (Mesh::FaceIter f_it = cut_mesh.faces_begin(); f_it != cut_mesh.faces_end(); f_it++)
	{
		int id = f_it.handle().idx();
		cut_face_area[id] = comput_mesh_area_face(cut_mesh, f_it.handle());
		for (Mesh::FaceVertexIter fv_it = cut_mesh.fv_begin(*f_it); fv_it != cut_mesh.fv_end(*f_it); fv_it++)
		{
			int idv = fv_it.handle().idx();
			cut_vert_type_list[idv] = cut_face_type_list[id];
		}
	}

	int min_type = -(type + 1);
	std::vector< double > cut_part_area(min_type, 0);
	for (int i = 0; i < cut_face_type_list.size(); i++)
	{
		int f_type = -(cut_face_type_list[i] + 1);
		cut_part_area[f_type] += cut_face_area[i];
	}
	if (cut_mesh_mg.size() > 0)
	{
		for (int i = 0; i < cut_mesh_mg.size(); i++)
		{
			delete cut_mesh_mg[i];
		}
	}
	cut_mesh_mg.clear();

	Mesh c_mesh_a;
	int c_overload_a;
	for (int i = 0; i < min_type; i++)
	{
		if (abs(cut_part_area[i]) > 0.05)
		{
			c_mesh_a = cut_mesh;
			for (Mesh::VertexIter v_it = c_mesh_a.vertices_begin(); v_it != c_mesh_a.vertices_end(); v_it++)
			{
				int id = v_it.handle().idx();
				if (cut_vert_type_list[id] != (-1 - i))
				{
					c_mesh_a.delete_vertex(v_it.handle());
				}
			}
			c_mesh_a.garbage_collection();
			if (cut_part_area[i] > 0) c_overload_a = 1;
			if (cut_part_area[i] < 0) c_overload_a = -1;

			CutMesh_Manager* cm = new CutMesh_Manager();
			cm->cut_mesh_a = c_mesh_a;
			cm->is_over_load = c_overload_a;
			cut_mesh_mg.push_back(cm);
		}		
	}
}
void InteractiveViewerWidget::FillArea(Mesh::FaceHandle& fh, std::vector<int>& cut_face_type_list, int type)
{
	int id = fh.idx();
	if (cut_face_type_list[id] > 0)
	{
		cut_face_type_list[id] = type;
		for (Mesh::FaceFaceIter ff_it = cut_mesh.ff_begin(fh); ff_it != cut_mesh.ff_end(fh); ff_it++)
		{
			FillArea(ff_it.handle(), cut_face_type_list, type);
		}
	}
	return;
}

void InteractiveViewerWidget::comput_cover_area(Mesh& my_mesh, double &area_r, double& area_g)
{
	double my_area_red = 0;
	double my_area_green = 0;
	for (Mesh::FaceIter f_it = my_mesh.faces_begin(); f_it != my_mesh.faces_end(); f_it++)
	{
		Mesh::FaceVertexIter fv_it = my_mesh.fv_begin(*f_it); Mesh::Point p0 = my_mesh.point(*fv_it);
		++fv_it;   Mesh::Point p1 = my_mesh.point(*fv_it);
		++fv_it;   Mesh::Point p2 = my_mesh.point(*fv_it);

		OpenMesh::Vec3d v1 = p1 - p0;
		OpenMesh::Vec3d v2 = p2 - p0;
		double area_f = (OpenMesh::cross(v1, v2)).norm() / 2;

		if (!my_mesh.property(is_singular, *f_it))
		{
			my_area_red = my_area_red + area_f;
		}
		else
		{
			my_area_green = my_area_green + area_f;
		}
	}
	std::cout << "the area is:" << std::endl;
	std::cout << my_area_red << std::endl;

	area_r = my_area_red;
	area_g = my_area_green;
}

double InteractiveViewerWidget::comput_mesh_area_face(Mesh& my_mesh, Mesh::FaceHandle face)
{
	Mesh::Normal f_normal = my_mesh.normal(face);

	Mesh::FaceVertexIter fv_it = my_mesh.fv_begin(face); Mesh::Point p0 = my_mesh.point(*fv_it);
	++fv_it;   Mesh::Point p1 = my_mesh.point(*fv_it);
	++fv_it;   Mesh::Point p2 = my_mesh.point(*fv_it);

	Mesh::Point f_center = (p0 + p1 + p2) / 3;

	OpenMesh::Vec3d v1 = p1 - p0;
	OpenMesh::Vec3d v2 = p2 - p0;
	double area_f = (OpenMesh::cross(v1, v2)).norm() / 2;
	double sig = dot(f_normal, f_center) > 0 ? 1 : -1;
	return sig*area_f;
}

void  InteractiveViewerWidget::My_Deformation_preapre(int handle_id_in_fix)
{
	std::set<int> my_fix_id;
	my_fix_id.clear();
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		int i = v_it.handle().idx();
		if (mesh.is_boundary(v_it))
		{
			my_fix_id.insert(i);
		}
	}
	if (handle_point.size())
	{
		for (int i = 0; i < handle_point.size(); i++)
		{
			if (i != handle_id_in_fix)
			{
				my_fix_id.insert(handle_point[i]);
			}			
		}
	}

	std::vector<int> my_handle_id;
	my_handle_id.clear();
	my_handle_id.push_back(handle_point[handle_id_in_fix]);

	std::vector<int> my_face_in_out;
	my_face_in_out.clear();	
	vector<double> s_p00(mesh.n_faces()); vector<double> s_p01(mesh.n_faces()); vector<double> s_p10(mesh.n_faces()); vector<double> s_p11(mesh.n_faces());
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		int f_id = f_it->idx();
		if (mesh.property(is_singular, *f_it))
		{
			Mesh::FaceVertexIter fv_it = mesh.fv_begin(*f_it); int f0 = fv_it->idx();
			++fv_it;   int f1 = fv_it->idx();
			++fv_it;   int f2 = fv_it->idx();

			OpenMesh::Vec3d x_ = (mesh.point(mesh.vertex_handle(f1)) - mesh.point(mesh.vertex_handle(f0)));
			double x1_0 = x_.length();
			OpenMesh::Vec3d l_ = mesh.point(mesh.vertex_handle(f2)) - mesh.point(mesh.vertex_handle(f0));
			OpenMesh::Vec3d n = OpenMesh::cross(x_, l_).normalize();
			OpenMesh::Vec3d y_ = n % (1 / x1_0 * x_);

			double x2_0 = 1 / x1_0*l_ | x_;
			double y2_0 = l_ | y_;

			s_p00[f_id] = 1 / x1_0;
			s_p01[f_id] = -x2_0 / (x1_0*y2_0);
			s_p10[f_id] = 0;
			s_p11[f_id] = 1 / y2_0;

			my_face_in_out.push_back(-1);
		}
		else
		{
			Mesh::FaceVertexIter fv_it = mesh.fv_begin(*f_it);
			Eigen::Vector3f p0 = mesh.property(vert_ref, *fv_it);
			++fv_it;  Eigen::Vector3f p1 = mesh.property(vert_ref, *fv_it);
			++fv_it;  Eigen::Vector3f p2 = mesh.property(vert_ref, *fv_it);
			
			OpenMesh::Vec3d x_(p1(0) - p0(0), p1(1) - p0(1), p1(2)-p0(2)); // (mesh.point(mesh.vertex_handle(f1)) - mesh.point(mesh.vertex_handle(f0)));
			double x1_0 = x_.length();
			OpenMesh::Vec3d l_(p2(0) - p0(0), p2(1) - p0(1), p2(2) - p0(2)); //mesh.point(mesh.vertex_handle(f2)) - mesh.point(mesh.vertex_handle(f0));
			
			OpenMesh::Vec3d n = OpenMesh::cross(x_, l_).normalize();
			OpenMesh::Vec3d y_ = n % (1 / x1_0 * x_);
			double x2_0 = 1 / x1_0*l_ | x_;
			double y2_0 = l_ | y_;

			s_p00[f_id] = 1 / x1_0;
			s_p01[f_id] = -x2_0 / (x1_0*y2_0);
			s_p10[f_id] = 0;
			s_p11[f_id] = 1 / y2_0;

			my_face_in_out.push_back(1);
		}
	}

	parafun_m = Parafun(mesh, 2, s_p00,s_p01,s_p10,s_p11);
	parafun_m.init_Deformation(my_fix_id, my_handle_id, my_face_in_out);

	is_deformation_prepare = true;
}

void  InteractiveViewerWidget::My_Deformation(Mesh::Point P)
{
	std::vector<double> h_x, h_y;
	h_x.clear();
	h_y.clear();
	h_x.push_back(P[0]);
	h_y.push_back(P[1]);

	parafun_m.init_handle(h_x, h_y);
	bool isok = parafun_m.run_cm();
	if (isok)
	{
		parafun_m.ChangeVertPoint(mesh);
	}
	else
	{
		std::cout << "deformation is not ok" << std::endl;
	}
	
	updateGL();
}

void  InteractiveViewerWidget::My_Deformotion_finally()
{
	int n = 120;
	for (Mesh::VertexIter v_it = ref_mesh.vertices_begin(); v_it != ref_mesh.vertices_end(); v_it++)
	{
		int id = v_it.handle().idx();
		Mesh::VertexHandle mvh = mesh.vertex_handle(id + n);
		Mesh::Point p = mesh.point(mvh);
		ref_mesh.set_point(v_it.handle(), p);
	}

	int pn = ref_boundary_->polyline_x.size() - 1;
	for (int i = 0; i < pn; i++)
	{
		Mesh::Point p = ref_mesh.point(ref_mesh.vertex_handle(i));
		ref_boundary_->polyline_x[i] = p[0];
		ref_boundary_->polyline_y[i] = p[1];
	}
	ref_boundary_->polyline_x[pn] = ref_boundary_->polyline_x[0];
	ref_boundary_->polyline_y[pn] = ref_boundary_->polyline_y[0];

	def_initial_mg->mesh_scaffold_Inline(ref_boundary_, ref_mesh, mesh);

	buildIndex();

	is_deformation_prepare = false;
}

void InteractiveViewerWidget::get_ref_boundary_poly(Mesh& my_mesh)
{
	vector<Vector3f> boundary;
	boundary.clear();
	auto heit = my_mesh.halfedges_begin();
	while (!my_mesh.is_boundary(*heit))
		heit++;
	auto he_start = *heit;
	auto he_it = he_start;
	do
	{
		he_it = my_mesh.next_halfedge_handle(he_it);
		boundary.push_back(my_mesh.property(vert_ref, my_mesh.to_vertex_handle(he_it)));
	} while (he_it != he_start);

	ref_boundary_->IntPolyline(boundary, true);
}
void InteractiveViewerWidget::get_curved_boundary_poly(Mesh& my_mesh)
{
	vector<Vector3f> boundary;
	boundary.clear();
	auto heit = my_mesh.halfedges_begin();
	while (!my_mesh.is_boundary(*heit))
		heit++;
	auto he_start = *heit;
	auto he_it = he_start;
	do
	{
		he_it = my_mesh.next_halfedge_handle(he_it);
		boundary.push_back(my_mesh.property(vert_ref, my_mesh.to_vertex_handle(he_it)));
	} while (he_it != he_start);

	ref_boundary_->IntPolyline(boundary, true);
}

void InteractiveViewerWidget::comput_in_out_time(Mesh& my_mesh, int& in_time, int& out_time)
{
	in_time = 0;
	out_time = 0;
	int nv = my_mesh.n_vertices();
	std::vector<Eigen::Vector3f> vert_list(nv);
	std::vector<int> in_out_list(nv);
	for (Mesh::VertexIter v_it = my_mesh.vertices_begin(); v_it != my_mesh.vertices_end(); v_it++)
	{
		vert_list[v_it.handle().idx()] = my_mesh.property(vert_ref, *v_it);
	}

#pragma omp parallel for
	for (int i = 0; i < nv; i++)
	{
		in_out_list[i] = ref_boundary_->PointInCurve(vert_list[i]);
	}

	for (int i = 0; i < nv; i++)
	{
		if (in_out_list[i] > 0) in_time++;
		if (in_out_list[i] < 0) out_time++;
	} 
}

void InteractiveViewerWidget::setChangeMesh(Mesh mesh_, bool All_change_ok)
{
	draw_plane_mode = true;
	is_show_papamaterization = false;
	User_Interaction_Time++;

	if (def_initial_mg != NULL)
	{
		delete def_initial_mg;
	}

	mesh = mesh_;
	ref_mesh = mesh_;
	
	def_initial_mg = new MeshGeneration();
	def_initial_mg->surface_manager = &surface_manager;
	def_initial_mg->mesh_scaffold_Inline(ref_boundary_, ref_mesh, mesh);

	is_finally_change_shape = All_change_ok;

	OpenMesh::Vec3d c(0, 0, 0);
	set_scene_pos(c, 3);
	updateGL();
}

void InteractiveViewerWidget::SetEye()
{
	draw_plane_mode = false;

	def_initial_mg->mesh = mesh;
	def_initial_mg->PlaneMesh_To_Surface(plane_Eye_id);
	mesh = def_initial_mg->mesh;

	if (!is_finally_change_shape)
	{
		area_lambda = 1;		
		papamaterization_times = 0;		
	}
	else
	{
		papamaterization_times++;
	}

	OpenMesh::Vec3d c(0, 0, 0);
	set_scene_pos(c, 1.25);
	handle_point.clear();
	updateGL();
}

void InteractiveViewerWidget::GetCut()
{
	QString meshname = "liuhaoFin";

	QString meshpath = meshname + ".obj";

	QString postfix = "v1s";

	ZipUp zip_tool;
	bool read_OK = OpenMesh::IO::read_mesh(zip_tool.mesh_, "liuhaoFin.obj");

	zip_tool.getVertex_list();
	zip_tool.getBoundaryIndex();
	zip_tool.getBoundaryAngle();

	zip_tool.Prim();

	zip_tool.evaluateComponent();
	zip_tool.cleanComponent(0.1f);
//	zip_tool.cleanComponent(0.05f);

	cout << zip_tool.tree_root_->degree_ << ' ' << zip_tool.tree_root_->sons_.size() << endl;


	zip_tool.writeComponent("Interaction\\" + meshname.toStdString() + "_component_" + postfix.toStdString() + ".txt");
	zip_tool.writeBoundaryStructure("Interaction\\" + meshname.toStdString() + "_structure_" + postfix.toStdString() + ".txt");
	zip_tool.writeTreeEdge("Interaction\\" + meshname.toStdString() + "_tedge_" + postfix.toStdString() + ".txt");
	zip_tool.writeBoundaryIndex("Interaction\\" + meshname.toStdString() + "_bindex_" + postfix.toStdString() + ".txt");
}

void InteractiveViewerWidget::Parameterization_3D(Mesh& plane_mesh)
{
	ARAPPara_3D* solve_height = new ARAPPara_3D();
	solve_height->mesh_surface = cut_mesh;
	solve_height->mesh_plane = plane_mesh;
	solve_height->optimizeMesh();
	curved_mesh = solve_height->mesh_plane;
	for (Mesh::VertexIter v_it = curved_mesh.vertices_begin(); v_it != curved_mesh.vertices_end(); v_it++)
	{
		Mesh::Point p = curved_mesh.point(*v_it);
		curved_mesh.property(vert_ref, *v_it) = Vector3f(p[0], p[1], 0);
	}

	delete solve_height;
}

void InteractiveViewerWidget::Parameterzation_to_OBJ()
{
	if (!NO_Papamaterization)
	{
		for (Mesh::VertexIter v_it = ref_mesh.vertices_begin(); v_it != ref_mesh.vertices_end(); v_it++)
		{
			Mesh::Point p = Mesh::Point(ref_mesh.property(vert_ref, *v_it).x(), ref_mesh.property(vert_ref, *v_it).y(), 0);
			ref_mesh.set_point(*v_it, p);
		}

		OpenMesh::IO::write_mesh(ref_mesh, "Interaction\\Acut.obj");

		for (Mesh::VertexIter v_it = curved_mesh.vertices_begin(); v_it != curved_mesh.vertices_end(); v_it++)
		{
			Mesh::Point p = Mesh::Point(curved_mesh.property(vert_ref, *v_it).x(), curved_mesh.property(vert_ref, *v_it).y(), 0);
			curved_mesh.set_point(*v_it, p);
		}

		OpenMesh::IO::write_mesh(curved_mesh, "Interaction\\Ccut.obj");
	}
}
void InteractiveViewerWidget::ARAPresult_to_OBJ()
{
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		if (mesh.property(is_singular, *f_it) == true)
		{
			mesh.delete_face(f_it.handle());
		}
	}
	mesh.garbage_collection();
	OpenMesh::IO::write_mesh(mesh, "liuhaoFin.obj");
	OpenMesh::IO::write_mesh(mesh, "Interaction\\liuhaoFin.obj");

	boundary_poly_to_file();
}

void InteractiveViewerWidget::boundary_poly_to_file()
{
	vector<Mesh::Point> boundary;
	boundary.clear();
	auto heit = mesh.halfedges_begin();
	while (!mesh.is_boundary(*heit))
		heit++;
	auto he_start = *heit;
	auto he_it = he_start;
	do
	{
		he_it = mesh.next_halfedge_handle(he_it);
		boundary.push_back(mesh.point( mesh.to_vertex_handle(he_it)));
	} while (he_it != he_start);

	std::ofstream file("Aboundary.txt");
	file << boundary.size() << std::endl;
	for (int i = 0; i < boundary.size(); i++)
	{
		file << boundary[i][0] << " " << boundary[i][1] << " " << boundary[i][2] << std::endl;
	}
}

#pragma region Not Used
void InteractiveViewerWidget::Transition_Ref_To_Mesh()
{
	if (!NO_Papamaterization)
	{
		for (Mesh::VertexIter v_it = ref_mesh.vertices_begin(); v_it != ref_mesh.vertices_end(); v_it++)
		{
			Mesh::VertexHandle vh = mesh.vertex_handle(ref_mesh.property(ori_vert_idx, *v_it));
			mesh.property(vert_ref, vh) = ref_mesh.property(vert_ref, *v_it);
		}
	}
}
#pragma endregion
