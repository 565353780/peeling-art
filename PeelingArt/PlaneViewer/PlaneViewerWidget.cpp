#include "PlaneViewerWidget.h"
#include "Algorithm/BPoly/BPolyline.h"
#include "Algorithm/BPoly/BSpline.h"
#include <QPainter>
#include <QMouseEvent>
#include <QPushButton>
#include <fstream>


#include "Mesh/CutMesh_Manager.h"
#include "Mesh/MeshGeneration.h"

extern OpenMesh::FPropHandleT<bool> is_singular;
extern OpenMesh::FPropHandleT<int> face_type;

PlaneViewerWidget::PlaneViewerWidget(QWidget *parent /* = 0 */)
	:Childwidget(parent)
{
	is_user_select_part = false;
	have_spline = false;
	show_spline = NULL;
	choose_spline_control_id = -1;
	nodo_spline_ = true;
	deformation_fix_id.clear();
	handel_id = -1;
	Weight_bond_Point.clear();
	Weight_poly_list.clear();
	all_curve_reduct_cut.clear();
	planeview_draw_model = PLANE;

	All_Change_OK = false;
	double maxh = 0;
}


PlaneViewerWidget::~PlaneViewerWidget()
{
}

void PlaneViewerWidget::draw_scene(int drawmode)
{
	Childwidget::draw_scene(drawmode);

	glDepthFunc(GL_LEQUAL);
	glClear(GL_DEPTH_BUFFER_BIT);
	draw_Rectangel();

	if (planeview_draw_model == CURVED)
	{		
		draw_CurvedMesh_And_Boundary();
		draw_All_Cut();
	}
	else
	{
		draw_Refmesh_And_Boundary();
	}
	

	if (is_user_select_part)
	{
		glColor3f(1.0, 0.0, 0.0);
		glPointSize(4);
		BPolyline* bp = user_select_part->cut_b_boundary[user_select_part->user_select];
		for (int k = 0; k < bp->polyline_x.size() - 1; k++)
		{
			glBegin(GL_LINES);
			glVertex3d(bp->polyline_x[k], bp->polyline_y[k], 1e-5);
			glVertex3d(bp->polyline_x[k + 1], bp->polyline_y[k + 1], 1e-5);
			glEnd();
		}
	}
	else if (mouse_mode_!= WEIGHT)
	{
		draw_Cut_Mesh();
	}

	draw_Spline();

	draw_Handel();

	draw_weight_bound();
}

void PlaneViewerWidget::mousePressEvent(QMouseEvent *_event)
{
		if (mouse_mode_ == TRANS)
		{
			Childwidget::mousePressEvent(_event);
		}
		else if (mouse_mode_ == CATSPLINE && is_user_select_part)
		{
			if (_event->button() == Qt::RightButton)
			{
				Eigen::Vector3f p1(boundary_poly->polyline_x[index_in_ori[0]], boundary_poly->polyline_y[index_in_ori[0]], 0);
				Eigen::Vector3f p2(boundary_poly->polyline_x[index_in_ori[2]], boundary_poly->polyline_y[index_in_ori[2]], 0);
				ComputerCutSpline(p1, p2);
				ChangeShowSpline();
				
				nodo_spline_ = true;
				setMouseMode(TRANS);
			}
			else
			{
				Eigen::Vector3f np;
				pick_point(_event->x(), _event->y(), np);
				if (nodo_spline_)
				{
					if (!have_spline)
					{
						Inline_Cutspline_First(np);
						Spline_ControlPoint.push_back(np);
						ChangeShowSpline();
					}					
				}
				else
				{
					Spline_ControlPoint.push_back(np);
					ChangeShowSpline();
				}
			}		
		}
		else if (mouse_mode_ == CATSPLINE && !is_user_select_part)
		{
			if (_event->button() == Qt::RightButton)
			{
				swapSequence(Spline_ControlPoint);
				Eigen::Vector3f p1(boundary_poly->polyline_x[index_in_ori[0]], boundary_poly->polyline_y[index_in_ori[0]], 0);
				Eigen::Vector3f p2(boundary_poly->polyline_x[index_in_ori[2]], boundary_poly->polyline_y[index_in_ori[2]], 0);
				ComputerSplitSpline(p1, p2);
				ChangeShowSpline();

				nodo_spline_ = true;
				setMouseMode(TRANS);
			}
			else
			{
				Eigen::Vector3f np;
				pick_point(_event->x(), _event->y(), np);
				if (nodo_spline_)
				{
					if (!have_spline)
					{
						Inline_Splitspline_First(np);
						Spline_ControlPoint.push_back(np);
						ChangeShowSpline();
					}
				}
				else
				{
					Spline_ControlPoint.push_back(np);
					ChangeShowSpline();
				}
			}
		}
		else if (mouse_mode_ == FSPLINE && is_user_select_part)
		{
			if (_event->button() == Qt::RightButton)
			{
				Change_Spline_FinallyPoint(Spline_ControlPoint[Spline_ControlPoint.size() - 1]);
				ChangeCutPart();
				ChangeShowSpline();

				nodo_spline_ = true;
				setMouseMode(TRANS);
			}
			else
			{
				Eigen::Vector3f np;
				pick_point(_event->x(), _event->y(), np);
				if (nodo_spline_)
				{
					if (!have_spline)
					{
						Inline_Freespline_First(np);
						Spline_ControlPoint.push_back(np);
						ChangeShowSpline();
					}
				}
				else
				{
					Spline_ControlPoint.push_back(np);
					ChangeShowSpline();
				}
			}
		}
		else if (mouse_mode_ == CHSPLINE)
		{
			if (_event->button() == Qt::RightButton)
			{
				nodo_spline_ = true;
				mouse_mode_ = TRANS;
			}
			else
			{
				if (Spline_ControlPoint.size() == 0)
				{
					mouse_mode_ = TRANS;
				}
				else
				{
					Eigen::Vector3f np;
					pick_point(_event->x(), _event->y(), np);
					get_show_spline_choosed_control(np);
				}		
			}
		}
		else if (mouse_mode_ == PICK)
		{
			Eigen::Vector3f np;
			pick_point(_event->x(), _event->y(), np);
			select_interation_Part(np);
			mouse_mode_ = TRANS;
		}
		else if (mouse_mode_ == DEFHANDLE)
		{
			if (_event->button() == Qt::RightButton)
			{
				mouse_mode_ = TRANS;			
			}
			else
			{
				Eigen::Vector3f np;
				pick_point(_event->x(), _event->y(), np);
				int id = FindMeshClosetid(np);
				bool is_in = false;
				for (std::set< int>::iterator iter = deformation_fix_id.begin(); iter != deformation_fix_id.end(); iter++)
				{
					if (*iter == id)
					{
						deformation_fix_id.erase(iter);
						is_in = true;
						break;
					}
				}
				if (!is_in && id >= 0)
				{
					deformation_fix_id.insert(id);
				}
			}			
		}
		else if (mouse_mode_ == DEFMOVE)
		{
			if (_event->button() == Qt::RightButton)
			{
				mouse_mode_ = TRANS;
				deformation_fix_id.clear();
			}
			else
			{
				Eigen::Vector3f np;
				pick_point(_event->x(), _event->y(), np);
				handel_id = FindMeshClosetid(np);
				if (handel_id >= 0) My_Deformation_prepare(handel_id);
			}
		}
		else if (mouse_mode_ == WEIGHT)
		{
			if (_event->button() == Qt::RightButton)
			{
				mouse_mode_ = TRANS;
				SetWeightPolyline();
			}
			else
			{
				Eigen::Vector3f np;
				pick_point(_event->x(), _event->y(), np);
				Weight_bond_Point.push_back(np);
			}
		}

		updateGL();
}

void PlaneViewerWidget::mouseMoveEvent(QMouseEvent *_event)
{
	if (mouse_mode_ == TRANS)
	{
		Childwidget::mouseMoveEvent(_event);
	}
	else if ((mouse_mode_ == CATSPLINE || mouse_mode_ == FSPLINE) && _event->button() == Qt::RightButton)
	{
		Childwidget::mouseMoveEvent(_event);
	}
	else if (mouse_mode_ == CHSPLINE)
	{
		Eigen::Vector3f np;
		pick_point(_event->x(), _event->y(), np);

		Spline_ControlPoint[choose_spline_control_id] = np;

		if (spline_type == CUT_SPLINE)
		{
			Change_Spline_FirstPoint(Spline_ControlPoint[0], spline_type);
			Eigen::Vector3f p1(boundary_poly->polyline_x[index_in_ori[0]], boundary_poly->polyline_y[index_in_ori[0]], 0);
			Eigen::Vector3f p2(boundary_poly->polyline_x[index_in_ori[2]], boundary_poly->polyline_y[index_in_ori[2]], 0);
			ComputerCutSpline(p1, p2);
			ChangeShowSpline();
		}
		else if (spline_type == FREE_SPLINE)
		{
			Change_Spline_FirstPoint(Spline_ControlPoint[0], spline_type);
			Change_Spline_FinallyPoint(Spline_ControlPoint[Spline_ControlPoint.size() - 1]);
			ChangeCutPart();
			ChangeShowSpline();
		}
		else if (spline_type == SPLIT_SPLINE)
		{
			Eigen::Vector3f p1(boundary_poly->polyline_x[index_in_ori[0]], boundary_poly->polyline_y[index_in_ori[0]], 0);
			Eigen::Vector3f p2(boundary_poly->polyline_x[index_in_ori[2]], boundary_poly->polyline_y[index_in_ori[2]], 0);
			ComputerSplitSpline(p1, p2);
			ChangeShowSpline();
		}
	}
	//is not used in program
	else if (mouse_mode_ == DEFMOVE)
	{
		Eigen::Vector3f np;
		pick_point(_event->x(), _event->y(), np);
		Mesh::Point P(np.x(), np.y(), 0);
		My_Deformotion(P);
	}

	updateGL();
}

void PlaneViewerWidget::mouseReleaseEvent(QMouseEvent *_event)
{
	if (mouse_mode_ == TRANS)
	{

	}
	else if (mouse_mode_ == CHSPLINE)
	{
		if (spline_type == CUT_SPLINE)
		{
			Change_Spline_FirstPoint(Spline_ControlPoint[0], spline_type);
			Eigen::Vector3f p1(boundary_poly->polyline_x[index_in_ori[0]], boundary_poly->polyline_y[index_in_ori[0]], 0);
			Eigen::Vector3f p2(boundary_poly->polyline_x[index_in_ori[2]], boundary_poly->polyline_y[index_in_ori[2]], 0);
			ComputerCutSpline(p1, p2);
			ChangeShowSpline();
		}
		else if (spline_type == FREE_SPLINE)
		{
			Change_Spline_FirstPoint(Spline_ControlPoint[0], spline_type);
			Change_Spline_FinallyPoint(Spline_ControlPoint[Spline_ControlPoint.size() - 1]);
			ChangeCutPart();
			ChangeShowSpline();
		}
		else if(spline_type == SPLIT_SPLINE)
		{
			Eigen::Vector3f p1(boundary_poly->polyline_x[index_in_ori[0]], boundary_poly->polyline_y[index_in_ori[0]], 0);
			Eigen::Vector3f p2(boundary_poly->polyline_x[index_in_ori[2]], boundary_poly->polyline_y[index_in_ori[2]], 0);
			ComputerSplitSpline(p1, p2);
			ChangeShowSpline();
		}
	}
	else if (mouse_mode_ == DEFMOVE)
	{
		My_Deformotion_finally();
		handel_id = -1;
		updateGL();
	}

	updateGL();

}

void PlaneViewerWidget::wheelEvent(QWheelEvent* _event)
{
	if (mouse_mode_ != N_MODE && mouse_mode_ != T2_MODE)
	{
		Childwidget::wheelEvent(_event);
	}
}

void PlaneViewerWidget::pick_point(int x, int y, Eigen::Vector3f& spoint)
{
	GLdouble object_x, object_y, object_z;
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	GLdouble winX = double(x);
	GLdouble winY = double(height() - y);
	GLfloat winZ = 0.0;
	glReadBuffer(GL_BACK);
	glReadPixels((int)winX, (int)winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);

	gluUnProject((GLdouble)winX, (GLdouble)winY, (GLdouble)winZ, &ModelViewMatrix[0], &ProjectionMatrix[0], viewport, &object_x, &object_y, &object_z);
	spoint = Eigen::Vector3f(object_x, object_y, object_z);
}

void PlaneViewerWidget::draw_All_Cut()
{
	glColor3f(0.0, 0.0, 0.0);
	glLineWidth(3);
	for (int a = 0; a < all_curve_reduct_cut.size(); a++)
	{
		for (int i = 0; i < all_curve_reduct_cut[a].size() - 1; i++)
		{
			glBegin(GL_LINES);
			glVertex3d(all_curve_reduct_cut[a][i].x(), all_curve_reduct_cut[a][i].y(), 3e-4);
			glVertex3d(all_curve_reduct_cut[a][i + 1].x(), all_curve_reduct_cut[a][i + 1].y(), 3e-4);
			glEnd();
		}
	}
}

void PlaneViewerWidget::draw_CurvedMesh_And_Boundary()
{
	glColor3f(0.0, 0.0, 0.0);
	glLineWidth(3);
	glPointSize(4);
	for (int i = 0; i < curved_bopundary_poly->polyline_x.size() - 1; i++)
	{
		glBegin(GL_LINES);
		glVertex3d(curved_bopundary_poly->polyline_x[i], curved_bopundary_poly->polyline_y[i], 3e-4);
		glVertex3d(curved_bopundary_poly->polyline_x[i + 1], curved_bopundary_poly->polyline_y[i + 1], 3e-4);
		glEnd();
	}

	double cn = 255;
	Eigen::Vector3f p_color = Vector3f(184 / cn, 131 / cn, 186 / cn);
	Eigen::Vector3f w_color = Vector3f(1.0, 1.0, 1.0);
	Eigen::Vector3f m_color;
	Mesh::ConstFaceVertexIter fv_it;
	Mesh::FaceHandle f_handle;
	for (Mesh::FaceIter f_it = curved_mesh.faces_begin(); f_it != curved_mesh.faces_end(); f_it++)
	{
		f_handle = f_it.handle();
		fv_it = curved_mesh.fv_iter(f_handle);
		glBegin(GL_POLYGON);
		for (fv_it; fv_it; ++fv_it)
		{
			Mesh::Point p = curved_mesh.point(fv_it);
			double h = p[2] / (max_curved_h);
			m_color = h*p_color + (1 - h)*w_color;
			glColor3f(m_color.x(), m_color.y(), m_color.z());
			glVertex3d(p[0], p[1], 5e-5);
		}
		glEnd();
	}
}
void PlaneViewerWidget::draw_Refmesh_And_Boundary()
{
	glColor3f(0.0, 0.0, 0.0);
	glPointSize(4);
	for (int i = 0; i < boundary_poly->polyline_x.size() - 1; i++)
	{
		glBegin(GL_LINES);
		glVertex3d(boundary_poly->polyline_x[i], boundary_poly->polyline_y[i], 2e-5);
		glVertex3d(boundary_poly->polyline_x[i + 1], boundary_poly->polyline_y[i + 1], 2e-5);
		glEnd();
	}

	double cn = 255;
	glColor3f(171 / cn, 197 / cn, 38 / cn);
	Mesh::ConstFaceVertexIter fv_it;
	Mesh::FaceHandle f_handle;
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		f_handle = f_it.handle();
		fv_it = mesh.fv_iter(f_handle);
		glBegin(GL_POLYGON);
		for (fv_it; fv_it; ++fv_it)
		{
			Mesh::Point p = mesh.point(fv_it);
			glVertex3d(p[0], p[1], 2e-5);
		}
		glEnd();
	}
}

void PlaneViewerWidget::draw_Cut_Mesh()
{
	for (int i = 0; i < cut_Manager.size(); i++)
	{
		for (int j = 0; j < cut_Manager[i]->cut_mesh_b_list.size(); j++)
		{
			if (cut_Manager[i]->is_show_list[j])
			{
				glColor3f(1.0, 0.0, 0.0);
				glPointSize(4);
				BPolyline* bp = cut_Manager[i]->cut_b_boundary[j];
				for (int k = 0; k < bp->polyline_x.size() - 1; k++)
				{
					glBegin(GL_LINES);
					glVertex3d(bp->polyline_x[k], bp->polyline_y[k], 0);
					glVertex3d(bp->polyline_x[k + 1], bp->polyline_y[k + 1], 0);
					glEnd();
				}

				//glColor3f(1.0, 1.0, 0.0);
				//glPointSize(5);
				//glBegin(GL_LINES);
				//Eigen::Vector3f end1 = cut_Manager[i]->cut_b_endpoints[j].first;
				//glVertex3d(end1.x(), end1.y(), 1e-4);
				//Eigen::Vector3f end2 = cut_Manager[i]->cut_b_endpoints[j].second;
				//glVertex3d(end2.x(), end2.y(), 1e-4);
				//glEnd();

				double r=0, g=0, b=0;
				cut_Manager[i]->getColor(i + 1, r, g, b);
				//if (cut_Manager[i]->is_over_load > 0)
				//	cut_Manager[i]->HSV_to_RGB(r, g, b, 120, 0.15 + 0.85*i / cn, 0.85);
				//else
				//	cut_Manager[i]->HSV_to_RGB(r, g, b, 240, 0.15 + 0.85*i / cn, 0.85);
				glColor3f(r, g, b);
				Mesh::ConstFaceVertexIter fv_it;
				Mesh::FaceHandle f_handle;
				Mesh& mymesh = cut_Manager[i]->cut_mesh_b_list[j];
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
}

void PlaneViewerWidget::draw_Spline()
{
	glColor3f(1.0, 0, 1.0);
	glPointSize(10);
	glBegin(GL_POINTS);
	for (int i = 0; i < Spline_ControlPoint.size(); i++)
	{
		glVertex3d(Spline_ControlPoint[i].x(), Spline_ControlPoint[i].y(), 2e-4);
	}
	glEnd();
	glBegin(GL_LINES);
	if (show_spline!= NULL)
	{
		for (int j = 0; j < show_spline->draw_Contral_Point.size() - 1; j++)
		{			
			Eigen::Vector3f s = show_spline->draw_Contral_Point[j];
			glVertex3d(s.x(), s.y(), 1e-4);
			s = show_spline->draw_Contral_Point[j + 1];
			glVertex3d(s.x(), s.y(), 1e-4);
		}
	}
	glEnd();
	glColor3f(0.0, 1.0, 1.0);

//	if (have_spline && (spline_type == CUT_SPLINE || spline_type == SPLIT_SPLINE) )
	if (have_spline && (spline_type == CUT_SPLINE))
	{
		glColor3f(0.0, 1.0, 1.0);
		glPointSize(10);

		glBegin(GL_POINTS);
		Eigen::Vector3f s = show_spline->finally_contral_Point[show_spline->finally_contral_Point.size() - 1];
		glVertex3d(s.x(), s.y(), 1e-4);

		for (int j = 0; j < show_spline->finally_first_Point.size(); j++)
		{
			Eigen::Vector3f s = show_spline->finally_first_Point[j];
			glVertex3d(s.x(), s.y(), 1e-4);
		}
		for (int j = 0; j < show_spline->finally_second_Point.size(); j++)
		{
			Eigen::Vector3f s = show_spline->finally_second_Point[j];
			glVertex3d(s.x(), s.y(), 1e-4);
		}
		glEnd();

		glBegin(GL_LINES);
		for (int j = 0; j < show_spline->draw_first_Point.size() - 1; j++)
		{
			Eigen::Vector3f s = show_spline->draw_first_Point[j];
			glVertex3d(s.x(), s.y(), 1e-4);
			s = show_spline->draw_first_Point[j + 1];
			glVertex3d(s.x(), s.y(), 1e-4);
		}
		glEnd();
		glBegin(GL_LINES);
		for (int j = 0; j < show_spline->draw_second_Point.size() - 1; j++)
		{
			Eigen::Vector3f s = show_spline->draw_second_Point[j];
			glVertex3d(s.x(), s.y(), 1e-4);
			s = show_spline->draw_second_Point[j + 1];
			glVertex3d(s.x(), s.y(), 1e-4);
		}
		glEnd();
	}
}

void PlaneViewerWidget::draw_Rectangel()
{
	glColor3f(1.0, 1.0, 1.0);
	glPointSize(1);
	glBegin(GL_QUADS);
	glVertex3d(-10.0, 10.0, -1e-5);
	glVertex3d(10.0, 10.0, -1e-5);
	glVertex3d(10.0, -10.0, -1e-5);
	glVertex3d(-10.0, -10.0, -1e-5);
	glEnd();
}

void PlaneViewerWidget::draw_Handel()
{
	if (deformation_fix_id.size() > 0)
	{
		Mesh::Point p;
		glColor3f(0.0, 0.0, 0.0);
		glPointSize(10);
		glBegin(GL_POINTS);

		std::cout << deformation_fix_id.size() << std::endl;

		for (std::set< int>::iterator iter = deformation_fix_id.begin(); iter != deformation_fix_id.end(); iter++)
		{
			p = mesh.point(mesh.vertex_handle(*iter));
			glVertex3d(p[0], p[1], 0.001);
		}

		glEnd();
	}
}

void PlaneViewerWidget::draw_weight_bound()
{
	glPointSize(5);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < Weight_bond_Point.size(); i++)
	{
		Eigen::Vector3f s = Weight_bond_Point[i];
		glVertex3d(s.x(), s.y(), 2e-4);
	}
	glEnd();

	if (Weight_bond_Point.size() > 1)
	{
		glBegin(GL_LINES);
		for (int i = 0; i < Weight_bond_Point.size() - 1; i++)
		{
			Eigen::Vector3f s = Weight_bond_Point[i];
			glVertex3d(s.x(), s.y(), 2e-4);
			s = Weight_bond_Point[i+1];
			glVertex3d(s.x(), s.y(), 2e-4);
		}
		glEnd();
	}

	glColor3f(0.0, 0.0, 0.0);
	glPointSize(6);
	for (int i = 0; i < Weight_poly_list.size(); i++)
	{
		for (int j = 0; j < Weight_poly_list[i]->polyline_x.size() - 1; j++)
		{
			glBegin(GL_LINES);
			glVertex3d(Weight_poly_list[i]->polyline_x[j], Weight_poly_list[i]->polyline_y[j], 1e-4);
			glVertex3d(Weight_poly_list[i]->polyline_x[j + 1], Weight_poly_list[i]->polyline_y[j + 1], 1e-4);
			glEnd();
		}
	}

}


void PlaneViewerWidget::ChangeShowSpline()
{
	if (have_spline)
	{
		show_spline->GetSpline_from_data(Spline_ControlPoint);
		show_spline->CalculatePointlist();
	}
	else if (Spline_ControlPoint.size() > 1 )
	{
		show_spline = new BSpline();
		show_spline->GetSpline_from_data(Spline_ControlPoint);
		show_spline->CalculatePointlist();
	}

}

void PlaneViewerWidget::ComputerCutSpline(Eigen::Vector3f& bound_b, Eigen::Vector3f& bound_e)
{
	std::vector<double> knock;
	knock.push_back(0);
	knock.push_back(0);
	knock.push_back(0);

	std::vector<Eigen::Vector3f> control_p;
	control_p.clear();
	double n_fin = 0;
	for (int i = 0; i < Spline_ControlPoint.size() - 1; i++)
	{
		double lambda = compute_insert(Spline_ControlPoint[i], Spline_ControlPoint[i + 1], bound_b, bound_e);
		if (lambda < 0)
		{
			control_p.push_back(Spline_ControlPoint[i]);
		}
		else
		{
			n_fin = 1-lambda;
			
			control_p.push_back(Spline_ControlPoint[i]);
			control_p.push_back(lambda*Spline_ControlPoint[i] + (1 - lambda)*Spline_ControlPoint[i + 1]);
			break;
		}
	}
	if (control_p.size() == 2)
	{
		control_p.resize(4);
		control_p[3] = control_p[1];
		control_p[1] = 0.66*control_p[0] + 0.33*control_p[3];
		control_p[2] = 0.33*control_p[0] + 0.66*control_p[3];
	}
	if (control_p.size() == 3)
	{
		control_p.resize(4);
		control_p[3] = control_p[2];
		control_p[2] = 0.5*control_p[1] + 0.5*control_p[3];
		control_p[1] = 0.5*control_p[0] + 0.5*control_p[1];
	}

	for (int i = 0; i <= control_p.size() - 4; i++)
	{
		knock.push_back(i);
	}

	n_fin = n_fin + knock[knock.size() - 1];
	knock.push_back(n_fin);
	knock.push_back(n_fin);
	knock.push_back(n_fin);
	knock.push_back(n_fin);
	
	show_spline->GetSpline_from_data(control_p, knock);
	show_spline->get_sharp_spline(0.13, bound_b, bound_e);
	show_spline->CalculatePointlist_finally();

	ChangeCutPart();
}

void PlaneViewerWidget::ComputerSplitSpline(Eigen::Vector3f& bound_b, Eigen::Vector3f& bound_e)
{
	std::vector<double> knock;
	knock.push_back(0);
	knock.push_back(0);
	knock.push_back(0);

	std::vector<Eigen::Vector3f> control_p;
	control_p.clear();
	double n_fin = 0;
	for (int i = 0; i < Spline_ControlPoint.size() - 1; i++)
	{
		double lambda = compute_insert(Spline_ControlPoint[i], Spline_ControlPoint[i + 1], bound_b, bound_e);
		if (lambda < 0)
		{
			control_p.push_back(Spline_ControlPoint[i]);
		}
		else
		{
			n_fin = 1 - lambda;

			control_p.push_back(Spline_ControlPoint[i]);
			control_p.push_back(lambda*Spline_ControlPoint[i] + (1 - lambda)*Spline_ControlPoint[i + 1]);
			break;
		}
	}
	if (control_p.size() == 2)
	{
		control_p.resize(4);
		control_p[3] = control_p[1];
		control_p[1] = 0.66*control_p[0] + control_p[3];
		control_p[2] = 0.33*control_p[0] + control_p[3];
	}
	if (control_p.size() == 3)
	{
		control_p.resize(4);
		control_p[3] = control_p[2];
		control_p[2] = 0.5*control_p[1] + 0.5*control_p[3];
		control_p[1] = 0.5*control_p[0] + 0.5*control_p[1];
	}

	for (int i = 0; i <= control_p.size() - 4; i++)
	{
		knock.push_back(i);
	}

	n_fin = n_fin + knock[knock.size() - 1];
	knock.push_back(n_fin);
	knock.push_back(n_fin);
	knock.push_back(n_fin);
	knock.push_back(n_fin);

	show_spline->GetSpline_from_data(control_p, knock);
	show_spline->get_sharp_spline(0.05, bound_b, bound_e);
	show_spline->CalculatePointlist_finally();

	have_spline = true;
}

double PlaneViewerWidget::compute_insert(Eigen::Vector3f& line11, Eigen::Vector3f& line12, Eigen::Vector3f& line21, Eigen::Vector3f& line22)
{
	Eigen::Matrix2d A;
	A << (line11.x() - line12.x()), -(line21.x() - line22.x()), (line11.y() - line12.y()), -(line21.y() - line22.y());
	Eigen::Vector2d b;
	b << line22.x() - line12.x(), line22.y() - line12.y();
	Eigen::Vector2d lambda;
	lambda = A.lu().solve(b);
	if (lambda(0) >= 0 && lambda(0) <= 1)
	{
		return lambda(0);
	}
	else
	{
		return -1;
	}

}

void PlaneViewerWidget::Inline_Cutspline_First(Eigen::Vector3f& np)
{
	Eigen::Vector3f end1 = user_select_part->cut_b_endpoints[user_select_part->user_select].first;
	Eigen::Vector3f end2 = user_select_part->cut_b_endpoints[user_select_part->user_select].second;

	double mdistance = 100, fdistance = 100, sdistance=100;
	int midx = -1, fidx = -1, sidx = -1;
	for (int i = 0; i < boundary_poly->polyline_x.size()-1; i++)
	{
		double dis = (np.x() - boundary_poly->polyline_x[i])*(np.x() - boundary_poly->polyline_x[i]) + (np.y() - boundary_poly->polyline_y[i])*(np.y() - boundary_poly->polyline_y[i]);
		if (dis < mdistance)
		{
			mdistance = dis;
			midx = i;
		}
		dis = (end1.x() - boundary_poly->polyline_x[i])*(end1.x() - boundary_poly->polyline_x[i]) + (end1.y() - boundary_poly->polyline_y[i])*(end1.y() - boundary_poly->polyline_y[i]);
		if (dis < fdistance)
		{
			fdistance = dis;
			fidx = i;
		}
		dis = (end2.x() - boundary_poly->polyline_x[i])*(end2.x() - boundary_poly->polyline_x[i]) + (end2.y() - boundary_poly->polyline_y[i])*(end2.y() - boundary_poly->polyline_y[i]);
		if (dis < sdistance)
		{
			sdistance = dis;
			sidx = i;
		}
	}

	np = Eigen::Vector3f(boundary_poly->polyline_x[midx], boundary_poly->polyline_y[midx], 0);

	int nn = boundary_poly->polyline_x.size()-1;
	int d1 = ((sidx - fidx) + nn) % nn;
	int d2 = ((fidx - sidx) + nn) % nn;

	if (d1 > d2)
	{
		double a = fidx;
		fidx = sidx;
		sidx = a;
	}

	index_in_ori.resize(3);
	index_in_ori[0] = fidx;
	index_in_ori[1] = midx;
	index_in_ori[2] = sidx;

	spline_type = CUT_SPLINE;
	nodo_spline_ = false;
}

void PlaneViewerWidget::Inline_Splitspline_First(Eigen::Vector3f& np)
{
	int sp_id = -1;
	get_boundary_closed_id(np, sp_id);

	int nn = boundary_poly->polyline_x.size() - 1;
	
	index_in_ori.resize(3);
	index_in_ori[0] = (sp_id-1)%nn;
	index_in_ori[1] = sp_id;
	index_in_ori[2] = (sp_id+1)%nn;

	spline_type = SPLIT_SPLINE;
	nodo_spline_ = false;
}

void PlaneViewerWidget::Inline_Freespline_First(Eigen::Vector3f& np)
{
	double mdistance = 100;
	int midx = -1;
	for (int i = 0; i < boundary_poly->polyline_x.size() - 1; i++)
	{
		double dis = (np.x() - boundary_poly->polyline_x[i])*(np.x() - boundary_poly->polyline_x[i]) + (np.y() - boundary_poly->polyline_y[i])*(np.y() - boundary_poly->polyline_y[i]);
		if (dis < mdistance)
		{
			mdistance = dis;
			midx = i;
		}
	}
	np = Eigen::Vector3f(boundary_poly->polyline_x[midx], boundary_poly->polyline_y[midx], 0);

	index_in_ori.resize(3);	
	index_in_ori[0] = midx;

	spline_type = FREE_SPLINE;
	nodo_spline_ = false;
}

void PlaneViewerWidget::Change_Spline_FinallyPoint(Eigen::Vector3f& np)
{
	double mdistance = 100;
	double midx = -1;
	for (int i = 0; i < boundary_poly->polyline_x.size() - 1; i++)
	{
		double dis = (np.x() - boundary_poly->polyline_x[i])*(np.x() - boundary_poly->polyline_x[i]) + (np.y() - boundary_poly->polyline_y[i])*(np.y() - boundary_poly->polyline_y[i]);
		if (dis < mdistance)
		{
			mdistance = dis;
			midx = i;
		}
	}
	np = Eigen::Vector3f(boundary_poly->polyline_x[midx], boundary_poly->polyline_y[midx], 0);

	if (spline_type == FREE_SPLINE )
	{
		index_in_ori[2] = midx;

		int fidx = index_in_ori[0];
		int sidx = index_in_ori[2];
		int nn = boundary_poly->polyline_x.size() - 1;
		int d1 = ((sidx - fidx) + nn) % nn;
		int d2 = ((fidx - sidx) + nn) % nn;

		if (d1 > d2)
		{
			index_in_ori[2] = fidx;
			index_in_ori[0] = sidx;
			swapSequence(Spline_ControlPoint);
		}
	}


}

void PlaneViewerWidget::Change_Spline_FirstPoint(Eigen::Vector3f& np, int s_type)
{
	double mdistance = 100;
	int midx = -1;
	for (int i = 0; i < boundary_poly->polyline_x.size() - 1; i++)
	{
		double dis = (np.x() - boundary_poly->polyline_x[i])*(np.x() - boundary_poly->polyline_x[i]) + (np.y() - boundary_poly->polyline_y[i])*(np.y() - boundary_poly->polyline_y[i]);
		if (dis < mdistance)
		{
			mdistance = dis;
			midx = i;
		}
	}
	np = Eigen::Vector3f(boundary_poly->polyline_x[midx], boundary_poly->polyline_y[midx], 0);
	if (s_type == FREE_SPLINE)
	{
		index_in_ori[0] = midx;
	}
	if (s_type == CUT_SPLINE || s_type == SPLIT_SPLINE)
	{
		index_in_ori[1] = midx;
	}

}


void PlaneViewerWidget::ChangeCutPart()
{
	for (int j = 0; j < user_select_part->is_show_list.size(); j++)
	{
		user_select_part->is_show_list[j] = false;
	}
	is_user_select_part = false;
	have_spline = true;
}

void PlaneViewerWidget::select_interation_Part(Eigen::Vector3f& np)
{
	for (int i = 0; i < cut_Manager.size(); i++)
	{
		for (int j = 0; j < cut_Manager[i]->cut_mesh_b_list.size(); j++)
		{
			if (cut_Manager[i]->is_show_list[j])
			{
				BPolyline* bp = cut_Manager[i]->cut_b_boundary[j];
				int a = bp->PointInCurve(np);
				if (a > -1)
					select_interation_Part(i, j);
			}
		}
	}
}
void PlaneViewerWidget::select_interation_Part(int i, int j)
{
	user_select_part = cut_Manager[i];
	user_select_part->user_select = j;
	is_user_select_part = true;
}
void PlaneViewerWidget::ClearAllPart()
{
	for (int i = 0; i < cut_Manager.size(); i++)
	{
		for (int j = 0; j < cut_Manager[i]->is_show_list.size(); j++)
		{
			cut_Manager[i]->is_show_list[j] = false;
		}
	}

	cut_Manager.clear();

	updateGL();
}

void PlaneViewerWidget::swapSequence(std::vector<Eigen::Vector3f>& Pointlist)
{
	std::vector<Eigen::Vector3f> newlist;
	newlist.clear();
	for (std::vector<Eigen::Vector3f>::const_reverse_iterator i = Pointlist.rbegin(); i != Pointlist.rend(); ++i)
	{
		newlist.push_back(*i);
	}
	for (int i = 0; i < newlist.size(); i++)
	{
		Pointlist[i] = newlist[i];
	}
}

void PlaneViewerWidget::get_show_spline_choosed_control(Eigen::Vector3f& np)
{
	double dis = 1000;
	for (int i = 0; i < Spline_ControlPoint.size(); i++)
	{
		if ((np - Spline_ControlPoint[i]).norm() < dis)
		{
			dis = (np - Spline_ControlPoint[i]).norm();
			choose_spline_control_id = i;
		}
	}
	if (dis > 0.1)
	{
		choose_spline_control_id = -1;
	}
}

void PlaneViewerWidget::get_boundary_closed_id(Eigen::Vector3f& np, int& id)
{
	double square_mindis = 1000;
	double min_id = -1;
	for (int i = 0; i < boundary_poly->polyline_x.size() - 1; i++)
	{
		double dis = (np.x() - boundary_poly->polyline_x[i])*(np.x() - boundary_poly->polyline_x[i]) + (np.y() - boundary_poly->polyline_y[i])*(np.y() - boundary_poly->polyline_y[i]);
		if (dis < square_mindis)
		{
			square_mindis = dis;
			min_id = i;
		}
	}
	id = min_id;
}


//Slot
void PlaneViewerWidget::setCutManager(std::vector<CutMesh_Manager*>& cut_list)
{
	cut_Manager = cut_list;
}

void PlaneViewerWidget::setBoundary(BPolyline* bp)
{
	boundary_poly = bp;
	setMouseMode(TRANS);
	planeview_draw_model = NOIC_PLANE;

	OpenMesh::Vec3d bbMax(bp->rectangle_max.x(), bp->rectangle_max.y(), 0);
	OpenMesh::Vec3d bbMin(bp->rectangle_min.x(), bp->rectangle_min.y(), 0);
	set_scene_pos((bbMin + bbMax)*0.5, (bbMin - bbMax).norm()*0.35);
}

void PlaneViewerWidget::setMouseMode(int mm)
{
	if (mouse_mode_ != T2_MODE)
	{
		mouse_mode_ = mm;
	}
}

void PlaneViewerWidget::changeRefMesh()
{

	std::cout << spline_type << std::endl;

	if (show_spline != NULL)
	{
		if (spline_type == CUT_SPLINE || spline_type == SPLIT_SPLINE)
		{
			std::vector<Vector3f> one_cut = show_spline->draw_Contral_Point;
			all_curve_reduct_cut.push_back(one_cut);
			boundary_poly->ChangePolylineAndReduce(index_in_ori, show_spline->draw_first_Point, show_spline->draw_second_Point);
		}
		else
		{		
			boundary_poly->ChangePolylineAndReduce(index_in_ori, show_spline->draw_Contral_Point);
		}
	}

	MeshGeneration* mg = new MeshGeneration();
	mg->meshFromChangePolyline(boundary_poly, mesh);
	delete mg;


	Spline_ControlPoint.clear();
	delete show_spline;
	show_spline = NULL;
	have_spline = false;
	is_user_select_part = false;
	nodo_spline_ = true;

	setMouseMode(TRANS);

	updateGL();
}





void PlaneViewerWidget::SetMeshWeight()
{
	for (int wi = 0; wi < Weight_poly_list.size(); wi++)
	{
		std::vector<Vector3f> center_list;
		BPolyline* bp = Weight_poly_list[wi];
		std::vector<Mesh::Point> f_p(3);
		Mesh::Point f_center;
		for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
		{
			Mesh::FaceVertexIter fv_it = mesh.fv_begin(*f_it);
			f_p[0] = mesh.point(*fv_it); fv_it++;
			f_p[1] = mesh.point(*fv_it); fv_it++;
			f_p[2] = mesh.point(*fv_it);
			f_center = (f_p[0] + f_p[1] + f_p[2]) / 3;
			center_list.push_back(Eigen::Vector3f(f_center[0], f_center[1], 0));
		}

		int nf = mesh.n_faces();
		std::vector<int> in_out_list(nf);
#pragma omp parallel for
		for (int i = 0; i < nf; i++)
		{
			in_out_list[i] = bp->PointInCurve(center_list[i]);
		}
		
		for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
		{
			int id = f_it->idx();
			if (in_out_list[id] >= 0)
			{
				mesh.property(face_type, *f_it) = 2;
			}
		}

		delete bp;
	}
	Weight_poly_list.clear();

}

void PlaneViewerWidget::SetWeightPolyline()
{
	if (Weight_bond_Point.size() >= 3)
	{
		BPolyline* bp = new BPolyline();
		bp->IntPolyline(Weight_bond_Point, true);
		Weight_poly_list.push_back(bp);
		Weight_bond_Point.clear();
	}

}

void PlaneViewerWidget::ChangeRefOK()
{
	double center_x, center_y;
	boundary_poly->FindCenter(center_x, center_y);
	double saler;
	if (All_Change_OK)
	{
		saler = sqrt(11.05 / computer_mesh_area(mesh));
	}
	else
	{
		saler = 1;
	}
	boundary_poly->Centralization(saler);
	MeshGeneration* mg = new MeshGeneration();
	mg->meshFromChangePolyline(boundary_poly, mesh);

	for (int wi = 0; wi < Weight_poly_list.size(); wi++)
	{
		Weight_poly_list[wi]->Centralization(center_x, center_y, saler);
	}

	SetMeshWeight();

	deformation_fix_id.clear();
	cut_Manager.clear();
	all_curve_reduct_cut.clear();

	planeview_draw_model = NOIC_PLANE;

	updateGL();
}
double PlaneViewerWidget::computer_mesh_area(Mesh& my_mesh)
{
	if (my_mesh.n_vertices() < 1)
	{
		return 0;
	}
	double my_area = 0;
	for (Mesh::FaceIter f_it = my_mesh.faces_begin(); f_it != my_mesh.faces_end(); f_it++)
	{
		std::vector<Mesh::Point> f_p;
		f_p.clear();
		for (Mesh::FaceVertexIter fv_it = my_mesh.fv_begin(*f_it); fv_it != my_mesh.fv_end(*f_it); fv_it++)
		{
			Mesh::Point p = my_mesh.point(*fv_it);
			f_p.push_back(p);
		}
		Mesh::Point v1 = f_p[1] - f_p[0];
		Mesh::Point v2 = f_p[2] - f_p[0];
		double area_f = (cross(v1, v2)).norm() / 2;
		my_area += area_f;
	}
	return my_area;
}

void PlaneViewerWidget::ComputeMaxH()
{
	double true_maxh = 0;

	for (Mesh::VertexIter v_it = curved_mesh.vertices_begin(); v_it != curved_mesh.vertices_end(); v_it++)
	{
		Mesh::Point p = curved_mesh.point(*v_it);
		if (max_curved_h < p[2])
		{
			max_curved_h = p[2];
		}
		if (true_maxh < p[2])
		{
			true_maxh = p[2];
		}
	}

	curved_bopundary_poly = new BPolyline();
	vector<Vector3f> boundary;
	boundary.clear();
	auto heit = curved_mesh.halfedges_begin();
	while (!curved_mesh.is_boundary(*heit))
		heit++;
	auto he_start = *heit;
	auto he_it = he_start;
	do
	{
		he_it = curved_mesh.next_halfedge_handle(he_it);
		Mesh::Point p = curved_mesh.point(curved_mesh.to_vertex_handle(he_it));
		boundary.push_back(Vector3f(p[0], p[1], 0));
	} while (he_it != he_start);

	curved_bopundary_poly->IntPolyline(boundary, true);

	if (max_curved_h > 3 * true_maxh && true_maxh > 0.1)
	{
		max_curved_h = 2.0*true_maxh;
	}

	std::cout << "max h is====" << max_curved_h << "====" << std::endl;
}


#pragma region TestDeformationFunction
int  PlaneViewerWidget::FindMeshClosetid(Eigen::Vector3f np)
{
	Mesh::Point p(np.x(), np.y(), 0);
	Mesh::Point v;
	double dis = 10000;
	int id = -1;
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		v = mesh.point(*v_it);
		if ((v - p).norm() < dis)
		{
			dis = (v - p).norm();
			id = v_it.handle().idx();
		}
	}
	if (dis < 0.1)
	{
		return id;
	}
	else
	{
		return -1;
	}
}
void PlaneViewerWidget::My_Deformation_prepare(int handle_id)
{
	std::vector<int> my_face_in_out;
	my_face_in_out.clear();
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		if (mesh.property(is_singular, *f_it))
		{
			my_face_in_out.push_back(-1);
		}
		else
		{
			my_face_in_out.push_back(1);
		}
	}

	std::vector<int> my_handle_id;
	my_handle_id.clear();
	my_handle_id.push_back(handle_id);

	parafun_m = Parafun(mesh, 4);
	parafun_m.init_Deformation(deformation_fix_id, my_handle_id, my_face_in_out);

}
void PlaneViewerWidget::My_Deformotion(Mesh::Point P)
{
	std::vector<double> h_x, h_y;
	h_x.clear();
	h_y.clear();
	h_x.push_back(P[0]);
	h_y.push_back(P[1]);

	parafun_m.init_handle(h_x, h_y);
	parafun_m.run_cm();
	parafun_m.ChangeVertPoint(mesh);

	updateGL();
}
void  PlaneViewerWidget::My_Deformotion_finally()
{
	int pn = boundary_poly->polyline_x.size() - 1;
	for (int i = 0; i < pn; i++)
	{
		Mesh::Point p = mesh.point(mesh.vertex_handle(i));
		boundary_poly->polyline_x[i] = p[0];
		boundary_poly->polyline_y[i] = p[1];
	}
	boundary_poly->polyline_x[pn] = boundary_poly->polyline_x[0];
	boundary_poly->polyline_y[pn] = boundary_poly->polyline_y[0];
}
#pragma endregion
