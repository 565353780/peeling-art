#ifndef MESHPROCESSING_MAIN_VIEWWE_WIDGET_H
#define MESHPROCESSING_MAIN_VIEWWE_WIDGET_H

#include <QtGui>
#include <QString>
#include <QMessageBox>
#include <QFileDialog>
//main widget
#include "InteractiveViewerWidget.h"
#include "PlaneViewer/PlaneViewerWidget.h"
#include "PlaneViewer/mainchildwidget.h"

extern bool _DEBUG_MODE;

class MainViewerWidget : public QDialog
{
	Q_OBJECT
public:
	MainViewerWidget(QWidget* _parent=0);
	~MainViewerWidget();

	void setMouseMode(int dm)
	{
		MeshViewer->setMouseMode(dm);
	}
	void openMesh_fromMain(char* filename)
	{
		QString str(filename);
		open_mesh_gui(filename);
	}

public slots:
	void open_mesh_query()
	{
		QString fileName = QFileDialog::getOpenFileName(this,
			tr("Open mesh file"),
			tr("../models/"),
			tr("OBJ Files (*.obj);;"
			"OFF Files (*.off);;"
			"PLY Files (*.ply);;"
			"STL Files (*.stl);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
		{
			open_mesh_gui(fileName);
		}
	}
	void save_mesh_query() 
	{
		QString fileName = QFileDialog::getSaveFileName(this,
			tr("Save mesh file"),
			tr("../models/untitled.obj"),
			tr("OFF Files (*.obj);;"
			"OBJ Files (*.off);;"
			"PLY Files (*.ply);;"
			"STL Files (*.stl);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
		{
			save_mesh_gui(fileName);
		}
	}
	void saveOpenGLScreen()
	{
		QString fileName = QFileDialog::getSaveFileName(this,
			("Save screen as image file"),
			("../Results/untitled.png"),
			("PNG Files (*.png);;BMP Files (*.bmp);;JPG Files (*.jpg);;"
			"All Files (*)"));
		if (!fileName.isEmpty())
		{
			save_screen_gui(fileName);
		}
	}
	void save_opengl_screen(QString str)
	{
		MeshViewer->saveScreen(str.toLocal8Bit());
	}

	void LoadMeshFromInner(bool OK,QString fname)
	{
		LoadMeshSuccess = OK;
		if(LoadMeshSuccess)
		{
			SetMeshForALL();
		}
		emit( haveLoadMesh(fname) );
	};

signals:
	void haveLoadMesh(QString filePath);
	void setMouseMode_signal_main(int);
	void setDrawMode_signal_main(int);

protected:
	virtual void initViewerWindow();
	virtual void createViewerDialog();
	virtual void save_mesh_gui(QString fname);
	virtual void open_mesh_gui(QString fname);
	virtual void save_screen_gui(QString fname);

protected:
	bool LoadMeshSuccess;
	void SetMeshForALL()
	{
	}

private:
	
	InteractiveViewerWidget* MeshViewer;
	MainChildWidget* child_viewer;
	

public:
	void MyARAP()
	{
		if (MeshViewer->plane_Eye_id == -2)
		{
			MeshViewer->SetEye();
		}


		MeshViewer->MyARAP3d();
        if(_DEBUG_MODE)
        {
            if (!child_viewer->isVisible() && !MeshViewer->is_finally_change_shape)
            {
                child_viewer->show();
                child_viewer->PlaneViewer->setBoundary(MeshViewer->ref_boundary_);
                child_viewer->PlaneViewer->setCutManager(MeshViewer->cut_mesh_mg);
    //			child_viewer->PlaneViewer->mesh = MeshViewer->getRefMesh();
                child_viewer->PlaneViewer->mesh = MeshViewer->getCurvedMesh();
            }
            else if (!child_viewer->isVisible() && MeshViewer->is_finally_change_shape)
            {
                child_viewer->show();
                child_viewer->PlaneViewer->setBoundary(MeshViewer->ref_boundary_);
                child_viewer->PlaneViewer->cut_Manager.clear();
    //			child_viewer->PlaneViewer->mesh = MeshViewer->getRefMesh();
                child_viewer->PlaneViewer->mesh = MeshViewer->getCurvedMesh();
            }
        }
	}
	void LoadFace()
	{
        MeshViewer->png_idx = png_idx;

        if(!LoadMeshSuccess)
        {
            LoadMeshSuccess = true;
            MeshViewer->setDrawMode(InteractiveViewerWidget::FLAT_POINTS);
            MeshViewer->setMouseMode(InteractiveViewerWidget::TRANS);
            MeshViewer->LoadFace(false);
        }
        else
        {
            MeshViewer->clearAllMesh();

            MeshViewer->have_first_load = false;
            MeshViewer->is_show_papamaterization = false;
            MeshViewer->NO_Papamaterization = true;
            MeshViewer->papamaterization_times = 0;
            MeshViewer->is_finally_change_shape = false;
            MeshViewer->is_deformation_prepare = false;
            MeshViewer->plane_Eye_id = -2;
            MeshViewer->draw_plane_mode = false;
            MeshViewer->graph_init = false;

            MeshViewer->setDrawMode(InteractiveViewerWidget::FLAT_POINTS);
            MeshViewer->setMouseMode(InteractiveViewerWidget::TRANS);
            MeshViewer->LoadFace(true);
        }

	}

	void RefChangeOK()
	{
		child_viewer->PlaneViewer->ChangeRefOK();
		bool a = child_viewer->PlaneViewer->All_Change_OK;
		MeshViewer->setChangeMesh(child_viewer->PlaneViewer->mesh, a);
		child_viewer->hide();
	}

public:
    int png_idx = 0;

};


#endif
