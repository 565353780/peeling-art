#ifndef MESHPROCESSING_MESHVIEWERWIDGET_H
#define MESHPROCESSING_MESHVIEWERWIDGET_H

#include <QString>
#include <QMessageBox>
#include <QFileDialog>

#include "QGLViewerWidget.h"
#include "MeshDefinition.h"

class MeshViewerWidget : public QGLViewerWidget 
{
	Q_OBJECT
public:
	MeshViewerWidget(QWidget* parent = 0);
	MeshViewerWidget(QGLFormat& _fmt, QWidget* _parent);
	~MeshViewerWidget();
public:
	bool openMesh(const char* filename);
	void initMesh();
	bool saveMesh(const char* filename);
	bool saveScreen(const char* filePath);
	void updateMesh()
	{
		updateMeshNormals();
		updateMeshCenter();
	};
	virtual void clearAllMesh()
    {
		mesh_vector.clear(); mesh_vector_index = -1;
        mesh.clear();
		updateGL();
	};
	void printBasicMeshInfo();

signals:
	void loadMeshOK(bool,QString);

protected:
	void updateMeshCenter(); // used by update_mesh().
	void updateMeshNormals(); // used by update_mesh().

protected:
	virtual void draw_scene(int drawmode);
	void draw_scene_mesh(int drawmode);
	
private:
	void draw_mesh_wireframe();
	void draw_mesh_hidden_lines() const;
	void draw_mesh_solidflat() const;
	void draw_mesh_solidsmooth() const;
	void draw_mesh_pointset() const;

protected:
	bool first_init;
	OpenMesh::Vec3d bbMin;
	OpenMesh::Vec3d bbMax;

	Mesh mesh;
	Mesh ref_mesh;
	Mesh cut_mesh; //储存没被覆盖的部分，在分割之后储存覆盖了的部分。

	Mesh curved_mesh;

	bool draw_BBox_OK;
	bool draw_mesh_boundary_ok;
	std::vector<Mesh> mesh_vector;
	int mesh_vector_index;

public:
	// mesh modes.
	enum { TRIANGLE = 0, QUAD, N_MESH_MODES };
	void setMeshMode(int mm) { mesh_mode_ = mm;}
	int meshMode() const { return mesh_mode_; }
	void checkMeshMode();
private:
	int mesh_mode_;
};

#endif // MESHPROCESSING_MESHVIEWERWIDGET_H
