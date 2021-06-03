#include <iostream>
#include <algorithm>
#include <cmath>

#include <qapplication.h>

#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "MeshViewerWidget.h"
#include "Common/CommonDefinitions.h"

using namespace Qt;

MeshViewerWidget::MeshViewerWidget(QWidget* parent)
 : QGLViewerWidget(parent)
{
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();

	mesh.request_face_normals();
	mesh.request_vertex_normals();

	mesh_vector.clear(); mesh_vector_index = -1;

	draw_BBox_OK = false;
	draw_mesh_boundary_ok = false;
	first_init = true;
}

MeshViewerWidget::MeshViewerWidget(QGLFormat& _fmt, QWidget* _parent)
: QGLViewerWidget(_fmt, _parent)
{
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();

	mesh.request_face_normals();
	mesh.request_vertex_normals();

	mesh_vector.clear(); mesh_vector_index = -1;

	draw_BBox_OK = false;
	draw_mesh_boundary_ok = false;
	first_init = true;
}

MeshViewerWidget::~MeshViewerWidget()
{
}

void MeshViewerWidget::updateMeshCenter()
{
	typedef Mesh::Point Point;
	Mesh::VertexIter vIt = mesh.vertices_begin();
	Mesh::VertexIter vEnd = mesh.vertices_end();
	bbMin = bbMax = OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh.point(vIt));

	size_t count = 0;
	for(; vIt != vEnd; ++vIt, ++count)
	{
		bbMin.minimize( OpenMesh::vector_cast<OpenMesh::Vec3d>( mesh.point(vIt) ) );
		bbMax.maximize( OpenMesh::vector_cast<OpenMesh::Vec3d>( mesh.point(vIt) ) );
	}

	Mesh::EdgeIter e_it = mesh.edges_begin();
	Mesh::EdgeIter e_end = mesh.edges_end();
	double aveLen = 0.0; double maxLen = 0.0; double minLen = mesh.calc_edge_length(e_it);
	double e_len = 0.0;
	for(; e_it != e_end; ++e_it)
	{
		double e_len = mesh.calc_edge_length(e_it);
		if( e_len > maxLen )
		{
			maxLen = e_len;
		}
		else if(e_len < minLen )
		{
			minLen = e_len;
		}
		aveLen += e_len;
	}

	if( first_init )
	{
		set_scene_pos( (bbMin+bbMax)*0.5, (bbMin-bbMax).norm()*0.5 );
		first_init = false;
	}
	else
	{
		set_scene_pos((bbMin + bbMax)*0.5, (bbMin - bbMax).norm()*0.5);
	}
	
	
	printf("BoundingBox:\nX : [ %f , %f ]\n",bbMin[0],bbMax[0]);
	printf("Y : [ %f , %f ]\n",bbMin[1],bbMax[1]);
	printf("Z : [ %f , %f ]\n",bbMin[2],bbMax[2]);
	printf("Diag length of BBox : %f\n", (bbMax-bbMin).norm() );
	printf("Edge Length : Max : %f; Min : %f; AVG : %f\n",maxLen,minLen, aveLen/mesh.n_edges() );
}

void MeshViewerWidget::updateMeshNormals()
{
	mesh.update_face_normals();
	mesh.update_vertex_normals();
}

bool MeshViewerWidget::openMesh(const char* filename)
{
	clearAllMesh();
	bool read_OK = OpenMesh::IO::read_mesh( mesh, filename );

	printf("%s\n", filename);
	if ( read_OK )
	{
		initMesh();
		mesh_vector.push_back(mesh); mesh_vector_index = 0;
		return true;
	}
	return false;
}

void MeshViewerWidget::initMesh()
{
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();

	mesh.request_face_normals();
	mesh.request_vertex_normals();
	printBasicMeshInfo();
	updateMesh();
}

void MeshViewerWidget::printBasicMeshInfo()
{
	if (mesh.n_vertices() == 0)
		printf("No Mesh\n");

	checkMeshMode();

	QString meshType;
	switch(meshMode())
	{
	case TRIANGLE:
		printf("Triangle Mesh.\n");
		break;
	case QUAD:
		printf("Quadrilateral Mesh.\n");
		break;
	default:
		printf("General Mesh.\n");
		break;
	}

	printf("Information of the input mesh:\nVertex : %d;\nFace : %d;\nEdge : %d, HalfEdge : %d\n",
		mesh.n_vertices(),mesh.n_faces(),mesh.n_edges(),mesh.n_halfedges());;
}


bool MeshViewerWidget::saveMesh(const char* filename)
{
	return OpenMesh::IO::write_mesh(mesh, filename);
}
bool MeshViewerWidget::saveScreen(const char* filePath)
{
	QImage image = grabFrameBuffer();
	image.save(filePath);
	return true;
}

void MeshViewerWidget::checkMeshMode()
{
	Mesh::FaceIter fIt = mesh.faces_begin();
	Mesh::FaceIter fEnd = mesh.faces_end();
	Mesh::FaceEdgeIter fe_it;
	int count = 1;
	int meshType[3] = {0};
	for(fIt; fIt != fEnd; ++fIt)
	{
		fe_it = mesh.fe_iter(fIt);
		while(--fe_it)
		{
			++count;
		}
		if(count == 4)
		{
			meshType[1]++;
		}
		else if(count == 3)
		{
			meshType[0]++;
		}
		else
		{
			meshType[2]++;
		}
		count = 1;
	}
	int faceNum = mesh.n_faces();
	if(meshType[0] == faceNum)//triangle
	{
		setMeshMode(TRIANGLE);
	}
	else if(meshType[1] == faceNum)//no 
	{
		setMeshMode(QUAD);
	}
	else
	{
		setMeshMode(N_MESH_MODES);
	}
}

void MeshViewerWidget::draw_scene(int drawmode)
{
	QFont Text_Font("Courier", 12);
	glViewport ( 0,0, width(),height());
	glMatrixMode( GL_PROJECTION );
	glLoadMatrixd( &ProjectionMatrix[0] );
	glMatrixMode( GL_MODELVIEW );
	glLoadMatrixd( &ModelViewMatrix[0] );

	if(draw_BBox_OK)
	{
		OpenMesh::Vec3d temp0 = bbMin;
		OpenMesh::Vec3d temp1;
		glLineWidth(2.0);
		glColor3f(1.0, 1.0, 0.0);
		glBegin(GL_LINES);
		temp1 = bbMin; temp1[0] = bbMax[0];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );
		temp1 = bbMin; temp1[1] = bbMax[1];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );
		temp1 = bbMin; temp1[2] = bbMax[2];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[0] = bbMax[0];
		temp1 = bbMax; temp1[1] = bbMin[1];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[0] = bbMax[0];
		temp1 = bbMax; temp1[2] = bbMin[2];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[1] = bbMax[1];
		temp1 = bbMax; temp1[2] = bbMin[2];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[1] = bbMax[1];
		temp1 = bbMax; temp1[0] = bbMin[0];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[2] = bbMax[2];
		temp1 = bbMax; temp1[1] = bbMin[1];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[2] = bbMax[2];
		temp1 = bbMax; temp1[0] = bbMin[0];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMax;
		temp1 = bbMax; temp1[0] = bbMin[0];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );
		temp1 = bbMax; temp1[1] = bbMin[1];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );
		temp1 = bbMax; temp1[2] = bbMin[2];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );
		glEnd();
	}

	if(draw_mesh_boundary_ok)
	{
		glLineWidth(2.0);
		glColor3f(1.0, 0.5, 0.0);
		glBegin(GL_LINES);
		for(Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
		{
			if( mesh.is_boundary(e_it ) )
			{
				Mesh::HalfedgeHandle heh0 = mesh.halfedge_handle(e_it, 0);
				glVertex3dv( mesh.point(mesh.to_vertex_handle(heh0)).data() );
				glVertex3dv( mesh.point(mesh.from_vertex_handle(heh0)).data() );
			}
		}
		glEnd();

		/*FILE* f_bde = fopen("bde.de", "w");
		for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
		{
			if (mesh.is_boundary(v_it))
			{
				OpenMesh::Vec3d& np = mesh.point(v_it);
				fprintf(f_bde, "%d %20.19f %20.19f 0.0\n", v_it.handle().idx(), np[0], np[1]);
			}
		}
		fclose(f_bde);*/
	}

	draw_scene_mesh(drawmode);
}

void MeshViewerWidget::draw_scene_mesh(int drawmode)
{
	if(mesh.n_vertices() == 0) { return; }

	if (drawmode == -100) {
		return;
	}

	switch (drawmode)
	{
	case WIRE_FRAME:
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		draw_mesh_wireframe();
		//draw_meshpointset();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		break;
	case HIDDEN_LINES:
		glDisable(GL_LIGHTING);
		draw_mesh_hidden_lines();
		break;
	case SOLID_FLAT:
		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		draw_mesh_solidflat();
		//draw_meshpointset();
		glDisable(GL_LIGHTING);
		break;
	case FLAT_POINTS:
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.5f,2.0f);
		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		draw_mesh_solidflat();
		glDisable(GL_POLYGON_OFFSET_FILL);
		//draw_meshpointset();
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		draw_mesh_wireframe();
		//draw_meshpointset();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		break;
	case SOLID_SMOOTH:
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);
		draw_mesh_solidsmooth();
		//draw_meshpointset();
		glDisable(GL_LIGHTING);
		break;
	case POINT_SET:
		glDisable(GL_LIGHTING);
		draw_mesh_pointset();
		break;
	case CURVATURE:
		
		break;
	default:
		break;
	}
}


void MeshViewerWidget::draw_mesh_wireframe()
{
	glLineWidth(1);
	//glColor3f(0.753, 0.753, 0.753);
	//glColor3f(0.0, 0.0, 0.0);
	glColor3f(0.0, 0.0, 0.25);
	//if(meshMode() != TRIANGLE && meshMode() != QUAD)
	{
		Mesh::ConstFaceIter fIt(mesh.faces_begin()),
			fEnd(mesh.faces_end());
		Mesh::ConstFaceVertexIter fvIt;
		for (; fIt != fEnd; ++fIt)
		{
			//GL::glNormal(dualMesh.normal(f_it));
			fvIt = mesh.cfv_iter(fIt); 
			glBegin(GL_POLYGON);
			for( fvIt; fvIt; ++fvIt )
			{
				glVertex3dv( mesh.point(fvIt).data() );
			}
			glEnd();
		}
	}

	/*OpenMesh::Vec3d pos = mesh.point(mesh.vertex_handle(0));
	GLdouble  winX, winY, winZ;
	GLint     viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	gluProject(pos[0],pos[1],pos[2],&ModelViewMatrix[0][0],&ProjectionMatrix[0][0],viewport,&winX,&winY,&winZ);
	int x = (long)winX;
	int y = viewport[3]-(long)winY;
	QString str = "fxum";
	render_text(x,y,str);*/
	/*else
	{
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_DOUBLE, 0, mesh.points());
		switch(meshMode())
		{
		case TRIANGLE:
			glDrawElements(GL_TRIANGLES,Indices.size(),GL_UNSIGNED_INT,&Indices[0]);
			break;
		case QUAD:
			glDrawElements(GL_QUADS,Indices.size(),GL_UNSIGNED_INT,&Indices[0]);
			break;
		default:
			break;
		}

		glDisableClientState(GL_VERTEX_ARRAY);
	}*/

}


void MeshViewerWidget::draw_mesh_hidden_lines() const
{
	Mesh::ConstFaceIter f_it(mesh.faces_begin());
	Mesh::ConstFaceIter f_end(mesh.faces_end());
	Mesh::ConstFaceVertexIter fv_it;
	glLineWidth(2.0);
	glColor3f(0.0, 1.0, 1.0);

	glDrawBuffer(GL_NONE);
	glDepthRange(0.01, 1.0);
	switch(meshMode())
	{
	case TRIANGLE:
		glBegin(GL_TRIANGLES);
		for (; f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}

		}
		glEnd();
		break;
	case QUAD:
		glBegin(GL_QUADS);
		for (; f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
		}
		glEnd();
		break;
	default:
		for (; f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			glBegin(GL_POLYGON);
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
			glEnd();
		}
		break;
	}

	glDrawBuffer(GL_BACK);
	glDepthRange(0.0, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glDepthFunc(GL_LEQUAL);

	switch(meshMode())
	{
	case TRIANGLE:
		glBegin(GL_TRIANGLES);
		for (f_it = mesh.faces_begin(); f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
		}
		glEnd();
		break;
	case QUAD:
		glBegin(GL_QUADS);
		for (f_it = mesh.faces_begin(); f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
		}
		glEnd();
		break;
	default:
		for (f_it = mesh.faces_begin(); f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			glBegin(GL_POLYGON);
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
			glEnd();
		}
		break;
	}

	glDepthFunc(GL_LESS);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void MeshViewerWidget::draw_mesh_solidflat() const
{
	Mesh::ConstFaceIter fIt(mesh.faces_begin()),
		fEnd(mesh.faces_end());
	Mesh::ConstFaceVertexIter fvIt;

	GLfloat mat_a[] = { 0.7f, 0.7f, 0.7f, 1.0f };
	GLfloat mat_d[] = { 0.88f, 0.84f, 0.76f, 1.0f };
	GLfloat mat_s[] = { 0.4f, 0.4f, 0.4f, 1.0f };
	GLfloat shine[] = { 120.0f };

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_a);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_d);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_s);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);

	switch(meshMode())
	{
	case TRIANGLE:
		glBegin(GL_TRIANGLES);
		for (fIt; fIt != fEnd; ++fIt) 
		{
			glNormal3dv(mesh.normal(fIt).data());
			fvIt = mesh.cfv_iter(fIt.handle());
			for(fvIt; fvIt; ++fvIt)
			{
				glVertex3dv(mesh.point(fvIt).data());
			}
		}
		glEnd();
		break;
	case QUAD:
		glBegin(GL_QUADS);
		for (;fIt != fEnd; ++fIt) 
		{
			glNormal3dv(&mesh.normal(fIt)[0]);
			fvIt = mesh.cfv_iter(fIt.handle());
			for(fvIt; fvIt; ++fvIt)
			{
				glVertex3dv(&mesh.point(fvIt)[0]);
			}
		}
		glEnd();
		break;
	default:
		for (; fIt != fEnd; ++fIt) 
		{
			glBegin(GL_POLYGON);
			glNormal3dv(&mesh.normal(fIt)[0]);
			fvIt = mesh.cfv_iter(fIt.handle());
			for(fvIt; fvIt; ++fvIt)
			{
				glVertex3dv(&mesh.point(fvIt)[0]);
			}
			glEnd();
		}
		break;
	}
}

void MeshViewerWidget::draw_mesh_solidsmooth() const 
{
	bool drawOK = false;
	glLoadName(mesh.n_vertices());

	Mesh::ConstFaceIter fIt(mesh.faces_begin()),
		fEnd(mesh.faces_end());
	Mesh::ConstFaceVertexIter fvIt;

	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, mesh.points());
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_DOUBLE, 0, mesh.vertex_normals());

	for (; fIt != fEnd; ++fIt)
	{
		glBegin(GL_POLYGON);
		fvIt = mesh.cfv_iter(fIt.handle());
		for(fvIt; fvIt ; --fvIt)
		{
			glArrayElement(fvIt.handle().idx());
		}
		glEnd();
	}

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
}

void MeshViewerWidget::draw_mesh_pointset() const 
{
	glPolygonMode(GL_FRONT_AND_BACK, GL_POINTS);

	glColor3f(0.0, 1.0, 1.0);
	glPointSize(10);
	Mesh::VertexIter v_it = mesh.vertices_begin();
	glBegin(GL_POINTS);
	for(v_it;v_it != mesh.vertices_end();++v_it)
	{
		glVertex3dv(mesh.point(v_it).data());
	}
	glEnd();

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

}
