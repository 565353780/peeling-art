#include "MainViewerWidget.h"
#include <QHBoxLayout>

MainViewerWidget::MainViewerWidget(QWidget* _parent/* =0 */)
{
	initViewerWindow();
	LoadMeshSuccess = false;

	child_viewer = new MainChildWidget();
}
MainViewerWidget::~MainViewerWidget()
{
};

void MainViewerWidget::initViewerWindow()
{
	createViewerDialog();

	QHBoxLayout* main_layout = new QHBoxLayout();
	main_layout->addWidget(MeshViewer, 1);
	this->setLayout(main_layout);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 

	connect(MeshViewer,SIGNAL(setMouseMode_signal(int)),SIGNAL(setMouseMode_signal_main(int)));
	connect(MeshViewer,SIGNAL(setDrawMode_signal(int)),SIGNAL(setDrawMode_signal_main(int)));
}

void MainViewerWidget::createViewerDialog()
{
	QGLFormat glFormat;
	glFormat.setSampleBuffers(true);
	glFormat.setSamples(16);

	MeshViewer = new InteractiveViewerWidget(glFormat, NULL);
	MeshViewer->setAcceptDrops(true);
	connect(MeshViewer,SIGNAL(loadMeshOK(bool,QString)), this, SLOT(LoadMeshFromInner(bool,QString)) );
}

void MainViewerWidget::open_mesh_gui(QString fname)
{
	if (fname.isEmpty() || !MeshViewer->openMesh(fname.toLocal8Bit())) 
	{
		QString msg = "Cannot read mesh from file:\n '";
		msg += fname;
		msg += "'";
		QMessageBox::critical(NULL, windowTitle(), msg);
	}
	else
	{
		LoadMeshSuccess = true;
		MeshViewer->setDrawMode(InteractiveViewerWidget::FLAT_POINTS);
		MeshViewer->setMouseMode(InteractiveViewerWidget::TRANS);
		if(LoadMeshSuccess)
		{
			SetMeshForALL();
		}
		emit(haveLoadMesh(fname));
	}
}

void MainViewerWidget::save_mesh_gui(QString fname)
{
	if (fname.isEmpty() || !MeshViewer->saveMesh(fname.toLocal8Bit()))
	{
		QString msg = "Cannot read mesh from file:\n '";
		msg += fname;
		msg += "'";
		QMessageBox::critical(NULL, windowTitle(), msg);
	}
}

void MainViewerWidget::save_screen_gui(QString fname)
{
	if (fname.isEmpty() || !MeshViewer->saveScreen(fname.toLocal8Bit()))
	{
		QString msg = "Cannot save image to file:\n '";
		msg += fname;
		msg += "'";
		QMessageBox::critical(NULL, windowTitle(), msg);
	}
}
