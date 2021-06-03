#include "SurfaceMeshProcessingWindow.h"
#include "MeshViewer/MainViewerWidget.h"

SurfaceMeshProcessingWindow::SurfaceMeshProcessingWindow(QWidget *parent)
	: QMainWindow(parent)
{
	setWindowTitle("Peeling");

	viewer = new MainViewerWidget(this);
	setCentralWidget(viewer);
	createActions();
	createToolBars();
	createStatusBar();
}

SurfaceMeshProcessingWindow::~SurfaceMeshProcessingWindow()
{

}

void SurfaceMeshProcessingWindow::createActions()
{
	//按钮_打开文件
	openAction = new QAction(tr("&Open"), this);
	openAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/Open.png"));
	openAction->setShortcut(QKeySequence::Open);
	openAction->setStatusTip(tr("Open a mesh file"));
	connect(openAction, SIGNAL(triggered()), viewer, SLOT(open_mesh_query()));

	//按钮_保存文件
	saveAction = new QAction(tr("&Save"), this);
	saveAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/Save.png"));
	saveAction->setShortcut(QKeySequence::Save);
	saveAction->setStatusTip(tr("Save the mesh to file"));
	connect(saveAction, SIGNAL(triggered()), viewer, SLOT(save_mesh_query()));

	//按钮_退出程序
	exitAction = new QAction(tr("&Exit"), this);
	exitAction->setShortcut(tr("Ctrl+Q"));
	exitAction->setStatusTip(tr("Exit the application"));
	connect(exitAction, SIGNAL(triggered()), this, SLOT(close()));

	//按钮_保存截图
	saveScreenAction = new QAction("Save Screen",this);
	saveScreenAction->setIcon(QIcon(":/SurfaceMeshProcessing/Images/saveScreen.png"));
	saveScreenAction->setStatusTip(tr("Save Screen"));
	connect(saveScreenAction, SIGNAL(triggered()), viewer, SLOT(saveOpenGLScreen()));

	//按钮_
	InlineFixAction = new QAction("Set Handle", this);
	InlineFixAction->setStatusTip(tr("Set Deformation Handle"));
	connect(InlineFixAction, SIGNAL(triggered()), this, SLOT(SetHandle()));

	moveAction = new QAction("Move Handle", this);
	moveAction->setStatusTip(tr("Move Deformation handle"));
	connect(moveAction, SIGNAL(triggered()), this, SLOT(moveHandle()));

	MyARAPAction = new QAction("Begin", this);
	MyARAPAction->setStatusTip(tr("Begin to get a Citrus"));
	connect(MyARAPAction, SIGNAL(triggered()), this, SLOT(MyARAP3D()));

	MyLoadAction = new QAction("Load ", this);
	MyLoadAction->setStatusTip(tr("Load Curve"));
	connect(MyLoadAction, SIGNAL(triggered()), this, SLOT(LoadFace()));

	MyRefOKAction = new QAction("Change Ok", this);
	MyRefOKAction->setStatusTip(tr("Ref Change Ok"));
	connect(MyRefOKAction, SIGNAL(triggered()), this, SLOT(RefChangeOk()));

	SetEyeAction = new QAction("Set Root", this);
	SetEyeAction->setStatusTip(tr("set eye position int Citrus"));
	connect(SetEyeAction, SIGNAL(triggered()), this, SLOT(SetEye()));

}

void SurfaceMeshProcessingWindow::createMenus()
{
	fileMenu = menuBar()->addMenu(tr("&File"));
	fileMenu->addAction(openAction);
	fileMenu->addAction(saveAction);
	fileMenu->addSeparator();
	fileMenu->addAction(saveScreenAction);
	fileMenu->addSeparator();
	fileMenu->addAction(exitAction);

	helpMenu = menuBar()->addMenu(tr("&Help"));
}

void SurfaceMeshProcessingWindow::createToolBars()
{
	Peeling = addToolBar(tr("Peeling Art Design"));
	
	Peeling->addAction(MyLoadAction);
	Peeling->addAction(SetEyeAction);
	Peeling->addAction(MyARAPAction);	

	Peeling->addAction(InlineFixAction);
	Peeling->addAction(moveAction);
	
	Peeling->addAction(MyRefOKAction);
}

void SurfaceMeshProcessingWindow::createStatusBar()
{
	statusLabel = new QLabel(tr("no mesh"));
	statusLabel->setAlignment(Qt::AlignHCenter);

	connect(viewer,SIGNAL( haveLoadMesh(QString) ),statusLabel,SLOT( setText(QString) ) );

	statusBar()->addWidget(statusLabel);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void SurfaceMeshProcessingWindow::moveHandle()
{
	moveAction->setChecked(true);
	viewer->setMouseMode(InteractiveViewerWidget::MOVE);
}
void SurfaceMeshProcessingWindow::SetHandle()
{
	InlineFixAction->setChecked(true);
	viewer->setMouseMode(InteractiveViewerWidget::HANDLE);
}
void SurfaceMeshProcessingWindow::MyARAP3D()
{
	MyARAPAction->setChecked(true);
	viewer->MyARAP();
}
void SurfaceMeshProcessingWindow::LoadFace()
{
	MyARAPAction->setChecked(true);
	viewer->LoadFace();
}

void SurfaceMeshProcessingWindow::RefChangeOk()
{
	MyRefOKAction->setChecked(true);
	viewer->RefChangeOK();
}
void SurfaceMeshProcessingWindow::SetEye()
{
	SetEyeAction->setChecked(true);
	viewer->setMouseMode(InteractiveViewerWidget::EYE);
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

void SurfaceMeshProcessingWindow::open_mesh_from_main(char* filename)
{
	viewer->openMesh_fromMain(filename);
}

