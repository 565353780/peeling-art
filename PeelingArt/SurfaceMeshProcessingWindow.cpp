#include "SurfaceMeshProcessingWindow.h"
#include "MeshViewer/MainViewerWidget.h"

SurfaceMeshProcessingWindow::SurfaceMeshProcessingWindow(QWidget *parent)
	: QMainWindow(parent)
{
    setWindowTitle("Peeling Art v1.0");

	viewer = new MainViewerWidget(this);
	setCentralWidget(viewer);
	createActions();
    createMenus();
	createToolBars();
	createStatusBar();

    if(!_DEBUG_MODE)
    {
        createButtons();
        ButtonVisiable(false);
    }
}

SurfaceMeshProcessingWindow::~SurfaceMeshProcessingWindow()
{

}

void SurfaceMeshProcessingWindow::createActions()
{
	openAction = new QAction(tr("&Open"), this);
    openAction->setIcon(QIcon("./Images/Open.png"));
	openAction->setShortcut(QKeySequence::Open);
	openAction->setStatusTip(tr("Open a mesh file"));
	connect(openAction, SIGNAL(triggered()), viewer, SLOT(open_mesh_query()));

	saveAction = new QAction(tr("&Save"), this);
    saveAction->setIcon(QIcon("./Images/Save.png"));
	saveAction->setShortcut(QKeySequence::Save);
	saveAction->setStatusTip(tr("Save the mesh to file"));
	connect(saveAction, SIGNAL(triggered()), viewer, SLOT(save_mesh_query()));

	exitAction = new QAction(tr("E&xit"), this);
	exitAction->setShortcut(tr("Ctrl+Q"));
	exitAction->setStatusTip(tr("Exit the application"));
	connect(exitAction, SIGNAL(triggered()), this, SLOT(close()));

	saveScreenAction = new QAction("Save Screen",this);
    saveScreenAction->setIcon(QIcon("./Images/saveScreen.png"));
	saveScreenAction->setStatusTip(tr("Save Screen"));
	connect(saveScreenAction, SIGNAL(triggered()), viewer, SLOT(saveOpenGLScreen()));

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
    MyLoadAction->setIcon(QIcon("./Images/Open.png"));
    MyLoadAction->setStatusTip(tr("Load Curve"));
    connect(MyLoadAction, SIGNAL(triggered()), this, SLOT(LoadFace()));

	MyRefOKAction = new QAction("Change Ok", this);
	MyRefOKAction->setStatusTip(tr("Ref Change Ok"));
	connect(MyRefOKAction, SIGNAL(triggered()), this, SLOT(RefChangeOk()));

	SetEyeAction = new QAction("Set Root", this);
    SetEyeAction->setStatusTip(tr("Set eye position int Citrus"));
	connect(SetEyeAction, SIGNAL(triggered()), this, SLOT(SetEye()));

    /*
    PngHeAction = new QAction("He", this);
    PngHeAction->setIcon(QIcon("./Pngs/he.png"));
    connect(PngHeAction, SIGNAL(triggered()), this, SLOT(ChoosePng1()));

    PngMaAction = new QAction("Ma", this);
    PngMaAction->setIcon(QIcon("./Pngs/ma.png"));
    connect(PngMaAction, SIGNAL(triggered()), this, SLOT(ChoosePng2()));

    PngNiaoAction = new QAction("Niao", this);
    PngNiaoAction->setIcon(QIcon("./Pngs/niao.png"));
    connect(PngNiaoAction, SIGNAL(triggered()), this, SLOT(ChoosePng3()));

    PngShuAction = new QAction("Shu", this);
    PngShuAction->setIcon(QIcon("./Pngs/shu.png"));
    connect(PngShuAction, SIGNAL(triggered()), this, SLOT(ChoosePng4()));

    PngTuAction = new QAction("Tu", this);
    PngTuAction->setIcon(QIcon("./Pngs/tu.png"));
    connect(PngTuAction, SIGNAL(triggered()), this, SLOT(ChoosePng5()));

    PngWaAction = new QAction("Wa", this);
    PngWaAction->setIcon(QIcon("./Pngs/wa.png"));
    connect(PngWaAction, SIGNAL(triggered()), this, SLOT(ChoosePng6()));

    PngXiaAction = new QAction("Xia", this);
    PngXiaAction->setIcon(QIcon("./Pngs/xia.png"));
    connect(PngXiaAction, SIGNAL(triggered()), this, SLOT(ChoosePng7()));

    PngXieAction = new QAction("Xie", this);
    PngXieAction->setIcon(QIcon("./Pngs/xie.png"));
    connect(PngXieAction, SIGNAL(triggered()), this, SLOT(ChoosePng8()));

    PngYingAction = new QAction("Ying", this);
    PngYingAction->setIcon(QIcon("./Pngs/ying.png"));
    connect(PngYingAction, SIGNAL(triggered()), this, SLOT(ChoosePng9()));
    */
}

void SurfaceMeshProcessingWindow::createMenus()
{
    if(_DEBUG_MODE)
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
}

void SurfaceMeshProcessingWindow::createToolBars()
{
    Peeling = addToolBar(tr("Peeling Art Design"));

    Peeling->addAction(MyLoadAction);

    if(_DEBUG_MODE)
    {
        Peeling->addAction(SetEyeAction);
        Peeling->addAction(MyARAPAction);

        Peeling->addAction(InlineFixAction);
        Peeling->addAction(moveAction);

        Peeling->addAction(MyRefOKAction);
    }
}

void SurfaceMeshProcessingWindow::createStatusBar()
{
    if(_DEBUG_MODE)
    {
        statusLabel = new QLabel(tr("no mesh"));
        statusLabel->setAlignment(Qt::AlignHCenter);

        connect(viewer,SIGNAL( haveLoadMesh(QString) ),statusLabel,SLOT( setText(QString) ) );

        statusBar()->addWidget(statusLabel);
    }
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
    if(_DEBUG_MODE)
    {
        MyARAPAction->setChecked(true);
        viewer->LoadFace();
    }
    else
    {
        ButtonVisiable(true);
    }
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

void SurfaceMeshProcessingWindow::createButtons()
{
    int xSpacing = 200;
    int ySpacing = 200;
    int pngWidth = 150;
    int pngHigh = 150;
    int xStart = 270;
    int yStart = 70;

    PngHeButton = new QPushButton(this);
    PngHeButton->setGeometry(xStart, yStart, pngWidth, pngHigh);
    PngHeButton->setIcon((QIcon("./Pngs/he.png")));
    PngHeButton->setIconSize(QSize(pngWidth, pngHigh));
    connect(PngHeButton, SIGNAL(clicked()), this, SLOT(ChoosePng1()));

    PngMaButton = new QPushButton(this);
    PngMaButton->setGeometry(xStart + xSpacing, yStart, pngWidth, pngHigh);
    PngMaButton->setIcon((QIcon("./Pngs/ma.png")));
    PngMaButton->setIconSize(QSize(pngWidth, pngHigh));
    connect(PngMaButton, SIGNAL(clicked()), this, SLOT(ChoosePng2()));

    PngNiaoButton = new QPushButton(this);
    PngNiaoButton->setGeometry(xStart + 2*xSpacing, yStart, pngWidth, pngHigh);
    PngNiaoButton->setIcon((QIcon("./Pngs/niao.png")));
    PngNiaoButton->setIconSize(QSize(pngWidth, pngHigh));
    connect(PngNiaoButton, SIGNAL(clicked()), this, SLOT(ChoosePng3()));

    PngShuButton = new QPushButton(this);
    PngShuButton->setGeometry(xStart, yStart + ySpacing, pngWidth, pngHigh);
    PngShuButton->setIcon((QIcon("./Pngs/shu.png")));
    PngShuButton->setIconSize(QSize(pngWidth, pngHigh));
    connect(PngShuButton, SIGNAL(clicked()), this, SLOT(ChoosePng4()));

    PngTuButton = new QPushButton(this);
    PngTuButton->setGeometry(xStart + xSpacing, yStart + ySpacing, pngWidth, pngHigh);
    PngTuButton->setIcon((QIcon("./Pngs/tu.png")));
    PngTuButton->setIconSize(QSize(pngWidth, pngHigh));
    connect(PngTuButton, SIGNAL(clicked()), this, SLOT(ChoosePng5()));

    PngWaButton = new QPushButton(this);
    PngWaButton->setGeometry(xStart + 2*xSpacing, yStart + ySpacing, pngWidth, pngHigh);
    PngWaButton->setIcon((QIcon("./Pngs/wa.png")));
    PngWaButton->setIconSize(QSize(pngWidth, pngHigh));
    connect(PngWaButton, SIGNAL(clicked()), this, SLOT(ChoosePng6()));

    PngXiaButton = new QPushButton(this);
    PngXiaButton->setGeometry(xStart, yStart + 2*ySpacing, pngWidth, pngHigh);
    PngXiaButton->setIcon((QIcon("./Pngs/xia.png")));
    PngXiaButton->setIconSize(QSize(pngWidth, pngHigh));
    connect(PngXiaButton, SIGNAL(clicked()), this, SLOT(ChoosePng7()));

    PngXieButton = new QPushButton(this);
    PngXieButton->setGeometry(xStart + xSpacing, yStart + 2*ySpacing, pngWidth, pngHigh);
    PngXieButton->setIcon((QIcon("./Pngs/xie.png")));
    PngXieButton->setIconSize(QSize(pngWidth, pngHigh));
    connect(PngXieButton, SIGNAL(clicked()), this, SLOT(ChoosePng8()));

    PngYingButton = new QPushButton(this);
    PngYingButton->setGeometry(xStart + 2*xSpacing, yStart + 2*ySpacing, pngWidth, pngHigh);
    PngYingButton->setIcon((QIcon("./Pngs/ying.png")));
    PngYingButton->setIconSize(QSize(pngWidth, pngHigh));
    connect(PngYingButton, SIGNAL(clicked()), this, SLOT(ChoosePng9()));
}

void SurfaceMeshProcessingWindow::ChoosePng1()
{
    ButtonVisiable(false);

    viewer->png_idx = 1;

    MyARAPAction->setChecked(true);
    viewer->LoadFace();
    MyARAP3D();
}
void SurfaceMeshProcessingWindow::ChoosePng2()
{
    ButtonVisiable(false);

    viewer->png_idx = 2;

    MyARAPAction->setChecked(true);
    viewer->LoadFace();
    MyARAP3D();
}
void SurfaceMeshProcessingWindow::ChoosePng3()
{
    ButtonVisiable(false);

    viewer->png_idx = 3;

    MyARAPAction->setChecked(true);
    viewer->LoadFace();
    MyARAP3D();
}
void SurfaceMeshProcessingWindow::ChoosePng4()
{
    ButtonVisiable(false);

    viewer->png_idx = 4;

    MyARAPAction->setChecked(true);
    viewer->LoadFace();
    MyARAP3D();
}
void SurfaceMeshProcessingWindow::ChoosePng5()
{
    ButtonVisiable(false);

    viewer->png_idx = 5;

    MyARAPAction->setChecked(true);
    viewer->LoadFace();
    MyARAP3D();
}
void SurfaceMeshProcessingWindow::ChoosePng6()
{
    ButtonVisiable(false);

    viewer->png_idx = 6;

    MyARAPAction->setChecked(true);
    viewer->LoadFace();
    MyARAP3D();
}
void SurfaceMeshProcessingWindow::ChoosePng7()
{
    ButtonVisiable(false);

    viewer->png_idx = 7;

    MyARAPAction->setChecked(true);
    viewer->LoadFace();
    MyARAP3D();
}
void SurfaceMeshProcessingWindow::ChoosePng8()
{
    ButtonVisiable(false);

    viewer->png_idx = 8;

    MyARAPAction->setChecked(true);
    viewer->LoadFace();
    MyARAP3D();
}
void SurfaceMeshProcessingWindow::ChoosePng9()
{
    ButtonVisiable(false);

    viewer->png_idx = 9;

    MyARAPAction->setChecked(true);
    viewer->LoadFace();
    MyARAP3D();
}

void SurfaceMeshProcessingWindow::ButtonVisiable(bool visiable)
{
    PngHeButton->setVisible(visiable);
    PngMaButton->setVisible(visiable);
    PngNiaoButton->setVisible(visiable);
    PngShuButton->setVisible(visiable);
    PngTuButton->setVisible(visiable);
    PngWaButton->setVisible(visiable);
    PngXiaButton->setVisible(visiable);
    PngXieButton->setVisible(visiable);
    PngYingButton->setVisible(visiable);
}
