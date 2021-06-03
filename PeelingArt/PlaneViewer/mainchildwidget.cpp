#include "mainchildwidget.h"
#include <QHBoxLayout>
#include <QScrollArea>
#include <QTabWidget>
#include <QDialog>
#include <QLabel>
#include <QRadioButton>
#include <QButtonGroup>
#include <QSlider>
#include <QPushButton>
#include "Algorithm/BPoly/BPolyline.h"


MainChildWidget::MainChildWidget(QWidget *parent)
	: QWidget(parent)
{
	setWindowTitle("Peeling");
	initViewerWindow();
}

MainChildWidget::~MainChildWidget()
{

}

void MainChildWidget::initViewerWindow()
{
	creatButtonDialog();
	createViewerDialog();

	QVBoxLayout* main_layout = new QVBoxLayout();
	main_layout->addLayout(layout_up);
	main_layout->addWidget(PlaneViewer, 5);
	this->setLayout(main_layout);
}

void MainChildWidget::creatButtonDialog()
{
	PartButton = new QPushButton(tr("Select Part"));
	FreeSplineButton = new QPushButton(tr("Free Spline"));
	CutSplineButton = new QPushButton(tr("Cut Spline"));
	CSplineButton = new QPushButton(tr("Change Spline"));
	CmeshButton = new QPushButton(tr("Change Mesh"));
	SetHandleButton = new QPushButton(tr("Set handle"));
	DeformationButton = new QPushButton(tr("Move"));

	SetWeightButton = new QPushButton(tr("Weight Set"));

	ISOk_change = new QPushButton(tr("All_Change_OK"));
	ISOk_change->setFixedSize(150, 30);

	ShowCurvedButton = new QPushButton(tr("Curved"));
	ShowCurvedButton->setFixedSize(70, 30);

	connect(PartButton, SIGNAL(clicked()), this, SLOT(SelectPart()));
	connect(FreeSplineButton, SIGNAL(clicked()), this, SLOT(DrawFreeSpline()));
	connect(CutSplineButton, SIGNAL(clicked()), this, SLOT(DrawCutSpline()));
	connect(CSplineButton, SIGNAL(clicked()), this, SLOT(ChangeSpline()));
	connect(CmeshButton, SIGNAL(clicked()), this, SLOT(ChangeMesh()));

	connect(SetHandleButton, SIGNAL(clicked()), this, SLOT(Deformationhandle()));
	connect(DeformationButton, SIGNAL(clicked()), this, SLOT(DeformationMove()));

	connect(SetWeightButton, SIGNAL(clicked()), this, SLOT(Weight_Set()));

	connect(ISOk_change, SIGNAL(clicked()), this, SLOT(SetAllChangOK()));

	connect(ShowCurvedButton, SIGNAL(clicked()), this, SLOT(ShowCurvedMesh()));

	layout_up = new QHBoxLayout();
	layout_up->addWidget(PartButton);
	layout_up->addWidget(FreeSplineButton);
	layout_up->addWidget(CutSplineButton);
	layout_up->addWidget(CSplineButton);
	layout_up->addWidget(CmeshButton);

//	layout_up->addWidget(SetHandleButton);
//	layout_up->addWidget(DeformationButton);

	layout_up->addWidget(SetWeightButton);
	layout_up->addWidget(ISOk_change);

	layout_up->addWidget(ShowCurvedButton);
}

void MainChildWidget::createViewerDialog()
{
	PlaneViewer = new PlaneViewerWidget();
}

void MainChildWidget::SelectPart()
{
	PlaneViewer->setMouseMode(PlaneViewer->PICK);
}
void MainChildWidget::DrawFreeSpline()
{
	PlaneViewer->setMouseMode(PlaneViewer->FSPLINE);
}
void MainChildWidget::DrawCutSpline()
{
	PlaneViewer->setMouseMode(PlaneViewer->CATSPLINE);
}
void MainChildWidget::ChangeSpline()
{
	PlaneViewer->setMouseMode(PlaneViewer->CHSPLINE);
}
void MainChildWidget::ChangeMesh()
{
	PlaneViewer->changeRefMesh();
}

void MainChildWidget::Deformationhandle()
{
	PlaneViewer->ClearAllPart();
	PlaneViewer->setMouseMode(PlaneViewer->DEFHANDLE);	
}
void MainChildWidget::DeformationMove()
{
	PlaneViewer->setMouseMode(PlaneViewer->DEFMOVE);
}
void MainChildWidget::Weight_Set()
{
	PlaneViewer->setMouseMode(PlaneViewer->WEIGHT);
	PlaneViewer->updateGL();
}

void MainChildWidget::SetAllChangOK()
{
	PlaneViewer->All_Change_OK = true;
}

void MainChildWidget::ShowCurvedMesh()
{
	PlaneViewer->ClearAllPart();
	if (PlaneViewer->planeview_draw_model == PlaneViewer->NOIC_PLANE)
	{
		PlaneViewer->curved_mesh = PlaneViewer->mesh;
		PlaneViewer->planeview_draw_model = PlaneViewer->CURVED;
		PlaneViewer->ComputeMaxH();
	}
	else
	{
		if (PlaneViewer->planeview_draw_model == PlaneViewer->CURVED)
		{
			PlaneViewer->planeview_draw_model = PlaneViewer->PLANE;
		}
		else
		{
			PlaneViewer->planeview_draw_model = PlaneViewer->CURVED;
		}
	}

	PlaneViewer->updateGL();

}
