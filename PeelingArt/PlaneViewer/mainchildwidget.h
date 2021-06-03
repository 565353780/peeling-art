#ifndef MAINCHILDWIDGET_H
#define MAINCHILDWIDGET_H

#include <QWidget>
#include "PlaneViewerWidget.h"

class QHBoxLayout;
class QPushButton;

class MainChildWidget : public QWidget
{
	Q_OBJECT

public:
	MainChildWidget(QWidget *parent = 0);
	~MainChildWidget();

protected:
	void initViewerWindow();
	void creatButtonDialog();
	void createViewerDialog();

public:
	PlaneViewerWidget* PlaneViewer;

private slots:
	void SelectPart();
	void DrawFreeSpline();
	void DrawCutSpline();
	void ChangeSpline();
	void ChangeMesh();

	void Deformationhandle();
	void DeformationMove();

	void Weight_Set();

	void SetAllChangOK();

	void ShowCurvedMesh();


private:
	QHBoxLayout* layout_up;

	QPushButton* PartButton;
	QPushButton* FreeSplineButton;
	QPushButton* CutSplineButton;
	QPushButton* CSplineButton;
	QPushButton* CmeshButton;

	QPushButton* SetHandleButton;
	QPushButton* DeformationButton;
	QPushButton* SetWeightButton;

	QPushButton* ISOk_change;

	QPushButton* ShowCurvedButton;

};

#endif // MAINCHILDWIDGET_H
