#ifndef SURFACEMESHPROCESSING_H
#define SURFACEMESHPROCESSING_H

#include <QMainWindow>
#include <QtGui>
#include <QtWidgets>

extern bool _DEBUG_MODE;

class QAction;
class MainViewerWidget;

class SurfaceMeshProcessingWindow : public QMainWindow
{
	Q_OBJECT
public:
    SurfaceMeshProcessingWindow(QWidget *parent = 0);
    ~SurfaceMeshProcessingWindow();

	void open_mesh_from_main(char* filemane);

private:
	void createActions();
	void createMenus();
	void createToolBars();
	void createStatusBar();

private slots:
	void SetHandle();
	void moveHandle();
	void MyARAP3D();
	void LoadFace();
	void RefChangeOk();
	void SetEye();

private:
	// File Actions.
	QAction* openAction;
	QAction* saveAction;
	QAction* exitAction;
	QAction* saveScreenAction;

	// Menus.
	QMenu* fileMenu;
	QMenu* helpMenu;

	//label
	QLabel* statusLabel;

	//my action and toolbar------------------------------------------------------------------
	QToolBar* Peeling;
	QAction* InlineFixAction;
	QAction* moveAction;
	QAction* MyARAPAction;
	QAction* MyLoadAction;
	QAction* MyRefOKAction;
	QAction* SetEyeAction;
	//my action and toolbar------------------------------------------------------------------

private:
	MainViewerWidget* viewer;

private:
    void createButtons();
    void ButtonVisiable(bool visiable);

    QAction* PngHeAction;
    QAction* PngMaAction;
    QAction* PngNiaoAction;
    QAction* PngShuAction;
    QAction* PngTuAction;
    QAction* PngWaAction;
    QAction* PngXiaAction;
    QAction* PngXieAction;
    QAction* PngYingAction;

    QPushButton* PngHeButton;
    QPushButton* PngMaButton;
    QPushButton* PngNiaoButton;
    QPushButton* PngShuButton;
    QPushButton* PngTuButton;
    QPushButton* PngWaButton;
    QPushButton* PngXiaButton;
    QPushButton* PngXieButton;
    QPushButton* PngYingButton;

private slots:
    void ChoosePng1();
    void ChoosePng2();
    void ChoosePng3();
    void ChoosePng4();
    void ChoosePng5();
    void ChoosePng6();
    void ChoosePng7();
    void ChoosePng8();
    void ChoosePng9();

};

#endif // SURFACEMESHPROCESSING_H
