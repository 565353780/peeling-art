#ifndef SURFACEMESHPROCESSING_H
#define SURFACEMESHPROCESSING_H

#include <QMainWindow>
#include <QtGui>
#include <QtWidgets>

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

};

#endif // SURFACEMESHPROCESSING_H
