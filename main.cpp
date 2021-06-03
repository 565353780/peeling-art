#include "surfacemeshprocessingwindow.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);

	SurfaceMeshProcessingWindow mainWin;
	/*mainWin.setGeometry(100,100,mainWin.sizeHint().width(),mainWin.sizeHint().height());
	mainWin.resize( mainWin.sizeHint() );*/
//	mainWin.showMaximized();
    mainWin.showNormal();
	if( argc > 1 )
	{
		mainWin.open_mesh_from_main(argv[1]);
	}

	QFont font = app.font();
	font.setPointSize(15);
	app.setFont(font);

	return app.exec();
}

