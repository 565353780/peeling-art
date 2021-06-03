#include "surfacemeshprocessingwindow.h"
#include <QtWidgets/QApplication>

bool _DEBUG_MODE = 0;

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

    QDir dir("./output");
    if(!dir.exists())
    {
        dir.mkdir("../output");
        dir.mkdir("boundary");
        dir.mkdir("obj");
    }

	return app.exec();
}

//暂时屏蔽掉的输出
//可以通过全局搜索查找来恢复
//No Fix id    2
//begin_sharp_optimize    1
//end_sharp_optimize1    1
//the area is    1
