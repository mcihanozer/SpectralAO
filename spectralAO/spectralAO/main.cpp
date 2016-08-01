#include "UserInterface.h"

//#pragma comment(lib, "embree.lib")

int main(int argc, char *argv[])
{
	// Evoke the king of hell, strike the death knell!
	UserInterface ui(new igl::viewer::Viewer());
	return ui.launch();
}