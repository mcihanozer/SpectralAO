// User interface class
// Handles user inputs and visualize the data

// TODO Introduce Controller to control and update parts

#ifndef USERINTERFACE_H
#define USERINTERFACE_H

#include "Controller.h"

#include <igl\viewer\ViewerPlugin.h>
#include <igl\viewer\Viewer.h>

class UserInterface : igl::viewer::ViewerPlugin
{
public:
	// CONSTRUCTORS
	UserInterface(igl::viewer::Viewer* _viewer);	// Set the viewer, invoke Controller to read mesh files

	// METHODS

	// Viewer related
	// This function is called when a keyboard key is down
	// - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
	bool key_down(int key, int modifiers);

	// Launch the viewer
	int launch()	{ return viewer->launch(); }

	void updateDisplay();	// Update the mesh in display
	void updateDisplay(const Eigen::VectorXd& color);
	
	bool takeScreenShot();

private:
	RENDERING_TYPE mRenderType;

	// Info for the viewer
	bool isOriginalMesh;

	Eigen::MatrixXd mVertInUse;
	Eigen::MatrixXi mFaceInUse;

	// Controls, updates the data, and respond UserInterface's queries
	Controller mController;
};

#endif