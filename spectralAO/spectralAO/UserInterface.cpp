#include "UserInterface.h"

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;

// Set the viewer, invoke Controller to read mesh files
UserInterface::UserInterface(igl::viewer::Viewer* _viewer)
{
	mRenderType = RENDERING_TYPE::USUAL;

	// Read and set rendering type

	ifstream listFile("settings.txt");
	if (listFile)
	{
		std::string setting;

		while (getline(listFile, setting))
		{
			if (setting == "MH_BASES" || setting == "MH" || setting == "BASES" || setting == "MHB" || setting == "HARMONICS")
			{
				mRenderType = (RENDERING_TYPE)(mRenderType | RENDERING_TYPE::MH_BASES);
			}
			else if (setting == "AO")
			{
				mRenderType = (RENDERING_TYPE)(mRenderType | RENDERING_TYPE::SPECTRAL_AO);
			}
		}

		listFile.close();
	}
	else
	{
		std::cout << "\nCANNOT OPEN SETTINGS.TXT!\nRENDERING TYPE SET USUAL MESH RENDERING AUTOMATICALLY!\n";
	}

	mController.changeRenderingType(mRenderType);

	// Dumb work around...
	viewer = _viewer;
	viewer->plugins.push_back(this);
	// If it works, it's not stupid

	// Read mesh files
	if (mController.readMeshes())
	{
		isOriginalMesh = true;

		mController.getMeshAt(0, mVertInUse, mFaceInUse);

		cout << "\n\nPress [N] to get the next mesh." << endl;
		cout << "Press [B] to get the previous mesh." << endl;
		cout << "Press [M] to get the mesh with given index" << endl;
		cout << "Press [S] to take screenshots of all bases" << endl;
		cout << "Press [C] to take screenshot of current basis" << endl;

		if (mRenderType == (mRenderType | RENDERING_TYPE::MH_BASES))
		{
			cout << "Press [UP-DOWN ARROWS] to increase/decrease # eigens in use by 1." << endl;
			cout << "Press [E] to enter # eigens you want to use." << endl;
			cout << "Press [R] to reload harmonic basis given in updateMH.txt" << endl;
			cout << "Press [Q] to toggle between original and modified mesh.\n\n" << endl;
		}
		else if (mRenderType == (mRenderType | RENDERING_TYPE::SPECTRAL_AO))
		{
			viewer->core.lighting_factor = 0.f;

			cout << "Press [UP-DOWN ARROWS] to increase/decrease # eigens for AO by 1." << endl;
			cout << "Press [E] to enter # eigens you want to use." << endl;
			cout << "Press [Q] to toggle between original and modified mesh.\n\n" << endl;
		}

		updateDisplay();
	}
	else
	{
		cout << "\nMESH FILES CANNOT BE READ!" << endl;
	}
}

// Press [N] to get the next mesh
// Press [B] to get the previous mesh
// Press [M] to get the mesh with given index
// Press [UP-DOWN ARROWS] to increase/decrease # eigens in use by 1
// Press [E] to enter # eigens you want to use
// Press [Q] to toggle between original and modified mesh   // TODO 1 -> ORIGINAL, 2-> MHB, 3-> SPECTRAL AO?

// This function is called when a keyboard key is down
// - modifiers is a bitfield that might one or more of the following bits Preview3D::NO_KEY, Preview3D::SHIFT, Preview3D::CTRL, Preview3D::ALT;
bool UserInterface::key_down(int key, int modifiers)
{
	//std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;

	//if (mRenderType == (mRenderType | RENDERING_TYPE::MH_BASES))
	//{
		if ((unsigned int)key == 265)	// UP: Increase Eigenfunctions by 1		-	For this project it is 265, not 9. Find a better way
		{
			if (mRenderType == (mRenderType | RENDERING_TYPE::SPECTRAL_AO))
			{
				isOriginalMesh = false;

				Eigen::VectorXd newAo;
				mController.increaseEigenInUse(newAo);

				updateDisplay(newAo);
			}
			else
			{
				isOriginalMesh = false;

				mController.increaseEigenInUse(mVertInUse);
				updateDisplay();
			}
			
		}
		else if ((unsigned int)key == 264)	// DOWN: Decrease Eigenfunctions by 1		-	For this project it is 264, not 8. Find a better way
		{
			if (mRenderType == (mRenderType | RENDERING_TYPE::SPECTRAL_AO))
			{
				isOriginalMesh = false;

				Eigen::VectorXd newAo;
				mController.decreaseEigenInUse(newAo);

				updateDisplay(newAo);
			}
			else
			{
				isOriginalMesh = false;

				mController.decreaseEigenInUse(mVertInUse);
				updateDisplay();
			}
		}
		else if (key == 'E')	// Specific eigenfunction
		{
			if (mRenderType == (mRenderType | RENDERING_TYPE::SPECTRAL_AO))
			{
				isOriginalMesh = false;

				int newEigenNum = 0;
				std::cout << "Enter number of eigens you want to use: ";
				std::cin >> newEigenNum;

				Eigen::VectorXd newAo;
				if (mController.setSpecificEigen(newEigenNum, newAo))
				{
					updateDisplay(newAo);
				}
			}
			else
			{
				isOriginalMesh = false;

				int newEigenNum = 0;
				std::cout << "Enter number of eigens you want to use: ";
				std::cin >> newEigenNum;

				mController.setSpecificEigen(newEigenNum, mVertInUse);
				updateDisplay();
			}
		}
		else if (key == 'Q')	// Toggle between original and modified mesh
		{
			if (isOriginalMesh)	// Toggle to MHB mesh
			{
				isOriginalMesh = false;

				mController.getCurrentEigenInUse(mVertInUse);
				updateDisplay();
			}
			else	// Toggle to original mesh
			{
				isOriginalMesh = true;

				mController.getCurrentMesh(mVertInUse);
				updateDisplay();
			}
		}
		else if (key == 'R')	// Reload spectral bases
		{
			if (mController.reloadSpectralBases())
			{
				isOriginalMesh = false;

				unsigned int newEigenNum = mController.getVertexNumber();

				Eigen::VectorXd newAo;
				if (mController.setSpecificEigen(newEigenNum, newAo))
				{
					updateDisplay(newAo);
				}
			}
			else
			{
				cout << "\n\nCONTROLLER RETURNS FALSE!!!!\n\n";
			}
			
		}
	//}
	else if (key == 'N')	// Press [N] to get the next mesh
	{
		isOriginalMesh = true;
		mController.getNextMesh(mVertInUse, mFaceInUse);
		updateDisplay();
	}
	else if (key == 'B')	// Press [B] to get the previous mesh
	{
		isOriginalMesh = true;
		mController.getPreviousMesh(mVertInUse, mFaceInUse);
		updateDisplay();
	}
	else if (key == 'M')	// Press [M] to get the mesh with given index
	{
		unsigned int size = mController.getMeshSize();
		std::cout << "Enter the mesh in the range [1, " << size << "]: ";
		
		unsigned int userInput;
		std::cin >> userInput;

		userInput -= 1;	// Because C++ indices start at 0

		if (mController.getMeshAt(userInput, mVertInUse, mFaceInUse))
		{
			isOriginalMesh = true;
			updateDisplay();
		}
	}
	else if (key == 'S')	// Take screenshot
	{
		takeScreenShot();
	}
	else if (key == 'C')	// Single screenshot
	{
		std::cout << "\n\nSCREENSHOT IS TAKING!";

		int BASES[23] = {3485,3450,3250,3115,2879,1356,758,541,223,50,33,19,11,10,9,8,7,6,5,4,3,2,1};


		std::string path = "Data\\Screenshots\\New\\";
		std::string fileName = "snap";
		std::string extension = ".tga";

		int w = 1200, h = 800;

		for (int bi = 0; bi < 23; bi++)
		{
			Eigen::VectorXd newAo;
			mController.setSpecificEigen(BASES[bi], newAo);
			updateDisplay(newAo);

			viewer->core.clear_framebuffers();
			viewer->core.draw(viewer->data, viewer->opengl, false);

			std::string file = path + fileName + std::to_string(BASES[bi]) + extension;

			// SCREEN-SHOT CODE FROM: http://www.david-amador.com/2012/09/how-to-take-screenshot-in-opengl/

			//This prevents the images getting padded 
			// when the width multiplied by 3 is not a multiple of 4
			glPixelStorei(GL_PACK_ALIGNMENT, 1);

			int nSize = w*h * 3;
			// First let's create our buffer, 3 channels per Pixel
			char* dataBuffer = (char*)malloc(nSize*sizeof(char));

			if (!dataBuffer)
			{
				std::cout << "\n\nERROR IN DATA BUFFER!!!!";
				return false;
			}

			// Let's fetch them from the backbuffer	
			// We request the pixels in GL_BGR format, thanks to Berzeger for the tip
			glReadPixels((GLint)0, (GLint)0,
				(GLint)w, (GLint)h,
				GL_BGR, GL_UNSIGNED_BYTE, dataBuffer);

			//Now the file creation
			FILE *filePtr = fopen(file.c_str(), "wb");
			if (!filePtr)
			{
				std::cout << "\n\nCANNOT OPEN THE FILE!!!!";
				return false;
			}

			unsigned char TGAheader[12] = { 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
			unsigned char header[6] = { w % 256, w / 256,
				h % 256, h / 256,
				24, 0 };
			// We write the headers
			fwrite(TGAheader, sizeof(unsigned char), 12, filePtr);
			fwrite(header, sizeof(unsigned char), 6, filePtr);
			// And finally our image data
			fwrite(dataBuffer, sizeof(GLubyte), nSize, filePtr);
			fclose(filePtr);

			free(dataBuffer);

		}	// End of for

		std::cout << "\nDONE!\n";
	}

	return false;
}

bool UserInterface::takeScreenShot()
{
	std::string path = "Data\\Screenshots\\New\\";
	std::string fileName = "snap";
	std::string extension = ".tga";

	int w = 1200, h = 800;

	// SCREEN-SHOT CODE FROM: http://www.david-amador.com/2012/09/how-to-take-screenshot-in-opengl/

	//This prevents the images getting padded 
	// when the width multiplied by 3 is not a multiple of 4
	glPixelStorei(GL_PACK_ALIGNMENT, 1);

	int nSize = w*h * 3;
	// First let's create our buffer, 3 channels per Pixel
	char* dataBuffer = (char*)malloc(nSize*sizeof(char));

	if (!dataBuffer)
	{
		return false;
	}

	// Let's fetch them from the backbuffer	
	// We request the pixels in GL_BGR format, thanks to Berzeger for the tip
	glReadPixels((GLint)0, (GLint)0,
		(GLint)w, (GLint)h,
		GL_BGR, GL_UNSIGNED_BYTE, dataBuffer);

	//Now the file creation
	std::string original = path + fileName + "0" + extension;
	FILE *filePtr = fopen(original.c_str(), "wb");
	if (!filePtr)
	{
		return false;
	}

	unsigned char TGAheader[12] = { 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	unsigned char header[6] = { w % 256, w / 256,
		h % 256, h / 256,
		24, 0 };
	// We write the headers
	fwrite(TGAheader, sizeof(unsigned char), 12, filePtr);
	fwrite(header, sizeof(unsigned char), 6, filePtr);
	// And finally our image data
	fwrite(dataBuffer, sizeof(GLubyte), nSize, filePtr);
	fclose(filePtr);

	free(dataBuffer);

	std::cout << "\n\nSCREENSHOT IS TAKING!";

	int vertexNumber = mController.getVertexNumber();
	for (int si = vertexNumber; si > 0; si--)
	{
		Eigen::VectorXd newAo;
		mController.decreaseEigenInUse(newAo);
		updateDisplay(newAo);

		viewer->core.clear_framebuffers();
		viewer->core.draw(viewer->data, viewer->opengl, false);

		std::string file = path + fileName + std::to_string(vertexNumber - si + 1) + extension;

		// SCREEN-SHOT CODE FROM: http://www.david-amador.com/2012/09/how-to-take-screenshot-in-opengl/

		//This prevents the images getting padded 
		// when the width multiplied by 3 is not a multiple of 4
		glPixelStorei(GL_PACK_ALIGNMENT, 1);

		int nSize = w*h * 3;
		// First let's create our buffer, 3 channels per Pixel
		char* dataBuffer = (char*)malloc(nSize*sizeof(char));

		if (!dataBuffer)
		{
			std::cout << "\n\nERROR IN DATA BUFFER!!!!";
			return false;
		}

		// Let's fetch them from the backbuffer	
		// We request the pixels in GL_BGR format, thanks to Berzeger for the tip
		glReadPixels((GLint)0, (GLint)0,
			(GLint)w, (GLint)h,
			GL_BGR, GL_UNSIGNED_BYTE, dataBuffer);

		//Now the file creation
		FILE *filePtr = fopen(file.c_str(), "wb");
		if (!filePtr)
		{
			std::cout << "\n\nCANNOT OPEN THE FILE!!!!";
			return false;
		}

		unsigned char TGAheader[12] = { 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		unsigned char header[6] = { w % 256, w / 256,
			h % 256, h / 256,
			24, 0 };
		// We write the headers
		fwrite(TGAheader, sizeof(unsigned char), 12, filePtr);
		fwrite(header, sizeof(unsigned char), 6, filePtr);
		// And finally our image data
		fwrite(dataBuffer, sizeof(GLubyte), nSize, filePtr);
		fclose(filePtr);

		free(dataBuffer);

	}	// End of for

	std::cout << "\n\nSCREENSHOT IS DONE!";
}

// Update the mesh in display
void UserInterface::updateDisplay()
{
	viewer->data.clear();
	viewer->data.set_mesh(mVertInUse, mFaceInUse);

	if (mController.getFlag() == (RENDERING_TYPE::USUAL | RENDERING_TYPE::SPECTRAL_AO))
	{
		Eigen::VectorXd c = mController.returnAO();
		Eigen::MatrixXd color(mController.getVertexNumber(), 3);

		color.col(0) = c;
		color.col(1) = c;
		color.col(2) = c;

		viewer->data.set_colors(color);
	}

	

	viewer->core.align_camera_center(mVertInUse, mFaceInUse);
}

void UserInterface::updateDisplay(const Eigen::VectorXd& _color)
{
	viewer->data.clear();
	viewer->data.set_mesh(mVertInUse, mFaceInUse);

	//if (mController.getFlag() == (RENDERING_TYPE::USUAL | RENDERING_TYPE::SPECTRAL_AO))
	//{
		//Eigen::VectorXd c = mController.returnAO();
		Eigen::MatrixXd color(mController.getVertexNumber(), 3);

		color.col(0) = _color;
		color.col(1) = _color;
		color.col(2) = _color;

		viewer->data.set_colors(color);
	//}



	viewer->core.align_camera_center(mVertInUse, mFaceInUse);
}

//int UserInterface::launch()
//{
//	return viewer->launch();
//}