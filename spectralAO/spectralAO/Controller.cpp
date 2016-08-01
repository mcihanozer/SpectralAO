#include "Controller.h"

#include <fstream>
using std::ifstream;

Controller::Controller()
{ 
	mIndexInUse = 0;
	mEigenInUse = 0;
	mFlag = RENDERING_TYPE::USUAL;
}


// Operation based on rendering type: I.E. if it's just rendering MH bases and AO values will not be loaded
Controller::Controller(RENDERING_TYPE flag)
{
	mIndexInUse = 0;

	// TODO NOT SURE ABOUT THIS USAGE
	if (flag == RENDERING_TYPE::USUAL)
	{
		Mesh m;
		mMeshList.push_back(m);

		mEigenInUse = m.getVertexNumber();
	}
	else if ((flag == RENDERING_TYPE::MH_BASES) || (flag == (RENDERING_TYPE::USUAL | RENDERING_TYPE::MH_BASES)))
	{

	}
}

// Change rendering type:  I.E. set it MH bases to load MH bases
bool Controller::changeRenderingType(RENDERING_TYPE newFlag)
{
	mFlag = newFlag;

	return true;
}

// Read meshes whose pathes are given in the file
bool Controller::readMeshes()
{
	// Test scene
	//Mesh m;
	//mMeshList.push_back(m);
	//return true;

	//Mesh m;
	//m.readMCOFiles();
	//m.readSpectralBases();
	//mMeshList.push_back(m);
	//return true;

	ifstream listFile("meshList.txt");
	if (listFile)
	{
		std::string filePath;

		while (getline(listFile, filePath))
		{
			Mesh nextMesh(filePath);

			if (mFlag == (RENDERING_TYPE::USUAL | RENDERING_TYPE::MH_BASES))
			{
				nextMesh.readSpectralBases();
			}
			else if (mFlag == (RENDERING_TYPE::USUAL | RENDERING_TYPE::SPECTRAL_AO))
			{
				nextMesh.getPerVertexAO();
				nextMesh.readSpectralBases();
			}

			mMeshList.push_back(nextMesh);
			
		}

		listFile.close();

		std::cout << "\nNUMBER OF MESHES: "<<mMeshList.size()<<std::endl;

		return true;
	}
	else
	{
		std::cout << "\nCANNOT OPEN MESTLIST.TXT!\n";
	}

	return false;
}

// Returns vertices
void Controller::getCurrentMesh(Eigen::MatrixXd& destVerts)
{
	mMeshList[mIndexInUse].setVertexInfo(destVerts);
}

// Returns vertices and faces of the next mesh
void Controller::getCurrentMesh(Eigen::MatrixXd& destVerts, Eigen::MatrixXi& destFaces)
{
	mMeshList[mIndexInUse].setYourInfo(destVerts, destFaces);
}

// Returns vertices and faces of the next mesh
bool Controller::getNextMesh(Eigen::MatrixXd& destVerts, Eigen::MatrixXi& destFaces)
{
	mIndexInUse = (mIndexInUse + 1) % mMeshList.size();	// Keep index in the range
	// For getting rid of "if (index >= 0 && index < mMeshList.size())" at getMeshAt() ?
	// getMeshAt(mIndexInUse, destVerts, destFaces);
	mMeshList[mIndexInUse].setYourInfo(destVerts, destFaces);
	mEigenInUse = mMeshList[mIndexInUse].getVertexNumber();

	return true;
}

// Returns vertices and faces of the previous mesh
bool Controller::getPreviousMesh(Eigen::MatrixXd& destVerts, Eigen::MatrixXi& destFaces)
{
	mIndexInUse = (mIndexInUse + mMeshList.size() - 1) % mMeshList.size();
	// For getting rid of "if (index >= 0 && index < mMeshList.size())" at getMeshAt() ?
	// getMeshAt(mIndexInUse, destVerts, destFaces);
	mMeshList[mIndexInUse].setYourInfo(destVerts, destFaces);
	mEigenInUse = mMeshList[mIndexInUse].getVertexNumber();
	
	return true;
}

// Returns vertices and faces of the mesh with the given index
bool Controller::getMeshAt(int index, Eigen::MatrixXd& destVerts, Eigen::MatrixXi& destFaces)
{
	if (index >= 0 && index < mMeshList.size())
	{
		mIndexInUse = index;

		mMeshList[mIndexInUse].setYourInfo(destVerts, destFaces);
		mEigenInUse = mMeshList[mIndexInUse].getVertexNumber();

		return true;
	}

	std::cout << "\nTHERE IS NO SUCH MESH WITH THE GIVEN INDEX \" " << index + 1 << "\"!\n";
	return false;
}

void Controller::getCurrentEigenInUse(Eigen::MatrixXd& newVertexPos)
{
	mMeshList[mIndexInUse].smooth(mEigenInUse, newVertexPos);
}

void Controller::increaseEigenInUse(Eigen::MatrixXd& newVertexPos)
{
	if (mEigenInUse < mMeshList[mIndexInUse].getVertexNumber())
	{
		mEigenInUse += 1;
		mMeshList[mIndexInUse].smooth(mEigenInUse, newVertexPos);

		std::cout << mEigenInUse << " eigens are in use.\n";
	}
}

void Controller::decreaseEigenInUse(Eigen::MatrixXd& newVertexPos)
{
	if (mEigenInUse > 1)
	{
		mEigenInUse -= 1;
		mMeshList[mIndexInUse].smooth(mEigenInUse, newVertexPos);

		std::cout << mEigenInUse << " eigens are in use.\n";
	}
}

void Controller::setSpecificEigen(unsigned int eigenNumber, Eigen::MatrixXd& newVertexPos)
{
	if (eigenNumber > 0 && eigenNumber <= mMeshList[mIndexInUse].getVertexNumber())
	{
		mEigenInUse = eigenNumber;
		mMeshList[mIndexInUse].smooth(mEigenInUse, newVertexPos);

		std::cout << mEigenInUse << " eigens are in use.\n";
	}
	else
	{
		std::cout << "Chosen eigen should be in [1, " << mMeshList[mIndexInUse].getVertexNumber() << "] range!\n";
	}
}

void Controller::increaseEigenInUse(Eigen::VectorXd& newAo)
{
	if (mEigenInUse < mMeshList[mIndexInUse].getVertexNumber())
	{
		mEigenInUse += 1;
	}

	mMeshList[mIndexInUse].smoothAO(mEigenInUse, newAo);
	std::cout << mEigenInUse << " eigens are in use.\n";

}

void Controller::decreaseEigenInUse(Eigen::VectorXd& newAo)
{
	if (mEigenInUse > 1)
	{
		mEigenInUse -= 1;
	}

	mMeshList[mIndexInUse].smoothAO(mEigenInUse, newAo);
	std::cout << mEigenInUse << " eigens are in use.\n";
}

bool Controller::setSpecificEigen(unsigned int eigenNumber, Eigen::VectorXd& newVertexPos)
{
	if (eigenNumber > 0 && eigenNumber <= mMeshList[mIndexInUse].getVertexNumber())
	{
		mEigenInUse = eigenNumber;
		mMeshList[mIndexInUse].smoothAO(mEigenInUse, newVertexPos);

		std::cout << mEigenInUse << " eigens are in use.\n";

		return true;
	}
	else
	{
		std::cout << "Chosen eigen should be in [1, " << mMeshList[mIndexInUse].getVertexNumber() << "] range!\n";
	}

	return false;
}