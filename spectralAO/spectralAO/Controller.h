// Controller class
// Controls, updates the data, and respond UserInterface's queries

#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "Mesh.h"

#include "Common.h"


class Controller
{
public:
	// CONSTRUCTOR
	Controller();	// Okay people! Move Along, there's nothing to see here!
	Controller(RENDERING_TYPE flag);	// Operation based on rendering type: I.E. if it's just rendering MH bases and AO values will not be loaded

	// METHODS
	bool changeRenderingType(RENDERING_TYPE newFlag);	// Change rendering type:  I.E. set it MH bases to load MH bases
	const RENDERING_TYPE getFlag()	{ return mFlag; }

	// MHB Related
	void getCurrentEigenInUse(Eigen::MatrixXd& newVertexPos);
	void increaseEigenInUse(Eigen::MatrixXd& newVertexPos);
	void decreaseEigenInUse(Eigen::MatrixXd& newVertexPos);

	bool reloadSpectralBases()	{ return mMeshList[mIndexInUse].reloadSpectralBases(); }	// Read new spectral bases (Added for the 4D Laplacian experiment)


	void increaseEigenInUse(Eigen::VectorXd& newVertexPos);
	void decreaseEigenInUse(Eigen::VectorXd& newVertexPos);
	bool setSpecificEigen(unsigned int eigenNumber, Eigen::VectorXd& newVertexPos);

	void setSpecificEigen(unsigned int eigenNumber, Eigen::MatrixXd& newVertexPos);

	// Mesh related
	bool readMeshes();	// Read meshes whose pathes are given in the file
	
	void getCurrentMesh(Eigen::MatrixXd& destVerts);	// Returns vertices
	void getCurrentMesh(Eigen::MatrixXd& destVerts, Eigen::MatrixXi& destFaces);	// Returns vertices and faces of the next mesh
	bool getNextMesh(Eigen::MatrixXd& destVerts, Eigen::MatrixXi& destFaces);	// Returns vertices and faces of the next mesh
	bool getPreviousMesh(Eigen::MatrixXd& destVerts, Eigen::MatrixXi& destFaces);	// Returns vertices and faces of the previous mesh
	bool getMeshAt(int index, Eigen::MatrixXd& destVerts, Eigen::MatrixXi& destFaces);	// Returns vertices and faces of the mesh with the given index

	unsigned int getMeshSize()	{ return mMeshList.size(); }
	unsigned int getVertexNumber()	{ return mMeshList[mIndexInUse].getVertexNumber(); }
	unsigned int getNumberOfBasisInUse()	{ return mEigenInUse; }


	const Eigen::VectorXd& returnAO(){ return mMeshList[mIndexInUse].returnAO(); }

private:
	RENDERING_TYPE mFlag;

	unsigned int mIndexInUse;
	unsigned int mEigenInUse;	// Number of eigenfunctions (Spectral bases) we are currently using

	std::vector<Mesh> mMeshList;	// Mesh list
};

#endif