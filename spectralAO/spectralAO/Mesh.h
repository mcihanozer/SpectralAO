// Mesh class
// Keeps mesh related info such as vertex position, face, MHB, and AO info

#ifndef MESH_H
#define MESH_H


#include <igl\read_triangle_mesh.h>	// For reading a mesh obviously
// (FUNNY ERROR: Compile keeps complaining when you use any Eigen related stuff unless you add an igl related header

// For spectral bases calculations
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>

#include <embree2\rtcore.h>

//using namespace Eigen;

class Mesh
{
public:
	// CONSTRUCTORS
	Mesh();	// Put your own mesh info such as test square
	Mesh(const std::string& filePath);	// Create a Mesh object directly via reading the file
	Mesh(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, std::string name = "dumbass");	// Set the mesh with vertex and face info

	// GETTER - SETTER
	unsigned int getVertexNumber(){ return mVertexNumber; }

	// METHODS
	bool readYourself(const std::string& filePath);	// Uses igl::read_triangle_mesh(). Read either OBJ or OFF
	bool readMCOFiles();

	void setYourInfo(Eigen::MatrixXd& destVerts, Eigen::MatrixXi& destFaces);	// Set given vertices and faces with your own: For visualization
	void setVertexInfo(Eigen::MatrixXd& destVerts);	// Set given vertices

	//void computerVertNormals();	// Compute per-vertex normals

	// MHB related
	void calculateSpectralBases();	// Calculate spectral bases (Manifold Harmonics Bases calculation)
	bool writeSpectralBases();	// Save spectral bases into a file	TODO Should it be private?
	bool readSpectralBases();	// Read saved spectral bases TODO Should it be private?
	bool reloadSpectralBases();	// Read new spectral bases (Added for the 4D Laplacian experiment)

	void smooth(const unsigned int givenEigenNumber, Eigen::MatrixXd& newVertexPos);

	void smoothAO(const unsigned int givenEigenNumber, Eigen::VectorXd& newAOs);

	// AO related
	bool getPerVertexAO();	// Load or ray trace per-vertex AO

	const Eigen::VectorXd& returnAO();

private:

	// METHODS

	// MHB related
	//void smooth(Eigen::MatrixXd& newVertexPos);
	bool calculateSpectralCoefficients(const Eigen::SparseMatrix<double>& M);	// Calculate coefficients in frequency space before hand (MHT operation)
	bool loadSpectralBases(const char* filePath);	// Read spectral bases (Added after adding  reloadSpectralBases(), old readSpectralBases())

	// AO related
	unsigned int addToScene(RTCScene scene, const RTCDevice device);
	void rayTracePerVertexAO(const RTCScene scene);	// Ray trace per-vertex AO
	bool loadAO();	// Read AO values from the file
	bool writeAO();	// Write AO values to the file

	// MEMBERS

	// Mesh ralted
	unsigned int mVertexNumber;	// Vertex number of the mesh

	std::string mMeshName;

	Eigen::MatrixXd mVertices;	// Vertex information of the mesh
	Eigen::MatrixXd mVertNormals;	// Per-vertex normals for AO
	Eigen::MatrixXi mFaces;	// Face information of the mesh

	// MHB related
	//unsigned int mEigenInUse;	// Number of eigenfunctions (Spectral bases) we are currently using

	// Coefficients in frequency space
	Eigen::VectorXd Xk;
	Eigen::VectorXd Yk;
	Eigen::VectorXd Zk;

	Eigen::VectorXd AOk;

	Eigen::MatrixXd Hks; // Spectral bases

	// AO related
	Eigen::VectorXd mAOs;
};

#endif