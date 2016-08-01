#ifndef MESH_H
#define MESH_H

#include <string>

#include <igl\read_triangle_mesh.h>	// For reading a mesh obviously
#include <igl/massmatrix.h>

class Mesh
{
public:
	//Mesh(const std::string& filePath);	// Create a Mesh object directly via reading the file
	Mesh(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, const Eigen::VectorXd& AOs, const std::string& name, const std::string& MHPath);	// Set the mesh with vertex and face info

	bool readYourself(const std::string& filePath);	// Uses igl::read_triangle_mesh(). Read either OBJ or OFF
	bool readSpectralBases(const std::string& MHPath);	// Read saved spectral bases TODO Should it be private?

	void smoothAO(const unsigned int givenEigenNumber);

	const Eigen::VectorXd& returnSmoothAO();

	std::string mMeshName;

private:

	bool calculateSpectralCoefficients(const Eigen::SparseMatrix<double>& M);	// Calculate coefficients in frequency space before hand (MHT operation)

	unsigned int mVertexNumber;	// Vertex number of the mesh

	Eigen::MatrixXd mVertices;	// Vertex information of the mesh
	Eigen::MatrixXi mFaces;	// Face information of the mesh

	// Coefficients in frequency space
	Eigen::VectorXd Xk;
	Eigen::VectorXd Yk;
	Eigen::VectorXd Zk;

	Eigen::VectorXd AOk;

	Eigen::MatrixXd Hks; // Spectral bases

	// AO related
	Eigen::VectorXd mAOs;
	Eigen::VectorXd smoothedAO;
};

#endif