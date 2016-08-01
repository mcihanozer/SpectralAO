#include "Mesh.h"

#include <vector>
#include <iostream>

// Takes the line and return each token
//void tokenizer(const char* delimeters, const char* line, std::vector<std::string>& tokens)
//{
//	char token[255] = "\0";
//	int tokenIndex = 0;
//	int index = 0;
//	bool isToken = false;
//
//	// Remove delimeter above of the line
//	for (int deli = 0; deli < strlen(delimeters); deli++)
//	{
//		while (line[index] == delimeters[deli])
//		{
//			index++;
//		}
//	}
//
//
//	// Check the each character of the line
//	for (index; index < strlen(line); index++)
//	{
//		for (int deli = 0; deli < strlen(delimeters); deli++)
//		{
//			if (line[index] != delimeters[deli])	// Add characters to the token until find the delimeter
//			{
//				token[tokenIndex++] = line[index];
//
//				isToken = true;
//			}
//			else if (isToken)	// Token ends, add into the std::vector
//			{
//				token[tokenIndex + 1] = '\0';
//				tokens.push_back(std::string(token));
//
//				tokenIndex = 0;
//				memset(token, 0, 255);
//
//				isToken = false;
//			}
//		}	// End of for strlen(delimeters)
//
//	}	// End of for strlen(line)
//
//	if (isToken)
//	{
//		token[tokenIndex + 1] = '\0';
//		tokens.push_back(token);
//	}
//}

// Insertion sort for sorting the eigenvalues in ascending order
// TODO Is quicksort needed?
void insertionSort(const int range, Eigen::VectorXd& eigenValues, std::vector<int>& sortedIndices)
{
	// Original algorithm: http://www.sorting-algorithms.com/insertion-sort
	for (int i = 1; i < range; i++)
	{
		for (int j = i; j > 0 && (eigenValues[j] < eigenValues[j - 1]); j--)	// errors[j] > errors[j - 1]: Descending - errors[j] < errors[j - 1]: Ascending
		{
			double tempEig = eigenValues[j];
			int tempIndex = sortedIndices[j];

			eigenValues[j] = eigenValues[j - 1];
			sortedIndices[j] = sortedIndices[j - 1];

			eigenValues[j - 1] = tempEig;
			sortedIndices[j - 1] = tempIndex;
		}
	}
}

// Set the mesh with vertex and face info
Mesh::Mesh(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, const Eigen::VectorXd& AOs, const std::string& name, const std::string& MHPath)
{
	mVertices = vertices;
	mFaces = faces;
	mAOs = AOs;

	mMeshName = name;

	mVertexNumber = mVertices.rows();

	smoothedAO.resize(mVertexNumber);
	
	readSpectralBases(MHPath);
}

// Create a Mesh object directly via reading the file
//Mesh::Mesh(const std::string& filePath)
//{
//
//	if (readYourself(filePath))
//	{
//		// Get mesh name from the path
//		std::vector<std::string> tokens;
//		std::vector<std::string> nameTokens;
//
//		tokenizer("\\", filePath.c_str(), tokens);
//		tokenizer(".", tokens[tokens.size() - 1].c_str(), nameTokens);
//
//		mMeshName = nameTokens[0];
//
//		nameTokens.clear();
//		tokens.clear();
//	}
//}

// Uses igl::read_triangle_mesh(). Read either OBJ or OFF
bool Mesh::readYourself(const std::string& filePath)
{
	bool result = igl::read_triangle_mesh(filePath, mVertices, mFaces);

	if (result)
	{
		mVertexNumber = mVertices.rows();

		// Get per-vertex normals
		//igl::per_vertex_normals(mVertices, mFaces, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, mVertNormals);
	}
	else
	{
		std::cout << "\nERROR AT READING MESH!\n";
	}

	return result;
}

// Read saved spectral bases
// TODO Should it be private?
bool Mesh::readSpectralBases(const std::string& MHPath)
{
	FILE *filePtr = fopen(MHPath.c_str(), "rb");

	if (filePtr)
	{
		Hks.resize(mVertexNumber, mVertexNumber);

		for (int col = 0; col < mVertexNumber; col++)
		{
			for (int row = 0; row < mVertexNumber; row++)
			{
				double eigen;
				fread(&eigen, sizeof(double), 1, filePtr);

				Hks(row, col) = eigen;
			}
		}

		fclose(filePtr);

		Eigen::SparseMatrix<double> M;
		igl::massmatrix(mVertices, mFaces, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);	// Compute dual area of the vertices (Hodge star 0)
		calculateSpectralCoefficients(M);
	}
	else
	{
		std::cout << "\nCANNOT OPEN " << MHPath << " FOR READING SPECTRAL BASES!\nMHB calculation is staring for " << mMeshName << " mesh.\n";
		return false;
	}

	return true;

}

// Calculate coefficients in frequency space before hand (MHT operation)
bool Mesh::calculateSpectralCoefficients(const Eigen::SparseMatrix<double>& M)
{
	// MHT
	Xk = mVertices.col(0).transpose() * M * Hks;
	Yk = mVertices.col(1).transpose() * M * Hks;
	Zk = mVertices.col(2).transpose() * M * Hks;

	AOk = mAOs.transpose() * M * Hks;

	return true;
}

void Mesh::smoothAO(const unsigned int givenEigenNumber)
{
	smoothedAO = Eigen::VectorXd::Zero(mVertexNumber);

	for (int k = 0; k < givenEigenNumber; k++)
	{
		smoothedAO = smoothedAO + (AOk(k) * Hks.col(k));
	}
}

const Eigen::VectorXd& Mesh::returnSmoothAO()
{
	return smoothedAO;
}