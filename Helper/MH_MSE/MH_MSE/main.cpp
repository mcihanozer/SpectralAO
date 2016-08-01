#include "Mesh.h"

#include <fstream>
using std::ifstream;

using std::cout;

std::vector<Mesh*> meshList;

int vertNum;
Eigen::MatrixXd mVertices;	// Vertex information of the mesh
Eigen::MatrixXi mFaces;	// Face information of the mesh
Eigen::VectorXd mAOs;
std::string meshName;

// Takes the line and return each token
void tokenizer(const char* delimeters, const char* line, std::vector<std::string>& tokens)
{
	char token[255] = "\0";
	int tokenIndex = 0;
	int index = 0;
	bool isToken = false;

	// Remove delimeter above of the line
	for (int deli = 0; deli < strlen(delimeters); deli++)
	{
		while (line[index] == delimeters[deli])
		{
			index++;
		}
	}


	// Check the each character of the line
	for (index; index < strlen(line); index++)
	{
		for (int deli = 0; deli < strlen(delimeters); deli++)
		{
			if (line[index] != delimeters[deli])	// Add characters to the token until find the delimeter
			{
				token[tokenIndex++] = line[index];

				isToken = true;
			}
			else if (isToken)	// Token ends, add into the std::vector
			{
				token[tokenIndex + 1] = '\0';
				tokens.push_back(std::string(token));

				tokenIndex = 0;
				memset(token, 0, 255);

				isToken = false;
			}
		}	// End of for strlen(delimeters)

	}	// End of for strlen(line)

	if (isToken)
	{
		token[tokenIndex + 1] = '\0';
		tokens.push_back(token);
	}
}

bool loadAO(const std::string& filePath, Eigen::VectorXd& mAOs)
{
	FILE *filePtr = fopen(filePath.c_str(), "rb");

	if (filePtr)
	{
		mAOs.resize(vertNum);

		while (!feof(filePtr))
		{
			// VertexId (int), AOvalue (float)

			int vertexId;
			float aoValue;

			// Read vertex Id
			fread(&vertexId, sizeof(int), 1, filePtr);

			// Read AO value
			fread(&aoValue, sizeof(float), 1, filePtr);

			// Put the data into the vertex
			mAOs[vertexId] = aoValue;
		}

		fclose(filePtr);

		return true;
	}
	else
	{
		std::cout << "\nCANNOT OPEN " << filePath << " FOR READING AO VALUES!\nAO calculation is staring for " << filePath << " \n";
	}

	return false;
}

// Read mesh files and init. the system
void init()
{
	// 1. Open mesh file:
	// - Read mesh info (Vertices, faces etc.)
	// - Read AO information
	ifstream listFile("meshFile.txt");
	if (listFile)
	{
		std::string meshPath;
		getline(listFile, meshPath);

		cout << "\nREAD MESH";

		// Read mesh
		bool result = igl::read_triangle_mesh(meshPath, mVertices, mFaces);

		vertNum = mVertices.rows();

		// Get mesh name from the path
		std::vector<std::string> tokens;
		std::vector<std::string> nameTokens;
		tokenizer("\\", meshPath.c_str(), tokens);
		tokenizer(".", tokens[tokens.size() - 1].c_str(), nameTokens);

		meshName = nameTokens[0];


		cout << "\nREAD AO";
		// Read AO
		std::string aoPath;
		getline(listFile, aoPath);

		loadAO(aoPath, mAOs);

		listFile.close();

		// 2. Open MH files:
		// - Read bases
		// - Generate new mesh with realted bases
		// - Set the new mesh into meshList
		ifstream mhFile("mhFiles.txt");
		if (mhFile)
		{
			cout << "\nINIT MESH FILES";

			std::string mhPath;
			while (getline(mhFile, mhPath))
			{
				Mesh* newMesh = new Mesh(mVertices, mFaces, mAOs, meshName, mhPath);
				meshList.push_back(newMesh);
			}

			mhFile.close();
		}
		else
		{
			std::cout << "\nTHINGS WENT WRONG WHILE READING mhFiles.txt!";
		}

	}
	else
	{
		std::cout << "\nTHINGS WENT WRONG WHILE READING meshFile.txt!";
	}
}

//#pragma optimize( "", off )

void writeMSE(const std::vector<float>& errors, const unsigned int baseNumber)
{
	// Load MSE error info

	std::vector<std::string> coeffs;
	ifstream eiFile("errorInfo.txt");
	if (eiFile)
	{
		std::string nextCoeff;
		while (getline(eiFile, nextCoeff))
		{
			coeffs.push_back(nextCoeff);
		}

		eiFile.close();

		std::string fileName = meshList[0]->mMeshName + "_BaseNO_" + std::to_string(baseNumber) + ".txt";
		std::stringstream lineStream;

		std::ofstream filePtr(fileName);

		if (filePtr)
		{
			for (int ei = 0; ei < errors.size(); ei++)
			{
				lineStream << "ORGINAL VS COEFF. " << coeffs[ei] << ":\n" << errors[ei] << "\n\n";
			}

			// Write into the file
			std::string line = lineStream.str();
			filePtr << line;

			filePtr.close();
		}
	}
	else
	{
		std::cout << "\n\nTHINGS WENT WRONG WHILE READING errorInfo.txt!";
	}

	

}

//#pragma optimize( "", on )

void startMSE(const unsigned int baseNumber)
{
	std::vector<float> errors(meshList.size() - 1); // Because 0 is original and we don't need it
	Eigen::VectorXd orgAO = meshList[0]->returnSmoothAO();


	for (int mi = 1; mi < meshList.size(); mi++)
	{
		// 1. Calculate MSE

		float error = 0.f;
		Eigen::VectorXd nextAO = meshList[mi]->returnSmoothAO();

		for (int vi = 0; vi < vertNum; vi++)
		{
			float difference = orgAO(vi) - nextAO(vi);
			error += (difference * difference);
		}

		error /= vertNum;
		
		errors[mi-1] = error;
	}	// End of for of MSE

	// 2. Write into file (meshName_0vsMeshListPos_NumberOfBases.txt)
	writeMSE(errors, baseNumber);
}

//#pragma optimize( "", on )

void start()
{
	while (true)
	{
		int nextBasis;
		
		cout << "\n\nEnter number of vertices in [1," << vertNum << "] range: ";
		std::cin >> nextBasis;

		if (nextBasis > 0 && nextBasis <= vertNum)
		{
			// 1. Smooth out AO
			for (int mi = 0; mi < meshList.size(); mi++)
			{
				cout << "\nSMOOTHING OUT THE MESH "<<mi;
				meshList[mi]->smoothAO(nextBasis);
			}

			// 2. Calculate and write MSE
			cout << "\n\nSTART MSE CALCULATION";
			startMSE(nextBasis);
		}
	}
}

//#pragma optimize( "", off )

void main()
{
	// 1. Read mesh files and init. the system
	cout << "\nINIT the SYSTEM";
	init();

	// 2. Start main loop
	cout << "\nSTART LOOP";
	start();
}

//#pragma optimize( "", on )