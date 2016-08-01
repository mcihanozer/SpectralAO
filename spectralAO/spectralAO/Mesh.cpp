#include "Common.h"

#include "Mesh.h"

#include <igl\per_vertex_normals.h>



// Put your own mesh info such as test square
Mesh::Mesh()
{
	// Say something! Reggae, reggae! Say something! Rockers, rockers!
	mVertices.resize(9, 3);
	mFaces.resize(8, 3);

	// Test mesh
	mVertices(0, 0) = -1;	mVertices(0, 1) = 1;	mVertices(0, 2) = 0;
	mVertices(1, 0) = 0;	mVertices(1, 1) = 1;	mVertices(1, 2) = 0;
	mVertices(2, 0) = 1;	mVertices(2, 1) = 1;	mVertices(2, 2) = 0;
	mVertices(3, 0) = -1;	mVertices(3, 1) = 0;	mVertices(3, 2) = 0;
	mVertices(4, 0) = 0;	mVertices(4, 1) = 0;	mVertices(4, 2) = 0;
	mVertices(5, 0) = 1;	mVertices(5, 1) = 0;	mVertices(5, 2) = 0;
	mVertices(6, 0) = -1;	mVertices(6, 1) = -1;	mVertices(6, 2) = 0;
	mVertices(7, 0) = 0;	mVertices(7, 1) = -1;	mVertices(7, 2) = 0;
	mVertices(8, 0) = 1;	mVertices(8, 1) = -1;	mVertices(8, 2) = 0;


	mFaces(0, 0) = 0;	mFaces(0, 1) = 3;	mFaces(0, 2) = 4;
	mFaces(1, 0) = 3;	mFaces(1, 1) = 6;	mFaces(1, 2) = 4;
	mFaces(2, 0) = 6;	mFaces(2, 1) = 7;	mFaces(2, 2) = 4;
	mFaces(3, 0) = 7;	mFaces(3, 1) = 8;	mFaces(3, 2) = 4;
	mFaces(4, 0) = 8;	mFaces(4, 1) = 5;	mFaces(4, 2) = 4;
	mFaces(5, 0) = 5;	mFaces(5, 1) = 2;	mFaces(5, 2) = 4;
	mFaces(6, 0) = 2;	mFaces(6, 1) = 1;	mFaces(6, 2) = 4;
	mFaces(7, 0) = 1;	mFaces(7, 1) = 0;	mFaces(7, 2) = 4;

	mMeshName = "TestMesh";

	// DO NOT FORGET TO SET VERTEX NUMBER BY HAND!
	mVertexNumber = 9;
	//mEigenInUse = mVertexNumber;

	// Get per-vertex normals
	igl::per_vertex_normals(mVertices, mFaces, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, mVertNormals);
}

// Create a Mesh object directly via reading the file
Mesh::Mesh(const std::string& filePath)
{

	if (readYourself(filePath))
	{
		// Get mesh name from the path
		std::vector<std::string> tokens;
		std::vector<std::string> nameTokens;

		tokenizer("\\", filePath.c_str(), tokens);
		tokenizer(".", tokens[tokens.size() - 1].c_str(), nameTokens);

		mMeshName = nameTokens[0];

		nameTokens.clear();
		tokens.clear();
	}
}

// Set the mesh with vertex and face info
Mesh::Mesh(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, std::string name)
{
	mVertices = vertices;
	mFaces = faces;

	mMeshName = name;

	mVertexNumber = mVertices.rows();
	//mEigenInUse = mVertexNumber;
}

// Uses igl::read_triangle_mesh(). Read either OBJ or OFF
bool Mesh::readYourself(const std::string& filePath)
{
	bool result = igl::read_triangle_mesh(filePath, mVertices, mFaces);

	if (result)
	{
		mVertexNumber = mVertices.rows();

		// Get per-vertex normals
		igl::per_vertex_normals(mVertices, mFaces, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE, mVertNormals);
	}
	else
	{
		std::cout<<"\nERROR AT READING MESH!\n";
	}

	return result;
}

bool Mesh::readMCOFiles()
{
	// For leather_armor (fixed)
	mVertexNumber = 15140;
	int faceNumber = 22288;

	mVertices.resize(mVertexNumber, 3);
	mFaces.resize(faceNumber, 3);
	mVertNormals.resize(mVertexNumber, 3);
	mAOs.resize(mVertexNumber);

	mMeshName = "mco";


	// Open files
	FILE* faceFile = fopen("Data\\MCO\\faceInfo.fi", "rb");
	FILE* vertexFile = fopen("Data\\MCO\\vertexInfo.vi", "rb");
	FILE* normalsFile = fopen("Data\\MCO\\normalInfo.ni", "rb");
	FILE* aoFile = fopen("Data\\MCO\\aoInfo.ao", "rb");

	if (faceFile && vertexFile )//&& normalsFile && aoFile)
	{
		// Read face info
		for (int fi = 0; fi < faceNumber; fi++)
		{
			int f1, f2, f3;

			fread(&f1, sizeof(int), 1, faceFile);
			fread(&f2, sizeof(int), 1, faceFile);
			fread(&f3, sizeof(int), 1, faceFile);

			mFaces(fi, 0) = f1;
			mFaces(fi, 1) = f2;
			mFaces(fi, 2) = f3;
		}
		fclose(faceFile);

		// Read vertex, normal, AO info
		for (int vi = 0; vi < mVertexNumber; vi++)
		{
			float posX, posY, posZ;
			float norX, norY, norZ;
			float ao;

			// Vertex
			fread(&posX, sizeof(float), 1, vertexFile);
			fread(&posY, sizeof(float), 1, vertexFile);
			fread(&posZ, sizeof(float), 1, vertexFile);

			mVertices(vi, 0) = posX;
			mVertices(vi, 1) = posY;
			mVertices(vi, 2) = posZ;

			//// Normal
			fread(&norX, sizeof(float), 1, normalsFile);
			fread(&norY, sizeof(float), 1, normalsFile);
			fread(&norZ, sizeof(float), 1, normalsFile);

			mVertNormals(vi, 0) = norX;
			mVertNormals(vi, 1) = norY;
			mVertNormals(vi, 2) = norZ;

			//Ao
			fread(&ao, sizeof(float), 1, aoFile);
			mAOs[vi] = ao;
		}
		fclose(vertexFile);
		fclose(normalsFile);
		fclose(aoFile);

		return true;
	}
	else
	{
		std::cout << "\nFAILED WHILE OPEN FILES!\n";
	}

	return false;
}

// Set given vertices
void Mesh::setVertexInfo(Eigen::MatrixXd& destVerts)
{
	destVerts = mVertices;
}

// Set given vertices and faces with your own: For visualization
void Mesh::setYourInfo(Eigen::MatrixXd& destVerts, Eigen::MatrixXi& destFaces)
{
	destVerts = mVertices;
	destFaces = mFaces;
}

//#include <Eigen/Core>
//#include <Eigen/Geometry>
//
//class Vec3
//{
//public:
//
//	Vec3(){}
//
//	Vec3(double _x, double _y, double _z)
//	{
//		x = _x;
//		y = _y;
//		z = _z;
//	}
//
//	Vec3(const Eigen::Vector3d& vec)
//	{
//		x = vec.coeff(0, 0);
//		y = vec.coeff(0, 1);
//		z = vec.coeff(0, 2);
//	}
//
//	bool operator < (Vec3 const& p) const
//	{
//		return	(z != p.z) ? (z<p.z) :
//				(y != p.y) ? (y<p.y) :
//				(x<p.x);
//
//		//return	(_v[2] != p_v[2]) ? (_v[2]<p._v[2]) :
//		//	(_v[1] != p._v[1]) ? (_v[1]<p._v[1]) :
//		//	(_v[0]<p._v[0]);
//	}
//
//	double x, y, z;
//};

// Compute per-vertex normals
//void Mesh::computerVertNormals()
//{
//	//mVertNormals.resize(mVertexNumber, 3);
//
//	//// A map from vert positions to normals
//	//std::map<Vec3, Vec3> p2n;
//
//	//// Pass 1: Reset
//	//for (unsigned int vi = 0; vi < mVertexNumber; vi++)
//	//{
//	//	Vec3 dumdum(0, 0, 0);
//	//	Vec3 dummy(mVertices.row(vi));
//
//	//	p2n[dummy] = dumdum;
//	//}
//
//	//// Pass 2: Cycle over faces, cumulate normals
//	//for (unsigned int fi = 0; fi < mFaces.rows(); fi++)
//	//{
//	//	// Get vertex ids
//	//	int vi[3];
//
//	//	vi[0] = mFaces.coeff(fi, 0);
//	//	vi[1] = mFaces.coeff(fi, 1);
//	//	vi[2] = mFaces.coeff(fi, 2);
//
//	//	// Get vertex positions
//	//	Eigen::Vector3d p0 = mVertices.row(vi[0]);
//	//	Eigen::Vector3d p1 = mVertices.row(vi[1]);
//	//	Eigen::Vector3d p2 = mVertices.row(vi[2]);
//
//	//	p1 -= p0;
//	//	p2 -= p0;
//
//	//	// This includes area weighting
//	//	Eigen::Vector3d faceNorm = p2.cross(p1);
//
//	//	for (int z = 0; z < 3; z++)
//	//	{
//	//		Eigen::Vector3d e1 = mVertices.row(vi[(z + 1) % 3]) - mVertices.row(vi[z]);
//	//		Eigen::Vector3d e2 = mVertices.row(vi[(z + 2) % 3]) - mVertices.row(vi[z]);
//
//	//		double wedgeAngle = angle(e1, e2);
//
//	//		p2n[Vec3(mVertices.row(vi[z]))] += faceNorm * wedgeAngle;
//	//	}
//
//	//}	// End of pass 2 for
//
//	//// pass 3: Normalize
//	//for (unsigned int vi = 0; vi < mVertexNumber; vi++)
//	//{
//	//	//mVertNormals.row(vi) = p2n[mVertices.row(vi)].normalize();
//	//}
//
//
//
//}	// End of computerVertNormals()

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

// Calculate spectral bases (Manifold Harmonics Bases calculation)
void Mesh::calculateSpectralBases()
{
	// ALGORITHM:

	// 1. Assemble positive semi - definite discrete Laplacian
	// 2. Compute eigenvectors for the Laplacian for getting MHB
	// 3. Map the bases into canonical basis

	// Details:  http://www.mcihanozer.com/wp-content/uploads/MCO-Mesh_Smoothing_with_Manifold_Harmonics.pdf

	// Laplacian, Hodge star 0, inverse of Hodge star 0, cotangent laplacian
	Eigen::SparseMatrix<double> beltrami, M, Minv, L;
	Eigen::MatrixXd eVec;	// Eigen vectors of the Laplacian (beltrami)
	Hks.resize(mVertexNumber, mVertexNumber);

	// STEP 1: Positive semi - definite discrete Laplacian
	igl::cotmatrix(mVertices, mFaces, L);	// Compute cotangent Laplace operator: #V by #V
	igl::massmatrix(mVertices, mFaces, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);	// Compute dual area of the vertices (Hodge star 0)

	// Get M-1/2
	igl::invert_diag(M, Minv);
	Minv = Minv.cwiseSqrt();

	// Set up positive semi-definite discrete Laplacian
	beltrami = Minv * L * Minv;
	beltrami = beltrami * -1;

	//  Handle numerical precision issue: http://stackoverflow.com/a/33259074
	Eigen::SparseMatrix<double> dummy = beltrami.transpose();
	beltrami = (beltrami + dummy) * 0.5;

	//-------------------------------------------------------------------------------//




	//-------------------------------------------------------------------------------//
	
	std::cout << "Eigen decomposition is started...\n";

	// Apperantly, EigenSolver<> does not give orthogonal eigenvectors even if the Laplacian is positive semi-definite
	// http://www.mathworks.com/matlabcentral/newsreader/view_thread/29459 (First entry)
	// But this is slower...
	Eigen::JacobiSVD<Eigen::MatrixXd> es(beltrami, Eigen::ComputeFullV);
	Eigen::VectorXd eigenValues = es.singularValues();	// Get eigenvalues
	eVec = es.matrixV();	// Get eigenvectors (Left eigenvectors)

	// Sort eigenvectors by increasing eigenvalues (Ascending order)
	std::vector<int> sortedIndices(mVertexNumber);
	for (int i = 0; i < mVertexNumber; i++)
	{
		sortedIndices[i] = i;
	}

	insertionSort(mVertexNumber, eigenValues, sortedIndices);

	// Map the bases into canonical basis
	// TODO probably slow
	for (int vi = 0; vi < mVertexNumber; vi++)
	{
		int next = sortedIndices[vi];
		Hks.col(vi) = Minv * eVec.col(next);
	}
	
	calculateSpectralCoefficients(M);

	// Write the bases into the file
	writeSpectralBases();

}

// Save spectral bases into a file
// TODO Should it be private?
bool Mesh::writeSpectralBases()
{
	std::string filePath = "Data\\" + mMeshName + ".hb";

	FILE* filePtr = fopen(filePath.c_str(), "wb");
	if (filePtr)
	{
		std::cout << "Spectralel basis are being written to " << filePath << " path...";

		for (int col = 0; col < mVertexNumber; col++)
		{
			for (int row = 0; row < mVertexNumber; row++)
			{
				double eigen = Hks.coeff(row, col);
				fwrite(&eigen, sizeof(double), 1, filePtr);
			}
		}

		fclose(filePtr);

		std::cout << "Done!";

		return true;
	}
	
	std::cout << "\nCANNOT OPEN "<<filePath<<" FOR WRITING SPECTRAL BASES!\n";

	return false;
}

// Read new spectral bases (Added for the 4D Laplacian experiment)
bool Mesh::reloadSpectralBases()
{
	std::cout << "\n\nReloading of Spectral Bases is starting...";

	std::ifstream newHbFile("updateMH.txt");
	if (newHbFile)
	{
		std::string path;
		std::getline(newHbFile, path);
		loadSpectralBases(path.c_str());
		return true;
	}
	else
	{
		std::cout << "\nCANNOT OPEN updateMH.txt!\n";
	}

	return false;
}

// Read saved spectral bases
// TODO Should it be private?
bool Mesh::readSpectralBases()
{
	std::string filePath = "Data\\" + mMeshName + ".hb";

	return loadSpectralBases(filePath.c_str());

	//FILE *filePtr = fopen(filePath.c_str(), "rb");

	//if (filePtr)
	//{
	//	Hks.resize(mVertexNumber, mVertexNumber);

	//	for (int col = 0; col < mVertexNumber; col++)
	//	{
	//		for (int row = 0; row < mVertexNumber; row++)
	//		{
	//			double eigen;
	//			fread(&eigen, sizeof(double), 1, filePtr);

	//			Hks(row, col) = eigen;
	//		}
	//	}

	//	fclose(filePtr);

	//	Eigen::SparseMatrix<double> M;
	//	igl::massmatrix(mVertices, mFaces, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);	// Compute dual area of the vertices (Hodge star 0)
	//	calculateSpectralCoefficients(M);
	//}
	//else
	//{
	//	std::cout << "\nCANNOT OPEN " << filePath << " FOR READING SPECTRAL BASES!\nMHB calculation is staring for "<<mMeshName<<" mesh.\n";
	//	calculateSpectralBases();
	//}

	//return true;
	
}

// Read spectral bases (Added after adding  reloadSpectralBases(), old readSpectralBases())
bool Mesh::loadSpectralBases(const char* filePath)
{
	FILE *filePtr = fopen(filePath, "rb");

	if (filePtr)
	{
		std::cout << "\nRead spectral bases located at" << filePath << std::endl;

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

		return true;
	}
	else
	{
		std::cout << "\nCANNOT OPEN " << filePath << " FOR READING SPECTRAL BASES!\nMHB calculation is staring for " << mMeshName << " mesh.\n";
		calculateSpectralBases();
	}

	return false;
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

void Mesh::smooth(const unsigned int givenEigenNumber, Eigen::MatrixXd& newVertexPos)
{
	newVertexPos.resize(mVertexNumber, 3);

	// Smoothing and MHT-1
	Eigen::VectorXd newX = Eigen::VectorXd::Zero(mVertexNumber);
	Eigen::VectorXd newY = Eigen::VectorXd::Zero(mVertexNumber);
	Eigen::VectorXd newZ = Eigen::VectorXd::Zero(mVertexNumber);

	for (int k = 0; k < givenEigenNumber; k++)
	{
		newX = newX + (Xk(k) * Hks.col(k));
		newY = newY + (Yk(k) * Hks.col(k));
		newZ = newZ + (Zk(k) * Hks.col(k));
	}

	// Map new vertices
	for (int vi = 0; vi < mVertexNumber; vi++)
	{
		newVertexPos(vi, 0) = newX(vi);
		newVertexPos(vi, 1) = newY(vi);
		newVertexPos(vi, 2) = newZ(vi);
	}
}

void Mesh::smoothAO(const unsigned int givenEigenNumber, Eigen::VectorXd& newAOs)
{
	//newAOs.resize(mVertexNumber);

	//Eigen::VectorXd newAO = Eigen::VectorXd::Zero(mVertexNumber);

	newAOs = Eigen::VectorXd::Zero(mVertexNumber);

	for (int k = 0; k < givenEigenNumber; k++)
	{
		newAOs = newAOs + (AOk(k) * Hks.col(k));
	}

	//// Map new AOs
	//for (int vi = 0; vi < mVertexNumber; vi++)
	//{
	//	newAOs[vi] = newAO[vi];
	//}
}

#include <embree2\rtcore_ray.h>

// Data structures for sending vertex, and face information to Embree
struct Vertex	{ float x, y, z, a; };	// For vertices
struct Triangle { int v0, v1, v2; };	// For faces


unsigned int Mesh::addToScene(RTCScene scene, const RTCDevice device)
{
	// Create RTC mesh

	int faces = mFaces.rows();
	int verts = mVertices.rows();

	unsigned int geomID = rtcNewTriangleMesh(scene, RTC_GEOMETRY_STATIC, faces, verts, 1);

	// Set vertices
	Vertex* vertices = (Vertex*)rtcMapBuffer(scene, geomID, RTC_VERTEX_BUFFER);
	for (int i = 0; i < mVertices.rows(); i++)
	{
		vertices[i].x = mVertices.coeff(i, 0);
		vertices[i].y = mVertices.coeff(i, 1);
		vertices[i].z = mVertices.coeff(i, 2);
	}
	rtcUnmapBuffer(scene, geomID, RTC_VERTEX_BUFFER);

	// Set faces/triangles
	Triangle* triangles = (Triangle*)rtcMapBuffer(scene, geomID, RTC_INDEX_BUFFER);

	for (int i = 0; i < mFaces.rows(); i++)
	{
		triangles[i].v0 = mFaces.coeff(i, 0);
		triangles[i].v1 = mFaces.coeff(i, 1);
		triangles[i].v2 = mFaces.coeff(i, 2);
	}
	rtcUnmapBuffer(scene, geomID, RTC_INDEX_BUFFER);

	return geomID;
}

bool Mesh::loadAO()
{
	std::string filePath = "Data\\" + mMeshName + ".ao";

	FILE *filePtr = fopen(filePath.c_str(), "rb");

	if (filePtr)
	{

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
		std::cout << "\nCANNOT OPEN " << filePath << " FOR READING AO VALUES!\nAO calculation is staring for " << mMeshName << " mesh.\n";
	}

	return false;
}

void error_handler(const RTCError code, const char *str) {
	printf("Embree: ");
	switch (code) {
	case RTC_UNKNOWN_ERROR:
		printf("RTC_UNKNOWN_ERROR");
		break;
	case RTC_INVALID_ARGUMENT:
		printf("RTC_INVALID_ARGUMENT");
		break;
	case RTC_INVALID_OPERATION:
		printf("RTC_INVALID_OPERATION");
		break;
	case RTC_OUT_OF_MEMORY:
		printf("RTC_OUT_OF_MEMORY");
		break;
	case RTC_UNSUPPORTED_CPU:
		printf("RTC_UNSUPPORTED_CPU");
		break;
	case RTC_CANCELLED:
		printf("RTC_CANCELLED");
		break;
	default:
		printf("invalid error code");
		break;
	}
	if (str) {
		printf(" (");
		while (*str) putchar(*str++);
		printf(")\n");
	}
	exit(1);
}

// Load or ray trace per-vertex AO
bool Mesh::getPerVertexAO()
{
	mAOs.resize(mVertexNumber);

	if (!loadAO())
	{
		//rtcInit();

		RTCDevice device = rtcNewDevice(NULL);
		rtcDeviceSetErrorFunction(device, error_handler);

		RTCScene scene = rtcDeviceNewScene(device, RTC_SCENE_STATIC, RTC_INTERSECT1);

		addToScene(scene, device);

		rtcCommit(scene);

		// Calculate per vertex AO for the current pose
		rayTracePerVertexAO(scene);

		rtcDeleteScene(scene);

		rtcExit();

		writeAO();	// Write AO factor (just a float) values into "modelName_animationName.ao" file
	}

	return true;
}

#define EPSILON 1e-4f
#define SAMPLE_SIZE 10000	// Number of ray samples

// Generate a random number between 0 and 1
// return an uniform number in [0,1].
double unifRand()
{
	return rand() / double(RAND_MAX);
}

// Return direction cosine-weighted hemisphere
Eigen::Vector3d getDirectionCosine()
{
	float r1 = unifRand();
	float r2 = unifRand();

	float theta = acosf(sqrt(r1));
	float phi = 2 * M_PI * r2;

	Eigen::Vector3d direction;
	direction[0] = sinf(theta) * cosf(phi);		// X
	direction[1] = sinf(theta) * sinf(phi);		// Y
	direction[2] = cosf(theta);					// Z

	direction.normalize();

	return direction;
}

//	Converts the direction into the world coordinates
// Generates a new coordinate system with using "a" parameter
// In the project "a" is vertex normals
void coordinateSystem(const Eigen::Vector3d& a, Eigen::Vector3d& b, Eigen::Vector3d& c)
{
	// Developed by refering Nori Ray Tracer source code

	if (std::abs(a[0]) > std::abs(a[1]))
	{
		float invLen = 1.0f / std::sqrt(a[0] * a[0] + a[2] * a[2]);

		c[0] = a[2] * invLen;
		c[1] = 0.f;
		c[2] = -a[0] * invLen;
	}
	else
	{
		float invLen = 1.f / std::sqrt(a[1] * a[1] + a[2] * a[2]);

		c[0] = 0.f;
		c[1] = a[2] * invLen;
		c[2] = -a[1] * invLen;
	}

	b = c.cross(a); //cross(c, a);
}

//	Converts the direction into the world coordinates
// Creates a new coordinate system by using vertex normal,
// and aligns the direction picked by "getDirectionCosine()" 
// within this new coordinate syste
Eigen::Vector3d toWorld(const Eigen::Vector3d& direction, const Eigen::Vector3d& n)
{
	// Developed by refering Nori Ray Tracer source code

	// 1. Calculate tangent vectors (Nori)
	Eigen::Vector3d s, t;
	coordinateSystem(n, s, t);

	// 2. Convert from local coordinates to world coordinates
	Eigen::Vector3d newDir = s * direction[0] + t * direction[1] + n * direction[2];
	newDir.normalize();

	return newDir;

}

void Mesh::rayTracePerVertexAO(const RTCScene scene)
{
	//	 ALOGIRTH	\\
	// For each vertices
		// Get vertex position

		// For each sample
			// Pick a direction (Cosine-weighted)
			// Convert direction into the world space
			// Generate the ray
			// Trace the ray
			// With respect to visibility add 1 or continue with the next sample
		// End of for each sample

		// Assign the result with respect to the current pose
	// End of for each vertex

	for (int i = 0; i < mVertexNumber; i++)
	{
		float aoFactor = 0.f;
		Eigen::Vector3d vertN = mVertNormals.row(i);
		Eigen::Vector3d rayOrigin = mVertices.row(i) + (mVertNormals.row(i) * EPSILON);

		for (int j = 0; j < SAMPLE_SIZE; j++)
		{
			// Pick a direction (Cosine-weighted). and convert it into the world space
			Eigen::Vector3d rayDir = toWorld(getDirectionCosine(), vertN);

			// Generate the ray
			RTCRay ray;

			ray.org[0] = rayOrigin[0];
			ray.org[1] = rayOrigin[1];
			ray.org[2] = rayOrigin[2];

			ray.dir[0] = rayDir[0];
			ray.dir[1] = rayDir[1];
			ray.dir[2] = rayDir[2];

			ray.tnear = EPSILON;
			ray.tfar = INFINITY;	//rayLength;
			ray.geomID = RTC_INVALID_GEOMETRY_ID;
			ray.primID = RTC_INVALID_GEOMETRY_ID;
			ray.mask = -1;
			ray.time = 0;

			rtcOccluded(scene, ray);	// Trace the ray

			if (ray.geomID != 0)	// No hit
			{
				aoFactor++;
			}

		}	// End of for each sample
		
		aoFactor /= SAMPLE_SIZE;
		mAOs[i] = aoFactor;

	}	// End of for each vertices

}	// End of rayTracePerVertexAO()

bool Mesh::writeAO()
{
	std::string filePath = "Data\\" + mMeshName + ".ao";

	FILE* filePtr = fopen(filePath.c_str(), "wb");

	if (filePtr)
	{
		for (int vertexId = 0; vertexId < mVertexNumber; vertexId++)	// For each vertices
		{
			// Write vertex id
			fwrite(&vertexId, sizeof(int), 1, filePtr);

			// Write AO value
			float ao = mAOs[vertexId];
			fwrite(&ao, sizeof(float), 1, filePtr);
		}

		fclose(filePtr);
	}
	else
	{
		std::cout << "\nCANNOT OPEN " << filePath << " FOR WRITING AO VALUES! for " << mMeshName << " MESH!!!\n";
	}

	return false;
}

const Eigen::VectorXd& Mesh::returnAO()
{
	return mAOs;
}