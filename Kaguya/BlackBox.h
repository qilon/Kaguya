#pragma once

#include "residual.h"
#include "FileIO.h"
#include "Parameters.h"
#include <Eigen\Dense>
//=============================================================================
using namespace std;
using namespace Eigen;
using namespace OpenMesh;
//=============================================================================
class BlackBox
{
private:
	static parameters::Parameters params;

	static MyMesh mesh;

	static unsigned int n_vertices;
	static unsigned int n_faces;

	static vector<Vertex> vertices;
	static vector<Normal> normals;
	static vector<Color> colors;
	static vector<Intensity> grays;

	static vector<Face> faces;
	static vector< vector<FaceIndex> > adj_faces;
	static vector< vector<VertexIndex> > adj_vertices;

	static unsigned int sh_order;

	static vector<Intensity> sh_coeff;
	static vector<Color> albedos;
	static vector<Color> lighting_variations;

	static OpenMesh::IO::Options mesh_read_opt;
	static OpenMesh::IO::Options mesh_write_opt;
	//=========================================================================

	/**** PRIVATE FUNCTIONS ****/

	static void initFromMesh(MyMesh& _mesh);

	static void computeSHFunctions(const MatrixXf &_normals, const int _sh_order,
		MatrixXf &_sh_functions);
	static VectorXf getSHCoeff(const VectorXf &_shading, const MatrixXf &_normals, 
		const int _sh_order);
	static VectorXf getShading(const MatrixXf &_normals, const VectorXf &_sh_coeff);
	static void setMeshColors(const vector<Color> &_colors, MyMesh &_mesh);
	static void setMeshColors(const vector<Intensity> &_intensities, MyMesh &_mesh);

	static Intensity rgb2gray(const Color &_color);
	static void sortAdjacentVerticesAndFaces();

public:
	static void initialize(int *argc, char **argv);
	static void run();
	static void destroy();
};
//=============================================================================
