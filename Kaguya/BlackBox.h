#pragma once

#include "residual.h"
#include "FileIO.h"
#include "Parameters.h"
#include <Eigen\Dense>
#include "ceres\ceres.h"
//=============================================================================
using namespace std;
using namespace Eigen;
using namespace OpenMesh;
//=============================================================================
class BlackBox
{
private:
	static Parameters params;

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
	static unsigned int n_sh_basis;

	static vector<Intensity> sh_coeff;
	static vector<Color> albedos;
	static vector<Color> local_lightings;
	static vector< vector<Intensity> > diff_weights;

	static OpenMesh::IO::Options mesh_read_opt;
	static OpenMesh::IO::Options mesh_write_opt;
	//=========================================================================

	/**** PRIVATE FUNCTIONS ****/

	static void initFromMesh(MyMesh& _mesh);

	static Intensity rgb2gray(const Color &_color);
	static void sortAdjacentVerticesAndFaces();

	static void initDiffWeights();
	static Intensity computeDiffWeight(Color &_color1, Color &_color2, 
		Normal &_normal1, Normal &_normal2, Vertex &_vertex1, Vertex &_vertex2, 
		Intensity _color_diff_threshold, Intensity _color_diff_var, 
		Coordinate _normal_diff_var);

	static void estimateSHCoeff(const ceres::Solver::Options &_options);
	static void estimateAlbedo(const ceres::Solver::Options &_options);
	static void estimateLocalLighting(const ceres::Solver::Options &_options);
	static void estimateShape(const ceres::Solver::Options &_options);

	static void computeShading(const MyMesh &_mesh, 
		const vector<Intensity> &_sh_coeff, const unsigned int &_sh_order,
		vector<Intensity> &_shading);
	static Intensity computeShading(const Normal &_normal, 
		const vector<Intensity> &_sh_coeff, const unsigned int &_sh_order);
	static void computeEstIntensity(const vector<Color> &_albedos,
		const vector<Intensity> &_shadings, const vector<Color> &_local_lightings,
		vector<Color> &_est_intensities);

	static void setMeshVertices(const vector<Vertex> &_vertices, MyMesh &_mesh);
	static void setMeshColors(const vector<Color> &_colors, MyMesh &_mesh);
	static void setMeshColors(const vector<Intensity> &_intensities, MyMesh &_mesh);

public:
	static void initialize(int *argc, char **argv);
	static void run();
	static void save();
	static void destroy();
};
//=============================================================================
