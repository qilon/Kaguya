#pragma once

#include "residual.h"
#include "Parameters.h"
#include "ceres\ceres.h"
#include <algorithm>
#include "utils.h"
#include "MeshData.h"
#include "qx_highlight_removal/qx_basic.h"
#include "qx_highlight_removal/qx_highlight_removal_bf.h"
//=============================================================================
using namespace std;
using namespace Eigen;
using namespace OpenMesh;
//=============================================================================
class BlackBox
{
private:
	static const string DEFAULT_PARAMS_FILENAME;

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
	static vector<Intensity> shadings;
	static vector<Intensity> lighting_weights;
	static vector< vector<Intensity> > diff_weights;

	static vector<bool> visibility;

	static OpenMesh::IO::Options mesh_read_opt;
	static OpenMesh::IO::Options mesh_write_opt;
	//=========================================================================

	/**** PRIVATE FUNCTIONS ****/

	static void initFromMesh(MyMesh& _mesh);
	static void initLocalLightings(MyMesh &_mesh, const string _image_filename, 
		const string _intrinsics_filename);

	static Intensity rgb2gray(const Color &_color);
	static void sortAdjacentVerticesAndFaces();

	static void initAlbedos(const string &mesh_filename);
	static void initDiffWeights();
	static void initLightingWeights();
	static Intensity computeDiffWeight(const unsigned int _v_idx, 
		const unsigned int _adj_v_idx);

	static void estimateSHCoeff(const ceres::Solver::Options &_options);
	static void estimateAlbedo(const ceres::Solver::Options &_options);
	static void estimateLocalLighting(const ceres::Solver::Options &_options);
	static void estimateAlbedoLocalLighting(const ceres::Solver::Options &_options);
	static void estimateShape(const ceres::Solver::Options &_options);

	static void refine();
	static void refineSHCoeff(const ceres::Solver::Options &_options);

	static void computeShading(const MyMesh &_mesh, 
		const vector<Intensity> &_sh_coeff, const unsigned int &_sh_order,
		vector<Intensity> &_shading);
	static Intensity computeShading(const Normal &_normal,
		const vector<Intensity> &_sh_coeff, const unsigned int &_sh_order);
	static void computeEstDiffuse(const vector<Color> &_albedos,
		const vector<Intensity> &_shadings, vector<Color> &_est_diffuse);
	static void computeIntensityEstDiffuseDiff(const vector<Color> &_intensity,
		const vector<Color> &_diffuse, vector<Color> &_diff);
	static void computeEstIntensity(const vector<Color> &_albedos,
		const vector<Intensity> &_shadings, const vector<Color> &_local_lightings,
		vector<Color> &_est_intensities);

	static void setMeshVertices(const vector<Vertex> &_vertices, MyMesh &_mesh);
	static void setMeshColors(const vector<Color> &_colors, MyMesh &_mesh);
	static void setMeshColors(const vector<Intensity> &_intensities, MyMesh &_mesh);

	static void updateVertices(vector<double> &_vertices_disp);

	static void updateShadings();
	static void updateAlbedos();
	static void updateLocalLightings();

	static void writeToPLY(std::string& filename, MyMesh& meshData);

	static void estimateDiffuse(const cv::Mat src_image, cv::Mat &diffuse_image);

public:
	static void initialize(int *argc, char **argv);
	static void run();
	static void save();
	static void destroy();
};
//=============================================================================
