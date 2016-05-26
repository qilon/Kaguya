#include "BlackBox.h"
//=============================================================================
const string BlackBox::DEFAULT_PARAMS_FILENAME =
	"C:/Users/Qi/Documents/GitHub/config/Kaguya/KaguyaConfig.yml";

Parameters BlackBox::params;

MyMesh BlackBox::mesh;

MeshData BlackBox::extra_mesh;

unsigned int BlackBox::n_vertices;
unsigned int BlackBox::n_faces;

vector<Vertex> BlackBox::vertices;
vector<Normal> BlackBox::normals;
vector<Color> BlackBox::colors;
vector<Intensity> BlackBox::grays;

vector<Face> BlackBox::faces;
vector< vector<FaceIndex> > BlackBox::adj_faces;
vector< vector<VertexIndex> > BlackBox::adj_vertices;

unsigned int BlackBox::sh_order;
unsigned int BlackBox::n_sh_basis;

vector<Intensity> BlackBox::sh_coeff;
vector<Color> BlackBox::albedos;
vector<Color> BlackBox::local_lightings;
vector<Intensity> BlackBox::shadings;
vector<Intensity> BlackBox::lighting_weights;
vector< vector<Intensity> > BlackBox::diff_weights;

vector<bool> BlackBox::visibility;
cv::Mat BlackBox::color_image;
Eigen::Matrix3d BlackBox::K;

OpenMesh::IO::Options BlackBox::mesh_read_opt;
OpenMesh::IO::Options BlackBox::mesh_write_opt;
//=============================================================================
void BlackBox::initialize(int *argc, char **argv)
{
	string params_filename;

	if (argc[0] > 1)
	{
		params_filename = argv[1];
	}
	else
	{
		params_filename = DEFAULT_PARAMS_FILENAME;
	}

	params.load(params_filename);

	mesh.request_vertex_colors();
	mesh.request_face_normals();
	mesh.request_vertex_normals();

	mesh_read_opt += OpenMesh::IO::Options::VertexColor;
	mesh_read_opt += OpenMesh::IO::Options::VertexNormal;

	if (!readMesh(mesh, params.input_mesh_filename.c_str(), mesh_read_opt))
	{
		exit(0);
	}

	mesh.update_normals();

	mesh_write_opt += OpenMesh::IO::Options::VertexColor;
	//mesh_write_opt += OpenMesh::IO::Options::VertexNormal;
	//mesh_write_opt += OpenMesh::IO::Options::Custom;


	if (params.save_mesh_binary)
	{
		mesh_write_opt += OpenMesh::IO::Options::Binary;
		mesh_write_opt += OpenMesh::IO::Options::LSB;
	}

	sh_order = params.sh_order;
	n_sh_basis = (unsigned int)pow(sh_order + 1, 2);

	initFromMesh(mesh);

	if (!params.input_extra_mesh_filename.empty())
	{
		extra_mesh.readPLY(params.input_extra_mesh_filename);
	}

	initAlbedos();

	shadings.resize(n_vertices, 0);

	Color black = { 0.0, 0.0, 0.0 };
	local_lightings.resize(n_vertices, black);

	initVisibility(mesh, params.input_image_filename,
		params.input_intrinsics_filename);

	if (params.use_specular_weight)
	{
		initLocalLightings(mesh, params.input_image_filename,
			params.input_intrinsics_filename);

		initLightingWeights();
	}
}
//=============================================================================
void BlackBox::run()
{
	if (params.estimate_sh)
	{
		estimateSHCoeff(params.ceres_solver.options[0],
			params.input_extra_mesh_filename.empty());
	}

	updateShadings();
	initDiffWeights();

	if (params.input_extra_mesh_filename.empty())
	{
		updateAlbedos();
	}

	if (params.combine_albedo_lighting)
	{
		if (params.estimate_local_lighting)
		{
			updateLocalLightings();
			estimateAlbedoLocalLighting(params.ceres_solver.options[1]);
		}
	}
	else
	{
		if (params.input_extra_mesh_filename.empty())
		{
			estimateAlbedo(params.ceres_solver.options[1]);
		}

		if (params.reestimate_sh_coeff_with_albedo)
		{
			estimateSHCoeff(params.ceres_solver.options[0], false);
		}

		if (params.estimate_local_lighting)
		{
			if (params.specular_from_image)
			{
				updateColorsFromImage();
			}

			updateLocalLightings();
			estimateLocalLighting(params.ceres_solver.options[2]);
		}
	}

	if (params.estimate_shape)
	{
		estimateShape(params.ceres_solver.options[3]);
	}

	if (params.refine_solution)
	{
		refine();
	}
}
//=============================================================================
void BlackBox::save()
{
	if (params.save_sh_coeff)
	{
		cout << "Writing sh coeff: " << params.output_sh_coeff_filename << endl;

		writeVector(sh_coeff, params.output_sh_coeff_filename.c_str());
	}

	setMeshVertices(vertices, mesh);
	mesh.update_normals();

	if (params.save_orig_intensity)
	{
		cout << "Writing mesh: " << params.output_orig_intensity_mesh_filename
			<< endl;

		writeMesh(mesh, params.output_orig_intensity_mesh_filename, mesh_write_opt);
	}

	if (params.save_albedo)
	{
		setMeshColors(albedos, mesh);

		cout << "Writing mesh: " << params.output_albedo_mesh_filename << endl;

		writeMesh(mesh, params.output_albedo_mesh_filename, mesh_write_opt);
	}

	if (params.save_local_lighting)
	{
		setMeshColors(local_lightings, mesh);

		cout << "Writing mesh: " << params.output_local_lighting_mesh_filename
			<< endl;

		writeMesh(mesh, params.output_local_lighting_mesh_filename, mesh_write_opt);
	}

	vector<Intensity> shadings;
	computeShading(mesh, sh_coeff, sh_order, shadings);

	if (params.save_shading)
	{
		setMeshColors(shadings, mesh);

		cout << "Writing mesh: " << params.output_shading_mesh_filename << endl;

		writeMesh(mesh, params.output_shading_mesh_filename, mesh_write_opt);
	}

	vector<Color> est_diffuse;
	computeEstDiffuse(albedos, shadings, est_diffuse);

	if (params.save_est_diffuse)
	{
		setMeshColors(est_diffuse, mesh);

		cout << "Writing mesh: " << params.output_est_diffuse_mesh_filename
			<< endl;

		writeMesh(mesh, params.output_est_diffuse_mesh_filename, mesh_write_opt);
	}

	if (params.save_diffuse_diff)
	{
		vector<Color> diff;
		computeIntensityEstDiffuseDiff(colors, est_diffuse, diff);
		setMeshColors(diff, mesh);

		cout << "Writing mesh: " << params.output_diffuse_diff_mesh_filename
			<< endl;

		writeMesh(mesh, params.output_diffuse_diff_mesh_filename, mesh_write_opt);
	}

	if (params.save_est_intensity)
	{
		vector<Color> est_intensities;
		computeEstIntensity(albedos, shadings, local_lightings, est_intensities);

		setMeshColors(est_intensities, mesh);

		cout << "Writing mesh: " << params.output_est_intensity_mesh_filename
			<< endl;

		writeMesh(mesh, params.output_est_intensity_mesh_filename, mesh_write_opt);
	}

	if (params.save_combined)
	{
		writeToPLY(params.output_combined_mesh_filename, mesh);
	}
}
//=============================================================================
void BlackBox::destroy()
{

}
//=============================================================================
void BlackBox::initFromMesh(MyMesh& _mesh)
{
	n_vertices = (unsigned int)_mesh.n_vertices();
	n_faces = (unsigned int)_mesh.n_faces();

	vertices.resize(n_vertices);
	normals.resize(n_vertices);
	colors.resize(n_vertices);
	grays.resize(n_vertices);

	faces.resize(n_faces);
	adj_faces.resize(n_vertices);
	adj_vertices.resize(n_vertices);

	// Copy vertex coordinates, normal, color, gray, adjacent vertices and 
	// adjacent faces
	MyMesh::ConstVertexIter v_it;
	MyMesh::ConstVertexIter v_end(_mesh.vertices_end());
	for (v_it = _mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		int v_idx = v_it->idx();

		Vertex vertex(3);
		MyMesh::Point v = vector_cast<MyMesh::Point>(_mesh.point(*v_it));
		vertex[0] = v[0];
		vertex[1] = v[1];
		vertex[2] = v[2];
		vertices[v_idx] = std::move(vertex);

		Normal normal(3);
		MyMesh::Normal n = vector_cast<MyMesh::Normal>(_mesh.normal(*v_it));
		normal[0] = n[0];
		normal[1] = n[1];
		normal[2] = n[2];
		normals[v_idx] = std::move(normal);

		Color color(3);
		MyMesh::Color c = vector_cast<MyMesh::Color>(_mesh.color(*v_it));
		color[0] = (double)c[0] / 255.0;
		color[1] = (double)c[1] / 255.0;
		color[2] = (double)c[2] / 255.0;

		Intensity gray = rgb2gray(color);
		grays[v_idx] = gray;

		colors[v_idx] = std::move(color);

		vector<VertexIndex> v_idxs;
		MyMesh::ConstVertexVertexIter vv_it;
		MyMesh::ConstVertexVertexIter vv_end = _mesh.vv_end(*v_it);
		for (vv_it = _mesh.vv_begin(*v_it); vv_it != vv_end; ++vv_it)
		{
			VertexIndex vv_idx = vv_it->idx();
			v_idxs.push_back(vv_idx);
		}
		adj_vertices[v_idx] = std::move(v_idxs);

		vector<FaceIndex> f_idxs;
		MyMesh::ConstVertexFaceIter vf_it;
		MyMesh::ConstVertexFaceIter vf_end = _mesh.vf_end(*v_it);
		for (vf_it = _mesh.vf_begin(*v_it); vf_it != vf_end; ++vf_it)
		{
			FaceIndex f_idx = vf_it->idx();
			f_idxs.push_back(f_idx);
		}
		adj_faces[v_idx] = std::move(f_idxs);
	}

	// Copy face vertices
	MyMesh::ConstFaceIter f_it;
	MyMesh::ConstFaceIter f_end = _mesh.faces_end();
	for (f_it = _mesh.faces_begin(); f_it != f_end; ++f_it)
	{
		FaceIndex f_idx = f_it->idx();

		Face face;
		face.reserve(3);
		MyMesh::ConstFaceVertexIter fv_it;
		MyMesh::ConstFaceVertexIter fv_end = _mesh.fv_end(*f_it);
		for (fv_it = _mesh.fv_begin(*f_it); fv_it != fv_end; ++fv_it)
		{
			VertexIndex v_idx = fv_it->idx();
			face.push_back(v_idx);
		}
		faces[f_idx] = std::move(face);
	}

	// Sort adjacent vertices and faces in the same order as the faces
	sortAdjacentVerticesAndFaces();
}
//=============================================================================
void BlackBox::computeShading(const MyMesh &_mesh,
	const vector<Intensity> &_sh_coeff, const unsigned int &_sh_order,
	vector<Intensity> &_shading)
{
	_shading.resize(_mesh.n_vertices());
	MyMesh::ConstVertexIter v_it;
	MyMesh::ConstVertexIter v_end(_mesh.vertices_end());
	for (v_it = _mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		int v_idx = v_it->idx();

		Normal normal(3);
		MyMesh::Normal n = vector_cast<MyMesh::Normal>(_mesh.normal(*v_it));
		normal[0] = n[0];
		normal[1] = n[1];
		normal[2] = n[2];

		Intensity shading = computeShading(normal, _sh_coeff, _sh_order);

		_shading[v_idx] = shading;
	}
}
//=============================================================================
Intensity BlackBox::computeShading(const Normal &_normal,
	const vector<Intensity> &_sh_coeff, const unsigned int &_sh_order)
{
	Coordinate x = _normal[0];	// x
	Coordinate y = _normal[1];	// y
	Coordinate z = _normal[2];	// z
	Coordinate x2 = x * x;		// x^2
	Coordinate y2 = y * y;		// y^2
	Coordinate z2 = z * z;		// z^2
	Coordinate xy = x * y;		// x * y
	Coordinate xz = x * z;		// x * z
	Coordinate yz = y * z;		// y * z
	Coordinate x2_y2 = x2 - y2;	// x^2 - y^2

	vector<Intensity> sh_functions;
	sh_functions.reserve(_sh_coeff.size());
	if (_sh_order >= 0)
	{
		sh_functions.push_back(1);
	}

	if (_sh_order > 0)
	{
		sh_functions.push_back(x);	// x
		sh_functions.push_back(y);	// y
		sh_functions.push_back(z);	// z
	}

	if (_sh_order > 1)
	{
		sh_functions.push_back(xy);			// x * y
		sh_functions.push_back(xz);			// x * z
		sh_functions.push_back(yz);			// y * z
		sh_functions.push_back(x2_y2);		// x^2 - y^2
		sh_functions.push_back(3 * z2 - 1);	// 3 * z^2 - 1
	}

	if (_sh_order > 2)
	{
		sh_functions.push_back((3 * x2 - y2) * y);	// (3 * x^2 - y^2) * y 
		sh_functions.push_back(xy * z);				// x * y * z
		sh_functions.push_back((5 * z2 - 1) * y);	// (5 * z^2 - 1) * y
		sh_functions.push_back((5 * z2 - 3) * z);	// (5 * z^2 - 3) * z
		sh_functions.push_back((5 * z2 - 1) * x);	// (5 * z^2 - 1) * x
		sh_functions.push_back(x2_y2 * z);			// (x^2 - y^2) * z
		sh_functions.push_back((x2 - 3 * y2) * x);	// (x^2 - 3 * y^2) * x
	}

	if (_sh_order > 3)
	{
		sh_functions.push_back(x2_y2 * xy);					// (x^2 - y^2) * x * y
		sh_functions.push_back((3 * x2 - y2) * yz);			// (3 * x^2 - y^2) * yz
		sh_functions.push_back((7 * z2 - 1) * xy);			// (7 * z^2 - 1) * x * y
		sh_functions.push_back((7 * z2 - 3) * yz);			// (7 * z^2 - 3) * y * z
		sh_functions.push_back(3 - 30 * z2 + 35 * z2 * z2);	// 3 - 30 * z^2 + 35 * z^4
		sh_functions.push_back((7 * z2 - 3) * xz);			// (7 * z^2 - 3) * x * z
		sh_functions.push_back((7 * z2 - 1) * x2_y2);		// (7 * z^2 - 1) * (x^2 - y^2)
		sh_functions.push_back((x2 - 3 * y2) * xz);			// (x^2 - 3 * y^2) * x * z
		sh_functions.push_back((x2 - 3 * y2) * x2			// (x^2 - 3 * y^2) * x^2 - (3 * x^2 - y^2) * y^2 
			- (3 * x2 - y2) * y2);
	}

	Intensity shading = 0;
	for (size_t i = 0; i < _sh_coeff.size(); i++)
	{
		shading += _sh_coeff[i] * sh_functions[i];
	}

	return max(shading, std::numeric_limits<Intensity>::epsilon());
}
//=============================================================================
void BlackBox::computeEstDiffuse(const vector<Color> &_albedos,
	const vector<Intensity> &_shadings, vector<Color> &_est_diffuse)
{
	_est_diffuse.reserve(_albedos.size());
	for (size_t i = 0; i < _albedos.size(); i++)
	{
		Color est_intensity(3);
		for (size_t j = 0; j < 3; j++)
		{
			est_intensity[j] = _albedos[i][j] * _shadings[i];
		}
		_est_diffuse.push_back(est_intensity);
	}
}
//=============================================================================
void BlackBox::computeIntensityEstDiffuseDiff(const vector<Color> &_intensity,
	const vector<Color> &_diffuse, vector<Color> &_diff)
{
	_diff.reserve(_intensity.size());
	for (size_t i = 0; i < _intensity.size(); i++)
	{
		Color diff_color(3);
		for (size_t j = 0; j < 3; j++)
		{
			diff_color[j] = _intensity[i][j] - _diffuse[i][j];
		}
		_diff.push_back(diff_color);
	}
}
//=============================================================================
void BlackBox::computeEstIntensity(const vector<Color> &_albedos,
	const vector<Intensity> &_shadings, const vector<Color> &_local_lightings,
	vector<Color> &_est_intensities)
{
	_est_intensities.reserve(_albedos.size());
	for (size_t i = 0; i < _albedos.size(); i++)
	{
		Color est_intensity(3);
		for (size_t j = 0; j < 3; j++)
		{
			est_intensity[j] = 
				_albedos[i][j] * _shadings[i] + _local_lightings[i][j];
		}
		_est_intensities.push_back(est_intensity);
	}
}
//=============================================================================
void BlackBox::setMeshVertices(const vector<Vertex> &_vertices, MyMesh &_mesh)
{
	MyMesh::ConstVertexIter v_it;
	MyMesh::ConstVertexIter v_end(_mesh.vertices_end());

	for (v_it = _mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		int i = v_it->idx();

		MyMesh::Point p(_vertices[i][0], _vertices[i][1], _vertices[i][2]);
		mesh.set_point(*v_it, p);
	}
}
//=============================================================================
void BlackBox::setMeshColors(const vector<Color> &_colors, MyMesh &_mesh)
{
	MyMesh::ConstVertexIter v_it;
	MyMesh::ConstVertexIter v_end(_mesh.vertices_end());
	for (v_it = _mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		int i = v_it->idx();

		MyMesh::Color c(
			(unsigned char)(_colors[i][0] * 255.0), 
			(unsigned char)(_colors[i][1] * 255.0),
			(unsigned char)(_colors[i][2] * 255.0)
			);
		mesh.set_color(*v_it, c);
	}
}
//=============================================================================
void BlackBox::setMeshColors(const vector<Intensity> &_intensities, MyMesh &_mesh)
{
	MyMesh::ConstVertexIter v_it;
	MyMesh::ConstVertexIter v_end(_mesh.vertices_end());

	for (v_it = _mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		int i = v_it->idx();

		MyMesh::Color c(
			(unsigned char)(_intensities[i] * 255.0), 
			(unsigned char)(_intensities[i] * 255.0), 
			(unsigned char)(_intensities[i] * 255.0)
			);

		if (_intensities[i] > 1)
		{
			c[0] = 255;
			c[1] = 0;
			c[2] = 0;
		}

		if (_intensities[i] < 0)
		{
			c[0] = 0;
			c[1] = 0;
			c[2] = 255;
		}

		mesh.set_color(*v_it, c);
	}
}
//=============================================================================
Intensity BlackBox::rgb2gray(const Color &_color)
{
	//Intensity gray = RGB2GRAY_R_FACTOR * _color[0] 
	//	+ RGB2GRAY_G_FACTOR * _color[1] 
	//	+ RGB2GRAY_B_FACTOR * _color[2];

	Color::const_iterator max_it = max_element(_color.begin(), _color.end());
	Intensity gray = *max_it;

	return gray;
}
//=============================================================================
void BlackBox::sortAdjacentVerticesAndFaces()
{
	// Sort faces and vertices clockwise or anti-clockwise depending on the 
	// definition of faces
	for (unsigned int i = 0; i < n_vertices; i++)
	{
		int n_adj_faces = (int)adj_faces[i].size();

		// List of local indexes of pair of adjacent vertices per face.
		// Stored as it was defined in the mesh (clockwise or anti-clockwise)
		vector< pair<VertexIndex, VertexIndex> > adjFacesVerticesInd(n_adj_faces);

		// For the current vertex, loop for each adjacent face
		for (int j = 0; j < n_adj_faces; j++)
		{
			int faceIdx = adj_faces[i][j];

			int fv_idx0 = faces[faceIdx][0];
			int fv_idx1 = faces[faceIdx][1];
			int fv_idx2 = faces[faceIdx][2];

			// Neighbours from current vertex in order
			pair<VertexIndex, VertexIndex> neigh_idxs;

			if (fv_idx0 == i)
			{
				neigh_idxs.first = fv_idx1;
				neigh_idxs.second = fv_idx2;
			}
			else
			{
				if (fv_idx1 == i)
				{
					neigh_idxs.first = fv_idx2;
					neigh_idxs.second = fv_idx0;
				}
				else // fv_idx2 == i
				{
					neigh_idxs.first = fv_idx0;
					neigh_idxs.second = fv_idx1;
				}
			}

			adjFacesVerticesInd[j] = neigh_idxs;
		}

		vector<VertexIndex> ordered_neigh_vertices;
		vector<FaceIndex> ordered_neigh_faces;
		VertexIndex next_idx = adjFacesVerticesInd[0].first;
		bool has_full_ring = true;
		for (int j = 0; j < adj_vertices[i].size(); j++)
		{
			// Add current neighbour vertex idx
			ordered_neigh_vertices.push_back(next_idx);

			// Find the neighbour edge that has the current vertex index as start
			std::vector< pair<VertexIndex, VertexIndex> >::iterator it =
				find_if(adjFacesVerticesInd.begin(),
				adjFacesVerticesInd.end(),
				[&next_idx](const pair<VertexIndex, VertexIndex>& element){
				return element.first == next_idx;
			}
			);

			// If it was found, add the face and update the next index to find
			if (it != adjFacesVerticesInd.end())
			{
				next_idx = (*it).second;

				int neigh_face_index = (int)(it - adjFacesVerticesInd.begin());
				ordered_neigh_faces.push_back(adj_faces[i][neigh_face_index]);
			}
			else // the current vertex has not a full ring and we stop 
			{
				has_full_ring = false;
				break;
			}
		}

		// If the vertex has not complete one-ring neighbours, there could be 
		// vertices that we have not added yet in the opposite direction to the
		// way that the faces are defined
		if (!has_full_ring)
		{
			// Number of vertex that have to be added yet 
			int n_neigh_left = (int)adj_vertices[i].size()
				- (int)ordered_neigh_vertices.size();
			VertexIndex next_idx = ordered_neigh_vertices[0];
			vector<VertexIndex> ordered_neigh_vertices_left;
			vector<VertexIndex> ordered_neigh_faces_left;
			for (int j = 0; j < n_neigh_left; ++j)
			{
				std::vector< pair<VertexIndex, VertexIndex> >::iterator it =
					find_if(adjFacesVerticesInd.begin(),
					adjFacesVerticesInd.end(),
					[&next_idx](const pair<VertexIndex, VertexIndex>& element){
					return element.second == next_idx;
				}
				);

				// If it was found, add the face and update the next index to find
				if (it != adjFacesVerticesInd.end())
				{
					next_idx = (*it).first;
					ordered_neigh_vertices_left.push_back(next_idx);

					int neigh_face_index = (int)(it - adjFacesVerticesInd.begin());
					ordered_neigh_faces_left.push_back(adj_faces[i][neigh_face_index]);
				}
				else // there are holes in the neighbour faces of the current vertex
				{
					cout << "Vertex " << i << " has holes in the neighbour faces" << endl;
					exit(0);
				}
			}

			// Reverse vertices left and add them to the beginning of the list
			reverse(ordered_neigh_vertices_left.begin(),
				ordered_neigh_vertices_left.end());
			ordered_neigh_vertices.insert(ordered_neigh_vertices.begin(),
				ordered_neigh_vertices_left.begin(), ordered_neigh_vertices_left.end());

			// Reverse faces left and add them to the beginning of the list
			reverse(ordered_neigh_faces_left.begin(),
				ordered_neigh_faces_left.end());
			ordered_neigh_faces.insert(ordered_neigh_faces.begin(),
				ordered_neigh_faces_left.begin(), ordered_neigh_faces_left.end());
		}

		// Overwrite ordered faces and vertices
		adj_vertices[i].assign(
			ordered_neigh_vertices.begin(), ordered_neigh_vertices.end());
		adj_faces[i].assign(
			ordered_neigh_faces.begin(), ordered_neigh_faces.end());
	}
}
//=============================================================================
void BlackBox::initDiffWeights()
{
	diff_weights.reserve(n_vertices);
	for (size_t i = 0; i < n_vertices; i++)
	{
		vector<Intensity> weights;
		weights.reserve(adj_vertices[i].size());
		for (size_t j = 0; j < adj_vertices[i].size(); j++)
		{
			unsigned int adj_v_idx = adj_vertices[i][j];
			Intensity weight = computeDiffWeight((unsigned int)i, adj_v_idx);
			weights.push_back(weight);
		}
		diff_weights.push_back(std::move(weights));
	}
}
//=============================================================================
void BlackBox::initLightingWeights()
{
	lighting_weights.resize(n_vertices, 1);
	if (params.specular_weight_var > 0)
	{
		for (size_t i = 0; i < n_vertices; i++)
		{
			Intensity weight = 
				exp(- local_lightings[i][0] / (2 * params.specular_weight_var));
			lighting_weights[i] = weight;
		}
	}
}
//=============================================================================
Intensity BlackBox::computeDiffWeight(const unsigned int _v_idx,
	const unsigned int _adj_v_idx)
{
	Intensity weight = 0.0;

	weight = 10 * max(local_lightings[_v_idx][0], local_lightings[_adj_v_idx][0]);

	//weight += 1 * (1 - min(shadings[_v_idx], shadings[_adj_v_idx]));

	//weight = 10 * exp(
	//	(max(local_lightings[_v_idx][0], local_lightings[_adj_v_idx][0]) - 1) 
	//	/ (2 * params.color_diff_var)
	//	);

	if (params.color_diff_var > 0)
	{
		Color &color1 = colors[_v_idx];
		Color &color2 = colors[_adj_v_idx];
		Color color_diff = {
			color1[0] - color2[0],
			color1[1] - color2[1],
			color1[2] - color2[2] 
		};

		Intensity color_diff_norm2 = pow(color_diff[0], 2) 
			+ pow(color_diff[1], 2)
			+ pow(color_diff[2], 2);

		if (color_diff_norm2 <= params.color_diff_threshold)
		{
			weight += exp(- color_diff_norm2 / (2 * params.color_diff_var));
		}
		//else
		//{
		//	weight = 0;
		//}
	}

	if (params.depth_diff_var > 0)
	{
		if (params.use_depth_weight)
		{
			Vertex &vertex1 = vertices[_v_idx];
			Vertex &vertex2 = vertices[_adj_v_idx];

			Coordinate depth_diff = vertex1[2] - vertex2[2];
			weight *= exp(- pow(depth_diff, 2) / (2 * params.depth_diff_var));
		}
		else
		{
			Normal &normal1 = normals[_v_idx];
			Normal &normal2 = normals[_adj_v_idx];

			Coordinate normal_diff_cos = normal1[0] * normal2[0]
				+ normal1[1] * normal2[1] 
				+ normal1[2] * normal2[2];
			weight *= exp(- (1 - normal_diff_cos) / (2 * params.depth_diff_var));
		}
	}

	return weight;
}
//=============================================================================
void BlackBox::initVisibility(MyMesh &_mesh,
	const string _image_filename, const string _intrinsics_filename)
{
	cv::Mat color_image_uchar = cv::imread(_image_filename.c_str(),
		CV_LOAD_IMAGE_COLOR);

	// Convert color image to double
	color_image_uchar.convertTo(color_image,
		CV_64FC3, 1. / 255);

	readMatrix3d(K, _intrinsics_filename);

	computeVisibility(color_image_uchar.cols, color_image_uchar.rows,
		K, _mesh, visibility);
}
//=============================================================================
void BlackBox::estimateSHCoeff(const ceres::Solver::Options &_options,
	bool _use_gray_color)
{
	Color black = { 0.0, 0.0, 0.0 };

	vector<Color> intensities = colors;

	if (_use_gray_color)
	{
		for (size_t v_idx = 0; v_idx < n_vertices; v_idx++)
		{
			intensities[v_idx][0] = grays[v_idx];
			intensities[v_idx][1] = grays[v_idx];
			intensities[v_idx][2] = grays[v_idx];
		}
	}

	ceres::Problem problem;

	const unsigned int N_CHANNELS = 3;

	for (size_t i = 0; i < n_vertices; i++)
	{
		//if (visibility[i])
		//{
			ceres::LossFunction *huber_data_loss = NULL;
			double huber_width = params.data_huber_width[0];
			if (huber_width > 0.0)
			{
				huber_data_loss = new ceres::HuberLoss(huber_width);
			}

			Intensity weight = params.use_specular_weight ? lighting_weights[i] : 1;

			ceres::ScaledLoss* scaled_data_loss = new ceres::ScaledLoss(
				huber_data_loss,
				weight,
				ceres::TAKE_OWNERSHIP);

			ResidualPhotometricErrorWithNormal *residual =
				new ResidualPhotometricErrorWithNormal(
				&intensities[i][0], &normals[i][0], N_CHANNELS, sh_order,
				params.use_lower_bound_shading,
				params.use_upper_bound_shading);

			ceres::DynamicAutoDiffCostFunction<ResidualPhotometricErrorWithNormal, 5>* dyn_cost_function
				= new ceres::DynamicAutoDiffCostFunction< ResidualPhotometricErrorWithNormal, 5 >(residual);

			// List of pointers to translations per vertex
			vector<Intensity*> v_parameter_blocks;

			// Albedo
			dyn_cost_function->AddParameterBlock(N_CHANNELS);
			v_parameter_blocks.push_back(&albedos[i][0]);
			// SH Coeff
			dyn_cost_function->AddParameterBlock(n_sh_basis);
			v_parameter_blocks.push_back(&sh_coeff[0]);
			// Local lighting variations
			dyn_cost_function->AddParameterBlock(N_CHANNELS);
			v_parameter_blocks.push_back(&black[0]);

			dyn_cost_function->SetNumResiduals(N_CHANNELS);

			ceres::ResidualBlockId residualBlockId = problem.AddResidualBlock(
				dyn_cost_function,
				scaled_data_loss,
				v_parameter_blocks);
		//}
	}

	problem.SetParameterBlockConstant(&black[0]);
	for (size_t i = 0; i < n_vertices; i++)
	{
		//if (visibility[i])
		//{
			problem.SetParameterBlockConstant(&albedos[i][0]);
		//}
	}

	ceres::Solver::Summary summary;

	ceres::Solve(_options, &problem, &summary);

	cout << summary.FullReport() << endl;
}
//=============================================================================
void BlackBox::estimateAlbedo(const ceres::Solver::Options &_options)
{
	double data_weight = params.data_weights[1];
	double huber_width = params.data_huber_width[1];

	Color black = { 0, 0, 0 };

	ceres::Problem problem;
	for (size_t i = 0; i < n_vertices; i++)
	{
		if (shadings[i] >= 0 /*&& visibility[i]*/)
		{
			ceres::LossFunction *huber_data_loss = NULL;
			if (huber_width > 0.0)
			{
				huber_data_loss = new ceres::HuberLoss(huber_width);
			}

			Intensity shading_weight = shadings[i];
			if (params.shading_weight_var > 0)
			{
				shading_weight = exp((abs(shadings[i]) - 1) / (2 * params.shading_weight_var));
			}

			Intensity weight = params.use_specular_weight ? 
				min(shadings[i], lighting_weights[i]) : shading_weight;

			ceres::ScaledLoss* lighting_scaled_loss = new ceres::ScaledLoss(
				huber_data_loss,
				weight,
				ceres::TAKE_OWNERSHIP);

			ceres::ScaledLoss* scaled_data_loss = new ceres::ScaledLoss(
				lighting_scaled_loss,
				data_weight,
				ceres::TAKE_OWNERSHIP);

			ResidualPhotometricErrorWithShading *residual = 
				new ResidualPhotometricErrorWithShading(
				&colors[i][0], shadings[i], 3);

			ceres::DynamicAutoDiffCostFunction<ResidualPhotometricErrorWithShading, 5>* dyn_cost_function
				= new ceres::DynamicAutoDiffCostFunction< ResidualPhotometricErrorWithShading, 5 >(residual);

			// List of pointers to translations per vertex
			vector<Intensity*> v_parameter_blocks;

			// Albedo
			dyn_cost_function->AddParameterBlock(3);
			v_parameter_blocks.push_back(&albedos[i][0]);
			// No local lighting variations
			dyn_cost_function->AddParameterBlock(3);
			v_parameter_blocks.push_back(&black[0]);

			dyn_cost_function->SetNumResiduals(3);

			ceres::ResidualBlockId residualBlockId = problem.AddResidualBlock(
				dyn_cost_function,
				scaled_data_loss,
				v_parameter_blocks);
		}
	}

	if (params.smooth_albedo_weight > 0.0)
	{
		ceres::LossFunction* huber_loss = NULL;

		if (params.smooth_albedo_huber_width > 0)
		{
			huber_loss = new ceres::HuberLoss(params.smooth_albedo_huber_width);
		}

		ceres::ScaledLoss* loss_function = new ceres::ScaledLoss(
			huber_loss, 
			params.smooth_albedo_weight, 
			ceres::TAKE_OWNERSHIP);

		for (size_t i = 0; i < n_vertices; i++)
		{
			if (shadings[i] > 0 /*&& visibility[i]*/)
			{
				for (size_t j = 0; j < adj_vertices[i].size(); j++)
				{
					if (shadings[j] > 0 /*&& visibility[j]*/)
					{
						unsigned int adj_v_idx = adj_vertices[i][j];

						Intensity weight = diff_weights[i][j];

						if (weight > 0.0)
						{
							ResidualWeightedDifference *residual =
								new ResidualWeightedDifference(weight, 3);

							ceres::AutoDiffCostFunction<ResidualWeightedDifference, 3, 3, 3>* cost_function =
								new ceres::AutoDiffCostFunction<ResidualWeightedDifference, 3, 3, 3>(residual);

							ceres::ResidualBlockId residualBlockId =
								problem.AddResidualBlock(
								cost_function,
								loss_function,
								&albedos[i][0],
								&albedos[adj_v_idx][0]);
						}
					}
				}
			}
		}
	}

	problem.SetParameterBlockConstant(&black[0]);
	for (size_t i = 0; i < n_vertices; ++i)
	{
		if (shadings[i] > 0 /*&& visibility[i]*/)
		{
			if (params.use_lower_bound_albedo)
			{
				problem.SetParameterLowerBound(&albedos[i][0], 0, 0.0);
				problem.SetParameterLowerBound(&albedos[i][0], 1, 0.0);
				problem.SetParameterLowerBound(&albedos[i][0], 2, 0.0);
			}

			if (params.use_upper_bound_albedo)
			{
				problem.SetParameterUpperBound(&albedos[i][0], 0, 1.0);
				problem.SetParameterUpperBound(&albedos[i][0], 1, 1.0);
				problem.SetParameterUpperBound(&albedos[i][0], 2, 1.0);
			}
		}
	}

	ceres::Solver::Summary summary;

	ceres::Solve(_options, &problem, &summary);

	cout << summary.FullReport() << endl;
}
//=============================================================================
void BlackBox::estimateLocalLighting(const ceres::Solver::Options &_options)
{
	ceres::Problem problem;

	ceres::LossFunction *data_loss_function = NULL;
	double huber_width = params.data_huber_width[2];
	if (huber_width > 0.0)
	{
		data_loss_function = new ceres::HuberLoss(huber_width);
	}

	for (size_t i = 0; i < n_vertices; i++)
	{
		if (true /*visibility[i]*/)
		{
			ResidualPhotometricErrorWithAlbedoShading *residual =
				new ResidualPhotometricErrorWithAlbedoShading(
				&colors[i][0], &albedos[i][0], shadings[i], 3);

			ceres::DynamicAutoDiffCostFunction<ResidualPhotometricErrorWithAlbedoShading, 5>* dyn_cost_function
				= new ceres::DynamicAutoDiffCostFunction< ResidualPhotometricErrorWithAlbedoShading, 5 >(residual);

			// List of pointers to translations per vertex
			vector<Intensity*> v_parameter_blocks;

			// Local lighting variations
			dyn_cost_function->AddParameterBlock(3);
			v_parameter_blocks.push_back(&local_lightings[i][0]);

			dyn_cost_function->SetNumResiduals(3);

			ceres::ResidualBlockId residualBlockId = problem.AddResidualBlock(
				dyn_cost_function,
				data_loss_function,
				v_parameter_blocks);
		}
	}

	if (params.smooth_local_lighting_weight > 0.0)
	{
		ceres::LossFunction* huber_loss = NULL;

		if (params.smooth_local_lighting_huber_width > 0)
		{
			huber_loss = new ceres::HuberLoss(params.smooth_local_lighting_huber_width);
		}

		ceres::ScaledLoss* loss_function = new ceres::ScaledLoss(
			huber_loss, 
			params.smooth_local_lighting_weight, 
			ceres::TAKE_OWNERSHIP);

		for (size_t i = 0; i < n_vertices; i++)
		{
			for (size_t j = 0; j < adj_vertices[i].size(); j++)
			{
				unsigned int adj_v_idx = adj_vertices[i][j];

				Intensity weight = 1;// diff_weights[i][j];

				if (weight > 0.0)
				{
					ResidualWeightedDifference *residual = new ResidualWeightedDifference(
						weight, 3);

					ceres::AutoDiffCostFunction<ResidualWeightedDifference, 3, 3, 3>* cost_function =
						new ceres::AutoDiffCostFunction<ResidualWeightedDifference, 3, 3, 3>(residual);

					ceres::ResidualBlockId residualBlockId = problem.AddResidualBlock(
						cost_function,
						loss_function,
						&local_lightings[i][0],
						&local_lightings[adj_v_idx][0]);
				}
			}
		}
	}

	if (params.local_lighting_magnitude_weight > 0.0)
	{
		ceres::LossFunction* huber_loss = NULL;

		if (params.local_lighting_magnitude_huber_width > 0)
		{
			huber_loss = new ceres::HuberLoss(params.local_lighting_magnitude_huber_width);
		}

		ceres::ScaledLoss* loss_function = new ceres::ScaledLoss(
			huber_loss,
			params.local_lighting_magnitude_weight, 
			ceres::TAKE_OWNERSHIP);

		for (size_t i = 0; i < n_vertices; i++)
		{
			ResidualValueMagnitude *residual = new ResidualValueMagnitude(3);

			ceres::AutoDiffCostFunction<ResidualValueMagnitude, 3, 3>* cost_function =
				new ceres::AutoDiffCostFunction<ResidualValueMagnitude, 3, 3>(residual);

			ceres::ResidualBlockId residualBlockId = problem.AddResidualBlock(
				cost_function,
				loss_function,
				&local_lightings[i][0]);
		}
	}

	for (size_t i = 0; i < n_vertices; ++i)
	{
		if (params.use_lower_bound_lighting)
		{
			problem.SetParameterLowerBound(&local_lightings[i][0], 0, 0.0);
			problem.SetParameterLowerBound(&local_lightings[i][0], 1, 0.0);
			problem.SetParameterLowerBound(&local_lightings[i][0], 2, 0.0);
		}

		if (params.use_upper_bound_lighting)
		{
			problem.SetParameterUpperBound(&local_lightings[i][0], 0, 1.0);
			problem.SetParameterUpperBound(&local_lightings[i][0], 1, 1.0);
			problem.SetParameterUpperBound(&local_lightings[i][0], 2, 1.0);
		}
	}

	ceres::Solver::Summary summary;

	ceres::Solve(_options, &problem, &summary);

	cout << summary.FullReport() << endl;
}
//=============================================================================
void BlackBox::estimateAlbedoLocalLighting(const ceres::Solver::Options &_options)
{
	ceres::Problem problem;

	ceres::LossFunction *huber_data_loss = NULL;
	double huber_width = params.data_huber_width[2];
	if (huber_width > 0.0)
	{
		huber_data_loss = new ceres::HuberLoss(huber_width);
	}

	for (size_t i = 0; i < n_vertices; i++)
	{
		ceres::ScaledLoss* scaled_data_loss = new ceres::ScaledLoss(
			huber_data_loss,
			shadings[i],
			ceres::TAKE_OWNERSHIP);

		ResidualPhotometricError *residual = new ResidualPhotometricError(
			&colors[i][0], 3, (unsigned int)adj_vertices[i].size(),
			(unsigned int)adj_faces[i].size(), sh_order);

		ceres::DynamicAutoDiffCostFunction<ResidualPhotometricError, 5>* dyn_cost_function
			= new ceres::DynamicAutoDiffCostFunction< ResidualPhotometricError, 5 >(residual);

		// List of pointers to translations per vertex
		vector<Intensity*> v_parameter_blocks;

		// Albedo
		dyn_cost_function->AddParameterBlock(3);
		v_parameter_blocks.push_back(&albedos[i][0]);
		// SH Coeff
		dyn_cost_function->AddParameterBlock(n_sh_basis);
		v_parameter_blocks.push_back(&sh_coeff[0]);
		// Local lighting variations
		dyn_cost_function->AddParameterBlock(3);
		v_parameter_blocks.push_back(&local_lightings[i][0]);
		// Vertex
		dyn_cost_function->AddParameterBlock(3);
		v_parameter_blocks.push_back(&vertices[i][0]);
		// Adjacent vertices
		for (size_t j = 0; j < adj_vertices[i].size(); j++)
		{
			int v_idx = adj_vertices[i][j];
			dyn_cost_function->AddParameterBlock(3);
			v_parameter_blocks.push_back(&vertices[v_idx][0]);
		}

		dyn_cost_function->SetNumResiduals(3);

		ceres::ResidualBlockId residualBlockId = problem.AddResidualBlock(
			dyn_cost_function,
			scaled_data_loss,
			v_parameter_blocks);
	}

	if (params.smooth_albedo_weight > 0.0)
	{
		ceres::ScaledLoss* loss_function = new ceres::ScaledLoss(
			NULL, params.smooth_albedo_weight, ceres::TAKE_OWNERSHIP);
		for (size_t i = 0; i < n_vertices; i++)
		{
			for (size_t j = 0; j < adj_vertices[i].size(); j++)
			{
				unsigned int adj_v_idx = adj_vertices[i][j];

				Intensity weight = diff_weights[i][j];

				if (weight > 0.0)
				{
					ResidualWeightedDifference *residual = new ResidualWeightedDifference(
						weight, 3);

					ceres::AutoDiffCostFunction<ResidualWeightedDifference, 3, 3, 3>* cost_function =
						new ceres::AutoDiffCostFunction<ResidualWeightedDifference, 3, 3, 3>(residual);

					ceres::ResidualBlockId residualBlockId = problem.AddResidualBlock(
						cost_function,
						loss_function,
						&albedos[i][0],
						&albedos[adj_v_idx][0]);
				}
			}
		}
	}

	if (params.smooth_local_lighting_weight > 0.0)
	{
		ceres::ScaledLoss* loss_function = new ceres::ScaledLoss(
			NULL, params.smooth_local_lighting_weight, ceres::TAKE_OWNERSHIP);
		for (size_t i = 0; i < n_vertices; i++)
		{
			for (size_t j = 0; j < adj_vertices[i].size(); j++)
			{
				unsigned int adj_v_idx = adj_vertices[i][j];

				Intensity weight = diff_weights[i][j];

				if (weight > 0.0)
				{
					ResidualWeightedDifference *residual = new ResidualWeightedDifference(
						weight, 3);

					ceres::AutoDiffCostFunction<ResidualWeightedDifference, 3, 3, 3>* cost_function =
						new ceres::AutoDiffCostFunction<ResidualWeightedDifference, 3, 3, 3>(residual);

					ceres::ResidualBlockId residualBlockId = problem.AddResidualBlock(
						cost_function,
						loss_function,
						&local_lightings[i][0],
						&local_lightings[adj_v_idx][0]);
				}
			}
		}
	}

	if (params.local_lighting_magnitude_weight > 0.0)
	{
		ceres::ScaledLoss* loss_function = new ceres::ScaledLoss(
			NULL, params.local_lighting_magnitude_weight, ceres::TAKE_OWNERSHIP);
		for (size_t i = 0; i < n_vertices; i++)
		{
			ResidualValueMagnitude *residual = new ResidualValueMagnitude(3);

			ceres::AutoDiffCostFunction<ResidualValueMagnitude, 3, 3>* cost_function =
				new ceres::AutoDiffCostFunction<ResidualValueMagnitude, 3, 3>(residual);

			ceres::ResidualBlockId residualBlockId = problem.AddResidualBlock(
				cost_function,
				loss_function,
				&local_lightings[i][0]);
		}
	}

	for (size_t i = 0; i < n_vertices; ++i)
	{
		problem.SetParameterBlockConstant(&vertices[i][0]);

		if (params.use_lower_bound_albedo)
		{
			problem.SetParameterLowerBound(&albedos[i][0], 0, 0.0);
			problem.SetParameterLowerBound(&albedos[i][0], 1, 0.0);
			problem.SetParameterLowerBound(&albedos[i][0], 2, 0.0);
		}

		if (params.use_upper_bound_albedo)
		{
			problem.SetParameterUpperBound(&albedos[i][0], 0, 1.0);
			problem.SetParameterUpperBound(&albedos[i][0], 1, 1.0);
			problem.SetParameterUpperBound(&albedos[i][0], 2, 1.0);
		}

		if (params.use_lower_bound_lighting)
		{
			problem.SetParameterLowerBound(&local_lightings[i][0], 0, 0.0);
			problem.SetParameterLowerBound(&local_lightings[i][0], 1, 0.0);
			problem.SetParameterLowerBound(&local_lightings[i][0], 2, 0.0);
		}

		if (params.use_upper_bound_lighting)
		{
			problem.SetParameterUpperBound(&local_lightings[i][0], 0, 1.0);
			problem.SetParameterUpperBound(&local_lightings[i][0], 1, 1.0);
			problem.SetParameterUpperBound(&local_lightings[i][0], 2, 1.0);
		}
	}

	ceres::Solver::Summary summary;

	ceres::Solve(_options, &problem, &summary);

	cout << summary.FullReport() << endl;
}
//=============================================================================
void BlackBox::estimateShape(const ceres::Solver::Options &_options)
{
	ceres::Problem problem;

	vector<Vertex> initial_vertices = vertices;

	vector<double> vertex_disp(n_vertices, 0.0);

	ceres::LossFunction *data_loss_function = NULL;
	double huber_width = params.data_huber_width[3];
	if (huber_width > 0.0)
	{
		data_loss_function = new ceres::HuberLoss(huber_width);
	}

	for (size_t i = 0; i < n_vertices; i++)
	{
		ResidualPhotometricError *residual = new ResidualPhotometricError(
			&colors[i][0], 3, (unsigned int)adj_vertices[i].size(),
			(unsigned int)adj_faces[i].size(), sh_order, 
			params.use_normal_for_shape, false, false);

		ceres::DynamicAutoDiffCostFunction<ResidualPhotometricError, 5>* dyn_cost_function
			= new ceres::DynamicAutoDiffCostFunction< ResidualPhotometricError, 5 >(residual);

		// List of pointers to translations per vertex
		vector<Intensity*> v_parameter_blocks;

		// Albedo
		dyn_cost_function->AddParameterBlock(3);
		v_parameter_blocks.push_back(&albedos[i][0]);
		// SH Coeff
		dyn_cost_function->AddParameterBlock(n_sh_basis);
		v_parameter_blocks.push_back(&sh_coeff[0]);
		// Local lighting variations
		dyn_cost_function->AddParameterBlock(3);
		v_parameter_blocks.push_back(&local_lightings[i][0]);
		// Vertex
		dyn_cost_function->AddParameterBlock(3);
		v_parameter_blocks.push_back(&vertices[i][0]);
		// Adjacent vertices
		for (size_t j = 0; j < adj_vertices[i].size(); j++)
		{
			int v_idx = adj_vertices[i][j];
			dyn_cost_function->AddParameterBlock(3);
			v_parameter_blocks.push_back(&vertices[v_idx][0]);
		}

		if (params.use_normal_for_shape)
		{
			// Vertex normal
			dyn_cost_function->AddParameterBlock(3);
			v_parameter_blocks.push_back(&normals[i][0]);
			// Adjacent vertices normal
			for (size_t j = 0; j < adj_vertices[i].size(); j++)
			{
				int v_idx = adj_vertices[i][j];
				dyn_cost_function->AddParameterBlock(3);
				v_parameter_blocks.push_back(&normals[v_idx][0]);
			}

			// Vertex disp
			dyn_cost_function->AddParameterBlock(1);
			v_parameter_blocks.push_back(&vertex_disp[i]);
			// Adjacent vertices disp
			for (size_t j = 0; j < adj_vertices[i].size(); j++)
			{
				int v_idx = adj_vertices[i][j];
				dyn_cost_function->AddParameterBlock(1);
				v_parameter_blocks.push_back(&vertex_disp[v_idx]);
			}
		}

		dyn_cost_function->SetNumResiduals(3);

		ceres::ResidualBlockId residualBlockId = problem.AddResidualBlock(
			dyn_cost_function,
			data_loss_function,
			v_parameter_blocks);
	}

	if (params.displacement_weight > 0.0)
	{
		ceres::ScaledLoss* loss_function = new ceres::ScaledLoss(
			NULL, params.displacement_weight, ceres::TAKE_OWNERSHIP);
		for (size_t i = 0; i < n_vertices; i++)
		{
			ResidualInitialVertexDiff *residual = 
				new ResidualInitialVertexDiff(initial_vertices[i]);

			ceres::ResidualBlockId residualBlockId;
			if (params.use_normal_for_shape)
			{
				ceres::AutoDiffCostFunction<ResidualInitialVertexDiff, 1, 3, 1>* cost_function =
					new ceres::AutoDiffCostFunction<ResidualInitialVertexDiff, 1, 3, 1>(residual);

				residualBlockId = problem.AddResidualBlock(
					cost_function,
					loss_function,
					&vertices[i][0],
					&vertex_disp[i]
					);
			}
			else
			{
				ceres::AutoDiffCostFunction<ResidualInitialVertexDiff, 3, 3>* cost_function =
					new ceres::AutoDiffCostFunction<ResidualInitialVertexDiff, 3, 3>(residual);

				residualBlockId = problem.AddResidualBlock(
					cost_function,
					loss_function,
					&vertices[i][0]);
			}
		}
	}

	if (params.laplacian_smooth_weight > 0.0)
	{
		ceres::ScaledLoss* loss_function = new ceres::ScaledLoss(
			NULL, params.laplacian_smooth_weight, ceres::TAKE_OWNERSHIP);
		for (size_t i = 0; i < n_vertices; i++)
		{
			ResidualLaplacianSmoothing *residual 
				= new ResidualLaplacianSmoothing(
				(unsigned int)adj_vertices[i].size(), 
				params.use_normal_for_shape);

			ceres::DynamicAutoDiffCostFunction<ResidualLaplacianSmoothing, 5>* dyn_cost_function
				= new ceres::DynamicAutoDiffCostFunction< ResidualLaplacianSmoothing, 5 >(residual);

			// List of pointers to translations per vertex
			vector<Coordinate*> v_parameter_blocks;

			// Vertex
			dyn_cost_function->AddParameterBlock(3);
			v_parameter_blocks.push_back(&vertices[i][0]);
			// Adjacent vertices
			for (size_t j = 0; j < adj_vertices[i].size(); j++)
			{
				int v_idx = adj_vertices[i][j];
				dyn_cost_function->AddParameterBlock(3);
				v_parameter_blocks.push_back(&vertices[v_idx][0]);
			}

			if (params.use_normal_for_shape)
			{
				// Vertex normal
				dyn_cost_function->AddParameterBlock(3);
				v_parameter_blocks.push_back(&normals[i][0]);
				// Adjacent vertices normal
				for (size_t j = 0; j < adj_vertices[i].size(); j++)
				{
					int v_idx = adj_vertices[i][j];
					dyn_cost_function->AddParameterBlock(3);
					v_parameter_blocks.push_back(&normals[v_idx][0]);
				}

				// Vertex disp
				dyn_cost_function->AddParameterBlock(1);
				v_parameter_blocks.push_back(&vertex_disp[i]);
				// Adjacent vertices disp
				for (size_t j = 0; j < adj_vertices[i].size(); j++)
				{
					int v_idx = adj_vertices[i][j];
					dyn_cost_function->AddParameterBlock(1);
					v_parameter_blocks.push_back(&vertex_disp[v_idx]);
				}
			}

			dyn_cost_function->SetNumResiduals(1);

			ceres::ResidualBlockId residualBlockId = problem.AddResidualBlock(
				dyn_cost_function,
				loss_function,
				v_parameter_blocks);
		}
	}

	if (params.tv_weight > 0.0)
	{
		ceres::ScaledLoss* loss_function = new ceres::ScaledLoss(
			NULL, params.tv_weight, ceres::TAKE_OWNERSHIP);
		for (size_t i = 0; i < n_vertices; i++)
		{
			for (size_t j = 0; j < adj_vertices[i].size(); j++)
			{
				int v_idx = adj_vertices[i][j];

				ResidualTV *residual = new ResidualTV(
					initial_vertices[i], initial_vertices[v_idx]);

				ceres::ResidualBlockId residualBlockId;
				if (params.use_normal_for_shape)
				{
					ceres::AutoDiffCostFunction<ResidualTV, 3, 3, 3, 3, 3, 1, 1>* cost_function =
						new ceres::AutoDiffCostFunction<ResidualTV, 3, 3, 3, 3, 3, 1, 1>(residual);

					residualBlockId = problem.AddResidualBlock(
						cost_function,
						loss_function,
						&vertices[i][0],
						&vertices[v_idx][0],
						&normals[i][0],
						&normals[v_idx][0],
						&vertex_disp[i],
						&vertex_disp[v_idx]);
				}
				else
				{
					ceres::AutoDiffCostFunction<ResidualTV, 3, 3, 3>* cost_function =
						new ceres::AutoDiffCostFunction<ResidualTV, 3, 3, 3>(residual);

					residualBlockId = problem.AddResidualBlock(
						cost_function,
						loss_function,
						&vertices[i][0],
						&vertices[v_idx][0]);
				}
			}
		}
	}

	problem.SetParameterBlockConstant(&sh_coeff[0]);
	for (size_t i = 0; i < n_vertices; ++i)
	{
		problem.SetParameterBlockConstant(&albedos[i][0]);
		problem.SetParameterBlockConstant(&local_lightings[i][0]);

		if (params.use_normal_for_shape)
		{
			problem.SetParameterBlockConstant(&vertices[i][0]);
			problem.SetParameterBlockConstant(&normals[i][0]);
		}
	}

	ceres::Solver::Summary summary;

	ceres::Solve(_options, &problem, &summary);

	cout << summary.FullReport() << endl;

	updateVertices(vertex_disp);
}
//=============================================================================
void BlackBox::refine()
{
	//refineSHCoeff(params.ceres_solver.options[0]);
	estimateSHCoeff(params.ceres_solver.options[0],
		params.input_extra_mesh_filename.empty());

	updateShadings();
	//updateAlbedos();

	if (params.combine_albedo_lighting)
	{
		updateLocalLightings();

		estimateAlbedoLocalLighting(params.ceres_solver.options[1]);
	}
	else
	{
		//estimateAlbedo(params.ceres_solver.options[1]);

		updateLocalLightings();
		estimateLocalLighting(params.ceres_solver.options[2]);
	}
	
	if (params.estimate_shape)
	{
		estimateShape(params.ceres_solver.options[3]);
	}
}
//=============================================================================
void BlackBox::refineSHCoeff(const ceres::Solver::Options &_options)
{
	ceres::Problem problem;

	ceres::LossFunction *data_loss_function = NULL;
	double huber_width = params.data_huber_width[0];
	if (huber_width > 0.0)
	{
		data_loss_function = new ceres::HuberLoss(huber_width);
	}

	for (size_t i = 0; i < n_vertices; i++)
	{
		ResidualPhotometricError *residual = new ResidualPhotometricError(
			&colors[i][0], 3, (unsigned int)adj_vertices[i].size(),
			(unsigned int)adj_faces[i].size(), sh_order, false,
			params.use_lower_bound_shading,
			params.use_upper_bound_shading);

		ceres::DynamicAutoDiffCostFunction<ResidualPhotometricError, 5>* dyn_cost_function
			= new ceres::DynamicAutoDiffCostFunction< ResidualPhotometricError, 5 >(residual);

		// List of pointers to translations per vertex
		vector<Intensity*> v_parameter_blocks;

		dyn_cost_function->AddParameterBlock(3);
		v_parameter_blocks.push_back(&albedos[i][0]);
		// SH Coeff
		dyn_cost_function->AddParameterBlock(n_sh_basis);
		v_parameter_blocks.push_back(&sh_coeff[0]);
		// Local lighting variations
		dyn_cost_function->AddParameterBlock(3);
		v_parameter_blocks.push_back(&local_lightings[i][0]);
		// Vertex
		dyn_cost_function->AddParameterBlock(3);
		v_parameter_blocks.push_back(&vertices[i][0]);
		// Adjacent vertices
		for (size_t j = 0; j < adj_vertices[i].size(); j++)
		{
			int v_idx = adj_vertices[i][j];
			dyn_cost_function->AddParameterBlock(3);
			v_parameter_blocks.push_back(&vertices[v_idx][0]);
		}

		dyn_cost_function->SetNumResiduals(3);

		ceres::ResidualBlockId residualBlockId = problem.AddResidualBlock(
			dyn_cost_function,
			data_loss_function,
			v_parameter_blocks);
	}

	for (size_t i = 0; i < n_vertices; ++i)
	{
		problem.SetParameterBlockConstant(&vertices[i][0]);
		problem.SetParameterBlockConstant(&albedos[i][0]);
		problem.SetParameterBlockConstant(&local_lightings[i][0]);
	}

	ceres::Solver::Summary summary;

	ceres::Solve(_options, &problem, &summary);

	cout << summary.FullReport() << endl;
}
//=============================================================================
void BlackBox::updateVertices(vector<double> &_vertices_disp)
{
	for (size_t i = 0; i < _vertices_disp.size(); i++)
	{
		vertices[i][0] += _vertices_disp[i] * normals[i][0];
		vertices[i][1] += _vertices_disp[i] * normals[i][1];
		vertices[i][2] += _vertices_disp[i] * normals[i][2];
	}
}
//=============================================================================
void BlackBox::updateShadings()
{
	for (size_t i = 0; i < n_vertices; i++)
	{
		shadings[i] = computeShading(normals[i], sh_coeff, sh_order);
	}
}
//=============================================================================
void BlackBox::updateAlbedos()
{
	for (size_t i = 0; i < n_vertices; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			switch (params.albedo_initialization)
			{
			case WHITE_ALBEDO:
				albedos[i][j] = 1;
				break;
			case COLOR_ALBEDO:
				albedos[i][j] = colors[i][j];
				break;
			case EST_ALBEDO:
				albedos[i][j] = colors[i][j] 
					/ max(shadings[i], std::numeric_limits<double>::epsilon());
				break;
			default:
				break;
			}
		}
	}
}
//=============================================================================
void BlackBox::updateLocalLightings()
{
	for (size_t i = 0; i < n_vertices; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			local_lightings[i][j] = colors[i][j] - albedos[i][j] * shadings[i];
		}
	}
}
//=============================================================================
void BlackBox::updateColorsFromImage()
{
	MyMesh::ConstVertexIter v_it;
	MyMesh::ConstVertexIter v_end = mesh.vertices_end();
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		int v_idx = v_it->idx();

		if (visibility[v_idx])
		{
			MyMesh::Point point = mesh.point(*v_it);

			Eigen::Vector3d mesh_point(point[0], point[1], point[2]);

			Eigen::Vector3d img_point = K * mesh_point;
			img_point /= img_point(2);

			//Bilinear interpolation
			double x = img_point(0);
			double y = img_point(1);

			double x1 = floor(x);
			double y1 = floor(y);
			double x2 = ceil(x);
			double y2 = ceil(y);

			cv::Vec3d top_left = color_image.at<cv::Vec3d>((int)y1, (int)x1);
			cv::Vec3d top_right = color_image.at<cv::Vec3d>((int)y1, (int)x2);
			cv::Vec3d bottom_left = color_image.at<cv::Vec3d>((int)y2, (int)x1);
			cv::Vec3d bottom_right = color_image.at<cv::Vec3d>((int)y2, (int)x2);

			double dx = x - x1;
			double dy = y - y1;
			double dx_1 = 1 - dx;
			double dy_1 = 1 - dy;

			Intensity* color = &colors[v_idx][0];

			for (int i = 0; i < 3; i++)
			{

				color[i] = dx * dy * bottom_right[2-i]
					+ dx * dy_1 * top_right[2-i]
					+ dx_1 * dy * bottom_left[2-i]
					+ dx_1 * dy_1 * top_left[2-i];
			}
		}
	}
}
//=============================================================================
void BlackBox::initLocalLightings(MyMesh &_mesh, const string _image_filename, 
	const string _intrinsics_filename)
{
	cv::Mat color_image_uchar = cv::imread(_image_filename.c_str(), 
		CV_LOAD_IMAGE_COLOR);

	// Estimate brightness
	std::vector<cv::Mat> planes(3);
	cv::split(color_image, planes);
	cv::Mat brightnessImage = cv::Mat(
		cv::max(planes[2],
		cv::max(planes[1], planes[0]))
		);

	// Estimate diffuse image
	cv::Mat diffuse_image_uchar;
	estimateDiffuse(color_image_uchar, diffuse_image_uchar);

	// Estimate diffuse brightness
	std::vector<cv::Mat> diffuse_planes(3);
	cv::split(diffuse_image_uchar, diffuse_planes);
	cv::Mat diffuseBrightnessImage_uchar = cv::Mat(
		cv::max(diffuse_planes[2],
		cv::max(diffuse_planes[1], diffuse_planes[0]))
		);

	// Convert diffuse brightness to double
	cv::Mat diffuseBrightnessImage;
	diffuseBrightnessImage_uchar.convertTo(diffuseBrightnessImage,
		cv::DataType<Intensity>::type, 1. / 255);

	// Compute specular brightness as the difference between full brightness 
	// and diffuse brightness
	cv::Mat specularImage;
	cv::subtract(brightnessImage, diffuseBrightnessImage, specularImage);

	MyMesh::ConstVertexIter v_it;
	MyMesh::ConstVertexIter v_end = _mesh.vertices_end();
	for (v_it = _mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		int v_idx = v_it->idx();

		if (visibility[v_idx])
		{
			MyMesh::Point point = _mesh.point(*v_it);

			Eigen::Vector3d mesh_point(point[0], point[1], point[2]);

			Eigen::Vector3d img_point = K * mesh_point;
			img_point /= img_point(2);

			//Bilinear interpolation
			double x = img_point(0);
			double y = img_point(1);

			double x1 = floor(x);
			double y1 = floor(y);
			double x2 = ceil(x);
			double y2 = ceil(y);

			Intensity top_left = specularImage.at<Intensity>((int)y1, (int)x1);
			Intensity top_right = specularImage.at<Intensity>((int)y1, (int)x2);
			Intensity bottom_left = specularImage.at<Intensity>((int)y2, (int)x1);
			Intensity bottom_right = specularImage.at<Intensity>((int)y2, (int)x2);

			double dx = x - x1;
			double dy = y - y1;
			double dx_1 = 1 - dx;
			double dy_1 = 1 - dy;

			Intensity* color = &local_lightings[v_idx][0];

			for (int i = 0; i < 3; i++)
			{
				color[i] = dx * dy * bottom_right
					+ dx * dy_1 * top_right
					+ dx_1 * dy * bottom_left
					+ dx_1 * dy_1 * top_left;
			}
		}
	}
}
//=============================================================================
void BlackBox::writeToPLY(std::string& filename, MyMesh& meshData)
{
	PlyFile *ply;
	float version;

	/* open either a binary or ascii PLY file for writing */
	/* (the file will be called "test.ply" because the routines */
	/*  enforce the .ply filename extension) */

	bool save_normal = false;
	bool save_alpha = false;
	bool save_specular = true;

	int n_elems = save_specular ? 3 : 2;

	if (params.save_mesh_binary)
	{
		ply = ply_open_for_writing(filename.c_str(), 3, ply::elem_names, PLY_BINARY_LE, &version);
	}
	else
	{
		ply = ply_open_for_writing(filename.c_str(), 3, ply::elem_names, PLY_ASCII, &version);
	}

	/* describe what properties go into the vertex and face elements */

	int num_vertices = (int)meshData.n_vertices();
	ply_element_count(ply, ply::elem_names[0], num_vertices);
	int vertex_type = 0;
	if (save_normal)
	{
		vertex_type |= 0x01;
	}
	if (save_alpha)
	{
		vertex_type |= 0x02;
	}
	if (save_specular)
	{
		vertex_type |= 0x04;
	}

	for (int i = 0; i < ply::n_vprops[vertex_type]; i++)
	{
		PlyProperty* prop = ply::get_vertex_property(vertex_type, i);
		ply_describe_property(ply, ply::elem_names[0], prop);
	}

	int num_faces = (int)meshData.n_faces();
	ply_element_count(ply, ply::elem_names[1], num_faces);
	for (int i = 0; i < ply::n_fprops; i++)
	{
		ply_describe_property(ply, ply::elem_names[1], &ply::face_props[i]);
	}

	if (save_specular)
	{
		int num_sh_coefficients = (int)sh_coeff.size();
		ply_element_count(ply, ply::elem_names[2], num_sh_coefficients);
		for (int i = 0; i < ply::n_sh_coeff_props; i++)
		{
			ply_describe_property(ply, ply::elem_names[2], &ply::sh_coeff_props[i]);
		}
	}

	/* we have described exactly what we will put in the file, so */
	/* we are now done with the header info */
	ply_header_complete(ply);

	/* set up and write the vertex elements */
	ply_put_element_setup(ply, ply::elem_names[0]);
	

	for (int i = 0; i < num_vertices; i++){
		switch (vertex_type)
		{
		case PLY_VERTEX_RGB:
		{
			ply::VertexColor v;
			v.x = (float)vertices[i][0];
			v.y = (float)vertices[i][1];
			v.z = (float)vertices[i][2];
			v.r = (unsigned char)(albedos[i][0] * 255.f);
			v.g = (unsigned char)(albedos[i][1] * 255.f);
			v.b = (unsigned char)(albedos[i][2] * 255.f);
			ply_put_element(ply, (void *)&v);
		}
		break;
		case PLY_VERTEX_NORMAL_RGB:
		{
			ply::VertexNormalColor v;
			v.x = (float)vertices[i][0];
			v.y = (float)vertices[i][1];
			v.z = (float)vertices[i][2];
			v.nx = (float)normals[i][0];
			v.ny = (float)normals[i][1];
			v.nz = (float)normals[i][2];
			v.r = (unsigned char)(albedos[i][0] * 255.f);
			v.g = (unsigned char)(albedos[i][1] * 255.f);
			v.b = (unsigned char)(albedos[i][2] * 255.f);
			ply_put_element(ply, (void *)&v);
		}
		break;
		case PLY_VERTEX_RGBA:
		{
			ply::VertexColorAlpha v;
			v.x = (float)vertices[i][0];
			v.y = (float)vertices[i][1];
			v.z = (float)vertices[i][2];
			v.r = (unsigned char)(albedos[i][0] * 255.f);
			v.g = (unsigned char)(albedos[i][1] * 255.f);
			v.b = (unsigned char)(albedos[i][2] * 255.f);
			v.a = 255;
			ply_put_element(ply, (void *)&v);
		}
		break;
		case PLY_VERTEX_NORMAL_RGBA:
		{
			ply::VertexNormalColorAlpha v;
			v.x = (float)vertices[i][0];
			v.y = (float)vertices[i][1];
			v.z = (float)vertices[i][2];
			v.nx = (float)normals[i][0];
			v.ny = (float)normals[i][1];
			v.nz = (float)normals[i][2];
			v.r = (unsigned char)(albedos[i][0] * 255.f);
			v.g = (unsigned char)(albedos[i][1] * 255.f);
			v.b = (unsigned char)(albedos[i][2] * 255.f);
			v.a = 255;
			ply_put_element(ply, (void *)&v);
		}
		break;
		case PLY_VERTEX_RGB_SPECULAR:
		{
			ply::VertexColorSpecular v;
			v.x = (float)vertices[i][0];
			v.y = (float)vertices[i][1];
			v.z = (float)vertices[i][2];
			v.r = (unsigned char)(albedos[i][0] * 255.f);
			v.g = (unsigned char)(albedos[i][1] * 255.f);
			v.b = (unsigned char)(albedos[i][2] * 255.f);
			v.specular_r = (unsigned char)(local_lightings[i][0] * 255.f);
			v.specular_g = (unsigned char)(local_lightings[i][1] * 255.f);
			v.specular_b = (unsigned char)(local_lightings[i][2] * 255.f);
			ply_put_element(ply, (void *)&v);
		}
		break;
		case PLY_VERTEX_NORMAL_RGB_SPECULAR:
		{
			ply::VertexNormalColorAlphaSpecular v;
			v.x = (float)vertices[i][0];
			v.y = (float)vertices[i][1];
			v.z = (float)vertices[i][2];
			v.nx = (float)normals[i][0];
			v.ny = (float)normals[i][1];
			v.nz = (float)normals[i][2];
			v.r = (unsigned char)(albedos[i][0] * 255.f);
			v.g = (unsigned char)(albedos[i][1] * 255.f);
			v.b = (unsigned char)(albedos[i][2] * 255.f);
			v.specular_r = (unsigned char)(local_lightings[i][0] * 255.f);
			v.specular_g = (unsigned char)(local_lightings[i][1] * 255.f);
			v.specular_b = (unsigned char)(local_lightings[i][2] * 255.f);
			ply_put_element(ply, (void *)&v);
		}
		break;
		case PLY_VERTEX_RGBA_SPECULAR:
		{
			ply::VertexColorAlphaSpecular v;
			v.x = (float)vertices[i][0];
			v.y = (float)vertices[i][1];
			v.z = (float)vertices[i][2];
			v.r = (unsigned char)(albedos[i][0] * 255.f);
			v.g = (unsigned char)(albedos[i][1] * 255.f);
			v.b = (unsigned char)(albedos[i][2] * 255.f);
			v.a = 255;
			v.specular_r = (unsigned char)(local_lightings[i][0] * 255.f);
			v.specular_g = (unsigned char)(local_lightings[i][1] * 255.f);
			v.specular_b = (unsigned char)(local_lightings[i][2] * 255.f);
			ply_put_element(ply, (void *)&v);
		}
		break;
		case PLY_VERTEX_NORMAL_RGBA_SPECULAR:
		{
			ply::VertexNormalColorAlphaSpecular v;
			v.x = (float)vertices[i][0];
			v.y = (float)vertices[i][1];
			v.z = (float)vertices[i][2];
			v.nx = (float)normals[i][0];
			v.ny = (float)normals[i][1];
			v.nz = (float)normals[i][2];
			v.r = (unsigned char)(albedos[i][0] * 255.f);
			v.g = (unsigned char)(albedos[i][1] * 255.f);
			v.b = (unsigned char)(albedos[i][2] * 255.f);
			v.a = 255;
			v.specular_r = (unsigned char)(local_lightings[i][0] * 255.f);
			v.specular_g = (unsigned char)(local_lightings[i][1] * 255.f);
			v.specular_b = (unsigned char)(local_lightings[i][2] * 255.f);
			ply_put_element(ply, (void *)&v);
		}
		break;
		default:
			break;
		}
	}

	/* set up and write the face elements */
	ply_put_element_setup(ply, ply::elem_names[1]);
	ply::Face f;
	f.verts = new int[10];

	MyMesh::ConstFaceIter f_it;
	MyMesh::ConstFaceIter f_end = meshData.faces_end();

	for (f_it = meshData.faces_begin(); f_it != f_end; ++f_it)
	{
		MyMesh::Face face = meshData.face(*f_it);
		f.nverts = (unsigned char) 3;

		MyMesh::ConstFaceVertexIter fv_it;
		MyMesh::ConstFaceVertexIter fv_end = meshData.fv_end(*f_it);
		int j = 0;
		for (fv_it = meshData.fv_begin(*f_it); fv_it != fv_end; ++fv_it)
		{
			int v_idx = fv_it->idx();
			f.verts[j] = v_idx;
			j++;
		}
		ply_put_element(ply, (void *)&f);
	}

	/* set up and write the face elements */
	ply_put_element_setup(ply, ply::elem_names[2]);
	for (int i = 0; i < sh_coeff.size(); i++)
	{
		ply::SH_Coefficient sh_c;
		sh_c.value = (float)sh_coeff[i];
		ply_put_element(ply, (void *)&sh_c);
	}

	/* close the PLY file */
	ply_close(ply);
}

void BlackBox::initAlbedos()
{
	sh_coeff.resize(n_sh_basis, 0);

	Color white = {1, 1, 1};
	albedos.resize(n_vertices, white);

	if (params.input_extra_mesh_filename.empty())
	{
		if (params.init_kmeans)
		{
			initAlbedoKmeans(params.kmeans_clusters);
		}
		else // Init using percentile
		{
			Intensity white = 1.0;

			if (params.albedo_percentile > 0.0)
			{
				vector<Intensity> aux_grays = grays;
				unsigned int nth = (unsigned int)(params.albedo_percentile * (double)aux_grays.size());
				nth_element(aux_grays.begin(), aux_grays.begin() + nth, aux_grays.end());
				white = aux_grays[nth - 1];
			}

			for (size_t v_idx = 0; v_idx < n_vertices; v_idx++)
			{
				albedos[v_idx][0] = white;
				albedos[v_idx][1] = white;
				albedos[v_idx][2] = white;
			}
		}
	}
	else
	{
		for (size_t v_idx = 0; v_idx < extra_mesh.n_vertices(); ++v_idx)
		{
			MeshData::ColorT color = extra_mesh.vertex_colors[v_idx];
			albedos[v_idx][0] = (Intensity)color[0];
			albedos[v_idx][1] = (Intensity)color[1];
			albedos[v_idx][2] = (Intensity)color[2];
		}

		sh_coeff.resize(extra_mesh.n_sh_coefficients());
		for (size_t i = 0; i < sh_coeff.size(); i++)
		{
			sh_coeff[i] = (Intensity)extra_mesh.sh_coefficients[i];
		}
	}
}

void BlackBox::estimateDiffuse(const cv::Mat src_image, cv::Mat &diffuse_image)
{
	qx_timer timer;

	int h = src_image.rows;
	int w = src_image.cols;
	int c = src_image.channels();

	uchar*** qx_image = qx_allocu_3(h, w, 3);
	uchar*** qx_diffuse_image = qx_allocu_3(h, w, 3);

	const uchar* p;
	for (size_t i = 0; i < h; i++)
	{
		p = src_image.ptr<uchar>(i);
		for (size_t j = 0; j < w; j++)
		{
			for (size_t k = 0; k < c; k++)
			{
				qx_image[i][j][k] = p[j * c + k];
			}
		}
	}

	qx_highlight_removal_bf m_highlight;
	m_highlight.init(h, w);/*initialization*/
	timer.start();

	int nr_iter_converge = m_highlight.diffuse(qx_diffuse_image, qx_image);/*extracting diffuse reflection*/

	timer.time_display("Highlight Removal");
	printf("# of iterations before convergence: [%02d]\n", nr_iter_converge);

	diffuse_image = cv::Mat(h, w, CV_8UC3);
	uchar* p2;
	for (size_t i = 0; i < h; i++)
	{
		p2 = diffuse_image.ptr<uchar>(i);
		for (size_t j = 0; j < w; j++)
		{
			for (size_t k = 0; k < c; k++)
			{
				p2[j * c + k] = qx_diffuse_image[i][j][k];
			}
		}
	}

	qx_freeu_3(qx_image);
	qx_image = NULL;

	qx_freeu_3(qx_diffuse_image);
	qx_diffuse_image = NULL;
}

void BlackBox::initAlbedoKmeans(const int _kmeans_clusters)
{
	int n_points = n_vertices;
	int n_channels = 3;
	int n_clusters = _kmeans_clusters;

	Mat cv_colors(n_points, 1, CV_32FC3);

	// K-means by chromaticity
	MatIterator_<cv::Vec3f> it;
	MatIterator_<cv::Vec3f> end = cv_colors.end<cv::Vec3f>();
	int i = 0;
	for (it = cv_colors.begin<cv::Vec3f>(); it != end; ++it)
	{
		float rgb = colors[i][0] + colors[i][1] + colors[i][2];

		if (rgb < std::numeric_limits<float>::epsilon())
		{
			rgb = std::numeric_limits<float>::epsilon();
		}

		(*it)[0] = (float)(colors[i][0] / rgb);
		(*it)[1] = (float)(colors[i][1] / rgb);
		(*it)[2] = (float)(colors[i][2] / rgb);

		i++;
	}

	Mat labels;
	Mat centers;
	TermCriteria crit = TermCriteria(CV_TERMCRIT_ITER + CV_TERMCRIT_EPS, 100, 0.5);

	kmeans(cv_colors, n_clusters, labels, crit, 3, KMEANS_PP_CENTERS, centers);

	// Compute intensity per cluster
	vector<Color> cluster_colors;
	cluster_colors.resize(_kmeans_clusters, Color({ 0, 0, 0 }));

	vector<vector<float>> red(_kmeans_clusters);
	vector<vector<float>> green(_kmeans_clusters);
	vector<vector<float>> blue(_kmeans_clusters);

	for (size_t i = 0; i < n_points; i++)
	{
		int cluster_idx = labels.at<int>(i);

		red[cluster_idx].push_back(colors[i][0]);
		green[cluster_idx].push_back(colors[i][1]);
		blue[cluster_idx].push_back(colors[i][2]);
	}

	for (size_t i = 0; i < _kmeans_clusters; i++)
	{
		float percentage = params.albedo_percentile;

		unsigned int percentile = (unsigned int)(red[i].size() * percentage);

		// Max colors
		std::nth_element(red[i].begin(), red[i].begin() + percentile, red[i].end());
		cluster_colors[i][0] = *(red[i].begin() + percentile);

		std::nth_element(green[i].begin(), green[i].begin() + percentile, green[i].end());
		cluster_colors[i][1] = *(green[i].begin() + percentile);

		std::nth_element(blue[i].begin(), blue[i].begin() + percentile, blue[i].end());
		cluster_colors[i][2] = *(blue[i].begin() + percentile);
	}

	for (size_t i = 0; i < n_points; i++)
	{	
		int cluster_idx = labels.at<int>(i);

		albedos[i][0] = (double)cluster_colors[cluster_idx][0];
		albedos[i][1] = (double)cluster_colors[cluster_idx][1];
		albedos[i][2] = (double)cluster_colors[cluster_idx][2];
	}
}