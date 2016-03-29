#include "BlackBox.h"
//=============================================================================
Parameters BlackBox::params;

MyMesh BlackBox::mesh;

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
vector< vector<Intensity> > BlackBox::diff_weights;

OpenMesh::IO::Options BlackBox::mesh_read_opt;
OpenMesh::IO::Options BlackBox::mesh_write_opt;
//=============================================================================
void BlackBox::initialize(int *argc, char **argv)
{
	params.load(argv[1]);

	mesh.request_vertex_colors();
	mesh.request_face_normals();
	mesh.request_vertex_normals();

	mesh_read_opt += OpenMesh::IO::Options::VertexColor;
	mesh_read_opt += OpenMesh::IO::Options::VertexNormal;

	readMesh(mesh, params.input_mesh_filename.c_str(), mesh_read_opt);

	mesh.update_normals();

	mesh_write_opt += OpenMesh::IO::Options::VertexColor;
	mesh_write_opt += OpenMesh::IO::Options::VertexNormal;

	if (params.save_mesh_binary)
	{
		mesh_write_opt += OpenMesh::IO::Options::Binary;
		mesh_write_opt += OpenMesh::IO::Options::LSB;
	}

	sh_order = params.sh_order;
	n_sh_basis = (unsigned int)pow(sh_order + 1, 2);

	initFromMesh(mesh);
}
//=============================================================================
void BlackBox::run()
{
	Color white = { 1.0, 1.0, 1.0 };
	Color black = { 0.0, 0.0, 0.0 };
	sh_coeff.resize(n_sh_basis, 0);
	albedos.resize(n_vertices, white);
	local_lightings.resize(n_vertices, black);

	estimateSHCoeff(params.ceres_solver.options[0]);
	estimateAlbedo(params.ceres_solver.options[1]);
	estimateLocalLighting(params.ceres_solver.options[2]);
	estimateShape(params.ceres_solver.options[3]);
}
//=============================================================================
void BlackBox::save()
{
	if (params.save_sh_coeff)
	{
		writeVector(sh_coeff, params.output_sh_coeff_filename.c_str());
	}

	setMeshVertices(vertices, mesh);
	mesh.update_normals();

	if (params.save_orig_intensity)
	{
		writeMesh(mesh, params.output_orig_intensity_mesh_filename, mesh_write_opt);
	}

	if (params.save_albedo)
	{
		setMeshColors(albedos, mesh);

		writeMesh(mesh, params.output_albedo_mesh_filename, mesh_write_opt);
	}

	if (params.save_local_lighting)
	{
		setMeshColors(local_lightings, mesh);

		writeMesh(mesh, params.output_local_lighting_mesh_filename, mesh_write_opt);
	}

	vector<Intensity> shadings;
	computeShading(mesh, sh_coeff, sh_order, shadings);

	if (params.save_shading)
	{
		setMeshColors(shadings, mesh);

		writeMesh(mesh, params.output_shading_mesh_filename, mesh_write_opt);
	}

	if (params.save_est_intensity)
	{
		vector<Color> est_intensities;
		computeEstIntensity(albedos, shadings, local_lightings, est_intensities);

		setMeshColors(est_intensities, mesh);

		writeMesh(mesh, params.output_est_intensity_mesh_filename, mesh_write_opt);
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
		color[0] = c[0];
		color[1] = c[1];
		color[2] = c[2];
		colors[v_idx] = std::move(color);

		Intensity gray = rgb2gray(color);
		grays[v_idx] = gray;

		vector<VertexIndex> v_idxs;
		MyMesh::ConstVertexVertexIter vv_it;
		MyMesh::ConstVertexVertexIter vv_end = _mesh.vv_end(*v_it);
		for (vv_it = _mesh.vv_begin(*v_it); vv_it != vv_end; ++vv_it)
		{
			VertexIndex v_idx = vv_it->idx();
			v_idxs.push_back(v_idx);
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
	for (f_it = _mesh.faces_begin(); f_it != f_end; ++f_end)
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
	if (_sh_order > -1)
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

	return shading;
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
			est_intensity[0] = 
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

		MyMesh::Color c(_colors[i][0], _colors[i][1], _colors[i][2]);
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

		MyMesh::Color c(_intensities[i], _intensities[i], _intensities[i]);
		mesh.set_color(*v_it, c);
	}
}
//=============================================================================
Intensity BlackBox::rgb2gray(const Color &_color)
{
	Intensity gray = RGB2GRAY_R_FACTOR * _color[0] 
		+ RGB2GRAY_G_FACTOR * _color[1] 
		+ RGB2GRAY_B_FACTOR * _color[2];

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

				int neigh_face_index = it - adjFacesVerticesInd.begin();
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

					int neigh_face_index = it - adjFacesVerticesInd.begin();
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
	diff_weights.resize(n_vertices);
	for (size_t i = 0; i < n_vertices; i++)
	{
		vector<Intensity> weights;
		weights.reserve(adj_vertices[i].size());
		for (size_t j = 0; j < adj_vertices[i].size(); j++)
		{
			unsigned int adj_v_idx = adj_vertices[i][j];
			Intensity weight = computeDiffWeight(colors[i], colors[adj_v_idx],
				normals[i], normals[adj_v_idx], params.color_diff_threshold,
				params.color_diff_std, params.normal_diff_std);
			weights.push_back(weight);
		}
		diff_weights[i] = std::move(weights);
	}
}
//=============================================================================
Intensity BlackBox::computeDiffWeight(Color &_color1, Color &_color2,
	Normal &_normal1, Normal &_normal2, Intensity _color_diff_threshold,
	Intensity _color_diff_std, Coordinate _normal_diff_std)
{
	Intensity weight = 0.0;

	Color color_diff = { _color1[0] - _color2[0],
		_color1[1] - _color2[1], _color1[2] - _color2[2] };
	Intensity color_diff_norm2 = pow(color_diff[0], 2) + pow(color_diff[1], 2)
		+ pow(color_diff[2], 2);

	if (color_diff_norm2 <= _color_diff_threshold)
	{
		weight = exp(color_diff_norm2 / (2 * pow(_color_diff_std, 2)));

		Coordinate normal_diff_cos = _normal1[0] * _normal2[0]
			+ _normal1[1] * _normal2[1] + _normal1[2] * _normal2[2];
		weight *= exp((1 - normal_diff_cos) / (2 * pow(_normal_diff_std, 2)));
	}

	return weight;
}
//=============================================================================
void BlackBox::estimateSHCoeff(const ceres::Solver::Options &_options)
{
	Color white = { 1.0, 1.0, 1.0 };
	Color black = { 0.0, 0.0, 0.0 };

	ceres::Problem problem;

	for (size_t i = 0; i < n_vertices; i++)
	{
		ResidualPhotometricError *residual = new ResidualPhotometricError(
			&grays[i], 1, (unsigned int)adj_vertices[i].size(),
			(unsigned int)adj_faces[i].size(), sh_order);

		ceres::DynamicAutoDiffCostFunction<ResidualPhotometricError, 5>* dyn_cost_function
			= new ceres::DynamicAutoDiffCostFunction< ResidualPhotometricError, 5 >(residual);

		// List of pointers to translations per vertex
		vector<Intensity*> v_parameter_blocks;

		// White albedo
		dyn_cost_function->AddParameterBlock(1);
		v_parameter_blocks.push_back(&white[0]);
		// SH Coeff
		dyn_cost_function->AddParameterBlock(n_sh_basis);
		v_parameter_blocks.push_back(&sh_coeff[0]);
		// Local lighting variations
		dyn_cost_function->AddParameterBlock(1);
		v_parameter_blocks.push_back(&black[0]);
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

		dyn_cost_function->SetNumResiduals(1);

		ceres::ResidualBlockId residualBlockId = problem.AddResidualBlock(
			dyn_cost_function,
			NULL,
			v_parameter_blocks);
	}

	problem.SetParameterBlockConstant(&white[0]);
	problem.SetParameterBlockConstant(&black[0]);
	for (size_t i = 0; i < n_vertices; ++i)
	{
		problem.SetParameterBlockConstant(&vertices[i][0]);
	}

	ceres::Solver::Summary summary;

	ceres::Solve(_options, &problem, &summary);

	cout << summary.FullReport() << endl;
}
//=============================================================================
void BlackBox::estimateAlbedo(const ceres::Solver::Options &_options)
{
	Color black = { 0.0, 0.0, 0.0 };

	ceres::Problem problem;

	for (size_t i = 0; i < n_vertices; i++)
	{
		ResidualPhotometricError *residual = new ResidualPhotometricError(
			&colors[i][0], 3, (unsigned int)adj_vertices[i].size(),
			(unsigned int)adj_faces[i].size(), sh_order);

		ceres::DynamicAutoDiffCostFunction<ResidualPhotometricError, 5>* dyn_cost_function
			= new ceres::DynamicAutoDiffCostFunction< ResidualPhotometricError, 5 >(residual);

		// List of pointers to translations per vertex
		vector<Intensity*> v_parameter_blocks;

		// White albedo
		dyn_cost_function->AddParameterBlock(3);
		v_parameter_blocks.push_back(&albedos[i][0]);
		// SH Coeff
		dyn_cost_function->AddParameterBlock(n_sh_basis);
		v_parameter_blocks.push_back(&sh_coeff[0]);
		// Local lighting variations
		dyn_cost_function->AddParameterBlock(3);
		v_parameter_blocks.push_back(&black[0]);
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

		dyn_cost_function->SetNumResiduals(1);

		ceres::ResidualBlockId residualBlockId = problem.AddResidualBlock(
			dyn_cost_function,
			NULL,
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

	problem.SetParameterBlockConstant(&sh_coeff[0]);
	problem.SetParameterBlockConstant(&black[0]);
	for (size_t i = 0; i < n_vertices; ++i)
	{
		problem.SetParameterBlockConstant(&vertices[i][0]);
	}

	ceres::Solver::Summary summary;

	ceres::Solve(_options, &problem, &summary);

	cout << summary.FullReport() << endl;
}
//=============================================================================
void BlackBox::estimateLocalLighting(const ceres::Solver::Options &_options)
{
	ceres::Problem problem;

	for (size_t i = 0; i < n_vertices; i++)
	{
		ResidualPhotometricError *residual = new ResidualPhotometricError(
			&colors[i][0], 3, (unsigned int)adj_vertices[i].size(),
			(unsigned int)adj_faces[i].size(), sh_order);

		ceres::DynamicAutoDiffCostFunction<ResidualPhotometricError, 5>* dyn_cost_function
			= new ceres::DynamicAutoDiffCostFunction< ResidualPhotometricError, 5 >(residual);

		// List of pointers to translations per vertex
		vector<Intensity*> v_parameter_blocks;

		// White albedo
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
			NULL,
			v_parameter_blocks);
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

	problem.SetParameterBlockConstant(&sh_coeff[0]);
	for (size_t i = 0; i < n_vertices; ++i)
	{
		problem.SetParameterBlockConstant(&vertices[i][0]);
		problem.SetParameterBlockConstant(&albedos[i][0]);
	}

	ceres::Solver::Summary summary;

	ceres::Solve(_options, &problem, &summary);

	cout << summary.FullReport() << endl;
}
//=============================================================================
void BlackBox::estimateShape(const ceres::Solver::Options &_options)
{
	ceres::Problem problem;

	for (size_t i = 0; i < n_vertices; i++)
	{
		ResidualPhotometricError *residual = new ResidualPhotometricError(
			&colors[i][0], 3, (unsigned int)adj_vertices[i].size(),
			(unsigned int)adj_faces[i].size(), sh_order);

		ceres::DynamicAutoDiffCostFunction<ResidualPhotometricError, 5>* dyn_cost_function
			= new ceres::DynamicAutoDiffCostFunction< ResidualPhotometricError, 5 >(residual);

		// List of pointers to translations per vertex
		vector<Intensity*> v_parameter_blocks;

		// White albedo
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
			NULL,
			v_parameter_blocks);
	}

	if (params.displacement_weight > 0.0)
	{
		ceres::ScaledLoss* loss_function = new ceres::ScaledLoss(
			NULL, params.displacement_weight, ceres::TAKE_OWNERSHIP);
		for (size_t i = 0; i < n_vertices; i++)
		{
			ResidualInitialVertexDiff *residual = new ResidualInitialVertexDiff(vertices[i]);

			ceres::AutoDiffCostFunction<ResidualInitialVertexDiff, 3, 3>* cost_function =
				new ceres::AutoDiffCostFunction<ResidualInitialVertexDiff, 3, 3>(residual);

			ceres::ResidualBlockId residualBlockId = problem.AddResidualBlock(
				cost_function,
				loss_function,
				&vertices[i][0]);
		}
	}

	if (params.laplacian_smoooth_weight > 0.0)
	{
		ceres::ScaledLoss* loss_function = new ceres::ScaledLoss(
			NULL, params.laplacian_smoooth_weight, ceres::TAKE_OWNERSHIP);
		for (size_t i = 0; i < n_vertices; i++)
		{
			ResidualLaplacianSmoothing *residual 
				= new ResidualLaplacianSmoothing((unsigned int)adj_vertices[i].size());

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

			dyn_cost_function->SetNumResiduals(3);

			ceres::ResidualBlockId residualBlockId = problem.AddResidualBlock(
				dyn_cost_function,
				loss_function,
				v_parameter_blocks);
		}
	}

	problem.SetParameterBlockConstant(&sh_coeff[0]);
	for (size_t i = 0; i < n_vertices; ++i)
	{
		problem.SetParameterBlockConstant(&albedos[i][0]);
		problem.SetParameterBlockConstant(&local_lightings[i][0]);
	}

	ceres::Solver::Summary summary;

	ceres::Solve(_options, &problem, &summary);

	cout << summary.FullReport() << endl;
}
//=============================================================================