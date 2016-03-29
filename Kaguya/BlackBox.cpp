#include "BlackBox.h"
//=============================================================================
parameters::Parameters BlackBox::params;

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

vector<Intensity> BlackBox::sh_coeff;
vector<Color> BlackBox::albedos;
vector<Color> BlackBox::lighting_variations;

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

	initFromMesh(mesh);
}
//=============================================================================
void BlackBox::run()
{


	//if (params.save_shading)
	//{
	//	setMeshColors(shading, mesh);

	//	writeMesh(mesh, params.output_shading_mesh_filename.c_str(), mesh_write_opt);
	//}

	//if (params.save_albedo)
	//{
	//	setMeshColors(albedo, mesh);

	//	writeMesh(mesh, params.output_albedo_mesh_filename.c_str(), mesh_write_opt);
	//}

	//if (params.save_est_intensity)
	//{
	//	MatrixXf intensities = albedo.cwiseProduct(shading.replicate(1, 3));

	//	setMeshColors(intensities, mesh);

	//	writeMesh(mesh, params.output_intensity_mesh_filename.c_str(), mesh_write_opt);
	//}

	if (params.save_sh_coeff)
	{
		writeVector(sh_coeff, params.output_sh_coeff_filename.c_str());
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
void BlackBox::computeSHFunctions(const MatrixXf &_normals, const int _sh_order,
	MatrixXf &_sh_functions)
{
	_sh_functions.resize(_normals.rows(), (_sh_order + 1) * (_sh_order + 1));

	VectorXf x = _normals.col(0);						// x
	VectorXf y = _normals.col(1);						// y
	VectorXf z = _normals.col(2);						// z
	VectorXf x2 = x.cwiseProduct(x);					// x^2
	VectorXf y2 = y.cwiseProduct(y);					// y^2
	VectorXf z2 = z.cwiseProduct(z);					// z^2
	VectorXf xy = x.cwiseProduct(y);					// x * y
	VectorXf xz = x.cwiseProduct(z);					// x * z
	VectorXf yz = y.cwiseProduct(z);					// y * z
	VectorXf x2_y2 = x2 - y2;							// x^2 - y^2
	VectorXf ones = VectorXf::Ones(_normals.rows());	// 1

	if (_sh_order > -1)
	{
		_sh_functions.col(0) = ones;
	}

	if (_sh_order > 0)
	{
		_sh_functions.block(0, 1, _normals.rows(), 3) = _normals;	// x, y, z
	}

	if (_sh_order > 1)
	{
		_sh_functions.col(4) = xy;				// x * y
		_sh_functions.col(5) = xz;				// x * z
		_sh_functions.col(6) = yz;				// y * z
		_sh_functions.col(7) = x2_y2;			// x^2 - y^2
		_sh_functions.col(8) = 3 * z2 - ones;	// 3 * z^2 - 1
	}

	if (_sh_order > 2)
	{
		_sh_functions.col(9) = (3 * x2 - y2).cwiseProduct(y);			// (3 * x^2 - y^2) * y 
		_sh_functions.col(10) = xy.cwiseProduct(z);						// x * y * z
		_sh_functions.col(11) = (5 * z2 - ones).cwiseProduct(y);		// (5 * z^2 - 1) * y
		_sh_functions.col(12) = (5 * z2 - 3 * ones).cwiseProduct(z);	// (5 * z^2 - 3) * z
		_sh_functions.col(13) = (5 * z2 - ones).cwiseProduct(x);		// (5 * z^2 - 1) * x
		_sh_functions.col(14) = x2_y2.cwiseProduct(z);					// (x^2 - y^2) * z
		_sh_functions.col(15) = (x2 - 3 * y2).cwiseProduct(x);			// (x^2 - 3 * y^2) * x
	}

	if (_sh_order > 3)
	{
		_sh_functions.col(16) = x2_y2.cwiseProduct(xy);							// (x^2 - y^2) * x * y
		_sh_functions.col(17) = (3 * x2 - y2).cwiseProduct(yz);					// (3 * x^2 - y^2) * yz
		_sh_functions.col(18) = (7 * z2 - ones).cwiseProduct(xy);				// (7 * z^2 - 1) * x * y
		_sh_functions.col(19) = (7 * z2 - 3 * ones).cwiseProduct(yz);			// (7 * z^2 - 3) * y * z
		_sh_functions.col(20) = 3 * ones - 30 * z2 + 35 * z2.cwiseProduct(z2);	// 3 - 30 * z^2 + 35 * z^4
		_sh_functions.col(21) = (7 * z2 - 3 * ones).cwiseProduct(xz);			// (7 * z^2 - 3) * x * z
		_sh_functions.col(22) = (7 * z2 - ones).cwiseProduct(x2_y2);			// (7 * z^2 - 1) * (x^2 - y^2)
		_sh_functions.col(23) = (x2 - 3 * y2).cwiseProduct(xz);					// (x^2 - 3 * y^2) * x * z
		_sh_functions.col(24) = (x2 - 3 * y2).cwiseProduct(x2)					// (x^2 - 3 * y^2) * x^2 - (3 * x^2 - y^2) * y^2 
			- (3 * x2 - y2).cwiseProduct(y2);
	}
}
//=============================================================================
VectorXf BlackBox::getSHCoeff(const VectorXf &_shading, const MatrixXf &_normals,
	const int _sh_order)
{
	MatrixXf sh_functions;
	computeSHFunctions(_normals, _sh_order, sh_functions);

	JacobiSVD<MatrixXf> sh_functions_svd(sh_functions, ComputeThinU | ComputeThinV);
	VectorXf sh_coeff = sh_functions_svd.solve(_shading);

	return sh_coeff;
}
//=============================================================================
VectorXf BlackBox::getShading(const MatrixXf &_normals, const VectorXf &_sh_coeff)
{
	int sh_order = (int)sqrt(_sh_coeff.size()) - 1;
	
	MatrixXf sh_functions;
	computeSHFunctions(_normals, sh_order, sh_functions);

	VectorXf shading = sh_functions * _sh_coeff;

	return shading;
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