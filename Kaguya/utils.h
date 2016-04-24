#include "FileIO.h"
#include <Eigen\Dense>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\highgui\highgui.hpp>
//=============================================================================
using namespace Eigen;
//=============================================================================
void readMatrix3f(Eigen::Matrix3f &_matrix, const string _filename);
void computeVisibility(const int _img_width, const int _img_height,
	const Matrix3f &_intrinsics, MyMesh &_mesh, vector<bool> &_visibility);
bool pointInTriangleTest2(const Vector2f &pointP,
	const Vector2f &pointA, const Vector2f &pointB, const Vector2f &pointC);
bool visibilityTest(const Vector3f &vertex, const Vector3f &center,
	const Vector3f & normal, const Vector3f & vertex1,
	const Vector3f & vertex2, const Vector3f & vertex3);
//=============================================================================
inline void readMatrix3f(Eigen::Matrix3f &_matrix, const string _filename)
{
	std::ifstream ifs(_filename);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)
		{
			ifs >> _matrix(i, j);
		}
	}
}
//=============================================================================
inline void computeVisibility(const int _img_width, const int _img_height,
	const Matrix3f &_intrinsics, MyMesh &_mesh, vector<bool> &_visibility)
{
	MatrixXf points3d(3, _mesh.n_vertices());

	MyMesh::ConstVertexIter v_it;
	MyMesh::ConstVertexIter v_end(_mesh.vertices_end());

	for (v_it = _mesh.vertices_begin(); v_it != v_end; ++v_it)
	{
		int i = v_it.handle().idx();

		MyMesh::Point v = _mesh.point(*v_it);
		points3d(0, i) = v[0];
		points3d(1, i) = v[1];
		points3d(2, i) = v[2];
	}

	MatrixXf points2d_hom = _intrinsics * points3d;
	MatrixXf points2d(2, _mesh.n_vertices());
	points2d.row(0) = points2d_hom.row(0).cwiseProduct(points2d_hom.row(2).cwiseInverse());
	points2d.row(1) = points2d_hom.row(1).cwiseProduct(points2d_hom.row(2).cwiseInverse());



	MatrixXi faces(3, _mesh.n_faces());
	MatrixXf face_normals(3, _mesh.n_faces());

	MyMesh::ConstFaceIter f_it;
	MyMesh::ConstFaceIter f_end = _mesh.faces_end();

	for (f_it = _mesh.faces_begin(); f_it != f_end; ++f_it)
	{
		int f_idx = f_it.handle().idx();

		MyMesh::Normal f_normal = _mesh.normal(*f_it);

		face_normals(0, f_idx) = f_normal[0];
		face_normals(1, f_idx) = f_normal[1];
		face_normals(2, f_idx) = f_normal[2];

		MyMesh::Face face = _mesh.face(*f_it);

		MyMesh::ConstFaceVertexIter fv_it;
		MyMesh::ConstFaceVertexIter fv_end = _mesh.fv_end(*f_it);

		int i = 0;
		for (fv_it = _mesh.fv_begin(*f_it); fv_it != fv_end; ++fv_it)
		{
			int v_idx = fv_it.handle().idx();
			faces(i, f_idx) = v_idx;

			i++;
		}
	}


	static vector< vector<unsigned int> > visibilityFacesTest;

	// visibility mask, used for speed up old rendering code
	// this will set up in the first call
	// will not do anything during the second time
	visibilityFacesTest.resize(_img_width*_img_height);
	for (int i = 0; i < _img_width*_img_height; ++i)
		visibilityFacesTest[i].reserve(5);   // 5 reserved faces for each pixel.

	// image width and height
	int m_nWidth, m_nHeight;
	m_nWidth = _img_width; m_nHeight = _img_height;

	// reset the visibility mask to true first.
	_visibility.resize(_mesh.n_vertices(), true);

	// update the visibility mask using previous mesh
	// assign faces to visibilityFacesTest according to 2d projections
	double x1, x2, x3, y1, y2, y3;
	double xmin, xmax, ymin, ymax;
	int xminI, xmaxI, yminI, ymaxI;
	int xx, yy;
	int vertex1, vertex2, vertex3;;
	MatrixXf faceCenters(3, _mesh.n_faces());

	for (int faceInd = 0; faceInd < _mesh.n_faces(); faceInd++)
	{
		vertex1 = faces(0, faceInd);
		vertex2 = faces(1, faceInd);
		vertex3 = faces(2, faceInd);

		x1 = points2d(0, vertex1);
		x2 = points2d(0, vertex2);
		x3 = points2d(0, vertex3);
		y1 = points2d(1, vertex1);
		y2 = points2d(1, vertex2);
		y3 = points2d(1, vertex3);
		xmin = min(min(x1, x2), x3);
		xmax = max(max(x1, x2), x3);
		ymin = min(min(y1, y2), y3);
		ymax = max(max(y1, y2), y3);
		// if the projection of the face is outside the image, we give up
		if (xmin > m_nWidth + 1 || ymin > m_nHeight + 1 ||
			xmax < 0 || ymax < 0)
			continue;

		xminI = max(1.0, floor(xmin));
		yminI = max(1.0, floor(ymin));
		xmaxI = min(ceil(xmax), double(m_nWidth));
		ymaxI = min(ceil(ymax), double(m_nHeight));
		for (xx = xminI - 1; xx < xmaxI; ++xx)
		{
			for (yy = yminI - 1; yy < ymaxI; ++yy)
			{
				visibilityFacesTest[yy*m_nWidth + xx].push_back(faceInd);
			}
		}

		Vector3f v1 = points3d.col(vertex1);
		Vector3f v2 = points3d.col(vertex2);
		Vector3f v3 = points3d.col(vertex3);

		// prepare the face center and face normals for the purpose of gettin
		// intersection between faces and back-projected ray.
		faceCenters.col(faceInd) = (v1 + v2 + v3) / 3;

		//compnorm(&meshVertices[vertex1][0], &meshVertices[vertex2][0],
		//	&meshVertices[vertex3][0], &faceNormals[3 * faceInd], 1);
	}

	// // check if the projection of all faces looks reasonable
	// // PixelType* pMaskBuffer = new PixelType[m_nHeight*m_nWidth];
	// // IntensityImageType maskImage(m_nHeight,m_nWidth,pMaskBuffer);
	// IntensityImageType maskImage(m_nHeight,m_nWidth);
	// for(int i = 0; i < m_nHeight; ++i)
	// {
	//     for(int j = 0; j < m_nWidth; ++j)
	//     {
	//         if(visibilityFacesTest[i*m_nWidth+j].size() > 0)
	//         maskImage(i,j) = 255;
	//         else
	//         maskImage(i,j) = 0;
	//     }
	// }
	// // cv::imshow("maskImage",maskImage);
	// // cv::waitKey(0);
	// cv::imwrite("visMaskImage.png",maskImage);
	// // delete[] pMaskBuffer;

	// loop over vertices to find the intersection between backprojected race
	// and the possible faces of the mesh
	double projX, projY;
	int projXminI, projYminI;
	int projXmaxI, projYmaxI;
	int position, faceInd;

	for (unsigned int vertexInd = 0; vertexInd < _mesh.n_vertices();
		++vertexInd)
	{
		// if the projection is outside the image, it is not visible

		if (points2d(0, vertexInd) < 0 || points2d(0, vertexInd) > m_nWidth ||
			points2d(1, vertexInd) < 0 || points2d(1, vertexInd) > m_nHeight)
		{
			_visibility[vertexInd] = false;
			continue;
		}

		projX = points2d(0, vertexInd);
		projY = points2d(1, vertexInd);

		projXminI = max(1.0, floor(projX));
		projYminI = max(1.0, floor(projY));

		projXmaxI = min(ceil(projX), double(m_nWidth));
		projYmaxI = min(ceil(projY), double(m_nHeight));

		bool done = false;

		for (xx = projXminI - 1; xx < projXmaxI; ++xx)
		{
			for (yy = projYminI - 1; yy < projYmaxI; ++yy)
			{
				// loop over all the faces here
				position = yy*m_nWidth + xx;
				int facesNum = visibilityFacesTest[position].size();
				for (int k = 0; k < facesNum; ++k)
				{
					faceInd = visibilityFacesTest[position][k];

					Vector3i face = faces.col(faceInd);

					vertex1 = face(0);
					vertex2 = face(1);
					vertex3 = face(2);
					// check if the vertex belongs to this face
					if (vertex1 == vertexInd ||
						vertex2 == vertexInd ||
						vertex3 == vertexInd)
						continue;
					// First test if the projection of the point is inside the triangle,
					// If not, we continue to the next face
					if (pointInTriangleTest2(points2d.col(vertexInd),
						points2d.col(vertex1),
						points2d.col(vertex2),
						points2d.col(vertex3)))
					{
						Vector3f v0 = points3d.col(vertexInd);
						Vector3f v1 = points3d.col(vertex1);
						Vector3f v2 = points3d.col(vertex2);
						Vector3f v3 = points3d.col(vertex3);
						Vector3f fn = face_normals.col(faceInd);

						// compute the intersection between the
						// backtraced ray of vertexInd and the face faceInd
						// point, center, normal, point1, point2, point3
						_visibility[vertexInd] = visibilityTest(v0,
							faceCenters.col(faceInd),
							fn,
							v1,
							v2,
							v3);
					}
					else
						continue;

					// check if the vertex is already occluded or not
					if (_visibility[vertexInd] == false)
					{
						done = true;
						break;
					}
				}

				if (done)
					break;

			}

			if (done)
				break;

		}

		//cout << "Vertex " << vertexInd << ": " << visibility[vertexInd] << endl;

	}

	for (int i = 0; i < visibilityFacesTest.size(); ++i)
		visibilityFacesTest[i].clear();   // 5 reserved faces for each pixel, capacity doesn't change

	// print out how many points are visible and how many points are occluded
	int n_visible = 0;
	int n_occluded = 0;
	for (int i = 0; i < _visibility.size(); ++i)
	{
		if (_visibility[i])
			n_visible++;
		else
			n_occluded++;
	}

	cout << "visible num: " << n_visible << endl;
	cout << "occluded num: " << n_occluded << endl << endl;
}
//=============================================================================
inline bool pointInTriangleTest2(const Vector2f &pointP,
	const Vector2f &pointA, const Vector2f &pointB, const Vector2f &pointC)
{
	float A, B, D, E, C, F;
	A = pointA(0) - pointC(0);
	B = pointB(0) - pointC(0);
	D = pointA(1) - pointC(1);
	E = pointB(1) - pointC(1);

	C = pointC(0) - pointP(0);
	F = pointC(1) - pointP(1);

	float alpha, beta, gamma;
	alpha = (B*(F)-C*(E)) / (A*(E)-B*(D));
	beta = (A*(F)-C*(D)) / (B*(D)-A*(E));
	gamma = 1 - alpha - beta;

	return (alpha >= 0 && beta >= 0 && gamma >= 0);

}
//=============================================================================
inline bool visibilityTest(const Vector3f &vertex, const Vector3f &center,
	const Vector3f & normal, const Vector3f & vertex1,
	const Vector3f & vertex2, const Vector3f & vertex3)
{
	// tell if a point is in front of a triangle face
	float faceDist, pointDist, scale;
	faceDist = center.dot(normal);
	pointDist = vertex.dot(normal);
	scale = faceDist / pointDist;

	Vector3f intersection = scale*vertex;

	// if the intersection is in front,
	// then the test vertex is occluded by
	// the front face, we give up
	if (intersection(2) < vertex(2))
		return false;
	else
		return true;

}
//=============================================================================