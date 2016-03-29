#pragma once
#include <OpenMesh\Core\IO\MeshIO.hh>
#include <OpenMesh\Core\Mesh\PolyMesh_ArrayKernelT.hh>
#include <Eigen\Dense>
//=============================================================================
struct Traits
{
	/// The default coordinate type is OpenMesh::Vec3d.
	typedef OpenMesh::Vec3d  Point;

	/// The default normal type is OpenMesh::Vec3d.
	typedef OpenMesh::Vec3d  Normal;

	/// The default 1D texture coordinate type is float.
	typedef double  TexCoord1D;
	/// The default 2D texture coordinate type is OpenMesh::Vec2d.
	typedef OpenMesh::Vec2d  TexCoord2D;
	/// The default 3D texture coordinate type is OpenMesh::Vec3d.
	typedef OpenMesh::Vec3d  TexCoord3D;

	/// The default texture index type
	typedef int TextureIndex;

	/// The default color type is OpenMesh::Vec3uc.
	typedef OpenMesh::Vec3d Color;

#ifndef DOXY_IGNORE_THIS
	VertexTraits{};
	HalfedgeTraits{};
	EdgeTraits{};
	FaceTraits{};
#endif

	VertexAttributes(OpenMesh::Attributes::Status);
	HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
	EdgeAttributes(OpenMesh::Attributes::Status);
	FaceAttributes(OpenMesh::Attributes::Status);
};
//=============================================================================
typedef OpenMesh::PolyMesh_ArrayKernelT<Traits> MyMesh;
//=============================================================================
