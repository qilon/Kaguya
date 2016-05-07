#pragma once

#include <fstream> 

#include "MyMesh.h"
//=============================================================================
bool writeMesh(const MyMesh& _mesh, const string &_filename,
	const OpenMesh::IO::Options _wopt);
//=============================================================================
bool readMesh(MyMesh& _mesh, const string &_filename,
	OpenMesh::IO::Options _ropt);
//=============================================================================
inline bool readMesh(MyMesh& _mesh, const string &_filename,
	OpenMesh::IO::Options _ropt)
{
	if (!OpenMesh::IO::read_mesh(_mesh, _filename, _ropt))
	{
		std::cout << "Could not read file: " << _filename << std::endl << std::endl;

		return false;
	}

	return true;
}
//=============================================================================
inline bool writeMesh(const MyMesh& _mesh, const string &_filename,
	const OpenMesh::IO::Options _wopt)
{
	if (!OpenMesh::IO::write_mesh(_mesh, _filename, _wopt))
	{
		std::cout << "Could not write file: " << _filename << std::endl << std::endl;

		return false;
	}

	return true;
}
//=============================================================================
template<typename T>
int writeVector(const std::vector<T> &_v, const string &_filename)
{
	std::ofstream ofs;
	ofs.open(_filename, std::ofstream::out | std::ofstream::trunc);
	for (int i = 0; i < _v.size(); i++)
	{
		ofs << _v[i] << std::endl;
	}
	ofs.close();

	return 0;
}
//=============================================================================
template<typename T>
int readVector(std::vector<T> &_v, const string &_filename)
{
	std::ifstream ifs(_filename);

	for (int i = 0; i < _v.size(); i++)
	{
		ifs >> _v[i];
	}

	ifs.close();

	return 0;
}
//=============================================================================