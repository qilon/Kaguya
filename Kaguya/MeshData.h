#include <vector>
#include "ply.h"

using namespace std;

struct MeshData
{
	typedef float VertexType;
	typedef vector<VertexType> VertexT;
	typedef vector<VertexType> NormalT;
	typedef unsigned int FaceType;
	typedef vector<FaceType> FaceT;
	typedef float ColorType;
	typedef vector<ColorType> ColorT;

	vector< VertexT > vertices;
	vector< VertexT > vertex_normals;
	vector< ColorT > vertex_colors;
	vector< ColorT > vertex_specular_colors;
	vector< FaceT > faces;
	vector< VertexT > face_normals;
	vector<VertexType> sh_coefficients;

	size_t n_vertices() const
	{
		return vertices.size();
	}

	size_t n_faces() const
	{
		return faces.size();
	}

	bool has_vertex_normals() const
	{
		return vertex_normals.size() == n_vertices();
	}

	bool has_vertex_colors() const
	{
		return vertex_colors.size() == n_vertices();
	}

	bool has_vertex_specular_colors() const
	{
		return vertex_specular_colors.size() == n_vertices();
	}

	size_t n_sh_coefficients() const
	{
		return sh_coefficients.size();
	}

	void readPLY(const std::string& filename)
	{
		PlyFile *ply;
		char **elist;
		int file_type;
		float version;
		int nelems, num_elems, nprops;
		char *elem_name;
		PlyProperty **plist;
		int num_comments;
		char **comments;
		int num_obj_info;
		char **obj_info;

		/* open a PLY file for reading */
		ply = ply_open_for_reading(filename.c_str(), &nelems, &elist, &file_type, &version);

		/* go through each kind of element that we learned is in the file */
		/* and read them */

		unsigned int vertex_type;

		for (int i = 0; i < nelems; i++) {

			/* get the description of the first element */
			elem_name = elist[i];
			plist = ply_get_element_description(ply, elem_name, &num_elems, &nprops);

			/* if we're on vertex elements, read them in */
			if (equal_strings(ply::elem_names[0], elem_name)) {
				int num_vertices = num_elems;

				vertex_type = ply::get_vertex_type(plist, nprops);

				switch (vertex_type)
				{
				case PLY_VERTEX_RGB:
					ply_read_vertices<ply::VertexColor>(ply, vertex_type, num_vertices,
						*this);
					break;
				case PLY_VERTEX_NORMAL_RGB:
					ply_read_vertices<ply::VertexNormalColor>(ply, vertex_type, num_vertices,
						*this);
					break;
				case PLY_VERTEX_RGBA:
					ply_read_vertices<ply::VertexColorAlpha>(ply, vertex_type, num_vertices,
						*this);
					break;
				case PLY_VERTEX_NORMAL_RGBA:
					ply_read_vertices<ply::VertexNormalColorAlpha>(ply, vertex_type, num_vertices,
						*this);
					break;
				case PLY_VERTEX_RGB_SPECULAR:
					ply_read_vertices<ply::VertexColorSpecular>(ply, vertex_type, num_vertices,
						*this);
					break;
				case PLY_VERTEX_NORMAL_RGB_SPECULAR:
					ply_read_vertices<ply::VertexNormalColorSpecular>(ply, vertex_type, num_vertices,
						*this);
					break;
				case PLY_VERTEX_RGBA_SPECULAR:
					ply_read_vertices<ply::VertexColorAlphaSpecular>(ply, vertex_type, num_vertices,
						*this);
					break;
				case PLY_VERTEX_NORMAL_RGBA_SPECULAR:
					ply_read_vertices<ply::VertexNormalColorAlphaSpecular>(ply, vertex_type, num_vertices,
						*this);
					break;
				default:
					break;
				}
			}

			/* if we're on face elements, read them in */
			if (equal_strings(ply::elem_names[1], elem_name)) {
				int num_faces = num_elems;

				/* set up for getting face elements */
				for (int i = 0; i < nprops; i++)
				{
					ply_get_property(ply, elem_name, &ply::face_props[i]);
				}

				faces.resize(num_faces);

				/* grab all the face elements */
				for (int j = 0; j < num_faces; j++) {
					ply::Face face;
					/* grab and element from the file */
					ply_get_element(ply, (void *)&face);

					vector<unsigned int> f;
					for (int k = 0; k < face.nverts; k++)
					{
						f.push_back(face.verts[k]);
					}
					faces[j] = f;
				}
			}

			/* if we're on sh coefficient elements, read them in */
			if (equal_strings(ply::elem_names[2], elem_name)) {
				int num_sh_coefficients = num_elems;

				/* set up for getting face elements */
				for (int i = 0; i < nprops; i++)
				{
					ply_get_property(ply, elem_name, &ply::sh_coeff_props[i]);
				}

				sh_coefficients.resize(num_sh_coefficients);

				/* grab all the sh coefficient elements */
				for (int j = 0; j < num_sh_coefficients; j++) {
					ply::SH_Coefficient sh_coeff;
					/* grab and element from the file */
					ply_get_element(ply, (void *)&sh_coeff);

					sh_coefficients[j] = sh_coeff.value;
				}
			}
		}

		/* grab and print out the comments in the file */
		comments = ply_get_comments(ply, &num_comments);

		/* grab and print out the object information */
		obj_info = ply_get_obj_info(ply, &num_obj_info);

		/* close the PLY file */
		ply_close(ply);

		for (int i = 0; i < nelems; i++)
		{
			delete elist[i];
		}
		delete[] elist;

		//if (vertex_type == PLY_VERTEX_RGBA || vertex_type == PLY_VERTEX_RGB
		//	|| vertex_type == PLY_VERTEX_RGB_SPECULAR
		//	|| vertex_type == PLY_VERTEX_RGBA_SPECULAR)
		//meshData.computeNormalsNeil();
	}

	//=============================================================================
	template<typename PlyVertexType>
	void ply_read_vertices(PlyFile *_ply, int _vertex_type, int _num_vertices,
		MeshData& meshData)
	{
		/* set up for getting vertex elements */
		for (int i = 0; i < ply::n_vprops[_vertex_type]; i++)
		{
			PlyProperty* prop = ply::get_vertex_property(_vertex_type, i);

			ply_get_property(_ply, ply::elem_names[0], prop);
		}

		meshData.vertices.resize(_num_vertices);
		meshData.vertex_colors.resize(_num_vertices);

		if (_vertex_type & PLY_WITH_NORMAL)
			meshData.vertex_normals.resize(_num_vertices);

		if (_vertex_type & PLY_WITH_SPECULAR)
			meshData.vertex_specular_colors.resize(_num_vertices);


		/* grab all the vertex elements */
		for (int i = 0; i < _num_vertices; i++) {
			PlyVertexType vertex;

			/* grab an element from the file */
			ply_get_element(_ply, (void *)&vertex);

			vector<VertexType> v = { (VertexType)vertex.x, (VertexType)vertex.y, (VertexType)vertex.z };
			meshData.vertices[i] = std::move(v);

			if (_vertex_type == PLY_VERTEX_NORMAL_RGB)
			{
				vector<VertexType> n = { (VertexType)((ply::VertexNormalColor*)&vertex)->nx,
					(VertexType)((ply::VertexNormalColor*)&vertex)->ny,
					(VertexType)((ply::VertexNormalColor*)&vertex)->nz };
				meshData.vertex_normals[i] = std::move(n);
			}

			if (_vertex_type == PLY_VERTEX_NORMAL_RGBA)
			{
				vector<VertexType> n = { (VertexType)((ply::VertexNormalColorAlpha*)&vertex)->nx,
					(VertexType)((ply::VertexNormalColorAlpha*)&vertex)->ny,
					(VertexType)((ply::VertexNormalColorAlpha*)&vertex)->nz };
				meshData.vertex_normals[i] = std::move(n);
			}

			if (_vertex_type == PLY_VERTEX_NORMAL_RGB_SPECULAR)
			{
				vector<VertexType> n = { (VertexType)((ply::VertexNormalColorSpecular*)&vertex)->nx,
					(VertexType)((ply::VertexNormalColorSpecular*)&vertex)->ny,
					(VertexType)((ply::VertexNormalColorSpecular*)&vertex)->nz };
				meshData.vertex_normals[i] = std::move(n);

				vector<ColorType> specular_c = {
					(ColorType)((ply::VertexNormalColorSpecular*)&vertex)->specular_r,
					(ColorType)((ply::VertexNormalColorSpecular*)&vertex)->specular_g,
					(ColorType)((ply::VertexNormalColorSpecular*)&vertex)->specular_b
				};
				specular_c[0] /= 255.f;
				specular_c[1] /= 255.f;
				specular_c[2] /= 255.f;
				meshData.vertex_specular_colors[i] = std::move(specular_c);
			}

			if (_vertex_type == PLY_VERTEX_NORMAL_RGBA_SPECULAR)
			{
				vector<VertexType> n = { (VertexType)((ply::VertexNormalColorAlphaSpecular*)&vertex)->nx,
					(VertexType)((ply::VertexNormalColorAlphaSpecular*)&vertex)->ny,
					(VertexType)((ply::VertexNormalColorAlphaSpecular*)&vertex)->nz };
				meshData.vertex_normals[i] = std::move(n);

				vector<ColorType> specular_c = {
					(ColorType)((ply::VertexNormalColorAlphaSpecular*)&vertex)->specular_r,
					(ColorType)((ply::VertexNormalColorAlphaSpecular*)&vertex)->specular_g,
					(ColorType)((ply::VertexNormalColorAlphaSpecular*)&vertex)->specular_b
				};
				specular_c[0] /= 255.f;
				specular_c[1] /= 255.f;
				specular_c[2] /= 255.f;
				meshData.vertex_specular_colors[i] = std::move(specular_c);
			}

			if (_vertex_type == PLY_VERTEX_RGB_SPECULAR)
			{
				vector<ColorType> specular_c = {
					(ColorType)((ply::VertexColorSpecular*)&vertex)->specular_r,
					(ColorType)((ply::VertexColorSpecular*)&vertex)->specular_g,
					(ColorType)((ply::VertexColorSpecular*)&vertex)->specular_b
				};
				specular_c[0] /= 255.f;
				specular_c[1] /= 255.f;
				specular_c[2] /= 255.f;
				meshData.vertex_specular_colors[i] = std::move(specular_c);
			}

			if (_vertex_type == PLY_VERTEX_RGBA_SPECULAR)
			{
				vector<ColorType> specular_c = {
					(ColorType)((ply::VertexColorAlphaSpecular*)&vertex)->specular_r,
					(ColorType)((ply::VertexColorAlphaSpecular*)&vertex)->specular_g,
					(ColorType)((ply::VertexColorAlphaSpecular*)&vertex)->specular_b
				};
				specular_c[0] /= 255.f;
				specular_c[1] /= 255.f;
				specular_c[2] /= 255.f;
				meshData.vertex_specular_colors[i] = std::move(specular_c);
			}

			vector<ColorType> c = { (ColorType)vertex.r, (ColorType)vertex.g, (ColorType)vertex.b };
			c[0] /= (ColorType)255;
			c[1] /= (ColorType)255;
			c[2] /= (ColorType)255;
			meshData.vertex_colors[i] = std::move(c);
		}
	}
	//=============================================================================

	void writeToPLY(std::string& filename, bool save_normal, bool save_alpha,
		bool save_specular, bool save_mesh_binary)
	{
		PlyFile *ply;
		float version;

		/* open either a binary or ascii PLY file for writing */
		/* (the file will be called "test.ply" because the routines */
		/*  enforce the .ply filename extension) */

		int n_elems = save_specular ? 3 : 2;

		if (save_mesh_binary)
		{
			ply = ply_open_for_writing(filename.c_str(), 3, ply::elem_names, 
				PLY_BINARY_LE, &version);
		}
		else
		{
			ply = ply_open_for_writing(filename.c_str(), 3, ply::elem_names, 
				PLY_ASCII, &version);
		}

		/* describe what properties go into the vertex and face elements */

		int num_vertices = (int)n_vertices();
		ply_element_count(ply, ply::elem_names[0], num_vertices);
		int vertex_type = 0;
		if (save_normal)
		{
			vertex_type |= PLY_WITH_NORMAL;
		}
		if (save_alpha)
		{
			vertex_type |= PLY_WITH_ALPHA;
		}
		if (save_specular)
		{
			vertex_type |= PLY_WITH_SPECULAR;
		}

		for (int i = 0; i < ply::n_vprops[vertex_type]; i++)
		{
			PlyProperty* prop = ply::get_vertex_property(vertex_type, i);
			ply_describe_property(ply, ply::elem_names[0], prop);
		}

		int num_faces = (int)n_faces();
		ply_element_count(ply, ply::elem_names[1], num_faces);
		for (int i = 0; i < ply::n_fprops; i++)
		{
			ply_describe_property(ply, ply::elem_names[1], &ply::face_props[i]);
		}

		if (save_specular)
		{
			int num_sh_coefficients = (int)n_sh_coefficients();
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
				v.r = (unsigned char)(vertex_colors[i][0] * 255.f);
				v.g = (unsigned char)(vertex_colors[i][1] * 255.f);
				v.b = (unsigned char)(vertex_colors[i][2] * 255.f);
				ply_put_element(ply, (void *)&v);
			}
			break;
			case PLY_VERTEX_NORMAL_RGB:
			{
				ply::VertexNormalColor v;
				v.x = (float)vertices[i][0];
				v.y = (float)vertices[i][1];
				v.z = (float)vertices[i][2];
				v.nx = (float)vertex_normals[i][0];
				v.ny = (float)vertex_normals[i][1];
				v.nz = (float)vertex_normals[i][2];
				v.r = (unsigned char)(vertex_colors[i][0] * 255.f);
				v.g = (unsigned char)(vertex_colors[i][1] * 255.f);
				v.b = (unsigned char)(vertex_colors[i][2] * 255.f);
				ply_put_element(ply, (void *)&v);
			}
			break;
			case PLY_VERTEX_RGBA:
			{
				ply::VertexColorAlpha v;
				v.x = (float)vertices[i][0];
				v.y = (float)vertices[i][1];
				v.z = (float)vertices[i][2];
				v.r = (unsigned char)(vertex_colors[i][0] * 255.f);
				v.g = (unsigned char)(vertex_colors[i][1] * 255.f);
				v.b = (unsigned char)(vertex_colors[i][2] * 255.f);
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
				v.nx = (float)vertex_normals[i][0];
				v.ny = (float)vertex_normals[i][1];
				v.nz = (float)vertex_normals[i][2];
				v.r = (unsigned char)(vertex_colors[i][0] * 255.f);
				v.g = (unsigned char)(vertex_colors[i][1] * 255.f);
				v.b = (unsigned char)(vertex_colors[i][2] * 255.f);
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
				v.r = (unsigned char)(vertex_colors[i][0] * 255.f);
				v.g = (unsigned char)(vertex_colors[i][1] * 255.f);
				v.b = (unsigned char)(vertex_colors[i][2] * 255.f);
				v.specular_r = (unsigned char)(vertex_specular_colors[i][0] * 255.f);
				v.specular_g = (unsigned char)(vertex_specular_colors[i][1] * 255.f);
				v.specular_b = (unsigned char)(vertex_specular_colors[i][2] * 255.f);
				ply_put_element(ply, (void *)&v);
			}
			break;
			case PLY_VERTEX_NORMAL_RGB_SPECULAR:
			{
				ply::VertexNormalColorAlphaSpecular v;
				v.x = (float)vertices[i][0];
				v.y = (float)vertices[i][1];
				v.z = (float)vertices[i][2];
				v.nx = (float)vertex_normals[i][0];
				v.ny = (float)vertex_normals[i][1];
				v.nz = (float)vertex_normals[i][2];
				v.r = (unsigned char)(vertex_colors[i][0] * 255.f);
				v.g = (unsigned char)(vertex_colors[i][1] * 255.f);
				v.b = (unsigned char)(vertex_colors[i][2] * 255.f);
				v.specular_r = (unsigned char)(vertex_specular_colors[i][0] * 255.f);
				v.specular_g = (unsigned char)(vertex_specular_colors[i][1] * 255.f);
				v.specular_b = (unsigned char)(vertex_specular_colors[i][2] * 255.f);
				ply_put_element(ply, (void *)&v);
			}
			break;
			case PLY_VERTEX_RGBA_SPECULAR:
			{
				ply::VertexColorAlphaSpecular v;
				v.x = (float)vertices[i][0];
				v.y = (float)vertices[i][1];
				v.z = (float)vertices[i][2];
				v.r = (unsigned char)(vertex_colors[i][0] * 255.f);
				v.g = (unsigned char)(vertex_colors[i][1] * 255.f);
				v.b = (unsigned char)(vertex_colors[i][2] * 255.f);
				v.a = 255;
				v.specular_r = (unsigned char)(vertex_specular_colors[i][0] * 255.f);
				v.specular_g = (unsigned char)(vertex_specular_colors[i][1] * 255.f);
				v.specular_b = (unsigned char)(vertex_specular_colors[i][2] * 255.f);
				ply_put_element(ply, (void *)&v);
			}
			break;
			case PLY_VERTEX_NORMAL_RGBA_SPECULAR:
			{
				ply::VertexNormalColorAlphaSpecular v;
				v.x = (float)vertices[i][0];
				v.y = (float)vertices[i][1];
				v.z = (float)vertices[i][2];
				v.nx = (float)vertex_normals[i][0];
				v.ny = (float)vertex_normals[i][1];
				v.nz = (float)vertex_normals[i][2];
				v.r = (unsigned char)(vertex_colors[i][0] * 255.f);
				v.g = (unsigned char)(vertex_colors[i][1] * 255.f);
				v.b = (unsigned char)(vertex_colors[i][2] * 255.f);
				v.a = 255;
				v.specular_r = (unsigned char)(vertex_specular_colors[i][0] * 255.f);
				v.specular_g = (unsigned char)(vertex_specular_colors[i][1] * 255.f);
				v.specular_b = (unsigned char)(vertex_specular_colors[i][2] * 255.f);
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

		vector<Face>::const_iterator f_it;
		vector<Face>::const_iterator f_end = faces.end();

		for (f_it = faces.begin(); f_it != f_end; ++f_it)
		{
			FaceT face = *f_it;
			f.nverts = (unsigned char)3;

			FaceT::const_iterator fv_it;
			FaceT::const_iterator fv_end = face.end();
			int j = 0;
			for (fv_it = face.begin(); fv_it != fv_end; ++fv_it)
			{
				unsigned int v_idx = *fv_it;
				f.verts[j] = v_idx;
				j++;
			}
			ply_put_element(ply, (void *)&f);
		}

		/* set up and write the face elements */
		ply_put_element_setup(ply, ply::elem_names[2]);
		for (int i = 0; i < n_sh_coefficients(); i++)
		{
			ply::SH_Coefficient sh_c;
			sh_c.value = (float)sh_coefficients[i];
			ply_put_element(ply, (void *)&sh_c);
		}

		/* close the PLY file */
		ply_close(ply);
	}
};
