#include "types.h"

using namespace std;

template<typename T>
void compnorm(const T* ver1, const T* ver2, const T* ver3, T* location,
	bool clockwise = false)
{
	// compute normals assume that the normal at each point
	//   is defined by the triangle consisting of the previous two
	//   points + current point.
	//   (p1-p3) x (p1-p2)

	T norm[3];
	T a[3];
	T b[3];

	if (clockwise)
	{
		a[0] = ver1[0] - ver3[0];
		a[1] = ver1[1] - ver3[1];
		a[2] = ver1[2] - ver3[2];

		b[0] = ver1[0] - ver2[0];
		b[1] = ver1[1] - ver2[1];
		b[2] = ver1[2] - ver2[2];
	}
	else	// Anti-clockwsie
	{
		a[0] = ver1[0] - ver2[0];
		a[1] = ver1[1] - ver2[1];
		a[2] = ver1[2] - ver2[2];

		b[0] = ver1[0] - ver3[0];
		b[1] = ver1[1] - ver3[1];
		b[2] = ver1[2] - ver3[2];
	}

	norm[0] = a[1] * b[2] - a[2] * b[1];
	norm[1] = a[2] * b[0] - a[0] * b[2];
	norm[2] = a[0] * b[1] - a[1] * b[0];

	T normalization = sqrt(norm[1] * norm[1] + norm[2] * norm[2] + norm[0] * norm[0]);
	if (normalization != T(0.0))
	{
		norm[0] /= normalization;
		norm[1] /= normalization;
		norm[2] /= normalization;
	}

	location[0] = norm[0];
	location[1] = norm[1];
	location[2] = norm[2];
}

// Computes vertex normal direction given its position, its one-ring neighbours 
// and the corresponding face indexes. Can handle clockwise and counter-clockwise
template <typename T>
void computeNormal(const T* p, const vector<T*> &adjP,
	const int n_faces, const bool clockwise, T* normal)
{
	normal[0] = T(0.0);
	normal[1] = T(0.0);
	normal[2] = T(0.0);

	for (int i = 0; i < n_faces; i++)
	{
		unsigned int vIdx1 = i;
		unsigned int vIdx2 = (i + 1) % adjP.size();

		T face_normal[3];
		compnorm(p, adjP[vIdx1], adjP[vIdx2], face_normal);

		normal[0] += face_normal[0];
		normal[1] += face_normal[1];
		normal[2] += face_normal[2];
	}

	T norm = sqrt(normal[1] * normal[1] + normal[2] * normal[2] + normal[0] * normal[0]);
	if (norm != T(0.0))
	{
		normal[0] /= norm;
		normal[1] /= norm;
		normal[2] /= norm;
	}
}

// Computes shading value given normal direction, spherical harmonic coefficients 
// and the SH order
template <typename T>
T computeShading(const T* _normal, const T* _sh_coeff,
	const unsigned int _sh_order)
{
	T n_x = _normal[0];
	T n_y = _normal[1];
	T n_z = _normal[2];

	T n_x2 = n_x * n_x;
	T n_y2 = n_y * n_y;
	T n_z2 = n_z * n_z;
	T n_xy = n_x * n_y;
	T n_xz = n_x * n_z;
	T n_yz = n_y * n_z;
	T n_x2_y2 = n_x2 - n_y2;

	T shading = T(_sh_coeff[0]);

	if (_sh_order > 0)
		shading = shading
		+ T(_sh_coeff[1]) * n_x					// x
		+ T(_sh_coeff[2]) * n_y					// y
		+ T(_sh_coeff[3]) * n_z;				// z

	if (_sh_order > 1)
		shading = shading
		+ T(_sh_coeff[4]) * n_xy						// x * y
		+ T(_sh_coeff[5]) * n_xz						// x * z
		+ T(_sh_coeff[6]) * n_yz						// y * z
		+ T(_sh_coeff[7]) * n_x2_y2						// x^2 - y^2
		+ T(_sh_coeff[8]) * (T(3.0) * n_z2 - T(1.0));	// 3 * z^2 - 1

	if (_sh_order > 2)
		shading = shading
		+ T(_sh_coeff[9]) * (T(3.0) * n_x2 - n_y2) * n_y		// (3 * x^2 - y^2) * y 
		+ T(_sh_coeff[10]) * n_x * n_y * n_z					// x * y * z
		+ T(_sh_coeff[11]) * (T(5.0) * n_z2 - T(1.0)) * n_y		// (5 * z^2 - 1) * y
		+ T(_sh_coeff[12]) * (T(5.0) * n_z2 - T(3.0)) * n_z		// (5 * z^2 - 3) * z
		+ T(_sh_coeff[13]) * (T(5.0) * n_z2 - T(1.0)) * n_x		// (5 * z^2 - 1) * x
		+ T(_sh_coeff[14]) * n_x2_y2 * n_z						// (x^2 - y^2) * z
		+ T(_sh_coeff[15]) * (n_x2 - T(3.0) * n_y2) * n_x;		// (x^2 - 3 * y^2) * x

	if (_sh_order > 3)
		shading = shading
		+ T(_sh_coeff[16]) * n_x2_y2 * n_x * n_y								// (x^2 - y^2) * x * y
		+ T(_sh_coeff[17]) * (T(3.0) * n_x2 - n_y2) * n_yz						// (3 * x^2 - y^2) * yz
		+ T(_sh_coeff[18]) * (T(7.0) * n_z2 - T(1.0)) * n_xy					// (7 * z^2 - 1) * x * y
		+ T(_sh_coeff[19]) * (T(7.0) * n_z2 - T(3.0)) * n_yz					// (7 * z^2 - 3) * y * z
		+ T(_sh_coeff[20]) * (T(3.0) - T(30.0) * n_z2 + T(35.0) * n_z2 * n_z2)	// 3 - 30 * z^2 + 35 * z^4
		+ T(_sh_coeff[21]) * (T(7.0) * n_z - T(3.0)) * n_xz						// (7 * z^2 - 3) * x * z
		+ T(_sh_coeff[22]) * (T(7.0) * n_z - T(1.0)) * n_x2_y2					// (7 * z^2 - 1) * (x^2 - y^2)
		+ T(_sh_coeff[23]) * (n_x2 - T(3.0) * n_y2) * n_xz						// (x^2 - 3 * y^2) * x * z
		+ T(_sh_coeff[24]) * ((n_x2 - T(3.0) * n_y2) * n_x2						// (x^2 - 3 * y^2) * x^2 - (3 * x^2 - y^2) * y^2 
		- (T(3.0) * n_x2 - n_y2) * n_y2);

	return shading;
}

// Photometric cost for the case of known albedo and illumination (sh representation)
class ResidualPhotometricError
{
public:
	ResidualPhotometricError(const Intensity* _intensity, 
		const unsigned int _n_channels,
		const unsigned int &_n_adj_vertices, 
		const unsigned int &_n_adj_faces,
		const unsigned int _sh_order,
		const bool _optimize_using_normal = false,
		const bool _use_lower_bound_shading = false,
		const bool _use_upper_bound_shading = false) :
		intensity(_intensity),
		n_channels(_n_channels),
		n_adj_vertices(_n_adj_vertices),
		n_adj_faces(_n_adj_faces),
		sh_order(_sh_order),
		optimize_using_normal(_optimize_using_normal),
		use_lower_bound_shading(_use_lower_bound_shading),
		use_upper_bound_shading(_use_upper_bound_shading)
	{

	}

	template<typename T>
	bool operator()(const T* const* const parameters, T* residuals) const
	{
		const T* albedo = parameters[0];
		const T* sh_coeff = parameters[1];
		const T* lighting_variation = parameters[2];
		const T* vertex = parameters[3];
		T* v = new T[3];
		v[0] = vertex[0];
		v[1] = vertex[1];
		v[2] = vertex[2];

		unsigned int param_idx = 4;

		vector<const T*> adj_vertices;
		adj_vertices.reserve(n_adj_vertices);

		vector<T*> adj_v;
		adj_v.reserve(n_adj_vertices);

		for (size_t i = 0; i < n_adj_vertices; i++)
		{
			const T* v_neighbour = parameters[param_idx];
			adj_vertices.push_back(v_neighbour);

			T* v_n = new T[3];
			v_n[0] = v_neighbour[0];
			v_n[1] = v_neighbour[1];
			v_n[2] = v_neighbour[2];

			adj_v.push_back(v_n);

			param_idx++;
		}

		if (optimize_using_normal)
		{
			const T *vertex_normal = parameters[param_idx];

			param_idx++;

			vector<const T*> adj_normals;
			adj_normals.reserve(n_adj_vertices);
			for (size_t i = 0; i < n_adj_vertices; i++)
			{
				const T* v_neighbour_normal = parameters[param_idx];
				adj_normals.push_back(v_neighbour_normal);

				param_idx++;
			}

			const T *vertex_disp = parameters[param_idx];
			v[0] += vertex_disp[0] * vertex_normal[0];
			v[1] += vertex_disp[0] * vertex_normal[1];
			v[2] += vertex_disp[0] * vertex_normal[2];

			param_idx++;

			for (size_t i = 0; i < n_adj_vertices; i++)
			{
				const T* v_neighbour_disp = parameters[param_idx];

				adj_v[i][0] += v_neighbour_disp[0] * adj_normals[i][0];
				adj_v[i][1] += v_neighbour_disp[0] * adj_normals[i][1];
				adj_v[i][2] += v_neighbour_disp[0] * adj_normals[i][2];

				param_idx++;
			}
		}

		T normal[3];
		computeNormal(&v[0], adj_v, n_adj_faces, false, normal);

		delete[] v;
		for (size_t i = 0; i < n_adj_vertices; i++)
		{
			delete[] adj_v[i];
		}

		T shading = computeShading(normal, sh_coeff, sh_order);

		for (size_t i = 0; i < n_channels; ++i)
		{
			T est_value = albedo[i] * shading + lighting_variation[i];
			residuals[i] = T(intensity[i]) - est_value;
		}

		bool shading_in_bounds = true;
		
		if (use_lower_bound_shading)
		{
			shading_in_bounds = shading_in_bounds && shading >= T(0.0);
		}

		if (use_upper_bound_shading)
		{
			shading_in_bounds = shading_in_bounds && shading <= T(1.0);
		}

		return shading_in_bounds;
	}

private:
	// Vertex intensity
	const Intensity* intensity;

	// Number of color channels
	const unsigned int n_channels;

	// Number of adjacent vertices
	const unsigned int n_adj_vertices;
	// Number of adjacent faces
	const unsigned int n_adj_faces;

	// SH order
	const unsigned int sh_order;

	// Optimize displacement in the direction of the normal
	const bool optimize_using_normal;

	// Constrain lower and/or upper shading bound
	const bool use_lower_bound_shading;
	const bool use_upper_bound_shading;
};

// Photometric cost for the case of known shape and albedo (sh coeff estimation)
class ResidualPhotometricErrorSHCoeff
{
public:
	ResidualPhotometricErrorSHCoeff(const Intensity* _intensity,
		const unsigned int _n_channels, const Intensity* _albedo,
		const Intensity* _lighting_variation, const Normal &_normal, 
		const unsigned int _sh_order,
		const bool _use_lower_bound_shading = false,
		const bool _use_upper_bound_shading = false) :
		intensity(_intensity),
		n_channels(_n_channels),
		albedo(_albedo),
		lighting_variation(_lighting_variation),
		normal(_normal),
		sh_order(_sh_order),
		use_lower_bound_shading(_use_lower_bound_shading),
		use_upper_bound_shading(_use_upper_bound_shading)
	{

	}

	template<typename T>
	bool operator()(const T* const _sh_coeff, T* residuals) const
	{
		T shading = computeShading(normal, _sh_coeff, sh_order);

		for (size_t i = 0; i < n_channels; ++i)
		{
			T est_value = albedo[i] * shading + T(lighting_variation[i]);
			residuals[i] = T(intensity[i]) - est_value;
		}

		bool shading_in_bounds = true;

		if (use_lower_bound_shading)
		{
			shading_in_bounds = shading_in_bounds && shading >= T(0.0);
		}

		if (use_upper_bound_shading)
		{
			shading_in_bounds = shading_in_bounds && shading <= T(1.0);
		}

		return shading_in_bounds;
	}

private:
	// Vertex intensity
	const Intensity* intensity;

	// Number of color channels
	const unsigned int n_channels;

	// Vertex albedo value
	const Intensity* albedo;

	// Vertex lighting variation
	const Intensity* lighting_variation;

	// Vertex normal vector
	const Normal &normal;

	// SH order
	const unsigned int sh_order;

	// Constrain lower and/or upper shading bound
	const bool use_lower_bound_shading;
	const bool use_upper_bound_shading;
};

// Photometric cost for the case of known shading (albedo estimation)
class ResidualPhotometricErrorAlbedo
{
public:
	ResidualPhotometricErrorAlbedo(const Intensity* _intensity,
		const unsigned int _n_channels, const Intensity _shading, 
		const Intensity* _lighting_variation) :
		intensity(_intensity),
		n_channels(_n_channels),
		shading(_shading),
		lighting_variation(_lighting_variation)
	{

	}

	template<typename T>
	bool operator()(const T* const _albedo, T* residuals) const
	{
		for (size_t i = 0; i < n_channels; ++i)
		{
			T est_value = _albedo[i] * T(shading) + lighting_variation[i];
			residuals[i] = T(intensity[i]) - est_value;
		}

		return true;
	}

private:
	// Vertex intensity
	const Intensity* intensity;

	// Number of color channels
	const unsigned int n_channels;

	// Vertex shading value
	const Intensity shading;

	// Vertex lighting variation
	const Intensity* lighting_variation;
};

// Photometric cost for the case of known albedo and shading 
// (lighting variation estimation)
class ResidualPhotometricErrorLightingVariation
{
public:
	ResidualPhotometricErrorLightingVariation(const Intensity* _intensity,
		const unsigned int _n_channels, const Intensity* _diffuse_intensity) :
		intensity(_intensity),
		n_channels(_n_channels),
		diffuse_intensity(_diffuse_intensity)
	{

	}

	template<typename T>
	bool operator()(const T* const _lighting_variation, T* residuals) const
	{
		for (size_t i = 0; i < n_channels; ++i)
		{
			T est_value = T(diffuse_intensity[i]) + _lighting_variation[i];
			residuals[i] = T(intensity[i]) - est_value;
		}

		return true;
	}

private:
	// Vertex intensity
	const Intensity* intensity;

	// Number of color channels
	const unsigned int n_channels;

	// Vertex estimated diffuse component
	const Intensity* diffuse_intensity;
};

// Photometric cost for the case of known shading 
// (albedo and lighting variation estimation)
class ResidualPhotometricErrorAlbedoLightingVariation
{
public:
	ResidualPhotometricErrorAlbedoLightingVariation(const Intensity* _intensity,
		const unsigned int _n_channels, const Intensity _shading) :
		intensity(_intensity),
		n_channels(_n_channels),
		shading(_shading)
	{

	}

	template<typename T>
	bool operator()(const T* const _albedo, 
		const T* const _lighting_variation, T* residuals) const
	{
		for (size_t i = 0; i < n_channels; ++i)
		{
			T est_value = _albedo[i] * T(shading) + _lighting_variation[i];
			residuals[i] = T(intensity[i]) - est_value;
		}

		return true;
	}

private:
	// Vertex intensity
	const Intensity* intensity;

	// Number of color channels
	const unsigned int n_channels;

	// Vertex shading
	const Intensity shading;
};

// Photometric cost for the case of known albedo and shading
// (shape estimation)
class ResidualPhotometricErrorShape
{
public:
	ResidualPhotometricErrorShape(const Intensity* _intensity,
		const unsigned int _n_channels,
		const Intensity* _albedo, const Intensity* _lighting_variation,
		const unsigned int &_n_adj_vertices,
		const unsigned int &_n_adj_faces,
		const Intensity *_sh_coeff,
		const unsigned int _sh_order,
		const bool _optimize_using_normal = false,
		const bool _use_lower_bound_shading = false,
		const bool _use_upper_bound_shading = false) :
		intensity(_intensity),
		n_channels(_n_channels),
		albedo(_albedo),
		lighting_variation(_lighting_variation),
		n_adj_vertices(_n_adj_vertices),
		n_adj_faces(_n_adj_faces),
		sh_coeff(_sh_coeff),
		sh_order(_sh_order),
		optimize_using_normal(_optimize_using_normal),
		use_lower_bound_shading(_use_lower_bound_shading),
		use_upper_bound_shading(_use_upper_bound_shading)
	{

	}

	template<typename T>
	bool operator()(const T* const* const parameters, T* residuals) const
	{
		const T* vertex = parameters[0];
		T* v = new T[3];
		v[0] = vertex[0];
		v[1] = vertex[1];
		v[2] = vertex[2];

		unsigned int param_idx = 1;

		vector<const T*> adj_vertices;
		adj_vertices.reserve(n_adj_vertices);

		vector<T*> adj_v;
		adj_v.reserve(n_adj_vertices);

		for (size_t i = 0; i < n_adj_vertices; i++)
		{
			const T* v_neighbour = parameters[param_idx];
			adj_vertices.push_back(v_neighbour);

			T* v_n = new T[3];
			v_n[0] = v_neighbour[0];
			v_n[1] = v_neighbour[1];
			v_n[2] = v_neighbour[2];

			adj_v.push_back(v_n);

			param_idx++;
		}

		if (optimize_using_normal)
		{
			const T *vertex_normal = parameters[param_idx];

			param_idx++;

			vector<const T*> adj_normals;
			adj_normals.reserve(n_adj_vertices);
			for (size_t i = 0; i < n_adj_vertices; i++)
			{
				const T* v_neighbour_normal = parameters[param_idx];
				adj_normals.push_back(v_neighbour_normal);

				param_idx++;
			}

			const T *vertex_disp = parameters[param_idx];
			v[0] += vertex_disp[0] * vertex_normal[0];
			v[1] += vertex_disp[0] * vertex_normal[1];
			v[2] += vertex_disp[0] * vertex_normal[2];

			param_idx++;

			for (size_t i = 0; i < n_adj_vertices; i++)
			{
				const T* v_neighbour_disp = parameters[param_idx];

				adj_v[i][0] += v_neighbour_disp[0] * adj_normals[i][0];
				adj_v[i][1] += v_neighbour_disp[0] * adj_normals[i][1];
				adj_v[i][2] += v_neighbour_disp[0] * adj_normals[i][2];

				param_idx++;
			}
		}

		T normal[3];
		computeNormal(&v[0], adj_v, n_adj_faces, false, normal);

		delete[] v;
		for (size_t i = 0; i < n_adj_vertices; i++)
		{
			delete[] adj_v[i];
		}

		T shading = computeShading(normal, sh_coeff, sh_order);

		for (size_t i = 0; i < n_channels; ++i)
		{
			T est_value = T(albedo[i]) * shading + T(lighting_variation[i]);
			residuals[i] = T(intensity[i]) - est_value;
		}

		bool shading_in_bounds = true;

		if (use_lower_bound_shading)
		{
			shading_in_bounds = shading_in_bounds && shading >= T(0.0);
		}

		if (use_upper_bound_shading)
		{
			shading_in_bounds = shading_in_bounds && shading <= T(1.0);
		}

		return shading_in_bounds;
	}

private:
	// Vertex intensity
	const Intensity* intensity;

	// Number of color channels
	const unsigned int n_channels;

	// Vertex albedo value
	const Intensity* albedo;

	// Vertex lighting variation
	const Intensity* lighting_variation;

	// Number of adjacent vertices
	const unsigned int n_adj_vertices;
	// Number of adjacent faces
	const unsigned int n_adj_faces;

	// SH coefficients
	const Intensity *sh_coeff;

	// SH order
	const unsigned int sh_order;

	// Optimize displacement in the direction of the normal
	const bool optimize_using_normal;

	// Constrain lower and/or upper shading bound
	const bool use_lower_bound_shading;
	const bool use_upper_bound_shading;
};

// Residual of difference between values of neighbour vertices
class ResidualWeightedDifference
{
public:
	ResidualWeightedDifference(
		const Intensity _weight,
		const unsigned int _n_channels) :
		weight(_weight),
		n_channels(_n_channels)
	{

	}

	template <typename T>
	bool operator()(const T* const _value1, const T* const _value2,
		T* residuals) const
	{
		for (size_t i = 0; i < n_channels; ++i)
		{
			residuals[i] = T(weight) * (_value1[i] - _value2[i]);
		}

		return true;
	}

private:
	// Weight
	const Intensity weight;

	// Number of color channels
	const unsigned int n_channels;
};

// Residual of magnitude of value
class ResidualValueMagnitude
{
public:
	ResidualValueMagnitude(const unsigned int _n_channels) :
		n_channels(_n_channels)
	{

	}

	template <typename T>
	bool operator()(const T* const _value,
		T* residuals) const
	{
		for (size_t i = 0; i < n_channels; ++i)
		{
			residuals[i] = _value[i];
		}

		return true;
	}

private:
	// Number of color channels
	const unsigned int n_channels;
};

// Residual of difference with initial value
class ResidualInitialVertexDiff
{
public:
	ResidualInitialVertexDiff(const Vertex &_vertex0) :
		vertex0(_vertex0)
	{

	}

	template <typename T>
	bool operator()(const T* const _vertex,
		T* residuals) const
	{
		for (int i = 0; i < 3; ++i)
		{
			residuals[i] = _vertex[i] - T(vertex0[i]);
		}

		return true;
	}

	template <typename T>
	bool operator()(const T* const _vertex, const T* const _vertex_disp,
		T* residuals) const
	{
		residuals[0] = _vertex_disp[0];

		return true;
	}

private:
	// Initial vertex
	const Vertex &vertex0;
};

// Laplacian smoothing residual
class ResidualLaplacianSmoothing
{
public:
public:
	ResidualLaplacianSmoothing(
		const unsigned int _n_adj_vertices, 
		const bool _optimize_using_normal = false) :
		n_adj_vertices(_n_adj_vertices), 
		optimize_using_normal(_optimize_using_normal)
	{

	}

	template <typename T>
	bool operator()(const T* const* const parameters, T* residuals) const
	{
		// Parameters:
		// 0 - Current vertex position or translation
		// >0 - Neighbour vertices positions or translations

		const T* vertex = parameters[0];
		T* p = new T[3];
		p[0] = vertex[0];
		p[1] = vertex[1];
		p[2] = vertex[2];

		vector<const T*> adj_vertices;
		adj_vertices.reserve(n_adj_vertices);

		vector<T*> adjP;
		adjP.reserve(n_adj_vertices);
		for (size_t i = 0; i < n_adj_vertices; i++)
		{
			const T* p_neighbour = parameters[i + 1];
			adj_vertices.push_back(p_neighbour);

			T* v_neigh = new T[3];
			v_neigh[0] = p_neighbour[0];
			v_neigh[1] = p_neighbour[1];
			v_neigh[2] = p_neighbour[2];
			adjP.push_back(v_neigh);
		}

		if (optimize_using_normal)
		{
			unsigned int param_idx = 1 + n_adj_vertices;

			const T *vertex_normal = parameters[param_idx];

			param_idx++;

			vector<const T*> adj_normals;
			adj_normals.reserve(n_adj_vertices);
			for (size_t i = 0; i < n_adj_vertices; i++)
			{
				const T* v_neighbour_normal = parameters[param_idx];
				adj_normals.push_back(v_neighbour_normal);

				param_idx++;
			}

			const T *vertex_disp = parameters[param_idx];
			p[0] += vertex_disp[0] * vertex_normal[0];
			p[1] += vertex_disp[0] * vertex_normal[1];
			p[2] += vertex_disp[0] * vertex_normal[2];

			param_idx++;

			for (size_t i = 0; i < n_adj_vertices; i++)
			{
				const T* v_neighbour_disp = parameters[param_idx];

				adjP[i][0] += v_neighbour_disp[0] * adj_normals[i][0];
				adjP[i][1] += v_neighbour_disp[0] * adj_normals[i][1];
				adjP[i][2] += v_neighbour_disp[0] * adj_normals[i][2];

				param_idx++;
			}
		}

		T laplacian_x = T(0.0);
		T laplacian_y = T(0.0);
		T laplacian_z = T(0.0);
		T weight_norm = T(0.0);

		for (size_t i = 1; i <= n_adj_vertices; i++)
		{
			T weight = T(0.0);
			T area = T(0.0);

			int prev_i = (i - 1) % n_adj_vertices;
			int curr_i = i % n_adj_vertices;
			int next_i = (i + 1) % n_adj_vertices;

			weight += computeCotangent(p, adjP[curr_i], adjP[prev_i]);
			weight += computeCotangent(p, adjP[curr_i], adjP[next_i]);

			laplacian_x += weight * (adjP[curr_i][0] - p[0]);
			laplacian_y += weight * (adjP[curr_i][1] - p[1]);
			laplacian_z += weight * (adjP[curr_i][2] - p[2]);

			weight_norm += weight;
		}

		if (weight_norm < T(std::numeric_limits<double>::epsilon()) && weight_norm > T(-std::numeric_limits<double>::epsilon()))
		{
			weight_norm = T(std::numeric_limits<double>::epsilon());
		}

		laplacian_x /= weight_norm;
		laplacian_y /= weight_norm;
		laplacian_z /= weight_norm;

		delete[] p;
		for (size_t i = 0; i < n_adj_vertices; i++)
		{
			delete[] adjP[i];
		}

		residuals[0] = ceres::sqrt(laplacian_x * laplacian_x
			+ laplacian_y * laplacian_y
			+ laplacian_z * laplacian_z);

		return true;
	}

	template <typename T>
	T computeCotangent(const T* const _p0, const T* const _p1,
		const T* const _p2) const
	{
		// Vector from the thrid vertex to the center of the ring
		T v1_x = _p0[0] - _p2[0];
		T v1_y = _p0[1] - _p2[1];
		T v1_z = _p0[2] - _p2[2];

		// Normalize first vector
		T v1_norm = ceres::sqrt(v1_x * v1_x + v1_y * v1_y + v1_z * v1_z);

		if (v1_norm < T(std::numeric_limits<double>::epsilon()) && v1_norm > T(-std::numeric_limits<double>::epsilon()))
		{
			v1_norm = T(std::numeric_limits<double>::epsilon());
		}

		v1_x /= v1_norm;
		v1_y /= v1_norm;
		v1_z /= v1_norm;

		// Vector from the third vertex to the current neighbour
		T v2_x = _p1[0] - _p2[0];
		T v2_y = _p1[1] - _p2[1];
		T v2_z = _p1[2] - _p2[2];

		// Normalize second vector
		T v2_norm = ceres::sqrt(v2_x * v2_x + v2_y * v2_y + v2_z * v2_z);

		if (v2_norm < T(std::numeric_limits<double>::epsilon()) && v2_norm > T(-std::numeric_limits<double>::epsilon()))
		{
			v2_norm = T(std::numeric_limits<double>::epsilon());
		}

		v2_x /= v2_norm;
		v2_y /= v2_norm;
		v2_z /= v2_norm;

		// Cosine of alpha as the dot product of both vectors
		T cos_a = v1_x * v2_x + v1_y * v2_y + v1_z * v2_z;

		// Projection of the first vector onto the orthogonal component of the 
		// second vector on the plane v1-v2
		T v1p_x = v1_x - cos_a * v2_x;
		T v1p_y = v1_y - cos_a * v2_y;
		T v1p_z = v1_z - cos_a * v2_z;

		// Sine of alpha as the norm of the projcted vector
		T sin_a = ceres::sqrt(v1p_x * v1p_x + v1p_y * v1p_y + v1p_z * v1p_z);

		//if (sin_a < T(1e-20) && sin_a > T(-1e-20))
		//{
		//	std::cout << "ERROR!!!!! sin(a) == 0" << endl;
		//}

		if (sin_a < T(std::numeric_limits<double>::epsilon()) && sin_a > T(-std::numeric_limits<double>::epsilon()))
		{
			sin_a = T(std::numeric_limits<double>::epsilon());
		}

		// cot(a) = cos(a) / sin(a)
		T cot_a = cos_a / sin_a;

		return cot_a;
	}

private:
	// Number of adjacent vertices
	const unsigned int n_adj_vertices;

	// Optimize displacement in the direction of the normal
	const bool optimize_using_normal;
};

// Residual Total Variation
class ResidualTV
{
public:
	ResidualTV(const Vertex &_vertex1, const Vertex &_vertex2) :
		vertex1(_vertex1),
		vertex2(_vertex2)
	{

	}

	template <typename T>
	bool operator()(const T* const _vertex1, const T* const _vertex2,
		T* residuals) const
	{
		for (size_t i = 0; i < 3; ++i)
		{
			residuals[i] = (_vertex1[i] - _vertex2[i]) 
				- (T(vertex1[i]) - T(vertex2[i]));
		}

		return true;
	}

	template <typename T>
	bool operator()(const T* const _vertex1, const T* const _vertex2,
		const T* const _normal1, const T* const _normal2,
		const T* const _vertex_disp1, const T* const _vertex_disp2,
		T* residuals) const
	{
		for (size_t i = 0; i < 3; ++i)
		{
			residuals[i] = (_vertex1[i] + _vertex_disp1[0] * _normal1[i] 
				- _vertex2[i] - _vertex_disp2[0] * _normal2[i])
				- (T(vertex1[i]) - T(vertex2[i]));
		}

		return true;
	}


private:
	const Vertex &vertex1, &vertex2;
};

// Photometric cost for the case of known normal
class ResidualPhotometricErrorWithNormal
{
public:
	ResidualPhotometricErrorWithNormal(const Intensity* _intensity,
		const Coordinate* _normal,
		const unsigned int _n_channels,
		const unsigned int _sh_order,
		const bool _use_lower_bound_shading = false,
		const bool _use_upper_bound_shading = false) :
		intensity(_intensity),
		normal(_normal),
		n_channels(_n_channels),
		sh_order(_sh_order),
		use_lower_bound_shading(_use_lower_bound_shading),
		use_upper_bound_shading(_use_upper_bound_shading)
	{

	}

	template<typename T>
	bool operator()(const T* const* const parameters, T* residuals) const
	{
		const T* albedo = parameters[0];
		const T* sh_coeff = parameters[1];
		const T* lighting_variation = parameters[2];

		T aux_normal[3];
		aux_normal[0] = T(normal[0]);
		aux_normal[1] = T(normal[1]);
		aux_normal[2] = T(normal[2]);
		T shading = computeShading(aux_normal, sh_coeff, sh_order);

		for (size_t i = 0; i < n_channels; ++i)
		{
			T est_value = albedo[i] * shading + lighting_variation[i];
			residuals[i] = T(intensity[i]) - est_value;
		}

		bool shading_in_bounds = true;

		if (use_lower_bound_shading)
		{
			shading_in_bounds = shading_in_bounds && shading >= T(0.0);
		}

		if (use_upper_bound_shading)
		{
			shading_in_bounds = shading_in_bounds && shading <= T(1.0);
		}

		return shading_in_bounds;
	}

private:
	// Vertex intensity
	const Intensity* intensity;

	// Vertex normal vector
	const Coordinate* normal;

	// Number of color channels
	const unsigned int n_channels;

	// SH order
	const unsigned int sh_order;

	// Constrain lower and/or upper shading bound
	const bool use_lower_bound_shading;
	const bool use_upper_bound_shading;
};

// Photometric cost for the case of known shading
class ResidualPhotometricErrorWithShading
{
public:
	ResidualPhotometricErrorWithShading(const Intensity* _intensity,
		const Intensity _shading,
		const unsigned int _n_channels) :
		intensity(_intensity),
		shading(_shading),
		n_channels(_n_channels)
	{

	}

	template<typename T>
	bool operator()(const T* const* const parameters, T* residuals) const
	{
		const T* albedo = parameters[0];
		const T* lighting_variation = parameters[1];

		for (size_t i = 0; i < n_channels; ++i)
		{
			T est_value = albedo[i] * T(shading) + lighting_variation[i];
			residuals[i] = T(intensity[i]) - est_value;
		}

		return true;
	}

private:
	// Vertex intensity
	const Intensity* intensity;

	// Vertex shading
	const Intensity shading;

	// Number of color channels
	const unsigned int n_channels;
};

// Photometric cost for the case of known shading
class ResidualPhotometricErrorWithAlbedoShading
{
public:
	ResidualPhotometricErrorWithAlbedoShading(const Intensity* _intensity,
		const Intensity* _albedo,
		const Intensity _shading,
		const unsigned int _n_channels) :
		intensity(_intensity),
		albedo(_albedo),
		shading(_shading),
		n_channels(_n_channels)
	{

	}

	template<typename T>
	bool operator()(const T* const* const parameters, T* residuals) const
	{
		const T* lighting_variation = parameters[0];

		for (size_t i = 0; i < n_channels; ++i)
		{
			T est_value = T(albedo[i]) * T(shading) + lighting_variation[i];
			residuals[i] = T(intensity[i]) - est_value;
		}

		return true;
	}

private:
	// Vertex intensity
	const Intensity* intensity;

	// Vertex albedo
	const Intensity* albedo;

	// Vertex shading
	const Intensity shading;

	// Number of color channels
	const unsigned int n_channels;
};