#ifndef PARAMETERS_H_
#define PARAMETERS_H_

// I/O parameters
#define MODE				"mode"
#define ESTIMATE_SH_COEFF	0
#define DECOMPOSE_INTRINSIC	1

#define INPUT_MESH_FILENAME				"input_mesh_filename"
#define INPUT_SH_COEFF_FILENAME			"input_sh_coeff_filename"

#define OUTPUT_ALBEDO_MESH_FILENAME		"output_albedo_mesh_filename"
#define OUTPUT_SHADING_MESH_FILENAME	"output_shading_mesh_filename"
#define OUTPUT_SH_COEFF_FILENAME		"output_sh_coeff_filename"
#define OUTPUT_INTENSITY_MESH_FILENAME	"output_intensity_mesh_filename"

#define SH_ORDER						"sh_order"

#define MIN_GRAY						"min_gray"
#define MAX_GRAY						"max_gray"

#define MAX_SHADING			"max_shading"

#define SAVE_SHADING		"save_shading"
#define SAVE_ALBEDO			"save_albedo"
#define SAVE_EST_INTENSITY	"save_est_intensity"
#define SAVE_SH_COEFF		"save_sh_coeff"
#define SAVE_MESH_BINARY	"save_mesh_binary"

#define USE_MASK_VALUES		"use_mask_values"

#include <string>
#include <vector>
#include "opencv2/core/core.hpp"

using namespace std;

namespace parameters {
	struct Parameters {

		// Default constructor that sets up a generic sparse problem.
		Parameters() {
			//input_mesh_filename = "C:/Users/Qi/Desktop/generated/meshes/cropped/mesh_0180.ply";

			mode = ESTIMATE_SH_COEFF;

			input_mesh_filename = "C:/Users/Qi/Desktop/generated/mesh_0180_albedo_0001.ply";

			output_albedo_mesh_filename = "C:/Users/Qi/Desktop/generated/mesh_0180_albedo.ply";
			output_shading_mesh_filename = "C:/Users/Qi/Desktop/generated/mesh_0180_shading.ply";
			output_intensity_mesh_filename = "C:/Users/Qi/Desktop/generated/mesh_0180_intensity.ply";
			output_sh_coeff_filename = "C:/Users/Qi/Desktop/generated/sh_coeff.txt";
			sh_order = 2;

			//min_gray = 0.3f;
			//max_gray = 0.598f;

			min_gray = 0.4f;
			max_gray = 0.7f;

			max_shading = 1.0f;

			save_shading = true;
			save_albedo = true;
			save_est_intensity = true;
			save_sh_coeff = true;
			save_mesh_binary = true;

			use_mask_values = true;
		}


		int mode;

		string input_mesh_filename;
		string input_sh_coeff_filename;
		string output_albedo_mesh_filename;
		string output_shading_mesh_filename;
		string output_intensity_mesh_filename;
		string output_sh_coeff_filename;
		int sh_order;
		float min_gray, max_gray;

		float max_shading;
		
		bool save_shading;
		bool save_albedo;
		bool save_est_intensity;
		bool save_sh_coeff;
		bool save_mesh_binary;

		bool use_mask_values;

		// Load parameters from an ini file
		inline void Parameters::load(const std::string &_filename)
		{
			cv::FileStorage fs(_filename, cv::FileStorage::READ);

			// I/O parameters
			if (!fs[MODE].empty())
			{
				fs[MODE] >> mode;
			}

			if (!fs[INPUT_MESH_FILENAME].empty())
			{
				fs[INPUT_MESH_FILENAME] >> input_mesh_filename;
			}

			if (!fs[INPUT_SH_COEFF_FILENAME].empty())
			{
				fs[INPUT_SH_COEFF_FILENAME] >> input_sh_coeff_filename;
			}

			if (!fs[OUTPUT_ALBEDO_MESH_FILENAME].empty())
			{
				fs[OUTPUT_ALBEDO_MESH_FILENAME] >> output_albedo_mesh_filename;
			}

			if (!fs[OUTPUT_SHADING_MESH_FILENAME].empty())
			{
				fs[OUTPUT_SHADING_MESH_FILENAME] >> output_shading_mesh_filename;
			}

			if (!fs[OUTPUT_INTENSITY_MESH_FILENAME].empty())
			{
				fs[OUTPUT_INTENSITY_MESH_FILENAME] >> output_intensity_mesh_filename;
			}

			if (!fs[OUTPUT_SH_COEFF_FILENAME].empty())
			{
				fs[OUTPUT_SH_COEFF_FILENAME] >> output_sh_coeff_filename;
			}

			if (!fs[SH_ORDER].empty())
			{
				fs[SH_ORDER] >> sh_order;
			}

			if (!fs[MIN_GRAY].empty())
			{
				fs[MIN_GRAY] >> min_gray;
			}

			if (!fs[MAX_GRAY].empty())
			{
				fs[MAX_GRAY] >> max_gray;
			}

			if (!fs[SAVE_SHADING].empty())
			{
				fs[SAVE_SHADING] >> save_shading;
			}

			if (!fs[SAVE_ALBEDO].empty())
			{
				fs[SAVE_ALBEDO] >> save_albedo;
			}

			if (!fs[SAVE_EST_INTENSITY].empty())
			{
				fs[SAVE_EST_INTENSITY] >> save_est_intensity;
			}

			if (!fs[SAVE_SH_COEFF].empty())
			{
				fs[SAVE_SH_COEFF] >> save_sh_coeff;
			}

			if (!fs[MAX_SHADING].empty())
			{
				fs[MAX_SHADING] >> max_shading;
			}

			if (!fs[SAVE_MESH_BINARY].empty())
			{
				fs[SAVE_MESH_BINARY] >> save_mesh_binary;
			}

			if (!fs[USE_MASK_VALUES].empty())
			{
				fs[USE_MASK_VALUES] >> use_mask_values;
			}
		}

		inline void Parameters::save(const std::string &_filename)
		{

		}
	};
}

#endif  // PARAMETERS_H_
