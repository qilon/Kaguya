#ifndef PARAMETERS_H_
#define PARAMETERS_H_

// Data terms
#define DATA_WEIGHTS		"data_weights"
#define DATA_HUBER_WIDTH	"data_huber_width"

// Lighting confidence
#define SPECULAR_WEIGHT_VAR	"specular_weight_var"

// SH Coefficients estimation
#define SH_ORDER				"sh_order"
#define ALBEDO_PERCENTILE		"albedo_percentile"
#define USE_LOWER_BOUND_SHADING	"use_lower_bound_shading"
#define USE_UPPER_BOUND_SHADING	"use_upper_bound_shading"

// Smoothing Weight
#define COLOR_DIFF_THRESHOLD	"color_diff_threshold"
#define COLOR_DIFF_VAR			"color_diff_var"
#define DEPTH_DIFF_VAR			"depth_diff_var"
#define USE_DEPTH_WEIGHT		"use_depth_weight"
#define SMOOTH_SPECULAR_WEIGHT_VAR	"smooth_specular_weight_var"

// Albedo estimation
#define SMOOTH_ALBEDO_WEIGHT		"smooth_albedo_weight"
#define SMOOTH_ALBEDO_HUBER_WIDTH	"smooth_albedo_huber_width"
#define USE_LOWER_BOUND_ALBEDO		"use_lower_bound_albedo"
#define USE_UPPER_BOUND_ALBEDO		"use_upper_bound_albedo"
#define ALBEDO_INITIALIZATION		"albedo_initialization"
#define WHITE_ALBEDO	0
#define COLOR_ALBEDO	1
#define EST_ALBEDO		2

// Local lighting variation estimation
#define SMOOTH_LOCAL_LIGHTING_WEIGHT		"smooth_local_lighting_weight"
#define LOCAL_LIGHTING_MAGNITUDE_WEIGHT		"local_lighting_magnitude_weight"
#define SMOOTH_LOCAL_LIGHTING_HUBER_WIDTH	"smooth_local_lighting_huber_width"
#define LOCAL_LIGHTING_MAGNITUDE_HUBER_WIDTH "local_lighting_magnitude_huber_width"
#define USE_LOWER_BOUND_LIGHTING			"use_lower_bound_lighting"
#define USE_UPPER_BOUND_LIGHTING			"use_upper_bound_lighting"

// Join estimation of albedo and local lighting variation
#define COMBINE_ALBEDO_LIGHTING "combine_albedo_lighting"

// Shape from shading
#define ESTIMATE_SHAPE			"estimate_shape"
#define DISPLACEMENT_WEIGHT		"displacement_weight"
#define LAPLACIAN_SMOOTH_WEIGHT	"laplacian_smooth_weight"
#define TV_WEIGHT				"tv_weight"
#define USE_NORMAL_FOR_SHAPE	"use_normal_for_shape"

// Solution refinement
#define REFINE_SOLUTION	"refine_solution"

// I/O parameters
#define INPUT_MESH_FILENAME				"input_mesh_filename"
#define INPUT_SPECULAR_IMAGE_FILENAME	"input_specular_image_filename"
#define INPUT_INTRINSICS_FILENAME		"input_intrinsics_filename"

#define OUTPUT_SH_COEFF_FILENAME			"output_sh_coeff_filename"
#define OUTPUT_ALBEDO_MESH_FILENAME			"output_albedo_mesh_filename"
#define OUTPUT_SHADING_MESH_FILENAME		"output_shading_mesh_filename"
#define OUTPUT_LOCAL_LIGHTING_MESH_FILENAME	"output_local_lighting_mesh_filename"
#define OUTPUT_EST_DIFFUSE_MESH_FILENAME	"output_est_diffuse_mesh_filename"
#define OUTPUT_DIFFUSE_DIFF_MESH_FILENAME	"output_diffuse_diff_mesh_filename"
#define OUTPUT_EST_INTENSITY_MESH_FILENAME	"output_est_intensity_mesh_filename"
#define OUTPUT_ORIG_INTENSITY_MESH_FILENAME	"output_orig_intensity_mesh_filename"

#define SAVE_SH_COEFF		"save_sh_coeff"
#define SAVE_ALBEDO			"save_albedo"
#define SAVE_SHADING		"save_shading"
#define SAVE_LOCAL_LIGHTING	"save_local_lighting"
#define SAVE_EST_DIFFUSE	"save_est_diffuse"
#define SAVE_DIFFUSE_DIFF	"save_diffuse_diff"
#define SAVE_EST_INTENSITY	"save_est_intensity"
#define SAVE_ORIG_INTENSITY	"save_orig_intensity"
#define SAVE_MESH_BINARY	"save_mesh_binary"

// Ceres options root
#define CERES_SOLVER "ceres_solver"

// Ceres options
#define MINIMIZER_TYPE "minimizer_type"
#define LINE_SEARCH_DIRECTION_TYPE "line_search_direction_type"
#define LINE_SEARCH_TYPE "line_search_type"
#define NONLINEAR_CONJUGATE_GRADIENT_TYPE "nonlinear_conjugate_gradient_type"
#define MAX_LBFGS_RANK "max_lbfgs_rank"
#define USE_APPROXIMATE_EIGENVALUE_BFGS_SCALING "use_approximate_eigenvalue_bfgs_scaling"
#define LINE_SEARCH_INTERPOLATION_TYPE "line_search_interpolation_type"
#define MIN_LINE_SEARCH_STEP_SIZE "min_line_search_step_size"
#define LINE_SEARCH_SUFFICIENT_FUNCTION_DECREASE "line_search_sufficient_function_decrease"
#define MAX_LINE_SEARCH_STEP_CONTRACTION "max_line_search_step_contraction"
#define MIN_LINE_SEARCH_STEP_CONTRACTION "min_line_search_step_contraction"
#define MAX_NUM_LINE_SEARCH_STEP_SIZE_ITERATIONS "max_num_line_search_step_size_iterations"
#define MAX_NUM_LINE_SEARCH_DIRECTION_RESTARTS "max_num_line_search_direction_restarts"
#define LINE_SEARCH_SUFFICIENT_CURVATURE_DECREASE "line_search_sufficient_curvature_decrease"
#define MAX_LINE_SEARCH_STEP_EXPANSION "max_line_search_step_expansion"
#define TRUST_REGION_STRATEGY_TYPE "trust_region_strategy_type"
#define DOGLEG_TYPE "dogleg_type"
#define USE_NONMONOTONIC_STEPS "use_nonmonotonic_steps"
#define MAX_CONSECUTIVE_NONMONOTONIC_STEPS "max_consecutive_nonmonotonic_steps"
#define MAX_NUM_ITERATIONS "max_num_iterations"
#define MAX_SOLVER_TIME_IN_SECONDS "max_solver_time_in_seconds"
#define NUM_THREADS "num_threads"
#define INITIAL_TRUST_REGION_RADIUS "initial_trust_region_radius"
#define MAX_TRUST_REGION_RADIUS "max_trust_region_radius"
#define MIN_TRUST_REGION_RADIUS "min_trust_region_radius"
#define MIN_RELATIVE_DECREASE "min_relative_decrease"
#define MIN_LM_DIAGONAL "min_lm_diagonal"
#define MAX_LM_DIAGONAL "max_lm_diagonal"
#define MAX_NUM_CONSECUTIVE_INVALID_STEPS "max_num_consecutive_invalid_steps"
#define FUNCTION_TOLERANCE "function_tolerance"
#define GRADIENT_TOLERANCE "gradient_tolerance"
#define PARAMETER_TOLERANCE "parameter_tolerance"
#define LINEAR_SOLVER_TYPE "linear_solver_type"
#define PRECONDITIONER_TYPE "preconditioner_type"
#define VISIBILITY_CLUSTERING_TYPE "visibility_clustering_type"
#define DENSE_LINEAR_ALGEBRA_LIBRARY_TYPE "dense_linear_algebra_library_type"
#define SPARSE_LINEAR_ALGEBRA_LIBRARY_TYPE "sparse_linear_algebra_library_type"
#define NUM_LINEAR_SOLVER_THREADS "num_linear_solver_threads"
#define USE_EXPLICIT_SCHUR_COMPLEMENT "use_explicit_schur_complement"
#define USE_POSTORDERING "use_postordering"
#define DYNAMIC_SPARSITY "dynamic_sparsity"
#define MIN_LINEAR_SOLVER_ITERATIONS "min_linear_solver_iterations"
#define MAX_LINEAR_SOLVER_ITERATIONS "max_linear_solver_iterations"
#define ETA "eta"
#define JACOBI_SCALING "jacobi_scaling"
#define USE_INNER_ITERATIONS "use_inner_iterations"
#define INNER_ITERATION_TOLERANCE "inner_iteration_tolerance"
#define LOGGING_TYPE "logging_type"
#define MINIMIZER_PROGRESS_TO_STDOUT "minimizer_progress_to_stdout"
#define TRUST_REGION_PROBLEM_DUMP_DIRECTORY "trust_region_problem_dump_directory"
#define TRUST_REGION_PROBLEM_DUMP_FORMAT_TYPE "trust_region_problem_dump_format_type"
#define CHECK_GRADIENTS "check_gradients"
#define GRADIENT_CHECK_RELATIVE_PRECISION "gradient_check_relative_precision"
#define NUMERIC_DERIVATIVE_RELATIVE_STEP_SIZE "numeric_derivative_relative_step_size"
#define UPDATE_STATE_EVERY_ITERATION "update_state_every_iteration"


#define NUM_OPTIMIZATIONS	4

#include <string>
#include <vector>
#include "opencv2/core/core.hpp"
#include <ceres/ceres.h>

using namespace std;

	struct SolverOptions
	{
		SolverOptions()
		{
			options = vector<ceres::Solver::Options>(NUM_OPTIMIZATIONS);
		}

		vector<ceres::Solver::Options> options;

		void read(const cv::FileNode& node)
		{
			if (!node[MINIMIZER_TYPE].empty())
			{
				vector<int> aux_vector;
				node[MINIMIZER_TYPE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].minimizer_type 
						= ceres::MinimizerType(aux_vector[i]);
				}
			}

			if (!node[LINE_SEARCH_DIRECTION_TYPE].empty())
			{
				vector<int> aux_vector;
				node[LINE_SEARCH_DIRECTION_TYPE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].line_search_direction_type 
						= ceres::LineSearchDirectionType(aux_vector[i]);
				}
			}

			if (!node[LINE_SEARCH_TYPE].empty())
			{
				vector<int> aux_vector;
				node[LINE_SEARCH_TYPE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].line_search_type 
						= ceres::LineSearchType(aux_vector[i]);
				}
			}

			if (!node[NONLINEAR_CONJUGATE_GRADIENT_TYPE].empty())
			{
				vector<int> aux_vector;
				node[NONLINEAR_CONJUGATE_GRADIENT_TYPE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].nonlinear_conjugate_gradient_type 
						= ceres::NonlinearConjugateGradientType(aux_vector[i]);
				}
			}

			if (!node[MAX_LBFGS_RANK].empty())
			{
				vector<int> aux_vector;
				node[MAX_LBFGS_RANK] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].max_lbfgs_rank = aux_vector[i];
				}
			}

			if (!node[USE_APPROXIMATE_EIGENVALUE_BFGS_SCALING].empty())
			{
				vector<int> aux_vector;
				node[USE_APPROXIMATE_EIGENVALUE_BFGS_SCALING] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].use_approximate_eigenvalue_bfgs_scaling 
						= aux_vector[i] != 0;
				}
			}

			if (!node[LINE_SEARCH_INTERPOLATION_TYPE].empty())
			{
				vector<int> aux_vector;
				node[LINE_SEARCH_INTERPOLATION_TYPE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].line_search_interpolation_type 
						= ceres::LineSearchInterpolationType(aux_vector[i]);
				}
			}

			if (!node[MIN_LINE_SEARCH_STEP_SIZE].empty())
			{
				vector<double> aux_vector;
				node[MIN_LINE_SEARCH_STEP_SIZE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].min_line_search_step_size = aux_vector[i];
				}
			}

			if (!node[LINE_SEARCH_SUFFICIENT_FUNCTION_DECREASE].empty())
			{
				vector<double> aux_vector;
				node[LINE_SEARCH_SUFFICIENT_FUNCTION_DECREASE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].line_search_sufficient_function_decrease 
						= aux_vector[i];
				}
			}

			if (!node[MAX_LINE_SEARCH_STEP_CONTRACTION].empty())
			{
				vector<double> aux_vector;
				node[MAX_LINE_SEARCH_STEP_CONTRACTION] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].max_line_search_step_contraction 
						= aux_vector[i];
				}
			}

			if (!node[MIN_LINE_SEARCH_STEP_CONTRACTION].empty())
			{
				vector<double> aux_vector;
				node[MIN_LINE_SEARCH_STEP_CONTRACTION] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].min_line_search_step_contraction 
						= aux_vector[i];
				}
			}

			if (!node[MAX_NUM_LINE_SEARCH_STEP_SIZE_ITERATIONS].empty())
			{
				vector<int> aux_vector;
				node[MAX_NUM_LINE_SEARCH_STEP_SIZE_ITERATIONS] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].max_num_line_search_step_size_iterations 
						= aux_vector[i];
				}
			}

			if (!node[MAX_NUM_LINE_SEARCH_DIRECTION_RESTARTS].empty())
			{
				vector<int> aux_vector;
				node[MAX_NUM_LINE_SEARCH_DIRECTION_RESTARTS] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].max_num_line_search_direction_restarts 
						= aux_vector[i];
				}
			}

			if (!node[LINE_SEARCH_SUFFICIENT_CURVATURE_DECREASE].empty())
			{
				vector<double> aux_vector;
				node[LINE_SEARCH_SUFFICIENT_CURVATURE_DECREASE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].line_search_sufficient_curvature_decrease 
						= aux_vector[i];
				}
			}

			if (!node[MAX_LINE_SEARCH_STEP_EXPANSION].empty())
			{
				vector<double> aux_vector;
				node[MAX_LINE_SEARCH_STEP_EXPANSION] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].max_line_search_step_expansion = aux_vector[i];
				}
			}

			if (!node[TRUST_REGION_STRATEGY_TYPE].empty())
			{
				vector<int> aux_vector;
				node[TRUST_REGION_STRATEGY_TYPE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].trust_region_strategy_type 
						= ceres::TrustRegionStrategyType(aux_vector[i]);
				}
			}

			if (!node[DOGLEG_TYPE].empty())
			{
				vector<int> aux_vector;
				node[DOGLEG_TYPE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].dogleg_type 
						= ceres::DoglegType(aux_vector[i]);
				}
			}

			if (!node[USE_NONMONOTONIC_STEPS].empty())
			{
				vector<int> aux_vector;
				node[USE_NONMONOTONIC_STEPS] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].use_nonmonotonic_steps = aux_vector[i] != 0;
				}
			}

			if (!node[MAX_CONSECUTIVE_NONMONOTONIC_STEPS].empty())
			{
				vector<int> aux_vector;
				node[MAX_CONSECUTIVE_NONMONOTONIC_STEPS] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].max_consecutive_nonmonotonic_steps 
						= aux_vector[i];
				}
			}

			if (!node[MAX_NUM_ITERATIONS].empty())
			{
				vector<int> aux_vector;
				node[MAX_NUM_ITERATIONS] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].max_num_iterations = aux_vector[i];
				}
			}

			if (!node[MAX_SOLVER_TIME_IN_SECONDS].empty())
			{
				vector<double> aux_vector;
				node[MAX_SOLVER_TIME_IN_SECONDS] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].max_solver_time_in_seconds = aux_vector[i];
				}
			}

			if (!node[NUM_THREADS].empty())
			{
				vector<int> aux_vector;
				node[NUM_THREADS] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].num_threads = aux_vector[i];
				}
			}

			if (!node[INITIAL_TRUST_REGION_RADIUS].empty())
			{
				vector<double> aux_vector;
				node[INITIAL_TRUST_REGION_RADIUS] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].initial_trust_region_radius = aux_vector[i];
				}
			}

			if (!node[MAX_TRUST_REGION_RADIUS].empty())
			{
				vector<double> aux_vector;
				node[MAX_TRUST_REGION_RADIUS] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].max_trust_region_radius = aux_vector[i];
				}
			}

			if (!node[MIN_TRUST_REGION_RADIUS].empty())
			{
				vector<double> aux_vector;
				node[MIN_TRUST_REGION_RADIUS] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].min_trust_region_radius = aux_vector[i];
				}
			}

			if (!node[MIN_RELATIVE_DECREASE].empty())
			{
				vector<double> aux_vector;
				node[MIN_RELATIVE_DECREASE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].min_relative_decrease = aux_vector[i];
				}
			}

			if (!node[MIN_LM_DIAGONAL].empty())
			{
				vector<double> aux_vector;
				node[MIN_LM_DIAGONAL] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].min_lm_diagonal = aux_vector[i];
				}
			}

			if (!node[MAX_LM_DIAGONAL].empty())
			{
				vector<double> aux_vector;
				node[MAX_LM_DIAGONAL] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].max_lm_diagonal = aux_vector[i];
				}
			}

			if (!node[MAX_NUM_CONSECUTIVE_INVALID_STEPS].empty())
			{
				vector<int> aux_vector;
				node[MAX_NUM_CONSECUTIVE_INVALID_STEPS] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].max_num_consecutive_invalid_steps 
						= aux_vector[i];
				}
			}

			if (!node[FUNCTION_TOLERANCE].empty())
			{
				vector<double> aux_vector;
				node[FUNCTION_TOLERANCE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].function_tolerance = aux_vector[i];
				}
			}

			if (!node[GRADIENT_TOLERANCE].empty())
			{
				vector<double> aux_vector;
				node[GRADIENT_TOLERANCE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].gradient_tolerance = aux_vector[i];
				}
			}

			if (!node[PARAMETER_TOLERANCE].empty())
			{
				vector<double> aux_vector;
				node[PARAMETER_TOLERANCE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].parameter_tolerance = aux_vector[i];
				}
			}

			if (!node[LINEAR_SOLVER_TYPE].empty())
			{
				vector<int> aux_vector;
				node[LINEAR_SOLVER_TYPE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].linear_solver_type 
						= ceres::LinearSolverType(aux_vector[i]);
				}
			}

			if (!node[PRECONDITIONER_TYPE].empty())
			{
				vector<int> aux_vector;
				node[PRECONDITIONER_TYPE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].preconditioner_type 
						= ceres::PreconditionerType(aux_vector[i]);
				}
			}

			if (!node[VISIBILITY_CLUSTERING_TYPE].empty())
			{
				vector<int> aux_vector;
				node[VISIBILITY_CLUSTERING_TYPE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].visibility_clustering_type 
						= ceres::VisibilityClusteringType(aux_vector[i]);
				}
			}

			if (!node[DENSE_LINEAR_ALGEBRA_LIBRARY_TYPE].empty())
			{
				vector<int> aux_vector;
				node[DENSE_LINEAR_ALGEBRA_LIBRARY_TYPE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].dense_linear_algebra_library_type 
						= ceres::DenseLinearAlgebraLibraryType(aux_vector[i]);
				}
			}

			if (!node[SPARSE_LINEAR_ALGEBRA_LIBRARY_TYPE].empty())
			{
				vector<int> aux_vector;
				node[SPARSE_LINEAR_ALGEBRA_LIBRARY_TYPE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].sparse_linear_algebra_library_type 
						= ceres::SparseLinearAlgebraLibraryType(aux_vector[i]);
				}
			}

			if (!node[NUM_LINEAR_SOLVER_THREADS].empty())
			{
				vector<int> aux_vector;
				node[NUM_LINEAR_SOLVER_THREADS] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].num_linear_solver_threads = aux_vector[i];
				}
			}

			if (!node[USE_EXPLICIT_SCHUR_COMPLEMENT].empty())
			{
				vector<int> aux_vector;
				node[USE_EXPLICIT_SCHUR_COMPLEMENT] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].use_explicit_schur_complement = aux_vector[i] != 0;
				}
			}

			if (!node[USE_POSTORDERING].empty())
			{
				vector<int> aux_vector;
				node[USE_POSTORDERING] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].use_postordering = aux_vector[i] != 0;
				}
			}

			if (!node[DYNAMIC_SPARSITY].empty())
			{
				vector<int> aux_vector;
				node[DYNAMIC_SPARSITY] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].dynamic_sparsity = aux_vector[i] != 0;
				}
			}

			if (!node[MIN_LINEAR_SOLVER_ITERATIONS].empty())
			{
				vector<int> aux_vector;
				node[MIN_LINEAR_SOLVER_ITERATIONS] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].min_linear_solver_iterations = aux_vector[i];
				}
			}

			if (!node[MAX_LINEAR_SOLVER_ITERATIONS].empty())
			{
				vector<int> aux_vector;
				node[MAX_LINEAR_SOLVER_ITERATIONS] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].max_linear_solver_iterations = aux_vector[i];
				}
			}

			if (!node[ETA].empty())
			{
				vector<double> aux_vector;
				node[ETA] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].eta = aux_vector[i];
				}
			}

			if (!node[JACOBI_SCALING].empty())
			{
				vector<int> aux_vector;
				node[JACOBI_SCALING] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].jacobi_scaling = aux_vector[i] != 0;
				}
			}

			if (!node[USE_INNER_ITERATIONS].empty())
			{
				vector<int> aux_vector;
				node[USE_INNER_ITERATIONS] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].use_inner_iterations = aux_vector[i] != 0;
				}
			}

			if (!node[INNER_ITERATION_TOLERANCE].empty())
			{
				vector<double> aux_vector;
				node[INNER_ITERATION_TOLERANCE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].inner_iteration_tolerance = aux_vector[i];
				}
			}

			if (!node[LOGGING_TYPE].empty())
			{
				vector<int> aux_vector;
				node[LOGGING_TYPE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].logging_type 
						= ceres::LoggingType(aux_vector[i]);
				}
			}

			if (!node[MINIMIZER_PROGRESS_TO_STDOUT].empty())
			{
				vector<int> aux_vector;
				node[MINIMIZER_PROGRESS_TO_STDOUT] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].minimizer_progress_to_stdout = aux_vector[i] != 0;
				}
			}

			//if (!node[TRUST_REGION_PROBLEM_DUMP_DIRECTORY].empty())
			//{
			//	vector<string> aux_vector;
			//	node[TRUST_REGION_PROBLEM_DUMP_DIRECTORY] >> aux_vector;
			//	for (size_t i = 0; i < aux_vector.size(); ++i)
			//	{
			//		options[i].trust_region_problem_dump_directory 
			//			= aux_vector[i];
			//	}
			//}

			if (!node[TRUST_REGION_PROBLEM_DUMP_FORMAT_TYPE].empty())
			{
				vector<int> aux_vector;
				node[TRUST_REGION_PROBLEM_DUMP_FORMAT_TYPE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].trust_region_problem_dump_format_type 
						= ceres::DumpFormatType(aux_vector[i]);
				}
			}

			if (!node[CHECK_GRADIENTS].empty())
			{
				vector<int> aux_vector;
				node[CHECK_GRADIENTS] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].check_gradients = aux_vector[i] != 0;
				}
			}

			if (!node[GRADIENT_CHECK_RELATIVE_PRECISION].empty())
			{
				vector<double> aux_vector;
				node[GRADIENT_CHECK_RELATIVE_PRECISION] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].gradient_check_relative_precision 
						= aux_vector[i];
				}
			}

			if (!node[NUMERIC_DERIVATIVE_RELATIVE_STEP_SIZE].empty())
			{
				vector<double> aux_vector;
				node[NUMERIC_DERIVATIVE_RELATIVE_STEP_SIZE] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].numeric_derivative_relative_step_size 
						= aux_vector[i];
				}
			}

			if (!node[UPDATE_STATE_EVERY_ITERATION].empty())
			{
				vector<int> aux_vector;
				node[UPDATE_STATE_EVERY_ITERATION] >> aux_vector;
				for (size_t i = 0; i < aux_vector.size(); ++i)
				{
					options[i].update_state_every_iteration = aux_vector[i] != 0;
				}
			}
		}

	};


	struct Parameters {

		// Default constructor that sets up a generic sparse problem.
		Parameters() {
			// Data term
			data_weights.resize(NUM_OPTIMIZATIONS, 0);
			data_huber_width.resize(NUM_OPTIMIZATIONS, 0);

			// Lighting confidence
			specular_weight_var = 0.1;

			// SH Coefficients estimation
			sh_order = 2;
			albedo_percentile = 0.98;
			use_lower_bound_shading = false;
			use_upper_bound_shading = false;

			// Smoothing Weight
			color_diff_threshold = 0.5;
			color_diff_var = 0.05;
			depth_diff_var = 50;
			use_depth_weight = true;
			smooth_specular_weight_var = 0.1;

			// Albedo estimation
			smooth_albedo_weight = 0.1;
			smooth_albedo_huber_width = 0;
			use_lower_bound_albedo = false;
			use_upper_bound_albedo = false;
			albedo_initialization = WHITE_ALBEDO;

			// Local lighting variation estimation
			smooth_local_lighting_weight = 1.0;
			local_lighting_magnitude_weight = 1.0;
			smooth_local_lighting_huber_width = 0;
			local_lighting_magnitude_huber_width = 0;
			use_lower_bound_lighting = false;
			use_upper_bound_lighting = false;

			// Join estimation of albedo and local lighting variation
			combine_albedo_lighting = false;

			// Shape from shading
			estimate_shape = true;
			displacement_weight = 4e-3;
			laplacian_smooth_weight = 7.5e-3;
			tv_weight = 7.5e-3;
			use_normal_for_shape = false;

			// Solution refinement
			refine_solution = false;

			// I/O parameters
			input_mesh_filename = "C:/template.ply";
			input_specular_image_filename = "C:/specular.png";
			input_intrinsics_filename = "C:/intrinsics.txt";

			output_sh_coeff_filename = "C:/sh_coeff.txt";
			output_albedo_mesh_filename = "C:/albedo.ply";
			output_shading_mesh_filename = "C:/shading.ply";
			output_local_lighting_mesh_filename = "C:/local_lighting.ply";
			output_est_diffuse_mesh_filename = "C:/diffuse.ply";
			output_diffuse_diff_mesh_filename = "C:/diffuse_diff.ply";
			output_est_intensity_mesh_filename = "C:/est_intensity.ply";
			output_orig_intensity_mesh_filename = "C:/orig_intensity.ply";

			save_sh_coeff = true;
			save_albedo = true;
			save_shading = true;
			save_local_lighting = true;
			save_est_diffuse = true;
			save_diffuse_diff = true;
			save_est_intensity = true;
			save_orig_intensity = true;
			save_mesh_binary = true;

			ceres_solver = SolverOptions();
		}

		// Data term weight and Huber width
		vector<double> data_weights;
		vector<double> data_huber_width;

		// Lighting confidence
		double specular_weight_var;

		// SH Coefficients estimation
		int sh_order;
		double albedo_percentile;
		bool use_lower_bound_shading;
		bool use_upper_bound_shading;

		// Smoothing Weight
		double color_diff_threshold;
		double color_diff_var;
		double depth_diff_var;
		bool use_depth_weight;
		double smooth_specular_weight_var;

		// Albedo estimation
		double smooth_albedo_weight;
		double smooth_albedo_huber_width;
		bool use_lower_bound_albedo;
		bool use_upper_bound_albedo;
		int albedo_initialization;

		// Local lighting variation estimation
		double smooth_local_lighting_weight;
		double local_lighting_magnitude_weight;
		double smooth_local_lighting_huber_width;
		double local_lighting_magnitude_huber_width;
		bool use_lower_bound_lighting;
		bool use_upper_bound_lighting;

		// Join estimation of albedo and local lighting variation
		bool combine_albedo_lighting;

		// Shape from shading
		bool estimate_shape;
		double displacement_weight;
		double laplacian_smooth_weight;
		double tv_weight;
		bool use_normal_for_shape;

		// Solution refinement
		bool refine_solution;

		// I/O parameters
		string input_mesh_filename;
		string input_specular_image_filename;
		string input_intrinsics_filename;

		string output_sh_coeff_filename;
		string output_albedo_mesh_filename;
		string output_shading_mesh_filename;
		string output_local_lighting_mesh_filename;
		string output_est_diffuse_mesh_filename;
		string output_diffuse_diff_mesh_filename;
		string output_est_intensity_mesh_filename;
		string output_orig_intensity_mesh_filename;
		
		bool save_sh_coeff;
		bool save_albedo;
		bool save_shading;
		bool save_local_lighting;
		bool save_est_diffuse;
		bool save_diffuse_diff;
		bool save_est_intensity;
		bool save_orig_intensity;
		bool save_mesh_binary;

		// Ceres Solver options
		SolverOptions ceres_solver;

		// Load parameters from an YAML file
		inline void Parameters::load(const std::string &_filename)
		{
			cv::FileStorage fs(_filename, cv::FileStorage::READ);

			// Data terms
			if (!fs[DATA_WEIGHTS].empty())
			{
				fs[DATA_WEIGHTS] >> data_weights;
			}

			if (!fs[DATA_HUBER_WIDTH].empty())
			{
				fs[DATA_HUBER_WIDTH] >> data_huber_width;
			}

			// Lighting confidence
			if (!fs[SPECULAR_WEIGHT_VAR].empty())
			{
				fs[SPECULAR_WEIGHT_VAR] >> specular_weight_var;
			}

			// SH Coefficients estimation
			if (!fs[SH_ORDER].empty())
			{
				fs[SH_ORDER] >> sh_order;
			}

			if (!fs[ALBEDO_PERCENTILE].empty())
			{
				fs[ALBEDO_PERCENTILE] >> albedo_percentile;
			}

			if (!fs[USE_LOWER_BOUND_SHADING].empty())
			{
				fs[USE_LOWER_BOUND_SHADING] >> use_lower_bound_shading;
			}

			if (!fs[USE_UPPER_BOUND_SHADING].empty())
			{
				fs[USE_UPPER_BOUND_SHADING] >> use_upper_bound_shading;
			}


			// Smoothing Weight
			if (!fs[COLOR_DIFF_THRESHOLD].empty())
			{
				fs[COLOR_DIFF_THRESHOLD] >> color_diff_threshold;
			}

			if (!fs[COLOR_DIFF_VAR].empty())
			{
				fs[COLOR_DIFF_VAR] >> color_diff_var;
			}

			if (!fs[DEPTH_DIFF_VAR].empty())
			{
				fs[DEPTH_DIFF_VAR] >> depth_diff_var;
			}

			if (!fs[USE_DEPTH_WEIGHT].empty())
			{
				fs[USE_DEPTH_WEIGHT] >> use_depth_weight;
			}

			if (!fs[SMOOTH_SPECULAR_WEIGHT_VAR].empty())
			{
				fs[SMOOTH_SPECULAR_WEIGHT_VAR] >> smooth_specular_weight_var;
			}

			// Albedo estimation
			if (!fs[SMOOTH_ALBEDO_WEIGHT].empty())
			{
				fs[SMOOTH_ALBEDO_WEIGHT] >> smooth_albedo_weight;
			}

			if (!fs[SMOOTH_ALBEDO_HUBER_WIDTH].empty())
			{
				fs[SMOOTH_ALBEDO_HUBER_WIDTH] >> smooth_albedo_huber_width;
			}

			if (!fs[USE_LOWER_BOUND_ALBEDO].empty())
			{
				fs[USE_LOWER_BOUND_ALBEDO] >> use_lower_bound_albedo;
			}

			if (!fs[USE_UPPER_BOUND_ALBEDO].empty())
			{
				fs[USE_UPPER_BOUND_ALBEDO] >> use_upper_bound_albedo;
			}

			if (!fs[ALBEDO_INITIALIZATION].empty())
			{
				fs[ALBEDO_INITIALIZATION] >> albedo_initialization;
			}


			// Local lighting variation estimation
			if (!fs[SMOOTH_LOCAL_LIGHTING_WEIGHT].empty())
			{
				fs[SMOOTH_LOCAL_LIGHTING_WEIGHT] >> smooth_local_lighting_weight;
			}

			if (!fs[LOCAL_LIGHTING_MAGNITUDE_WEIGHT].empty())
			{
				fs[LOCAL_LIGHTING_MAGNITUDE_WEIGHT] >> local_lighting_magnitude_weight;
			}

			if (!fs[SMOOTH_LOCAL_LIGHTING_HUBER_WIDTH].empty())
			{
				fs[SMOOTH_LOCAL_LIGHTING_HUBER_WIDTH] >> smooth_local_lighting_huber_width;
			}

			if (!fs[LOCAL_LIGHTING_MAGNITUDE_HUBER_WIDTH].empty())
			{
				fs[LOCAL_LIGHTING_MAGNITUDE_HUBER_WIDTH] >> local_lighting_magnitude_huber_width;
			}

			if (!fs[USE_LOWER_BOUND_LIGHTING].empty())
			{
				fs[USE_LOWER_BOUND_LIGHTING] >> use_lower_bound_lighting;
			}

			if (!fs[USE_UPPER_BOUND_LIGHTING].empty())
			{
				fs[USE_UPPER_BOUND_LIGHTING] >> use_upper_bound_lighting;
			}


			// Join estimation of albedo and local lighting variation
			if (!fs[COMBINE_ALBEDO_LIGHTING].empty())
			{
				fs[COMBINE_ALBEDO_LIGHTING] >> combine_albedo_lighting;
			}


			// Shape from shading
			if (!fs[ESTIMATE_SHAPE].empty())
			{
				fs[ESTIMATE_SHAPE] >> estimate_shape;
			}

			if (!fs[DISPLACEMENT_WEIGHT].empty())
			{
				fs[DISPLACEMENT_WEIGHT] >> displacement_weight;
			}

			if (!fs[LAPLACIAN_SMOOTH_WEIGHT].empty())
			{
				fs[LAPLACIAN_SMOOTH_WEIGHT] >> laplacian_smooth_weight;
			}

			if (!fs[TV_WEIGHT].empty())
			{
				fs[TV_WEIGHT] >> tv_weight;
			}

			if (!fs[USE_NORMAL_FOR_SHAPE].empty())
			{
				fs[USE_NORMAL_FOR_SHAPE] >> use_normal_for_shape;
			}

			// Solution refinement
			if (!fs[REFINE_SOLUTION].empty())
			{
				fs[REFINE_SOLUTION] >> refine_solution;
			}

			// I/O parameters
			if (!fs[INPUT_MESH_FILENAME].empty())
			{
				fs[INPUT_MESH_FILENAME] >> input_mesh_filename;
			}

			if (!fs[INPUT_SPECULAR_IMAGE_FILENAME].empty())
			{
				fs[INPUT_SPECULAR_IMAGE_FILENAME] >> input_specular_image_filename;
			}

			if (!fs[INPUT_INTRINSICS_FILENAME].empty())
			{
				fs[INPUT_INTRINSICS_FILENAME] >> input_intrinsics_filename;
			}


			if (!fs[OUTPUT_SH_COEFF_FILENAME].empty())
			{
				fs[OUTPUT_SH_COEFF_FILENAME] >> output_sh_coeff_filename;
			}

			if (!fs[OUTPUT_ALBEDO_MESH_FILENAME].empty())
			{
				fs[OUTPUT_ALBEDO_MESH_FILENAME] >> output_albedo_mesh_filename;
			}

			if (!fs[OUTPUT_SHADING_MESH_FILENAME].empty())
			{
				fs[OUTPUT_SHADING_MESH_FILENAME] >> output_shading_mesh_filename;
			}

			if (!fs[OUTPUT_LOCAL_LIGHTING_MESH_FILENAME].empty())
			{
				fs[OUTPUT_LOCAL_LIGHTING_MESH_FILENAME] >> output_local_lighting_mesh_filename;
			}

			if (!fs[OUTPUT_EST_DIFFUSE_MESH_FILENAME].empty())
			{
				fs[OUTPUT_EST_DIFFUSE_MESH_FILENAME] >> output_est_diffuse_mesh_filename;
			}

			if (!fs[OUTPUT_DIFFUSE_DIFF_MESH_FILENAME].empty())
			{
				fs[OUTPUT_DIFFUSE_DIFF_MESH_FILENAME] >> output_diffuse_diff_mesh_filename;
			}

			if (!fs[OUTPUT_EST_INTENSITY_MESH_FILENAME].empty())
			{
				fs[OUTPUT_EST_INTENSITY_MESH_FILENAME] >> output_est_intensity_mesh_filename;
			}

			if (!fs[OUTPUT_ORIG_INTENSITY_MESH_FILENAME].empty())
			{
				fs[OUTPUT_ORIG_INTENSITY_MESH_FILENAME] >> output_orig_intensity_mesh_filename;
			}


			if (!fs[SAVE_SH_COEFF].empty())
			{
				fs[SAVE_SH_COEFF] >> save_sh_coeff;
			}

			if (!fs[SAVE_ALBEDO].empty())
			{
				fs[SAVE_ALBEDO] >> save_albedo;
			}

			if (!fs[SAVE_SHADING].empty())
			{
				fs[SAVE_SHADING] >> save_shading;
			}

			if (!fs[SAVE_LOCAL_LIGHTING].empty())
			{
				fs[SAVE_LOCAL_LIGHTING] >> save_local_lighting;
			}

			if (!fs[SAVE_EST_DIFFUSE].empty())
			{
				fs[SAVE_EST_DIFFUSE] >> save_est_diffuse;
			}

			if (!fs[SAVE_DIFFUSE_DIFF].empty())
			{
				fs[SAVE_DIFFUSE_DIFF] >> save_diffuse_diff;
			}

			if (!fs[SAVE_EST_INTENSITY].empty())
			{
				fs[SAVE_EST_INTENSITY] >> save_est_intensity;
			}

			if (!fs[SAVE_ORIG_INTENSITY].empty())
			{
				fs[SAVE_ORIG_INTENSITY] >> save_orig_intensity;
			}

			if (!fs[SAVE_MESH_BINARY].empty())
			{
				fs[SAVE_MESH_BINARY] >> save_mesh_binary;
			}

			// Ceres Solver options
			if (!fs[CERES_SOLVER].empty())
			{
				fs[CERES_SOLVER] >> ceres_solver;
			}
		}

		inline void Parameters::save(const std::string &_filename)
		{

		}
	};

void read(const cv::FileNode& node, SolverOptions& options,
	const SolverOptions& default_options = SolverOptions());

inline void read(const cv::FileNode& node, SolverOptions& options,
	const SolverOptions& default_options)
{
	if (node.empty())
		options = default_options;
	else
		options.read(node);
}

#endif  // PARAMETERS_H_
