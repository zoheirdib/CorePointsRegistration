// Copyright Zoheir DIB - ZDEngine, Inc. All Rights Reserved.
// Mail : Zohier.dib@gmail.com

#pragma once

//CoreMinimal
#include <ZDCoreTypes.h>
#include <chrono>

#include "Geom\PointsRegistration\ZDCPDRegistration.h"

namespace ZDEngine
{
	template< typename source_matrix_exp_type, typename target_matrix_exp_type>
	class ZDCPDNonRigidRegistration
	{

		/*!

		REQUIREMENTS ON source_matrix_exp_type and target_matrix_exp_type
		source_matrix_exp_type and target_matrix_exp_type should have the same NC
		source_matrix_exp_type and target_matrix_exp_type should be of the same type T
		the type T should be double or float


		WHAT THIS OBJECT REPRESENTS
		this object implements the non-rigid registration algorithm as discribed if the following paper
		A.Myronenko and X.Song - Point Set Registration: Coherent Point Drift - 2009 - IEEE Trans
		!*/

	public:

		const static long source_NR = source_matrix_exp_type::NR;
		const static long target_NR = target_matrix_exp_type::NR;

		// You have supplied an invalid type of matrix_exp_type for source and target.  You have
		// to use this object with matrices the same number of columns (DIM).
		COMPILE_TIME_ASSERT(source_matrix_exp_type::NC == target_matrix_exp_type::NC);

		const static long DIM = source_matrix_exp_type::NC;

		typedef typename source_matrix_exp_type::type source_type;
		typedef typename target_matrix_exp_type::type target_type;

		// You have supplied an invalid type of matrix_exp_type for source and target.  You have
		// to use this object with matrices of the same type
		COMPILE_TIME_ASSERT((Is_Same_Type<source_type, target_type>::value));

		// You have supplied an invalid type of matrix_exp_type.  You have
		// to use this object with matrices that contain float or double type data.
		COMPILE_TIME_ASSERT((Is_Same_Type<float, source_type>::value ||
			Is_Same_Type<double, source_type>::value));

		typedef typename source_matrix_exp_type::type type;

		typedef typename source_matrix_exp_type::mem_manager_type mem_manager_type;
		typedef typename source_matrix_exp_type::layout_type layout_type;
		typedef matrix<type, source_NR, DIM, mem_manager_type, layout_type>  source_matrix_type;
		typedef matrix<type, target_NR, DIM, mem_manager_type, layout_type>  target_matrix_type;


		typedef matrix<uint32, source_NR, 1, mem_manager_type, layout_type> correspondance_column_vector_type;

		typedef ZDPoint<type, DIM> TPoint;

		typedef matrix<type, source_NR, 1, mem_manager_type, layout_type>  p1_matrix_type;
		typedef matrix<type, target_NR, 1, mem_manager_type, layout_type>  pt1_matrix_type;
		typedef typename source_matrix_type px_matrix_type;

		typedef CPDImpl::Probabilities<p1_matrix_type, pt1_matrix_type, px_matrix_type, correspondance_column_vector_type> probabilities_type;

		struct NonRigidResult
		{
			source_matrix_type points;
			correspondance_column_vector_type correspondence;
			double sigma2;
		};

		/*!
		ensures
		- # sets whether the correspondence vector will be computed.
		!*/
		void correspondence(bool correspondence)
		{
			m_correspondence = correspondence;
		}

		/*!
		ensures
		- # sets the max iterations for this transform.
		!*/
		void max_iterations(double max_iterations)
		{
			m_max_iterations = static_cast<size_t>(max_iterations);
		}

		/*!
		ensures
		- # sets whether to normalize the points before running cpd.
		!*/
		void normalize(bool normalize)
		{
			m_normalize = normalize;
		}

		/*!
		ensures
		- # sets the outlier tolerance.
		!*/
		void outliers(double outliers)
		{
			m_outliers = outliers;
		}

		/*!
		ensures
		- # sets the sigma2 value for this transform.
		!*/
		void sigma2(double sigma2)
		{
			m_sigma2 = sigma2;
		}

		/*!
		ensures
		- # sets the final tolerance.
		!*/
		void tolerance(double tolerance)
		{
			m_tolerance = tolerance;
		}

		/*!
		ensures
		- # sets whether this rigid transform allows reflections.
		!*/
		void reflections(bool reflections)
		{
			m_reflections = reflections;
		}

		/*!
		ensures
		- # sets whether this rigid transform allows scaling.
		!*/
		void scale(bool scale)
		{
			m_scale = scale;
		}

		/*!
		ensures
		- # sets the beta.
		!*/
		void beta(double beta)
		{
			m_beta = beta;
		}

		/*!
		ensures
		- # sets the lamda.
		!*/
		void lambda(double lambda)
		{
			m_lambda = lambda;
		}

		/*!
		ensures
		- # get the result
		- # points ==  result.points
		!*/
		template<typename matrix_type>
		void getTransformedPoints(matrix_type& points)
		{
			points = result.points;
		}

		/*!
		ensures
		- # get the correspondance
		- # correspondance ==  result.correspondance
		!*/
		void getCorrespondance(correspondance_column_vector_type& correspondence)
		{
			correspondence = result.correspondence;
		}

		/*!
		ensures
		- # get the runtime
		!*/
		const std::chrono::microseconds& getRunTime()
		{
			return m_runtime;
		}

		/*!
		ensures
		- # get the number of iterations
		!*/
		uint32 getIterations()
		{
			return m_iterations;
		}

		/*!
		ensures
		- # run the registration process
		!*/
		void run(const source_matrix_exp_type& source, const target_matrix_exp_type& target)
		{
			auto tic = std::chrono::high_resolution_clock::now();

			target_matrix_type fixed;
			source_matrix_type moving;

			TPoint fixed_mean, moving_mean;
			double fixed_scale, moving_scale;

			if (m_normalize)
			{
				normalizePoints(target, fixed, fixed_mean.asBase(), fixed_scale);
				normalizePoints(source, moving, moving_mean.asBase(), moving_scale);
			}
			else
			{
				fixed = target;
				moving = source;
			}

			this->init(fixed, moving);

			result.points = moving;

			if (m_sigma2 == 0.0)
			{
				result.sigma2 = CPDImpl::default_sigma2(fixed, moving);
			}
			else
			{
				result.sigma2 = m_sigma2;
			}


			uint32 iter = 0;
			double ntol = m_tolerance + 10.0;
			double l = 0.;

			probabilities_type probabilities;
			probabilities.sigmaInit = result.sigma2;

			while (iter < m_max_iterations && ntol > m_tolerance &&
				result.sigma2 > 10 * std::numeric_limits<double>::epsilon())
			{
				CPDImpl::ZDGaussTransformFGT(fixed, result.points, result.sigma2, m_outliers, probabilities);
			
				ntol = std::abs((probabilities.l - l) / probabilities.l);
				l = probabilities.l;

				this->compute_one(fixed, moving, probabilities, result.sigma2, result);

				++iter;
			}

			if (m_correspondence)
			{
				probabilities_type probabilities;
				CPDImpl::ZDGaussTransform(
					fixed, result.points, result.sigma2, m_outliers, probabilities);
				result.correspondence = probabilities.correspondence;
				ZD_ASSERT(result.correspondence.nr() > 0);
			}

			if (m_normalize)
			{
				denormalizePoints(result.points, result.points, fixed_mean.asBase(), fixed_scale);
			}

			auto toc = std::chrono::high_resolution_clock::now();
			m_runtime = std::chrono::duration_cast<std::chrono::microseconds>(toc - tic);
			m_iterations = iter;
		}

		/*!
		ensures
		- # initialize this transform for the provided matrices
		!*/
		void init(const target_matrix_type& fixed, const source_matrix_type& moving) 
		{
			typedef typename matrix_traits<source_matrix_type>::type type;
			CPDImpl::affinity(moving, moving, m_beta, m_g);
			m_w.set_size(moving.nr(), moving.nc());
			ZDMemory::Memset(m_w.asPtr(), 0, m_w.size() * sizeof(type));
		}

		/*!
		ensures
		- # computes one iteration of the transform.
		!*/
		void compute_one(const target_matrix_type& fixed, const source_matrix_type& moving,
			const probabilities_type& probabilities,
			double sigma2, NonRigidResult& result)
		{
			typedef typename matrix_traits<source_matrix_type>::type type;

			uint32 cols = fixed.nc();
			auto dp = diagm(probabilities.p1);

			qr_decomposition<matrix<type>> QRSolver(dp * m_g + m_lambda * sigma2 * identity_matrix<type>(moving.nr()));
			matrix<type>  w = QRSolver.solve(probabilities.px - dp * moving);

			result.points = moving + m_g * w;
			double np = sum(probabilities.p1);
			result.sigma2 = std::abs(
				(sum(pointwise_multiply(pow(fixed, 2), repmat(probabilities.pt1, 1, cols)))
					+
					sum(pointwise_multiply(pow(result.points, 2),
						repmat(probabilities.p1, 1, cols)))
					-
					2 * trace(trans(probabilities.px) * result.points)) /
					(np * cols));
		}

	private:
		size_t m_max_iterations = 150;
		NonRigidResult result;
		bool m_normalize = true;
		bool m_correspondence = false;
		bool m_linked = true;
		double m_outliers = 0.1;
		double m_sigma2 = 0.0;
		double m_tolerance = 1e-5;
		bool m_reflections = false;
		bool m_scale = false;
		matrix<type> m_g;
		matrix<type> m_w;
		double m_lambda = 3.0;
		double m_beta = 3.0;
		std::chrono::microseconds m_runtime;
		uint32 m_iterations;
	};
}
