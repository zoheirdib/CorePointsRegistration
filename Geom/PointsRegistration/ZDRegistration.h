// Copyright Zoheir DIB - ZDEngine, Inc. All Rights Reserved.
// Mail : Zohier.dib@gmail.com

#pragma once

//STD
#include <chrono>

//CoreComon
#include <ZDCoreTypes.h>
#include <Threads\ZDParallelFor.h>

//CoreMaths
#include "Maths\FGT\FGTCModel.c"
#include "Maths\FGT\FGTCPredict.c"

//CoreGeom
#include "Geom\ZDPointsTransform.h"


namespace ZDEngine
{

	namespace CPDImpl
	{	
		/*!
		ensures
		- # computes the default sigma2 for the given matrices.
		!*/
		template<typename source_matrix_exp_type, typename target_matrix_exp_type>
		static double default_sigma2(const target_matrix_exp_type& fixed, const source_matrix_exp_type& moving)
		{
			return (
				(moving.nr() * trace(trans(fixed) * fixed)) +
				(fixed.nr() * trace(trans(moving) * moving)) -
				2 * (sum_rows(fixed) * trans(sum_rows(moving)))) /
				(fixed.nr() * moving.nr() * fixed.nc());
		}

		namespace impl
		{
			class ZDCPDRegistration_affinity_error : public ZDError 
			{
				public: ZDCPDRegistration_affinity_error(const std::string& str) : ZDError(EFATAL, str) {}
			};
		}

		/*!
		ensures
		- # computes the affinity matrix between the two matrices
		!*/
		template<typename source_matrix_exp_type, typename target_matrix_exp_type, typename g_matrix_exp_type>
		void affinity(const source_matrix_exp_type& x, const target_matrix_exp_type& y, double beta, g_matrix_exp_type &g) 
		{
			double k = -2.0 * beta * beta;
			uint32 x_rows = static_cast<uint32>(x.nr());
			uint32 y_rows = static_cast<uint32>(y.nr());

			try 
			{
				g.set_size(x_rows, y_rows);
			}
			catch (const std::bad_alloc& err)
			{
				std::stringstream msg;
				msg << "Unable to allocate " << x_rows << " by " << y_rows
					<< " affinity matrix, try again with fewer points" << " Error from std::bad_alloc "<< err.what();
				throw impl::ZDCPDRegistration_affinity_error(msg.str());
			}

			ZDParallelFor(0, y_rows, [&](uint32 i)
			{
				set_colm(g, i) = exp(sum_cols(pow((x - repmat(rowm(y, i), x_rows, 1)), 2)) / k);
			});
		}

		/*!
		ensures
		- # denormalize the given rigid transform
		!*/
		template<typename T, long DIM, typename MM, typename L, typename point_type>
		void denormalizeTransform(ZDPointTransformRigide<T,DIM, MM,L> &transform, const point_type& fixed_mean, const point_type& moving_mean, double fixed_scale, double moving_scale)
		{
			transform.scale = transform.scale * fixed_scale / moving_scale;
			transform.translation = fixed_scale * transform.translation + fixed_mean.asBase() -
			transform.scale * transform.rotation * moving_mean.asBase();
		}

		/*!
		ensures
		- # denormalize the given affine transform
		!*/
		template<typename T, long DIM, typename MM, typename L, typename point_type>
		void denormalizeTransform(ZDPointTransformAffine<T, DIM, MM, L> &transform, const point_type& fixed_mean, const point_type& moving_mean, double fixed_scale, double moving_scale)
		{
			transform.translation = fixed_scale * transform.translation + fixed_mean.asBase() -
			transform.rotation *moving_mean.asBase();
			transform.rotation = transform.rotation * fixed_scale / moving_scale;
		}

		/*!
		ensures
		- # computes the probability matrices by camparing two data sets with a GaussTransform
		!*/
		template<typename p1_matrix_type, typename pt1_matrix_type, typename px_matrix_type, typename correspondence_type>
		struct Probabilities
		{
			// The probability matrix, multiplied by the identity matrix.
			p1_matrix_type p1;
			// The probability matrix, transposes, multiplied by the identity matrix.
			pt1_matrix_type pt1;
			// The probability matrix multiplied by the fixed points.
			px_matrix_type px;
			// The total error.
			double l;
			// The correspondence vector between the two datasets.
			correspondence_type correspondence;
			//the initial sigma for non regid registration
			double sigmaInit;
		};

		/*!
		ensures
		- # computes the gaussian transform
		!*/
		template<typename target_matrix_type, typename source_matrix_type, typename probabilities_type>
		void ZDGaussTransform(const target_matrix_type& fixed, const source_matrix_type& moving,
			double sigma2, double outliers, probabilities_type& prob)
		{

			typedef typename matrix_traits<source_matrix_type>::type type;
			typedef typename matrix_traits<source_matrix_type>::const_ret_type const_ret_type;
			typedef typename matrix_traits<source_matrix_type>::mem_manager_type mem_manager_type;
			const static long Source_NR = matrix_traits<source_matrix_type>::NR;
			const static long Target_NR = matrix_traits<target_matrix_type>::NR;

			const static long NC = matrix_traits<source_matrix_type>::NC;
			typedef typename matrix_traits<source_matrix_type>::layout_type layout_type;

			double ksig = -2.0 * sigma2;
			long cols = fixed.nc();
			double outlier_tmp =
				(outliers * moving.nr() * std::pow(-ksig * ZDMath::pi, 0.5 * cols)) /
				((1 - outliers) * fixed.nr());

			matrix<type, Source_NR, 1> p; p.set_size(moving.nr());
			matrix<type, Source_NR, 1> p1_max;  p1_max.set_size(moving.nr());
			ZDMemory::Memzero(p1_max.asPtr(), p1_max.size() * sizeof(type));

			prob.p1.set_size(moving.nr());
			ZDMemory::Memzero(prob.p1.asPtr(), prob.p1.size() * sizeof(type));

			prob.pt1.set_size(fixed.nr());

			prob.px.set_size(moving.nr(), cols);
			ZDMemory::Memzero(prob.px.asPtr(), prob.px.size() * sizeof(type));

			prob.correspondence.set_size(moving.nr());
			set_all_elements(prob.correspondence, -1);
			
			prob.l = 0.0;

			for (uint32 i = 0; i < static_cast<uint32>(fixed.nr()); ++i)
			{
				double sp = 0;

				for (uint32 j = 0; j < static_cast<uint32>(moving.nr()); ++j)
				{
					double razn = sum(pow(rowm(fixed, i) - rowm(moving, j), 2));
					p(j) = std::exp(razn / ksig);
					sp += p(j);
				}
				sp += outlier_tmp;
				prob.pt1(i) = 1 - outlier_tmp / sp;

				for (uint32 j = 0; j < static_cast<uint32>(moving.nr()); ++j)
				{
					prob.p1(j) += p(j) / sp;
					set_rowm(prob.px, j) += rowm(fixed, i) * p(j) / sp;
					if (p(j) / sp > p1_max(j))
					{
						prob.correspondence(j) = i;
						p1_max(j) = p(j) / sp;
					}
				}
				prob.l += -std::log(sp);
			}
			prob.l += cols * fixed.nr() * std::log(sigma2) / 2;

		}

	}

}