// Copyright Zoheir DIB - ZDEngine, Inc. All Rights Reserved.
// Mail : Zohier.dib@gmail.com

#pragma once

//CoreCommon
#include <ZDCoreTypes.h>
#include <chrono>

//CoreGeom
#include "Geom\ZDPointsTransform.h"

//CorePointsRegistration
#include "Geom\PointsRegistration\ZDRegistration.h"
#include "Geom\PointsRegistration\ZDCPDRigidRegistration.h"
#include "Geom\PointsRegistration\ZDCPDAffineRegistration.h"
#include "Geom\PointsRegistration\ZDCPDNonRigidRegistration.h"

namespace ZDEngine
{

	/*!
	ensures
	- # computes the transfomration that maps the target to the source points with the given registration algorithm
	- # result is sets with the transformed points
	!*/
	template< typename Registrator, typename PointsType>
	void ZDRegistration(const PointsType& source, const PointsType& target, PointsType& result)
	{
		Registrator registrator;
		registrator.run(source, target);
		registrator.getResult(result);
	}

	/*!
	ensures
	- # computes the rigid transfomration that maps the target to the source points
	- # result is sets with the transformed points
	- # transform is sets with the computed rigid transformation
	!*/

	template<typename PointsType, typename Registrator, typename TransformType>
	typename std::enable_if<std::is_same<Registrator, ZDCPDRigidRegistration<PointsType, PointsType>>::value>::type ZDRigidRegistration(const PointsType& source, const PointsType& target, PointsType& result, TransformType& transform)
	{
		Registrator registrator;
		registrator.run(source, target);
		registrator.getResult(result);
		registrator.getTransform(transform);
	}
}