/**
* @file    NearestNeighbor.h
* @author  T.Sato
* @date    2015.05.05
* @version 1.0
*/

#ifndef __NEAREST_NEIGHBOR__
#define __NEAREST_NEIGHBOR__

#include <float.h>
#include "define.h"
#include "point.h"

/**
* @brief distance calculation
* @return distance
*/
template<typename data_t, typename query_t>
inline double Distance(
	int dim,	   //!< [in] dimension
	data_t* sample,//!< [in] data
	query_t* query //!< [in] query
	)
{
	double dist = 0.0;

	for (int d = 0; d < dim; ++d){
		dist += NORM(query[d] - sample[d]);
	}

	return dist;
}

/**
* @brief distance calculation
* @return distance
*/
template<typename data_t, typename query_t>
inline double Distance(
	const int& dim,
	data_t* sample,
	query_t* query,
	const double& Limit
	)
{

	double dist = 0.0;
	for (int d = 0; d < dim && dist < Limit; ++d)
	{
		dist += NORM((*query++) - (*sample++));
	}

	return dist;
}

/**
* @brief distance calculation
* @return distance
*/
//template<typename data_t, typename query_t>
//inline double Distance(
//	const int& dim,
//	data_t* const & sample,
//	query_t* const & query,
//	const double& Limit,
//	int& d
//	)
//{
//
//	double dist = 0.0;
//
//	for (; d < dim && dist < Limit; ++d)
//	{
//		dist += NORM(query[d] - sample[d]);
//	}
//
//	return dist;
//}

/**
* @brief k-means
*/
template<typename data_t, typename query_t>
inline int NearestNeighbor(
	int dim,
	int num,
	data_t** sample,
	query_t* query
	)
{

	int NNidx = 0;
	double distance;
	double NNdis = Distance(dim, *sample++, query);
	for (int n = 1; n < num; n++){

		distance = Distance(dim, *sample++, query, NNdis);
		if (distance < NNdis)
		{
			NNdis = distance;
			NNidx = n;
		}
	}
	return NNidx;
}

/**
* @brief k-means
*/
template<typename data_t, typename query_t>
inline void NearestNeighbor(
	int dim,
	int num,
	data_t** sample,
	query_t* query,
	point_t<data_t>& NNpoint
	)
{
	double NNdist = DBL_MAX;
	int NNindex = -1;
	double distance;
	for (int n = 0; n < num; n++)
	{

		distance = Distance(dim, *sample++, query, NNdist);
		if (distance < NNdist)
		{
			NNdist = distance;
			NNindex = n;
		}
	}

	NNpoint.index = NNindex;
	NNpoint.distance = NNdist;
}

#endif// __NEAREST_NEIGHBOR__
