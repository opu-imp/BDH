/**
* @file Subspace.h
* @author Tomokazu Sato
* @date 2015/05/06
*/

#ifndef __SUBSPACE__
#define __SUBSPACE__

#include <algorithm>

#include "define.h"
#include "baseset.h"
#include <NearestNeighbor.h>

/**
* @brief manage subspace 
*/
class Subspace
{

public:
	static int dim;		 //!< dimension
	int subDim;			 //!< dimension of subspace
	int subHashSize;	 //!< hash size at subspace = 1<<bit
	double bit;			 //!< information volume
	double variance;	 //!< sum of variance
	double** base;		 //!< base direction[P][dim]
	size_t* hashKey;	 //!< hash value of bin corespond to centroid[subHashSize]
	double* cellVariance;//!< variance in cell[subHashSize]
	double** centroid;	 //!< centroid[subHashSize][subDim]

public:

	/**
	* @brief default constructor
	*/
	Subspace()
		: subDim(0)
		, subHashSize(0)
		, bit(0.0)
		, variance(0.0)
		, base(nullptr)
		, hashKey(nullptr)
		, cellVariance(nullptr)
		, centroid(nullptr)
	{}

	~Subspace()
	{
		if (base != nullptr)
		{
			for (int p = 0; p < subDim; ++p)
			{
				delete[] base[p];
			}
			delete[] base;
		}

		delete[] hashKey;

		delete[] cellVariance;

		if (centroid != nullptr)
		{
			for (int b = 0; b < subHashSize; ++b)
			{
				delete[] centroid[b];
			}
			delete[] centroid;
		}
	}

	/**
	* @brief initialize all member variables
	*/
	void clear();

	/**
	* @brief set training paramters
	*/
	void setParameters(
		const baseset_t& baseSet
		);

	/**
	* @brief inner product
	*/
	template<typename data_t>
	double innerProduct(
		double* base, 
		const data_t* data
		) const;

	/**
	* @brief project data into Principal Component space
	*/
	template<typename data_t>
	void getPCAdata(
		const data_t* data, 
		double* PCAdata) const;

	/**
	* @brief get sub hash value
	*/
	template<typename data_t>
	size_t getSubHashValue(
		const data_t* data
		) const;

	/**
	* @brief set node param
	*/
	template<typename data_t>
	void setNodeParam(
		node_t* node, 
		data_t* query
		);

	/**
	* @brief get distance to centroid
	*/
	double getDistanceToCentroid (
		double* PCAquery, 
		int centroidIndex
		)const;
};

template<typename data_t>
double Subspace::innerProduct(double* base, const data_t* data) const
{
	double val = 0.0;
	for (int d = 0; d < dim; ++d){
		val += base[d] * data[d];
	}
	return val;
}

template<typename data_t>
void Subspace::getPCAdata(const data_t* data, double* PCAdata) const
{
	double** base_p = base;
	double* PCAdata_end = PCAdata + subDim;
	for (; PCAdata != PCAdata_end; ++PCAdata)
	{
		*PCAdata = innerProduct(*base_p++, data);
	}
}


template<typename data_t>
size_t Subspace::getSubHashValue(
	const data_t* data
	) const
{

	//work space
	double* PCAdata = new double[subDim];

	getPCAdata(data, PCAdata);
	int idx = NearestNeighbor(subDim, subHashSize, centroid, PCAdata);
	delete[] PCAdata;

	return hashKey[idx];
}

template<typename data_t>
void Subspace::setNodeParam(
	node_t* node,
	data_t* query
	)
{

	double* PCAquery = new double[subDim];
	getPCAdata(query, PCAquery);

	for (int b = 0; b < subHashSize; ++b)
	{
		node[b].distance = getDistanceToCentroid(PCAquery, b);
		node[b].hashKey = hashKey[b];
	}
	sort(node, node + subHashSize);

	delete[] PCAquery;
}


#endif