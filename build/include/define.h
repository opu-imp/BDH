/**
* @file		define.h
* @author	Tomokazu Sato
* @date		2015/05/05
*/

#ifndef __DEFINE__
#define __DEFINE__

#include <iostream>
using namespace std;

/**
* @brief norm for distance
* @return scalar distance
*/
static inline double NORM(const double& x)
{
	return x*x;
}

//#define NORM(x) ((x)*(x))

/**
* @brief type of node. a node coresponde to a centroid
*/
struct node_t{
	size_t hashKey;//!< hash value
	double distance; //!< sub bucket distance

	/**
	* @brief compare the distance
	* @return is which lesser ?
	*/
	bool operator < (
		const node_t& obj //!< compared object
		)
	{
		return distance < obj.distance;
	}
};

/**
* @brief type of layer. a node coresponde to a subspace
*/
struct layer_t
{
	int k;			//!< number of nodes
	double restMin;	//!< the minimam rest distance from this layer
	double restMax;	//!< the maximam rest distance from this layer
	double gap;		//!< the gap of distance between max and min
	node_t* node;	//!< nodes.

	void calc_gap()
	{
		gap	= node[k - 1].distance - node[0].distance;
	}

	/**
	* @brief compare the gap
	* @return is which gap larger ?
	*/
	bool operator < (const layer_t& obj)
	{
		return gap > obj.gap;
	}
};

/**
* @brief status of neare bucket search
*/
struct status_t{

	int m;			//!< index of layer
	int nodeIdx;	//!< index of nodes
	double dist;	//!< distance
	size_t hashKey;	//!< hash value

	/**
	* @brief default constructor
	*/
	status_t()
		: m(0)
		, nodeIdx(0)
		, dist(0.0)
		, hashKey(0)
	{}

	/**
	* @brief constructor
	*/
	status_t(
		const int& m,			//!< index of layer
		const int& nodeIdx,		//!< index of nodes
		const size_t& hashKey,	//!< hash value
		const double& dist)		//!< distance
		: m(m)
		, nodeIdx(nodeIdx)
		, dist(dist)
		, hashKey(hashKey)
	{}

	/**
	* @brief constructor
	*/
	status_t(
		const int& m			//!< index of layer
		)
		: m(m)
		, nodeIdx(0)
		, dist(0.0)
		, hashKey(0)
	{}
};

/**
* @brief status of neare bucket search
*/
struct hashKey_t{

	size_t hashKey;	//!< hash value
	double dist;	//!< distance

	/**
	* @brief default constructor
	*/
	hashKey_t()
		:hashKey(0)
		, dist(0)
	{}

	/**
	* @brief constructor. allocate memory of nodeIdx and deep copied.
	*/
	hashKey_t(
		size_t hashKey,	//!< hash value
		double dist)	//!< distance
		: hashKey(hashKey)
		, dist(dist)
	{}

	void setVariable(size_t hashKey, double dist)
	{
		this->hashKey = hashKey;
		this->dist = dist;
	}
};


enum search_mode
{
	Radius,
	NumPoints,
	NumPoints2
};

struct dist_index_pair_t
{
	unsigned index;
	float dist;

	dist_index_pair_t()
	{
		this->index = (unsigned)(-1);
		this->dist = FLT_MAX;
	}

	dist_index_pair_t(const unsigned& index, const float& dist)
	{
		this->index = index;
		this->dist = dist;
	}


	/**
	* @brief compare the distance
	* @return is my distance lessor than e's ?
	*/
	bool operator <(
		const dist_index_pair_t &e	//!< [in] compare object
		) const
	{
		return dist < e.dist;
	}

	/**
	* @brief compare the distance
	* @return is my distance equal to e's ?
	*/
	bool operator ==(
		const dist_index_pair_t &e	//!< [in] compare object
		) const
	{
		return index == e.index;
	}
};

template<typename data_t>
struct param_for_incremental_search
{
	//static parameters
	const size_t nData;
	const unsigned nQuery;					//!< [in] number of query vecotr set
	data_t** const query;					//!< [in] query point set
	const search_mode searchMode;
	const double epsilon;		//!< [in] search points near than epsilon
	int M;
	layer_t** layer;
	double* lestSpaceDist;

	//incremental parameters
	dist_index_pair_t* NNquery;
	unsigned* NNC;
	double* preRadius;

	param_for_incremental_search(
		size_t nData,
		int nQuery,					//!< [in] number of query vecotr set
		data_t** const query,					//!< [in] query point set
		search_mode searchMode,
		float epsilon = FLT_MAX		//!< [in] search points near than epsilon
		) :
		nData(nData),
		nQuery(nQuery),
		query(query),
		searchMode(searchMode),
		epsilon(epsilon)
	{
		NNquery = new dist_index_pair_t[nData];
		preRadius = new double[nQuery];
		NNC = new unsigned[nQuery];
		layer = new layer_t*[nQuery];
		lestSpaceDist = new double[nQuery];

		memset(NNC, 0, sizeof(unsigned)*nQuery);
		for (int q = 0; q < nQuery; ++q)
		{
			NNquery[q].dist = epsilon;
			NNquery[q].index = UINT_MAX;
			preRadius[q] = 0.0;
		}
	}

	~param_for_incremental_search()
	{
		delete[] NNquery;
		delete[] preRadius;
		delete[] NNC;
		delete[] lestSpaceDist;
		
		for (unsigned q = 0; q < nQuery; ++q)
		{
			for (int m = 0; m < M; ++m)
			{
				delete[] layer[q][m].node;
			}
			delete[] layer[q];
		}
		delete[] layer;

	}
};



#endif//__define__