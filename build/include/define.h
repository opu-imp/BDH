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
inline double NORM(double x)
{
	return x*x;
}

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
	double dist = 0.0;//!< distance
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

#endif//__define__