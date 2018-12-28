/**
* @file    baseset.h
* @author  T.Sato
* @date    2015.05.04
* @version 1.0
*/

#ifndef __BASE_SET__
#define __BASE_SET__

/**
* @brief base infomation, mainly PCA base.
* @details destructor and copy constructor are not defined for fast sort.
*/
struct base_t
{
	static int dim;		//!< dimension of data space (static member)

	int idx;			//!< index of base
	double mean;		//!< mean at base direction
	double variance;	//!< variance at base direction
	double* direction;	//!< direction of base [dim]

	/**
	* @brief default constructor
	*/
	base_t()
		: idx(0)
		, mean(0.0)
		, variance(0.0)
		, direction(nullptr)
	{}

	/**
	* @brief initialize all attributes default and release memory
	*/
	void clear()
	{
		idx = 0;
		mean = 0.0;
		variance = 0.0;

		if (direction != nullptr)
		{
			delete[] direction;
		}

	}

	/**
	* @brief compare the variance
	*/
	bool operator < (
		const base_t& obj //!< [in] object
		)
	{
		return variance > obj.variance;
	}

};

/**
* @brief set of base_t
* @details destructor and copy constructor are not defined for fast sort
*/
struct baseset_t{

	int idx;				//!< index of base set
	int subDim;				//!< number of base
	double variance;		//!< sum of variance of base

	base_t* base;			//!< set of base_t[subDim]

	int k;					//!< number of centorid = 2^bit
	double bit;				//!< amount of infomation for this base set = log2(k)
	double error;			//!< error of quantization = sum of cellVariance
	double score;			//!< efficiency for incrementint the number of centroid
	double** centroid;		//!< centroid[k][subDim]
	double* cellVariance;	//!< variance in cell[k]

	/**
	* @brief default constructor
	*/
	baseset_t()
		: idx(0)
		, subDim(0)
		, variance(0.0)
		, base(nullptr)
		, k(0)
		, bit(0)
		, error(0.0)
		, score(0.0)
		, centroid(nullptr)
		, cellVariance(nullptr)
	{}

	/**
	* @brief initialize all attributes default and release memory
	*/
	void clear()
	{
		for (int sd = 0; sd < subDim; ++sd)
		{
			base[sd].clear();
		}

		if (centroid != nullptr)
		{
			for (int i = 0; i < k; ++i)
			{
				delete[] centroid[i];
			}
			delete[] centroid;
			centroid = nullptr;
		}

		if (cellVariance != nullptr)
		{
			delete[] cellVariance;
			cellVariance = nullptr;
		}

	}

	/**
	* @brief compare the variance
	*/
	bool operator < (
		const baseset_t& obj //!< [in] object
		)
	{
		return variance > obj.variance;
	}

};

#endif