/**
* @file BDHtraining.h
* @author Tomokazu Sato
* @date 2015/05/04
*/

#ifndef __BDH_TRAINING__
#define __BDH_TRAINING__

#include <baseset.h>
#include <string>
using namespace std;

#include "k_means.h"

/**
* @brief training BDH parameters
*/
template <typename data_t>
class BDHtraining
{
	//for data sample
	int dim;	//!< dimension of dataspace
	unsigned num;	//!< number of data samples

	//parameters
	int M;	//!< number of subspace
	int P;	//!< dimension of subspace
	int U;	//!< numer of base used for hashing = P*M
	int bit;//!< bits num of hash table

	//result values
	size_t hashSize;	//!< hash size
	baseset_t* baseSet;	//!< baseSet[M]
	baseset_t lestSet;	//!< baseSet not used for hashing

public:

	/**
	* @brief default constructor
	*/
	BDHtraining()
		: dim(0)
		, num(0)
		, M(0)
		, P(0)
		, U(0)
		, bit(0)
		, hashSize(0)
		, baseSet(nullptr)
		, lestSet()
	{}

	/**
	* @brief destructor
	*/
	~BDHtraining()
	{
		if (baseSet != nullptr)
		{
			for (int m = 0; m < M; ++m)
			{
				baseSet[m].clear();
			}
			delete[] baseSet;
			baseSet = nullptr;
		}
		lestSet.clear();
	}

	int getM()
	{
		return M;
	}

	/**
	* @brief training BDH parameters
	*/
	void training(
		int dim,						//!< [in] dimension
		unsigned num,						//!< [in] number of sample
		data_t** data,					//!< [in] sample data set
		const base_t* const baseInput,	//!< [in] base
		int M,							//!< [in] number of subspace
		int P,							//!< [in] dimension of subspace
		int bit,						//!< [in] bits num of hash table
		double bit_step = 1.0			//!< [in] training parameter.
		);

	/**
	* @brief training BDH parameters
	*/
	void training_ICCV2013(
		int dim,						//!< [in] dimension
		unsigned num,					//!< [in] number of sample
		data_t** data,					//!< [in] sample data set
		const base_t* const baseInput,	//!< [in] base
		int M,							//!< [in] number of subspace
		int bit,						//!< [in] bits num of hash table
		double bit_step = 1.0			//!< [in] training parameter.
		);

	/**
	* @brief getter
	* @return baseSet
	*/
	const baseset_t* const getBaseSet()
	{
		return baseSet; 
	}

	/**
	* @brief training BDH parameters
	* @return lsetSet
	*/
	const baseset_t& getLestSet()
	{
		return lestSet;
	}

	/**
	* @brief training BDH parameters
	* @return is file open
	*/
	bool saveParameters(
		const string& path//!< [in] file path
		);

private:

	//データ空間をM個の部分空間に分割する
	void partitioningDataspace(
		const base_t* const base
		);

	/**
	*上位の基底セットに上位の主成分をまとめて登録する
	*
	*/
	void partitioningDataspace_ICCV2013(
		const base_t* const base
		);


	//各部分空間のセントロイドを求める
	//セントロイドの数は自動的に求められる
	void calclateCentroid(
		float*** subPrjData,
		double bit_step
		);

	//各部分空間のセントロイドを求める
	//セントロイドの数は自動的に求められる
	void calclateCentroid_ICCV2013(
		float*** subPrjData,
		double bit_step
		);

	void updateCentroid(
		double bit_step,
		float*** subPrjData,
		K_Means<float, double>*& k_means);

	//推定距離の計算用に各セル内の分散を求める
	void calculateCellVariance(
		float*** subPrjData
		);

	//内積
	double innerProduct(const double* base, const data_t* data){
		double val = 0.0;
		for (int d = 0; d < dim; ++d){
			val += base[d] * data[d];
		}
		return val;
	}


};

#endif