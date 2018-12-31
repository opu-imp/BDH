/**
* @file BDH.h
* @author Tomokazu Sato
* @date 2015/01/13
*/

#ifndef __BDH__
#define __BDH__

#include <list>
#include <map>
#include <queue>
#include <vector>
#include <float.h>
using namespace std;

#include <Subspace.h>
#include <HashTable.h>
#include <point.h>
#include <unordered_map> 
enum search_mode
{
	Radius,
	NumPoints,
	NumPoints2
};

/**
* @brief Bucket Distance Hashing
*/
template <typename data_t>
class BDH
{
public:
	typedef unsigned index_t;//!< type of index for point

private:

	int dim;			//!< dimension of dataspace
	int M;				//!< number of subspace
	int P;				//!< dimension of subspace
	int bit;			//!< bits num of hash table
	double delta;		//!< increment step of search radius for C search
	int subHashSizeMax; //!< max of sub hash size
	size_t pointSize;	//!< number of data points
	size_t entrySize;	//!< size of entory = sum of size of index and data point
	size_t hashSize;	//!< hash size = 2^bit
	double variance;

	Subspace* subspace;	//!< classes handling parameters of subspace
	Subspace  lestspace;//!< classe handling parameters of subspace not which construct the hash table

	HashTable hashTable;	//!< hash table

public:

	/**
	* @brief default constructor
	*/
	BDH()
		: dim(0)
		, M(0)
		, bit(0)
		, delta(0.0)
		, pointSize(0)
		, entrySize(0)
		, hashSize(0)
		, subspace(nullptr)
		, hashTable()
	{}

	~BDH()
	{
		delete[] subspace;
	}

	index_t get_nDdataPoints() const
	{
		return hashTable.get_nEntry();
	}

	double get_variance() const
	{
		return variance;
	}

	////////////////parameterTuning///////////////

	/**
	* @brief parameterTuning parameters for hashing
	*/
	void parameterTuning(
		int dim,				//!< [in] dimension of data space
		index_t num,			//!< [in] number of data points
		data_t** const data,	//!< [in] sample points for training
		base_t* const base,		//!< [in] base for projectToPCspace
		int M,					//!< [in] number of subspace
		int P,					//!< [in] dimension of subspace
		int bit,				//!< [in] bits num of hash table
		double bit_step = 1.0,	//!< [in] training parameter. 0 < bit_step <= bit.
		double sampling_rate = 1.0		//!< [in] training parameter.  0 < rate <= 1.
		);

	/**
	* @brief parameterTuning parameters for hashing
	*/
	void parameterTuning_ICCV2013(
		int dim,				//!< [in] dimension of data space
		index_t num,			//!< [in] number of data points
		data_t** const data,	//!< [in] sample points for training
		base_t* const base,		//!< [in] base for projectToPCspace
		int P,					//!< [in] dimension of subspace
		int bit,				//!< [in] bits num of hash table
		double bit_step = 1.0,	//!< [in] training parameter. 0 < bit_step <= bit.
		double sampling_rate = 1.0		//!< [in] training parameter.  0 < rate <= 1.
		);

	////////// store data points /////////////////

	/**
	* @brief store point set into hash table
	*/
	void storePoint(
		index_t num,	//!< [in] number of data points. 
		data_t**data	//!< [in] data point set. 
		);

	/**
	* @brief hash function
	* @return hash value
	*/
	size_t hashFunction(
		data_t* data	//!< [in] a point 
		);

	///////////////////// file handle ///////////////////////////////////

	/**
	* @brief save hash table
	* @return is file open ?
	*/
	bool saveTable(
		const string& path	//!< [in] file path
		);

	/**
	* @brief load hash table
	* @return is file open ?
	*/
	bool loadTable(
		const string& path	//!< [in] file path
		);

	/**
	* @brief save pareameters for hashing
	* @return is file open ?
	*/
	bool saveParameters(
		const string& path	//!< [in] file path
		);

	/**
	* @brief load pareameters for hashing
	* @return is file open ?
	*/
	bool loadParameters(
		const string& path	//!< [in] file path
		);


	///////////// Search Function ////////////////////////

	/**
	* @brief search in Bucket Distance R from query
	* @return number of points in search area
	*/
	int NearestNeighbor(
		data_t* query,
		point_t<data_t>* point,
		double searchParam,
		search_mode searchMode = NumPoints,
		int K = 1,
		double epsilon = DBL_MAX
		)const;

	/**
	* @brief search in Bucket Distance R from query
	* @return number of points in search area
	*/
	int BicromaticReverseNearestNeighbor(
		size_t nQuery,					//!< [in] number of query vecotr set
		data_t** query,					//!< [in] query point set
		vector<point_t<data_t>>* point,	//!< [out] result reverse nearest neighbors 
		double searchParam,
		search_mode searchMode,
		int K = 1,						//!< [in] number of nearest neighbors
		double epsilon = DBL_MAX		//!< [in] search points near than epsilon
		)const;

private:

	/**
	* @brief set Indexing Prameters
	*/
	void setParameters(
		const baseset_t* const baseSet, //!< [in] 
		const baseset_t& lestSet		//!< [in] 
		);

	void setLayerParam(
		layer_t* layer,	//!< [in]
		data_t* query	//!< [in]
		) const;

	int NearBucket_R(
		const double Radius,//探索半径
		layer_t* const layer,//クエリから求・ﾟたレイヤごとの部分距離情報
		const status_t& status,//ノードの状態を表す
		vector<hashKey_t>& bucketList //![out] collect hash key of buckets near than Radius from query
		) const;

	int NearBucket_C(
		const double& Lbound,//探索下限
		const double& Ubound,//探索上限
		layer_t* const layer,//クエリから求めたレイヤごとの部分距離情報
		const status_t& status,
		vector<hashKey_t>& bucketList
		) const;

	int NearBucket_C_list(
		const double Rbound,//探索半径
		layer_t* const layer,//クエリから求めたレイヤごとの部分距離情報
		list<status_t>& statusQue,//探索途中のノードを保持
		list<status_t>::iterator* itr,//ノードの状態を表す
		vector<hashKey_t>& bucketList
		) const;

	void linearSearchInNNcandidates(
		data_t* query,
		point_t<data_t>* point,
		int K,
		double epsilon,
		vector<hashKey_t>& bucketList
		) const;



	int searchInBucket(
		data_t* query,//クエリ
		size_t hashKey,
		priority_queue<point_t<data_t>>& NNpointQue
		) const;

	int getBucketList(
		data_t* query,
		double searchParam,
		search_mode searchMode,
		vector<hashKey_t>& bucketList
		)const;

	void getQueryListFromBucket(
		size_t nQuery,
		vector<hashKey_t>* bucketList,
		unordered_map<size_t, vector<size_t>>& queryListFromBucket
		)const;

	void getSubHashkey(size_t hashKey, int* subHashkey)const
	{
		for (int m = 0; m < M; ++m)
		{
			subHashkey[m] = hashKey % subspace[m].subHashSize;
			hashKey /= subspace[m].subHashSize;
		}
	}

	void extractRNNpoints(
		int K,
		vector<point_t<data_t>>* RKNNpoint,
		unordered_map<size_t, vector<size_t>>& bucketToQuery,
		data_t** query,
		double epsilon = DBL_MAX
		)const
	{

		address_t address;
		address_t address_end;
		size_t coll;
		bin_t bin;
		double dist, KNNdist;

		priority_queue<point_t<data_t>> NNqueries;

		point_t<data_t> dummy(-1, epsilon);
		point_t<data_t> pushPoint;

		//bucketToQuery内の各バケットにアクセス．
		for (auto btq_itr : bucketToQuery )
		{

			//バケット内のデータを取り出す
			hashTable.getBin(btq_itr.first, bin);
			coll = bin.collision;
			address = bin.addressOfChainList;
			address_end = address + entrySize*coll;

			for (; address != address_end; address += entrySize)
			{
				point_t<data_t> dataPoint(
					*reinterpret_cast<index_t*>(address + pointSize),
					reinterpret_cast<data_t*>(address));

				//K近傍点の初期化
				for (int i = 0; i < K; ++i)
				{
					NNqueries.push(dummy);
				}
				KNNdist = epsilon;

				for (auto ql_itr : btq_itr.second)
				{

					//距離計算
					dist = Distance(
						dim,
						dataPoint.addressOfpoint,
						query[ql_itr],
						KNNdist);

					//K近傍点の更新
					if (dist <= KNNdist)
					{
						NNqueries.pop();
						pushPoint.index = ql_itr;
						pushPoint.distance = dist;
						NNqueries.push(pushPoint);
						KNNdist = NNqueries.top().distance;
					}

				}

				//extract KNN queries
				while (NNqueries.empty() == false)
				{
					// is not top dummy query ?
					if (NNqueries.top().index != index_t(-1))
					{
						break;
					}
					NNqueries.pop();
				}
				while (NNqueries.empty() == false)
				{
					dataPoint.distance = NNqueries.top().distance;

					RKNNpoint[NNqueries.top().index].push_back(dataPoint);
					NNqueries.pop();
				}
			}
		}

	}
};





#endif //__BDH__