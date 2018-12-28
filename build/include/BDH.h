/**
* @file BDH.h
* @author Tomokazu Sato
* @date 2015/01/13
*/

#ifndef __BDH__
#define __BDH__

//#define NO_DISCAL

#include <list>
#include <map>
#include <queue>
#include <vector>
#include <unordered_map> 
#include <float.h>
using namespace std;

#include <Subspace.h>
#include <HashTable.h>
#include <point.h>

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

	/**
	* @brief destructor
	*/
	~BDH()
	{
		delete[] subspace;
	}

	size_t get_nDdataPoints() const
	{
		return hashTable.get_nEntry();
	}

	double get_variance() const
	{
		return variance;
	}

	void setCollisionMax(collision_t collisionMax)
	{
		hashTable.setCollisionMax(collisionMax);
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
	void storePoints(
		index_t num,	//!< [in] number of data points. 
		data_t**data	//!< [in] data point set. 
		);

	/**
	* @brief store point set into hash table
	*/
	void storePoints(
		index_t num,	//!< [in] number of data points. 
		index_t* index,	//!< [in] number of data points. 
		data_t**data	//!< [in] data point set. 
		);

	///////////////////// file handle ///////////////////////////////////

	/**
	* @brief データ登録済みのテーブルをファイルに保存する
	* @return is file open ?
	*/
	bool saveTable(
		const string& path	//!< [in] file path
		)const
	{
		return hashTable.writeTable(path);
	}

	/**
	* @brief dataをテーブルに登録した形でファイルに保存する
	* @return is file open ?
	*/
	bool saveTable(
		index_t  num,	//!< [in] number of data points.
		data_t** data,	//!< [in] data point set. size[num][dim].
		const string& path
		);

	/**
	* @brief dataをテーブルに登録した形でファイルに保存する．インデックスはindexの中のものを使う
	* @return is file open ?
	*/
	bool saveTable(
		index_t  num,	//!< [in] number of data points.
		index_t* index, //!< [in] index of data point. size[num].
		data_t** data,	//!< [in] data point set. size[num][dim].
		const string& path
		);

	/**
	* @brief load hash table
	* @return is file open ?
	*/
	bool loadTable(
		const string& path	//!< [in] file path
		)
	{
		return hashTable.readTable(path);
	}

	bool loadTableSkipData(
		const string& path	//!< [in] file path
		)
	{
		return hashTable.readTableSkipData(path);
	}

	/**
	* @brief save pareameters for hashing
	* @return is file open ?
	*/
	bool saveParameters(
		const string& path	//!< [in] file path
		) const;

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
	size_t BicromaticReverseNearestNeighbor(
		int nQuery,					//!< [in] number of query vecotr set
		data_t** query,					//!< [in] query point set
		vector<point_t<data_t>>* point,	//!< [out] result reverse nearest neighbors 
		double searchParam,
		search_mode searchMode,
		int K = 1,						//!< [in] number of nearest neighbors
		double epsilon = DBL_MAX		//!< [in] search points near than epsilon
		)const;

	//////////////////// Incremental Search  ///////////////////////////////

	/**
	* @brief インクリメンタルサーチに必要なパラメータの初期化．
	*/
	void InitializeForBRNN_Incremental(param_for_incremental_search<data_t>& pfis);

	/**
	* @brief インクリメンタルバイクロマティック逆最近傍探索．
	*        初期化済みのpfisを入力すること．
	*        searchParamを徐々に大きくしていく 
	*/
	size_t BicromaticReverseNearestNeighbor_Incremental(
		param_for_incremental_search<data_t>& pfis,
		double searchParam,
		vector<point_t<data_t>>* RNNpoint
		)const;

private:

	/**
	* @brief hash function
	* @return hash value
	*/
	size_t hashFunction(
		data_t* data	//!< [in] a point 
		)const;



	/**
	* @brief set Indexing Prameters
	*/
	void setParameters(
		const baseset_t* const baseSet, //!< [in] 
		const baseset_t& lestSet		//!< [in] 
		);

	void calcLayerParam(
		layer_t*& layer,	//!< [in]
		data_t* const query	//!< [in]
		) const;

	int NearBucket_R(
		const double& Radius,//探索半径
		layer_t* const& layer,//クエリから求めたレイヤごとの部分距離情報
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

	int calcBucketList(
		data_t* query,
		double searchParam,
		search_mode searchMode,
		vector<hashKey_t>& bucketList
		)const;

	void addQueryListFromBucket(
		int queryIndex,
		const vector<hashKey_t>& bucketList,
		unordered_map<size_t, vector<size_t>*>& queryListFromBucket
		)const;

	void getSubHashkey(size_t hashKey, int* subHashkey)const;

	int extractRNNpoints(
		int K,
		vector<point_t<data_t>>* RKNNpoint,
		data_t** query,
		unordered_map<size_t, vector<size_t>*>& bucketToQuery,
		double epsilon = DBL_MAX
		)const ;

	int extractNNqueries(
		param_for_incremental_search<data_t>& pfis,
		unordered_map<size_t, vector<size_t>*>& bucketToQuery
		)const;

	int NearBucket_R(
		const double Radius,//! 探索半径
		layer_t* const layer,//! クエリから求めたレイヤごとの部分距離情報
		const status_t& status,//! ノードの状態を表す
		data_t* query,//! クエリ
		priority_queue<point_t<data_t>> NNpointQue //! K近傍点を保持するコンテナ 
		)const;

	int NearBucket_C(
		const double& Lbound,//探索下限
		const double& Ubound,//探索上限
		layer_t* const& layer,//クエリから求めたレイヤごとの部分距離情報
		const status_t& status,
		data_t* const& query,//! クエリ
		priority_queue<point_t<data_t>>& NNpointQue //! K近傍点を保持するコンテナ 
		) const;


	unsigned nearBucket(
		data_t* query,
		layer_t* layer,
		double lestSpaceDist,
		double searchParam,
		search_mode searchMode,
		unsigned& NNC,
		double& preRadius,
		vector<hashKey_t>& bucketList
		)const;
};

template <typename data_t>
inline int BDH<data_t>::NearBucket_R(
	const double& Radius,//探索半径
	layer_t* const& layer,//クエリから求めたレイヤごとの部分距離情報
	const status_t& status,//ノードの状態を表す
	vector<hashKey_t>& bucketList //![out] collect hash key of buckets near than Radius from query
	)const
{
	const int m_plus1 = status.m + 1;

	int count = 0;

	if (m_plus1 == M)
	{
		size_t hashKey;
		index_t collision;
		node_t* node = layer[status.m].node;
		const double layerBound = Radius - status.dist;
		for (; node->distance <= layerBound; ++node)
		{
			hashKey = status.hashKey + node->hashKey;

			collision = hashTable.getCollision(hashKey);
			if (collision > 0)
			{
				count += collision;
				bucketList.push_back(
					hashKey_t(hashKey, status.dist + node->distance)
					);
			}
		}
	}
	else
	{
		status_t statusNext(m_plus1);
		node_t* node = layer[status.m].node;

		const double layerBound = Radius - status.dist - layer[m_plus1].restMin;
		for (; node->distance <= layerBound; ++node)
		{
			statusNext.hashKey = status.hashKey + node->hashKey;
			statusNext.dist = status.dist + node->distance;
			count += NearBucket_R(
				Radius,
				layer,
				statusNext,
				bucketList
				);
		}
	}

	return count;
}

template <typename data_t>
inline int BDH<data_t>::NearBucket_C(
	const double& Lbound,//探索下限
	const double& Ubound,//探索上限
	layer_t* const layer,//クエリから求めたレイヤごとの部分距離情報
	const status_t& status,
	vector<hashKey_t>& bucketList
	) const
{

	const int m_plus1 = status.m + 1;
	int count = 0;

	if (m_plus1 == M)
	{
		size_t hashKey;
		hashKey_t Key;
		address_t bucket_p;

		const double layerLowerBound = Lbound - status.dist;
		const double layerUpperBound = Ubound - status.dist;
		node_t* node = layer[status.m].node;
		for (; node->distance <= layerLowerBound; ++node){}
		for (; node->distance <= layerUpperBound; ++node)
		{
			hashKey = status.hashKey + node->hashKey;
			bucket_p = hashTable.getPointer(hashKey);
			if (bucket_p)
			{
				Key.setVariable(hashKey, status.dist + node->distance);
				bucketList.push_back(Key);
				count += *reinterpret_cast<collision_t*>(bucket_p);
			}
		}
	}
	else
	{
		const double layerLowerBound = Lbound - status.dist - layer[m_plus1].restMax;
		const double layerUpperBound = Ubound - status.dist - layer[m_plus1].restMin;

		status_t statusNext(m_plus1);

		node_t* node = layer[status.m].node;
		for (; node->distance <= layerLowerBound; ++node){}
		for (; node->distance <= layerUpperBound; ++node)
		{
			statusNext.hashKey = status.hashKey + node->hashKey;
			statusNext.dist = status.dist + node->distance;

			count += NearBucket_C(
				Lbound, Ubound,
				layer, statusNext, bucketList
				);

		}

	}

	return count;
}

template <typename data_t>
inline int BDH<data_t>::NearBucket_C_list(
	const double Rbound,//探索半径
	layer_t* const layer,//クエリから求めたレイヤごとの部分距離情報
	list<status_t>& statusQue,//探索途中のノードを保持
	list<status_t>::iterator* itr,//ノードの状態を表す
	vector<hashKey_t>& bucketList
	) const
{

	const int m = (*itr)->m;
	const int m_plus1 = m + 1;
	node_t* const node = layer[m].node;

	int count = 0;
	double layerBound = Rbound - (*itr)->dist;
	int i = (*itr)->nodeIdx;
	if (m_plus1 == M)
	{
		size_t hashKey;
		index_t collision;

		for (; node[i].distance <= layerBound; ++i)
		{
			hashKey = (*itr)->hashKey + node[i].hashKey;
			collision = hashTable.getCollision(hashKey);
			if (collision > 0)
			{
				bucketList.push_back(
					hashKey_t(hashKey, (*itr)->dist + node[i].distance)
					);
				count += collision;
			}
		}
	}
	else
	{
		layerBound -= layer[m_plus1].restMin;

		status_t statusNext(m_plus1);
		list<status_t>::iterator itr_next;
		for (; node[i].distance <= layerBound; ++i){

			statusNext.hashKey = (*itr)->hashKey + node[i].hashKey;
			statusNext.dist = (*itr)->dist + node[i].distance;
			statusQue.push_front(statusNext);

			itr_next = statusQue.begin();
			count += NearBucket_C_list(
				Rbound, layer, statusQue, &itr_next, bucketList
				);
		}
	}

	//すべてのノードにアクセスしたか
	if (i == layer[m].k)
	{
		//ノードを消して次へ
		statusQue.erase((*itr)++);
	}
	else
	{
		//ノードの状態を更新して次へ
		(*itr)->nodeIdx = i;
		++(*itr);
	}

	return count;
}

template <typename data_t>
inline int BDH<data_t>::searchInBucket(
	data_t* query,//クエリ
	size_t hashKey,
	priority_queue<point_t<data_t>>& NNpointQue
	) const
{

	bin_t bin;
	hashTable.getBin(hashKey, bin);
	collision_t coll = bin.collision;
	if (coll == 0)
	{
		return 0;
	}

	/*set point's pointer*/
	address_t addr = bin.addressOfChainList;
	address_t addr_end = addr + coll*entrySize;
	double KNNdist = NNpointQue.top().distance;

	double dist;
	data_t* data;
	point_t<data_t> pushData;

	while (addr != addr_end)
	{
		/*Distance Calculation*/
		data = reinterpret_cast<data_t*>(addr);
		dist = Distance(dim, data, query, KNNdist);

		if (dist < KNNdist)
		{
			NNpointQue.pop();
			pushData.setMemberVariable(*reinterpret_cast<index_t*>(addr + pointSize), data, dist);
			NNpointQue.push(pushData);

			KNNdist = NNpointQue.top().distance;
		}
		addr += entrySize;
	}

	return coll;
}

template <typename data_t>
inline int BDH<data_t>::NearBucket_R(
	const double Radius,//! 探索半径
	layer_t* const layer,//! クエリから求めたレイヤごとの部分距離情報
	const status_t& status,//! ノードの状態を表す
	data_t* query,//! クエリ
	priority_queue<point_t<data_t>> NNpointQue //! K近傍点を保持するコンテナ 
	)const
{
	const int m_plus1 = status.m + 1;
	int count = 0;

	if (m_plus1 == M)
	{//!最終レイヤでの処理
		size_t hashKey;
		node_t* node = layer[status.m].node;
		const double layerBound = Radius - status.dist;
		for (; node->distance <= layerBound; ++node)
		{//バケット距離が指定したパラメータ以下である限り探索を継続

			//ハッシュキーを取得
			hashKey = status.hashKey + node->hashKey;
			//ハッシュキーに対応するバケットに登録されている点と距離計算を行い，NNpointQueに格納されている最近傍点を更新する
			searchInBucket(query, hashKey, NNpointQue);
		}
	}
	else
	{
		status_t statusNext(m_plus1);
		node_t* node = layer[status.m].node;

		const double layerBound = Radius - status.dist - layer[m_plus1].restMin;
		for (; node->distance <= layerBound; ++node)
		{
			statusNext.hashKey = status.hashKey + node->hashKey;
			statusNext.dist = status.dist + node->distance;
			count += NearBucket_R(
				Radius, layer, statusNext,
				query, NNpointQue
				);
		}
	}

	return count;
}

template <typename data_t>
inline int BDH<data_t>::NearBucket_C(
	const double& Lbound,//探索下限
	const double& Ubound,//探索上限
	layer_t* const& layer,//クエリから求めたレイヤごとの部分距離情報
	const status_t& status,
	data_t* const& query,//! クエリ
	priority_queue<point_t<data_t>>& NNpointQue //! K近傍点を保持するコンテナ 
	) const
{

	const int m_plus1 = status.m + 1;
	int count = 0;

	if (m_plus1 == M)
	{
		size_t hashKey;

		const double layerLowerBound = Lbound - status.dist;
		const double layerUpperBound = Ubound - status.dist;
		node_t* node = layer[status.m].node;
		for (; node->distance <= layerLowerBound; ++node){}
		for (; node->distance <= layerUpperBound; ++node)
		{
			hashKey = status.hashKey + node->hashKey;
			count += searchInBucket(query, hashKey, NNpointQue);
		}
	}
	else
	{
		const double layerLowerBound = Lbound - status.dist - layer[m_plus1].restMax;
		const double layerUpperBound = Ubound - status.dist - layer[m_plus1].restMin;

		status_t statusNext(m_plus1);

		node_t* node = layer[status.m].node;
		for (; node->distance <= layerLowerBound; ++node){}
		for (; node->distance <= layerUpperBound; ++node)
		{
			statusNext.hashKey = status.hashKey + node->hashKey;
			statusNext.dist = status.dist + node->distance;

			count += NearBucket_C(
				Lbound, Ubound,
				layer, statusNext, query, NNpointQue
				);

		}

	}

	return count;
}

#endif //__BDH__