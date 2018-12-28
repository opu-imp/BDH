#include <BDHtraining.h>
#include <BDH.h>
const double deltaRate = 50;
#include <omp.h>


/************** Indexing *******************/
template <typename data_t>
void BDH<data_t>::parameterTuning(
	int dim,
	index_t num,
	data_t** const data,
	base_t* const base,
	int M,
	int P,
	int bit,
	double bit_step,
	double sampling_rate
	)
{

	this->dim = dim;
	this->M = M;
	this->P = P;
	this->bit = bit;
	hashSize = (size_t(1) << bit);
	Subspace::dim = dim;

	pointSize = sizeof(data_t)*dim;
	entrySize = pointSize + sizeof(index_t);

	variance = 0;
	for (int d = 0; d < dim; ++d)
	{
		variance += base[d].variance;
	}
	delta = variance / deltaRate;

	hashTable.initialize(entrySize, hashSize);

	BDHtraining<data_t> bdh_learning;

	if (sampling_rate < 1.0)
	{
		vector<data_t*> data_train;
		double step = 1.0 / sampling_rate;
		for (double n = 0; n < num; n += step)
		{
			data_train.push_back(data[(size_t)n]);
		}

		bdh_learning.training(
			dim, (index_t)data_train.size(), data_train.data(),
			base, M, P, bit, bit_step
			);

	}
	else
	{

		bdh_learning.training(
			dim, num, data,
			base, M, P, bit, bit_step
			);

	}

	setParameters(
		bdh_learning.getBaseSet(), bdh_learning.getLestSet()
		);

}

/************** Indexing *******************/
template <typename data_t>
void BDH<data_t>::parameterTuning_ICCV2013(
	int dim,
	index_t num,
	data_t** const data,
	base_t* const base,
	int P,
	int bit,
	double bit_step,
	double sampling_rate
	)
{

	this->dim = dim;
	this->P = P;
	this->bit = bit;
	hashSize = (size_t(1) << bit);
	Subspace::dim = dim;

	pointSize = sizeof(data_t)*dim;
	entrySize = pointSize + sizeof(index_t);

	variance = 0;
	for (int d = 0; d < dim; ++d)
	{
		variance += base[d].variance;
	}
	delta = variance / deltaRate;

	hashTable.initialize(entrySize, hashSize);

	BDHtraining<data_t> BDHtrainer;

	if (sampling_rate < 1.0)
	{
		data_t** l_data = new data_t*[num];
		index_t l_num = 0;

		double tmp = 0.0;
		for (size_t n = 0; n < num; ++n)
		{
			tmp += sampling_rate;
			if (sampling_rate >= 1.0)
			{
				l_data[l_num++] = data[n];
				sampling_rate -= 1.0;
			}
		}

		BDHtrainer.training_ICCV2013(
			dim, l_num, l_data,
			base, P, bit, bit_step
			);

		delete[] l_data;

	}
	else
	{

		BDHtrainer.training_ICCV2013(
			dim, num, data,
			base, P, bit, bit_step
			);

	}

	M = BDHtrainer.getM();
	setParameters(
		BDHtrainer.getBaseSet(), BDHtrainer.getLestSet()
		);

}


template <typename data_t>
void BDH<data_t>::setParameters(
	const baseset_t* const baseSet, const baseset_t& lestSet
	)
{
	subspace = new Subspace[M];
	size_t rank = 1;

	for (int m = 0; m < M; ++m)
	{
		subspace[m].setParameters(baseSet[m]);

		subspace[m].hashKey = new size_t[subspace[m].subHashSize];

		for (int i = 0; i < subspace[m].subHashSize; ++i)
		{
			subspace[m].hashKey[i] = rank*i;

		}

		rank *= subspace[m].subHashSize;
	}

	lestspace.subDim = lestSet.subDim;
	lestspace.variance = lestSet.variance;
	lestspace.centroid = new double*[1];
	lestspace.centroid[0] = new double[lestSet.subDim];

	lestspace.base = new double*[lestspace.subDim];
	for (int d = 0; d < lestspace.subDim; ++d)
	{
		lestspace.centroid[0][d] = lestSet.base[d].mean;
		lestspace.base[d] = new double[dim];
		memcpy(lestspace.base[d], lestSet.base[d].direction, sizeof(double)*dim);
	}
	lestspace.cellVariance = new double[1];
	lestspace.cellVariance[0] = lestspace.variance;

}

template <typename data_t>
bool BDH<data_t>::saveParameters(const string& path) const
{

	ofstream ofs(path);
	if (ofs.is_open() == false)
	{
		return false;
	}

	ofs << dim << "\t"
		<< M << "\t"
		<< P << "\t"
		<< bit << "\t"
		<< delta << "\t"
		<< pointSize << "\t"
		<< entrySize << "\t"
		<< hashSize << "\t"
		<< variance << endl;

	for (int m = 0; m < M; ++m)
	{
		ofs << subspace[m].subHashSize << "\t"
			<< subspace[m].variance << endl;

		for (int sd = 0; sd < P; ++sd)
		{
			for (int d = 0; d < dim; ++d)
			{
				ofs << subspace[m].base[sd][d] << "\t";
			}
			ofs << endl;
		}

		for (int i = 0; i < subspace[m].subHashSize; ++i)
		{

			ofs << subspace[m].cellVariance[i] << "\t" << subspace[m].hashKey[i] << endl;

			for (int d = 0; d < subspace[m].subDim; ++d)
			{
				ofs << subspace[m].centroid[i][d] << "\t";
			}

			ofs << endl;
		}
	}

	ofs << lestspace.subDim << "\t"
		<< lestspace.variance << endl;
	for (int sd = 0; sd < lestspace.subDim; ++sd)
	{
		ofs << lestspace.centroid[0][sd] << "\t";
	}
	ofs << endl;

	for (int sd = 0; sd < lestspace.subDim; ++sd)
	{
		for (int d = 0; d < dim; ++d)
		{
			ofs << lestspace.base[sd][d] << "\t";
		}
		ofs << endl;
	}


	ofs.close();
	return true;
}

template <typename data_t>
bool BDH<data_t>::loadParameters(
	const string& path
	)
{

	ifstream ifs(path);
	if (ifs.is_open() == false)
	{
		cerr << "failed to open " << path << endl;
		return false;
	}

	ifs >> dim
		>> M
		>> P
		>> bit
		>> delta
		>> pointSize
		>> entrySize
		>> hashSize
		>> variance;

	Subspace::dim = dim;

	subHashSizeMax = 0;
	subspace = new Subspace[M];
	for (int m = 0; m < M; ++m)
	{
		subspace[m].subDim = P;

		ifs >> subspace[m].subHashSize
			>> subspace[m].variance;

		if (subspace[m].subHashSize > subHashSizeMax)
		{
			subHashSizeMax = subspace[m].subHashSize;
		}

		subspace[m].base = new double*[P];
		for (int sd = 0; sd < P; ++sd)
		{
			subspace[m].base[sd] = new double[dim];
			for (int d = 0; d < dim; ++d)
			{
				ifs >> subspace[m].base[sd][d];
			}
		}

		subspace[m].cellVariance = new double[subspace[m].subHashSize];
		subspace[m].hashKey = new size_t[subspace[m].subHashSize];
		subspace[m].centroid = new double*[subspace[m].subHashSize];
		for (int i = 0; i < subspace[m].subHashSize; ++i)
		{
			ifs >> subspace[m].cellVariance[i] >> subspace[m].hashKey[i];

			subspace[m].centroid[i] = new double[P];
			for (int d = 0; d < P; ++d)
			{
				ifs >> subspace[m].centroid[i][d];
			}
		}
	}

	ifs >> lestspace.subDim
		>> lestspace.variance;

	lestspace.centroid = new double*[1];
	lestspace.centroid[0] = new double[lestspace.subDim];
	for (int sd = 0; sd < lestspace.subDim; ++sd)
	{
		ifs >> lestspace.centroid[0][sd];
	}
	lestspace.cellVariance = new double[1];
	lestspace.cellVariance[0] = lestspace.variance;

	lestspace.base = new double*[lestspace.subDim];
	for (int sd = 0; sd < lestspace.subDim; ++sd)
	{
		lestspace.base[sd] = new double[dim];
		for (int d = 0; d < dim; ++d)
		{
			ifs >> lestspace.base[sd][d];
		}
	}

	ifs.close();

	hashTable.initialize(entrySize, hashSize);
	return true;
}


template <typename data_t>
void BDH<data_t>::storePoints(index_t num, data_t** data)
{
	//alloc workspace
	collision_t* collision = new collision_t[hashSize];
	memset(collision, 0, sizeof(collision_t)*hashSize);

	size_t* hashKey = new size_t[num];
	for (index_t n = 0; n < num; ++n)
	{
		//get hash value
		hashKey[n] = hashFunction(data[n]);
		//increment collision
		++collision[hashKey[n]];
	}

	hashTable.allocTable(collision);
	delete[] collision;

	char* entry = new char[entrySize];
	if (entry == nullptr)
	{
		exit(__LINE__);
	}

	for (index_t n = 0; n < num; ++n)
	{
		memcpy(entry, data[n], pointSize);
		*reinterpret_cast<index_t*>(entry + pointSize) = n;

		hashTable.storeEntryWithoutAlloc(hashKey[n], entry);
	}
	delete[] entry;

	delete[] hashKey;
}

template <typename data_t>
void BDH<data_t>::storePoints(
	index_t  num,	//!< [in] number of data points.
	index_t* index, //!< [in] index of data point
	data_t** data	//!< [in] data point set. 
	)
{
	//alloc workspace
	collision_t* collision;
	try
	{
		collision = new collision_t[hashSize];
	}
	catch (bad_alloc)
	{ // 例外 bad_alloc をここで受け取る
		cerr << "bad_alloc for collision" << endl;
		abort();
	}
	catch (...)
	{ // それ以外の例外はここで受け取る
		cerr << "iregular for collison" << endl;
		abort();
	}
	//初期化
	memset(collision, 0, sizeof(collision_t)*hashSize);

	cout << "calculate hash values." << endl;
	const int int_num = num;
	size_t* hashKey;
	try
	{
		hashKey = new size_t[num];
	}
	catch (bad_alloc)
	{ // 例外 bad_alloc をここで受け取る
		cerr << "bad_alloc for collision" << endl;
		abort();
	}
	catch (...)
	{ // それ以外の例外はここで受け取る
		cerr << "iregular for collison" << endl;
		abort();
	}

	//各データ点のハッシュ値を計算しコリジョンを取得
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (int n = 0; n < int_num; ++n)
	{
		//get hash value
		hashKey[n] = hashFunction(data[n]);

		//increment collision
#ifdef _OPENMP
#pragma omp atomic 
#endif
		++collision[hashKey[n]];

	}

	hashTable.allocTable(collision);
	delete[] collision;

	cout << "store data points." << endl;
	char* entry = new char[entrySize];
	for (index_t n = 0; n < num; ++n)
	{
		memcpy(entry, data[n], pointSize);
		*reinterpret_cast<index_t*>(entry + pointSize) = index[n];
		hashTable.storeEntryWithoutAlloc(hashKey[n], entry);
	}
	delete[] entry;

	delete[] hashKey;
}


template <typename data_t>
inline size_t BDH<data_t>::hashFunction(data_t* data)const
{

	size_t hashKey = 0;
	for (int m = 0; m < M; ++m)
	{
		hashKey += subspace[m].getSubHashValue(data);
	}
	return hashKey;
}

/************* file handle ***************/

template <typename data_t>
bool BDH<data_t>::saveTable(
	index_t  num,	//!< [in] number of data points.
	data_t** data,	//!< [in] data point set. 
	const string& path
	)
{
	ofstream ofs(path, ios::binary);
	if (ofs.fail())
	{
		ofs.close();
		cerr << "failed to open " << path << endl;
		return false;
	}

	//alloc workspace
	size_t* hashValue = new size_t[num];

	long long longNum = num;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (long long n = 0; n < longNum; ++n)
	{
		hashValue[n] = hashFunction(data[n]);
	}

	cout << "calculate hash values." << endl;
	//各データ点のハッシュ値を計算しコリジョンを取得
	vector<vector<index_t>> IDtable(hashSize);
	for (index_t n = 0; n < num; ++n)
	{
		IDtable[hashValue[n]].push_back(n);
	}
	delete[] hashValue;

	cout << "store data points." << endl;
	collision_t* collision = new collision_t[hashSize];
	long long longHashSize = hashSize;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (long long h = 0; h < longHashSize; ++h)
	{
		collision[h] = (collision_t)IDtable[h].size();
	}

	ofs << entrySize << "\t" << hashSize << "\t" << num;
	ofs.write((char*)collision, sizeof(collision_t)*hashSize);
	delete[] collision;

	vector<char> entry;
	char* entry_p;
	size_t progressStep = hashSize / 20;
	size_t progressPoint = progressStep;
	for (size_t h = 0; h < hashSize; ++h)
	{
		if (IDtable[h].empty() == false)
		{
			entry.resize(entrySize*IDtable[h].size());
			entry_p = entry.data();

			// write hash chain structure on RAM
			for (auto& obj : IDtable[h])
			{
				memcpy(entry_p, data[obj], pointSize);
				*reinterpret_cast<index_t*>(entry_p + pointSize) = obj;
				entry_p += entrySize;
			}

			// write on HD
			ofs.write(entry.data(), entry.size());
		}
		else if (h > progressPoint)
		{
			cout << h * 100 / hashSize << " % done." << endl;
			progressPoint += progressStep;
		}
	}
	ofs.close();

	return true;
}

template <typename data_t>
bool BDH<data_t>::saveTable(
	index_t  num,	//!< [in] number of data points.
	index_t* index, //!< [in] index of data point
	data_t** data,	//!< [in] data point set. 
	const string& path
	)
{
	ofstream ofs(path, ios::binary);
	if (ofs.fail())
	{
		ofs.close();
		cerr << "failed to open " << path << endl;
		return false;
	}

	//alloc workspace
	size_t* hashValue = new size_t[num];

	long long longNum = num;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (long long n = 0; n < longNum; ++n)
	{
		hashValue[n] = hashFunction(data[n]);
	}

	cout << "calculate hash values." << endl;
	//各データ点のハッシュ値を計算しコリジョンを取得
	vector<vector<index_t>> IDtable(hashSize);
	for (index_t n = 0; n < num; ++n)
	{
		IDtable[hashValue[n]].push_back(n);
	}
	delete[] hashValue;

	cout << "store data points." << endl;
	collision_t* collision = new collision_t[hashSize];
	long long longHashSize = hashSize;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (long long h = 0; h < longHashSize; ++h)
	{
		collision[h] = (collision_t)IDtable[h].size();
	}

	ofs << entrySize << "\t" << hashSize << "\t" << num;
	ofs.write((char*)collision, sizeof(collision_t)*hashSize);
	delete[] collision;

	vector<char> entry;
	char* entry_p;
	size_t progressStep = hashSize / 20;
	size_t progressPoint = progressStep;
	for (size_t h = 0; h < hashSize; ++h)
	{
		if (IDtable[h].empty() == false)
		{
			entry.resize(entrySize*IDtable[h].size());
			entry_p = entry.data();

			// write hash chain structure on RAM
			for (auto& obj : IDtable[h])
			{
				memcpy(entry_p, data[obj], pointSize);
				*reinterpret_cast<index_t*>(entry_p + pointSize) = index[obj];
				entry_p += entrySize;
			}

			// write on HD
			ofs.write(entry.data(), entry.size());
		}
		else if (h > progressPoint)
		{
			cout << h * 100 / hashSize << " % done." << endl;
			progressPoint += progressStep;
		}
	}
	ofs.close();

	return true;
}


/*4 SearchFunction*/
/*4(A) Estimated distance<=Rmax,K Nearest NeighborSearchFunction*/
template <typename data_t>
size_t BDH<data_t>::BicromaticReverseNearestNeighbor(
	int nQuery,
	data_t** query,
	vector<point_t<data_t>>* RKNNpoint,
	double searchParam,
	search_mode searchMode,
	int K,
	double epsilon
	)const
{

	//探索対象となるバケットのハッシュキーを生成 start//
	vector<hashKey_t> bucketList;
	unordered_map<size_t, vector<size_t>*> queryListFromBucket;//first: hash key, second: vector of query index 
	size_t NNC = 0;
	for (int q = 0; q < nQuery; ++q)
	{
		bucketList.clear();
		NNC += calcBucketList(query[q], searchParam, searchMode, bucketList);

		//生成されたクエリ→ハッシュキーのリストからハッシュキー→クエリのmapを作成
		addQueryListFromBucket(q, bucketList, queryListFromBucket);
	}

	//クエリリストから逆最近傍店を求める
	int count = extractRNNpoints(K, RKNNpoint, query, queryListFromBucket, epsilon);

	return NNC;
}

template <typename data_t>
int BDH<data_t>::extractRNNpoints(
	int K,
	vector<point_t<data_t>>* RKNNpoint,
	data_t** query,
	unordered_map<size_t, vector<size_t>*>& bucketToQuery,
	double epsilon
	)const
{

	address_t address;
	collision_t coll;
	bin_t bin;
	double dist, KNNdist;

	priority_queue<point_t<data_t>> NNqueries;

	point_t<data_t> dummy(index_t(-1), epsilon);
	point_t<data_t> pushPoint;
	point_t<data_t> dataPoint;

	//bucketToQuery内の各バケットにアクセス．
	for (auto& btq_itr : bucketToQuery)
	{
		//バケット内のデータを取り出す
		hashTable.getBin(btq_itr.first, bin);
		coll = bin.collision;
		address = bin.addressOfChainList;
		for (collision_t c = 0; c < coll; ++c, address += entrySize)
		{
			dataPoint.setMemberVariable(
				*reinterpret_cast<index_t*>(address + pointSize),
				reinterpret_cast<data_t*>(address),
				0.0
				);
			//K近傍点の初期化
			for (int i = 0; i < K; ++i)
			{
				NNqueries.push(dummy);
			}
			KNNdist = epsilon;

			for (auto& ql_itr : *btq_itr.second)
			{

				//距離計算
				dist = Distance(dim, dataPoint.addressOfpoint, query[ql_itr], KNNdist);
				if (dist < KNNdist)
				{
					NNqueries.pop();
					pushPoint.index = ql_itr;
					pushPoint.distance = dist;
					NNqueries.push(pushPoint);
					KNNdist = NNqueries.top().distance;
				}

			}

			while (NNqueries.empty() == false)
			{
				dataPoint.distance = NNqueries.top().distance;

				RKNNpoint[NNqueries.top().index].push_back(dataPoint);
				NNqueries.pop();
			}
		}

		delete btq_itr.second;
	}

	return 0;
}


template <typename data_t>
void BDH<data_t>::calcLayerParam(
	layer_t*& layer,
	data_t* const query
	) const {

	layer_t* layer_p = layer;
	layer_t* layer_p_end = layer + M;
	for (int m = 0; layer_p != layer_p_end; ++m, ++layer_p)
	{
		layer_p->k = subspace[m].subHashSize;
		subspace[m].setNodeParam(layer_p->node, query);
		layer_p->calc_gap();
	}
	sort(layer, layer + M);

	//m番目の部分空間からの残り距離の最大と最小を計算
	double distRestMin = 0;
	double distRestMax = 0;
	layer_p_end = layer - 1;
	for (layer_p = layer + M - 1; layer_p != layer_p_end; --layer_p)
	{
		layer_p->restMin = distRestMin += layer_p->node[0].distance;
		layer_p->restMax = distRestMax += layer_p->node[layer_p->k - 1].distance;
	}

}

template <typename data_t>
void BDH<data_t>::linearSearchInNNcandidates(
	data_t* query,
	point_t<data_t>* point,
	int K,
	double epsilon,
	vector<hashKey_t>& bucketList
	) const
{
	//生成されたハッシュキーを元にバケットを参照して最近傍点を探索 start//
	priority_queue<point_t<data_t>> NNpointQue;
	//最近傍点保持用のヒープ木を初期化
	index_t nullIndex = index_t(-1);
	point_t<data_t> dummy(nullIndex, nullptr, epsilon);
	for (int i = 0; i < K; ++i)
	{
		NNpointQue.push(dummy);
	}

	//見つけてきたハッシュキーを参照して最近傍点を探索する
	vector<hashKey_t>::iterator keyList_itr = bucketList.begin();
	vector<hashKey_t>::iterator keyList_itr_end = bucketList.end();
	for (; keyList_itr != keyList_itr_end; ++keyList_itr)
	{
		searchInBucket(query, (*keyList_itr).hashKey, NNpointQue);
	}
	//生成されたハッシュキーを元にバケットを参照して最近傍点を探索 end//

	//extract KNN point
	for (int i = K - 1; NNpointQue.empty() == false; --i)
	{
		point[i] = NNpointQue.top();
		NNpointQue.pop();
	}
}



template <typename data_t>
int BDH<data_t>::calcBucketList(
	data_t* query,
	double searchParam,
	search_mode searchMode,
	vector<hashKey_t>& bucketList
	)const
{

	//部分距離を計算し，優先度の高い順にソート
	layer_t* layer = new layer_t[M];
	for (int m = 0; m < M; ++m)
	{
		layer[m].node = new node_t[subspace[m].subHashSize + 1];
		layer[m].node[subspace[m].subHashSize].distance = DBL_MAX;
	}
	calcLayerParam(layer, query);

	double* lestSpaceVal = new double[lestspace.dim];
	lestspace.getPCAdata(query, lestSpaceVal);
	const double lestSpaceDist 
		= lestspace.getDistanceToCentroid(lestSpaceVal, 0);
	delete[] lestSpaceVal;

	unsigned NNC = 0;
	double preRadius = 0.0;
	nearBucket(
		query, layer, lestSpaceDist, searchParam, 
		searchMode, NNC, preRadius, bucketList
		);

	//探索が終わったのでデリート
	for (int m = 0; m < M; ++m)
	{
		delete[] layer[m].node;
	}
	delete[] layer;

	return NNC;
}

template <typename data_t>
unsigned BDH<data_t>::nearBucket(
	data_t* query,
	layer_t* layer,
	double lestSpaceDist,
	double searchParam,
	search_mode searchMode,
	unsigned& NNC,
	double& preRadius,
	vector<hashKey_t>& bucketList
	)const
{

	status_t status;
	status.dist = lestSpaceDist;

	switch (searchMode)
	{
	case Radius:
		//Radius以下の距離を探索
		NNC += NearBucket_C(preRadius, searchParam, layer, status, bucketList);
		break;

	case NumPoints:
	{
		unsigned C = static_cast<unsigned>(searchParam);
		bucketList.reserve(C);

		double Ubound = preRadius + delta;
		double tmpUbound = lestSpaceDist + layer[0].restMin;
		if (Ubound <= tmpUbound)
		{
			Ubound = tmpUbound + 1.0e-10;
		}
		for (double Lbound = preRadius
			; NNC < C
			; Ubound += delta)
		{
			NNC += NearBucket_C(Lbound, Ubound, layer, status, bucketList);
			Lbound = Ubound;
		}

	}
	break;

	case NumPoints2:
	{
		unsigned C = static_cast<unsigned>(searchParam);
		bucketList.reserve(C);

		//前回の探索で打ち切られたルートを再探索
		list<status_t> statusQue;
		statusQue.push_front(status);//ルートノード
		size_t loop, l;
		list<status_t>::iterator itr;

		double Ubound = preRadius + delta;
		double tmpUbound = lestSpaceDist + layer[0].restMin;
		if (Ubound <= tmpUbound)
		{
			Ubound = tmpUbound + 1.0e-10;
		}
		for (double Lbound = preRadius
			; NNC < C
			; Ubound += delta)
		{
			loop = statusQue.size();
			itr = statusQue.begin();
			for (l = 0; l < loop; ++l)
			{
				NNC += NearBucket_C_list(Ubound, layer, statusQue, &itr, bucketList);
			}
		}
	}
	break;

	default:
		break;
	}

	return NNC;
}

template <typename data_t>
size_t BDH<data_t>::BicromaticReverseNearestNeighbor_Incremental(
	param_for_incremental_search<data_t>& pfis,
	double searchParam,
	vector<point_t<data_t>>* RNNpoint
	)const
{

	//探索対象となるバケットのハッシュキーを生成 start//
	vector<hashKey_t> bucketList;
	unordered_map<size_t, vector<size_t>*> queryListFromBucket;//first: hash key, second: vector of query index 
	size_t NNC = 0;
	for (unsigned q = 0; q < pfis.nQuery; ++q)
	{
		bucketList.clear();
		NNC += nearBucket(
			pfis.query[q],
			pfis.layer[q],
			pfis.lestSpaceDist[q],
			searchParam,
			pfis.searchMode,
			pfis.NNC[q],
			pfis.preRadius[q],
			bucketList);

		//生成されたクエリ→ハッシュキーのリストからハッシュキー→クエリのmapを作成
		addQueryListFromBucket(q, bucketList, queryListFromBucket);
	}

	extractNNqueries(pfis, queryListFromBucket);
	for (unsigned q = 0; q < pfis.nQuery; ++q)
	{
		RNNpoint[q].clear();
	}
	for (size_t n = 0; n < pfis.nData; ++n)
	{
		if (pfis.NNquery[n].index != UINT_MAX)
		{
			RNNpoint[pfis.NNquery[n].index].push_back(point_t<data_t>(n, NULL, pfis.NNquery[n].dist));
		}
	}

	return NNC;
}

template <typename data_t>
void BDH<data_t>::InitializeForBRNN_Incremental(param_for_incremental_search<data_t>& pfis)
{
	pfis.M = M;
	double* lestSpaceVal = new double[lestspace.dim];
	for (unsigned q = 0; q < pfis.nQuery; ++q)
	{
		//部分距離を計算し，優先度の高い順にソート
		layer_t* layer = new layer_t[M];
		for (int m = 0; m < M; ++m)
		{
			layer[m].node = new node_t[subspace[m].subHashSize + 1];
			layer[m].node[subspace[m].subHashSize].distance = DBL_MAX;
		}
		calcLayerParam(layer, pfis.query[q]);

		pfis.layer[q] = layer;

		lestspace.getPCAdata(pfis.query[q], lestSpaceVal);
		pfis.lestSpaceDist[q] = lestspace.getDistanceToCentroid(lestSpaceVal, 0);
	}
	delete[] lestSpaceVal;
}

template <typename data_t>
int BDH<data_t>::NearestNeighbor(
	data_t* query,
	point_t<data_t>* point,
	double searchParam,
	search_mode searchMode,
	int K,
	double epsilon
	)const
{

	//layer（部分空間のバケット距離などのパラメータ）のメモリ確保
	layer_t* layer = new layer_t[M];
	for (int m = 0; m < M; ++m)
	{
		layer[m].node = new node_t[subspace[m].subHashSize + 1];
		layer[m].node[subspace[m].subHashSize].distance = DBL_MAX;
	}
	//layerのパラメータを計算する
	calcLayerParam(layer, query);

	//バケット距離の計算に使われない基底空間(lestspace)における推定距離を計算
	double* lestSpaceVal = new double[lestspace.dim];
	lestspace.getPCAdata(query, lestSpaceVal);
	const double lestSpaceDist = lestspace.getDistanceToCentroid(lestSpaceVal, 0);
	delete[] lestSpaceVal;

	unsigned NNC = 0;
	status_t status;
	status.dist = lestSpaceDist;

	//最近傍点を管理するキューにダミーのデータを登録しておく
	//epsilonより距離の近い点が見つかればダミーは排除される．
	point_t<data_t> dummy(UINT_MAX, nullptr, epsilon);
	priority_queue<point_t<data_t>> NNpointQue;
	for (int i = 0; i < K; ++i)
	{
		NNpointQue.push(dummy);
	}

	//searchModeに探索を行う
	switch (searchMode)
	{
	case Radius:
		//Radius以下の距離を探索
		NNC = NearBucket_R(searchParam, layer, status, query, NNpointQue);

		break;

	case NumPoints:
	{
		unsigned C = static_cast<unsigned>(searchParam);
		for (double Lbound = 0, Ubound = lestSpaceDist + layer[0].restMin + 1.0e-10
			; NNC < C
			; Ubound += delta)
		{
			NNC += NearBucket_C(Lbound, Ubound, layer, status, query, NNpointQue);
			Lbound = Ubound;
		}

		break;
	}

	}

	//探索が終わったのでデリート
	for (int m = 0; m < M; ++m)
	{
		delete[] layer[m].node;
	}
	delete[] layer;

	//発見した最近傍点を答え用のコンテナ(point)にセットする
	for (int i = K - 1; NNpointQue.empty() == false; --i)
	{
		point[i] = NNpointQue.top();
		NNpointQue.pop();
	}

	return NNC;//距離計算した点数
}

template <typename data_t>
void BDH<data_t>::addQueryListFromBucket(
	int queryIndex,
	const vector<hashKey_t>& bucketList,
	unordered_map<size_t, vector<size_t>*>& queryListFromBucket
	)const
{
	pair< unordered_map<size_t, vector<size_t>*>::iterator, bool> map_insert_ret;
	pair<size_t, vector<size_t>*> map_pair;

	for (auto itr : bucketList)
	{
		map_pair.first = itr.hashKey;
		map_pair.second = new vector<size_t>(1, queryIndex);

		map_insert_ret = queryListFromBucket.insert(map_pair);

		if (map_insert_ret.second == false)
		{//既にhashKeyが登録されていた場合は

			//後に挿入
			map_insert_ret.first->second->push_back(queryIndex);
		}
	}

}

template <typename data_t>
void BDH<data_t>::getSubHashkey(size_t hashKey, int* subHashkey)const
{
	for (int m = 0; m < M; ++m)
	{
		subHashkey[m] = hashKey % subspace[m].subHashSize;
		hashKey /= subspace[m].subHashSize;
	}
}

template <typename data_t>
int BDH<data_t>::extractNNqueries(
	param_for_incremental_search<data_t>& pfis,
	unordered_map<size_t, vector<size_t>*>& bucketToQuery
	)const
{

	address_t address;
	collision_t coll;
	bin_t bin;
	double dist;
	point_t<data_t> pushPoint;
	point_t<data_t> dataPoint;

	//bucketToQuery内の各バケットにアクセス．
	for (auto& btq_itr : bucketToQuery)
	{
		//バケット内のデータを取り出す
		hashTable.getBin(btq_itr.first, bin);
		coll = bin.collision;
		address = bin.addressOfChainList;
		for (collision_t c = 0; c < coll; ++c, address += entrySize)
		{
			dataPoint.setMemberVariable(
				*reinterpret_cast<index_t*>(address + pointSize),
				reinterpret_cast<data_t*>(address),
				0.0
				);

			pushPoint.distance = pfis.NNquery[dataPoint.index].dist;
			pushPoint.index = pfis.NNquery[dataPoint.index].index;

			for (auto& ql_itr : *btq_itr.second)
			{

				//距離計算
				dist = Distance(dim, dataPoint.addressOfpoint, pfis.query[ql_itr], pushPoint.distance);
				if (dist < pushPoint.distance)
				{
					pushPoint.index = ql_itr;
					pushPoint.distance = dist;
				}
			}

			//extract KNN queries
			pfis.NNquery[dataPoint.index].dist = static_cast<float>(pushPoint.distance);
			pfis.NNquery[dataPoint.index].index = static_cast<unsigned>(pushPoint.index);
		}

		delete btq_itr.second;
	}

	return 0;
}



template class BDH < unsigned char >;
//template class BDH < unsigned short >;
template class BDH < float >;
template class BDH < double >;
