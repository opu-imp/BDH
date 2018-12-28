/**
* @brief test sample code
* @date  2015/05/04 new
* @author Tomokazu Sato
*/

#include <measure.h>
#include <FileHandle.h>
#include <sstream>
#include <BDH.h>
#include <NNpointIO.h>

#ifdef __unix__
#include <octavePCA.h>
#else
#include <opencvPCA.h>
#endif

typedef unsigned char data_t;
typedef unsigned index_t;

struct argumentsSet
{
	//ハッシュテーブルに登録される点
	string dataPointPath;
	string PCApath;

	//粗量子化パラメータ
	string parameterPath;
	int B;//インデクシングのビット数
	int M; //ハッシュテーブル構築に利用する部分空間の数
	int P;//一つの部分空間内の次元数
	int K; //K近傍探索

	//パラメータチューニングの精度パラメータ
	double bit_step;
	double sampling_rate;
	string tablePath;

	//クエリ点
	string queryPointPath;

	//探索パラメータ
	search_mode searchMode;//探索アルゴリズム
	double paramRate;//ハッシュに登録した点数×paramRate点と距離計さんを行う
	double Epsilon;//Epsilon以内の点のみを対象とする（K近傍でもEpsilonより遠い点は解とならない）

	//探索結果
	string resultPath;
};


void setOptions(int argc, char** argv, argumentsSet& args);

void Indexing(
	int dim, unsigned num, data_t** data,
	BDH<data_t>& bdh, const argumentsSet& args);

void outputSearchResult(const argumentsSet& args, index_t nQuery, int* NNC, point_t<data_t>** KNNpoint, index_t nDataPoints, double queryTime);

int main(int argc, char** argv)
{
	argumentsSet args;
	args.dataPointPath = "./../../testSamples/sift10K.ucdat";
	args.parameterPath = "parameter.bdh";
	args.B = 13;
	args.M = 3;
	args.P = 10;
	args.K = 1;
	args.bit_step = 0.1;
	args.sampling_rate = 1.0;
	args.tablePath = "BDHtable.tbl";
	args.queryPointPath = "./../../testSamples/sift1K.ucquery";
	args.searchMode = NumPoints;
	args.paramRate = 0.001;
	args.Epsilon = DBL_MAX;
	args.resultPath = "result.txt";

	setOptions(argc, argv, args);

	//read data point set
	cout << "read data point set." << endl;
	int dim;
	unsigned num;
	data_t** data;
	readBinary(args.dataPointPath, dim, num, data);

	// build BDH hash tabe
	BDH<data_t> bdh;
	Indexing(dim, num, data, bdh, args);

	//delete data point
	for (index_t n = 0; n < num; ++n)
	{
		delete[] data[n];
	}
	delete[] data;


	// 探索パラメータ
	double searchParam;
	switch (args.searchMode)
	{
	case Radius:
		searchParam = bdh.get_variance()*args.paramRate;
		break;
	default:
		searchParam = static_cast<index_t>(bdh.get_nDdataPoints()*args.paramRate);
		break;
	}

	//read query point set
	cout << "read query point set." << endl;
	unsigned nQuery;
	data_t** query;
	fileErrorHandle(
		readBinary(args.queryPointPath, dim, nQuery, query),
		args.queryPointPath, __FILE__, __LINE__
		);

	cout << "start to test search." << endl;
	int* NNC = new int[nQuery];
	point_t<data_t>** KNNpoint = new point_t<data_t>*[nQuery];
	for (unsigned q = 0; q < nQuery; ++q)
	{
		KNNpoint[q] = new point_t<data_t>[args.K];
	}

	// test perform for each Cnum by pointsnum saerch
	double startTime = GetCPUTime();
	for (unsigned q = 0; q < nQuery; ++q)
	{
		NNC[q] = bdh.NearestNeighbor(query[q], KNNpoint[q], searchParam, args.searchMode, args.K, args.Epsilon);
	}
	double endTime = GetCPUTime();

	outputSearchResult(args, nQuery, NNC, KNNpoint, bdh.get_nDdataPoints(), (endTime - startTime) / nQuery);

	return 0;
}

void setOptions(int argc, char** argv, argumentsSet& args)
{
	for (int i = 1; i < argc; ++i)
	{
		if (argv[i][0] == '-')
		{
			switch (argv[i][1])
			{
			case 'd':
				args.dataPointPath = argv[++i];
				break;
			case 'a':
				args.PCApath = argv[++i];
				break;
			case 'p':
				args.parameterPath = argv[++i];
				break;
			case 't':
				args.tablePath = argv[++i];
				break;
			case 'q':
				args.queryPointPath = argv[++i];
				break;
			case 'r':
				args.resultPath = argv[++i];
				break;
			case 'B':
				args.B = atoi(argv[++i]);
				break;
			case 'M':
				args.M = atoi(argv[++i]);
				break;
			case 'P':
				args.P = atoi(argv[++i]);
				break;
			case 'K':
				args.K = atoi(argv[++i]);
				break;
			case 'S':
				args.bit_step = atof(argv[++i]);
				break;
			case 'L':
				args.sampling_rate = atof(argv[++i]);
				break;
			case 'R':
				args.paramRate = atof(argv[++i]);
				break;
			case 'E':
				args.Epsilon = atof(argv[++i]);
				break;
			case 'A':
				args.searchMode = static_cast<search_mode>(atoi(argv[++i]));
				break;
			}
		}
	}
}

void Indexing(
	int dim, unsigned num, data_t** data,
	BDH<data_t>& bdh, const argumentsSet& args)
{

	//Principal Component Analysis
	cout << "calculate PCA ." << endl;
	PrincipalComponentAnalysis pca;

	if (pca.loadPCA(args.PCApath) == false)
	{
		pca.executePCA(dim, num, data);
		pca.savePCA(args.PCApath);
	}

	// copy PCA direction to base_t for BDH
	const PC_t* pcDir = pca.getPCdir();
	base_t* base = new base_t[dim];
	for (int d = 0; d < dim; ++d)
	{
		base[d].mean = pcDir[d].mean;
		base[d].variance = pcDir[d].variance;
		base[d].direction = new double[dim];
		memcpy(base[d].direction, pcDir[d].direction, sizeof(double)*dim);
	}

	cout << "training Start ." << endl;
	// train parameters
	if (bdh.loadParameters(args.parameterPath) == false)
	{
		bdh.parameterTuning_ICCV2013(dim, num, data, base, args.P, args.B, args.bit_step, args.sampling_rate);
		bdh.saveParameters(args.parameterPath);
	}

	//delete base
	for (int d = 0; d < dim; ++d)
	{
		delete[] base[d].direction;
	}
	delete[] base;

	// entory data points into hash table
	if (bdh.loadTable(args.tablePath) == false)
	{
		bdh.storePoints(num, data);
		bdh.saveTable(args.tablePath);
	}

}

void outputSearchResult(const argumentsSet& args, index_t nQuery, int* NNC, point_t<data_t>** KNNpoint, index_t nDataPoints, double queryTime)
{
	ofstream ofst(args.resultPath);
	openErrorHandle(ofst, args.resultPath, __FILE__, __LINE__);
	ofst << "average query time : " << queryTime << " ms" << endl;
	cout << "average query time : " << queryTime << " ms" << endl;
	ofst << "memory usage : " << static_cast<double>(Memory()) / (1024 * 1024) << " MB" << endl;
	cout << "memory usage : " << static_cast<double>(Memory()) / (1024 * 1024) << " MB" << endl;

	for (index_t qv = 0; qv < nQuery; ++qv)
	{
		ofst << "query index : " << qv << "\tNNC : " << NNC[qv] << endl;
		for (int i = 0; i < args.K && (KNNpoint[qv][i].index < nDataPoints); ++i)
		{
			ofst << KNNpoint[qv][i].index << "\t" << KNNpoint[qv][i].distance << "\t";
		}
		ofst << endl;
	}
	ofst.close();
}