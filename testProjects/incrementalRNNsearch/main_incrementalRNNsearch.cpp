
#include <measure.h>
#include <sstream>

#include <FileHandle.h>
#include <BDH.h>
#include <NNpointIO.h>

typedef unsigned char data_t;
typedef unsigned index_t;

#define FIXED_PARAM

int main(int argc, char** argv)
{
	
	//引数の設定
#ifdef FIXED_PARAM
	//関数べた書き
	const string bdhParamPath = "../../testSamples/sift10K_bit14_M5_P10_BS0.1_SR1_ICCV2013.bdh";
	const string bdhTablePath = "../../testSamples/sift10K_bit14_M5_P10_BS0.1_SR1_ICCV2013_sift10K.tbl";
	const string dataPath = "../../testSamples/sift10K.ucdat";
	const string answerPath = "../../testSamples/sift10K.ucdat_sift1K.ucquery_K1.brnn";
	const string queryPath = "../../testSamples/sift1K.ucquery";
	const string resultPath = "ICCV2013_sift10K_bit14_M5_P10_BS0.1_SR1_sift10K_sift1K.result";
	const int K = 1;
	const int CYCLE = 1;
	search_mode searchMode = NumPoints;
#else
	//mainの引数から取得
	int i = 1;
	const string bdhParamPath = argv[i++];
	const string bdhTablePath = argv[i++];
	const string dataPath = argv[i++];
	const string answerPath = argv[i++];
	const string queryPath = argv[i++];
	const string resultPath = argv[i++];
	const int K = atoi(argv[i++]);
	const int CYCLE = atoi(argv[i++]);
#endif

	//BDHをインスタンスする
	BDH<data_t> bdh;

	//BDHのパラメータを取得する
	if (bdh.loadParameters(bdhParamPath) == false)
	{
		return __LINE__;
	}

	//ハッシュテーブルにデータを登録する
	if (bdh.loadTable(bdhTablePath) == false)
	{
		return __LINE__;
	}

	//クエリを読み込む
	int dim;
	unsigned nQuery;
	data_t** query;
	readBinary(queryPath, dim, nQuery, query);

	//正解データを読み込む
	vector<point_t<data_t>>* RNNpointAns = new vector<point_t<data_t>>[nQuery];
	size_t memoryUsage;
	double queryTime;
	int answerK;
	fileErrorHandle(
		LoadReverseNearestNeighbor(
		answerPath,
		nQuery,
		answerK,
		RNNpointAns,
		queryTime,
		memoryUsage
		),
		answerPath,
		__FILE__,
		__LINE__
		);


	//探索の途中状態を保持するための構造体
	param_for_incremental_search<data_t> pfis(bdh.get_nDdataPoints(), nQuery, query, searchMode);
	
	//BDHによる検索のための初期化。探索を行う前に実行することが必須。
	bdh.InitializeForBRNN_Incremental(pfis);

	double recall = 0;
	double startTime, endTime;
	vector<point_t<data_t>>* RNNpoint = new vector<point_t<data_t>>[nQuery];
	vector<point_t<data_t>>::iterator itr;
	vector<point_t<data_t>>::iterator itrAns;
	ofstream result(resultPath);
	fileErrorHandle(result.is_open(), resultPath, __FILE__, __LINE__);
	//探索点数を２倍ずつ増やしていく
	for (unsigned searchParam = 10; recall < 0.99; searchParam *= 2)
	{
		//探索を実行する
		startTime = GetCPUTime();
		bdh.BicromaticReverseNearestNeighbor_Incremental(pfis, searchParam,RNNpoint);
		endTime = GetCPUTime();

		//探索結果のPrecisionとRecallを計算する
		size_t resultSum = 0;
		size_t ansSum = 0;
		int Correct = 0;
		for (index_t qv = 0; qv < nQuery; qv++)
		{
			resultSum += RNNpoint[qv].size();
			ansSum += RNNpointAns[qv].size();
			for (itr = RNNpoint[qv].begin();
				itr != RNNpoint[qv].end();
				++itr)
			{
				for (itrAns = RNNpointAns[qv].begin();
					itrAns != RNNpointAns[qv].end();
					++itrAns)
				{
					if (itr->index == itrAns->index)
					{
						++Correct;
					}
				}
			}
		}
		recall = static_cast<double>(Correct) / ansSum;
		double precision = static_cast<double>(Correct) / resultSum;
		double queryTime = endTime - startTime;

		//結果の出力
		cout << searchParam << "\t" << precision << "\t" << recall << "\t" << queryTime << "\t" << Correct << endl;
		result << searchParam << "\t" << precision << "\t" << recall << "\t" << queryTime << "\t" << Correct << endl;
	}
	result.close();


	return 0;
}