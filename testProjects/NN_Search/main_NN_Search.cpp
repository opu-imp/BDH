/**
* @file main_BDH.cpp
* @brief * @author Tomokazu Sato
*/

#include <sstream>

#include <measure.h>
#include <FileHandle.h>
#include <BDH.h>
#include <NNpointIO.h>

typedef unsigned char data_t;
typedef unsigned index_t;

#define FIXED_PARAM

int main(int argc, char** argv)
{

#ifdef FIXED_PARAM
	const string bdhParamPath = "../../testSamples/sift10K_bit14_M5_P10_BS0.1_SR1_ICCV2013.bdh";
	const string bdhTablePath = "../../testSamples/sift10K_bit14_M5_P10_BS0.1_SR1_ICCV2013_sift10K.tbl";
	const string dataPath = "../../testSamples/sift10K.ucdat";
	const string answerPath = "../../testSamples/sift10K.ucdat_sift1K.ucquery_K1.nn";
	const string queryPath = "../../testSamples/sift1K.ucquery";
	const string resultPath = "ICCV2013_sift10K_bit14_M5_P10_BS0.1_SR1_sift10K_sift1K.result";
	const int K = 1;
	const int CYCLE = 1;

#else
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

	BDH<data_t> bdh;
	cout << "load parameters." << endl;
	if (bdh.loadParameters(bdhParamPath) == false)
	{
		fileErrorHandle(false, bdhParamPath, __FILE__, __LINE__);
		return __LINE__;
	}

	//ハッシュテーブルにデータを登録する
	cout << "load table." << endl;
	if (bdh.loadTable(bdhTablePath) == false)
	{
		//read Indexing samples	
		int dim;
		index_t num;
		data_t** data;
		if (readBinary(dataPath, dim, num, data) == false)
		{
			fileErrorHandle(false, dataPath, __FILE__, __LINE__);
			return 1;
		}

		bdh.storePoints(num, data);

		deleteDataArray(num, data);

		fileErrorHandle(
			bdh.saveTable(bdhTablePath),
			bdhTablePath,
			__FILE__,
			__LINE__
			);
	}

	cout << "load query." << endl;
	int dim;
	unsigned nQuery;
	data_t** query;
	fileErrorHandle(
		readBinary(queryPath, dim, nQuery, query),
		queryPath, __FILE__, __LINE__);

	cout << "load answer." << endl;
	point_t<data_t>** answerPoint;
	size_t memoryUsage;
	double queryTime;
	int answerK;

	if (LoadNearestNeighbor(
		answerPath,
		nQuery, answerK, answerPoint,
		queryTime, memoryUsage) == false)
	{
		fileErrorHandle(false, answerPath, __FILE__, __LINE__);
	}
	double meanNNC;
	double recall = 0;
	double startTime, endTime;
	int* NNC = new int[nQuery];
	point_t<data_t>** KNNpoint = new point_t<data_t>*[nQuery];
	for (unsigned q = 0; q < nQuery; ++q)
	{
		KNNpoint[q] = new point_t<data_t>[K];
	}

	ofstream result(resultPath);
	openErrorHandle(result, resultPath, __FILE__, __LINE__);

	cout << "Cnum\t recall\t queryTime\t meanNNC" << endl;
	for (double a = 2; recall < 0.99; a += 0.5)
	{
		index_t Cnum = static_cast<int>(pow(2, a));
		startTime = GetCPUTime();

		for (int cy = 0; cy < CYCLE; ++cy)
		{
			for (unsigned q = 0; q < nQuery; ++q)
			{
				NNC[q] = bdh.NearestNeighbor(query[q], KNNpoint[q], Cnum, NumPoints,K);
			}
		}
		endTime = GetCPUTime();

		meanNNC = 0;
		recall = 0;

		for (index_t qv = 0; qv < nQuery; qv++)
		{
			int Correct = 0;

			for (int k = 0; k < K; k++)
			{
				for (int k2 = 0; k2 <= k; k2++)
				{
					//cout << KNNpoint[qv][k2].index << "\t" << answerPoint[qv][k].index << endl;
					if (KNNpoint[qv][k2].index == answerPoint[qv][k].index)
					{
						Correct++;
						break;
					}
				}
			}

			recall += static_cast<double>(Correct) / K;
			meanNNC += NNC[qv];
		}

		recall /= nQuery;
		meanNNC /= nQuery;
		result << Cnum << "\t" << recall << "\t" << (endTime - startTime) / (nQuery*CYCLE) << "\t" << meanNNC << endl;
		cout << Cnum << "\t" << recall << "\t" << (endTime - startTime) / (nQuery*CYCLE) << "\t" << meanNNC << endl;
	}
	result.close();
}