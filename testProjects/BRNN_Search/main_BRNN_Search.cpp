/**
* @file main_BDH.cpp
* @brief main関数を定義しています。<br>
* @author Tomokazu Sato
*/

#include <measure.h>
#include <sstream>

#include <FileHandle.h>
#include <BDH.h>
#include <NNpointIO.h>

typedef unsigned char data_t;
typedef unsigned index_t;

//for training
const string datDir = "./../../testSamples";
const string trainingDir = "./../../testSamples";
const int bit = 14;
const int M = 5;
const int P = 10;
const double bit_step = 0.1;

//for brute force nearest neighbor search
const int K = 1;
const string datName = "sift10K.ucdat";
const string queryName = "sift1K.ucquery";
const string answerPath = datDir + "/sift10K.ucdat_sift1K.ucquery_K1.brnn";

int main(void)
{

	const string samplePath = datDir + "/" + datName;

	stringstream trainingFileNameStrm;
	trainingFileNameStrm << datName << "_bit" << bit << "_M" << M << "_P" << P << "_BS" << int(bit_step * 10);
	const string trainingName = trainingFileNameStrm.str();
	const string bdhParamPath = trainingDir + "/" + trainingName + ".bdh";
	const string bdhTablePath = trainingDir + "/" + trainingName + ".tbl";

	const string queryPath = datDir + "/" + queryName;
	stringstream ssResult;
	ssResult << datName << "_bit" << bit << "_M" << M << "_P" << P << "_BS" << int(bit_step * 10) << "_" << queryName << "_K" << K << ".result";
	const string resultPath = ssResult.str();

	//read Indexing samples
	int dim;
	unsigned num;
	data_t** data;
	readBinary(samplePath, dim, num, data);


	BDH<data_t> bdh;
	if (bdh.loadParameters(bdhParamPath) == false)
	{
		return __LINE__;
	}

	//ハッシュテーブルにデータを登録する
	if (bdh.loadTable(bdhTablePath) == false)
	{
		bdh.storePoint(num, data);

		fileErrorHandle(
			bdh.saveTable(bdhTablePath),
			bdhTablePath,
			__FILE__,
			__LINE__
			);
	}

	unsigned nQuery;
	data_t** query;
	readBinary(queryPath, dim, nQuery, query);

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

	double recall = 0;
	double startTime, endTime;
	vector<point_t<data_t>>* RNNpoint = new vector<point_t<data_t>>[nQuery];

	ofstream result(resultPath);
	fileErrorHandle(result.is_open(), resultPath, __FILE__, __LINE__);

	vector<point_t<data_t>>::iterator itr;
	vector<point_t<data_t>>::iterator itrAns;
	for (unsigned R = 10000; recall < 0.99; R += 10000)
	{
		for (index_t qv = 0; qv < nQuery; qv++)
		{
			RNNpoint[qv].clear();
		}
		startTime = GetCPUTime();
		bdh.BicromaticReverseNearestNeighbor(nQuery, query, RNNpoint, R, Radius, K);
		endTime = GetCPUTime();

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
		cout << M << "\t" << bit << "\t" << R << "\t" << precision << "\t" << recall << "\t" << queryTime << endl;
		result << M << "\t" << bit << "\t" << R << "\t" << precision << "\t" << recall << "\t" << queryTime << endl;
	}
	result.close();
}