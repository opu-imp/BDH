
#include <measure.h>
#include <FileHandle.h>
#include <sstream>
#include <BDH.h>

#ifdef __unix__
#include <octavePCA.h>
#else
#include <opencvPCA.h>
#endif

//#define FIXED_PARAM
#define ICCV2013

typedef unsigned char data_t;
typedef unsigned index_t;

int main(int argc,char** argv)
{

#ifdef FIXED_PARAM
	const string samplePath   = "../../testSamples/sift10K.ucdat";
	const string pcaPath      = "../../testSamples/sift10K.pca";
	const string bdhParamPath = "../../testSamples/sift10K_bit14_M5_P10_BS0.1_SR1_ICCV2013.bdh";
	const int bit = 14;
	const int M = 5;
	const int P = 10;
	const double bit_step = 0.1;
	const double sampling_rate = 1;

#else
	//set training parameters
	int i = 1;
	const int bit = atoi(argv[i++]);
	const int M = atoi(argv[i++]);
	const int P = atoi(argv[i++]);
	const double bit_step = atof(argv[i++]);
	const double sampling_rate = atof(argv[i++]);
	const string samplePath = argv[i++];
	const string pcaPath = argv[i++];
	const string bdhParamPath = argv[i++];
#endif

	if (isFileExist(bdhParamPath))
	{
		cerr << bdhParamPath << "already exists.";
		return __LINE__;
	}

	//read Indexing samples
	cout << "read training sample" << endl;
	int dim;
	index_t num;
	data_t** data = nullptr;
	readBinary(samplePath, dim, num, data);

	//principal component analysis
	base_t* base;
	{
		cout << "principal component analysis" << endl;
		PrincipalComponentAnalysis CVpca;
		if (CVpca.loadPCA(pcaPath) == false)
		{
			CVpca.executePCA(dim, num, data);
			fileErrorHandle(CVpca.savePCA(pcaPath), pcaPath, __FILE__, __LINE__);
		}

		//主成分分析で求めた情報をbase_tにセット
		base = new base_t[dim];
		const PC_t* pcDir = CVpca.getPCdir();
		for (int d = 0; d < dim; ++d)
		{
			base[d].mean = pcDir[d].mean;
			base[d].variance = pcDir[d].variance;
			base[d].direction = new double[dim];
			memcpy(
				base[d].direction,
				pcDir[d].direction,
				sizeof(double)*dim
				);
		}
	}

	cout << "train BDH paramters" << endl;
	BDH<data_t> bdh;
	//パラメータ学習
#ifdef ICCV2013
	bdh.parameterTuning_ICCV2013(dim, num, data, base, P, bit, bit_step, sampling_rate);
#else
	bdh.parameterTuning(dim, num, data, base, M, P, bit, bit_step, sampling_rate);
#endif

	fileErrorHandle(
		bdh.saveParameters(bdhParamPath),
		bdhParamPath,
		__FILE__,
		__LINE__
		);

	//基底のメモリ解放
	for (int d = 0; d < dim; ++d)
	{
		delete[] base[d].direction;
	}
	delete[] base;

	//サンプルメモリ解放
	deleteDataArray(num, data);

	return 0;
}