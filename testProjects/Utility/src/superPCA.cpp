#include <superPCA.h>

#include <fstream>
using namespace std;

/**
* @brief save result of PCA
*/
bool superPCA::savePCA(const string& path)
{
	ofstream ofs(path);
	if (ofs.is_open() == false)
	{
		return false;
	}

	ofs << dim << endl;
	for (int d = 0; d < dim; ++d)
	{
		ofs << pcDir[d].mean << "\t" << pcDir[d].variance << endl;
		for (int d2 = 0; d2 < dim; ++d2)
		{
			ofs << pcDir[d].direction[d2] << "\t";
		}
		ofs << endl;
	}

	return true;
}

/**
* @brief load result of PCA
*/
bool superPCA::loadPCA(
	const string& path
	)
{
	ifstream ifs(path);
	if (ifs.is_open() == false)
	{
		return false;
	}

	ifs >> dim;
	pcDir = new PC_t[dim];
	for (int d = 0; d < dim; ++d)
	{
		ifs >> pcDir[d].mean >> pcDir[d].variance;
		pcDir[d].direction = new double[dim];
		for (int d2 = 0; d2 < dim; ++d2)
		{
			ifs >> pcDir[d].direction[d2];
		}
	}

	return true;
}

//ŽŸŒ³”•Ï‚í‚Á‚½‚çƒƒ‚ƒŠŠm•Û‚µ‚È‚¨‚µ
void superPCA::resetDimension(int dim)
{
	if (this->dim != dim)
	{
		//ƒƒ‚ƒŠ‰ð•ú
		this->~superPCA();

		//ƒƒ‚ƒŠ‰Šú‰»
		this->dim = dim;
		pcDir = new PC_t[dim];
		
		for (int d = 0; d < dim; ++d)
		{
			pcDir[d].direction = new double[dim];
		}
	}
}
