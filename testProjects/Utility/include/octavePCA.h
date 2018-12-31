/**
* @file    opencvPCA.h
* @author  T.Sato
* @date    2015.05.06
* @version 1.0
*/

#ifndef __OCTAVE_PCA__
#define __OCTAVE_PCA__

#include <superPCA.h>

#include <octave/config.h>
#include <octave/Matrix.h>
#include <stdio.h>
#include <string>
using namespace std;

/**
* @brief Pricipal Component Analysis
*/
class PrincipalComponentAnalysis : public superPCA
{

public:
	////コンストラクタ
	PrincipalComponentAnalysis()
		: superPCA()
	{}

	/**
	* @brief PCA (covariance)
	*/
	template <typename data_t>
	void executePCA(int dim, unsigned num, data_t** data);

	/**
	* @brief PCA (coefficient)
	*/
	template <typename data_t>
	bool executePCAcorelationCoefficient(
		int dim, int num, data_t** data);
};

template <typename data_t>
void PrincipalComponentAnalysis::executePCA(
	int dim, 
	unsigned num, 
	data_t** data)
{
	//次元数変わったらメモリ確保しなおし
	resetDimension(dim);

	//各基底の平均と分散共分散行列を得る
	double* Mean = new double[dim];
	double** covariance = new double*[dim];
	for (int d = 0; d < dim; ++d)
	{
		covariance[d] = new double[dim];
	}
	calculateCovarianceMatrix(dim, num, data, Mean, covariance);

	Matrix Cov(dim, dim);
	for (int d = 0; d < dim; ++d)
	{
		Cov(d, d) = covariance[d][d];
		for (int d2 = d + 1; d2 < dim; ++d2)
		{
			Cov(d, d2) = Cov(d2, d) = covariance[d][d2];
		}
	}

	/*Compute EigenVector*/
	EIG eig(Cov);

	//Copy EigenVector to nomal pointer Matrix
	for (int d = 0; d<dim; d++)
	{
		for (int d2 = 0; d2<dim; d2++)
		{
			pcDir[d].direction[d2] = real(eig.eigenvectors()(d2, dim - 1 - d));
		}
		pcDir[d].variance = real(eig.eigenvalues()(dim - 1 - d));
		pcDir[d].mean = superPCA::innerProduct(dim, pcDir[d].direction, Mean);
	}

	for (int d = dim - ZeroCount; d<dim; d++)
	{
		pcDir[d].variance = 0;
	}

	//sort int Descending order of variance
	sort(pcDir, pcDir + dim);

	delete[] Mean;
	for (int d = 0; d < dim; ++d)
	{
		delete[] covariance[d];
	}
	delete[] covariance;
}

#endif