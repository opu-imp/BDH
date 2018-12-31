/**
* @file    opencvPCA.h
* @author  T.Sato
* @date    2015.05.06
* @version 1.0
*/

#ifndef __opencvPCA__
#define __opencvPCA__
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <superPCA.h>
#include <string>
#include <iostream>
#include <fstream>
#include <measure.h>

using namespace std;
using namespace cv;

#define SHOW_PROGRESS

/**
* @brief Pricipal Component Analysis
*/
class PrincipalComponentAnalysis : public superPCA
{

public:
	//コンストラクタ
	PrincipalComponentAnalysis()
		: superPCA()
	{}

	/**
	* @brief PCA (covariance)
	*/
	template <typename data_t>
	void executePCA(int dim, size_t num, data_t** data);

	/**
	* @brief executePCA (coefficient)
	*/
	template <typename data_t>
	bool executePCAcorelationCoefficient(
		int dim, int num, data_t** data);
};

template <typename data_t>
void PrincipalComponentAnalysis::executePCA(
	int dim, 
	size_t num, 
	data_t** data)
{
	//次元数変わったらメモリ確保しなおし
	resetDimension(dim);

	//各基底の平均と分散共分散行列を得る
#ifdef SHOW_PROGRESS
	cout << "calclate covariance matrix" << endl;
	double timeStart = GetCPUTime();
#endif
	double* Mean = new double[dim];
	double** covariance = new double*[dim];
	for (int d = 0; d < dim; ++d)
	{
		covariance[d] = new double[dim];
	}
	calculateCovarianceMatrix(dim, num, data, Mean, covariance);

	Mat Cov(dim, dim, CV_64FC1);
	for (int d = 0; d < dim; ++d)
	{
		Cov.at<double>(d, d) = covariance[d][d];
		for (int d2 = d + 1; d2 < dim; ++d2)
		{
			Cov.at<double>(d, d2) = Cov.at<double>(d2, d) = covariance[d][d2];
		}
	}
#ifdef SHOW_PROGRESS
	double timeEnd = GetCPUTime();
	cout << static_cast<int>(timeEnd-timeStart)/1000 << " sec\t" << endl;
	cout << "Eigenvalue decomposition" << endl;
	timeStart = GetCPUTime();
#endif

	//calc eigen vector
	Mat EigVal, EigVec;
	eigen(Cov, EigVal, EigVec);

	/*copy the Eigen values and eigen vectors*/
	for (int d = dim - ZeroCount; d<dim; d++)
	{
		EigVal.at<double>(d) = 0.0;
	}

	for (int d2, d = 0; d < dim; d++)
	{
		for (d2 = 0; d2 < dim; d2++)
		{
			pcDir[d].direction[d2] = EigVec.at<double>(d, d2);
		}

		pcDir[d].variance = EigVal.at<double>(d);
		pcDir[d].mean = innerProduct(dim, Mean, pcDir[d].direction);
	}

	//sort int Descending order of variance
	sort(pcDir, pcDir + dim);

#ifdef SHOW_PROGRESS
	timeEnd = GetCPUTime();
	cout << static_cast<int>(timeEnd - timeStart) / 1000 << " sec\t" << endl;
#endif

	//delete
	Cov.release();
	EigVal.release();
	EigVec.release();
	delete[] Mean;
	for (int d = 0; d < dim; ++d)
	{
		delete[] covariance[d];
	}
	delete[] covariance;
}

template <typename data_t>
bool PrincipalComponentAnalysis::executePCAcorelationCoefficient(
	int dim, 
	int num, 
	data_t** data
	)
{
	//次元数変わったらメモリ確保しなおし
	resetDimension(dim);

	//各基底の平均と分散共分散行列を得る
	double* Mean = new double[dim];
	Mat CorCoe(dim, dim, CV_64FC1);
	calCovarianceMatrix(dim, num, data, Mean, CorCoe);

	//正規化相関係数を得る
	calCorelationCoefficient(CorCoe);

	/*Calculate Eigen values and eigen vectors*/
	Mat EigVal, EigVec;
	EigenFunction(dim, Cov, EigVal, EigVec);
	EigVal.release();
	CorCoe.release();

	/*calculate Mean at Principal Component Space if need*/
	if (!PCAmean){
		PCAmean = new double[dim];
	}
	for (d = 0; d < dim; d++){
		PCAmean[d] = InnerProduct(&EigVec.at<double>(d, 0), Mean, dim);
	}
	delete[] Mean;

	double t;
	if (!Variance){
		Variance = new double[dim];
	}
	SortStructure_CV* SS = new SortStructure_CV[dim];
	for (d = 0; d < dim; d++){
		Variance[d] = 0;
		for (n = 0; n < num; n++){
			t = InnerProduct(&EigVec.at<double>(d, 0), data[n], dim) - PCAmean[d];
			Variance[d] += t*t;
		}
		SS[d].idx = d;
		SS[d].value = Variance[d] / num;
	}

	/*pass the Eigen values and eigen vectors*/
	qsort(SS, dim, sizeof(SortStructure), DesendIdx);
	if (!EigenVector){
		EigenVector = new double*[dim];
		for (d = 0; d < dim; d++){
			EigenVector[d] = new double[dim];
		}
	}
	for (d = 0; d < dim; d++){
		Variance[d] = SS[d].value;
		for (d2 = 0; d2 < dim; d2++){
			EigenVector[d][d2] = EigVec.at<double>(SS[d].idx, d2);
		}
	}
	delete[] SS;
	EigVec.release();

	return true;
}
#endif