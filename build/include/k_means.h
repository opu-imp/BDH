/**
* @file    k_means.h
* @author  T.Sato
* @date    2015.05.05
* @version 1.0
*/

#ifndef __k_means__
#define __k_means__

#include <iostream>
#include <limits.h>
#include <float.h>
#include <memory.h>

#include "NearestNeighbor.h"

#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

/**
* @brief k-means
*/
template<typename data_t, typename centroid_t>
class K_Means
{

private:
	int dim;	 //!< dimension
	int K;		 //!< number of centroid
	double error;//!< quantization error

	unsigned* centroid_size;//!< number of point included in corespond cluster
	double* centroid_error; //!< quantization error of correspond cluster
	centroid_t** centroid;  //!< centroid

	double CONVERGENCE_RATE;  //!< convergence rate 
	unsigned CONVERGENCE_LOOP;//!< convergence loop

public:

	/**
	* @brief default constructor
	*/
	K_Means();

	/**
	* @brief destructor
	*/
	~K_Means();

	/**
	* @brief error getter
	* @return error
	*/
	double getError()
	{
		return error;
	}

	/**
	* @brief centroid getter
	* @return centroid
	*/
	double** const getCentroid()const
	{
		return centroid;
	}

	/**
	* @brief calculate centroid
	*/
	void calclateCentroid(
		int dim,	   //!< [in] dimension
		int num,	   //!< [in] number of sample
		data_t** point,//!< [in] sample point set
		int K		   //!< [in] number of centroid
		);

private:

	void setParam(int dim, int K);

	void IniCentScala(int num, data_t** point);

	void updateMinMax(double& Min, double& Max, const double& val);

	/**
	* @brief find maximam value
	*/
	unsigned findMaxIndex(
		int num, 
		double* MinDisTable
		);

	void IniCent_PlaPla(
		int num,
		data_t** point);

	void updateMinDisTable(
		centroid_t* centroid,
		int num,
		data_t** point, 
		double* MinDisTable
		);

	double CalCentScala(
		int num, 
		data_t* point
		);

	double CalCent(
		int num,
		data_t** point
		);

	int ResetCent(
		long num, 
		data_t** point
		);

};

template<typename data_t, typename centroid_t>
K_Means<data_t, centroid_t>::K_Means()
	: dim(0)
	, K(0)
	, error(0.0)
	, centroid_size(nullptr)
	, centroid_error(nullptr)
	, centroid(nullptr)
	, CONVERGENCE_RATE(1.0e-5)
	, CONVERGENCE_LOOP(UINT_MAX)
{}

template<typename data_t, typename centroid_t>
K_Means<data_t, centroid_t>::~K_Means()
{
	if (centroid != nullptr)
	{
		for (int i = 0; i < K; ++i)
		{
			delete[] centroid[i];
		}
		delete[] centroid;
		centroid = nullptr;
	}

	delete[] centroid_error;
	centroid_error = nullptr;

	delete[] centroid_size;
	centroid_size = nullptr;
}

template<typename data_t, typename centroid_t>
void K_Means<data_t, centroid_t>::calclateCentroid(
	int dim,	   //!< dimension
	int num,	   //!< number of sample
	data_t** point,//!< sample point set
	int K		   //!< number of centroid
	)
{

	setParam(dim, K);

	if (dim == 1)
	{
		IniCentScala(num, point);
	}
	else
	{
		IniCent_PlaPla(num, point);
	}

	CalCent(num, point);
}

template<typename data_t, typename centroid_t>
void K_Means<data_t, centroid_t>::setParam(int dim, int K)
{

	if (this->dim == dim &&
		this->K == K)
	{
		return;
	}
	this->~K_Means();


	this->dim = dim;
	this->K = K;
	centroid_size = new unsigned[K];
	centroid_error = new double[K];
	centroid = new double*[K];
	for (int i = 0; i < K; ++i){
		centroid[i] = new double[dim];
	}
}

template<typename data_t, typename centroid_t>
void K_Means<data_t, centroid_t>::IniCentScala(int num, data_t** point)
{

	double Min = DBL_MAX;
	double Max = -DBL_MAX;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (int n = 0; n < num; n++)
	{
#ifdef _OPENMP
#pragma omp critical
#endif
			{
				updateMinMax(Min, Max, *point[n]);
			}
	}

	double delta = (Max - Min) / K;
	for (int c = 0; c < K; c++)
	{
		*centroid[c] = (centroid_t)(Min + delta*(c + 0.5));
	}
}

template<typename data_t, typename centroid_t>
void K_Means<data_t, centroid_t>::updateMinMax(double& Min, double& Max, const double& val)
{
	if (val < Min)
	{
		Min = val;
	}

	if (val > Max)
	{
		Max = val;
	}
}

template<typename data_t, typename centroid_t>
void K_Means<data_t, centroid_t>::IniCent_PlaPla(
	int num,
	data_t** point)
{

	/* 重心を求める */
	double* Mean = new double[dim]; /* 重心 */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (int d = 0; d < dim; d++)
	{
		Mean[d] = 0;
		for (int n = 0; n < num; n++)
		{
			Mean[d] += point[n][d];
		}
		Mean[d] /= num;
	}

	//最も近いセントロイドまでの距離
	double* MinDisTable = new double[num];

	//各サンプルの重心までの距離をもとめる
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (int n = 0; n < num; n++)
	{
		MinDisTable[n] = Distance(dim, Mean, point[n]);
	}

	/* 重心から最も遠い点を求め、
	* それを１つ目のセントロイドとする
	*/
	unsigned MinDisIndex_max = findMaxIndex(num, MinDisTable);

	for (int d = 0; d < dim; d++)
	{
		centroid[0][d] = point[MinDisIndex_max][d];
	}

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (int n = 0; n < num; n++)
	{
		MinDisTable[n] = DBL_MAX;
	}

	/* その他のセントロイドを求める */
	for (int c = 1; c < K; c++)
	{

		updateMinDisTable(centroid[c - 1], num, point, MinDisTable);

		/* choose the furthest point(except CentIdx)
		* it is a next centroid
		*/
		const int MinDisIndex_max = findMaxIndex(num, MinDisTable);

		for (int d = 0; d < dim; d++)
		{
			centroid[c][d] = point[MinDisIndex_max][d];
		}
	}

	/* メモリ解放 */
	delete[] Mean;
	delete[] MinDisTable;
}

/* 最も遠いインデックスを探す
*/
template<typename data_t, typename centroid_t>
unsigned K_Means<data_t, centroid_t>::findMaxIndex(int num, double* MinDisTable)
{

	unsigned MaxIndex = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (int n = 1; n<num; n++)
	{

#ifdef _OPENMP
#pragma omp critical
#endif
			{
				if (MinDisTable[n] > MinDisTable[MaxIndex]){
					MaxIndex = n;
				}
			}
	}

	return MaxIndex;
}

template<typename data_t, typename centroid_t>
void K_Means<data_t, centroid_t>::updateMinDisTable(centroid_t* centroid, int num, data_t** point, double* MinDisTable){

	double dist;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(dist)
#endif
	for (int n = 0; n < num; n++)
	{
		//check the distance to the nearest centroid
		dist = Distance(dim, centroid, point[n], MinDisTable[n]);

		if (dist < MinDisTable[n])
		{
			MinDisTable[n] = dist;
		}
	}

}

template<typename data_t, typename centroid_t>
double K_Means<data_t, centroid_t>::CalCentScala(
	int num, 
	data_t* point)
{
	int loop = 0;
	int* nCP = new int[K];
	centroid_t* centroid2 = new centroid_t[K];
	for (int c = 0; c < K; c++)
	{
		centroid2[c] = 0;
		nCP[c] = 0;
	}

	error = DBL_MAX;
	double error0 = 0.0;

	while (error0 != error)
	{

		loop++;
		error0 = error;
		error = 0;

		//calculate new centroid
		for (int n = 0; n < num; n++)
		{

			int index = 0;
			double dist;
			double terror = fabs(point[n] - centroid[0]);
			for (int i = 1; i < K; ++i)
			{
				dist = fabs(point[n] - centroid[i]);
				if (dist < terror)
				{
					terror = dist;
				}
			}

			error += terror;
			centroid2[index] += point[n];
			nCP[index]++;
		}

		error /= num;
		for (int c = 0; c < K; c++){
			centroid[c] = centroid2[c] / nCP[c];
			centroid2[c] = 0;
			nCP[c] = 0;
		}
		if (error <= 1.0e-5){
			break;
		}
	}

	delete[] centroid2;

	return error;
}

template<typename data_t, typename centroid_t>
double K_Means<data_t, centroid_t>::CalCent(
	int num,
	data_t** point
	){

	centroid_t** tmpCent;
	centroid_t** lastCentroid = new centroid_t*[K];
	for (int i = 0; i < K; ++i){
		lastCentroid[i] = new centroid_t[dim];
	}
	//エラーの管理を実装せにゃ
	double rate = DBL_MAX;
	error = DBL_MAX;

	for (unsigned loop = 0;
		rate > CONVERGENCE_RATE &&
		loop < CONVERGENCE_LOOP &&
		error > 1.0e-5
		; ++loop)
	{

		//swap
		tmpCent = lastCentroid;
		lastCentroid = centroid;
		centroid = tmpCent;

		const double error0 = error;
		for (int i = 0; i < K; ++i)
		{
			centroid_size[i] = 0;
			centroid_error[i] = 0.0;

			for (int d = 0; d < dim; ++d)
			{
				centroid[i][d] = 0.0;
			}
		}

		//calculate new centroid
		point_t<centroid_t> NNcent;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(NNcent)
#endif
		for (int n = 0; n < num; ++n)
		{
			data_t* point_n_itr = point[n];
			data_t* point_n_itr_end = point_n_itr + dim;
			NearestNeighbor(dim, K, lastCentroid, point_n_itr, NNcent);
			size_t index = NNcent.index;
			double* centroid_index_itr = centroid[index];

#ifdef _OPENMP
#pragma omp critical
#endif
			{
				for (; point_n_itr != point_n_itr_end; ++point_n_itr)
				{
					*centroid_index_itr += *point_n_itr;
					++centroid_index_itr;
				}
				++centroid_size[index];
				centroid_error[index] += NNcent.distance;
			}
		}

		error = 0;
		for (int c = 0; c < K; c++)
		{
			error += centroid_error[c];
			centroid_error[c] /= centroid_size[c];

			if (centroid_size[c] == 0){

				int CentID = ResetCent(num, point);

				for (int d = 0; d < dim; d++)
				{
					centroid[c][d] = point[CentID][d];
				}
			}
			else{
				for (int d = 0; d < dim; d++)
				{
					centroid[c][d] /= centroid_size[c];
				}
			}
		}

		error /= num;
		rate = (error0 - error) / error0;
	}

	for (int i = 0; i < K; ++i){
		delete[] lastCentroid[i];
	}
	delete[] lastCentroid;


	return error;
}

template<typename data_t, typename centroid_t>
int K_Means<data_t, centroid_t>::ResetCent(long num, data_t** point){

	double distanceMax = DBL_MAX;
	long indexMax = 0;

	point_t<centroid_t> NNcent;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(NNcent)
#endif
	for (long n = 0; n < num; n++){
		NearestNeighbor(dim, K, centroid, point[n], NNcent);

#ifdef _OPENMP
#pragma omp critical
#endif
		{
			if (NNcent.distance > distanceMax){
				distanceMax = NNcent.distance;
				indexMax = n;
			}
		}
	}

	return indexMax;
}

#endif