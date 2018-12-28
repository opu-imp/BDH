/**
* @file    PCA.h
* @author  T.Sato
* @date    2015.09.05
* @version 1.0
*/

#ifndef __SUPER_PCA__
#define __SUPER_PCA__

#include <string>
#include <math.h>
using namespace std;

/**
* @brief handle eigen vector and correlation information
*/
struct PC_t{
	int		dim;	  //!< dimension
	double  mean;	  //!< mean at base direction
	double  variance; //!< variance at base direction = eigen value
	double* direction;//!< eign vector

	bool operator < (const PC_t& obj)
	{
		return this->variance > obj.variance;
	}
};

/**
* @brief Pricipal Component Analysis
*/
class superPCA
{
protected:
	int dim;
	int ZeroCount;
	PC_t* pcDir;

public:

	//コンストラクタ
	superPCA()
		: dim(0)
		, ZeroCount(0)
		, pcDir(nullptr)
	{}

	//デストラクタ
	virtual ~superPCA()
	{
		if (pcDir != nullptr)
		{
			for (int d = 0; d < dim; ++d)
			{
				delete[] pcDir[d].direction;
			}
			delete[] pcDir;
			pcDir = nullptr;
		}
	}

	/**
	* @brief pcDir getter
	* @return pcDir
	*/
	const PC_t* getPCdir()const
	{
		return pcDir;
	}

	/**
	* @brief project src vector to ret vector
	*/
	template<typename srcType, typename retType>
	void projectToPCspace(
		srcType* src,
		retType* ret
		);

	/**
	* @brief save result of PCA
	*/
	bool savePCA(
		const string& path
		);

	/**
	* @brief load result of PCA
	*/
	bool loadPCA(
		const string& path
		);

protected:

	template<typename data_t>
	double innerProduct(
		int dim,
		data_t* data,
		double* dir
		);

	/**
	* @brief realloc memory if dimension is changed
	*/
	void resetDimension(int dim);

	/**
	*/

	template<typename data_t>
	void calculateCovarianceMatrix(
		int dim, size_t num, data_t** data,
		double* mean, double** covariance);

	template<typename data_t>
	void calculateCoefficientMatrix(
		int dim, size_t num, data_t** data,
		double* mean, double** coefficient);
};

/**
* @brief inner product
*/
template<typename data_t>
double superPCA::innerProduct(
	int dim,
	data_t* data,
	double* dir
	)
{
	double val = 0.0;
	for (int d = 0; d < dim; ++d)
	{
		val += data[d] * dir[d];
	}
	return val;
}

/**
* @brief project src vector to ret vector
*/
template<typename srcType, typename retType>
void superPCA::projectToPCspace(
	srcType* src,
	retType* ret
	)
{
	for (int d = 0; d < dim; ++d)
	{
		ret[d] = static_cast<retType>(innerProduct(dim, src, pcDir[d].direction));
	}
}

template<typename data_t>
void superPCA::calculateCovarianceMatrix(
	int dim, size_t num, data_t** data,
	double* mean, double** covariance)
{

#ifdef _OPENMP
#pragma omp parallel
#endif
		{
			/* calculate Mean */
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
			for (int d = 0; d < dim; d++)
			{
				mean[d] = 0;
				for (unsigned n = 0; n < num; n++)
				{
					mean[d] += data[n][d];
				}
				mean[d] /= num;
			}

			/*Covariance*/
#ifdef _OPENMP
#pragma omp for schedule(guided)
#endif
			for (int d = 0; d<dim; d++)
			{
				double cov;
				double mean_d2;
				double mean_d = mean[d];
				data_t** data_n;
				data_t** data_n_End;
				for (int d2 = d; d2<dim; d2++)
				{
					mean_d2 = mean[d2];
					cov = 0;

					data_n_End = data + num;
					for (data_n = data; data_n != data_n_End; ++data_n)
					{
						cov += ((*data_n)[d] - mean_d)*((*data_n)[d2] - mean_d2);
					}
					covariance[d][d2] = covariance[d2][d] = cov / num;
				}
			}
		}

		bool* flag = new bool[dim];
		for (int d = 0; d < dim; d++)
		{
			if (covariance[d][d] < 1.0e-10)
			{
				flag[d] = false;
			}
			else
			{
				flag[d] = true;
			}
		}

		for (int d = 0; d < dim; d++)
		{
			for (int d2 = d + 1; d2<dim; d2++)
			{
				if (flag[d] && flag[d2]){
					if (covariance[d2][d] * covariance[d2][d]>(1 - 1.0e-10)*(covariance[d][d] * covariance[d2][d2]))
					{
						flag[d2] = false;
					}
				}
			}
		}

		ZeroCount = 0;
		for (int d = 0; d<dim; d++){
			if (!flag[d])
			{
				++ZeroCount;
			}
		}
		delete[] flag;
}

template<typename data_t>
void superPCA::calculateCoefficientMatrix(
	int dim, size_t num, data_t** data,
	double* mean, double** coefficient)
{
	calculateCovarianceMatrix(dim, num, data, mean, coefficient);

	for (int d = 0; d < dim; ++d)
	{
		for (int d2 = d + 1; d2 < dim; ++d2)
		{
			coefficient[d][d2] /= sqrt(coefficient[d][d] * coefficient[d2][d2]);
			coefficient[d2][d] = coefficient[d][d2];
		}

		coefficient[d][d] = 1.0;
	}
}
#endif