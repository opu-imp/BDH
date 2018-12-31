
#include <Subspace.h>
#include <NearestNeighbor.h>
#include <k_means.h>

int Subspace::dim;

void Subspace::setParameters(const baseset_t& baseSet)
{
	bit = baseSet.bit;
	subDim = baseSet.subDim;
	subHashSize = baseSet.k;
	variance = baseSet.variance;

	cellVariance = new double[subHashSize];
	memcpy(cellVariance, baseSet.cellVariance, sizeof(double)*subHashSize);

	centroid = new double*[subHashSize];
	for (int h = 0; h < subHashSize; ++h)
	{
		centroid[h] = new double[subDim];
		memcpy(centroid[h], baseSet.centroid[h], sizeof(double)*subDim);
	}

	base = new double*[subDim];
	for (int d = 0; d < subDim; ++d){

		base[d] = new double[dim];
		memcpy(base[d], baseSet.base[d].direction, sizeof(double)*dim);
	}
}

void Subspace::clear(){
	if (base != nullptr){
		for (int sd = 0; sd < subDim; ++sd){
			if (base[sd] != nullptr){
				delete[] base[sd];
			}
		}
		delete[] base;
		base = nullptr;
	}

	if (hashKey != nullptr){
		delete[] hashKey;
		hashKey = nullptr;
	}

	if (cellVariance != nullptr){
		delete[] cellVariance;
		cellVariance = nullptr;
	}

	if (centroid != nullptr)
	{
		for (int h = 0; h < subHashSize; ++h)
		{
			if (centroid[h] != nullptr)
			{
				delete[] centroid[h];
			}
		}
		delete[] centroid;
		centroid = nullptr;
	}
}

double Subspace::getDistanceToCentroid(double* PCAquery, int centroidIndex) const
{
	double t = cellVariance[centroidIndex];
	double* centroid_p = centroid[centroidIndex];
	double* centroid_p_end = centroid_p + subDim;
	for (; centroid_p != centroid_p_end;)
	{
		t += NORM((*PCAquery++) - (*centroid_p++));
	}

	return t;
}
