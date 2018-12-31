/**
* @file    NNpointIO.h
* @author  T.Sato
* @date    2015.05.03
* @version 1.0
*/

#include "point.h"

/**
* @brief save the result of reverse nearest neighbor search.
* @return is file open ?
*/
template <typename data_t>
bool SaveReverseNearestNeighbor(
	const string& path,				  //!< [in] file path
	unsigned Qnum,					  //!< [in] number of query point set
	int K,						  //!< [in] number of nearest neighbors
	vector<point_t<data_t>>* RNNpoint,//!< [in] revers nearest neighbors
	double time = 0.0,				  //!< [in] query time
	size_t memory = 0				  //!< [in] work memory
	);

/**
* @brief load the result of reverse nearest neighbor search.
* @return is file open ?
*/
template <typename data_t>
bool LoadReverseNearestNeighbor(
	const string& path,				  //!< [in] file path
	unsigned& Qnum,					  //!< [out] number of query point set
	int& K,						  //!< [out] number of nearest neighbors
	vector<point_t<data_t>>* RNNpoint,//!< [out] revers nearest neighbors
	double& time,					  //!< [out] query time
	size_t& memory					  //!< [out] work memory
	);

/**
* @brief save the result of reverse nearest neighbor search.
* @return is file open ?
*/
template <typename data_t>
bool SaveNearestNeighbor(
	const string& path, //!< [in]
	unsigned Qnum, //!< [in] number of query point set
	int K, //!< [in] number of nearest neighbors
	point_t<data_t>** NNpoint, //!< [in] nearest neighbors
	double time = 0.0, //!< [in]
	size_t memory = 0 //!< [in]
	);

/**
* @brief load the result of reverse nearest neighbor search.
* @return is file open ?
*/
template <typename data_t>
bool LoadNearestNeighbor(
	const string& path, //!< [in]
	unsigned& Qnum, //!< [out] number of query point set
	int& K, //!< [out] number of nearest neighbors
	point_t<data_t>**& NNpoint, //!< [out] nearest neighbors
	double& time, //!< [out]
	size_t& memory //!< [out]
	);

///////////////// implementation after here /////////////////

template <typename data_t>
bool SaveNearestNeighbor(
	const string& path,
	unsigned Qnum,
	int K,
	point_t<data_t>** NNpoint,
	double time,
	size_t memory){

	ofstream ofs(path);
	if (ofs.is_open() == false){
		cerr << "failed to open " << path << endl;
	}

	ofs << K << "\t" << Qnum << endl;
	for (unsigned q = 0; q < Qnum; ++q){

		for (int s = 0; s < K; ++s){
			ofs << NNpoint[q][s].index << "\t" << NNpoint[q][s].distance << "\t";
		}
	}
	ofs << time << "\t" << memory << endl;

	ofs.close();

	return true;
}

template <typename data_t>
bool LoadNearestNeighbor(
	const string& path,
	unsigned& Qnum,
	int& K,
	point_t<data_t>**& NNpoint,
	double& time,
	size_t& memory){

	ifstream ifs(path);
	if (ifs.is_open() == false){
		cerr << "failed to open " << path << endl;
		return false;
	}

	ifs >> K >> Qnum;
	NNpoint = new point_t < data_t >*[Qnum];
	for (unsigned q = 0; q < Qnum; ++q){

		NNpoint[q] = new point_t<data_t>[K];
		for (int s = 0; s < K; ++s){
			ifs >> NNpoint[q][s].index >> NNpoint[q][s].distance;
		}
	}
	ifs >> time >> memory;

	ifs.close();

	return true;
}


template <typename data_t>
bool SaveReverseNearestNeighbor(
	const string& path,
	unsigned Qnum,
	int K,
	vector<point_t<data_t>>* RNNpoint,
	double time,
	size_t memory){

	ofstream ofs(path);
	if (ofs.is_open() == false){
		cerr << "failed to open " << path << endl;

		return false;
	}

	size_t size;
	ofs << K << "\t" << Qnum << endl;
	for (unsigned q = 0; q < Qnum; ++q){

		size = RNNpoint[q].size();
		ofs << size << "\t";
		for (size_t s = 0; s < size; ++s){
			ofs << RNNpoint[q][s].index << "\t" << RNNpoint[q][s].distance << "\t";
		}
		ofs << endl;
	}
	ofs << time << "\t" << memory << endl;

	ofs.close();

	return true;
}

template <typename data_t>
bool LoadReverseNearestNeighbor(
	const string& path,
	unsigned& Qnum,
	int& K,
	vector<point_t<data_t>>* RNNpoint,
	double& time,
	size_t& memory){

	ifstream ifs(path);
	if (ifs.is_open() == false){
		cerr << "failed to open " << path << endl;

		return false;
	}

	size_t size;
	ifs >> K >> Qnum;
	for (unsigned q = 0; q < Qnum; ++q){

		size = RNNpoint[q].size();
		ifs >> size;
		RNNpoint[q].resize(size);
		for (size_t s = 0; s < size; ++s){
			ifs >> RNNpoint[q][s].index >> RNNpoint[q][s].distance;
		}
	}
	ifs >> time >> memory;

	ifs.close();

	return true;
}
