/**
* @file    FileHandle.h
* @author  T.Sato
* @date    2015.05.03
* @version 1.0
*/

#ifndef __FILE_HANDLE__
#define __FILE_HANDLE__

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

/**
* @brief check statement of ofstream. if opening file is faile, exit.
*/
void openErrorHandle(
	const ofstream& fs, //!< [in] file stream
	const string& path, //!< [in] file path of file stream
	const string& file, //!< [in] cpp file name which open the file
	int line //!< [in] the line of the cpp file where call this function 
	);

/**
* @brief check statement of ofstream. if opening file is faile, exit.
*/
void openErrorHandle(
	const ifstream& fs, //!< [in] file stream
	const string& path, //!< [in] file path of file stream
	const string& file, //!< [in] cpp file name which open the file
	int line //!< [in] the line of the cpp file where call this function
	);

/**
* @brief check statement of ofstream. if opening file is faile, exit.
*/
void fileErrorHandle(
	bool flag, //!< [in] is file open ?
	const string& path, //!< [in] file path of file stream
	const string& file, //!< [in] cpp file name which open the file
	int line //!< [in] the line of the cpp file where call this function
	);

/**
* @brief read point set
* @return is file open ?
*/
template<typename data_t>
bool readBinary(
	const string& path, //!< [in] file stream
	int& dim, //!< [out] dimension of point
	unsigned& num, //!< [out] number of points
	data_t**& data //!< [out] point set
	);

bool isFileExist(const string& path);

/* implementation of template function after here */

template<typename data_t>
bool readBinary(
	const string& path, 
	int& dim, 
	unsigned& num, 
	data_t**& data)
{

	ifstream ifs(path, ios::in | ios::binary);
	if (ifs.is_open() == false)
	{
		return false;
	}

	ifs.read((char*)&dim, sizeof(dim));
	ifs.read((char*)&num, sizeof(num));
	data = new data_t*[num];

	const int dSize = sizeof(data_t)*dim;
	for (unsigned n = 0; n < num; n++){

		data[n] = new data_t[dim];
		ifs.read((char*)data[n], dSize);
	}
	ifs.close();

	return true;
}

template<typename data_t>
void deleteDataArray(unsigned num, data_t**& data)
{
	for (unsigned n = 0; n < num; ++n)
	{
		delete[] data[n];
	}
	delete[] data;
}

#endif