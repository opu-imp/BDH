/**
* @file		HashTable.h
* @author	Tomokazu Sato
* @date		2015/05/05
*/

#ifndef __HashTable__
#define __HashTable__

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <memory.h>
using namespace std;

typedef unsigned collision_t;//!< type of collision
typedef char* address_t;	 //!< type of address

/**
* @brief a chain list
*/
struct bin_t
{
	collision_t collision;		 //!< collision
	address_t addressOfChainList;//!< head address of chain list
};

/**
* @brief hash table
*/
class HashTable
{

private:
	static const size_t collisionSize = sizeof(collision_t); //!< bytes of collision

	size_t entrySize;	 //!< bytes of entry
	size_t hashSize;	 //!< size of hash table
	size_t nEntry;
	address_t* hashTable;//!< hash table

public:

	/**
	* @brief default constructor
	*/
	HashTable()
		: entrySize(0)
		, hashSize(0)
		, nEntry(0)
		, hashTable(nullptr)
	{}

	size_t get_nEntry() const
	{
		return nEntry;
	}

	/**
	* @brief constructor
	*/
	HashTable(size_t entrySize, size_t hashSize)
		: entrySize(entrySize)
		, hashSize(hashSize)
	{
		initialize(entrySize, hashSize);
	}

	/**
	* @brief destructor
	*/
	~HashTable()
	{
		if (hashTable == nullptr)
		{
			return;
		}

		address_t* table_p = hashTable;
		for (size_t hv = 0; hv < hashSize; ++hv)
		{
			if (*table_p != nullptr)
			{
				free(*table_p);
			}
			++table_p;
		}

		delete[] hashTable;
		hashTable = nullptr;
	}

	/**
	* @brief initializer
	*/
	void initialize(
		size_t pointSize,//!< [in] sizeof(data_t)*dim;
		size_t hashSize	 //!< [in] hash size = 1<<bit;
		);

	/**
	* @brief bin getter
	*/
	void getBin(size_t hashKey, bin_t& bin) const
	{
		address_t table_p = hashTable[hashKey];
		if (!table_p)
		{
			bin.collision = 0;
			return;
		}
		bin.collision = *reinterpret_cast<collision_t*>(table_p);
		bin.addressOfChainList = table_p + collisionSize;
	}


	/**
	* @brief collision getter
	*/
	collision_t getCollision(
		size_t hashKey //!< [in] hash value
		) const
	{
		address_t table_p = hashTable[hashKey];
		if (table_p)
		{
			return (*reinterpret_cast<collision_t*>(table_p));
		}
		else
		{
			return 0;
		}
	}

	/**
	* @brief check is this hash value
	*/
	address_t getPointer(
		size_t hashKey//!< [in] hash value
		) const
	{
		return hashTable[hashKey];
	}

	/**
	* @brief check is hash value available
	* @return is hashTable[hashKey] active ?
	*/
	bool isEntried(
		size_t hashKey //!< [in] hash vlue
		)const 
	{
		if (hashTable[hashKey])
		{
			return true;
		}
		else
		{
			return false;
		}
	};


	/**
	* @brief Storing a point into hash table
	* @return is allocation complete ?
	*/
	bool storeEntry(
		size_t HashValue,	 //!< [in] hash value of point data
		const address_t point//!< [in] the point stored into hash table
		);

	/**
	* @brief allocation of hash table
	* @return is allocation complete ?
	*/
	bool allocTable(
		unsigned* collision//!< [in] collision list of the hash table
		);

	/**
	* @brief you have to take care of whether allocation of hash table has already completed
	*/
	void storeEntryWithoutAlloc(
		size_t HashValue,	//!< [in] hash value of point data
		const address_t data//!< [in] the point stored into hash table 
		);

	/**
	* @brief read hash table from binary file
	* @return is file open ?
	*/
	bool readTable(
		const string& tblFile//!< [in] file path
		);

	/**
	* @brief write hash table from binary file
	* @return is file open ?
	*/
	bool writeTable(
		const string& tblFile//!< [in] file path
		);

};






#endif