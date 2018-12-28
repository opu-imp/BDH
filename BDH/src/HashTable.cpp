#include "HashTable.h"
#include <iostream>
using namespace std;

void HashTable::initialize(size_t entrySize, size_t hashSize)
{
	//ƒƒ‚ƒŠ‰Šú‰»
	this->~HashTable();

	this->entrySize = entrySize;
	this->hashSize = hashSize;

	try
	{
		hashTable = new address_t[hashSize];
	}
	catch (bad_alloc)
	{ // —áŠO bad_alloc ‚ð‚±‚±‚ÅŽó‚¯Žæ‚é
		cerr << " Out of memory. size:" << hashSize << endl;
		cerr << __FILE__ <<"\t" << __LINE__ << endl;
		abort();
	}
	catch (...){ // ‚»‚êˆÈŠO‚Ì—áŠO‚Í‚±‚±‚ÅŽó‚¯Žæ‚é
		cerr << "error something." << endl;
		cerr << __FILE__ << "\t" << __LINE__ << endl;
		abort();
	}

	memset(hashTable, 0, sizeof(address_t)*hashSize);
}

/* allocation of hash table */
bool HashTable::allocTable(
	collision_t* collision)
{
	initialize(entrySize, hashSize);

	for (size_t hv = 0; hv < hashSize; ++hv)
	{
		if (collision[hv] == 0)
		{
			hashTable[hv] = nullptr;

		}
		else
		{
			hashTable[hv] = static_cast<address_t>(
				malloc(collisionSize + entrySize*collision[hv])
				);
			
			if (hashTable[hv] == nullptr)
			{
				cerr << "failed to malloc memory of hashTable[" << hv << "]" << endl;
				cerr << collisionSize << "\t" << entrySize << "\t" << collision[hv] << endl;
				cerr << "malloc size is " <<
					collisionSize + entrySize*collision[hv] << endl;
				return false;
			}

			*reinterpret_cast<collision_t*>(hashTable[hv]) = 0;
		}
	}

	return true;
}

/* read hash table from binary file */
bool HashTable::readTable(const string& tblFile)
{
	ifstream ifs(tblFile, ios::in | ios::binary);
	if (ifs.is_open() == false)
	{
		cerr << "failed to open " << tblFile.c_str() << endl;
		return false;
	}
	ifs >> entrySize >> hashSize >> nEntry;
	cout << "entrySize = " << entrySize << endl;
	cout << "hashSize = " << hashSize << endl;
	cout << "number of Entry = " << nEntry << endl;

	//initialize
	initialize(entrySize, hashSize);

	cout << "read collision." << endl;
	collision_t* collision = new collision_t[hashSize];
	ifs.read((char*)collision, collisionSize*hashSize);

	cout << "read table." << endl;
	size_t progressStep = hashSize / 20;
	size_t progressPoint = progressStep;
	nEntry2 = 0;
	collision_t* collision_p = collision;
	for (size_t h = 0; h < hashSize; ++h,++collision_p)
	{
		if (*collision_p > 0)
		{
			if (*collision_p <= collisionMax)
			{
				hashTable[h] = reinterpret_cast<address_t>(malloc(collisionSize + entrySize*(*collision_p)));
				if (hashTable[h] == nullptr)
				{
					cerr << "failed to malloc memory of hashTable[" << h << "]" << endl;
					cerr << collisionSize << "\t" << entrySize << "\t" << (*collision_p) << endl;
					cerr << "malloc size is " << collisionSize + entrySize*(*collision_p) << endl;
					cerr << __FILE__ << "\t" << __LINE__ << endl;
					exit(__LINE__);
				}
				*reinterpret_cast<collision_t*>(hashTable[h]) = (*collision_p);
				ifs.read(hashTable[h] + collisionSize, entrySize*(*collision_p));
				nEntry2 += (*collision_p);
			}
			else
			{
				ifs.seekg(entrySize*(*collision_p), ios::cur);
				hashTable[h] = nullptr;
			}
			/*
			cout << "hashValue " << h << "\t" << collision[h] << endl;
			char* datac = (char*)(hashTable[h] + collisionSize);
			for (collision_t c = 0; c < collision[h]; ++c)
			{
				float* data = (float*)datac;
				double leng = 0;
				for (int d = 0; d < 60; ++d)
				{
					leng += (data[d])*(data[d]);
				}

				//cout << leng << endl;
				
				if (leng<0.99 || leng>1.01)
				{
					cout << h << "\t" << c << "\t" << leng << "\t" << *(unsigned*)(datac + entrySize - sizeof(unsigned)) <<endl;
					for (int d = 0; d < 60; ++d)
					{
						cout << data[d] << "\t";
					}
					cout << endl;
				}

				datac += entrySize;
			}
			*/

			if (h > progressPoint)
			{
				cout << h * 100 / hashSize << " % done." << endl;
				progressPoint += progressStep;
			}
		}
		else
		{
			hashTable[h] = nullptr;
		}
	}

	delete[] collision;

	ifs.close();
	return true;
}

/* read hash table from binary file */
bool HashTable::readTableSkipData(const string& tblFile)
{
	ifstream ifs(tblFile, ios::in | ios::binary);
	if (ifs.is_open() == false)
	{
		cerr << "failed to open " << tblFile.c_str() << endl;
		return false;
	}
	ifs >> entrySize >> hashSize >> nEntry;
	cout << "entrySize = " << entrySize << endl;
	cout << "hashSize = " << hashSize << endl;
	cout << "number of Entry = " << nEntry << endl;

	//initialize
	initialize(entrySize, hashSize);

	cout << "read collision." << endl;
	collision_t* collision = new collision_t[hashSize];
	ifs.read((char*)collision, collisionSize*hashSize);
	ifs.close();

	cout << "read table." << endl;
	size_t progressStep = hashSize / 20;
	size_t progressPoint = progressStep;
	nEntry2 = 0;
	collision_t* collision_p = collision;
	for (size_t h = 0; h < hashSize; ++h, ++collision_p)
	{
		if (*collision_p > 0)
		{
			hashTable[h] = reinterpret_cast<address_t>(malloc(collisionSize));
			*reinterpret_cast<collision_t*>(hashTable[h]) = (*collision_p);
			if (h > progressPoint)
			{
				cout << h * 100 / hashSize << " % done." << endl;
				progressPoint += progressStep;
			}
		}
		else
		{
			hashTable[h] = nullptr;
		}
	}

	delete[] collision;

	return true;
}

/* write hash table into binary file */
bool HashTable::writeTable(const string& tblFile) const
{
	long long long_hashSize = hashSize;
	collision_t* collision = new collision_t[hashSize];
	for (long long h = 0; h < long_hashSize; ++h)
	{
		if (hashTable[h] == nullptr)
		{
			collision[h] = 0;
		}
		else
		{
			collision[h] = *reinterpret_cast<collision_t*>(hashTable[h]);
		}
	}

	ofstream ofs(tblFile, ios::out | ios::binary);
	if (ofs.is_open() == false)
	{
		cerr << "failed to open " << tblFile.c_str() << endl;
		return false;
	}
	ofs << entrySize << "\t" << hashSize << "\t" << nEntry;
	ofs.write((char*)collision, collisionSize*hashSize);
	for (size_t h = 0; h < hashSize; ++h)
	{
		if (collision[h] > 0)
		{
			ofs.write(
				reinterpret_cast<char*>(hashTable[h] + collisionSize),
				entrySize*(collision[h]));
		}
	}

	delete[] collision;

	ofs.close();

	return true;
}

/* you have to take care of whether allocation of hash table has already completed */
void HashTable::storeEntryWithoutAlloc(
	size_t HashValue,	/* hash value of point data */
	const address_t entry		/* the point stored into hash table */
	)
{
	//update collision
	address_t table_p = hashTable[HashValue];
	collision_t collision = *reinterpret_cast<collision_t*>(table_p);
	++(*reinterpret_cast<collision_t*>(table_p));

	//write point on memory
	table_p += collisionSize + entrySize*collision;
	memcpy(table_p, entry, entrySize);

	++nEntry;
}

/* Storing a point into hash table */
inline bool HashTable::storeEntry(
	size_t HashValue,	/* hash value of point data */
	const address_t point		/* the point stored into hash table */
	)
{
	//getCollision
	collision_t collision = 0;
	address_t table_p = hashTable[HashValue];
	if (table_p != nullptr){
		collision = *reinterpret_cast<collision_t*>(table_p);
	}

	//realloc heap memory
	table_p = static_cast<address_t>(
		realloc(table_p, collisionSize + entrySize*(collision + 1)));
	if (!table_p)
	{
		return false;
	}

	hashTable[HashValue] = table_p;

	//update collision
	++(*reinterpret_cast<collision_t*>(table_p));
	//write point on memory
	table_p += collisionSize + entrySize*collision;
	memcpy(table_p, point, entrySize);

	++nEntry;
	return true;
}

