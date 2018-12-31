#include "HashTable.h"

void HashTable::initialize(size_t entrySize, size_t hashSize)
{
	//ƒƒ‚ƒŠ‰Šú‰»
	this->~HashTable();

	this->entrySize = entrySize;
	this->hashSize = hashSize;

	try{
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
		return false;
	}
	ifs >> entrySize >> hashSize >> nEntry;

	initialize(entrySize, hashSize);

	address_t* table_p = hashTable;
	collision_t collision;
	for (size_t h = 0; h < hashSize; ++h)
	{
		ifs.read(reinterpret_cast<char*>(&collision), collisionSize);
		if (collision > 0)
		{
			*table_p = reinterpret_cast<address_t>(malloc(collisionSize + entrySize*collision));
			if (*table_p == nullptr)
			{
				cerr << "failed to malloc memory of hashTable[" << h << "]" << endl;
				cerr << collisionSize << "\t" << entrySize << "\t" << collision << endl;
				cerr << "malloc size is " << collisionSize + entrySize*collision << endl;
				cerr << __FILE__ << "\t" << __LINE__ << endl;
				exit(__LINE__);
			}
			*reinterpret_cast<collision_t*>(*table_p) = collision;
			ifs.read(*table_p + collisionSize, entrySize*collision);
		}
		else
		{
			*table_p = nullptr;
		}
		++table_p;
	}

	return true;
}

/* write hash table into binary file */
bool HashTable::writeTable(const string& tblFile) 
{
	ofstream ofs(tblFile, ios::out | ios::binary);
	if (ofs.is_open() == false)
	{
		return false;
	}
	ofs << entrySize << "\t" << hashSize << "\t" << nEntry;

	address_t* table_p = hashTable;
	collision_t collision;
	for (size_t h = 0; h < hashSize; ++h)
	{
		if (*table_p == nullptr)
		{
			collision = 0;
			ofs.write(
				reinterpret_cast<char*>(&collision),
				collisionSize);
		}
		else{
			collision = *reinterpret_cast<collision_t*>(*table_p);
			ofs.write(
				reinterpret_cast<char*>(*table_p),
				collisionSize + entrySize*collision);
		}
		++table_p;
	}
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

