
#include <stdlib.h>

#include <FileHandle.h>

void openErrorHandle(
	const ofstream& fs,
	const string& path,
	const string& file,
	int line)
{
	if (fs.is_open() == false)
	{
		cerr << "failed to open " << path << endl;
		cerr << "at file : " << file << endl;
		cerr << "at line : " << line << endl;
		exit(line);
	}

}

void openErrorHandle(
	const ifstream& fs,
	const string& path,
	const string& file,
	int line)
{
	if (fs.is_open() == false)
	{
		cerr << "failed to open " << path << endl;
		cerr << "at file : " << file << endl;
		cerr << "at line : " << line << endl;
		exit(line);
	}

}

void fileErrorHandle(
	bool flag,
	const string& path,
	const string& file,
	int line){
	
	if (flag == false)
	{
		cerr << "failed to open " << path << endl;
		cerr << "at file : " << file << endl;
		cerr << "at line : " << line << endl;
		exit(line);
	}
}

bool isFileExist(const string& path)
{
	ifstream ifs(path, ios::in);
	return ifs.is_open();
}
