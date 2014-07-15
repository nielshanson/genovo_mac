#ifndef _jvFastaParser_H_
#define _jvFastaParser_H_

#include <string>
#include "jvState.h"

using namespace std;

class jvStringReadData
{
public:
	string read;
	string meta;

	jvStringReadData() {};
	jvStringReadData(jvReadData);
	jvStringReadData(string inMeta,string inRead) {meta = inMeta; read = inRead;};
	jvReadData convertToReadData();
};

class jvFastaParser
{	
public:	
	static vector<jvStringReadData>* parseStringReadData(const char* name);		
	static void parseReadData(const char *name, vector<jvReadData>& readData);	
	static vector<jvStringReadData>* convertFromReadData(vector<jvReadData>& readData);
	static void dump(vector<jvStringReadData>* data,const char * fname);
        static void dump_fastq(vector<jvStringReadData>* data, const char * fname);	
};

#endif /*_jvFastaParser_H_*/
