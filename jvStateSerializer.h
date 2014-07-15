#ifndef _jvStateSerializer_H_
#define _jvStateSerializer_H_

#include "jvState.h"

class jvStateSerializer
{
public:
	static void dump_fasta(jvState *jvs, const char *fname);
	static void dump(jvState *jvs, const char *fname, bool dumpDataToFile=true, const char *fclust=NULL);
	static jvState* read(const char *fname, vector<jvReadData>& readData, bool dataInFile=true);
	static jvState* read_compressed(const char *fname, vector<jvReadData>& readData, bool dataInFile=true);
	
};

#endif /*_jvStateSerializer_H_*/
