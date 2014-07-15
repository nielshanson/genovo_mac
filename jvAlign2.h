#ifndef JVALIGN2_H_
#define JVALIGN2_H_
#include "jvCommon.h"

struct jvAlignCosts
{
    void print() {cout<<"indel="<<indel<<" miss="<<miss<<" unobserved="<<unobserved_match<<" pad="<<pad_insert<<endl; }
    jvAlignCosts(){}
    jvAlignCosts(double _indel, double _miss,	double _hit, double _unobserved_match, double _pad_insert):
        indel(_indel), miss(_miss), hit(_hit), unobserved_match(_unobserved_match), pad_insert(_pad_insert){}
    
    double indel;
    double miss;
    double hit;
    double unobserved_match;
    double pad_insert;
    
};

class jvAlign2
{
	int maxalen,maxrlen,maxeditlen;
	int m_band;
	int m_lenA, m_lenR;
	jvBase *m_A, *m_R;
	double p[5][5];  	// match/mismatch
	double sx,dx,dy; // insertion/deletion
	int si[3], sj[3];	
	double* fm;	
	unsigned char *bm;
	char* edit;
	int editLen;
	int local_start;
	bool m_scoreOnly;
		
	double score;		
	double alignFast();
	int backtraceFast();
public:
	void init(jvAlignCosts costs);
	jvAlign2(int maxalen, int maxrlen, jvAlignCosts costs);		
	int align(jvBase* A, int lenA, jvBase* R, int lenR, int band=MAX_INT);
	int align(jvBaseVec& A, jvBaseVec& R, int band=MAX_INT)
        {return align(&A[0], A.size(), &R[0], R.size(), band);}
	~jvAlign2();
	char* getEdit(int* outEditLen, bool local=false);
	double getScore() {return score;}
	void applyEdit(jvBase* inStr, int inLen, jvBase* outStr, int *outStrLen, bool local=false);
	int getLocalStart(){return local_start;}
};
#endif /*JVALIGN2_H_*/
