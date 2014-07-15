#ifndef JVALIGN_H_
#define JVALIGN_H_
#define HEIGHT 612
#define MAXALEN 1500
#define MAXRLEN 1500
#define MAXALIGNEDLEN (MAXALEN+MAXRLEN)+1
#define MAXEDITLEN    (MAXALEN+MAXRLEN)+1

#include "jvCommon.h"
#include "jvAlign2.h"

#include <set>



class jvAlignCache {
public:
	uint m_match[100];
	short m_len;
	int last;
	double m_score;
	vector<char> m_edit;
	jvBaseVec m_bases;
	int localStart;

	void setLast(int iter) { last = iter; }
	jvAlignCache(jvBase* first, int len);
	bool operator==(const jvAlignCache &other);
	bool operator!=(const jvAlignCache &other) { return !(*this == other); }

};
extern bool operator<(const jvAlignCache& left,const jvAlignCache& right);

static jvAlignCosts def_costs     (0, 0, 0, 0, 0);


class jvAlign
{
    jvAlign2 inner_aligner;
    int localStart;
    
    int maxalen,maxrlen,maxr1len;
    jvBase *buffA;
    int buffAlen;
    jvBase *aligned;
    char *edit;
    int alignedLen;
    int editLen;
    double score;
    int m_cache_hits;
    int m_cache_misses;

public:
    static jvAlign aligner;
    static jvAlign merge_aligner;
    static jvAlign edit_aligner;
    static int extend;
    static bool with_cache;
    
    jvAlign(int maxalen, int maxrlen, jvAlignCosts acosts=def_costs);
    void setCosts(jvAlignCosts costs){inner_aligner.init(costs);}
    double align(jvBaseVec *A, jvBaseVec *R, bool print = false);
    double align(jvBaseVec*A, int start, int end, jvBaseVec *R, set<jvAlignCache>* alignmentCache, bool print = false);
    double align(jvBase*A,int lenA,jvBase*R,int lenR,jvBase*R1, int *outLen,char*editR=NULL, int*outLenEditR=NULL);
    double alignFast(jvBase*A,int lenA,jvBase*R,int lenR,jvBase*R1, int *outLen,char*editR=NULL, int*outLenEditR=NULL);
    jvBase* getAligned();
    int getAlignedLen();
    char* getEdit();
    int getEditLen();
    double getScore();
    
    double getCacheHits(bool reset = false);
    double getAvgCacheSize(bool reset=false);
    
    int getLocalStart(){return localStart;}
    
    ~jvAlign();
};

#endif /*JVALIGN_H_*/
