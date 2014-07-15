#ifndef _jvDataGenerator_h_
#define _jvDataGenerator_h_

#include <iostream>
#include <vector>
#include <set>
#include "jvCommon.h"
#include "jvAlign.h"
//#include "values.h"
#include "limits.h"

using namespace std;

static const int DG_RANDOM = 0;
static const int DG_FIXED = 1;

class jvReadData{
private:
    jvBaseVec m_bases;
    jvBaseVec m_raw_bases;
    jvBaseVec m_rev_bases;
    bool m_isReversed;
    jvBaseVec m_correct;
    vector<char> m_edit;
    set<jvAlignCache> m_alignment_cache;
    set<jvAlignCache> m_alignment_cache_rev;
    double m_alignment_score;
    int m_seqId;
    int m_offset;
    int m_siteId;
    bool m_orientation;
public:
    bool modified;
    bool touched;
    int getSeqId(){return m_seqId;}
    int getOffset(){return m_offset;}
    int getSiteId(){return m_siteId;}
    jvBaseVec& getBases() {return m_bases;}
    int getBasesLength() const {return m_bases.size();}

    jvBaseVec* getRawBases() {return (m_isReversed ? &m_rev_bases : &m_raw_bases);}
    vector<char>* getEdit() {return &m_edit;}
    void resetAlignment(void);
    set<jvAlignCache>* getAlignmentCache() { return (m_isReversed?  &m_alignment_cache_rev : &m_alignment_cache); }
    int cleanAlignmentCache(int cutoff);

    void updateAlignment(bool local = true, bool print = false, jvAlign* aln = &jvAlign::aligner);
    jvReadData(jvBaseVec bases, int seqId=0, int offset=0, int siteId=0, bool orientation=false, jvBaseVec correct=emptyBaseVec);
    jvReadData(jvBaseVec bases, jvBaseVec rawBases, vector<char> edit, double score);
    jvReadData complement();
    void complementSelf();
    void print();
    bool printBoth();
    void setReversed(bool isReversed){m_isReversed = isReversed;}
    bool getReversed(){return m_isReversed;}

    void setCorrect(jvBaseVec correct) {m_correct = correct;}
    jvBaseVec getCorrect() { return m_correct;}
    bool isCorrect();

    double getScore() {return m_alignment_score;}
    void setScore(double score) {m_alignment_score = score;}

    void insertAtOffset(int pos);
    void deleteAtOffset(int pos);
    long getMem();
};

class jvDataGenerator{
private:
    jvBaseVec m_bases;
    int m_seqId;
    int m_siteId;
    vector<jvReadData> m_reads;
    int m_readLength;
    uint m_var;
    int m_interval;
    int m_rand_interval;
    int m_readsPerSite;
    double m_noise;
    double m_insProb;
    double m_delProb;
    double m_missing;
    double m_revProb;
    jvBaseVec m_backbone;
public:
    jvDataGenerator(int random_seed);
    void setReadLength(int length, uint var=0);
    void setInterval(int interval, int isRand){ m_interval = interval; m_rand_interval = isRand;}
    void setReadsPerSite(int readsPerSite, int isRand) { m_readsPerSite = readsPerSite;}
    void setNoise(double noiseProb=0, double insProb=0, double delProb=0, double missing=0) { m_noise = noiseProb; m_insProb = insProb; m_delProb = delProb; m_missing = missing;}
    void randomizeSequence(int nBases);
    void randomizeBackbone(int nBases, int nDistance);
    void randomizeVariant(jvBaseVec mut); //  insers the pattern 'mut' every nDistance letters
    void setSequence(jvBaseVec& input);
    void setSequence(char* strSeq);
    void setOrientation(double reverseProb){m_revProb = reverseProb;}
    void addToReads();
    int nSeq(){return m_seqId;}
    int nSites(){return m_siteId;}
    void print();
    void printSequence();
    jvBaseVec* getSeq() {return &m_bases;}
    void saveSequence(char* name, bool append = false);
    void shuffleReads();
    jvBaseVec* getBackbone() {return &m_backbone;}
    vector<jvReadData> getReads(){ return m_reads; }
};



#endif
