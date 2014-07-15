#ifndef _jvState_h_
#define _jvState_h_

#include "jvCommon.h"
#include "jvDataGenerator.h"
#include "jvAlign.h"
#include <set>
#include <map>
#include <assert.h>

using namespace std;

class jvInstance;
class jvReadAss;
class jvSequence;
class jvLikelihood;

typedef pair<jvSequence*, int> jvLocation;
typedef int jvKmer;
typedef vector<set<jvLocation> > jvHash;
static const set<jvLocation> emptyHashEntry;
static const set<jvInstance*> emptyInsSet;

enum GIBBS_POLICY{
  GIBBS_ALL = 0x00,
  GIBBS_FREEZE = 0x01,
  GIBBS_INCLUDE = 0x02,
  GIBBS_EXCLUDE = 0x04,
  GIBBS_NEW_SEQS_ONLY = 0x08,
  GIBBS_SNM_ONLY = 0x10,
  GIBBS_BOTH_ORIENTATIONS = 0x40,
  GIBBS_NO_CONSTRUCT = 0x80,
};

enum SPIKE_POLICY{
  SPIKE_HASH = 0x01,    
  SPIKE_NEW_SEQ = 0x02,
  SPIKE_LAST_TWO_SEQS = 0x04,  
  SPIKE_ALL_SEQS = 0x08,
  SPIKE_SPARSE = 0x10,
  SPIKE_DEFAULT = 0x01+0x02,
  SPIKE_SNM = 0x01+0x04
};

enum PRINT_OPTIONS{
    PRINT_INS = 1,
    PRINT_BASES = 2,
    PRINT_LL = 4,
    PRINT_ALL = 7,
    PRINT_DEBUG = 8,
    PRINT_ALL_DEBUG = 15,
    PRINT_READS = 16,
    PRINT_SUMMARY = 32,
    PRINT_SKELETON = 64
};

struct jvLogProb
{
    double partition;
    double choice;
    jvLogProb(double Z=0, double p=0):partition(Z),choice(p){}
    jvLogProb operator+(jvLogProb x){x.partition += partition; x.choice +=choice; return x;}
};

class jvParticle{
public:
	int o_t;
	int o_i;
	jvInstance* ins;
	int o_spike;
	//  alignment details
	set<jvAlignCache>* alignmentCache;
	jvParticle(){ alignmentCache = NULL; }
};

class jvLikelihood{
public:
    double s;
    double t;
    double o;
    double w;
    double y;
    double a;
    double ll(){return o+w+y+s+t+a;}
    jvLikelihood(){o=0; w=0; y=0; s=0; t=0; a=0;}
    void print(){
        cout.precision(2);
        cout<<fixed;
        cout<<"ll: "<<ll()<<" | o: "<<o<<", w: "<<w<<", y: "<<y<<", s: "<<s<<", t: "<<t<<", a: "<<a<<endl;}
    jvLikelihood operator+(jvLikelihood x){x.o += o; x.w +=w; x.y += y; x.s+=s; x.t +=t; x.a += a; return x;}
};


class jvSnmUndo{
public:
    bool isSplit;
    bool rightMoves;
    int offset;
    int back1;
    int back2;
    jvSequence* s;
    jvSequence* s2;
    double p_nst;
    double l_nst;
    void undo();
    jvSnmUndo():back1(0), back2(0), bases(0){}
    vector<int> bases;
    int ov_left;
    int ov_right;
};


struct jvSpikePolicy
{
    SPIKE_POLICY code;
    int bin;
    jvSpikePolicy(SPIKE_POLICY code_=SPIKE_DEFAULT) : code(code_),bin(0){}
};


//  READ

class jvReadAss{
    jvReadData* m_data;
    jvInstance* m_ins;
    int m_offset;
    int m_new_offset;  // this member is used by sampleInstance
    bool m_isAssigned;
    bool m_oneSpike;
    bool m_skip;

public:
    jvInstance* getIns(){return m_ins;}
    jvReadData* getData() {return m_data;}
    void setIns(jvInstance* ins, int offset=0);
    int getOffset() {return m_offset;}
    bool assign(bool isAssigned);
    bool isAssigned(){return m_isAssigned;}
    jvReadAss(jvReadData* data);
    void print(){}
    void revComplement();
    long getMem();
private:
//    void reg();
    friend class jvInstance;
    friend class jvState;
    friend class jvStateSerializer;
    friend class jvSequence;
};


//  INSTANCE

class jvInstance{
    int m_N; // #reads in instance
    int m_offset;
    jvSequence* m_seq;
    bool m_isAssigned;
    static uint WIDTH;
    friend class jvReadAss;
    friend class jvSequence;
    friend class jvState;
    friend class jvStateSerializer;
    set<jvReadAss*> m_reads;
public:
    static uint DEFAULT_READ_LENGTH;
    ~jvInstance(){}
    jvInstance(jvSequence* seq = NULL, int offset=0);
    int getN() {return m_N;}
    int getOffset() {return m_offset;}
    bool isAssigned(){return m_isAssigned;}
    jvSequence* getSeq(){return m_seq;}
    void setSeq(jvSequence* seq, int offset);
    void assign(bool isAssigned);
    void print();
    bool operator<(jvInstance& other){return (m_offset < other.m_offset);}
    int getLeft();
    int getRight();
    void revComplement();
    long getMem();

private:
    void updateCounts(int action, int count=1);
    void reg();
};

class jvEdit
{
public:
	int in;
	int dl;
	jvEdit(int i,int d):in(i),dl(d) {}
	jvEdit() {in=0;dl=0;}

	jvEdit& operator+=(jvEdit a) {in += a.in;dl += a.dl; return *this;}
	jvEdit& operator-=(jvEdit a) {assert(in>=a.in && dl>=a.dl);in -= a.in;dl -= a.dl; return *this;}
	void reset() {in=0;dl=0;}
};


class jvAlignmentJob
{
public:
	bool empty;
	jvReadData *read;
	jvSequence *seq;
	int offs;
	jvAlignmentJob() { empty = true; }
};

//  SEQUENCE


class jvSequence{
    int m_M;
    int m_left;
    int m_right;
    double m_p;
    double m_logP;
    double m_logInvP;
    long int m_sumO;
    jvBaseVec m_bases;
    vector<int> m_map;
    vector<jvEdit> m_editSummary;
    int m_insCount;
    int m_delCount;
    jvBaseCountsVec m_counts;
    bool m_isInstanceFrozen;
    vector<jvInstance*> m_ins;
    void pad(int left, int right);
    int pos(int offset) {return offset-m_left;}
    void updateP();
    jvLikelihood m_ll;
    jvDVec m_logPy;  // byproduct of running sampleO
    vector<int> m_spikes; // byproduct of running sampleInstance
    vector<int> m_spikes_bu; // byproduct of running sampleInstance
    vector<int> m_problematic;
    bool m_never_clean;
    bool m_inprocess;  // used in jvState::computeSpikes to tell if it was hit.
    bool m_getRight_dirty;
    int m_getRight;
    int m_bin;


    vector<set<jvInstance*> > m_ins2;
    double logPo(int offset, bool tail=false);
    jvHash* m_hash;
    void buildHash(jvHash* hash){m_hash = hash; updateHash();}
    void updateHash(int left=MIN_INT, int right=MAX_INT, int action=1);
    void computeSpikedInstances(vector<jvParticle>& spikeIns, jvDVec& logPty);
    bool getOverlap(jvSequence* seq, int offset,
    		int& shift_offset, int& ov_left, int&ov_right);
    friend class jvInstance;
    friend class jvState;
    friend class jvReadAss;
    friend class jvStateSerializer;
    friend class jvSnmUndo;
    static double pseudoA;
    static double pseudoB;
public:
    ~jvSequence(){}
    jvSequence();
    jvSequence(jvBaseVec* bb_vec);
    int getM() {return m_M;}
    jvBase& getBases(int offset){return m_bases[pos(offset)];}
    void print(PRINT_OPTIONS opt=PRINT_ALL);
    int getLeft();
    int getRight();
    int center();
    int shift(int shiftOffset);
    double getP() { return m_p;}
    jvLikelihood ll();
    jvSnmUndo split(int offset, bool rightMoves=true, jvSequence* new_seq=NULL);
    jvSnmUndo merge(jvSequence* seq, int offset);
    void merge_simple(jvSequence* seq, int offset);
    jvLikelihood mergedLl(jvSequence* seq, int offset, double* p);
    jvLikelihood splitLl(int offset, double gamma);
    static double ll_o(vector<int>& os, double* p=NULL);
    bool isEqual(jvSequence* s);
    bool isSkeleton(){return m_never_clean;}

    void revComplement();

    void addEdit(int ofs, vector<char>* edit, bool print=false);
    void summarizeEdits();
    void resetEditSummary();
    bool checkAllOffsets();
    void applyEdit(vector<char> edit, int st);

    void printCounts(); // for debug
    bool insertAtOffset(int ofs); // bool for debug purposes
    void deleteAtOffset(int ofs);
    jvBaseVec buildAlignBuff(int L, int o);
    jvAlignmentJob getNextAlignmentJob(int &m, const jvAlignmentJob templateJob);
    pair<jvInstance*, jvInstance*> getGapInstances(bool all = false);

    long getMem();
    int getBin(){return m_bin;}
    void setBin(int bin){m_bin = bin;}
    vector<jvReadAss*> getReads();    
private:
    int pruneSpikes(jvSpikePolicy policy);
    void updateReadCounts(jvBaseCountsVec& newCounts, int offset, int action);
    void updateReadCounts(jvBaseVec& newCounts, int offset, int action);
    void updateInstanceCounts(jvInstance* ins, int offset, int action, int count=1);

    void setFreeze(bool freeze);
    void sampleBases(bool greedy_override=false);
    double logPy(jvBaseVec& read, int offset);
    double logPy(int offset);
    jvLogProb sampleO(jvInstance* ins, jvBaseVec* read, jvParticle* particle, bool debug_print = false);
    double sampleInternalO(int offset, int* res);
    double logPy_specific_bb(jvBaseCountsVec& counts, int o);
    double logPy_specific_bb(jvBaseVec& read, int o);
    double logPy_align(jvBaseVec& read, int* o, set<jvAlignCache>* alignmentCache=NULL, bool print=false);
    double logPy_align_ins(jvInstance* ins, int o, bool change_reads=false, bool print=false);
    void init();  // called by c-tor

    void propagateDeleteToRead(int ofs, jvReadAss* ra);
    void propagateInsertToRead(int ofs, jvReadAss* ra);
};



// STATE

class jvState{
protected:
    vector<jvReadData>& m_data;
    vector<jvReadAss> m_reads;
    vector<jvInstance*> m_ins;
    vector<jvSequence*> m_seq;
    double m_gamma;
    double m_alpha;
    jvHash m_hash;
    bool m_withInstances;
    bool m_batch;
    int m_sum_spikes;
    int m_count_reads;
    vector<jvSequence*> m_spiked_seqs;
    int m_nSpikes; // internal counter updated by sampleInstance to measure no. of spikes
    clock_t m_start;
    jvSnmUndo snm(uint t1, uint t2, bool changeState=true);
    void setTemperature();
    vector<int> m_read_stack;
    int m_total_splits;
    int m_total_merges;
    int m_total_failures;

public:
    vector<jvReadAss*> popReads(int n=1);
    bool isThereAGap();
    int candI,candJ; // a hack to force a particular merge
    int seq_choice, cutoff; // a hack to force a particular gibbs choice in SNM
    double prop;
    void saveStatus(const char* filename, uint snm_count);
    void buildHash();
    void checkHash();
    jvHash& getHash(){return m_hash;}
    void printHash();

    jvState(vector<jvReadData>& data, bool empty=false, int seed=0);
    ~jvState();
    jvState* copy();
    vector<jvReadAss>& getReads(){return m_reads;}
    vector<jvInstance*>& getIns(){return m_ins;}
    vector<jvSequence*>& getSeq(){return m_seq;}
    jvInstance* newIns(int seqId = -1, int offset = 0);
    jvInstance* newIns(jvSequence* p_seq = NULL, int offset = 0);
    jvSequence* newSeq(jvSequence* seq=NULL);
    void print(PRINT_OPTIONS opt=PRINT_ALL);
    void cleanInstances();
    int cleanSequences();
    bool isEqual(jvState* st);
    void center();
    void flipSequences();
    void unassignEnds();
    void sampleBases(bool greedy_override=false);
    void freezeInstances(bool freeze);
    jvLikelihood ll();
    jvLogProb sampleInstance(jvInstance* ins,  jvReadAss* read,
                          jvSpikePolicy policy = jvSpikePolicy(), bool debug_print=false);
    uint computeSpikes(jvBaseVec* ptr, jvSpikePolicy policy);
    uint setSpikesAtZero();
    jvLogProb gibbsReads(vector<jvReadAss*> reads=vector<jvReadAss*>(), GIBBS_POLICY gibbs_policy = GIBBS_ALL);
    void gibbsInstances();
    jvDVec getM(bool sparse = false, bool with_last=false);
    jvDVec getN();
    void sampleGamma();
    void sampleAlpha();
    bool splitAndMerge3(bool out_print=false);

    static bool splitAndMerge(jvState* st, bool out_print=false);
    static uint K;
    static int endAl;
    static bool alignment_only;
    static double garbage_thresh;
    static bool construct;
    static bool both_orientations;
    static bool all_reads;
    static bool with_temperature;
    static double temperature;
    static bool with_geometric;

    void summarizeEdits();
    void checkAllOffsets();
    static double computeCountsAlignmentScore(vector<char> edit, int local_start, jvBaseVec seqI, jvBaseVec seqJ, jvBaseCountsVec countsI, jvBaseCountsVec countsJ, vector<char> &editI, vector<char>  &editJ);
    bool alignEnds(jvSequence *seqI, jvSequence *seqJ, int dirJ, bool out_print = false);
    bool updateHashHits(jvKmer &Kmer, jvSequence* seqI,map<jvSequence*,int> &hits,set<jvSequence*> &seen,map<jvSequence*,bool> &touched, bool rev);
    bool findAndAlignEnds();
    bool findAndSplitEnds(bool all = false);
    double cleanAlignmentCache(int cutoff);
    void initAlignmentFromCache();


    //  get/set
    double getGamma(){return m_gamma;}
    void setGamma(double gamma){cout<<"Setting gamma="<<gamma<<endl; m_gamma = gamma;}
    double getAlpha(){return m_alpha;}
    void setAlpha(double alpha){m_alpha = alpha;}
    void setGeomPseudos(double pseudoA, double pseudoB)
    {jvSequence::pseudoA = pseudoA; jvSequence::pseudoB = pseudoB;}
    void setInstanceWidth(uint width){cout<<"Setting instance width = "<<width<<endl; jvInstance::WIDTH = width;}
    void setWithInstances(bool withInstances){
        cout<<"Setting with instances "<<(withInstances ? "true" : "false")<<endl;
        m_withInstances = withInstances;   if (!withInstances) setInstanceWidth(1);}
    void setEndAl(int _endAl){cout<<"Setting endAL="<<_endAl<<endl; endAl = _endAl;}
    void setBothOrientations(bool val){both_orientations = val; cout<<"will "<<(val? "":"not ")<<"check both orientations\n";}
    void setAlignmentOnly(bool _construct = false, double _garbage_thresh = log(0.25));
    void setAuto();
    void setBatch(bool batch){m_batch = batch;}
    void setAllReads(bool val){all_reads = val; cout<<"will "<<(val? "":"not ")<<"gibbs all the reads\n";}
    void setAnnealing(bool ann){with_temperature = ann; if (ann) cout<<"will "<<(ann? "":"not ")<<"perform annealing\n";}
    void setGeometric(bool geom){with_geometric = geom;  cout<<"will "<<(geom? "":"not ")<<"use geometric.\n";}
    void setErrorModel(double noise, double indels);
    long getMem();
    void dumpHashCounts();
    friend class jvStateSerializer;
    bool sanityEditString();
    bool sanityGetRight();

};

#endif
