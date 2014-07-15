#ifndef _jvSampler_h_
#define _jvSampler_h_

#include <vector>
#include <gsl/gsl_rng.h>
#include "jvCommon.h"
#include "jvAlign.h"
#include <limits.h>

using namespace std;

typedef unsigned int uint;

enum CONV_FLAGS{
    NO_BORDER = 0,
    LEFT_BORDER = 1,
    RIGHT_BORDER = 2,
    LEFT_RIGHT_BORDER = 3    
};

struct jvMergePair{
    int i1;
    int i2;
    double sim;
};

class jvSampler{
    static const gsl_rng_type * T;
    static gsl_rng * r;
    static bool m_isMAP;
public:
    static void setMAP(bool map){m_isMAP = map;}
    static bool getMAP(){return m_isMAP;}
    static void init(unsigned long int s = 0);
    static void release();
    static uint geornd(double p);
    static uint intrnd(uint N);
    static double unirnd();
    static uint lmnrnd(const jvDVec& l, bool MAP=false);
    static double DPll(const jvDVec& v, double alpha);
    static vector<int> randomPermutation(uint N, bool identity=false);

//  functions to find the likelihood of a set of counts given the bases    
    static double logPy_base_observed(const jvBaseCounts& counts, jvBase base);
    static void   logPy_all_base_observed(const jvBaseCounts& counts, jvDVec& res);
    static double logPy_base_unobserved(const jvBaseCounts& counts, const jvDVec& bb_pssm=jvDVec(0));

    static double logPy_specific_lr(jvBaseCountsVec& counts, jvBaseVec& bases, uint left=0, uint right=UINT_MAX);

    static double logPy_specific(jvBaseCountsVec& counts, jvBaseVec& bases, int o=0, uint left=0, uint right= UINT_MAX);
    
    //static double logPy_specific(jvBaseCountsVec& counts, jvBaseVec& bases, int o);

    static double logPy_specific(jvBaseVec& read, jvBaseVec& bases, int o=0, uint left=0, uint right=UINT_MAX);
    
    static void sampleFromCounts(jvBaseCountsVec& counts, jvBaseVec& bases, const vector<int>& map, int left = 0, vector<int>* problematic = NULL);
    static void readToCounts(jvBaseVec& read, jvBaseCountsVec& counts);
    static bool isActive(){return (r != NULL);}
    static double sampleConcentration(jvDVec v, double alpha);
    
};



double log_sum_exp(const jvDVec& v);
jvDVec log(const jvDVec& v);
jvDVec norm(const jvDVec& v);        


#endif
