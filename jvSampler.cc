#include "jvSampler.h"
#include <gsl/gsl_randist.h>
#include<gsl/gsl_sf_gamma.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <assert.h>
#include <iterator>
#include <vector>
#include <algorithm>

//#include "DisjointSets.h"
#include "jvAlign.h"


using namespace std;
jvDVec norm(const jvDVec& v)
{
    jvDVec r(v.size());
    double s = sum(v);
    if (s == 0) return r;
    for (uint i=0; i<v.size(); i++)
        r[i] = v[i] / s;
    return r;
}


double log_sum_exp(const jvDVec& v)
{
//    double v_max = *max_element(v.begin(), v.end());
    double v_max = max(v);
    double res = 0;

//  Approximation for faster computation: exp(-10) ~= 0
    for (uint i=0; i<v.size(); i++)
    {
    	double z =  v[i]-v_max;
        res += (z < -10 ? 0 : exp(z) );
    }
    return log(res) + v_max;

}


jvDVec log(const jvDVec& v)
{
    jvDVec r(v.size());
    if (v.size() == 0) return r;
    double s = *min_element(v.begin(), v.end());
    if (s < 0) return r;
    for (uint i=0; i<v.size(); i++)
        r[i] = log(v[i]+1e-100);
    return r;
}


const gsl_rng_type* jvSampler::T = NULL;
gsl_rng * jvSampler::r;
bool jvSampler::m_isMAP = false;

void jvSampler::init(unsigned long int s)
{

    /* create a generator chosen by the
       environment variable GSL_RNG_TYPE */

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set (r, s);
}

void jvSampler::release()
{
    gsl_rng_free(r);
    r = NULL;
}


uint jvSampler::lmnrnd(const jvDVec& l, bool MAP)
{
    assert(l.size() != 0);
    if (l.size() == 1) return 0;

    jvDVec::const_iterator it = std::max_element(l.begin(), l.end());

    // create weights so that max(l_norm) = 1.
    double l_max = *it;
    vector <double> l_norm(l.size());
    for (uint i=0; i<l.size(); i++)
        l_norm[i] = exp(l[i]-l_max);
    
    if (m_isMAP || MAP)
    {
        int candidates = 0; int winner = -1;
        for (uint i=0; i<l.size(); i++)
            if (l_norm[i]< 0.99)
                l_norm[i] = 0;
            else {candidates++; winner = i;}
        if (candidates == 1)
            return winner;
        else assert(candidates>0);
    }
    
    gsl_ran_discrete_t* g = gsl_ran_discrete_preproc (l_norm.size(), &l_norm[0]);
    uint res = gsl_ran_discrete (r, g);
    gsl_ran_discrete_free (g);
    return res;
}




double jvSampler::DPll(const jvDVec& v, double alpha)
{
    double sum_v = 0;
    double sum_lnGamma_v = 0;
    double count_v = 0;
    for (uint i=0; i<v.size(); i++)
    {
        if (v[i] == 0) continue;
        sum_v += v[i];
        sum_lnGamma_v += gsl_sf_lngamma(v[i]);
        count_v++;
    }
    if (alpha > 1e+10)
    {
    	return (count_v - sum_v )* log(alpha) + sum_lnGamma_v;
    }
    else
    {
    	return count_v * log(alpha) + sum_lnGamma_v + gsl_sf_lngamma(alpha) - gsl_sf_lngamma(alpha + sum_v);
    }
}


// returns P(counts|base)
double jvSampler::logPy_base_observed(const jvBaseCounts& counts, jvBase base)
{
    double res = 0;
    for (jvBase j=0; j<(jvBase)nBases; j++)
        res += counts[j]*(j == base ? errorModel.logInvNoise() : errorModel.logNoise());
    return res;
}


//  Given a set of counts for *one* base location, return the (un-normalized) log-posterior
//  for the base in that location
//  in other words:  P(counts|base) for every one of the four bases
void jvSampler::logPy_all_base_observed(const jvBaseCounts& counts, jvDVec& res)
{
    res.resize(nBases);
    for (jvBase j=0; j<(jvBase) nBases; j++)
        res[j] = logPy_base_observed(counts, j);
}


// returns P(counts)
double jvSampler::logPy_base_unobserved(const jvBaseCounts& counts, const jvDVec& bb_pssm)
{
//    bool empty = (*max_element(counts.begin(), counts.end()) == 0);
    double N = (double) nBases;
    int max_counts = max(counts);
    int sum_counts = sum(counts);

    if (max_counts == 0){
        return 0;
    }

    // Optimize for "pure" counts (when only one entry has non-zero counts)
//  (1-3noise)^counts + 3*noise^counts /4
    if ((bb_pssm.size() == 0) && (max_counts == sum_counts))
        return log(  pow ( (1-(N-1)*errorModel.noise()), max_counts) + ((N-1)*pow(errorModel.noise(),max_counts)))+logBasePrior;

    jvDVec res(counts.size());
    logPy_all_base_observed(counts, res);

    if (bb_pssm.size() == 0)
    	return  log_sum_exp(res)+logBasePrior;   // = 0 if counts are 0
    else
    	return  log_sum_exp(bb_pssm+res);

}



void jvSampler::sampleFromCounts(jvBaseCountsVec& counts, jvBaseVec& bases, const vector<int>& map, int left, vector<int>* problematic)
{
    jvDVec res(nBases);
    for (uint i=0; i<counts.size(); i++)
    {
        if (map[i] == 0)
            bases[i] = 4;
        else
        {
            if (getMAP() == true)
                for (int j=0; j<nBases; j++)
                    res[j] = counts[i][j];
            else
                logPy_all_base_observed(counts[i], res);
            
            bases[i] = lmnrnd(res);
            if (problematic != NULL)
            {
                double conf = 1;
                if ( (conf < 0.99) || (map[i] <= 2) || (i+1 == counts.size())
                     || (map[i+1] == 0) || (i == 0) || (map[i-1] == 0))
                    problematic->push_back(i+left);                    
            }
        }
    }

}



//  likelihood counts given that they were taken from offset o in bases,
//  taking into account only the entries in [left, right).
double jvSampler::logPy_specific(jvBaseCountsVec& counts, jvBaseVec& bases,  int o, uint left, uint right)
{

    assert(left < right);
    assert(left < bases.size());
    right = min(right, (uint)bases.size());
//    left = max(left, 0);  // not necessary since left is uint

    double res = 0;
    for (uint i=0; i<counts.size(); i++)
    {
        if ( (o+i < left) || (o+i >= right) || (bases[o+i] == (int)nBases) )
            res += logPy_base_unobserved(counts[i]);
        else
            res += logPy_base_observed(counts[i], bases[o+i]);
    }
    return res;

}

//  likelihood counts taking into account only the entries (in bases and in counts) in [left, right).
double jvSampler::logPy_specific_lr(jvBaseCountsVec& counts, jvBaseVec& bases, uint left, uint right)
{
    right = min((uint)bases.size(), right);
    assert(left < right);
    assert(right <= bases.size());
    assert(right <= counts.size());

    double res = 0;
    for (uint i=left; i<right; i++)
    {
        if (bases[i] == (int)nBases)
            res += logPy_base_unobserved(counts[i]);
        else
            res += logPy_base_observed(counts[i], bases[i]);
    }
    return res;
}



void jvSampler::readToCounts(jvBaseVec& read, jvBaseCountsVec& counts)
{
    counts.resize(read.size(), emptyBaseCount);
    for (uint i=0; i<read.size(); i++)
        for (jvBase j=0; j<4; j++)
            counts[i][j] = (j == read[i] ? 1 : 0);
}



double jvSampler::logPy_specific(jvBaseVec& read, jvBaseVec& bases, int o, uint left, uint right)
{
    jvBaseCountsVec counts;
    readToCounts(read, counts);
    return logPy_specific(counts, bases, o,left, right);
}


double jvSampler::unirnd()
{
    return gsl_rng_uniform(r);
}

uint jvSampler::intrnd(uint N)
{
    return gsl_rng_uniform_int(r, N);
}


uint jvSampler::geornd(double p)
{
    return gsl_ran_geometric(r, p)-1;
}


vector<int> jvSampler::randomPermutation(uint N, bool identity)
{
    vector<int> res(N);
    for(uint i=0; i<N; i++) res[i] = i;

    if (!identity)
        random_shuffle(res.begin(), res.end());

    return res;
}


double jvSampler::sampleConcentration(jvDVec v, double alpha)
{
    uint k = v.size(); //  total number of tables
    double n = sum(v); //  number of guests
    double a=1, b=1;   //  Gamma distribution prior

// 1) sample mu from Beta(alpha+1, n)
// 2) sample alpha from pi G(a+k, 1/b-log mu) + (1-pi) G(a+k-1, 1/b-log mu)
//  pi/(1-pi) = ( a+k-1 ) / ( n(1/b-log mu) )

    assert( !((n==0) || (a==0) || (b==0)));

    double mu = gsl_ran_beta (r, alpha+1, n);

    double j = 1/b-log(mu);
    double tmp =  ( a+k-1 ) / (n * j) ;
    double Pi = tmp/(tmp+1);

    double coin = gsl_rng_uniform (r);

    if (coin < Pi)
        return gsl_ran_gamma(r, a+k, 1/j);
    else
        return gsl_ran_gamma(r, a+k-1, 1/j);

}


int argmax(const jvDVec& v, int exclude, const vector<int>& I)
{
    int N = v.size();
    assert(N>0);
    assert(N == (int)I.size());
    assert( exclude >= -1 && exclude < N);

    int ix = 0;
    while ( ix<N && (/*ix == exclude ||*/ I[ix] != ix))
        ix++;

    assert(ix != N);
    for (uint i=ix+1; i<N; i++)
    {
        if ( I[i] != i)
            continue;
        if (v[i] > v[ix])
            ix = i;
    }
    return ix;    
}

