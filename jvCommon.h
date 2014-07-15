#ifndef _jvCommon_h_
#define _jvCommon_h_

#include <vector>
#include <assert.h>
#include <iostream>
#include <math.h>
//#include <values.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>

using namespace std;

typedef unsigned int uint;
const uint nBases = 4;

const double logBasePrior = log(1.0/nBases);
const double logTwo = log(2);

typedef unsigned char jvBase;
typedef vector<jvBase> jvBaseVec;
typedef vector<int> jvBaseCounts; // always length 4
typedef vector<jvBaseCounts> jvBaseCountsVec;
typedef vector<double> jvDVec;
static const jvBaseCounts emptyBaseCount(nBases,0);
static const jvBaseVec emptyBaseVec(0);

//const double iT = 1;
const int MAX_INT = 327670000;
const int MIN_INT = -327680000;
const double PI = 3.141592653589793;

class jvErrorModel
{
   double m_noise;
   double m_indel;
   double m_log_noise;
   double m_log_indel;
   double m_log_inv_noise;
public:
   void setParams(double noise, double indel)
   {
	m_noise = noise; 
	m_indel = indel;
   	m_log_noise = log(noise);
   	m_log_indel = log(indel);
   	m_log_inv_noise = log(1-(nBases-1)*noise);
	cout<<"base insertion prob: "<<indel<<"\tbase deletion prob: "<<indel<<"\tmismatch prob: "<<noise<<endl;
   }
   double noise(){return m_noise;}
   double indel(){return m_indel;}
   double logIndel(){return m_log_indel;}
   double logNoise(){return m_log_noise;}
   double logInvNoise(){return m_log_inv_noise;}
   jvErrorModel(){cout<<"init jvError Model\n"; setParams(0.01, 0.01);}
};

extern jvErrorModel errorModel;

template <class T>
ostream& operator<<(ostream& os, vector<T> vec)
{
    for (uint i=0; i<vec.size(); i++)
        os<<vec[i]<<" ";
    os<<endl;
    return os;
}


template <class T>
T sum(const vector<T> v)
{
    T r = 0;
    for (uint i=0; i<v.size(); i++)
        r += v[i];
    return r;
}

template <class T>
T max(const vector<T>& v)
{
    T r = v[0];
    for (uint i=1; i<v.size(); i++)
        if (v[i]>r) r = v[i];
    return r;
}


template <class T>
vector<T>& operator+=(vector<T>& v1, const vector<T>& v2)
{
    assert(v1.size() == v2.size());
    vector<T> r(v1.size());
    for (uint i=0; i<v1.size(); i++)
        v1[i]+=v2[i];
    return v1;
}


template <class T>
vector<T> operator+(const vector<T>& v1, const vector<T>& v2)
{
    assert(v1.size() == v2.size());
    vector<T> r(v1.size());
    for (uint i=0; i<v1.size(); i++)
        r[i] = v1[i]+v2[i];
    return r;
}

template <class T>
vector<T> operator-(const vector<T>& v1, const vector<T>& v2)
{
    assert(v1.size() == v2.size());
    vector<T> r(v1.size());
    for (uint i=0; i<v1.size(); i++)
        r[i] = v1[i]-v2[i];
    return r;
}

template <class T>
vector<T> operator+(const vector<T>& v1, const T c)
{
    vector<T> r(v1.size());
    for (uint i=0; i<v1.size(); i++)
        r[i] = v1[i]+c;
    return r;
}

template <class T>
vector<T> operator*(const vector<T>& v1, const T c)
{
    vector<T> r(v1.size());
    for (uint i=0; i<v1.size(); i++)
        r[i] = v1[i]*c;
    return r;
}


template <class T>
bool comp_ref(T* a, T* b)
{
    return (*a < *b);
}

void print_joni(jvBaseVec* vec);
void complement_bases(jvBaseVec& bases);
jvBase complement_base(jvBase base);
void print_edit(vector<char>& edit);
string get_edit(vector<char>& edit);
void set_edit(string str, vector<char>& edit);


jvBaseCounts complement_count(jvBaseCounts val);
void complement_counts(jvBaseCountsVec& bases);

#endif
