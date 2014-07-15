//#include <values.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <iostream>
#include <assert.h>


#include "jvAlign.h"
#include "jvState.h"
#include "jvCommon.h"

using namespace std;

extern int iteration;

jvAlign jvAlign::aligner(MAXALEN,MAXRLEN);
jvAlign jvAlign::merge_aligner(2*MAXALEN,MAXRLEN);
jvAlign jvAlign::edit_aligner(MAXALEN,MAXRLEN);

int jvAlign::extend = 5;
bool jvAlign::with_cache = true;


jvAlignCache::jvAlignCache(jvBase* first, int len)
{
	int val;
	int i;
	int j=0;
	for (int i=0;i<len;i++) assert(*(first+i) <= nBases);

	for (i=0;i<len-len%10;i+=10)
	{
		val = (((int)*(first+i)  )<<27) + (((int)*(first+i+1))<<24)  + (((int)*(first+i+2))<<21)  + (((int)*(first+i+3))<<18)  + (((int)*(first+i+4))<<15) +
		      (((int)*(first+i+5))<<12) + (((int)*(first+i+6))<< 9)  + (((int)*(first+i+7))<< 6)  + (((int)*(first+i+8))<< 3)  + (((int)*(first+i+9)));
		m_match[j] = val;
		j++;
	}
	val = 0;
	for (;i<len;i++)
	{
		val = (val<<3)+((int)*(first+i));
	}
	assert(j<100); // supports reads up to size 1000;
	m_match[j] = val;
	m_len = j+1;
}

bool jvAlignCache::operator==(const jvAlignCache& other)
{
	if (this->m_len != other.m_len)
	{
		return false;
	}
	for (int i=m_len-1;i>=0;i--)
	{
		if (this->m_match[i] != other.m_match[i])
		{
			return false;
		}
	}
	return true;
}
bool operator<(const jvAlignCache& left,const jvAlignCache& right)
{
	if (left.m_len != right.m_len)
	{
		return (left.m_len<right.m_len);
	}
	for (int i=left.m_len-1;i>=0;i--)
	{
		if (left.m_match[i] != right.m_match[i])
		{
			return left.m_match[i] < right.m_match[i];
		}
	}
	return false;
}



jvAlign::jvAlign(int imaxalen,int imaxrlen, jvAlignCosts acosts):inner_aligner(imaxalen, imaxrlen, acosts)
{
	string _method_name = "jvAlign::jvAlign";        
	maxalen = imaxalen+1; maxrlen = imaxrlen+1;
	aligned = new jvBase[MAXALIGNEDLEN];

//	aligned = (jvBase*)malloc(sizeof(jvBase)*MAXALIGNEDLEN);
//
	edit =  new char[MAXEDITLEN];
	//buffA = (jvBase*)malloc(sizeof(jvBase)*maxalen);
	//
	buffA = new jvBase[maxalen];
}

jvAlign::~jvAlign()
{
        if( aligned != 0) {
           delete [] aligned;
           aligned=0;
        }

        if( edit != 0) {
	   free(edit);edit = 0;
        }
        if(buffA!= 0 ) {
	  free(buffA);buffA=0;
        }
}

jvBase* jvAlign::getAligned(void)
{
	return aligned;
}

int jvAlign::getAlignedLen(void)
{
	return alignedLen;
}

char* jvAlign::getEdit(void)
{
	return edit;
}

int jvAlign::getEditLen(void)
{
	return editLen;
}


double jvAlign::getScore(void)
{
	return score;
}

double jvAlign::align(jvBaseVec *A, jvBaseVec *R, bool print)
{
	return align(A,0,A->size(),R,NULL,print);
}



double jvAlign::getCacheHits(bool reset)
{
    double res = (double)m_cache_hits/(m_cache_misses+m_cache_hits);
    if (reset)
    {
        m_cache_hits = 0;
        m_cache_misses = 0;
    }
    return res;
}
    


// 2-side align
// input: string A, string R, mask M
// process:  align A to R where M == 0 and R to M where M==1
// output: aligned A, aligned R, edit string 1, edit string 2


double jvAlign::align(jvBaseVec *A, int start, int end, jvBaseVec *R, set<jvAlignCache>* cache, bool print)
{
	int st = max(start,0);       int leftPad = max(0,0-start);
	int en = min(end,(int)A->size()); int rightPad = max(0,end-(int)A->size());
	buffAlen = end-start;

	if (st<=en)
	{
		fill(buffA,buffA+leftPad,nBases);
		copy(A->begin()+st,A->begin()+en,buffA+leftPad);
		fill(buffA+leftPad+en-st,buffA+leftPad+en-st+rightPad,nBases);
		assert(leftPad+en-st+rightPad == buffAlen);
	}
	else
		fill(buffA,buffA+buffAlen,nBases);

	//  toggle cache on/off
	if (!with_cache) cache = NULL;

	set<jvAlignCache>::iterator it;
	if (print)
	{
            cout << endl << " st:" << st << " en: " << en << " leftPad: " << leftPad << " right: " << rightPad << " end-start " << end-start  << endl;
            cout << "buffA: ";
            for (int i=0;i<end-start;i++)
            {
                cout << buffA[i] + 1;
            }
            cout << "\nbuffR: ";
            for (uint i=0;i<R->size();i++)
            {
                cout  << ((*R)[i]+1);
            }
            cout << endl;
	}
        
	if (cache)
	{
		it = cache->find(jvAlignCache(buffA,buffAlen));
    }
    
	if (cache && (it != cache->end()))
	{
            assert(it->m_bases.size()<MAXALIGNEDLEN);
            assert(it->m_edit.size()<MAXEDITLEN);
            alignedLen = it->m_bases.size();
            copy(it->m_bases.begin(),it->m_bases.end(),aligned);
            editLen = it->m_edit.size();
            copy(it->m_edit.begin(),it->m_edit.end(),edit);
            score = it->m_score;
            localStart = it->localStart;
            if (iteration - it->last > 5)
            {
                jvAlignCache ut(*it);
                cache->erase(it);
                ut.setLast(iteration);
                cache->insert(ut);
            }
            if (print) cout<<"Cache! ";
            m_cache_hits++;
	}
	else
	{
            alignedLen = MAXALIGNEDLEN;
            editLen = MAXEDITLEN;
            score = alignFast(buffA,  end-start,   &( (*R)[0] ), (int)R->size(),  aligned,  &alignedLen, edit,&editLen);
            if (print) cout<<"Cache miss! ";
            m_cache_misses++;
            if (cache && cache->size()<3 && (it == cache->end()))
            {
                assert(alignedLen<MAXALIGNEDLEN);
                assert(editLen<MAXEDITLEN);
                jvAlignCache ac(buffA,end-start);
                ac.m_bases.resize(alignedLen);
                copy(aligned,aligned+alignedLen,ac.m_bases.begin());
                ac.m_edit.resize(editLen);
                copy(edit,edit+editLen,ac.m_edit.begin());
                ac.m_score = score;
                ac.localStart = localStart;
                ac.last = iteration;
                pair <set<jvAlignCache>::iterator,bool> pr = cache->insert(ac);
                assert(pr.second);
            }
	}

	if (print)
	{
            cout<<"score: "<<score<<endl;
            vector<char> str_edit(edit, edit+editLen);
            print_edit(str_edit);
            cout<<"local Start: "<<localStart<<endl;
	}

	return score;
}



double jvAlign::alignFast(jvBase*A,int lenA,jvBase*R, int lenR,jvBase*R1, int *outLen,char*editR,int*outLenEditR)
{
	assert( editR && outLenEditR && R1 && outLenEditR);
	inner_aligner.align(A, lenA, R, lenR, 9); // -->  THIS IS THE BAND PARAMETER FOR BANDED ALIGNMENT
	char* srcEditR = inner_aligner.getEdit(outLenEditR, true);
	for (int j=0; j<*outLenEditR; j++)
	{
		char newchar;
		switch (srcEditR[j])
		{
			case 0:editR[j]=-1; break;
	 		case 1:editR[j]=4; break;
 			case -1: editR[j] = -2;
		}
	}
	localStart = inner_aligner.getLocalStart();
	inner_aligner.applyEdit(R, lenR, R1, outLen, true);
	return inner_aligner.getScore();
}




