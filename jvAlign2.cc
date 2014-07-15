#include "jvAlign2.h"
#include <limits.h>
#include <float.h>
#define FM1(i,j) (*(fm+(i)*(m_lenR+1) + (j)))
#define BM1(i,j) (*(bm+(i)*(m_lenR+1) + (j)))


void inline max3(double *vals,double *mx,unsigned char* argmx)
{
	if (vals[0]>=vals[1] && vals[0]>=vals[2])
	{
		*mx = vals[0];
		*argmx = 0;
	}
	else
	{
		if (vals[1] >= vals[2])
		{
			*mx = vals[1];
			*argmx = 1;
		}
		else
		{
			*mx = vals[2];
			*argmx = 2;
		}
	}

}

jvAlign2::~jvAlign2()
{
        if( bm!=0) { free(bm);bm=0; }
	if( edit != 0 ) { free(edit);edit = 0; }

        if(fm!= 0) { free(fm);fm=0;	}
}




void jvAlign2::init(jvAlignCosts costs)
{
	string _method_name = "jvAlign2::init";

        cout<<"setting costs to: ";
        costs.print();

	dx = costs.indel; dy = costs.indel;
	sx = costs.pad_insert;

        // p[i][j] refers to matching letter i in the *sequence*
        //         with letter j in the *read*
	for (int i=0;i<5;i++)
	{
            for (int j=0;j<5;j++)
            {
                p[i][j] = costs.miss;		
            }
            p[i][i] = costs.hit;
	}
	for (int i=0;i<5;i++)
	{
            p[i][4] = 0; 
            p[4][i] = costs.unobserved_match; 
	}
	p[4][4] = 0;     
}


jvAlign2::jvAlign2(int imaxalen, int imaxrlen, jvAlignCosts costs)
{
	string _method_name = "jvAlign::jvAlign";

	maxalen = imaxalen+1; maxrlen = imaxrlen+1;
	maxeditlen = maxalen + maxrlen; 

	//  prepare shared memory structures

	//  edit strings and aligned strings	
	//  edit string can be compressed!
	edit = (char*)malloc(sizeof(char)*maxeditlen*2);

	//  prepare matrix for forward and backward	
	fm = (double*)malloc(sizeof(double)*maxalen*maxrlen);
	bm = (unsigned char*)malloc(sizeof(unsigned char)*maxalen*maxrlen);


	if (!fm || !bm || !edit) 
	{ 
		cerr << _method_name << " : out of memory\n"; throw exception();
	}
	
	m_scoreOnly = true;
	si[0] = -1; si[1] =  0; si[2] = -1;
	sj[0] =  0; sj[1] = -1; sj[2] = -1;
			
	//  prepare matrix for easy access to penalty values	
	init(costs);
	 	
}

int jvAlign2::align(jvBase* A, int lenA,jvBase* R, int lenR, int band)//, bool scoreOnly)
{
	if (abs(lenA-lenR) > band)
	{
		band = (abs(lenA-lenR));
//		cout<<"Warning - band too short - increasing band\n";
	}
	
	m_lenA = lenA;  m_lenR = lenR;	
	m_A = A; m_R = R; m_band = band;
		
	alignFast();
	int misses = backtraceFast();
	
	return misses;
}



//  Aligns A ("the sequence") to R ("the read").  
//  meaning of FM(i,j) is:  Best score after aligning 1..i-1 to 1..j-1 (i and j were not handled!)
//   
double jvAlign2::alignFast()
{
    
    FM1(0,0) = 0;		
    string _method_name = "jvAlign::alignFast";	
    if (m_lenA>=maxalen || m_lenR>=maxrlen)
    {
        cerr << _method_name << " sequences too long\n";
        throw exception();
    }
    
    double vals[3];
    
    for (int i=0; i<=m_lenA; i++) // i goes over sequence
    {
        int left_margin = max(0, i-m_band);
        int right_margin = min(m_lenR, i+m_band);
        
        for (int j=left_margin; j<=right_margin; j++) // j goes over read		
        {
            if (i==0 && j==0)  // remove this if later
            {
                continue;
            }
            
            vals[0] = -DBL_MAX; vals[1] = -DBL_MAX; vals[2] = -DBL_MAX;
            
            // coming from above - deletion in the sequence/insertion in the read
            // ==> letter i-1 in the sequence is skipped.      
            if ( (i>0) && (j < i+m_band) )
            {
                vals[0] = FM1(i-1,j);
                if (j>0 && j<m_lenR)  // penalize insertions only if they are *inside* the read
                    vals[0] = vals[0] + dx;// + q[A[i-1]];  // insertion in the read
                else
                    vals[0] = vals[0] + sx;// + q[A[i-1]];  // insertion outside  the read								
            }
            
            // coming from the left - insertion in the sequence/deletion in the read
            // ==> letter j-1 in the read is skipped.      
            if (j>left_margin) //( (j>0) && (j > i-K) ) 
            {
                vals[1] =  dy + FM1(i,j-1);	// deletion in read			
            }
            
            if (i>0 && j>0)
            {				
                vals[2] = FM1(i-1,j-1) + p[m_A[i-1]][m_R[j-1]]; // match/mismatch
            }
            
            max3(vals,&FM1(i,j),&BM1(i,j)); // max and updates the matrices				
        }
    }
    score =  FM1(m_lenA, m_lenR);
    m_scoreOnly = true;
    return score;
}

int jvAlign2::backtraceFast()
{
		
	string _method_name = "jvAlign::backtraceFast";
	int misses = 0;  int indels = 0;
	
	if ((maxeditlen < m_lenA + m_lenR))
	{
		cerr << _method_name << " output string buffers too short \n";
		throw exception();
	}
	
	unsigned char id;
	int i=m_lenA; int j=m_lenR;
	id = BM1(i,j);
	int ei = m_lenA + m_lenR;
	local_start = 0;
	
	while (1)
	{				
		edit[--ei] = sj[id]-si[id];

		i = i+si[id];
		j = j+sj[id];

		if (i<0 || j<0)
		{
			cerr << _method_name << " negative indexes" << endl;			
			throw exception();		
		}

		id = BM1(i,j);

		if (edit[ei] != 1)
			local_start = ei;
				
		if ((edit[ei] == 0) && (m_A[i] != nBases) && (m_R[j] != nBases) && (m_A[i] != m_R[j]))
			misses++;
		else
			indels++;

		if (i==0 && j==0)
			break;
		
			
	}	

	//  shift back edit string
	editLen = m_lenA+m_lenR-ei;
	memmove(edit,edit+ei, sizeof(char)*(editLen));
	m_scoreOnly = false;
	local_start -= ei;
	return misses;
}	
	
char* jvAlign2::getEdit(int* outEditLen, bool local)
{
	assert(!m_scoreOnly);
	*outEditLen = editLen;			

	if (!local)
            return edit;
	
	char* edit_trimmed = edit;
	while (*edit_trimmed == 1)
	{
            edit_trimmed++;
            (*outEditLen)--;
	}
	assert(*outEditLen >0);
	assert (edit + local_start == edit_trimmed); 
	while (*(edit_trimmed + *outEditLen -1) == 1)
	{
            (*outEditLen)--;		
	}
	assert(*outEditLen >0);
	return edit_trimmed;

}

void jvAlign2::applyEdit(jvBase* inStr, int inLen, jvBase* outStr, int *outStrLen, bool local)
{
	string _method_name = "jvAlign::applyEdit";	
	assert(!m_scoreOnly);
	int in=0, out=0;
	for (int i=0; i<editLen; i++)
	{
		switch (edit[i])
		{
			case 0:						
				outStr[out++] = ( in < inLen ? inStr[in++] : nBases); break;
			case -1:
				in++; break;
			case 1:
				if (!(local && (in==0 || in == inLen)))
					outStr[out++] = nBases;
		}					
	}
	if (out >= *outStrLen)
	{
		cerr << _method_name << " output string too short\n";
		throw exception();
	}
	*outStrLen = out;
	
	return;
}

