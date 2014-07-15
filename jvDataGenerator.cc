#include "jvDataGenerator.h"
#include "jvAlign.h"
#include <string>
#include <stdlib.h>
#include <assert.h>
#include <iterator>
#include <algorithm>
#include <stdio.h>
#include <float.h>

void complement_bases(jvBaseVec& bases)
{
    reverse(bases.begin(), bases.end());
    for (int i=0; i<(int)bases.size(); i++)
        bases[i] = complement_base(bases[i]);
    return;
}

jvBase complement_base(jvBase val)
{
    if (val < nBases)
        val = nBases -1 - val;
    return val;
}

jvBaseCounts complement_count(jvBaseCounts val)
{
    jvBaseCounts nval(4);
    for (jvBase i=0;i<nBases;i++)
    {
        nval[i] = val[(int)complement_base(i)];
    }
    return nval;
}

void complement_counts(jvBaseCountsVec& counts)
{
    reverse(counts.begin(),counts.end());
    for (int i=0;i<(int)counts.size();i++)
    {
        counts[i] = complement_count(counts[i]);
    }

}



void jvReadData::print()
{
    cout<<"s: "<<m_seqId<<"  o: "<<m_offset<<"  data: ";
    for(uint i=0; i < m_bases.size(); i++) cout << m_bases[i];
    cout<<endl;
}

bool jvReadData::isCorrect()
{
	bool correct = true;
	if (m_bases.size()==m_correct.size())
	{
		for (uint i=0; i<m_bases.size();i++)
		{
			if (m_bases[i] != m_correct[i])
			{
				correct = false;
				break;
			}
		}
	}
	else
	{
		correct = false;
	}
	return correct;

}
bool jvReadData::printBoth()
{
    bool correct = isCorrect();
    cout<<(getReversed() ? "Read is reversed\n" : "Read is original\n");
    cout<<"Raw:     ";
    for(uint i=0; i < getRawBases()->size(); i++) if ( (*getRawBases())[i] == 4) cout << "-"; else cout << (*getRawBases())[i]+1;
    cout<<endl;
    cout<<"Aligned: ";
    for(uint i=0; i < m_bases.size(); i++) if (m_bases[i] == 4) cout << "-"; else cout << m_bases[i]+1;
    cout<<endl;

    cout<<"Correct: ";
    for(uint i=0; i < m_correct.size(); i++) if (m_correct[i] == 4) cout << "-"; else cout << m_correct[i]+1;
    cout<<endl;

    cout<<"Edit:    ";
    print_edit(m_edit);
    
    return correct;
}

void jvDataGenerator::saveSequence(char* name, bool append)
{
        FILE *fout;
	if (!append) fout = fopen(name,"w+");
        else fout = fopen(name,"a+");
        if (!fout)
        {
                cerr << "Could not open: " << name << " for writing\n";
                throw exception();
        }
        char map[] = "ACGTN";
        fprintf(fout, ">\n");
        for (uint l=0;l<m_bases.size();l++)
        {
                fprintf(fout,"%c",map[m_bases[l]]);
        }
        fprintf(fout,"\n");

        fclose(fout);
}

jvDataGenerator::jvDataGenerator(int random_seed)
{
	m_seqId = 0;
	m_siteId = 0;
	srand(random_seed);
	m_noise = 0;
	m_insProb = 0;
	m_delProb = 0;
	m_missing = 0;
	m_revProb = 0;
}


void jvDataGenerator::setSequence(jvBaseVec& input)
{
    m_bases = input;
}


void jvDataGenerator::setSequence(char* strSeq)
{
    int nBases = strlen(strSeq);
    m_bases.resize(nBases);
    for(int i=0; i<nBases; i++)
        m_bases[i] = strSeq[i]-'1';
}

void jvDataGenerator::randomizeBackbone(int nBases, int nDistance)
{
	m_backbone.clear();
    m_backbone.resize(nBases, 4);
    for(int i=0; i<nBases; i++)
    {
		if (i && (i%nDistance == 0)) continue;
		m_backbone[i] = rand() % 4;
    }
}

void jvDataGenerator::randomizeVariant(jvBaseVec mut)
{
	m_bases.clear();
	for (uint i=0;i<m_backbone.size();i++)
	{
		if (m_backbone[i] < 4)
			m_bases.push_back(m_backbone[i]);
		else
			for (int j=0; j<(int)mut.size(); j++)
				m_bases.push_back(mut[j]);
	}
}

void jvDataGenerator::setReadLength(int length, uint var)
{
	m_readLength = length;
	m_var = var;
}


void jvDataGenerator::randomizeSequence(int nBases)
{
    m_bases.resize(nBases);
    for(int i=0; i<nBases; i++)
        m_bases[i] = rand() % 4;
}

void jvDataGenerator::printSequence()
{
	cout << "Sequence:\n";
	for (uint i=0;i<m_bases.size();i++)
	{
		cout << m_bases[i]+1;
	}
	cout << "\n";
}

void jvDataGenerator::addToReads()
{
    int insCutOff = (int)(m_insProb*((double)RAND_MAX));
    int delCutOff = (int)(m_delProb*((double)RAND_MAX));
    int noiseCutOff = (int)(m_noise*((double)RAND_MAX));
    int missCutOff = (int)(m_missing*((double)RAND_MAX));
    int revCutOff = (int)(m_revProb*((double)RAND_MAX));
    cout<<"cutoffs: "<<noiseCutOff<<" "<<insCutOff<<endl;
    bool modified;

    long n_ins=0, n_del=0, n_noise=0, n_miss=0, n_fine=0, n_rev=0;  // counters

    //  First randomize offsets
    vector<int> offsets(1,0);
    for (uint i=0; i<m_bases.size(); )
    {
        int diff = m_interval;;
        if (m_rand_interval == DG_RANDOM)
            diff = rand() % (2*m_interval);
        i+= diff;
        offsets.push_back(i);
    }

    for (int o=0; o<offsets.size(); o++)
    {
        int i = offsets[o];

        for (int k=0; k<m_readsPerSite; k++)
        {
            jvBaseVec read;

            // set length of copied area (estimated length of read)
            int length = m_readLength;
            if (m_var)
                length += (2*(rand() % 2)-1)*(rand() % m_var);
            if (i+length > m_bases.size())
                continue;

            jvBaseVec correct(0);
            modified = false;

            for (int j=0; j<length;)  // for every base
            {
                // insertion
                if (rand()<insCutOff)
                {
                    n_ins++;
                    modified = true;
                    read.push_back(rand()%4);
                    continue;
                }

                // deletion
                if (rand()<delCutOff)
                {
                    n_del++;
                    correct.push_back(nBases);
                    modified = true;
                    j++;
                    continue;
                }

                // mismatches
            	if (rand()<noiseCutOff)
            	{
                    n_noise++;
                    read.push_back(rand() % (int)nBases);
                    correct.push_back(read.back());
                    j++;
                    continue;
            	}


                // missing letters
            	if (rand() < missCutOff)
            	{
                    n_miss++;
                    read.push_back((int)nBases);
                    correct.push_back(read.back());
                    j++;
                    continue;
            	}

                n_fine++;
                read.push_back(m_bases[i+j]);
                correct.push_back(read.back());
                j++;
            }

            bool reversed = false;
            if (rand() < revCutOff)
            {
                n_rev++;
                complement_bases(read);
                reversed = true;
            }
            m_reads.push_back(jvReadData(read, m_seqId, i, m_siteId, reversed, correct));
            m_reads.back().setCorrect(correct);

            if (modified)
                m_reads.back().modified = true;
        }
        m_siteId++;
    }
    m_seqId++;
    cout<<"mis-matches = "<<n_noise<<", w/o noise:  "<<n_fine<<", insertions = "<<n_ins;
    cout<<", deletions:  "<<n_del<<", missing letters: "<<n_miss<<", reversals: "<<n_rev<<endl;
}

int rand_n(int N)
{
    return rand() % N;
}

void jvDataGenerator::shuffleReads()
{
    random_shuffle(m_reads.begin(), m_reads.end());
}


void jvDataGenerator::print()
{
	cout<<"seq_length = "<<m_bases.size()<<" interval = "<<m_interval<<" reads_per_site = "<<m_readsPerSite;
    cout<<" n_seqs = "<<m_seqId<<" n_offsets = "<<m_siteId<<" n_reads = "<<m_reads.size()<<" noise = "<<m_noise<<endl;
}

jvReadData::jvReadData(jvBaseVec bases, int seqId, int offset, int siteId, bool orientation, jvBaseVec correct)
    : m_bases(bases), m_raw_bases(bases), m_rev_bases(bases), m_correct(correct), m_seqId(seqId), m_offset(offset), m_siteId(siteId)
{
    m_edit = vector<char>(bases.size()*2);
    m_alignment_score = -DBL_MAX;
    touched = false;
    modified = false;
    m_isReversed = false;
    complement_bases(m_rev_bases);
    m_alignment_cache.clear();
    m_alignment_cache_rev.clear();    
}

jvReadData::jvReadData(jvBaseVec bases, jvBaseVec rawBases, vector<char> edit, double score)
    : m_bases(bases), m_raw_bases(rawBases), m_rev_bases(rawBases), m_edit(edit), m_correct(0), m_seqId(0), m_offset(0), m_siteId(0)
{
    m_alignment_score = score;
    touched = false;
    modified = false;
    m_isReversed = false;
    complement_bases(m_rev_bases);
    m_alignment_cache.clear();
    m_alignment_cache_rev.clear();    
}


void jvReadData::updateAlignment(bool local, bool printAl, jvAlign* aln)
{
	touched = true;
	double score = aln->getScore();

        int ex = 0;
        if (!local)
            ex = aln->getLocalStart();

        int len = aln->getAlignedLen();
        jvBase *aligned = aln->getAligned();
        m_bases.resize(len+ex);
        fill(m_bases.begin(), m_bases.begin()+ex, nBases);
        memcpy(&m_bases[ex],aligned,sizeof(jvBase)*len);

        for (int i=0;i<m_bases.size();i++) assert(m_bases[i]<=(int)nBases);

        m_alignment_score = score;
        int editLen = aln->getEditLen();
        char* edit = aln->getEdit();        
        m_edit.resize(editLen+ex);
        fill(m_edit.begin(), m_edit.begin()+ex, nBases);
        memcpy(&m_edit[ex],edit,sizeof(char)*editLen);

        if (printAl)
        {
            cout << score << endl;
            printBoth();
        }

}

void jvReadData::resetAlignment(void)
{
	touched = false;
	m_alignment_score = -DBL_MAX;
	m_bases = m_raw_bases;
}

int jvReadData::cleanAlignmentCache(int cutoff)
{            
    set<jvAlignCache>::iterator it2;
    
    if (m_alignment_cache.size() >= 10)
    {
    	for (set<jvAlignCache>::iterator it=m_alignment_cache.begin();it != m_alignment_cache.end();)
    	{
            
            if (it->last < cutoff)
            {
                it2 = it;it++;
                m_alignment_cache.erase(it2);
            }
            else
            {
                it++;
            }
    	}
    }
    
    
    if (m_alignment_cache_rev.size() >= 10)
    {
    	for (set<jvAlignCache>::iterator it=m_alignment_cache_rev.begin();it != m_alignment_cache_rev.end();)
    	{
            
            if (it->last < cutoff)
            {
                it2 = it;it++;
                m_alignment_cache_rev.erase(it2);
            }
            else
            {
                it++;
            }
    	}
    }
    return m_alignment_cache_rev.size()+m_alignment_cache.size();
}

jvReadData jvReadData::complement()
{
    jvBaseVec bases = m_bases;
    complement_bases(bases);
    return jvReadData(bases, 0, 0, false);
}

void jvReadData::complementSelf()
{
    setReversed(!getReversed());
    complement_bases(getBases());
    reverse(getEdit()->begin(), getEdit()->end());
}


void jvReadData::insertAtOffset(int pos)
{
    m_bases.insert(m_bases.begin()+pos,nBases);
    jvAlign::edit_aligner.align(&m_bases, 0, m_bases.size(), getRawBases(), getAlignmentCache(), false);
    updateAlignment(false, false,  &jvAlign::edit_aligner);
}

void jvReadData::deleteAtOffset(int pos)
{
    m_bases.erase(m_bases.begin()+pos);
    jvAlign::edit_aligner.align(&m_bases, 0, m_bases.size(), getRawBases(), getAlignmentCache(), false);
    updateAlignment(false, false, &jvAlign::edit_aligner);
}


long jvReadData::getMem()
{  
  return 0 + 
    m_bases.capacity()*sizeof(jvBase) + //jvBaseVec m_bases;
    m_raw_bases.capacity()*sizeof(jvBase) + //jvBaseVec m_raw_bases;
    m_rev_bases.capacity()*sizeof(jvBase) + //jvBaseVec m_rev_bases;
    sizeof(bool) + //bool m_isReversed;
    m_correct.capacity()*sizeof(jvBase) + //jvBaseVec m_correct;
    m_edit.capacity()*sizeof(char) + //vector<char> m_edit;
    m_alignment_cache.size()*sizeof(jvAlignCache) + //set<jvAlignCache> m_alignment_cache;
    m_alignment_cache_rev.size()*sizeof(jvAlignCache) + //set<jvAlignCache> m_alignment_cache_rev;
    sizeof(double) +//double m_alignment_score;
    sizeof(int) +//int m_seqId;
    sizeof(int) +//int m_offset;
    sizeof(int) +//int m_siteId;
    sizeof(bool) +//bool m_orientation;
    sizeof(bool) +//bool modified;
    sizeof(bool);//bool touched;

}
