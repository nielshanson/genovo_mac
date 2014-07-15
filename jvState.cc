#include <algorithm>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <iomanip> // for setbase(4)
#include <fstream>
#include <sstream>

#include "jvState.h"
#include "jvSampler.h"
#include "jvAlign.h"

double jvSequence::pseudoA=1;
double jvSequence::pseudoB=1e2;


uint jvInstance::WIDTH = 1;
uint jvInstance::DEFAULT_READ_LENGTH;
uint jvState::K;
int jvState::endAl = 300;
bool jvState::alignment_only = false;
double jvState::garbage_thresh = log(0.25);
bool jvState::construct = false;
bool jvState::both_orientations = false;
bool jvState::all_reads = false;
bool jvState::with_temperature = false;
bool jvState::with_geometric = true;
double jvState::temperature = 1;
bool constantP = false;

int iteration;
int out_of_bounds = 0;
bool debugOn;
bool walign = true;
int change_seq = 0;

jvErrorModel errorModel;

void jvState::dumpHashCounts()
{
  buildHash();
  FILE *f = fopen("dump.m_hash","w+");
  map<jvSequence*,int> reverse; // reverse index not reverse sequence

  for (int i=0;i<m_seq.size();i++)
    {
         reverse.insert(pair<jvSequence*,int>(m_seq[i],i));
 
    }
  for (int i=0;i<m_hash.size();i++)
    {
      if (m_hash[i].size())
	{
	  fprintf(f,"%d ",i);

	  for (set<jvLocation>::iterator it=m_hash[i].begin(); it!=m_hash[i].end(); it++)
	    {
	      fprintf(f,"(%d,%d) ",reverse.find(it->first)->second,it->second);
	      
	    }
	  fprintf(f,"\n");
	}      
    }
  fclose(f);
  reverse.clear();
}

void jvState::setAlignmentOnly(bool _construct, double _garbage_thresh)
{
    cout<<"Setting alignment only.  Geometric and Dirichlet are off.\n";
    cout<<"Threshold "<<_garbage_thresh<<", alignment will "<<(_construct? "be dumped to mapping.dump" : "not be dumped")<<".\n";
    alignment_only = true;
    garbage_thresh = _garbage_thresh;
    construct = _construct;
    setWithInstances(false);
    setGeometric(false);
    setBatch(true);
    jvSampler::setMAP(true);
}



void prt(jvBaseVec buff)
{
	    for (int k=0; k<(int)buff.size(); k++)
	    	cout<<(int)buff[k]<<" ";
	    cout<<endl;
}


void print_bases(jvBaseVec* vec)
{
    for (int i=0; i<vec->size(); i++)
    {
        assert((*vec)[i] <= (int)nBases);
        if ((*vec)[i] == (int)nBases)
            cout<<".";
        else
            cout<<(int)(*vec)[i]+1;
    }
    cout<<endl;
}


void set_edit(string str, vector<char>& edit)
{
    char map[26] = {0};
    map['D'-'A'] = -2;
    map['N'-'A'] = -1;
    map['I'-'A'] = 4;
    
    edit.clear();
    stringstream ss(str);
    char action;
    int count;
    char check = '-';
    do
    {
        assert(check == '-');
        ss >> count>>action;
        if (ss.fail( ))
            throw (str +" has a problem\n");
        assert(count > 0);
        edit.insert ( edit.end(), count, map[action-'A'] );            
    }
    while(ss>>check);
}


string get_edit(vector<char>& edit)
{
    stringstream ss;
    
    char map[8] = {'D', 'N', 'A', 'C', 'G', 'T', 'I', 'E'};
    char curr_ch = edit[0];
    int count = 0;
    for (int i=0; i<edit.size(); i++)
    {
        if (edit[i] != curr_ch)
        {
            
            ss<<count<<map[curr_ch+2]<<'-';
            count = 0;
            curr_ch = edit[i];
        }
        count++;
    }
    if (count>0)
        ss<<count<<map[curr_ch+2];

    return ss.str();
}


void print_edit(vector<char>& edit)
{
    if (edit.size() == 0)
    {
        cout<<"Empty.\n";
        return;
    }
    cout<<get_edit(edit)<<endl;
}


void print_hist(vector<int>& hist)
{
    int sum_vec = sum(hist);
    for (int i=0; i<hist.size(); i++)
        cout<<100.0*hist[i]/sum_vec<<" ";
    cout<<endl;
}

//  State methods

void jvState::setAuto()
{
    assert(m_reads.size()>0);
    int sum_reads = 0;
    int N = min((int)m_reads.size(), 10000);
    for (int i=0; i<N ; i++)
    {
        sum_reads += m_reads[i].getData()->getRawBases()->size();
    }
    int avg_reads = sum_reads/N + 1;
    cout<<"Average read length: "<<avg_reads<<" bases.\n";

    setEndAl((int)(2.9*avg_reads));
    setBothOrientations(true);
    setAllReads(true);
    setAnnealing(false);
    setErrorModel(0.01,0.01);
    m_start = clock();
}


void jvState::setErrorModel(double noise, double indels)
{
    errorModel.setParams(noise, indels);
    jvAlignCosts costs     (errorModel.logIndel(), errorModel.logNoise(), errorModel.logInvNoise(), logBasePrior, 0);
    jvAlignCosts mergecosts(errorModel.logIndel(), errorModel.logNoise(), errorModel.logInvNoise(), logBasePrior, 0);
    jvAlignCosts editcosts (errorModel.logIndel(), errorModel.logNoise(), errorModel.logInvNoise(), logBasePrior, errorModel.logIndel());
    
    jvAlign::aligner.setCosts(costs);
    jvAlign::merge_aligner.setCosts(mergecosts);
    jvAlign::edit_aligner.setCosts(editcosts);
}


void jvState::print(PRINT_OPTIONS opt)
{
    assert(m_seq.back()->getM() == 0);
    if (opt & PRINT_BASES)
        for (uint i=0; i < m_seq.size()-1; i++)
	{
	    if ( (opt & PRINT_SKELETON) && !m_seq[i]->m_never_clean) continue;
            m_seq[i]->print((PRINT_OPTIONS)(opt & (PRINT_BASES | PRINT_SUMMARY)));
	}
    else
    {
        for (int i=0; i<(int)m_seq.size(); i++)
            if (m_seq[i]->getRight()-m_seq[i]->getLeft() > 400)
                cout<<m_seq[i]->getRight()-m_seq[i]->getLeft()<<" ";
        cout<<endl;
    }

    if (opt & PRINT_INS)
        for (uint i=0; i < m_seq.size()-1; i++)
	{
	    if ( (opt & PRINT_SKELETON) && !m_seq[i]->m_never_clean) continue;
            m_seq[i]->print(PRINT_INS);
	}

    cout<<"gamma = "<<m_gamma<<"  alpha = "<<m_alpha<<endl;
    if (opt & PRINT_LL)
        ll().print();

    if (opt & PRINT_DEBUG)
        for (uint i=0; i < m_seq.size()-1; i++)
            m_seq[i]->print(PRINT_DEBUG);

	if (opt & PRINT_READS)
	{
    	cout<<"\nReads:\n";
    	for (uint i=0; i < m_data.size(); i++)
        	m_data[i].print();
	}
    cout<<endl;
}

jvState::jvState(vector<jvReadData>& data, bool empty, int seed):m_data(data), m_reads()
{
        candI = -1;  candJ = -1;
        seq_choice = -1; cutoff = -1;
	if (!jvSampler::isActive())
		jvSampler::init(seed);

        endAl = 100; // max alignment length for align ends
	//  Allocating space for hash table
	K = 10;  //  Length of a hash entry
	m_hash.resize( (1<<2*K));


	//  Instance Parameters - length of a read
	jvInstance::DEFAULT_READ_LENGTH = m_data[0].getBases().size();

	//  Concentration parameters
	m_alpha = 1;
	m_gamma = 1;
	m_withInstances = false;
	m_batch = false;

    if (!empty)
	{
		for (uint i=0; i<data.size(); i++)
            m_reads.push_back(jvReadAss(&data[i]));

        m_seq.push_back(new jvSequence());
	}
        m_start = clock();
        m_total_splits = 0;
        m_total_merges = 0;
        m_total_failures = 0;
}

jvInstance* jvState::newIns(jvSequence* p_seq, int offset)
{
    jvInstance* pIns = new jvInstance(p_seq, offset);
    m_ins.push_back(pIns);
    return pIns;    
}

jvInstance* jvState::newIns(int seqId, int offset)
{
    jvSequence* p_seq = (seqId <0 ? NULL: m_seq[seqId]);
    return newIns(p_seq, offset);
}

jvState::~jvState()
{
    for (uint i=0; i<m_ins.size(); i++)
        delete m_ins[i];

    for (uint i=0; i<m_seq.size(); i++)
        delete m_seq[i];
}



jvSequence* jvState::newSeq(jvSequence* seq)
{
    if (seq)
    {
        delete m_seq.back();
        m_seq.back() = seq;
    }

    m_seq.push_back(new jvSequence());

    return m_seq.back();
}


void jvState::cleanInstances()
{
    for (int i=(int)m_ins.size()-1; i>=0; i--)
        if (m_ins[i]->getN() == 0)
        {
            m_ins[i]->assign(false);
            delete m_ins[i];
            m_ins.erase(m_ins.begin()+i);
        }

}

int jvState::cleanSequences()
{
    int ct = 0;
    assert(m_seq.back()->getM() == 0);
    for (int i=(int)m_seq.size()-2; i>=0; i--)
        if (!m_seq[i]->m_never_clean && (m_seq[i]->getM() == 0))
        {
            delete m_seq[i];
            m_seq.erase(m_seq.begin()+i);
            ct++;
        }
    assert(m_seq.back()->getM() == 0);
    return ct;
}

double jvState::cleanAlignmentCache(int cutoff)
{
    int sum_cache = 0;
    for (uint i=0;i<m_reads.size();i++)
    {
        sum_cache += m_reads[i].getData()->cleanAlignmentCache(cutoff);
    }
    double avg = (double)sum_cache/m_reads.size();
    cout<<"On average, cache size is "<<avg<<" per read.\n";
    return avg;
}

// greedy override - runs sampleBases in greedy mode regardless of the jvSampler::isMAP value.
void jvState::sampleBases(bool greedy_override)
{
    for (uint i=0; i<m_seq.size(); i++)
        m_seq[i]->sampleBases(greedy_override);
}


bool comp_seq(jvInstance* a, jvInstance* b)
{
    return (*a < *b);
}



void jvState::freezeInstances(bool freeze)
{
    for (uint i=0; i<m_seq.size(); i++)
        m_seq[i]->setFreeze(freeze);

}



jvLikelihood jvState::ll()
{
    jvLikelihood l;
    jvDVec M = getM();
    jvDVec N = getN();

    for (uint i=0; i<m_seq.size(); i++)
    {
        jvLikelihood l_other = m_seq[i]->ll();
        l = l + l_other;
    }

    l.s = jvSampler::DPll(M, m_gamma);
    l.t = jvSampler::DPll(N, m_alpha);

    return l;
}

vector<int> hist_hits_per_spike;

jvLogProb jvState::gibbsReads(vector<jvReadAss*> reads, GIBBS_POLICY gibbs_policy)
{
    assert( ! (gibbs_policy & GIBBS_BOTH_ORIENTATIONS) || both_orientations );
    assert( ! (gibbs_policy & GIBBS_INCLUDE & GIBBS_EXCLUDE) );
    freezeInstances(false);

    jvSpikePolicy spike_policy;
    if (gibbs_policy & GIBBS_NEW_SEQS_ONLY) //  force the first x reads into a singleton
    {
        spike_policy.code = SPIKE_NEW_SEQ;
        cout<<"spike only in new seq\n";
    }
    else if (gibbs_policy & GIBBS_SNM_ONLY)
    {
        spike_policy.code = SPIKE_SNM;
        cout<<"SNM - two spikes only\n";
    }
    else
    {
        spike_policy.code = SPIKE_DEFAULT;
        cout<<"Regular spikes.\n";
    }         
     

// construct is relevant only if m_batch is true.
    // construct = false:  don't map the reads, just compute likelihood.
    // construct = true:   map the reads and also compute the likelihood.
    //                     in order for the likelihood calculation to be consistent and not effected by the changes
    //                     to the sequences due to read assignments, the reads are assigned to instances, but the instances are assigned only at the end.
    //                     garbage reads are not assigned.
//    bool construct = false;  // REMARK:  Construct is now a global paramtere that can be set with a function

    buildHash();

    if (with_temperature) setTemperature();

    vector<int> hist_spikes(20,0);
    hist_hits_per_spike.clear();
    hist_hits_per_spike.resize(150,0);
    
    m_sum_spikes = 0;
    m_count_reads = 0;
    change_seq = 0;  // global variable
    int change_loc = 0;
    int nsingle = 0;
    jvAlign::aligner.getCacheHits(true);

    //  which set of reads are we working on
    int N = (gibbs_policy & GIBBS_INCLUDE ) ? reads.size() : m_reads.size();

    //  what is the read order ? if we move them to new seqs, keep fixed order.
    vector<int> perm = jvSampler::randomPermutation(N, gibbs_policy & GIBBS_NEW_SEQS_ONLY);

    if (gibbs_policy & GIBBS_EXCLUDE)
    	for (int i=0; i<reads.size(); i++)
       	    reads[i]->m_skip = true;
    
    uint j=0;
    jvLogProb whole_transition(0,0);

    for (vector<int>::iterator i=perm.begin(); i<perm.end(); i++)
    {
        bool debug = false;
        jvReadAss& read = ( (gibbs_policy & GIBBS_INCLUDE) ? (*reads[*i]) : m_reads[*i] );
        if ( (gibbs_policy & GIBBS_EXCLUDE) && read.m_skip)
        {
            read.m_skip = false;
            continue;
        }

        // check the current location of the read. count how many problematic spots it covers.
        // If there are less than X such spots - don't resample the read.
        if ( (gibbs_policy & GIBBS_FREEZE) &&(read.isAssigned()) )
        {
            if (read.m_oneSpike) continue;

            jvSequence* seq = read.getIns()->getSeq();
            assert(seq);
            vector<int>& v = seq->m_problematic;
            if (v.size() == 0)
                continue;
            
            int start = read.getIns()->getOffset() + read.getOffset();
            int end = start+read.getData()->getBases().size();
            vector<int>::iterator low,up;
            low = lower_bound (v.begin(), v.end(), start);
            up = lower_bound (v.begin(), v.end(), end);
            int bad = (int)(up-low);
            if (bad == 0) continue;
        }
                
        if (! ((++j) %1000)) cout<<j<<" "<<flush;

        read.assign(false);
        
        //  remove empty instances from sequences - they will not get any reads in the future.
        jvSequence* old_seq = NULL;
        int old_bin = 0;
        int old_o = 0;
        bool singelton = false;
        if (read.getIns())
        {
            old_seq = read.getIns()->getSeq();
            old_bin = read.getIns()->getSeq()->getBin();
            old_o = read.getIns()->getOffset();
            read.getIns()->assign(false);
            singelton = (old_seq->getM() == 0);
        }
                
        m_nSpikes = 0;
        spike_policy.bin = (singelton ? 0 : old_bin);

        //  hack to clamp the results of GIBBS READS according to split-and-merge requirement
        seq_choice = -1;
        if (cutoff > -1) seq_choice = ( ((*i)<cutoff) ? 0 : 1);
//        if (seq_choice == 1) cout<<(*i)<<" "<<flush;

        //  m_batch = false:  assign the reads incrementally
        //  m_batch = true:  first map *all* the reads, without assigning them.
        if ( (gibbs_policy & GIBBS_BOTH_ORIENTATIONS) || both_orientations)
        {
            //  try reverse direction
            bool curr_ort = read.getData()->getReversed();
            jvLogProb transition;
            jvDVec rev(2);
            read.getData()->setReversed(!curr_ort);
            transition = sampleInstance(NULL, &read, spike_policy, debug);
            rev[1] = transition.partition;
            // undo sample instance...  (it could potentially create new instance/sequence)
            read.getIns()->assign(false);
            read.setIns(NULL);

            //  try current direction
            read.getData()->setReversed(curr_ort);
            transition = sampleInstance(NULL, &read, spike_policy, debug);
            rev[0] = transition.partition;
            
            read.m_oneSpike = (m_nSpikes <= 1); // annotate reads with one spike
            hist_spikes[min((int)hist_spikes.size()-1, m_nSpikes)]++; //update histogram

            // choose best direction.
            int orientation = jvSampler::lmnrnd(rev);
            if (orientation == 1)
            {
                //  chose reverse.   undo last sampleInstance (original direction)
                read.getIns()->assign(false);
                read.setIns(NULL);
                //  now do sampleInstance again in the reverse direction.
                read.getData()->setReversed(!curr_ort);
                transition = sampleInstance(NULL, &read, spike_policy, debug);
            }
            else
                if (debug) cout<<"chose current\n";
            transition.partition = log_sum_exp(rev);
            whole_transition = whole_transition+transition;
        }
        else
        {
            jvLogProb transition = sampleInstance(NULL, &read, spike_policy, debug);
            whole_transition = whole_transition+transition;
            read.m_oneSpike = (m_nSpikes <= 1); // annotate reads with one spike
            hist_spikes[min((int)hist_spikes.size()-1, m_nSpikes)]++; //update histogram
        }


        jvSequence* new_seq = read.getIns()->getSeq();
        if (old_seq != new_seq)
        {
            change_seq++;
            change_loc++;
        }
        else if (read.getIns()->getOffset() != old_o)
            change_loc++;
        else //same location
            read.m_oneSpike = true;  // freeze
            
        if (singelton && (new_seq->getM() == 1))
        {
            nsingle++;
            read.m_oneSpike = true; // annotate singleton as if they have one spike
        }


	if (!m_batch)
            read.assign(true);
        else
        {
            if (construct)
            {
                read.getIns()->assign(false);

                if (read.getIns()->getSeq()->isSkeleton() == false)
                    read.setIns(NULL);
                else
                    read.assign(true);
            }
            else
                if (read.getIns()->getSeq()->isSkeleton() == false)
                    read.setIns(NULL);

            cleanInstances();
            cleanSequences();
        }
    }


    if (construct)
    {
        for (int i=0; i<m_ins.size(); i++)
        {
            if (!m_ins[i]->isAssigned()) m_ins[i]->assign(true);
        }
    }

    seq_choice = -1; 
    cleanInstances();
    int seq_clean = cleanSequences();
    
    cout<<"Cleaned "<<seq_clean<<" sequences.\n";
    cout<<endl<<"Total "<<j<<" moving reads (out of "<<m_reads.size()<<").  On average, there are "<< (double)m_sum_spikes/m_count_reads <<" spikes per read, and ";
    cout<<100 * jvAlign::aligner.getCacheHits()<<"% cache hits.\n";
    cout<<"Histogram for spikes (total = "<<sum(hist_spikes)<<"): ";
    print_hist(hist_spikes);
    cout<<endl<<"Total "<<change_loc<<" reads changed location.  There are "<<nsingle<<" singletons.  Temperature is "<<temperature<<".\n";
    cout<<"There are "<< m_seq.size()-1 <<" sequences.  They are (#reads): ";
    for (int k=0; k<min((int)m_seq.size()-3, 50); k++) cout<<m_seq[k]->getRight()-m_seq[k]->getLeft()<<"("<<m_seq[k]->getM()<<") ";
    cout<<" ... ";
    for (int k=max(0,(int)m_seq.size()-3); k<(int)m_seq.size()-1; k++) cout<<m_seq[k]->getRight()-m_seq[k]->getLeft()<<"("<<m_seq[k]->getM()<<") ";
    cout<<endl;

    return whole_transition;
}



void jvState::gibbsInstances()
{
//    freezeReads(true);
    freezeInstances(false);

    buildHash();

    vector<int> perm = jvSampler::randomPermutation(m_ins.size());
    uint j=0;

    for (vector<int>::iterator i=perm.begin(); i<perm.end(); i++, j++)
    {
        if (! ((j+1) %1000)) cout<<j+1<<" "<<flush;
        m_ins[*i]->assign(false);
        sampleInstance(m_ins[*i], NULL);
        m_ins[*i]->assign(true);
    }
    cout << endl;
    cleanSequences();

}

void jvState::sampleGamma()
{
    m_gamma = jvSampler::sampleConcentration(getM(), m_gamma);
    cout<<"New gamma: "<<m_gamma<<endl;
}

void jvState::sampleAlpha()
{
    m_alpha = jvSampler::sampleConcentration(getN(), m_alpha);
    cout<<"New alpha: "<<m_alpha<<endl;
}

void jvState::center()
{
    freezeInstances(true);
    for (uint i=0; i<m_seq.size(); i++)
        m_seq[i]->center();
}


int broken_instance = 0;

// Find if there is a gap in the sequence, and if yes, returns the first and last instances.
// If there is no gap return (NULL,NULL).
pair<jvInstance*, jvInstance*> jvSequence::getGapInstances(bool all)
{
    assert(false); //added by joni
    
    int GAP_THRESH = 3;
    int gap_length = -INT_MAX;
    pair<jvInstance*, jvInstance*> res(NULL, NULL);
    if (m_ins.size() == 0) return res;

    if (all && (m_ins.size() > 0))
    {
        res.first = m_ins.front();
        res.second = m_ins.back();
        return res;
    }


    for (int o=m_left; o<m_right; o++)
    {
        if (m_bases[pos(o)] != (int)nBases)
        {
            if (gap_length > GAP_THRESH)
            {
                //  found a gap.  The letter o-1 is where the gap ends.  The letter
                // o-gap_length is where the gap starts.
                int t = 0;
                for (; t<(int)m_ins.size(); t++)
                    if (m_ins[t]->getOffset() > o-gap_length)  break;
                if ( t == m_ins.size() )
                    t = m_ins.size() -1;

                assert( (t>=0) && (t<m_ins.size()));
                res.second = m_ins[t];
                res.first = m_ins.front();
                while (t > 0)
                {
                    int u = o-gap_length - m_ins[t-1]->getOffset();
                    int v = o-1          - m_ins[t-1]->getOffset();
                    jvReadAss* r = *(m_ins[t-1]->m_reads.begin());
                    int len_r = r->getData()->getBases().size();

                    if ( (u-1 >=0) && (u-1 <len_r))
                    {
                        if (r->getData()->getBases()[u-1] == nBases)
                        {
                            res.second = m_ins[t-1];
                            t = t-1;
                            continue;
                        }
                        else
                        {
                            if ((v+1<len_r) && (r->getData()->getBases()[v+1] == nBases))
                                res.first = m_ins[t-1];
                            else
                            {
                              cout<<"BAD INSTANCE to choose at offset "<< m_ins[t-1]->getOffset()<<endl;
                              broken_instance++;
                            }
                            break;
                        }
                    }
                    break;
                }
                cout<<res.first->getOffset()<<" "<<res.second->getOffset()<<endl;

                return res;
            }
            gap_length = 0;
        }
        else
            gap_length++;
    }
    return res;
}


bool jvState::isThereAGap()
{
    for (int i=0; i<m_seq.size(); i++)
    {
        pair<jvInstance*, jvInstance*> res = m_seq[i]->getGapInstances();
        if (res.first != res.second)
        {
            cout<<"###  Found A GAP! (seq "<<i<<")\n";
            return true;
        }
    }
    return false;
}



jvSnmUndo jvSequence::split(int offset, bool rightMoves, jvSequence* new_seq)
{
    assert(m_isInstanceFrozen);
    assert(new_seq->getM() == 0);

    jvBaseVec tempVec;

    // get the first instance at this offset (or after).
    uint i = 0;
    if (rightMoves)
    {
        for (; i<m_ins.size(); i++)
            if (m_ins[i]->getOffset() >= offset) break;
        offset = m_ins[i]->getOffset();

        // save bases
        tempVec.assign(&m_bases[pos(offset)], &m_bases[pos(getRight())]);

        // move instances
        new_seq->m_ins.insert(new_seq->m_ins.begin(), m_ins.begin()+i, m_ins.end());
        m_ins.erase(m_ins.begin()+i, m_ins.end());

    }
    else  // left subseq moves
    {
        for (; i<m_ins.size(); i++)
            if (m_ins[i]->getOffset() > offset) break;
        offset = m_ins[i-1]->getOffset();

        // save bases
        tempVec.assign(&m_bases[pos(m_ins.front()->getOffset())], &m_bases[pos(getRight())]);

        // move instances
        new_seq->m_ins.insert(new_seq->m_ins.begin(), m_ins.begin(), m_ins.begin()+i);
        m_ins.erase(m_ins.begin(), m_ins.begin()+i);
    }

    new_seq->m_isInstanceFrozen = false;
    m_isInstanceFrozen = false;
    for(uint j=0; j<new_seq->m_ins.size(); j++)
    {
        new_seq->m_ins[j]->assign(false);
        new_seq->m_ins[j]->m_seq = new_seq;
        new_seq->m_ins[j]->assign(true);
    }
    m_isInstanceFrozen = true;
    new_seq->m_isInstanceFrozen = true;

    int start = new_seq->getLeft() - new_seq->m_ins.front()->getOffset();
    assert(start >= 0);

    copy(tempVec.begin()+start, tempVec.begin()+start+new_seq->getRight()-new_seq->getLeft(),
    		&new_seq->m_bases[new_seq->pos(new_seq->getLeft())]);



    jvSnmUndo undo;
    undo.isSplit = true;
    undo.offset = offset;
    undo.s = this;
    undo.s2 = new_seq;

    return undo;
}


// merge seq into this sequence.  if offset on the right of this sequence,
// put the first instance of seq at offset.
// If offset to the left of this sequence, put there the last instance of seq.
jvSnmUndo jvSequence::merge(jvSequence* seq, int offset)
{
    assert(m_isInstanceFrozen);
    assert(seq != this);
    assert((uint)seq->m_M == seq->m_ins.size());


    jvSnmUndo undo;
    int shift_offset, ov_left, ov_right;

    bool mergeToRight = getOverlap(seq, offset, shift_offset, ov_left, ov_right);
    vector<jvInstance*>::iterator it;

    jvBaseVec tempVec2(seq->m_bases);
    int tempShift2 = shift_offset+seq->m_left;
    jvBaseVec tempVec1(m_bases);
    int tempShift1 = m_left;

    if (mergeToRight)
    {
        it = m_ins.end();
    }
    else
    {
        it = m_ins.begin();
    }


    // back-up overlap bases
    if (ov_left < ov_right)
    	undo.bases.insert(undo.bases.begin(), &seq->m_bases[seq->pos(ov_left)],
    			&seq->m_bases[seq->pos(ov_right)]);
    undo.ov_left = ov_left; undo.ov_right = ov_right;
	//

    seq->shift(shift_offset);



    m_ins.insert(it, seq->m_ins.begin(), seq->m_ins.end());

    m_isInstanceFrozen = false;
    seq->m_isInstanceFrozen = false;
    for(uint j=0; j<seq->m_ins.size(); j++)
    {
        seq->m_ins[j]->assign(false);
        seq->m_ins[j]->m_seq = this;
        seq->m_ins[j]->assign(true);
    }
    seq->m_ins.clear();
    assert((uint)seq->m_M == seq->m_ins.size());
    m_isInstanceFrozen = true;
    seq->m_isInstanceFrozen = true;

    assert(m_bases[pos(getRight()-1)] != (int)nBases);
    for (uint i=0; i<tempVec2.size(); i++)
    {
    	if (tempVec2[i] != (int)nBases)
    	{
    		// if the original m_bases was empty
        	int j= tempShift2+i-tempShift1;
        	if ((j<0) || (j>=(int)tempVec1.size()) || (tempVec1[j] == (int)nBases))
        		m_bases[pos(tempShift2 + i)] = tempVec2[i];

    	}

    }



    assert(m_bases[pos(getRight()-1)] != (int)nBases);

    undo.offset = offset;
    undo.isSplit = false;
    undo.back2 = -shift_offset;
    undo.s = this;
    undo.s2 = seq;
    undo.rightMoves = mergeToRight;
    return undo;
}


void jvSequence::merge_simple(jvSequence* seq, int offset)
{
    assert(!seq->m_isInstanceFrozen);
    assert(!m_isInstanceFrozen);
    assert(seq != this);
    vector<jvInstance*> instancesToMove;

    int left = seq->getLeft();

    for (uint i=0;i<seq->m_ins2.size();i++)
        if (seq->m_ins2[i].size())
            for (set<jvInstance*>::iterator insit = seq->m_ins2[i].begin(); insit != seq->m_ins2[i].end();insit++)
            {
                instancesToMove.push_back(*insit);
            }

    
    for (int i=0; i<(int)instancesToMove.size(); i++)
    {
        instancesToMove[i]->assign(false);
        instancesToMove[i]->setSeq(this, offset + instancesToMove[i]->getOffset()-left);
        instancesToMove[i]->assign(true);
    }

}


int jvSequence::pruneSpikes(jvSpikePolicy policy)
{
  if (m_spikes.size() == 0) return 0;

  int max_gap = 10; 
  int min_hits_per_spikes = (policy.code & SPIKE_SPARSE ? 3 : 6);
  
  vector<int> pruned;
  sort(m_spikes.begin(), m_spikes.end());
  
  m_spikes.push_back(MAX_INT);
  int last = m_spikes[0];
  int hits_per_spike = 0;
  for (uint i=0;i<m_spikes.size();i++)
  {
      assert(m_spikes[i] - last >= 0);
      if (m_spikes[i] - last > max_gap)
      {
          if (hits_per_spike >= min_hits_per_spikes)
              pruned.push_back(last);
          last = m_spikes[i];
          hist_hits_per_spike[min((int)hist_hits_per_spike.size()-1, hits_per_spike)]++;
          hits_per_spike = 1;
      }
      else
      {
	  hits_per_spike++;
      }
  }
  m_spikes = pruned;
  return m_spikes.size();
}

jvLogProb jvSequence::sampleO(jvInstance* ins, jvBaseVec* read, jvParticle* particle, bool debug_print)
{

    assert(read != NULL);

    if (debug_print)
    	cout<<"there are "<<m_spikes.size()<<" Spikes.\n";

    if ( m_spikes.empty())
    {
        m_logPy.clear();
        particle->o_i = 0;
        particle->o_t = 0;
        particle->o_spike = jvAlign::extend; //
        return -1e100;
    }

    uint L = read->size();
    m_logPy.reserve(m_spikes.size()+1);
    m_logPy.clear();
    m_spikes_bu = m_spikes;
    assert(!m_spikes.empty());

    for(vector<int>::iterator o=m_spikes.begin(); o!=m_spikes.end(); o++)
    {
        // notice that logPy_align *changes* the spike to the actual starting point
        if (read)
            m_logPy.push_back(logPy_align(*read, &(*o),particle->alignmentCache,debug_print));
    }
/////////////////////////////////////////////////
    uint K = m_logPy.size();

    if (debug_print)
    {
    	cout<<"log P(y|o): "<<m_logPy;
    }

    vector<double> logPyo(K);
    vector<double> logPo_(K);
    vector<int> o_(K, 0);  // important to init with 0, because instances rely on that


    for (uint i=0; i<K; i++)
    {
        if (read)
            logPo_[i] = sampleInternalO(m_spikes[i], &o_[i]);
        if (debug_print)
		    	cout<<"For spike "<<i<<", at location "<<m_spikes[i]<<", we sampled o_="<<o_[i]<<endl;
    }

    logPyo = m_logPy + logPo_;
    if (!jvState::with_geometric) 
    	logPyo = m_logPy;

    ///////
    if (jvState::with_temperature)
        logPyo = logPyo*jvState::temperature;
    ///////

    int o = jvSampler::lmnrnd(logPyo);

    jvLogProb logPartition(log_sum_exp(logPyo), logPyo[o]);
    

    int res = m_spikes[o];

    if (debug_print)
    {
	    cout<<"log P(o|s): "<<logPo_;
	    cout<<"log P(y,o|s): "<<logPyo;
	    cout<<"chose entry "<<o<<", at offset "<<res<<" (o_[o] = "<<o_[o]<<")"<<endl;
    }


    particle->o_t = res-o_[o];
    particle->o_i = o_[o];
    particle->o_spike = m_spikes_bu[o];
    return logPartition;

}

double jvSequence::sampleInternalO(int offset, int* res)
{
    jvDVec offsets(jvInstance::WIDTH);
    for (uint o=0; o<jvInstance::WIDTH; o++)
        offsets[o] = logPo(offset-o);

    *res = jvSampler::lmnrnd(offsets);
    return log_sum_exp(offsets)-log(jvInstance::WIDTH);
}


void jvState::setTemperature()
{
    if (with_temperature)
    {
        if (iteration <1000)
            temperature = max(1.0,0.97/(1-0.0009*iteration));
        else
            temperature = 10;
    }
    else
        cout<<"Error: no annealing!\n";
                   
}

//  this functions gets two output parameters (by reference).
//  it also uses two already cached parameters: m_spikes and m_logPy
//  for each spike in m_spikes the function checks which instances contain that spike.
//  If there are instances that contain it, an entry is added to spikeIns (pointer to the
//  instance) and to logPty (the weight to give this instance).
//void jvSequence::computeSpikedInstances(vector<jvInstance*>& spikeIns, jvDVec& logPty)
void jvSequence::computeSpikedInstances(vector<jvParticle>& spikeIns, jvDVec& logPty)
{
    assert(m_spikes.size() == m_logPy.size());
    for (uint i=0; i<m_spikes.size(); i++)
    {
        for (int o=0; o<(int)jvInstance::WIDTH; o++)
        {
            if (m_spikes[i]-o < m_left) break;
            if (m_spikes[i]-o >= m_right) continue;

            set<jvInstance*>& insAtPos = m_ins2[pos(m_spikes[i]-o)];

            for (set<jvInstance*>::iterator it = insAtPos.begin(); it != insAtPos.end(); it++)
            {
                double logPy_t = m_logPy[i];
                logPty.push_back(logPy_t + log((*it)->getN()));
                jvParticle particle;
                particle.ins = *it;
                particle.o_i = o;
                particle.o_spike = m_spikes_bu[i];
                spikeIns.push_back(particle);
            }
        }
    }
}


double jvSequence::logPy(int offset)
{
    if (m_spikes.size() != m_logPy.size())
    {
        assert(m_spikes.size() == 1);
        return m_logPy[offset-m_spikes.front()];
    }

    vector<int>::iterator it = lower_bound(m_spikes.begin(), m_spikes.end(), offset);
    if ((it != m_spikes.end()) && (*it == offset))
        return m_logPy[*it];
    else return -1e100;
}


bool jvSequence::isEqual(jvSequence* s)
{
    if (getRight()-getLeft() != s->getRight()-s->getLeft())
        return false;

    int i=getLeft(), j=s->getLeft();
    for (;i<getRight(); i++, j++)
        if (m_bases[pos(i)] != s->m_bases[s->pos(j)])
            return false;
    return true;
}

bool jvState::isEqual(jvState* st)
{
    if (m_seq.size() != st->m_seq.size())
        return false;

    for (uint s=0; s<m_seq.size(); s++)
    {
        bool tmp = false;
        for (uint s2=0; s2<st->m_seq.size(); s2++)
            if (m_seq[s]->isEqual(st->m_seq[s2]))
            {
                tmp = true;
                break;
            }
        if (!tmp) return false;
    }
    return true;
}



jvDVec jvState::getM(bool sparse, bool with_last)
{
	//  remark:  i should go until m_seq.size()-1, since other functions
	//           depends on that size of the returned vector

    //  if sparse = true, return values corresponding to the entries at m_spiked_seqs.
    //  if sparse = false, return values for all sequences.
    
    jvDVec M;
    vector<jvSequence*>& seq_array = (sparse? m_spiked_seqs : m_seq);
    assert(with_last || (seq_array.back() == m_seq.back()));
    int N = (with_last ? seq_array.size() : seq_array.size()-1);

    for (uint i=0; i<N; i++)
        if (seq_array[i]->m_never_clean)
            M.push_back((double)seq_array[i]->getM());
        else
            M.push_back((double)seq_array[i]->getM());

    if (with_last && (seq_array.back() == m_seq.back()))
    {
        assert(M.back() == 0);
        M.back() = m_gamma;
    }
        
    return M;
}

jvDVec jvState::getN()
{
    jvDVec N;
    for (uint i=0; i<m_ins.size(); i++)
        N.push_back((double)m_ins[i]->getN());
    return N;
}


uint jvState::setSpikesAtZero()
{
    assert(false);
    m_spiked_seqs = m_seq;
    for (int i=0; i<(int)m_spiked_seqs.size(); i++)
    {
    	m_spiked_seqs[i]->m_spikes.clear();
    	m_spiked_seqs[i]->m_spikes.push_back(m_spiked_seqs[i]->getLeft());
    }
    return m_seq.size();
}

// returns the number of spikes found.
// now supports also restricted spikes:  when s is not null it limits the possible spikes
// to be within 10 bases of (s,offset)
uint jvState::computeSpikes(jvBaseVec* ptr, jvSpikePolicy policy)
{
    bool debug_print = false;

    if (debug_print) print_bases(ptr);

    m_spiked_seqs.clear();
    int total_spikes = 0;    
    int N = m_seq.size(); 

    if (policy.code & SPIKE_LAST_TWO_SEQS)
    {
        // force the last two sequences to be spiked.
        assert(N >= 3);
        for (int s = N-3; s<=N-2; s++) // s= N-3, N-2.
        {
            jvSequence* seq = m_seq[s];
            seq->m_inprocess = true;
            seq->m_spikes.clear();
            m_spiked_seqs.push_back(seq);
        }
    }
            
    if (policy.code & SPIKE_HASH)
    {
        int L = ptr->size();
    
        //  sample/choose one or more K-mers from the read
        //  make a list of spikes;
        // UINT here is a bad idea
        int gap_between_kmers = (policy.code & SPIKE_SPARSE ? 10:1);

    
        jvKmer mask = (1<<2*K)-1;
        jvKmer stat = 0; 
        jvKmer Kmer = 0;
        jvBase base;

        for (int o=0; L>=(int)K && o<L; o++)
        {
            base = (jvBase)(*ptr)[o]; 
            stat = stat * nBases + (base == (int)nBases? 0 : nBases-1);
            Kmer = Kmer * nBases + base;
            stat &= mask;
            Kmer &= mask;
            if (stat != mask) continue;
            if (debug_print) cout<<"NEW: Chose Kmer: "<<setbase(nBases)<<Kmer<<endl<<dec;
            if (o % gap_between_kmers) continue;
            
            //  add to spike list
            for (set<jvLocation>::iterator i=m_hash[Kmer].begin(); i!=m_hash[Kmer].end(); i++)
            {
                jvSequence* seq = (*i).first;
                int spike_loc = (*i).second;
                if ( (policy.bin) && (policy.bin != seq->getBin()) )
                    continue;
                if (policy.code & SPIKE_LAST_TWO_SEQS)
                    if ( (seq != m_seq[N-3]) &&  (seq != m_seq[N-2]) )
                         continue;
                if (seq->m_inprocess == false)
                {
                    seq->m_inprocess = true;
                    seq->m_spikes.clear();
                    m_spiked_seqs.push_back(seq);
                }
                seq->m_spikes.push_back(spike_loc-o+K-1);
            }
            
        }
        for (int i=0; i<m_spiked_seqs.size(); i++)
        {
            total_spikes += m_spiked_seqs[i]->pruneSpikes(policy);
            m_spiked_seqs[i]->m_inprocess = false;
        }
    }

    
    if (policy.code & SPIKE_NEW_SEQ)
    {
        // set one spike at an empty seq
        m_spiked_seqs.push_back(m_seq.back());
        m_seq.back()->m_spikes.clear();
        m_seq.back()->m_spikes.push_back(jvAlign::extend);
        total_spikes++;
    }

    if (policy.code & SPIKE_LAST_TWO_SEQS)
    {
        // force the last two sequences to be spiked.
        assert(N >= 3);
        for (int s = N-3; s<=N-2; s++) // s= N-3, N-2.
        {
            jvSequence* seq = m_seq[s];
            seq->m_spikes.push_back(0);
            total_spikes++;
        }
    }

    assert(total_spikes>0);
        
    return total_spikes;
}



//////////////////////////////////

jvLogProb jvState::sampleInstance(jvInstance* ins,  jvReadAss* read, jvSpikePolicy policy, bool debug_print)
{

    assert(read != NULL);

    if (debug_print)
    	cout<<"-----    Sampling a location for "<<(ins ? "an instance" : "a read")<<endl;

    ///////////////  populate m_spikes /////////////////////
    jvBaseVec  ins2read;
    jvBaseVec* ptr;

    if (read)
    {
        ptr = read->getData()->getRawBases();
    }

    if (debug_print)
    {
        cout<<"sampling read/consensus: "<<scientific;
        for (uint i=0; i<(*ptr).size(); i++)
            cout<<(int)(*ptr)[i]+1;
        cout<<endl;
        cout<<"Collecting spikes...: ";
    }


    int nSpikes = computeSpikes(ptr, policy);
    assert(nSpikes > 0);
    m_sum_spikes += nSpikes-1;
    m_nSpikes += nSpikes-1;
    m_count_reads++;
    

    
    jvDVec logPy_s(m_spiked_seqs.size());
    vector<jvParticle> particle(m_spiked_seqs.size());
    
    if (read)
    {
    	for (uint s=0; s<particle.size(); s++)
    	{
            if (debug_print)
                cout<<"Sampling location in sequence "<<s<<": "<<endl;
            particle[s].alignmentCache = read->getData()->getAlignmentCache();
            jvLogProb logPyo_s = m_spiked_seqs[s]->sampleO(NULL, ptr, &particle[s], debug_print);
            logPy_s[s] = logPyo_s.partition;
    	}
    }

    //  add sequence weights
    jvDVec M = getM(true, true);
    
    jvDVec logPs = log(M);
    if (m_withInstances)
        logPs = logPs + (-log(m_ins.size())-1) + (-log(m_gamma));
    
    jvDVec logPsy = logPs + logPy_s;
    if ( alignment_only)
	logPsy = logPy_s;  // don't include DP during launch state computation

    ///////
    if (jvState::with_temperature)
        logPsy = logPsy*jvState::temperature;  // annealing
    //////

    if (debug_print)
    {
        cout<<"M: "<<M;
        cout<<"logP(s): "<<logPs;
        cout<<"logP(y|s): "<<logPy_s;
        cout<<"logP(s,y): "<<logPsy;
        cout<<"o(s): ";
        for (uint s=0; s<m_spiked_seqs.size(); s++)
            cout<<particle[s].o_t<<" ";
        cout<<endl;
    }
    
    uint s = jvSampler::lmnrnd(logPsy);
    if (seq_choice > -1)
    {
        assert(seq_choice < (int)logPsy.size());
        s = seq_choice;
    }
    // Collect transition probability
    
    if (alignment_only)
    {
        double garbage_per_base = logPsy.back()/ptr->size();
        double price_per_base = logPsy[s]/ptr->size();
	if (price_per_base < garbage_thresh)
            s = logPsy.size()-1;

    }

    jvLogProb log_partition(log_sum_exp(logPsy), logPsy[s]);


    if (read) //  Sampling an instance for a read
    {
        jvDVec logPty;
        vector<jvParticle> spikeIns;

        //  fill the values of spikeIns and logPty  ///////
        if (m_withInstances)
        {
            assert(false);
            for (uint j=0; j<m_spiked_seqs.size(); j++)
            {
                m_spiked_seqs[j]->computeSpikedInstances(spikeIns, logPty);
            }
        }

        logPty.push_back(log_partition.partition + log(m_alpha));
        //////////////////////////////////////////////////


        if (debug_print)
            cout<<"logP(t,y): "<<logPty;

        assert(logPty.size() == 1);
        uint t = 0;

        if (t == logPty.size()-1)
        {
            m_spiked_seqs[s]->setFreeze(false);
            ins = newIns(m_spiked_seqs[s], particle[s].o_t);
            ins->assign(true);  // this line changes the state!
            particle[s].ins = ins;
            spikeIns.push_back(particle[s]);
        }
    	read->getData()->resetAlignment();
    	jvSequence* seq = spikeIns[t].ins->getSeq();
    	int offset = spikeIns[t].ins->getOffset() + spikeIns[t].o_i;
    	int out_o = spikeIns[t].o_spike;
    	seq->logPy_align(*ptr, &out_o, read->getData()->getAlignmentCache(),false);
    	if ( (out_o != offset) && (logPsy[s] > -1e99) )
    	{
    	    print_bases(&read->getData()->getBases());
    	    cout<<"out_o="<<out_o<<"  offset="<<offset<<"  spikeIns[t].o_spike="<<spikeIns[t].o_spike<<endl;
    	    out_o = spikeIns[t].o_spike;
            seq->logPy_align(*ptr, &out_o, read->getData()->getAlignmentCache(),true);
            cout<<"total spikes: "<<nSpikes<<"  chose sequence "<<s<<endl;
            cout<<"M: "<<M;
            cout<<"logP(s): "<<logPs;
            cout<<"logP(y|s): "<<logPy_s;
            cout<<"logP(s,y): "<<logPsy;
            cout<<"o(s): ";
            for (uint s=0; s<m_spiked_seqs.size(); s++)
                cout<<particle[s].o_t<<" ";
            cout<<endl;
            assert(false);
    	}

    	read->getData()->updateAlignment(true, debug_print);
        read->setIns(spikeIns[t].ins, spikeIns[t].o_i);


    }

    if (m_seq.back()->getM() > 0)
    {
        m_seq.back()->setBin(policy.bin);
        newSeq()->buildHash(&m_hash);
    }

    assert(m_seq.back()->getM() == 0);

    return log_partition;

}


void jvSequence::resetEditSummary(void)
{
	m_editSummary.resize(m_map.size()+1);
	for (uint i=0;i<m_map.size();i++) m_editSummary[i].reset();
        m_insCount = 0;
        m_delCount = 0;
}

void jvState::summarizeEdits(void)
{
    jvInstance* ins;
    jvSequence* seq;
    int ofs;

    for (uint s=0;s<m_seq.size();s++)
    {
        m_seq[s]->resetEditSummary();
    }

    for (uint i=0;i<m_reads.size();i++)
    {
        if (!m_reads[i].m_isAssigned) continue;
        ins = m_reads[i].getIns();
        seq = ins->getSeq();
        ofs = m_reads[i].getOffset() + ins->getOffset();
        seq->addEdit(ofs,m_reads[i].getData()->getEdit());
    }
}

void jvSequence::addEdit(int ofs,vector<char>* edit, bool print)
{
    print = false;
    int curr = ofs-m_left;
    bool addDelete = true;

    for (uint i=0;i<edit->size();i++)
    {
        if (print) cout << (int)(*edit)[i];
        if (((*edit)[i] == -2) && addDelete)
        {
            if (m_editSummary.size() < curr+1) m_editSummary.resize(curr+1);
            m_editSummary[curr] += jvEdit(0,1);
            addDelete = false;
            m_delCount++;
        }

        if ((*edit)[i] == -1)
        {
            curr++;
            addDelete = true;
        }

        if ((*edit)[i] >= 0)
        {
            if (m_editSummary.size() < curr+1) m_editSummary.resize(curr+1);
            m_editSummary[curr] += jvEdit(1,0);
            curr++;
            addDelete = true;
            m_insCount++;
        }
    }
    if (print) cout << endl;
}

void jvSequence::propagateInsertToRead(int ofs, jvReadAss* ra)
{
	int st = ra->getOffset() + ra->getIns()->getOffset();
	int pos = ofs-st;
	ra->getData()->insertAtOffset(pos);
	ra->m_oneSpike = false;
}

void jvSequence::propagateDeleteToRead(int ofs, jvReadAss* ra)
{
	int st = ra->getOffset() + ra->getIns()->getOffset();
	int pos = ofs-st;
	ra->getData()->deleteAtOffset(pos);
	ra->m_oneSpike = false;
}

bool jvSequence::insertAtOffset(int ofs)
{
    if (pos(ofs) >= m_map.size())
    {
        cout<<"insert after end of sequence!\n";
        cout << "insert at offset: " << ofs  <<  "(total reads at offset: " << m_map[pos(ofs)] << ", total deletions at offset:" << m_editSummary[pos(ofs)].dl << ")\n";
        return false;
    }

    bool debug_print = false;
    if (debug_print)
        cout << "insert at offset: " << ofs  <<  "(total reads at offset: " << m_map[pos(ofs)] << ", total deletions at offset:" << m_editSummary[pos(ofs)].dl << ")\n";
    vector<jvInstance*> instancesToShift(0);
    vector<jvInstance*> instancesToUnfreeze(0);
    vector<jvReadAss*> readsToShift(0);
    vector<jvReadAss*> readsToUnassign(0);

    for (uint i=0;i<m_ins2.size();i++)
    {
        if (m_ins2[i].size())
        {
            for (set<jvInstance*>::iterator insit = m_ins2[i].begin(); insit != m_ins2[i].end();insit++)
            {
                assert((*insit)->getN() == (int)(*insit)->m_reads.size());
                if ((*insit)->getOffset() >= ofs)
                // shift right all the instances that start after ofs
                {
                    instancesToShift.push_back(*insit);
                }
                else
                {
                    for (set<jvReadAss*>::iterator it = (*insit)->m_reads.begin();it != (*insit)->m_reads.end();it++)
                    {
                        // for instances that start before ofs, unassign all reads that overlap with ofs
                        if (((*it)->getOffset() + (*it)->getIns()->getOffset()<=ofs) &&
                            ((*it)->getOffset() + (*it)->getIns()->getOffset() + ((int)(*it)->getData()->getBases().size())-1 >=ofs))
                        {
                            if ((*it)->m_isAssigned)
                            {
                                readsToUnassign.push_back((*it));
                            }
                        }
                        // shift right all reads that start after ofs
                        if ((*it)->getOffset() + (*it)->getIns()->getOffset()>ofs)
                        {
                            readsToShift.push_back((*it));
                        }
                    }
                    instancesToUnfreeze.push_back(*insit);
                }
            }

        }
    }

    if (false) //(m_map[pos(ofs)] > (int)readsToUnassign.size())
    {
        debug_print = true;
        cout<<"Offset: "<<ofs<<"  readsToUnasign.size()= "<< readsToUnassign.size()<< "m_map[pos(ofs)]="<< m_map[pos(ofs)]<<endl;
        print();
        printCounts();
        for (uint i=0;i<readsToUnassign.size();i++)
        {
            cout<<"read "<<i<<"  o= "<<readsToUnassign[i]->getOffset()+readsToUnassign[i]->getIns()->getOffset()<<endl;
            print_bases(&readsToUnassign[i]->getData()->getBases());
        }
    }

    if (debug_print)
    {
        cout << "instances to shift: " << instancesToShift.size() << " reads to reassign: " << readsToUnassign.size()<<"  ";
        cout << "instances to unfreeze: " << instancesToUnfreeze.size() << " reads to shift: " << readsToShift.size()<<"  "<<endl;
    }


    for (uint i=0;i<instancesToShift.size();i++)
    {
        instancesToShift[i]->assign(false);
        instancesToShift[i]->m_offset++;
        instancesToShift[i]->assign(true);
    }

    for (uint i=0;i<readsToUnassign.size();i++)
    {
        readsToUnassign[i]->assign(false);
        propagateInsertToRead(ofs,readsToUnassign[i]);
        readsToUnassign[i]->assign(true); // added!!
    }

    for (uint i=0;i<readsToShift.size();i++)
    {
        readsToShift[i]->assign(false);
        readsToShift[i]->setIns(readsToShift[i]->getIns(),readsToShift[i]->getOffset()+1);
        readsToShift[i]->assign(true);
    }

    sampleBases();
    return true;
}

void jvSequence::deleteAtOffset(int ofs)
{
    bool debug_print = false;
    if (debug_print)
    	cout << "delete at offset: " << ofs  <<  "(total reads at offset: " << m_map[pos(ofs)] << endl;
    vector<jvInstance*> instancesToUnfreeze(0);
    vector<jvInstance*> instancesToShift(0);
    vector<jvReadAss*> readsToShift(0);
    vector<jvReadAss*> readsToUnassign(0);

    for (uint i=0;i<m_ins2.size();i++)
    {
        if (m_ins2[i].size())
        {
            for (set<jvInstance*>::iterator insit = m_ins2[i].begin(); insit != m_ins2[i].end();insit++)
            {
                assert((*insit)->getN() == (int)(*insit)->m_reads.size());
                if ((*insit)->getOffset() > ofs)   //  shift back instances that start after ofs
                {
                    instancesToShift.push_back(*insit);
                }
                else
                {
                    //  unassign all reads that overlap this ofs
                    for (set<jvReadAss*>::iterator it = (*insit)->m_reads.begin();it != (*insit)->m_reads.end();it++)
                    {
                        if (((*it)->getOffset() + (*it)->getIns()->getOffset()<=ofs) &&
                            ((*it)->getOffset() + (*it)->getIns()->getOffset() + ((int)(*it)->getData()->getBases().size())-1 >=ofs))
                        {
                            assert((*insit)->getOffset() <= ofs);
                            if ((*it)->m_isAssigned)
                            {
                                readsToUnassign.push_back((*it));
                            }
                        }
                        //  shift back reads that start after ofs in instances that start before or at ofs
                        if ((*it)->getOffset() + (*it)->getIns()->getOffset()>ofs)
                        {
                            readsToShift.push_back((*it));
                        }

                    }
                    instancesToUnfreeze.push_back(*insit);
                }
            }

        }
    }
    assert(m_map[pos(ofs)] <= (int)readsToUnassign.size());
    if (debug_print)
    {
        cout << "instances to shift: " << instancesToShift.size() << " reads to reassign: " << readsToUnassign.size()<<"  ";
        cout << "instances to unfreeze: " << instancesToUnfreeze.size() << " reads to shift: " << readsToShift.size()<<"  "<<endl;
    }

    for (uint i=0;i<instancesToShift.size();i++)
    {
        instancesToShift[i]->assign(false);
        instancesToShift[i]->m_offset--;
        instancesToShift[i]->assign(true);
    }
    for (uint i=0;i<readsToUnassign.size();i++)
    {
        readsToUnassign[i]->assign(false);
        propagateDeleteToRead(ofs,readsToUnassign[i]);
        readsToUnassign[i]->assign(true);
    }
    for (uint i=0;i<readsToShift.size();i++)
    {
        readsToShift[i]->assign(false);
        assert(readsToShift[i]->getOffset()-1>=0);
        readsToShift[i]->setIns(readsToShift[i]->getIns(),readsToShift[i]->getOffset()-1);
        readsToShift[i]->assign(true);
    }

    sampleBases();
}


void jvState::checkAllOffsets(void)
{
	walign = false;
	freezeInstances(false);
        bool debug_print = false;

	summarizeEdits();
	int nseqs = m_seq.size();

	for (uint i=0;i<m_seq.size();i++)
	{
		if (m_seq.size() != nseqs)
			cout<<"stop here\n";

		if (m_seq[i]->getM())
		{
			int ofs = m_seq[i]->checkAllOffsets();
		        if (ofs && debug_print) cout << " -  seq " << i << "\n";
		}
	}
	walign = true;
}

bool jvSequence::checkAllOffsets(void)
{
    vector<int> edit;
    vector<bool> action;

    int left = getLeft();
    int right = getRight();

    for (int i=left; i<right; i++)
    {
        if (m_map[pos(i)]<=m_editSummary[pos(i)].in && m_editSummary[pos(i)].in && !( (m_map[pos(i)] == 1) && (m_editSummary[pos(i)].in ==1) ))
        {
            //  position for deletion
            edit.push_back(i);
            action.push_back(true);
        }
        else
            if (max(pos(i-1)>=0 ? m_map[pos(i-1)]:0, m_map[pos(i)])/2 < m_editSummary[pos(i)].dl)
            {
                //  position for insertion
                edit.push_back(i);
                action.push_back(false);
            }
    }

    int correction = 0;
    for (int j=0; j<edit.size(); j++)
    {
        if (action[j])
        {
            deleteAtOffset(edit[j] + correction);
            correction--;
        }
        else
        {
            if (!insertAtOffset(edit[j] + correction))
            {
                //  code for debugging purposes incase of assertion failure in insert at offset.
                cout<<"insert at offset failed!\n";
                for (int i=0; i<edit.size(); i++)
                    cout<<"loc: "<<edit[i]<<"  action: "<<(action[i]?"del":"ins")<<endl;
                assert(false);
            }

            correction++;
        }
    }

    return edit.size()>0;
}


int findProp = 0;
int findAcc = 0;


int problematicMerge = 0;

bool jvState::updateHashHits(jvKmer &Kmer, jvSequence* seqI,map<jvSequence*,int> &hits,set<jvSequence*> &seen,map<jvSequence*,bool> &touched, bool rev)
{
    jvSequence *aSeq;
    for (set<jvLocation>::iterator it = (m_hash)[Kmer].begin();it != (m_hash)[Kmer].end();it++)
    {
        aSeq = it->first;
        if (aSeq->getM() == 1) continue;  // don't bother with singletons
        bool is_edge = (rev ? it->second >= aSeq->getRight() - endAl : it->second < aSeq->getLeft() + endAl);
        // if the kmer appears in the first (last, if rev = true) endAl bases of aSeq
        if (aSeq != seqI && !touched.find(aSeq)->second && is_edge)
        {
            int curr = hits.find(aSeq)->second; // find aseq index
            hits.find(aSeq)->second = hits.find(aSeq)->second + 1;
            assert(hits.find(aSeq)->second );
            seen.insert(aSeq);
        }
    }
}

bool jvState::findAndAlignEnds()
{
    bool debug_print = false;
    int ct =0 ;
    vector<int> perm = jvSampler::randomPermutation(m_seq.size()-1);
    jvDVec prop(2*m_seq.size());

    vector<jvSequence*> seqs(m_seq);

    buildHash(); // based on seqs

    uint K = jvState::K;

    jvKmer mask = (1<<2*K)-1;
    jvKmer stat = 0, revStat = 0;
    jvKmer Kmer = 0, revKmer = 0;
    jvBase base, revBase;
    jvSequence *aSeq;

    set<jvLocation> locs;

    set<jvSequence*> seen;
    set<jvSequence*> revSeen;
    map<jvSequence*,int> hits;
    map<jvSequence*,int> revHits;
    map<jvSequence*,int> reverse; // reverse index not reverse sequence
    map<jvSequence*,bool> touched;
    bool gotOne = false;

    for (int i=0;i<m_seq.size();i++)
    {
        hits.insert(pair<jvSequence*,int>(seqs[i],0));
        revHits.insert(pair<jvSequence*,int>(seqs[i],0));
        reverse.insert(pair<jvSequence*,int>(seqs[i],i));
        touched.insert(pair<jvSequence*,bool>(seqs[i],false));
    }

    touched.find(seqs[seqs.size()-1])->second = true;

    // for each  sequence
    for (vector<int>::iterator pi=perm.begin(); pi<perm.end(); pi++, ct++)    
    {
    	int i = *pi;
        debug_print = false;

        if (debug_print) cout<<"Processing seq "<<i<<endl;

        if (seqs[i]->m_never_clean)
            continue;

        //  if it's been touched - continue.
    	if (touched.find(seqs[i])->second)
    	{
            if (debug_print) cout<<"already touched - abort\n";
            continue;
    	}

        jvSequence *seqI = seqs[i];

        //  if it's too short - continue
        if (seqI->getM() == 1)
        {
            if (debug_print) cout<<"singleton - abort\n";
            continue;
        }

        if (seqI->getRight()-seqI->getLeft() < endAl)
        {
            if (debug_print) cout<<"too short - abort\n";
            continue;
        }

        // reset score
        for (int ii=0;ii<2*seqs.size();ii++)
        {
            prop[ii] = -1e+150;
        }

        seen.clear();
        revSeen.clear();

        //reset hits to 0
        for (map<jvSequence*,int>::iterator hit = hits.begin();hit != hits.end();hit++)
        {
            hit->second = 0;
        }

        for (map<jvSequence*,int>::iterator hit = revHits.begin();hit != revHits.end();hit++)
        {
            hit->second = 0;
        }

        // prepare area to look for overlap (endI = the last endAl bases of seqI).
        int stI,enI;
        stI = seqI->getRight() - endAl; enI = seqI->getRight();
        jvBaseVec endI(seqI->m_bases.begin() + seqI->pos(stI),seqI->m_bases.begin() + seqI->pos(enI));

        mask = (1<<2*K)-1;
        stat = 0; revStat = 0;
        Kmer = 0; revKmer = 0;


        // get all kmers from the endI = last enaAl bases of the current sequence
        // and collect the hits they have on the rest of the sequences
        for (int l=0;l<endAl;l++)
        {
            base = endI[l];
            revBase = complement_base(endI[endAl - 1 - l]);

            stat = stat * nBases + (base == (int)nBases? 0 : nBases-1);
            Kmer = Kmer * nBases + base;
            stat &= mask;
            Kmer &= mask;

            revStat = revStat * nBases  + (revBase == (int)nBases? 0 :nBases-1);
            revKmer = revKmer * nBases + revBase;
            revStat &= mask;
            revKmer &= mask;

            if (stat == mask)
            {
                // a hit is Kmer appearing in the first endAl letters of another seq
                updateHashHits(Kmer,seqI,hits,seen,touched, false);
            }

            if (both_orientations && revStat == mask)
            {
                // a hit is revKmer appearing in the *last* endAl letters of another seq
                updateHashHits(revKmer,seqI,revHits,revSeen,touched, true);
            }
        }

        int best = 0; int bestSeq; int bestDir;
        int count = 0;
        int j=-1;

        //  for each sequence that has kmer hits calculate probability for merge proposal
        for (set<jvSequence*>::iterator sit=seen.begin();sit != seen.end();sit++)
        {
            aSeq = *sit;
            if (touched[aSeq]) continue;
            
            count = hits[aSeq];
            if (count>best)
                best = count;
            
            prop[reverse[aSeq]*2] = count;
        }

        for (set<jvSequence*>::iterator sit=revSeen.begin();sit != revSeen.end();sit++)
        {
            aSeq = *sit; 
            if (touched[aSeq]) continue;
            
            count = revHits[aSeq];
            if (count>best)
                best = count;

            prop[reverse[aSeq]*2+1] = count;
        }

        if (debug_print) cout<<"Best: "<<best<<endl;//<<"  prop: "<<prop;

        // propose up to 20 merges, or until one is accepted.
        if (best>5)
        {
            int tries = 0;
            int dirj = 0;
            bool gotOne = false;
            int lenI = seqs[i]->getRight()-seqs[i]->getLeft();
            
            while(tries<3 && !gotOne)
            {
                j = jvSampler::lmnrnd(prop, true);
                if (prop[j]>12)                
                    prop[j] = -1e+150;
                else
                    break;
                dirj = j % 2; j = (j-dirj)/2;
                tries++;
                
		if (debug_print)                
                    cout << "Trying align ends of :" << i << "," << j << (dirj ? " - reverse" : "")<<endl;

                int nreadsI = seqs[i]->getM();
                int nreadsJ = seqs[j]->getM();
		if (alignEnds(seqs[i],seqs[j],dirj,debug_print))
                {
                    if (debug_print)
                    {
                        cout << "merged " << i << "("<<nreadsI<<" reads)," << j << (dirj ? " - reverse" : "") << "("<<nreadsJ<<" reads).  ";
                        cout << "New sequence has length "<<seqs[i]->getRight()-seqs[i]->getLeft() << "(before "<<lenI<<")."<<endl;
                    }
                    touched[seqs[i]] = true;
                    touched[seqs[j]] = true;
                    gotOne = true;
		        }
            }
        }
    }

    hits.clear();revHits.clear();
    reverse.clear();
    seqs.clear();
    seen.clear();revSeen.clear();
    touched.clear();

    return true;
}

void jvInstance::revComplement()
{
    assert(m_reads.size() == 1);
    jvReadAss* rd = *m_reads.begin();
    assert(rd->getOffset() == 0);
    assign(false);
    vector<jvReadAss*> readsToFlip;
    for (set<jvReadAss*>::iterator it = m_reads.begin();it != m_reads.end();it++)
    {
        readsToFlip.push_back(*it);
    }
    for (int i=0; i<readsToFlip.size(); i++)
    {
        readsToFlip[i]->assign(false);;
        readsToFlip[i]->revComplement();
        readsToFlip[i]->assign(true);;
    }
    int L = rd->getData()->getBases().size();
    assert(getRight() == getOffset() + L);
    m_offset = -getRight();
    assign(true);
}

void jvReadAss::revComplement()
{
    assert(m_isAssigned == false);
    m_data->complementSelf();
}


void jvSequence::revComplement()
{
    setFreeze(false);
    vector<jvInstance*> instancesToMove;
    for (uint i=0;i<m_ins2.size();i++)
        if (m_ins2[i].size())
            for (set<jvInstance*>::iterator insit = m_ins2[i].begin(); insit != m_ins2[i].end();insit++)
            {
                instancesToMove.push_back(*insit);
            }

    for (int i=0; i<(int)instancesToMove.size(); i++)
    {
        instancesToMove[i]->revComplement();
    }
    sampleBases();  // TODO: the correct way is to restore the complement of the original
}

// gets an edit string for aligning seqJ to seqI,
// and coverts it to two edit strings, one for each sequence,
// and also provides an estimate for the likelihood of doing the merge
// according to these alignments.
double jvState::computeCountsAlignmentScore(vector<char> edit, int local_start,
                                            jvBaseVec seqI,
                                            jvBaseVec seqJ,
                                            jvBaseCountsVec countsI,
                                            jvBaseCountsVec countsJ,
                                            vector<char> &editI,
                                            vector<char>  &editJ)
{
    bool debug_print = false;
    // seqI is padded from the right; so is countsI
    int pos_j = 0;
    int pos_i = local_start;
    int edpos_i = local_start;
    int edpos_j = 0;
    int editLen = edit.size();
    double score = 0;
    double delPenalty = errorModel.logIndel();
    double baseSaving = -logBasePrior;
    int bases_saved = 0;
    
    editI.insert(editI.begin(),local_start,-1); // everything before local_start gets copied

    for (int i=0;i<editLen && pos_i < countsI.size();i++) // only work on overlap
    {
        int totI = sum(countsI[pos_i]);
        int totJ = sum(countsJ[pos_j]);

        if (edit[i]==-1)
        {
            
            if (seqI[pos_i] != seqJ[pos_j])
            {
                jvBase bJ = seqJ[pos_j];
                jvBase bI = seqI[pos_i];
                if ((bJ != nBases) &&  (bI != nBases))
                {
                    double score_old = score;
                    
                    if (totI>=totJ)
                    {
                        score = score - jvSampler::logPy_base_observed(countsJ[pos_j], bJ)
                                      + jvSampler::logPy_base_observed(countsJ[pos_j], bI);
                    }
                    else
                    {
                        score = score - jvSampler::logPy_base_observed(countsI[pos_i], bI)
                                      + jvSampler::logPy_base_observed(countsI[pos_i], bJ);
                    }

                    if (debug_print)
                    {
                        cout<<scientific<<"posi="<<pos_i<<" bI="<<(int)bI<<" totI="<<totI<<endl;
                        cout<<"posj="<<pos_j<<" bJ="<<(int)bJ<<" totJ="<<totJ<<endl;
                        cout<<"score += "<<score-score_old<<endl;
                    }
                }
            }
            if ((seqI[pos_i] != nBases) && (seqJ[pos_j] != nBases))
            {
                bases_saved++;
                score = score + baseSaving;
            }
            pos_i++; editI.push_back(-1); edpos_i++;
            pos_j++; editJ.push_back(-1); edpos_j++;
        }

        if (edit[i]==-2)
        {
            double score_old = score;

            if (debug_print)
            {
                cout<<"posi="<<pos_i<<" totI="<<totI<<endl;
                cout<<"posj="<<pos_j<<" totJ="<<totJ<<endl;
            }
            if (totI>=totJ) // delete in J
            {
                editJ.push_back(-2); edpos_j++;
                score = score + totJ*delPenalty;
                pos_j++;
            }
            else // insert into I
            {
                editI.push_back(seqJ[pos_j]); edpos_i++;
                pos_j++;
                score = score + totI*delPenalty;
            }
            if (debug_print)
                cout<<"score += "<<score-score_old<<endl;                      
        }

        if (edit[i]>=0)
        {
            double score_old = score;
            if (debug_print)
            {
                cout<<"posi="<<pos_i<<" totI="<<totI<<endl;
                cout<<"posj="<<pos_j<<" totJ="<<totJ<<endl;
            }
            if (totI>=totJ) // insert into J
            {
                editJ.push_back(seqI[pos_i]); edpos_j++;
                pos_i++;
                score = score + totJ*delPenalty;
            }
            else // delete in I
            {
                editI.push_back(-2); edpos_i++;
                score = score + totI*delPenalty;
                pos_i++;;
            }
            if (debug_print)
                cout<<"score += "<<score-score_old<<endl;

        }
    }

    for (int i=pos_j;i<seqJ.size();i++)
    {
        editJ.push_back(-1); edpos_j++;
    }

    if (debug_print)
    {
        cout<<"Total saved bases: "<<bases_saved<<endl;
        print_edit(edit);
        print_edit(editI);
        print_edit(editJ);
    }
    return score;
}

void jvSequence::applyEdit(vector<char> edit, int st)
{
    bool debug_print = false;
    bool gotChange = false;
    int pos = st;
    for (int i=0;i<edit.size();i++)
    {
         if (edit[i]==-1)
         {
             pos++;
         }

         if (edit[i]==-2)
         {
//             cout<<"del at "<<pos<<endl;
             deleteAtOffset(pos);
             gotChange = true;
             if (debug_print)
                 print();
         }

         if (edit[i]>=0)
         {
//             cout<<"ins at "<<pos<<endl;
             if (!insertAtOffset(pos))
             {
                 print();
                 assert(false);
             }
             gotChange = true;
             pos++;
             if (debug_print)
                 print();
         }
    }
}


//  propose to merge end of seqI with beginning of seqJ, (or with beginning of flipped seqJ if dirj=1);
bool jvState::alignEnds(jvSequence *seqI, jvSequence *seqJ, int dirj, bool out_print)
{
    assert(seqI != seqJ);
    bool flip_back = dirj;
    if (dirj)
    {
        assert(both_orientations);
        seqJ->revComplement();
        dirj = 0;        
    }

    if (false)//out_print)
    {
        cout << endl <<  "PRE: " << endl;
        seqI->print(PRINT_BASES);
        cout << endl <<  "PRE: " << endl;
        seqJ->print(PRINT_BASES);
    }
    

    bool gotOne = false;

    //  collect ending of seq i and beginning of seq j
    // dirj==0: the end of seqI is merged with beginning of seqJ
    // dirj==1: seqJ is flipped, so the end of seqI is merged with the flipped end of seqJ
    int stI,enI,stJ,enJ;

    // counts for both sequences, used to choose "optimal" changes to sequences
    jvBaseCountsVec countsI, countsJ;


    // always the same regardles of dirj
    stI = seqI->getRight() - endAl; enI = seqI->getRight();

    countsI.insert(countsI.begin(),
                   &seqI->m_counts[stI-seqI->m_left],
                   &seqI->m_counts[enI-seqI->m_left]);

    // depends on dirJ
    if (dirj == 0)
    {
        stJ = seqJ->getLeft();

        enJ = min(seqJ->getLeft() + endAl,
                  seqJ->getRight());

        countsJ.insert(countsJ.begin(),
		       &seqJ->m_counts[stJ - seqJ->m_left],
                       &seqJ->m_counts[enJ - seqJ->m_left]);
    }
    else
    {
        assert(false);
    }

    jvBaseVec endI(seqI->m_bases.begin() + seqI->pos(stI),
                   seqI->m_bases.begin() + seqI->pos(enI));

    jvBaseVec begJ(seqJ->m_bases.begin() + seqJ->pos(stJ),
                   seqJ->m_bases.begin() + seqJ->pos(enJ));


    endI.resize(2*endAl,nBases);
    countsI.resize(2*endAl,emptyBaseCount);

    //  We are trying to align seq j to the right of seq i
    jvAlign::merge_aligner.align(&endI,&begJ, out_print);
    int local_start = jvAlign::merge_aligner.getLocalStart();
    int overlap = min(enI-stI - local_start, enJ-stJ);    

    if (out_print)
        cout<<"Overlap:" <<overlap;

    // check score
    // if score is ok go through the editR and find first difference (insert or delete) and apply it

    int editLen = jvAlign::merge_aligner.getEditLen();
    vector<char> edit(jvAlign::merge_aligner.getEdit(), jvAlign::merge_aligner.getEdit()+editLen);
    vector<char> editJ, editI;


//    freezeInstances(false);
    double delta =  computeCountsAlignmentScore(edit,local_start, endI, begJ, countsI,countsJ, editI,editJ);
    delta = delta - log(pow(2,40)); // every merge colst you log(gamma)

    if (out_print)
        cout<<"Predicting this much improvement in score: "<<delta<<endl;
    
    if (delta > 0 && overlap >= 20)
    {
        gotOne = true;
        seqI->applyEdit(editI,stI);
        if (dirj)
        {
            seqJ->revComplement();
            seqJ->applyEdit(editJ,-(enJ-1));
        }
        else
        {
            seqJ->applyEdit(editJ,stJ);
        }

        seqI->merge_simple(seqJ, enI - endAl + local_start);
        seqI->sampleBases();

        if (false)
        {
            cout <<  "POST: " << endl;
            seqI->print(PRINT_BASES);
        }
    }
    else
        if (flip_back)
            seqJ->revComplement();  // flip back sequence in case align did not happen 


    return gotOne;
}


//////////////////////////////////////




//  Read methods

jvReadAss::jvReadAss(jvReadData* data): m_data(data), m_ins(NULL), m_offset(0), m_isAssigned(false)
{
    m_oneSpike = false;
    m_skip = false;
}

void jvReadAss::setIns(jvInstance* ins, int offset)
{
    assert(!m_isAssigned);
    assert( (offset >=0));
    if( ! ((offset >=0) && (offset < (int)jvInstance::WIDTH)))
		out_of_bounds++;

    m_ins = ins;
    m_offset = offset;
}

//  returns true iff something has changed.
bool jvReadAss::assign(bool isAssigned)
{
    if (m_isAssigned == isAssigned)
        return false;

    assert(m_ins != NULL);
//    assert(!m_ins->m_isReadFrozen);

    m_isAssigned = isAssigned;
    int action = (isAssigned ? 1 : -1);
    m_ins->updateCounts(action);

    if (action==-1)
    {
    	m_ins->m_reads.erase(this);
    }
    else
    {
    	m_ins->m_reads.insert(this);
    }

    if (m_ins->isAssigned())
        m_ins->getSeq()->updateReadCounts(m_data->getBases(), m_ins->getOffset()+m_offset, action);

    return true;

}


//  Instance methods


// C-tor
jvInstance::jvInstance(jvSequence* seq, int offset)
:m_N(0), m_offset(offset), m_seq(seq), m_isAssigned(false)
{
}


void jvInstance::setSeq(jvSequence* seq, int offset)
{
    if (m_isAssigned == 0)
    {
        m_seq = seq;
        m_offset = offset;
    }
    else cout<<"error!  setting an assigned instance\n";
};


void jvInstance::updateCounts(int action, int count)
{
    m_N += action * count;
}

void jvInstance::assign(bool isAssigned)
{
    if (m_isAssigned == isAssigned)
        return;

    assert(m_seq != NULL);
    assert(!m_seq->m_isInstanceFrozen);

    int action = (isAssigned ? 1 : -1);
    m_isAssigned = isAssigned;
    m_seq->updateInstanceCounts(this, m_offset, action);

//  read is still in the instance, but need to update the contribution to sequence stats.
    if (m_N > 0)
    {
        assert(m_N == 1);
        jvReadAss* r = *(m_reads.begin());
        m_seq->updateReadCounts(r->getData()->getBases(), getOffset()+r->getOffset(), action);
    }

}

void jvInstance::print()
{
    if (m_N > 1)
        cout<<m_reads.size()<<"*";
    jvReadAss* r = *(m_reads.begin());
    cout<<r->getData()->getBases().size()<<" ";
}


void jvInstance::reg()
{
    if ( (m_isAssigned) && (m_seq->m_isInstanceFrozen))
        m_seq->m_ins.push_back(this);
}


//  Sequence methods

void jvSequence::init()
{
    m_M = 0;
    m_sumO = 0;
    m_left = 0;
    m_right = 1;
    updateP();
    m_isInstanceFrozen = false;
    m_hash = NULL;
    m_editSummary.resize(0);
    m_never_clean = false;
    m_insCount = 0;
    m_delCount = 0;
    m_inprocess = false;
    m_getRight_dirty = true;
    m_bin = 0;
}


jvSequence::jvSequence()
    : m_bases(1, nBases), m_map(1,0), m_counts(1, emptyBaseCount), m_ins2(1, emptyInsSet)
{
    init();
}


jvSequence::jvSequence(jvBaseVec* bb_vec)
    : m_bases(1, nBases), m_map(1,0), m_counts(1, emptyBaseCount), m_ins2(1, emptyInsSet)
{
    init();
    int N = (*bb_vec).size();
    assert (N > 0);

    m_sumO = 100000;
    updateP();
    updateReadCounts(*bb_vec, -N/2, +1);

    setFreeze(true);
    center();
    m_never_clean = true;
}


void jvSequence::updateP()
{
    double a = pseudoA;
    double b = pseudoB;
    m_p = (m_M +a) / (m_M + m_sumO + a + b);
    m_logP = log(m_p);
    m_logInvP = log(1-m_p);
}

void jvSequence::updateInstanceCounts(jvInstance* ins, int offset, int action, int count)
{
    m_M += action * count;
    m_sumO += action * count * abs(offset);
    if (!constantP) updateP();
    pad(offset, offset+1);
    if (action == 1)
    	m_ins2[pos(offset)].insert(ins);
    if (action == -1)
    	m_ins2[pos(offset)].erase(ins);
    assert(  ( !m_ins2[pos(offset)].empty() ) || (m_ins2[pos(offset)].begin() == m_ins2[pos(offset)].end() ) );
}

//  Instance Version
void jvSequence::updateReadCounts(jvBaseCountsVec& newCounts, int offset, int action)
{
    assert(false);

     uint left = offset, right = offset + newCounts.size();
     pad(left, right);

     bool update = false;

     for (uint i=0; i<newCounts.size(); i++)
     {
    	if (sum(newCounts[i]) == 0) continue;
     	if (m_map[pos(offset + i)] == 0)
     	{
                	update = true;
                	break;
     	}

     	if (m_map[pos(offset + i)] == -action*sum(newCounts[i]))
     	{
                	update = true;
                	break;
        }
     }


     if (update)
     	updateHash(left, right, -1);

     //  Here we can assert that left, right are non-negative

     for (uint i=0; i<newCounts.size(); i++)
     {
         int sum_counts_i = sum(newCounts[i]);
         if (sum_counts_i == 0) continue;

         //  setting bases to maximal element
         if (m_map[pos(offset + i)] == 0)
         {
             jvBaseCounts::iterator max_base = max_element(newCounts[i].begin(), newCounts[i].end());
        	 m_bases[pos(offset + i)] = distance(newCounts[i].begin(), max_base);
         }


         m_counts[pos(offset + i)] = m_counts[pos(offset + i)] + newCounts[i]*action;
         m_map[pos(offset + i)] += action * sum_counts_i;

         if (m_map[pos(offset + i)] < 0)
         {
        	 cout<<"problem at position "<<offset<<"+"<<i<<endl;
        	 cout<<"instance counts: "<<newCounts;
        	 cout<<"sequence map at position offset+i: "<<m_map[pos(offset + i)]<<endl;
        	 cout<<"sequence full map: "<<m_map;
         }
         assert(m_map[pos(offset + i)] >= 0);


         if (m_map[pos(offset + i)] == 0)
             m_bases[pos(offset + i)] = nBases;
     }
     if (update)
     {
          m_getRight_dirty = true;
          updateHash(left, right, +1);
     }
}


//  Read Version
void jvSequence::updateReadCounts(jvBaseVec& newCounts, int offset, int action)
{
    uint left = offset, right = offset + newCounts.size();
    pad(left, right);

    bool update = false;

    for (uint i=0; i<newCounts.size(); i++)
    {
    	if (newCounts[i] == (int)nBases) continue;
    	if (m_map[pos(offset + i)] == 0)
    	{
               	update = true;
               	break;
    	}

    	if (m_map[pos(offset + i)] == -action)
    	{
               	update = true;
               	break;
		}
    }


    if (update)
    	updateHash(left, right, -1);

    //  Here we can assert that left, right are non-negative

    for (uint i=0; i<newCounts.size(); i++)
    {
    	if (newCounts[i] == (int)nBases) continue;

        if (m_map[pos(offset + i)] == 0)
            m_bases[pos(offset + i)] = newCounts[i];

        m_counts[pos(offset + i)][newCounts[i]] += action;
        m_map[pos(offset + i)] += action;

        if (m_map[pos(offset + i)] == 0)
            m_bases[pos(offset + i)] = nBases;
    }

    if (update)
    {
        m_getRight_dirty = true;
    	updateHash(left, right, +1);
    }
}

//  pads the sequence so that offsets left, left+1...right-1 are accessible.
//  returns true if there was no change.
void jvSequence::pad(int left, int right)
{

    right = pos(right);
    left = pos(left);

    if ( right > (int)m_bases.size())
    {
      right = max(right, (int)(2*(int)m_bases.size()));
        // postfix with don't-cares
        m_bases.resize(right, nBases);
        m_counts.resize(right, emptyBaseCount);
        m_map.resize(right,0);
        m_ins2.resize(right, emptyInsSet);

        m_right = right + m_left;
    }

    if ( left < 0)
    {
    	left = min(left, -(int)m_bases.size());
        // prefix with -left don't-cares
        m_bases.insert(m_bases.begin(), -left, nBases);
        m_counts.insert(m_counts.begin(), -left, emptyBaseCount);
        m_map.insert(m_map.begin(), -left, 0);
        m_ins2.insert(m_ins2.begin(), -left, emptyInsSet);

        m_left += left;
    }
}


void jvSequence::sampleBases(bool greedy_override)
{
    bool old_map = jvSampler::getMAP();    
    if (greedy_override)
        jvSampler::setMAP(true);

    m_problematic.clear();
    m_problematic.reserve(20);
    m_getRight_dirty = true;
    jvSampler::sampleFromCounts(m_counts, m_bases, m_map, m_left, &m_problematic);

    if (greedy_override)
        jvSampler::setMAP(old_map);
    
}

void jvSequence::printCounts()
{
    int left = getLeft();
    int right = getRight();
    for (uint j=0; j<nBases; j++)
    {
        for (int i=left; i<right; i++)
            cout<<m_counts[pos(i)][j]<<" ";
        cout<<endl;
    }
}

void jvSequence::print(PRINT_OPTIONS opt)
{
	int wdth = 1;
    if (opt & PRINT_BASES)
    {
        // print bases
        cout.precision(2);

        int left = getLeft();
        int right = getRight();
        cout<<m_p<<" | "<<left<<":";
        for (int i=left; i<right; i++)
        {
        	cout.width(wdth);
            if (m_bases[pos(i)] == (int)nBases)
                cout<<".";
            else
                cout<<(int)m_bases[pos(i)]+1;
        }
        cout<<":"<<getRight()<<endl;
        if (opt & PRINT_INS)
                  cout<<endl;
    }
    if (opt & PRINT_SUMMARY)
    {
    	cout.precision(2);
    	int left = getLeft();
    	int right = getRight();
        if (m_editSummary.size()>0)
        {
        	cout<<m_p<<" | "<<left<<":";
        	for (int i=left; i<right; i++)
        	{
        		cout.width(wdth);
        		cout << (int)m_editSummary[pos(i)].in;
        	}
        	cout<<":"<<getRight()<<endl;
        	cout<<m_p<<" | "<<left<<":";
        	for (int i=left; i<right; i++)
        	{
        		cout.width(wdth);
        		cout << (int)m_editSummary[pos(i)].dl;
        	}
        	cout<<":"<<getRight()<<endl;
        }

        cout<<endl;

    }

    // print instances
    if (opt & PRINT_INS )
    {
    	int left = max(m_left, getLeft()-(int)jvInstance::WIDTH);
		int right = getRight();
		cout<<m_M<<" | ";
		for (int i=left; i<right; i++)
		{
			uint S = m_ins2[pos(i)].size();
			if (S>0)
			{
				cout<<i<<":";
				cout<<S<<"(";
				for (set<jvInstance*>::iterator it=m_ins2[pos(i)].begin(); it != m_ins2[pos(i)].end(); it++)
					(*it)->print();
				cout<<") ";
			}
		}
		cout<<" | "<<m_sumO<<endl;
    }

    if (opt & PRINT_DEBUG)
    {
        cout<<m_left<<" "<<m_map;
    }
}


void jvState::printHash()
{
    assert(false);
	for (jvHash::iterator it=m_hash.begin(); it != m_hash.end(); it++)
	{
		if ((*it).empty()) continue;
		for (set<jvLocation>::iterator i=(*it).begin(); i != (*it).end(); i++)
			cout<<"("<<(*i).first<<", "<<(*i).second<<") ";
		cout<<endl;
	}
}


void jvState::checkHash()
{
    assert(false);
	for (jvHash::iterator it=m_hash.begin(); it != m_hash.end(); it++)
	{
		for (set<jvLocation>::iterator i=(*it).begin(); i != (*it).end(); i++)
		{

			int pos = i->second;
			jvSequence* seq = i->first;
			if (!((int)pos>=(int)seq->getLeft() - 100 -1 && (int)pos <= (int)seq->getRight() + 100 + 1))
			{
				cout << (int)pos << endl;
				cout << (int)seq->getLeft() << endl;
				cout << (int)seq->getRight() << endl;
				throw exception();
			}
		}
	}
}


void jvState::buildHash()
{
    clock_t toc, tic = clock();

    fill(m_hash.begin(), m_hash.end(), emptyHashEntry);
    
    for (uint i=0; i<m_seq.size(); i++)
        m_seq[i]->buildHash(&m_hash);

    toc = clock();
    cout<<"buildHash done in "<<(double)(toc-tic)/CLOCKS_PER_SEC<<" seconds.\n";
 
}

void jvSequence::updateHash(int left, int right, int action)
{
    if (!m_hash) return;

    uint K = jvState::K;

    left = max(getLeft(), left-(int)K+1);
    right = min(getRight(), right+(int)K-1);

    jvKmer mask = (1<<2*K)-1; // only works for nBases = 4
    jvKmer stat = 0;
    jvKmer Kmer = 0;

    for (int i=left; i<right; i++)
    {
    	jvBase base = (jvBase)nBases;

    	if ( (pos(i) >=0) && (pos(i) < (int)m_bases.size()))
    		base = m_bases[pos(i)];

    	stat = stat * nBases + (base == (int)nBases? 0 : nBases-1);
    	Kmer = Kmer * nBases + base;
    	stat &= mask;
    	Kmer &= mask;
    	if (stat == mask)
    		if (action == 1)
    			(*m_hash)[Kmer].insert(jvLocation(this, i-K+1));
    		else
    			(*m_hash)[Kmer].erase(jvLocation(this, i-K+1));
    }



}


bool isObserved(int x){return x!=(int)nBases;}

int jvSequence::getLeft()
{
    int i = m_left;

    for (;i<m_right; i++)
        if (m_map[pos(i)] != 0) break;
    if (i==m_right) return 0;
    return i;
}

// returns one slot after the last place
int jvSequence::getRight()
{
    if (!m_getRight_dirty)
    {
        assert( (pos(m_getRight) == m_map.size()) || (m_map[pos(m_getRight)] == 0));
        assert( (m_getRight == 1) || (m_map[pos(m_getRight-1)] > 0));        
        return m_getRight;
    }
    
    int i = m_right-1;
    
    for (; i>m_left; i--)
        if (m_map[pos(i)] != 0) break;
        
    if (i==m_left) m_getRight = 1;
    else m_getRight = i+1;

    m_getRight_dirty = false;
    return m_getRight;    
}

void jvSequence::setFreeze(bool freeze)
{
    m_ins.clear();
    m_isInstanceFrozen = freeze;

    if (freeze)
    {
        for (uint i=0;i<m_ins2.size();i++)
            if (m_ins2[i].size())
                for (set<jvInstance*>::iterator insit = m_ins2[i].begin(); insit != m_ins2[i].end();insit++)
                {
                    (*insit)->reg();
                }
    }
}


int jvSequence::shift(int shiftOffset)
{
    assert(m_isInstanceFrozen);

    pad(-shiftOffset, -shiftOffset+1);

    //  remove assertions later...
    assert(m_left + shiftOffset <= 0);
    assert(m_right + shiftOffset > 0);

    for (uint i=0; i<m_ins.size(); i++)
    {
        m_sumO -= abs(m_ins[i]->m_offset);
        m_ins[i]->m_offset += shiftOffset;
        m_sumO += abs(m_ins[i]->m_offset);
    }

    m_left += shiftOffset;
    m_right += shiftOffset;
    m_getRight += shiftOffset;
    if (!constantP) updateP();
    return -shiftOffset;
}



int jvSequence::center()
{
    assert(m_isInstanceFrozen);
    if (m_ins.empty()) return 0;
    int left = m_ins.front()->getOffset();
    int right = m_ins.back()->getOffset();
    int new_zero = (int) ceil(0.5*(left + right));
    shift(-new_zero);
    return new_zero;
}

double jvSequence::logPo(int offset, bool tail)
{
    if (!jvState::with_geometric)
	return 0; 
    return (tail ? 0 : m_logP) - (offset ? logTwo : 0) + abs(offset)*m_logInvP;
}


jvLikelihood jvSequence::ll()
{
    m_ll.o = 0; 
    if (jvState::with_geometric)
    	m_ll.o = m_sumO*m_logInvP + m_M*(m_logP-logTwo);
    int H = count_if(&m_bases[pos(getLeft())], &m_bases[pos(getRight())], isObserved);
    m_ll.w = H*logBasePrior;  
    m_ll.y = jvSampler::logPy_specific_lr(m_counts, m_bases, (uint)pos(getLeft()), (uint)pos(getRight()));
    m_ll.a = errorModel.logIndel()*(m_insCount + m_delCount);

    return m_ll;
}



bool jvState::splitAndMerge(jvState* st, bool out_print)
{
    assert(false);
    assert(st->m_seq.back()->getM() == 0);


    if (out_print)
        cout<<"\n----------------   SPLIT AND MERGE  -----------------\n";
    uint t1, t2;

    if (st->candI >-1 && st->candJ>-1)
    {
    	t1 = st->candI;
    	t2 = st->candJ;
    }
    else
    {

      // (make sure they are not located exactly at the same offset)
      int trials = 0;
      t1 = jvSampler::intrnd(st->m_ins.size());
      do
      {
	  t2 = jvSampler::intrnd(st->m_ins.size());
	  if (++trials == 5) return false;
      } while ( (st->m_ins[t1]->getOffset() == st->m_ins[t2]->getOffset()) &&
	        (st->m_ins[t1]->getSeq() == st->m_ins[t2]->getSeq()) );
    }

    bool isSplit = (st->m_ins[t1]->getSeq() == st->m_ins[t2]->getSeq());

    if (out_print)
    {
        cout<<"t1="<<t1<<" t2="<<t2<<endl;
        cout<<"proposing "<<(isSplit ? "split":"merge")<<"... ";
        if (!isSplit)
        {
            cout<<"Current:  split sequences are\n";
            st->m_ins[t1]->getSeq()->print();
            st->m_ins[t2]->getSeq()->print();
        }
	else
        {
            cout<<"Current:  merged sequence is\n";
            st->m_ins[t1]->getSeq()->print();
        }
        st->print(PRINT_LL);
    }

    jvSnmUndo undo1 = st->snm(t1,t2);
    double p_nst = undo1.p_nst;
    double l_nst = st->ll().ll();
    if (out_print)
    {
        if (!undo1.isSplit)
        {
            cout<<"merged sequence is\n";
            st->m_ins[t1]->getSeq()->print();
        }
	else
        {
            cout<<"split sequences are\n";
            st->m_ins[t1]->getSeq()->print();
            st->m_ins[t2]->getSeq()->print();
        }
        st->print(PRINT_LL);
        cout<<"forward transition prob = "<<p_nst<<", l_nst="<<l_nst<<endl;
    }

    assert(st->m_seq.back()->getM() == 0);

    jvSnmUndo undo2 = st->snm(t1, t2, true);
    double p_0 = undo2.p_nst;
    double l_0 = st->ll().ll();
    if (out_print)
    {
        if (!undo2.isSplit)
        {
            cout<<"backward: merged sequence is\n";
            st->m_ins[t1]->getSeq()->print();
        }
        else
        {
            cout<<"backward: split sequences are\n";
            st->m_ins[t1]->getSeq()->print();
            st->m_ins[t2]->getSeq()->print();
        }
	st->print(PRINT_LL);
    	cout<<"backward transition prob = "<<p_0<<", l_0="<<l_0<<endl;
    }
    undo2.undo();

    assert(st->m_seq.back()->getM() == 0);

    double accept = min(1.0, exp(l_nst-l_0)*p_0/p_nst);

    bool accepted = true;

    if (jvSampler::unirnd()>accept)
    {
    	undo1.undo();
    	accepted = false;
    }

    if (out_print)
    {
        cout<<"Acceptance probability:  "<<accept<<endl;
        cout<<"----------------->    ";
        cout<<( accepted ? "ACCEPTED!" : "REJECTED!");
        cout<<"    <-----------------\n\n";
    }

    assert(st->m_seq.back()->getM() == 0);

    st->cleanSequences();
    return accepted;

}




jvSnmUndo jvState::snm(uint t1, uint t2, bool changeState)
{
    assert(false);

    // Save original ll
    // Also, we need to call ll() to have the cached m_ll computation.
    jvLikelihood start_ll = ll();

    jvSequence* s  = m_ins[t1]->getSeq();
    jvSequence* s2 = m_ins[t2]->getSeq();
    int o1 = m_ins[t1]->getOffset();
    int o2 = m_ins[t2]->getOffset();

    set<int> set_offsets;
    vector<int> offsets;
    vector<double> lls;

    //  if the instances are from the same sequence, split it at some
    //  point that will separate t1 and t2.  sample from all possible such splits.
    if (s == s2)
    {
        assert(o1 != o2);
        if (o2<o1) {swap(o1,o2); swap(t1,t2);}

        vector<jvInstance*>::iterator it = find(s->m_ins.begin(), s->m_ins.end(), m_ins[t1]);
        for (; (it != s->m_ins.end()) && ((*it)->getOffset()<=o2); it++)
            set_offsets.insert((*it)->getOffset());
        set_offsets.erase(o1);

        for (set<int>::iterator ssp = set_offsets.begin(); ssp != set_offsets.end(); ssp++)
            offsets.push_back(*ssp);

        newSeq();
        s2 = m_seq[m_seq.size()-2];
        if (s2->getM() != 0)
            assert(false);

        for (uint i=0; i<offsets.size(); i++)
        {
        	jvLikelihood delta_ll = s->splitLl(offsets[i],m_gamma);
        	delta_ll = delta_ll + start_ll;
            lls.push_back(delta_ll.ll());
        }

        lls = lls;
        int o=jvSampler::lmnrnd(lls);

        jvSnmUndo undo;
        if (changeState)
        {
        	undo = s->split(offsets[o], true, s2);
        	undo.back1 = s->center();
        	undo.back2 = s2->center();
        }
        undo.l_nst = lls[o];
        undo.p_nst = exp(lls[o]-log_sum_exp(lls));
        return undo;
    }

    //  else (instances from different sequences), merge the sequences so that
    //  one is completely after the other.  sample from all possible such merges.
    else
    {

        //  merge sites from the left
        int s2_right = s2->getRight();
        int slack = s2_right - s2->m_ins.back()->getOffset(); // bases from last instances
        for (int o=s->m_ins.front()->getOffset()-slack; o<=s->m_ins.front()->getOffset()-1; o++)
            offsets.push_back(o);

        //  merge sites from the right
        int s_right = s->getRight();
        for (int o=s->m_ins.back()->getOffset()+1; o<=s_right; o++)
            offsets.push_back(o);


//        cout<<"offsets: "<<offsets;
        double p_merged = 0, p;
        assert(offsets.size() >= 3);

        //  this part can be removed later!
        //  It has no effect on the proposal dist.
        jvDVec mergedM(1), splitM(2);
        splitM[0] = s->getM();  splitM[1] = s2->getM();
        mergedM[0] = sum(splitM);
        double deltaLogPs = jvSampler::DPll(mergedM, m_gamma)
        - jvSampler::DPll(splitM, m_gamma);
        start_ll.s += deltaLogPs;
        ////////////////////////////////////

        for (uint i=0; i<offsets.size(); i++)
        {
        	jvLikelihood delta_ll = s->mergedLl(s2, offsets[i] ,&p_merged);
        	delta_ll = delta_ll + start_ll;
        	lls.push_back(delta_ll.ll());
        }

        //  Sample a merge in the tail (using the approximate last p_merged).
        p = 1-pow(1-p_merged, s2->m_ins.size());
        uint K = 5;
        offsets.front() -= K;
        offsets.back() += K;
        lls.front() -= log(p);  // increase ll for tails
        lls.back() -= log(p);

        lls = lls;
        cout.precision(5);

        // sample the merge
        int o=jvSampler::lmnrnd(lls);

        jvSnmUndo undo;
        if (changeState)
        {
        	undo = s->merge(s2, offsets[o]);
        	undo.back1 = s->center();
        }
        undo.l_nst = lls[o];
        undo.p_nst = exp(lls[o]-log_sum_exp(lls)); // winning prob. 
        return undo;

    }
}


void jvSnmUndo::undo()
{
	if (isSplit)
	{
        s->shift(back1);
        s2->shift(back2);
		s->merge(s2, offset);
	}
	else
	{
        s->shift(back1);  // from centering
        s->split(offset, rightMoves, s2);
        s2->shift(back2); // from aligning
        assert((bases.size() == 0) || (ov_right-ov_left == (int)bases.size()));
        for(int o=ov_left; o<ov_right; o++)
        	s2->m_bases[s2->pos(o)] = bases[o-ov_left];
	}

}





bool jvSequence::getOverlap(jvSequence* seq, int offset,
		int& shift_offset, int& ov_left, int& ov_right)
{
	if (m_ins.empty())
	{
		ov_left = seq->m_ins.front()->getOffset(); ov_right = ov_left;
		shift_offset = -ov_left + offset;
		return true;
	}


    int left = m_ins.front()->getOffset();
    int right = m_ins.back()->getOffset();

    if (offset > right)  // putting seq on the right
    {
        shift_offset = -seq->m_ins.front()->getOffset() + offset;
        ov_left = offset-shift_offset;  //assert(ov_left == seq->getLeft());
    }
    else if (offset < left) // putting seq on the left
    {
        shift_offset = -seq->m_ins.back()->getOffset() + offset;
        ov_left = left-shift_offset;
    }
    else {cout << "Sequence is "<<left<<" - "<<right<<", Can't merge into "<<offset<<"\n"; assert(false);}

    ov_right = min(getRight()-shift_offset, seq->getRight());
    ov_right = max(ov_right, ov_left);

    assert((ov_left == ov_right) || (seq->pos(ov_right)<= (int)seq->m_counts.size()));
	return (offset > right);
}







jvLikelihood jvSequence::mergedLl(jvSequence* seq, int offset, double* p_merged)
{
    //  find overlapping indices in both sequences
    //  store them in ov_src and ov_tgt
    assert(m_isInstanceFrozen);

    jvLikelihood ll;

    //  Calculate shift_offset - how much the offsets in seq should move.
    int shift_offset;
    int ov_left, ov_right;
    bool mergeToRight = getOverlap(seq, offset, shift_offset, ov_left, ov_right);

    //  log_Po
    vector<int> merged;
    merged.reserve(m_M + seq->m_M);
    if (mergeToRight)  // putting seq on the right
    {
        for (uint i=0; i<m_ins.size(); i++)
        	merged.push_back(m_ins[i]->getOffset());
        for (uint i=0; i<seq->m_ins.size(); i++)
        	merged.push_back(seq->m_ins[i]->getOffset()+shift_offset);
    }
    else // putting seq on the left
    {
        for (uint i=0; i<seq->m_ins.size(); i++)
        	merged.push_back(seq->m_ins[i]->getOffset()+shift_offset);
        for (uint i=0; i<m_ins.size(); i++)
        	merged.push_back(m_ins[i]->getOffset());
    }
    ll.o = -m_ll.o -seq->m_ll.o + ll_o(merged, p_merged);

    //  Log_Py

    //  prepare counts for overlapped area in seq
//    assert(seq->pos(ov_right)<<" "<<seq->m_counts.size()<<endl;
    jvBaseCountsVec ovCounts(&seq->m_counts[seq->pos(ov_left)],
    		&seq->m_counts[seq->pos(ov_right)]);


    //  subtract likelihood from its own, add likelihood with
    //  respect to this seq

    ll.y  =    -jvSampler::logPy_specific(ovCounts, seq->m_bases, seq->pos(ov_left))
          	+ jvSampler::logPy_specific(ovCounts, m_bases, pos(ov_left+shift_offset));

    if (false)
    {
        double d1 = jvSampler::logPy_specific(ovCounts, seq->m_bases, seq->pos(ov_left));
               ll.y=    -jvSampler::logPy_specific(ovCounts, seq->m_bases, seq->pos(ov_left));

        double d2 = jvSampler::logPy_specific(ovCounts, m_bases, pos(ov_left+shift_offset));

        cout<<"seq: "<<d1<<"  this: "<<d2<<endl;

            for (int a=0; a<ov_right-ov_left; a++)
            cout<<(int)seq->m_bases[seq->pos(ov_left+a)];
        cout<<endl;
        for (int a=0; a<ov_right-ov_left; a++)
            cout<<(int)m_bases[pos(shift_offset+ov_left+a)];
        cout<<endl;

        jvBaseVec  ins2read;
    	//  get "consensus read" from instance counts
        uint L = ovCounts.size();
        ins2read.resize(L, nBases);
        for (uint i=0; i<L; i++)
        {
            jvBaseCounts& counts = ovCounts[i];
            jvBaseCounts::iterator max_base = max_element(counts.begin(), counts.end());
            if (*max_base == 0) continue;
            ins2read[i] = distance(counts.begin(), max_base);
        }
        for (int a=0; a<ins2read.size(); a++)
            cout<<(int)ins2read[a];
        cout<<endl;
        cout<<ovCounts;

    }

    //  log_Pw
    int H = count_if(&seq->m_bases[seq->pos(ov_left)], &seq->m_bases[seq->pos(ov_right)], isObserved);
    ll.w = -H*logBasePrior;

    return ll;

}


double jvSequence::ll_o(vector<int>& os, double* p_merged)
{
    if (!jvState::with_geometric) 
	return 0;  
    int left = os.front();
    int right = os.back();
    int new_zero = (int) ceil(0.5*(left + right));
    os = os + (-new_zero);

    int O_sum = 0;
    for (uint i=0; i<os.size(); i++)
    	O_sum += abs(os[i]);

    double a = os.size()+pseudoA;
    double b = O_sum + pseudoB;
    double p = a/(a+b);
    if (constantP) p = pseudoA/(pseudoA+pseudoB); 
    if (p_merged) *p_merged = p;

    return O_sum*log(1-p) + os.size()*(log(p) - logTwo);

}



int jvInstance::getLeft()
{
    assert(m_N == 1);
    jvReadAss* r = *(m_reads.begin());
    return r->getOffset() + m_offset;    
}

int jvInstance::getRight()
{
    assert(m_N == 1);
    jvReadAss* r = *(m_reads.begin());
    return r->getOffset() + m_offset + r->getData()->getBases().size();
}


jvLikelihood jvSequence::splitLl(int offset, double gamma)
{

	//  find overlapping indices in both sequences
	//  store them in ov_src and ov_tgt
    assert(m_isInstanceFrozen);

    jvLikelihood ll;

	//  Log_Py is the same after a split
    // - assuming overlapping bases are duplicated.  !!!


    // logPo
    vector<int> os1, os2;
    uint i=0;
    for (; m_ins[i]->getOffset() < offset; i++)
    {
    	os1.push_back(m_ins[i]->getOffset());
    }
    for (uint j=i; j<m_ins.size(); j++)
    	os2.push_back(m_ins[j]->getOffset());

    ll.o = -m_ll.o + ll_o(os1) + ll_o(os2);

    // logPw
    int ov_left = m_ins[0]->getLeft();
    int ov_right = m_ins[0]->getRight();
    for (uint j=1; j<i; j++)
    {
    	ov_right = max(ov_right, m_ins[j]->getRight());
    	ov_left = min(ov_left, m_ins[j]->getLeft());
	}

    int ov_right2 = m_ins[i]->getRight();
    int ov_left2 = m_ins[i]->getLeft();
    for (uint j=i+1; j<m_ins.size(); j++)
    {
    	ov_right2 = max(ov_right2, m_ins[j]->getRight());
    	ov_left2 = min(ov_left2, m_ins[j]->getLeft());
	}

    ov_right = min(ov_right, ov_right2);
    ov_left = max(ov_left, ov_left2);

    int H=0;
    if (ov_left<ov_right)
    	H = count_if(&m_bases[pos(ov_left)], &m_bases[pos(ov_right)], isObserved);
    ll.w = H*logBasePrior;

    // logPs
    jvDVec mergedM(1), splitM(2);
    mergedM[0] = m_M;
    splitM[0] = i;  splitM[1] = m_M-i;
    double deltaLogPs = -jvSampler::DPll(mergedM, gamma)
    + jvSampler::DPll(splitM, gamma);
    ll.s += deltaLogPs;


    return ll;
}




//  likelihood counts given that they were taken from offset o in bases,
//  taking into account only the entries in [left, right).
double /*inline*/ jvSequence::logPy_specific_bb(jvBaseCountsVec& counts, int o)
{
	int left = 0;
	int right = (int)m_bases.size();
    assert(left < right);

    double res = 0;
    for (int i=0; i<(int)counts.size(); i++)
    {
    	//  if position in sequence is unobserved
        if ( (o+i < m_left) || (o+i >= m_right) || (m_bases[pos(o+i)] == (int)nBases))
            res += jvSampler::logPy_base_unobserved(counts[i]);
        else // position in sequence is observed
            res += jvSampler::logPy_base_observed(counts[i], m_bases[pos(o+i)]);
    }
    return res;

}

double /*inline*/ jvSequence::logPy_specific_bb(jvBaseVec& read, int o)
{
    jvBaseCountsVec counts;
    jvSampler::readToCounts(read, counts);
    return logPy_specific_bb(counts, o);
}

// builds the align buffer
// aligning the half-closed-interval [o-extend, o+read_length+extend)
jvBaseVec jvSequence::buildAlignBuff(int L, int o)
{
    assert(false);
	int start = o-jvAlign::extend;
	int end = o + L + 2*jvAlign::extend;
	jvBaseVec buff(end-start);

	for (int i=0; i<(int)buff.size(); i++)
	{
		//  if position in sequence is unobserved
		if ( (start+i < m_left) || (start+i >= m_right) || (m_bases[pos(start+i)] == (int)nBases))
                    buff[i] = nBases;
		else // position in sequence is observed
                    buff[i] = m_bases[pos(start+i)];
	}
	return buff;
}

//  align read to sequence at position o.
//  padding of jvAlign::extend bases is added so that we align [o-extend, o+read_length+2*extend)
//  o is also an out parameter.  It corrects the input o by adjusting it to trailing inserts.
double jvSequence::logPy_align(jvBaseVec& read, int *o, set<jvAlignCache>* cache, bool print)
{
	assert(walign);
	jvBaseVec* inSeq;
	int L  = read.size()+3*jvAlign::extend;
	int start = 0;
	jvBaseVec buff(0);

        inSeq = &m_bases;
        start = pos(*o-jvAlign::extend);

	double score = jvAlign::aligner.align(inSeq, start, start + L, &read, cache, print);

	*o += jvAlign::aligner.getLocalStart() - jvAlign::extend;

    return score;
}


//  align instance to sequence at position o.
//  padding of jvAlign::extend bases is added so that we align [o-extend, o+read_length+extend)
//  each read might change its position slightly.  The new locations for the read are set
//  in the m_new_offset member of each read.
double jvSequence::logPy_align_ins(jvInstance* ins, int o, bool change_reads, bool print)
{
	assert(walign);
    assert(false);
    return 0;
}




void jvState::saveStatus(const char* filename, uint iteration)
{
    summarizeEdits();
    
    int max_len = 0;
    int vals[] = {1,2,3,4,5,10,20,40,80,160,320,640};
    vector<int> hist_vals(vals, vals+12);
    hist_vals.push_back(INT_MAX);
    vector<int> n_reads_hist(hist_vals.size(),0);
    
    // collect stats
    for (int s=0; s<m_seq.size()-1; s++)
    {
        int M = m_seq[s]->getM(); 
        assert(M != 0);
        int len = m_seq[s]->getRight()-m_seq[s]->getLeft();
        if (len>max_len) max_len = len;
        for (int j=0; j<(int)hist_vals.size(); j++)
            if (M<=hist_vals[j])
            {
                n_reads_hist[j]++;
                break;
            }
    }

    clock_t toc = clock();
    int seconds = (toc-m_start)/CLOCKS_PER_SEC;
    m_start = toc;

            
    ofstream myfile;
    jvLikelihood l = ll();
    myfile.open(filename, ios::app);
    myfile << l.ll()<<" "<<l.s<<" "<<l.t;
    myfile <<" "<<l.o<<" "<<l.w<<" "<<l.y<<" "<<l.a;
    myfile<<" "<<max_len<<" "<<m_seq.size()-1;
    for (int i=0; i<n_reads_hist.size(); i++)
        myfile<<" "<<n_reads_hist[i];
    myfile<<" "<<m_total_splits<<" "<<m_total_merges<<" "<<m_total_failures;
    myfile<<" "<<change_seq<<" "<<seconds<<" "<<iteration<<endl;
    myfile.close();

    // reset split-and-merge stats    
    m_total_splits = 0;
    m_total_merges = 0;
    m_total_failures = 0;
    

    
}


int findPropSplit = 0;
int findAccSplit = 0;

bool jvState::findAndSplitEnds(bool all)
{
    bool acc = false;

    pair<jvInstance*, jvInstance*> ret;
    vector<int> perm = jvSampler::randomPermutation(m_seq.size());
    for (vector<int>::iterator pi=perm.begin(); pi<perm.end(); pi++)
    {
    	int i = *pi;

        //  don't split the skeleton sequences
        if (m_seq[i]->m_never_clean)
            continue;

        ret =  m_seq[i]->getGapInstances(all);

	if (ret.first == ret.second) continue;
	assert(ret.first); assert(ret.second);
        if (ret.first->getOffset() == ret.second->getOffset()) continue;

        cout<<"Found a place to split!\n";
        for (int ii=0; ii<m_ins.size(); ii++)
        {
            if ( m_ins[ii] == ret.first)
                candI = ii;
            if ( m_ins[ii] == ret.second)
                candJ = ii;
            if (candI>-1 && candJ>-1)
                break;
        }

        acc = splitAndMerge(this,true);
        candI = -1;
        candJ = -1;

        findPropSplit++;
        findAccSplit+=acc;
        if (acc) cout << "SNM move accepted ";
        else cout << "SNM move rejected ";
        cout << "\t\tAcceptance rate: " << ((double)findAccSplit)/((double)findPropSplit) << endl;
        cout << "So far we saw "<<broken_instance<<" broken instances.\n";
    }

    return acc;
}


void jvState::flipSequences()
{
    if (!both_orientations)
    {
        cout<<"Both orienations flag is false.\n";
        return;
    }
    
    double accept = 0.1;
    for (uint i=0; i<m_seq.size()-1; i++)
        if (jvSampler::unirnd()<accept)
            m_seq[i]->revComplement();

}

// don't remove these, need them to triage memory usage
// if you redefine classes update these
long jvInstance::getMem()
{
   return 0 +   
    sizeof(int) + //int m_N; 
    sizeof(int) + //int m_offset;
    sizeof(jvSequence*) + //jvSequence* m_seq;
    sizeof(bool) + //bool m_isAssigned;
    sizeof(bool) + //bool m_isReadFrozen;
//    m_readCounts.capacity()*sizeof(jvBaseCounts)+//jvBaseCountsVec m_readCounts;
    m_reads.size()*sizeof(jvReadAss*);//set<jvReadAss*> m_reads;
}

long jvSequence::getMem()
{
  return 0 + 
    sizeof(int) + //int m_M;
    sizeof(int) + //int m_left;
    sizeof(int) + //int m_right;
    sizeof(double) + //double m_p;
    sizeof(double) + //double m_logP;
    sizeof(double) + //double m_logInvP;
    sizeof(int) + //int m_sumO;
    sizeof(m_bases) + //jvBaseVec m_bases;
    m_map.capacity()*sizeof(int) + //vector<int> m_map;
    m_editSummary.capacity()*sizeof(jvEdit)+ //vector<jvEdit> m_editSummary;
    sizeof(int) + //int m_insCount;
    sizeof(int)+ //int m_delCount;
    m_counts.capacity()*sizeof(jvBaseCounts) + //jvBaseCountsVec m_counts;
    sizeof(bool) + //bool m_isInstanceFrozen;
    m_ins.capacity()*sizeof(jvInstance*) + //vector<jvInstance*> m_ins;
    sizeof(jvLikelihood) + //jvLikelihood m_ll;
    m_logPy.capacity()*sizeof(double) +//jvDVec m_logPy;  // byproduct of running sampleO
    m_spikes.capacity()*sizeof(int) + //vector<int> m_spikes; // byproduct of running sampleInstance
    m_spikes_bu.capacity()*sizeof(int) + //vector<int> m_spikes_bu; // byproduct of running sampleInstance
    m_problematic.capacity()*sizeof(int) + ///vector<int> m_problematic;
    sizeof(bool) + //bool m_never_clean;
    sizeof(bool) + // bool m_inprocess;  // used in jvState::computeSpikes to tell if it was hit.
    sizeof(jvHash*); //    jvHash* m_hash;

}

long jvReadAss::getMem()
{
  return 0 + 
    sizeof(jvReadData*) +//jvReadData* m_data;
    sizeof(jvInstance*) +//jvInstance* m_ins;
    sizeof(int) + //int m_offset;
    sizeof(int) + //int m_new_offset;  // this member is used by sampleInstance
    sizeof(bool); bool m_isAssigned;
    
}

long jvState::getMem()
{
  int seqTot = 0;
  int insTot = 0;
  int raTot = 0;
  int dataTot = 0;
  int hashTot = 0;
  
  for (int i=0;i < m_seq.size();i++) { seqTot = seqTot + m_seq[i]->getMem(); }  ; cout << "Seq Mem: " << seqTot << endl;
  for (int i=0;i < m_ins.size();i++) { insTot = insTot + m_ins[i]->getMem(); } ;  cout << "Ins Mem: " << insTot << endl;
  for (int i=0;i < m_reads.size();i++) {raTot = raTot + m_reads[i].getMem(); } ; cout << "RA Mem: " << raTot << endl;
  for (int i=0;i < m_data.size();i++) {dataTot = dataTot + m_data[i].getMem(); } ; cout << "Data Mem: " << dataTot << endl;
  hashTot += m_hash.capacity()*sizeof(set<jvLocation>);
  for (int i=0;i < m_hash.size(); i++) {hashTot += m_hash[i].size() * sizeof(jvLocation);} cout<< "Hash Mem: "<< hashTot << endl;
}


void jvState::unassignEnds()
{
    freezeInstances(true);
    // go over sequence
    for (int s=0; s<(int)m_seq.size()-1; s++)
    {
    // if sequence has more than X reads
        jvSequence* seq = m_seq[s];
        if (seq->getM() > 10)
        {            
            assert(seq->m_ins.front()->m_reads.size() == 1);
            //unassign first read
            jvReadAss* r = *(seq->m_ins.front()->m_reads.begin());
            r->assign(false);
            r->setIns(NULL);
            
            //unassign last (ending) read
            assert(seq->m_ins.back()->m_reads.size() == 1);
            r = *(seq->m_ins.back()->m_reads.begin());
            r->assign(false); 
            r->setIns(NULL);
       }
    }
    freezeInstances(false);
    cleanInstances();    
}

bool jvState::sanityGetRight()
{
    for (int s=0; s<m_seq.size(); s++)
    {
        m_seq[s]->getRight();
    }
    return true;
}


bool jvState::sanityEditString()
{
    int rd = 0;
    
    // Go over reads
    for (int i=0; i<(int)m_reads.size(); i++)
    {
        vector<char>* edt = m_reads[i].getData()->getEdit();
        // if the edit string is all zeros assert that rawBases is the same as m_bases.
        bool gotOne = true;
        for (int c=0; gotOne && c<edt->size(); c++)
            if ( (*edt)[c] != -1 ) gotOne = false;        

        if (gotOne)
        {
            rd++;
            jvBaseVec* v1 = m_reads[i].getData()->getRawBases();
            jvBaseVec* v2 = &m_reads[i].getData()->getBases();
            if(v1->size() != v2->size())
            {
                cout<<"read "<<i<<" has a problem.\n";
                print_bases(v1);
                print_bases(v2);
                print_edit(*edt);
                return false;

            }
            for (int c=0; c<v1->size(); c++)
                if( (*v1)[c] != (*v2)[c])
                {
                    cout<<"read "<<i<<" has a problem.\n";
                    print_bases(v1);
                    print_bases(v2);
                    print_edit(*edt);
                    return false;
                }
        }
    }
    cout<<rd<<" clean reads!\n";
    return true;
}
  

//  returns the first n reads out of a list, and then remove them from the list.
vector<jvReadAss*> jvState::popReads(int n)
{
    vector<jvReadAss*> popped(0);
    if (n == 0) return popped;

    //  pop the reads from the stuck and populate the 'popped' array
    while( popped.size() < n )
    {
        if (m_read_stack.size() == 0)
        {

            /////  Rebuild stack    /////
            vector<int> sigma(m_reads.size());
            vector< pair<double, int> > aln(m_reads.size());
            int cap1 = 0;
            for (int i=0; i<m_reads.size(); i++)
            {
                double score = m_reads[i].getData()->getScore();
                aln[i].first = score;
                aln[i].second = i;
                if (score < 3*errorModel.logNoise()) cap1++;
            }        
            
            cout<<"popReads: Total "<<cap1<<" reads with score < "<< 3*errorModel.logNoise()<<".\n";
            int cap2 = min(cap1, (int)aln.size());

            partial_sort(aln.begin(), aln.begin()+cap1, aln.end());
            
            sigma = jvSampler::randomPermutation(cap1);        
            m_read_stack.resize(cap2);
            for (int i=0; i<cap2; i++)
                m_read_stack[i] = aln[sigma[i]].second;
            /////////////////////////////
            
        }
        
        jvReadAss* read = &m_reads[m_read_stack.back()];
        m_read_stack.pop_back();
        if (read->getData()->getScore() < -10.0)
            popped.push_back(read);
    }
    return popped;    
}


bool jvState::splitAndMerge3(bool out_print)
{
    assert(m_seq.back()->getM() == 0);
    summarizeEdits();

    cutoff = -1;
    
    if (out_print)
        cout<<"\n----------------   SPLIT AND MERGE  -----------------\n";

    // choose x0 != x1:
    uint x0, x1;

    if (candI >-1 && candJ>-1)
    {
    	x0 = candI;
    	x1 = candJ;
        assert(max(x0,x1)<m_reads.size());
    }
    else
    {
    	int trials = 0;
        x0 = jvSampler::intrnd(m_reads.size());
	do
	{
            if (++trials == 10) break;
            x1 = jvSampler::intrnd(m_reads.size());           
	} while ( (x0 == x1) || (m_reads[x0].getIns()->getSeq()->getBin() != m_reads[x1].getIns()->getSeq()->getBin()) );
        
        if (trials == 10)  // strategy #2 (a little more computationally expensive)
        {
            trials = -1;
            vector<int> indexes = jvSampler::randomPermutation(m_reads.size());
            do
            {
                if (++trials == m_reads.size()) return false;
                x1 = indexes[trials];
                
            } while ( (x0 == x1) || (m_reads[x0].getIns()->getSeq()->getBin() != m_reads[x1].getIns()->getSeq()->getBin()) );                        
        }  
    }
    
    vector<jvReadAss*> x(2);
    x[0] = &m_reads[x0];
    x[1] = &m_reads[x1];

    jvSequence* seq0 = x[0]->getIns()->getSeq();
    jvSequence* seq1 = x[1]->getIns()->getSeq();

    bool isSplit = (seq0 == seq1);

    if (seq0->getBin() != seq1->getBin()) return false;

    if (out_print)
    {
        cout<<"x0="<<x0<<" x1="<<x1<<endl;
        cout<<"proposing "<<(isSplit ? "split":"merge")<<"... ";
        if (!isSplit)
        {
            cout<<"Current:  split sequences are\n";
            seq0->print(PRINT_BASES);
            seq1->print(PRINT_BASES);
        }
	else
        {
            cout<<"Current:  merged sequence is\n";
            seq0->print(PRINT_BASES);
        }
        print(PRINT_LL);
    }

   
    // save state and compute current likelihood p(x)
    jvLikelihood likelihood = ll();
    double ll_old = likelihood.ll()-likelihood.w;
    
    // assign reads x_0 and x_1 to two different variants
    gibbsReads(x, (GIBBS_POLICY)(GIBBS_INCLUDE | GIBBS_NEW_SEQS_ONLY) );
    jvSequence* new_seq0 = m_seq[m_seq.size()-3];
    jvSequence* new_seq1 = m_seq[m_seq.size()-2];
    new_seq0->m_never_clean = true;
    new_seq1->m_never_clean = true;
    assert(x[0]->getIns()->getSeq() == new_seq0);
    assert(x[1]->getIns()->getSeq() == new_seq1);
    
    //  build a list of reads:
    vector<jvReadAss*> reads = seq0->getReads();
    int n_seq0 = reads.size();
    if (!isSplit) // seq0 != seq1
    {
        vector<jvReadAss*> reads_ = seq1->getReads();
        //  append reads_ to reads
        reads.insert(reads.end(), reads_.begin(), reads_.end());
    }


    // run k restricted gibbs iterations on reads-(x_1,x_2)
    for (int k=0; k<2; k++)
    {
        gibbsReads(reads, (GIBBS_POLICY)(GIBBS_INCLUDE | GIBBS_SNM_ONLY) );
        new_seq0->sampleBases(true); // greedy 
        new_seq1->sampleBases(true); // greedy
    }

    if (out_print) cout<<"--->  At launch state!  <--- \n";

    print(PRINT_LL);
        
    
    //   run one restricted gibbs iteraiton on reads-(x_1,x_2).
    //   in a split->merge move force results to agree with original split
    //   report transition probability as p(x->x')
    cutoff = (isSplit ? -1 : n_seq0); // hidden parameter to gibbsReads
    jvLogProb p_m2s = gibbsReads(reads, (GIBBS_POLICY)(GIBBS_INCLUDE | GIBBS_SNM_ONLY) );
    
    if (!isSplit)  //going for a merge
    {        
    //   merge reads, run a gibbs iteration or two, compute likelihood p(x')
        cout<<"Merging...\n";
        reads.push_back(x[1]);
        cutoff = reads.size();
        for (int k=0; k<1; k++)
        {
            gibbsReads(reads, (GIBBS_POLICY)(GIBBS_INCLUDE | GIBBS_SNM_ONLY) );
            new_seq0->sampleBases();
        }

    }

    if (out_print) cout<<"--->  At new state!  <--- \n";
    summarizeEdits();
    print(PRINT_LL);
    likelihood = ll();
    double ll_new = likelihood.ll()-likelihood.w;
    
    //  compute acceptance prob p(x')p(x'->x)/[p(x)p(x->x')]
    cout<<"Merge->split transition: "<<p_m2s.choice-p_m2s.partition<<endl;
    cout<<"ll_new - ll_old: "<<ll_new -ll_old<<endl;
    double accept = exp(ll_new -ll_old + (isSplit? -1 :1)*(p_m2s.choice-p_m2s.partition));
    
    //  decide:  accept or reject.
    //  if reject - undo
    
    bool accepted = true;

    if (jvSampler::unirnd()>accept)
    	accepted = false;

    if (out_print)
    {
        cout<<"Acceptance probability:  "<<accept<<endl;
        cout<<"----------------->    ";
        cout<<( accepted ? "ACCEPTED!" : "REJECTED!");
        cout<<"    <-----------------\n\n";
    }
    
    if (!accepted) // undo if rejected
    {
        if (isSplit)
        {
            //   merge reads
            reads.push_back(x[1]);
            cutoff = reads.size();
        }            
        else
        {
            cutoff = n_seq0; 
        }
        for (int k=0; k<1; k++)
        {
            gibbsReads(reads, (GIBBS_POLICY)(GIBBS_INCLUDE | GIBBS_SNM_ONLY) );
            assert(isSplit == (new_seq1->getM() == 0));
            new_seq0->sampleBases();
            if (!isSplit)
                new_seq1->sampleBases();
        }
    }
        
    new_seq0->m_never_clean = false;
    new_seq1->m_never_clean = false;
    cleanSequences();
    cutoff = -1; 

    // split-and-merge stats:
    if (!accepted)
        m_total_failures++;
    else
        if (isSplit)
            m_total_splits++;
        else
            m_total_merges++;
    
    return accepted;

}

vector<jvReadAss*> jvSequence::getReads()
{
   vector<jvReadAss*> reads;
   for (uint i=0;i<m_ins2.size();i++)
            if (m_ins2[i].size())
                for (set<jvInstance*>::iterator insit = m_ins2[i].begin(); insit != m_ins2[i].end();insit++)
	            for (set<jvReadAss*>::iterator it = (*insit)->m_reads.begin();it != (*insit)->m_reads.end();it++)
              		reads.push_back(*it);
   return reads;

}

