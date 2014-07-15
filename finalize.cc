#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <assert.h>
#include <algorithm>

using namespace std;

#include "jvDataGenerator.h"
#include "jvState.h"
#include "jvStateSerializer.h"
#include "jvAlign.h"
#include "jvAlign2.h"
#include "jvFastaParser.h"


bool comp_length(const jvReadData& rd1, const jvReadData& rd2)
{
    return rd1.getBasesLength() >= rd2.getBasesLength();
}


bool split_seqs(vector<jvReadData>& seqs)
{
    int gap_thresh = 3;
    int gap_length = -1000000;
    for (int i=0; i<seqs.size(); i++)
        for (int j=0; j<seqs[i].getBases().size(); j++)           
            if (seqs[i].getBases()[j] != (int)nBases)
            {
                if (gap_length > gap_thresh)
                {
                    // found a gap.  the letter j-gap_length is where the
                    // gap starts.  the letter j-1 is where the gap ends.
                    jvBaseVec seq_right(&seqs[i].getBases()[j],
                                        &seqs[i].getBases()[seqs[i].getBasesLength()]);
                    jvBaseVec seq_left (&seqs[i].getBases()[0], &seqs[i].getBases()[j-gap_length]); 
                    seqs[i] = seq_left;
                    seqs.push_back(seq_right);
                    return true;
                }
                gap_length = 0;
            }
            else
                gap_length = gap_length+1;
    return false;
}
    
    
    
int main(int argc, char ** argv)
{
    if (argc < 4)
    {
        cout<<"will output to <fasta_file> all the contig sequences in <dump_file> that have length greater than <cutoff>.\n";
        cout<<"Example: \n";
        cout<<"finalize 0 output.fa all_reads.fa.dump.best\n";
        return 1;
    }
    
    string fasta_file = string(argv[2]);
    int B = atoi(argv[1]);
 
    vector<jvReadData> data;
    vector<jvReadData> output_seqs;

    for (int s=3; s<argc; s++)
    {
        vector<jvReadData> seqs;
        string state_file = string(argv[s]);
    
        cout<<"Loading State "<<state_file<<"..."<<flush;
        jvState* jvs = jvStateSerializer::read(state_file.c_str(), data);
        cout<<" Done!"<<endl;

        for (vector<jvSequence*>::iterator it=jvs->getSeq().begin(); it != jvs->getSeq().end()-1; it++)
        {
            int left = (*it)->getLeft();
            int right = (*it)->getRight();
            
            if (right - left < B) continue;
            if ((*it)->getM() <= 1) continue;

            jvBaseVec seq(&(*it)->getBases(left), &(*it)->getBases(right));
            seqs.push_back(jvReadData(seq));        
        }

        // split sequence when there is a large gap (>=3)
        while (split_seqs(seqs));
        
        // erase sequences with <=B bases
        for (int k=seqs.size()-1; k>=0; k--)
        {
            if (seqs[k].getBases().size() <= B)
                seqs.erase(seqs.begin()+k);
        }
                    
        
        // align and remove duplicates
        vector<int> thrown;
        for (int k=0; k<seqs.size(); k++)
        {
            if (find(thrown.begin(), thrown.end(), k) == thrown.end())
                output_seqs.push_back(seqs[k]);
        }
        delete jvs;

    }
    
    sort(output_seqs.begin(), output_seqs.end(), comp_length);    
    vector<jvStringReadData>* string_data = jvFastaParser::convertFromReadData(output_seqs);
    jvFastaParser::dump(string_data, fasta_file.c_str());
    
    return 0;
}



    


