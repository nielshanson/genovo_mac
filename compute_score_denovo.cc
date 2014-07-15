#include <iostream>
#include <vector>
#include <math.h>
#include <assert.h>
#include <fstream>
#include <string>
using namespace std;

#include "jvDataGenerator.h"
#include "jvState.h"
#include "jvSampler.h"
#include "jvStateSerializer.h"
#include "jvAlign.h"
#include "jvFastaParser.h"


vector<bool> success;

extern int iteration;
extern int out_of_bounds;

int main(int argc, char ** argv)
{
    if ( (argc < 3) || (argc > 4) )
    {
        cout<<"Usage:\t compute_score_denovo <read_fasta_file> <contig_fasta_file>\n\n";
        cout<<"Compute score_{denovo} and other statistics.\n\n";
        
        cout<<"compute_score_denovo <read_fasta_file> <contig_fasta_file> 1\n";
        cout<<"will output a state file mapping.dump with the best mapping of the reads to the contings.\n\n";
        return 1;
    }

    vector<jvReadData> data;
    vector<jvReadData> skeleton;
    jvState* st = NULL;
    uint j=0;    
    
    string pathToReads = string(argv[1]);
    string pathToContigs = string(argv[2]);
    bool construct_state = (argc >= 4) && (atoi(argv[3]) == 1);

    // load reads from file
    jvFastaParser::parseReadData(pathToReads.c_str(),data);

    // load skeletons from file
    jvFastaParser::parseReadData(pathToContigs.c_str(),skeleton);

    double gamma = pow(2,80);
       
    cout<<"Initialize an empty state...\n";
    st = new jvState(data);
    st->setAuto();
    
    st->setGamma(gamma);
    st->setAlignmentOnly(construct_state);
    st->setBothOrientations(true);

    int total_length = 0;
    for (int k=0; k<skeleton.size(); k++)
    {
        st->newSeq(new jvSequence(&skeleton[k].getBases()));
        total_length += (int)skeleton[k].getBases().size();
    }   

    clock_t toc, tic = clock();
    st->gibbsReads();
    if (construct_state)
    {
        jvStateSerializer::dump(st, "mapping.dump", true);
    }

    
    cout<<"Collecting read information...\n";
    int garbage_count = 0, total_garbage_length = 0, total_skeleton_length = 0;
    double total_score = 0;
    for (int i=0; i<st->getReads().size(); i++)
    {
        if (! st->getReads()[i].getIns()) 
        {
            garbage_count++;
            total_garbage_length += (int)st->getReads()[i].getData()->getRawBases()->size();
        }
        else
        {
            int read_len = (int)st->getReads()[i].getData()->getRawBases()->size();
            total_skeleton_length += read_len;
            double per_base_score = st->getReads()[i].getData()->getScore() /read_len;
            total_score += st->getReads()[i].getData()->getScore();
           }
    }
    gamma = 20*log(4);
    double skeleton_ll = total_score; // ll.a + ll.y        
    double garbage_ll = total_garbage_length*(logBasePrior+errorModel.logInvNoise());
    int N = st->getReads().size();
    int skeleton_count = N-garbage_count;
    double raw_score = total_length*logBasePrior + garbage_ll + skeleton_ll + gamma*(skeleton.size()+garbage_count);
    double naive_score = (total_garbage_length + total_skeleton_length)*(logBasePrior+errorModel.logInvNoise()) + gamma*N;
    double score = (raw_score-naive_score)/(total_garbage_length + total_skeleton_length);
    
    const char tab[] = "\t";    
    cout<<endl<<"===================================================================================="<<endl;   
    cout<<"* read file: "<<pathToReads.c_str()<<tab<<"contig file: "<<pathToContigs.c_str()<<endl<<endl;    
    cout<<"* no. reads: "<<N<<tab<<"Percent garbage: "<<100*garbage_count/N<<"%"<<" ("<<garbage_count<<" reads)"<<endl;
    cout<<"* no. contigs: "<<skeleton.size()<<tab<<"total contig length: "<<total_length<<endl;
    cout<<"* alignment score: "<<skeleton_ll<<tab<<"garbage score: "<<garbage_ll<<endl;
    cout<<"* Score: "<<score<<tab<<"(raw score: "<<raw_score<<" )"<<endl;
    cout<<endl<<"===================================================================================="<<endl;   
    toc = clock();
    cout<<"Done in "<<(double)(toc-tic)/CLOCKS_PER_SEC<<" seconds.\n";


    char file_str[100];
    sprintf(file_str, "jobs/%s.done", argv[2]);
    ofstream myfile;
    myfile.open(file_str, ios::app);
    myfile<<pathToReads.c_str()<<tab<<pathToContigs.c_str()<<tab;   
    myfile<<"#reads:"<<tab<<N<<tab<<"Percent_garbage:"<<tab<<100*garbage_count/N<<tab;
    myfile<<"#contigs:"<<tab<<skeleton.size()<<tab<<"total_contig_length:"<<tab<<total_length<<tab;
    myfile<<"alignment_score:"<<tab<<skeleton_ll<<tab<<"garbage_score:"<<tab<<garbage_ll<<tab;
    myfile<<"raw_Score:"<<tab<<raw_score<<tab<<"norm_score:"<<tab<<score<<endl;
    myfile.close();
    
    jvSampler::release();
    delete st;
    return 0;
}
