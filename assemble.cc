//  This test is for denovo sequencing of one sequence.
//  no backbone.


#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string>
#include <assert.h>
// #include <malloc.h>

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

#define BUFFLEN 10000


int main(int argc, char ** argv)
{
    if ( (argc < 2) || (argc > 5) )
    {
        cout<<"Usage:\n";
        cout<<"- assemble <fasta_file>\n";
        cout<<"will run Genovo for 10,000 iterations on the set of reads in <fasta_file>.\n\n";        

        cout<<"- assemble <fasta_file> N\n";        
        cout<<"will run Genovo for N iterations.\n\n";

        cout<<"- assemble <fasta_file> N <dump_file>\n";
        cout<<"will run Genovo for N iterations, loading the initial state from <dump_file>.\n\n";

        cout<<"Inputs:\n";
        cout<<"<fasta_file> is a fasta file of reads.\n";
        cout<<"<dump_file> is the .dump or .dump.best file generated as output by a previous Genovo run.\n\n";

        cout<<"Outputs:\n";
        cout<<"The output files are updated during the course of the algorithm.  It is ok to view them while the algorithm is running.  It is ok to kill the algorithm in which case the output files will represent the most updated results.\n\n";

        cout<<"Output all contigs with length > 500b to file genovo.fa by running:\n";
	cout<<"\tfinalize 500 genovo.fa <fasta_file>.dump.best\n\n"; 
        return 1;
    }
        
    string pathToDataset = string(argv[1]);
    string pathToDump0 = pathToDataset + ".dump0";
    string pathToDump1 = pathToDataset + ".dump1";
    string pathToClust = pathToDataset + ".clust";
    string pathToStatus = pathToDataset + ".status";
    string pathToBest = pathToDataset + ".dump.best";
    string pathToBestFasta = pathToDataset + ".dump.best.fa";

    const char *dump_str_even = pathToDump0.c_str();
    const char *dump_str_odd = pathToDump1.c_str();
    const char *dump_str_clust = pathToClust.c_str();
    const char *dump_str_best = pathToBest.c_str();
    const char *dump_str_best_fasta = pathToBestFasta.c_str();
    const char *status_str = pathToStatus.c_str();

    int max_iter = ( argc >= 3 ? atoi(argv[2]) : 10000 );
    int seed = 10;
    uint j = ( argc >= 5 ? atoi(argv[4]) : 0 );
    int pow_gamma = 35;

    vector<jvReadData> data;
    jvState* st = NULL;
    uint width = 100;
    int endAl = 450;
    double gamma = pow(2,pow_gamma);

    if (argc < 4)
    {
        cout<<"Loading fasta file...\n";
        jvFastaParser::parseReadData(pathToDataset.c_str(),data);
        cout<<"Initialize an empty state...\n";
        st = new jvState(data, false, seed);
    }
    else // argc == 4
    {
        cout<<"Loading dump file "<<argv[3]<<"..."<<flush;
        st = jvStateSerializer::read(argv[3], data);
        gamma = st->getGamma();
        cout<<" Done!"<<endl;
    }

    st->setAuto();
    st->setErrorModel(0.01, 0.01);
    st->setWithInstances(false);    
    st->setGamma(gamma);
    st->setBothOrientations(true);
    st->setGeometric(true);
    jvAlign::with_cache = true;

    uint snm_count = 0;
    double best_ll = -1e100;
    bool is_gap = false;

    cout<<"Running for "<<max_iter<<" iterations.\n";
    for (; j<max_iter; j++)
    {
    	iteration = j;
        clock_t toc, tic = clock();
        cout<<"\n----------------   Iteration "<<j<<"   ---------------------\n";

        jvSampler::setMAP((j+1)%100 == 0);  //every 100th iteration take the max    
        
        cout<<"Gibbs Reads...\n";
        if (j>1 && ( (j+1)%10 != 0))
            st->gibbsReads(vector<jvReadAss*>(), GIBBS_FREEZE);
        else
            st->gibbsReads();       
        st->center();

        toc = clock();
        cout<<"gibbsReads done in "<<(double)(toc-tic)/CLOCKS_PER_SEC<<" seconds.\n";
        
        st->sampleBases();
        cout<<"sampleBases done in "<<(double)(clock()-toc)/CLOCKS_PER_SEC<<" seconds.\n";
        toc = clock();
  
        st->cleanAlignmentCache(j-11);
        cout<<"cleaning cache done in "<<(double)(clock()-toc)/CLOCKS_PER_SEC<<" seconds.\n";
        toc = clock();
        
        st->summarizeEdits();
        cout<<"sumEdit done in "<<(double)(clock()-toc)/CLOCKS_PER_SEC<<" seconds.\n";
        toc = clock();

       	st->saveStatus(status_str, iteration);
        cout<<"saving status done in "<<(double)(clock()-toc)/CLOCKS_PER_SEC<<" seconds.\n";
        toc = clock();
        
        st->freezeInstances(true);
        cout<<"freezing instances done in "<<(double)(clock()-toc)/CLOCKS_PER_SEC<<" seconds.\n";
        toc = clock();
            
        if ( (j+1)%10 == 0)
        {
            if (st->ll().ll() > best_ll)
            {
                jvStateSerializer::dump(st, dump_str_best, true);
                best_ll = st->ll().ll();            
            }
        }
        cout<<"saving done in "<<(double)(clock()-toc)/CLOCKS_PER_SEC<<" seconds.\n";
	    toc = clock();
        
        if ( (j+2)%5 == 0)
		st->unassignEnds();
        cout<<"unAssignEnds done in "<<(double)(clock()-toc)/CLOCKS_PER_SEC<<" seconds.\n";
        toc = clock();
        
        st->sampleBases();
        st->freezeInstances(false);
	    cout<<"Align/Split ends...\n";
        if ((j>0) || (argc == 4))
            st->findAndAlignEnds();
        cout<<"alignEnds done in "<<(double)(clock()-toc)/CLOCKS_PER_SEC<<" seconds.\n";
        toc = clock();

        st->flipSequences();
        cout<<"flipSequences done in "<<(double)(clock()-toc)/CLOCKS_PER_SEC<<" seconds.\n";
	    toc = clock();
        
        if ((j>0) || (argc == 4))
            st->checkAllOffsets();
        
        cout<<"sampleBases and CAO done in "<<(double)(clock()-toc)/CLOCKS_PER_SEC<<" seconds.\n";
	    toc = clock();

        if ( (j+1)%10 == 0)
	        jvStateSerializer::dump(st, (j % 2 ? dump_str_odd : dump_str_even),true);
        cout<<"saving done in "<<(double)(clock()-toc)/CLOCKS_PER_SEC<<" seconds.\n";

	    toc = clock();
        cout<<"Done in "<<(double)(toc-tic)/CLOCKS_PER_SEC<<" seconds.\n";
        cout<<"\n--------------------------------------------------------\n";

}


    jvSampler::release();
    delete st;
    return 0;
}
