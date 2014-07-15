#include "jvStateSerializer.h"
#include "jvDataGenerator.h"
#include "jvFastaParser.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>



void print_bases(jvBaseVec* vec);

void jvStateSerializer::dump_fasta(jvState *jvs, const char *fname)
{
    vector<jvReadData> reads;
    for (vector<jvSequence*>::iterator it=jvs->getSeq().begin(); it != jvs->getSeq().end()-1; it++)
    {
        int left = (*it)->getLeft();
        int right = (*it)->getRight();
        jvBaseVec read(&(*it)->m_bases[(*it)->pos(left)], &(*it)->m_bases[(*it)->pos(right)]);
        reads.push_back(jvReadData(read));        
    }
    vector<jvStringReadData>* string_data = jvFastaParser::convertFromReadData(reads);
    jvFastaParser::dump(string_data, fname);
}


void jvStateSerializer::dump(jvState *jvs, const char *fname, bool dumpDataToFile, const char *fclust)
{
	string _method_name = "jvStateSerializer::dump";
	// sequence count
	// instance count
	// read count
	// alpha
	// gamma
	// left
	// for each sequence thetas
	// for each instance assignment
	// for each read assignment

	vector<jvReadAss> rdas = jvs->getReads();
	vector<jvInstance*> inss = jvs->getIns();
	vector<jvSequence*> seqs = jvs->getSeq();

	map<jvSequence*,int> seqMap;
	map<jvInstance*,int> insMap;

	ofstream ofs(fname,ifstream::out);
        
	if (!ofs.good())
	{
		cerr << _method_name << " file open failed" << endl;
		throw exception();
	}

	ofs << "#COMPRESSED" << endl;
        
	ofs << "#Sequence count " << endl << seqs.size() << endl;
	ofs << "#Instance count " << endl << inss.size() << endl;
	ofs << "#Read count " << endl << rdas.size() << endl;
	ofs << "#Alpha " << endl << jvs->getAlpha() << endl;
	ofs << "#Gamma " << endl << jvs->getGamma() << endl;       
	ofs << "#Sequences seqId,left,right, vector of bases " << endl;
	for (uint s=0;s<seqs.size();s++)
	{
		jvSequence* aSeq = seqs[s];
                int left = aSeq->getLeft();
                int right = aSeq->getRight();
                int bin = aSeq->getBin();
                ofs << "bin " << bin << endl;

		ofs << s << " " << left << " " << right << " ";
                
		for (int o=left; o<right; o++)
		{
                    ofs << ((int)aSeq->getBases(o));// << " ";
		}
		ofs << " " << endl;
		seqMap.insert( make_pair(aSeq,s) );	
	}
        
	ofs << "#Instances insId, seqId, offset" << endl;
	for (uint i=0;i<inss.size();i++)
	{
		jvInstance* anIns = inss[i];
		jvSequence* aSeq = anIns->getSeq();
                if (anIns->isAssigned() == false) continue;
		map<jvSequence*,int>::iterator iter = seqMap.find(aSeq);
		if (iter == seqMap.end())
		{
			cerr << _method_name << " sequence not in jvState" << endl;
			throw exception();
		}
		int sId = iter->second; 
		ofs << i << " " << sId << " " << anIns->getOffset() << endl;

		insMap.insert( make_pair(anIns,i) );		
	}

	ofs << "#Read readId, insId, offset,  length, vector of bases, raw_length, vector of raw bases, edit_length, vector for edit, alignment score" << endl;

        //////////////////////////////////////////////////
	ofstream cls_ofs;
        if (fclust != NULL) cls_ofs.open(fclust,ifstream::app);
        //////////////////////////////////////////////////
        
	for (uint r=0;r<rdas.size();r++)
	{
		jvReadAss ra = rdas[r];
		jvInstance* anIns = ra.getIns();
		jvReadData* rdata = ra.getData();

                int insId = -1;
                if (ra.isAssigned() == true)
                {
                    map<jvInstance*,int>::iterator iter = insMap.find(anIns);
                    if (iter == insMap.end())
                    {
			cerr << _method_name << " instance not in jvState" << endl;
			throw exception();
                    }                    
                    insId = iter->second; 
                }
                

		ofs << r << " " << insId << " " << ra.getOffset(); 
		if (dumpDataToFile)
		{
			ofs << " ";
			jvBaseVec aRead = rdata->getBases();
			ofs << aRead.size() << " ";
			for (uint ri=0;ri<aRead.size();ri++)
			{
                            ofs << ((int)aRead[ri]);// << " ";
			}
                        ofs << " ";

			jvBaseVec* aRawRead = rdata->getRawBases();
			ofs << aRawRead->size() << " ";
			for (uint ri=0;ri<aRawRead->size();ri++)
			{
                            ofs << ((int)((*aRawRead)[ri]));// << " ";
			}
                        ofs << " ";

			vector<char>* anEdit = rdata->getEdit();
			ofs << anEdit->size() << " ";
                        ofs << get_edit(*anEdit) << " ";

			ofs << rdata->getScore() << " ";
		}
		ofs << endl;

                //////////////////////////////////////////////////
                if (fclust != NULL)
                {
                    int seqId = -1;
                    if (insId != -1) 
                    {
                        jvInstance* anIns = inss[insId];
                        jvSequence* aSeq = anIns->getSeq();
                        if (anIns->isAssigned() == true)
                        {
                            map<jvSequence*,int>::iterator iter = seqMap.find(aSeq);
                            if (iter == seqMap.end())
                            {
                                cerr << _method_name << " sequence not in jvState" << endl;
                                throw exception();
                            }
                            seqId = iter->second; 
                        }
                    }
                    cls_ofs << seqId << " ";
                }
                //////////////////////////////////////////////////

	}

        //////////////////////////////////////////////////
        if (fclust != NULL)
        {
            cls_ofs << endl;
            cls_ofs.close();
        }
        //////////////////////////////////////////////////

}


template <class T> bool from_string(T& t, const string& s)
{
	istringstream iss(s);
	return !(iss >> std::dec >> t).fail();
}



void get_words(string str, vector<string*> *words)
{
	words->clear();

	size_t found,prevFound;
	str = " " + str + " ";
	found=str.find_first_of(" \t");
	prevFound = 0;
	while (found != string::npos)
	{
		if (found-prevFound>1)
		{			
			string *sub = new string(str.substr(prevFound+1,found-1-(prevFound+1)+1));
			words->push_back(sub);
		}
		prevFound = found;
		found=str.find_first_of(" \t",prevFound+1);
	}

}

enum parse_states {_seqct,_insct,_rdct,_alpha,_gamma,_eachseq,_eachins,_eachrd, _bb, _eachref };

jvState* jvStateSerializer::read(const char* fname, vector<jvReadData>& readData, bool dataInFile)
{
	string _method_name = "jvStateSerializer::read";

        cout<<"state serializer: Reading data"<<endl;
	int sequenceCount;
	int instanceCount;
	int readCount;
	double alpha;
	double gamma;
	bool done = false;
	if (dataInFile)
	{
		readData.clear();
	}

	parse_states currst = _seqct;

	vector<jvSequence*>* seqs = new vector<jvSequence*>;
	vector<jvInstance*>* inss = new vector<jvInstance*>;	
	vector<jvReadAss>  rdas;


	map<jvSequence*,int> seqMap;
	map<jvInstance*,int> insMap;
	vector<int>  readAssId;
	vector<int>  readAssOfs;
	vector<jvBaseVec> sequence_bases;   
    
	ifstream ifs(fname,ifstream::in);
	string line;
	vector<string*> words;

	int seqCt = -1;
	int insCt = -1;
	int readCt = -1;
    int bin = 0;

    cout<<"state serializer: going over file..."<<endl;
	while (!ifs.eof() && !done)
	{
		getline(ifs,line);
                
                if (line.compare("#COMPRESSED") == 0)
                {                
                    cout<<"compressed!\n";
                    ifs.close(); delete seqs; delete inss;
                    return read_compressed(fname, readData, dataInFile);
                }
                
		if (line[0] != '#')
		{
			get_words(line,&words);
			switch (currst)
			{
			case _seqct: 
			{
				if (words.size() != 1 || !from_string(sequenceCount,*words[0]))
				{
					cerr << _method_name << " sequence count parse failed" << endl;
					throw exception();
				}
				currst = _insct;
				break;
			}

			case _insct: 
			{
				if (words.size() != 1 || !from_string(instanceCount,*words[0]))
				{
					cerr << _method_name << " instance count parse failed" << endl;
					throw exception();
				}
				currst = _rdct;
				break;
			}	
			case _rdct:
			{
				if (words.size() != 1 || !from_string(readCount,*words[0]))
				{
					cerr << _method_name << " read count parse failed" << endl;
					throw exception();
				}
				currst = _alpha;
				break;
			}	
			case _alpha:
			{
				if (words.size() != 1 || !from_string(alpha,*words[0]))
				{
					cerr << _method_name << " alpha parse failed" << endl;
					throw exception();
				}
				currst = _gamma;
				break;			
			}	
			case _gamma:
			{
				if (words.size() != 1 || !from_string(gamma,*words[0]))
				{
					cerr << _method_name << " gamma parse failed" << endl;
					throw exception();
				}
				currst = _bb;
				break;
			}
			case _bb:
			{
                if (line[0] != 'b')
                    currst = _eachseq;
                else
                {
                    assert(words.size() == 2);
                    assert(from_string(bin,*words[1]));
                    currst = _eachseq;
                    break;
                }
			}
			
			case _eachseq:
			{
				if (words.size() <= 1 || seqCt+1>=sequenceCount)
				{
					cerr << "Seq no: " << seqCt << endl;
					cerr << _method_name << " each sequence parser failed" << endl;
					throw exception();
				}
				seqCt = seqCt+1;
				int tmp = 0;
				int left,right;
				if (!from_string(tmp,*words[0]) 	 || 
						tmp != seqCt 			   	 || 
						!from_string(left,*words[1])  || 
						!from_string(right,*words[2]) ||
						right - left != (int)words.size()-3) 
				{
					for (uint i=0;i<words.size();i++) 
					{
						cerr << *words[i] << " " ;
					}
					cerr << endl;
					cerr << "Header Seq no: " << seqCt << endl;
					cerr << _method_name << " each sequence parser failed" << endl;
					throw exception();
				}

				jvSequence* aSeq = new jvSequence();

				seqMap.insert( make_pair(aSeq,seqCt) );

				aSeq->m_left = 0;//left;
				aSeq->m_right = 1;//right;
				aSeq->pad(left,right);
                aSeq->m_bin = bin;  bin = 0;


				sequence_bases.push_back(jvBaseVec(0));
				for (uint wi=3;wi<words.size();wi++)
				{
					int base; 
					if (!from_string(base,*words[wi]))					
					{
						cerr << "Base Seq no: " << seqCt << endl;
						cerr << _method_name << " each sequence parser failed" << endl;
						throw exception();
					}
					sequence_bases[seqCt].push_back((jvBase)base);
				}
				seqs->push_back(aSeq);
				
                currst = _bb; // bin info
				if (seqCt == sequenceCount-1)
				{
					currst = _eachins;
				}
				break;
			}
			
			
			case _eachref:
			{
                assert(false);
                            
				int left,right, tmp;				
				if (!from_string(tmp,*words[0]) 	 || 
					tmp != seqCt 			   	 || 
				    !from_string(left,*words[1])  || 
					!from_string(right,*words[2]) ||
					right - left != (int)words.size()-3) 
				{
					for (uint i=0;i<words.size();i++) 
					{
						cerr << *words[i] << " " ;
					}
					cerr << endl;
					cerr << "ref " << seqCt << endl;
					cerr << _method_name << " ref parser failed" << endl;
					throw exception();
				}
				
				jvSequence* aSeq  = seqs->back(); 
				
				for (uint wi=3;wi<words.size();wi++)
				{
					int index; 
					if (!from_string(index,*words[wi]))					
					{
						cerr << "index ref for Sequnce " << seqCt << endl;
						cerr << _method_name << " ref parser failed" << endl;
						throw exception();
					}
				}
				assert( left == right);
			    currst = _eachseq;
				if (seqCt == sequenceCount-1)
				{
					currst = _eachins;
				}

				break;    			    
	        }
						
			case _eachins:
			{
				if (words.size() <= 1 || insCt+1>=instanceCount)
				{
					cerr << _method_name << " each instance parser failed" << endl;
					throw exception();
				}
				insCt = insCt+1;
				int tmp = 0;
				int seqId,o;
				if (!from_string(tmp,*words[0]) 	 || 
						tmp != insCt 			   	 || 
						!from_string(seqId,*words[1]) || 
						!from_string(o,*words[2])     ||
						seqId >= sequenceCount)						 
				{
					cerr << _method_name << " each instance parser failed" << endl;
					throw exception();
				}


				jvInstance* anIns = new jvInstance((*seqs)[seqId],o);					
				inss->push_back(anIns);					

				insMap.insert( make_pair(anIns,insCt));
				if (insCt == instanceCount-1)
				{
					currst = _eachrd;
				}
				break;
			}
			case _eachrd:
			{

				if (words.size() <= 1 || readCt+1>=readCount)
				{

					cerr << "ReadCt: " << readCt << " readCount " << endl;
					cerr << _method_name << "each read parser failed" << endl;
					throw exception();
				}
				readCt = readCt+1;
				int tmp;				
				int insId; 
				int offset;
				int len,rlen,elen;
				double score;


				if (!from_string(tmp,*words[0]) 	     ||	tmp != readCt 			   	 ||
						!from_string(insId,*words[1])   ||						    
						!from_string(offset,*words[2])    ||	insId >= instanceCount ||
						(!dataInFile && words.size() != 3))
				{
					cerr << "Header Read no: " << readCt << endl;

					cerr << _method_name << " each read parser failed" << endl;
					throw exception();
				}

				if (dataInFile)
				{
					int last = 3;
					if (!from_string(len,*words[last])   ||							
							len > (int)words.size() - last+1)
					{
						cerr << "Data Read no: " << readCt << endl;
						cerr << _method_name << " each read parser failed" << endl;
						throw exception();
					}

					last++;


					jvBaseVec read;
					jvBaseVec rawRead;
					vector<char> edit;

					for (int wi=last;wi<last+len;wi++)
					{
						int base; 
						if (!from_string(base,*words[wi]))					
						{
							cerr << "Datap Read no: " << readCt << endl;
							cerr << _method_name << " each read parser failed" << endl;
							throw exception();
						}
						read.push_back((jvBase)base);
					}

					last = last+len;


					if (!from_string(rlen,*words[last]) || 
							rlen > (int)words.size() - last+1)
					{
						cerr << "Raw Read no: " << readCt << endl;
						cerr << _method_name << " each read parser failed" << endl;
						throw exception();							
					}
					last++;

					for (int wi=last;wi < last+rlen;wi++)
					{
						int base; 
						if (!from_string(base,*words[wi]))					
						{
							cerr << "RawP Read no: " << readCt << endl;
							cerr << _method_name << " each read parser failed" << endl;
							throw exception();
						}
						rawRead.push_back((jvBase)base);
					}

					last = last + rlen;

					if (!from_string(elen,*words[last]) || 
							elen > (int)words.size() - last+1)
					{
						cerr << "Edit Read no: " << readCt << endl;
						cerr << _method_name << " each read parser failed" << endl;
						throw exception();							
					}

					last++;

					for (int wi=last;wi < last+elen;wi++)
					{
						int base; 
						if (!from_string(base,*words[wi]))					
						{
							cerr << "Editp Read no: " << readCt << endl;
							cerr << _method_name << " each read parser failed" << endl;
							throw exception();
						}
						edit.push_back((char)base);
					}

					last = last + elen;

					if (!from_string(score,*words[last]))					
					{
						cerr << "Score Read no: " << readCt << endl;
						cerr << _method_name << " each read parser failed" << endl;
						throw exception();
					}

					jvReadData rddata(read, rawRead,edit, score);
					readData.push_back(rddata);
				}

				readAssId.push_back(insId);
				readAssOfs.push_back(offset);

				if (readCt == readCount-1)
				{
					done = true;
				}	
				break;
			} // end case
			} // end switch
			for (uint wi=0;wi<words.size();wi++)
			{
				delete words[wi];
			}
			words.clear();						
		}
	}

        cout<<"State Serializer: Finished reading file"<<endl;
        
	for (uint i=0;i<readData.size();i++) 
	{
		jvReadAss rdass(&(readData[i]));
                if (readAssId[i] != -1)
                    rdass.setIns((*inss)[readAssId[i]],readAssOfs[i]);					
                else
                    rdass.setIns(NULL);
		rdas.push_back(rdass);
	}

	//  Createing a new state
	//  (The default for the new state has m_with_instances = true)
	jvState* jvs = new jvState(readData, true);

	jvs->m_alpha = alpha;
	jvs->m_gamma = gamma;
	jvs->m_ins = *inss;
	jvs->m_seq = *seqs;
	jvs->m_reads = rdas;

	//  An assumption: All the instances listed in the file are assigned!
    //      Reads may be unassigned if their insId is -1.

	for (uint i=0;i<inss->size();i++) 
	{
            jvInstance* jvi = (*inss)[i];
            jvi->assign(true);
	}

	for (uint i=0;i<jvs->m_reads.size();i++) 
	{
            if (jvs->m_reads[i].getIns() != NULL)
                jvs->m_reads[i].assign(true);		
	}

	jvs->freezeInstances(true);

    cout<<"State Serializer: All that is left is to recopy the bases"<<endl;

	for (uint i=0;i<jvs->m_seq.size();i++)
	{
            if (sequence_bases[i].size() == jvs->m_seq[i]->getRight()-jvs->m_seq[i]->getLeft())
		copy(sequence_bases[i].begin(),sequence_bases[i].end(), jvs->m_seq[i]->m_bases.begin()+jvs->m_seq[i]->pos(jvs->m_seq[i]->getLeft()));
	}

        delete inss;
        delete seqs;
        ifs.close();
        
	return jvs;
}




jvState* jvStateSerializer::read_compressed(const char* fname, vector<jvReadData>& readData, bool dataInFile)
{
	string _method_name = "jvStateSerializer::read";

        cout<<"state serializer: Reading data"<<endl;
	int sequenceCount;
	int instanceCount;
	int readCount;
	double alpha;
	double gamma;
	bool done = false;
	if (dataInFile)
	{
		readData.clear();
	}

	parse_states currst = _seqct;

	vector<jvSequence*>* seqs = new vector<jvSequence*>;
	vector<jvInstance*>* inss = new vector<jvInstance*>;	
	vector<jvReadAss>  rdas;


	map<jvSequence*,int> seqMap;
	map<jvInstance*,int> insMap;
	vector<int>  readAssId;
	vector<int>  readAssOfs;
	vector<jvBaseVec> sequence_bases;    
    
	ifstream ifs(fname,ifstream::in);
	string line;
	vector<string*> words;

	int seqCt = -1;
	int insCt = -1;
	int readCt = -1;

        int bin = 0;

        cout<<"state serializer: going over file..."<<endl;
	while (!ifs.eof() && !done)
	{
		getline(ifs,line);
		if (line[0] != '#')
		{
			get_words(line,&words);
			switch (currst)
			{
			case _seqct: 
			{
				if (words.size() != 1 || !from_string(sequenceCount,*words[0]))
				{
					cerr << _method_name << " sequence count parse failed" << endl;
					throw exception();
				}
				currst = _insct;
				break;
			}

			case _insct: 
			{
				if (words.size() != 1 || !from_string(instanceCount,*words[0]))
				{
					cerr << _method_name << " instance count parse failed" << endl;
					throw exception();
				}
				currst = _rdct;
				break;
			}	
			case _rdct:
			{
				if (words.size() != 1 || !from_string(readCount,*words[0]))
				{
					cerr << _method_name << " read count parse failed" << endl;
					throw exception();
				}
				currst = _alpha;
				break;
			}	
			case _alpha:
			{
				if (words.size() != 1 || !from_string(alpha,*words[0]))
				{
					cerr << _method_name << " alpha parse failed" << endl;
					throw exception();
				}
				currst = _gamma;
				break;			
			}	
			case _gamma:
			{
				if (words.size() != 1 || !from_string(gamma,*words[0]))
				{
					cerr << _method_name << " gamma parse failed" << endl;
					throw exception();
				}
				currst = _bb;
				break;
			}
			case _bb:
			{
                            if (line[0] != 'b')
                                currst = _eachseq;
                            else
                            {
                                assert(words.size() == 2);
                                assert(from_string(bin,*words[1]));
                                currst = _eachseq;
                                break;
                            }
			}
			
			case _eachseq:
			{
				if (words.size() <= 1 || seqCt+1>=sequenceCount)
				{
					cerr << "Seq no: " << seqCt << endl;
					cerr << _method_name << " each sequence parser failed" << endl;
					throw exception();
				}
				seqCt = seqCt+1;
				int tmp = 0;
				int left,right;
				if (!from_string(tmp,*words[0])   || 
                                    tmp != seqCt 		  || 
                                    !from_string(left,*words[1])  || 
                                    !from_string(right,*words[2]) ||
                                    (int)words.size() != 4 ||
                                    right - left != words[3]->length()) 
				{
					for (uint i=0;i<words.size();i++) 
					{
						cerr << *words[i] << " " ;
					}
					cerr << endl;
					cerr << "Header Seq no: " << seqCt << endl;
					cerr << _method_name << " each sequence parser failed" << endl;
					throw exception();
				}

				jvSequence* aSeq = new jvSequence();

				seqMap.insert( make_pair(aSeq,seqCt) );

				aSeq->m_left = 0;//left;
				aSeq->m_right = 1;//right;
				aSeq->pad(left,right);
                aSeq->m_bin = bin;  bin = 0;


				sequence_bases.push_back(jvBaseVec(0));
				for (uint i=0;i<words[3]->length();i++)
				{
                    char base = (*words[3])[i]-'0'; 
                    sequence_bases[seqCt].push_back((jvBase)base);
				}
				seqs->push_back(aSeq);
				
                currst = _bb; // bin info
				if (seqCt == sequenceCount-1)
				{
					currst = _eachins;
				}
				break;
			}
			
			
			case _eachref:
			{
                            assert(false);
                            
				int left,right, tmp;				
				if (!from_string(tmp,*words[0]) 	 || 
					tmp != seqCt 			   	 || 
				    !from_string(left,*words[1])  || 
					!from_string(right,*words[2]) ||
					right - left != (int)words.size()-3) 
				{
					for (uint i=0;i<words.size();i++) 
					{
						cerr << *words[i] << " " ;
					}
					cerr << endl;
					cerr << "ref " << seqCt << endl;
					cerr << _method_name << " ref parser failed" << endl;
					throw exception();
				}
				
				jvSequence* aSeq  = seqs->back(); 

				for (uint wi=3;wi<words.size();wi++)
				{
					int index; 
					if (!from_string(index,*words[wi]))					
					{
						cerr << "index ref for Sequnce " << seqCt << endl;
						cerr << _method_name << " ref parser failed" << endl;
						throw exception();
					}
				}
				assert( left == right);
			    currst = _eachseq;
				if (seqCt == sequenceCount-1)
				{
					currst = _eachins;
				}

				break;    			    
	        }
						
			case _eachins:
			{
				if (words.size() <= 1 || insCt+1>=instanceCount)
				{
					cerr << _method_name << " each instance parser failed" << endl;
					throw exception();
				}
				insCt = insCt+1;
				int tmp = 0;
				int seqId,o;
				if (!from_string(tmp,*words[0]) 	 || 
						tmp != insCt 			   	 || 
						!from_string(seqId,*words[1]) || 
						!from_string(o,*words[2])     ||
						seqId >= sequenceCount)						 
				{
					cerr << _method_name << " each instance parser failed" << endl;
					throw exception();
				}


				jvInstance* anIns = new jvInstance((*seqs)[seqId],o);					
				inss->push_back(anIns);					

				insMap.insert( make_pair(anIns,insCt));
				if (insCt == instanceCount-1)
				{
					currst = _eachrd;
				}
				break;
			}
			case _eachrd:
			{

				if (words.size() <= 1 || readCt+1>=readCount)
				{

					cerr << "ReadCt: " << readCt << " readCount " << endl;
					cerr << _method_name << "each read parser failed" << endl;
					throw exception();
				}
				readCt = readCt+1;
				int tmp;				
				int insId; 
				int offset;
				int len,rlen,elen;
				double score;


				if (!from_string(tmp,*words[0]) 	     ||	tmp != readCt 			   	 ||
						!from_string(insId,*words[1])   ||						    
						!from_string(offset,*words[2])    ||	insId >= instanceCount ||
						(!dataInFile && words.size() != 3))
				{
					cerr << "Header Read no: " << readCt << endl;

					cerr << _method_name << " each read parser failed" << endl;
					throw exception();
				}

				if (dataInFile)
				{
					int last = 3;
					if (!from_string(len,*words[last])   ||	 (int)words.size() == last + 1)
					{
						cerr << "Data Read no: " << readCt << endl;
						cerr << _method_name << " each read parser failed" << endl;
						throw exception();
					}

					last++;


					jvBaseVec read;
					jvBaseVec rawRead;
					vector<char> edit;


                                        assert(len == words[last]->length());
					for (int i=0;i<len;i++)
                                        {
                                            read.push_back((jvBase)( (*words[last])[i]-'0'));
                                        }
                                                           
					last++;


					if (!from_string(rlen,*words[last]) ||  (int)words.size() == last + 1)
					{
						cerr << "Raw Read no: " << readCt << endl;
						cerr << _method_name << " each read parser failed" << endl;
						throw exception();							
					}
					last++;

                                        assert(rlen == words[last]->length());
					for (int i=last;i < rlen;i++)
					{
                                            rawRead.push_back((jvBase)( (*words[last])[i]-'0'));
					}

					last++;

					if (!from_string(elen,*words[last]) || (int)words.size() == last + 1 )
					{
						cerr << "Edit Read no: " << readCt << endl;
						cerr << _method_name << " each read parser failed" << endl;
						throw exception();							
					}

					last++;

                                        set_edit(*words[last], edit);
                                        assert(elen == (int)edit.size());
					last = last++;

					if (!from_string(score,*words[last]))					
					{
						cerr << "Score Read no: " << readCt << endl;
						cerr << _method_name << " each read parser failed" << endl;
						throw exception();
					}

					jvReadData rddata(read, rawRead,edit, score);
					readData.push_back(rddata);
				}

				readAssId.push_back(insId);
				readAssOfs.push_back(offset);

				if (readCt == readCount-1)
				{
					done = true;
				}	
				break;
			} // end case
			} // end switch
			for (uint wi=0;wi<words.size();wi++)
			{
				delete words[wi];
			}
			words.clear();						
		}
	}

    cout<<"State Serializer: Finished reading file"<<endl;
        
	for (uint i=0;i<readData.size();i++) 
	{
		jvReadAss rdass(&(readData[i]));
                if (readAssId[i] != -1)
                    rdass.setIns((*inss)[readAssId[i]],readAssOfs[i]);					
                else
                    rdass.setIns(NULL);
		rdas.push_back(rdass);
	}

	//  Createing a new state
	//  (The default for the new state has m_with_instances = true)
	jvState* jvs = new jvState(readData, true);

	jvs->m_alpha = alpha;
	jvs->m_gamma = gamma;
	jvs->m_ins = *inss;
	jvs->m_seq = *seqs;
	jvs->m_reads = rdas;

	//  An assumption: All the instances listed in the file are assigned!
    //      Reads may be unassigned if their insId is -1.

	for (uint i=0;i<inss->size();i++) 
	{
            jvInstance* jvi = (*inss)[i];
            jvi->assign(true);
	}

	for (uint i=0;i<jvs->m_reads.size();i++) 
	{
            if (jvs->m_reads[i].getIns() != NULL)
                jvs->m_reads[i].assign(true);		
	}

	jvs->freezeInstances(true);

    cout<<"State Serializer: All that is left is to recopy the bases"<<endl;

	for (uint i=0;i<jvs->m_seq.size();i++)
	{
            if (sequence_bases[i].size() == jvs->m_seq[i]->getRight()-jvs->m_seq[i]->getLeft())
		copy(sequence_bases[i].begin(),sequence_bases[i].end(), jvs->m_seq[i]->m_bases.begin()+jvs->m_seq[i]->pos(jvs->m_seq[i]->getLeft()));
	}

        delete inss;
        delete seqs;
        
	return jvs;
}
