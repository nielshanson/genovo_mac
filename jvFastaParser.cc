#include "jvFastaParser.h"
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

void jvFastaParser::parseReadData(const char* name, vector<jvReadData>& readData)
{
	vector<jvStringReadData>* data = jvFastaParser::parseStringReadData(name);
	readData.clear();
	for (uint i=0;i<data->size();i++)
	{
		readData.push_back((*data)[i].convertToReadData());
	}
	delete data;
}

vector<jvStringReadData>* jvFastaParser::parseStringReadData(const char* name)
{
	vector<jvStringReadData> *data = new vector<jvStringReadData>(); 
	string line;
	string buff;
	string meta;
	bool hasData = false;
	
	int ct = 0;
	ifstream ifs (name, ifstream::in);
	if (!ifs.good())
	{
		cout << "Error: fasta file " << name << endl;
		throw exception();		
	}
	while (1==1)
	{
		getline(ifs,line);
		if (line[0] =='>' || ifs.eof())
		{
			if (hasData)
			{							
				data->push_back(jvStringReadData(meta,buff));
			}
			
			if (ifs.eof())
			{
				break;
			}
			else
			{
				meta = line.substr(1);
				ct = ct+1;			
				hasData = false;
				buff.clear();
				line.clear();
			}			
			
		}
		else
		{
			hasData = true;
			buff += line;
			line.clear();
		}				
	}
	ifs.close();
	return data;
}

vector<jvStringReadData>* jvFastaParser::convertFromReadData(vector<jvReadData>& readData)
{
	
	vector<jvStringReadData>* data = new vector<jvStringReadData>;
		
	for (uint i=0;i<readData.size();i++)
	{
		data->push_back(jvStringReadData(readData[i]));		
	}	
	
	return data;
}

void jvFastaParser::dump(vector<jvStringReadData>* data, const char * fname)
{
	ofstream ofs (fname, ifstream::out);
	if (!ofs.good())
	{
		cout << "Error: fasta file " << fname << endl;
		throw exception();			
	}
	for (uint i=0;i<data->size();i++)
	{		
		jvStringReadData jvrs = (*data)[i];
                if (jvrs.meta.empty())
                    ofs << '>' << i+1 << " length "<<jvrs.read.length()<<endl;
                else
                    ofs << '>' << jvrs.meta << endl;
                        
		for (uint l=0;l<jvrs.read.length();l+=60)
		{
			string tmp = jvrs.read.substr(l,60); 
			ofs << tmp << endl;		
		}
	}
	ofs.close();	
}


void jvFastaParser::dump_fastq(vector<jvStringReadData>* data, const char * fname)
{
	ofstream ofs (fname, ifstream::out);
	if (!ofs.good())
	{
		cout << "Error: fastq file " << fname << endl;
		throw exception();			
	}
	for (uint i=0;i<data->size();i++)
	{		
		jvStringReadData jvrs = (*data)[i];
                if (jvrs.meta.empty())
                    ofs << '@' << i+1 << " length "<<jvrs.read.length()<<endl;
                else
                    ofs << '@' << jvrs.meta << endl;

                ofs << jvrs.read << endl;
                ofs << '+' <<endl;
		for (uint l=0;l<jvrs.read.length();l++)
		{
			ofs << 'I';		
		}
                ofs<<endl;
	}
	ofs.close();	
}


jvStringReadData::jvStringReadData(jvReadData jvread)
{
	jvBaseVec iread = jvread.getBases();
	read.clear();
	
	
	for (uint i=0;i<iread.size();i++)
	{
		switch(iread[i]) 
		{
		case 0: read.append("A");break; 
		case 1: read.append("C");break;
		case 2: read.append("G");break;
		case 3: read.append("T");break;	
		case 4: read.append("N");break;
		default : cerr<< " Unknown nucleotide " << endl;throw exception(); break;
		}		
	}	
}

jvReadData jvStringReadData::convertToReadData()
{
	jvBaseVec iread;
	for (uint i=0;i<read.size();i++)
	{
		switch((read)[i]) 
		{
		case 'A':
                case 'a': iread.push_back(0);break; 
                case 'c':
		case 'C': iread.push_back(1);break;
                case 'g':
		case 'G': iread.push_back(2);break;
                case 't':
		case 'T': iread.push_back(3);break;
                case 'n':
		case 'N': iread.push_back(4);break;
		default : cerr<< " Unknown nucleotide " << endl;throw exception(); break;
		}		
	}
	return jvReadData(iread,0,0,0);
}

