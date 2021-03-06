// RJMcmc.cpp : Defines the entry point for the console application.
//
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include "Matrix.h"
#include "ClusterFlags.h"
#include "mcmc.h"
#include "Random.h"
#include "MCMCEnv.h"
#include <stdexcept>
#include <string>
#include <sstream>

using namespace std;


// With a new dataset, do not forget modify n_clusterhome!!!

void initializeParameters(long NumOfGenes, long n_time, const matrix &Data, double k1, double k0, bool bForest) {
	
	ClusterTree::setNumOfTimePoints(n_time);

	matrix mu1(1,NumOfGenes);
	matrix mu0(1,NumOfGenes);
	matrix sigma0(NumOfGenes,1.0,1);
	matrix sigma1(NumOfGenes,1.0,1);
	double beta0 = 1000.0;
	double beta1 = 1000.0;
	double deltar = 3.0;
	double alpha = 1.0;
	double shape = 1.0;
	double bk = 0.5; // For split_merge move and birth_death move, we give equal probability.
	double ranc  = 0.9; // Tail birth is given much more porbability.
	//double k1 = .3;
	//double k0 = .022;
	// K1 and deltar are both used to control data effect, it should be 
	// Note that K1 must be <= 1;
	// With a large or small K1, the algorithm both seems converge very slowly.
	mu1 = Data.matrix_mean();
	mu0 = mu1;
	sigma1 *= 1.0/k1;
	sigma0 *= 1.0/k0;

	// Beta controls the prior effect to the sampling probability
	// A large Beta value should be set, if you want to decrease the prior effect.
	// The data effect is controlled by deltar, if you want to increase data effect, a large deltar should be set.
	// However the minimum of 

	MCMCEnv::initializeParameters(mu1, mu0, sigma1, sigma0, beta1, beta0, deltar, alpha, shape, bk, ranc, bForest);

}

void getInputsFromFile(istream &is, long &n_time, long &n_gene, long &n_num,
	vector<unsigned long> &n_Nt, matrix &datamat, double &k1, double &k0,
	vector<pair<unsigned long, unsigned long> > &weightToMon, bool &bForest) {
		if(!is) {
			throw(logic_error("Cannot open file!"));
		}
		string varName;
		dvector datavec;
		while(is && varName != "END") {
			is >> varName;
			if(!varName.empty() 
				&& (*(varName.begin()) == '#' || *(varName.begin()) == '/')) {
					// is comment
					getline(is, varName);
			} else if(varName == "TIMEPOINTS") {
				// time points following
				is >> n_time >> ws;
			} else if(varName == "GENES") {
				is >> n_gene >> ws;
			} else if(varName == "TOTAL_SAMPLES") {
				is >> n_num >> ws;
			} else if(varName == "SAMPLES_BY_TIME") {
				n_Nt.clear();
				is.exceptions(istream::failbit | istream::eofbit | istream::badbit);
				is.clear();
				while(true) {
					try {
						unsigned long sampleByTime;
						is >> sampleByTime;
						n_Nt.push_back(sampleByTime);
					} catch(ios::failure &e) {
						is.clear();
						is.exceptions(istream::goodbit);
						is >> ws;
						if(((char) is.peek()) == '#' || ((char) is.peek()) == '/') {
							getline(is, varName);
						} else {
							break;
						}
					}
				}
			} else if(varName == "K0") {
				is >> k0 >> ws;
			} else if(varName == "K1") {
				is >> k1 >> ws;
			} else if(varName == "DATA") {
				datavec.clear();
				is.exceptions(istream::failbit | istream::eofbit | istream::badbit);
				is.clear();
				while(true) {
					try {
						double data;
						is >> data;
						datavec.push_back(data);
					} catch(ios::failure &e) {
						is.clear();
						is.exceptions(istream::goodbit);
						is >> ws;
						if(((char) is.peek()) == '#' || ((char) is.peek()) == '/') {
							getline(is, varName);
						} else {
							break;
						}
					}
				}
			} else if(varName == "WEIGHT_MONITOR") {
				weightToMon.clear();
				is.exceptions(istream::failbit | istream::eofbit | istream::badbit);
				is.clear();
				while(true) {
					try {
						unsigned long time, sample;
						is >> time >> sample;
						weightToMon.push_back(pair<unsigned long, unsigned long>(time, sample));
					} catch(ios::failure &e) {
						is.clear();
						is.exceptions(istream::goodbit);
						is >> ws;
						if(((char) is.peek()) == '#' || ((char) is.peek()) == '/') {
							getline(is, varName);
						} else {
							break;
						}
					}
				}
			} else if (varName == "BUILD_FOREST") {
				is >> varName >> ws;
				bForest = (((*varName.begin()) == 'T' || (*varName.begin()) == 't' 
					|| (*varName.begin()) == '1'));
			}
		}
		datamat = matrix(n_num, n_gene, datavec);
		for(vector<pair<unsigned long, unsigned long> >::iterator itor = weightToMon.begin();
			itor != weightToMon.end();) {
				if(itor->first >= n_time || itor->second >= n_Nt[itor->first]) {
					vector<pair<unsigned long, unsigned long> >::iterator itor_to_remove = itor;
					itor++;
					weightToMon.erase(itor_to_remove);
				} else {
					itor++;
				}
		}
}	

int main(int argc, char* argv[])
{   

	long i, j, k, l, n_time, n_gene, n_num, iterations = 100000;
	// The cell number in each time point
	vector <unsigned long> n_Nt;
	// Genes relating to clustering(1) and not relating to clustering(0).
	matrix n_Data;
	bool bForest = false;
	double k0, k1;
	if(argc < 2) {
		cerr << "Usage:" << endl << "./rjmcmc <data file name> [number of iterations]" << endl;
		return 1;
	}

	if(argc > 2) {
		istringstream istr(argv[2]);
		istr >> iterations;
	}

	cout << "Reading input file ... " << flush;
	ifstream fin(argv[1]);
	vector<pair<unsigned long, unsigned long> > weightToMon;
	getInputsFromFile(fin, n_time, n_gene, n_num, n_Nt, n_Data, k1, k0, weightToMon, bForest);
	fin.close();
	cout << "done." << endl;
	vector <bool> n_tao(n_gene, true);
	//vector <bool> newtao(n_gene, false);
	ClusterFlags n_z(n_time, n_Nt);

	// Note that n_z must be set according to tree value.
	for(i=0; i<n_time; i++){
		for(j=0; j<n_Nt.at(i); j++){
		    n_z.setFlag(i,j,1);
		}
	}

	cout << "Initializing parameters ... " << flush;
	initializeParameters(n_gene, n_time, n_Data, k1, k0, bForest);
	cout << "environment ... " << flush;
	MCMCEnv mainenv = MCMCEnv::initMCMCEnv(n_time, n_gene, n_num, n_Nt, n_tao, n_z, n_Data);
	cout << "agent ... " << flush;
	mcmc n_mcmc(mainenv);

	// This is to create a branch from tree 1 at time = 1
	/*ClusterTree & newbranch = n_mcmc.env.createTree(1, 1);
	ClusterTree & newbranch2 = n_mcmc.env.createTree(newbranch.ID, 3);*/

	// This is to create a new tree 
	// ClusterTree & newtree = n_mcmc.env.createTree();

	// to get the ID of the new created tree / branch
	// use newtree.ID / newbranch.ID

	// Edit flag
	/*n_mcmc.env.Flags.setFlag(0,0,1);

	n_mcmc.env.Flags.setFlag(1,0,1);
	n_mcmc.env.Flags.setFlag(1,1,newbranch.ID);

	n_mcmc.env.Flags.setFlag(2,0,1);
	n_mcmc.env.Flags.setFlag(2,1,newbranch.ID);

	n_mcmc.env.Flags.setFlag(3,0,1);
	n_mcmc.env.Flags.setFlag(3,1,newbranch.ID);
	n_mcmc.env.Flags.setFlag(3,2, newbranch2.ID);*/

	vector<long> cnum;
	cnum.push_back(10);
	cnum.push_back(15);
	cnum.push_back(20);
	cnum.push_back(20);

	for(i=0; i<4; i++){
		for(j=0; j<cnum.at(i); j++){
		    n_mcmc.env.Flags.setFlag(i, j, 1);
		}
	}
	
	/*for(int i = 0; i < 20; i++) {
		n_mcmc.env.Flags.setFlag(0, i, 1);
	}

	for(int i = 20; i < 30; i++) {
		n_mcmc.env.Flags.setFlag(1, i, newbranch.ID);
	}
	for(int i = 0; i < 20; i++) {
		n_mcmc.env.Flags.setFlag(1, i, 1);
	}

	for(int i = 20; i < 30; i++) {
		n_mcmc.env.Flags.setFlag(2, i, newbranch.ID);
	}
	for(int i = 30; i < 40; i++) {
		n_mcmc.env.Flags.setFlag(2, i, newbranch2.ID);
	}
	for(int i = 0; i < 20; i++) {
		n_mcmc.env.Flags.setFlag(2, i, 1);
	}

	for(int i = 20; i < 30; i++) {
		n_mcmc.env.Flags.setFlag(3, i, newbranch.ID);
	}
	for(int i = 30; i < 40; i++) {
		n_mcmc.env.Flags.setFlag(3, i, newbranch2.ID);
	}
	for(int i = 0; i < 20; i++) {
		n_mcmc.env.Flags.setFlag(3, i, 1);
	}*/

	// Edit weight
	ClusterTree &root = n_mcmc.env.getTreeFromID(1);

	root.weights[0] = 1;
	root.weights[1] = 1;
	root.weights[2] = 1;
	root.weights[3] = 1;

	/*newbranch.weights[0] = 0;
	newbranch.weights[1] = 0.5;
	newbranch.weights[2] = 0.5;
	newbranch.weights[3] = 0.33;

	newbranch2.weights[0] = 0;
	newbranch2.weights[1] = 0;
	newbranch2.weights[2] = 0;
	newbranch2.weights[3] = 0.34;*/

	// Edit sample number (number of samples in this tree at a certain time)
	root.samples[0] = 10;
	root.samples[1] = 15;
	root.samples[2] = 20;
	root.samples[3] = 20;

	/*newbranch.samples[0] = 0;
	newbranch.samples[1] = 1;
	newbranch.samples[2] = 1;
	newbranch.samples[3] = 1;

	newbranch2.samples[0] = 0;
	newbranch2.samples[1] = 0;
	newbranch2.samples[2] = 0;
	newbranch2.samples[3] = 1;*/
    
	
	n_mcmc.env.Tao[0] = true;
	n_mcmc.env.Tao[1] = true;
	n_mcmc.env.Tao[2] = true;
	n_mcmc.env.Tao[3] = false;
	n_mcmc.env.Tao[4] = false;
	n_mcmc.env.Tao[5] = false;
	n_mcmc.env.Tao[6] = false;
	n_mcmc.env.Tao[7] = false;
	n_mcmc.env.Tao[8] = false;
	n_mcmc.env.Tao[9] = false;
	

	/*cout<<"reData"<<endl;
	n_mcmc.env.reData(newbranch2.ID, 0, n_mcmc.env.Tao, cout);
	cout<<"datacov"<<endl;
	n_mcmc.env.datacov(newbranch2.ID, 0, n_mcmc.env.Tao,cout);*/

	cout << "done." << endl;
	
	

	// Write splitset
	MCMCEnv::TreeSet splitSet = n_mcmc.env.getSplitSet();
	n_mcmc.env.writeSet(cout, splitSet);

	// Write mergeset
	MCMCEnv::TreeMergeSet mergeSet = n_mcmc.env.getMergeSet();
	n_mcmc.env.writeSet(cout, mergeSet);
    
	
	
	
	
	MCMCEnv oldenv(n_mcmc.env);

	// Manually choose a split
	MCMCEnv::TreeSet::value_type chos_splitset = splitSet[3];
	double p_alloc = 1, f_ui = 1;
	bool NULLSet;
	unsigned long long newTreeID;

	// decide which sample will be put into the new tree branch
	typedef pair<unsigned long, unsigned long long> timeAndSample;
	vector<timeAndSample > samples;
	/*samples.push_back(timeAndSample(0, 0));
	samples.push_back(timeAndSample(0, 1));
	samples.push_back(timeAndSample(0, 2));
	samples.push_back(timeAndSample(0, 3));
	samples.push_back(timeAndSample(0, 4));
	samples.push_back(timeAndSample(0, 5));
	samples.push_back(timeAndSample(0, 6));
	samples.push_back(timeAndSample(0, 7));
	samples.push_back(timeAndSample(0, 8));
	samples.push_back(timeAndSample(0, 9));

	samples.push_back(timeAndSample(1, 0));
	samples.push_back(timeAndSample(1, 1));
	samples.push_back(timeAndSample(1, 2));
	samples.push_back(timeAndSample(1, 3));
	samples.push_back(timeAndSample(1, 4));
	samples.push_back(timeAndSample(1, 5));
	samples.push_back(timeAndSample(1, 6));
	samples.push_back(timeAndSample(1, 7));
	samples.push_back(timeAndSample(1, 8));
	samples.push_back(timeAndSample(1, 9));

	samples.push_back(timeAndSample(2, 0));
	samples.push_back(timeAndSample(2, 1));
	samples.push_back(timeAndSample(2, 2));
	samples.push_back(timeAndSample(2, 3));
	samples.push_back(timeAndSample(2, 4));
	samples.push_back(timeAndSample(2, 5));
	samples.push_back(timeAndSample(2, 6));
	samples.push_back(timeAndSample(2, 7));
	samples.push_back(timeAndSample(2, 8));
	samples.push_back(timeAndSample(2, 9));

	samples.push_back(timeAndSample(2, 15));
	samples.push_back(timeAndSample(2, 16));
	samples.push_back(timeAndSample(2, 17));
	samples.push_back(timeAndSample(2, 18));
	samples.push_back(timeAndSample(2, 19));*/
	
	samples.push_back(timeAndSample(3, 15));
	samples.push_back(timeAndSample(3, 16));
	samples.push_back(timeAndSample(3, 17));
	samples.push_back(timeAndSample(3, 18));
	samples.push_back(timeAndSample(3, 19));


	/*samples.push_back(timeAndSample(3, 5));
	samples.push_back(timeAndSample(3, 6));
	samples.push_back(timeAndSample(3, 7));
	samples.push_back(timeAndSample(3, 8));
	samples.push_back(timeAndSample(3, 9));*/

	


	// Notice flagSplitTest should be used to decide the ut and sample
	vector<double> uts;
	uts.push_back(1);
	uts.push_back(1);
	uts.push_back(0.75);
	uts.push_back(0.75);

	newTreeID = n_mcmc.env.flagSplitTest(p_alloc, NULLSet, f_ui, chos_splitset, samples, uts);
	MCMCEnv::TreeSet newSplitset = n_mcmc.env.getSplitSet();
	cout << "APSplit" << endl;
	cout << n_mcmc.env.apSplit(chos_splitset, NULLSet, newSplitset, newSplitset.size(), p_alloc,
		f_ui, oldenv, newTreeID) << endl;

	//n_mcmc.env.writeStatus(cout);

	string temp;
	cin >> temp;

	double probability = 1.0;

	//cerr << n_Data;

	cout<<"old parameters"<<endl;
	n_mcmc.env.writeStatus(cout);

    ofstream outfile( "outfile_sim1.txt" );
	  if ( ! outfile ) {
         cerr << "error: unable to open output file!\n";
         
      }

    ofstream outfile_tree("tree_sim1.txt");
		 if ( ! outfile_tree ) {
         cerr << "error: unable to open output file!\n";
         
      }

    ofstream outfile_treenum("tree_number.txt");
		 if ( ! outfile_treenum ) {
         cerr << "error: unable to open output file!\n";
         
      }

	vector<ofstream*> outfile_weights;
	i = 0;
	for(vector<pair<unsigned long, unsigned long> >::const_iterator itor = weightToMon.begin();
		itor != weightToMon.end(); itor++, i++) {
			ostringstream ostr;
			ostr << "w" << i << "_sim1.txt";
			ofstream *fout = new ofstream(ostr.str().c_str());
			if ( ! *fout ) {
				cerr << "error: unable to open output file!\n";
			}

			outfile_weights.push_back(fout);
	}


    ofstream outfile_branch("Branch_sim1.txt");
		if(! outfile_branch) {
		cerr << "error: unable to open output_branch file!";
		}

	ofstream outfile_taonum("Taonum_sim1.txt");
	    if(! outfile_taonum){
	    cerr << "error: unable to open output_taonum file!";
	    }
	bool move = false;

	for(j=0; j<iterations; j++){

		n_mcmc.env.sampleTao();
		n_mcmc.env.sampleWeight();
		n_mcmc.env.sampleZ();
		if(!(j % 10)) {
			cout<<"\rSTEP: "<<j+1<<' ' << n_mcmc.env.getClusterNumber() << ' '
				<< n_mcmc.env.treeSummary() << flush;
		}
		outfile<<"STEP: "<<j+1<<std::endl;
		n_mcmc.env.writeStatus(outfile);
		outfile_taonum<< n_mcmc.env.taoCount() <<std::endl;

		//cout<<"After SplitMerge_move!"<<std::endl;
		outfile<<"After SplitMerge_move!"<<std::endl;
		probability = 0.0;
		move = n_mcmc.SplitMerge_move(probability);
		i = 0;
		for(vector<pair<unsigned long, unsigned long> >::const_iterator itor = weightToMon.begin();
			itor != weightToMon.end(); itor++, i++) {
				*(outfile_weights[i]) 
				<< n_mcmc.env.getWeightFromSample(itor->first, itor->second) << std::endl;
		}


		outfile<<"move: "<<move<<std::endl;

		if(move){
			outfile<<"SplitSplitSplit:"<<endl;
			outfile<<"accepted probability:"<<endl;
		}else{
			outfile<<"MergeMergeMerge:"<<endl;
			outfile<<"accetped probability:"<<endl;
		}

		outfile<<probability<<endl;
		n_mcmc.env.writeStatus(outfile);

		//cout<<"After BirthDeath_move!"<<std::endl;
		outfile<<"After BirthDeath_move!"<<std::endl;

		probability = 0.0;
		move = n_mcmc.BirthDeath_move(probability);

		if(move){
			outfile<<"BirthBirthBirth:"<<std::endl;
			outfile<<"accepted probability:"<<endl;
		}else{
			outfile<<"DeathDeathDeath:"<<std::endl;
			outfile<<"accepted probability:"<<endl;
		}

		outfile<<probability<<endl;
		n_mcmc.env.writeStatus(outfile);

		outfile_branch << n_mcmc.env.getClusterNumber() << std::endl;
		outfile_treenum << n_mcmc.env.getTreeNumbers() << endl;
		outfile_tree<<n_mcmc.env.treeSummary()<<endl;

	}

	for(vector<ofstream*>::iterator itor = outfile_weights.begin();
		itor != outfile_weights.end(); itor++) {
			(*itor)->close();
			delete(*itor);
			*itor = NULL;
	}

	cout << endl;

	return 0;
}




