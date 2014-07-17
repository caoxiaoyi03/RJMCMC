#include "MCMCEnv.h"
#include "Random.h"
#include <stdexcept>
#include <deque>
#include <iomanip>
#include <cmath>
#include <sstream>

#define SMALLNUM 1e-200
#define TREE_ID_SIZE	4

unsigned long MCMCEnv::NumOfTP = 0;
unsigned long MCMCEnv::NumOfGenes = 0;
unsigned long MCMCEnv::NumOfData = 0;
ClusterTree::NumOfSampleContainer MCMCEnv::DataNumEachTP;	// number of data in each time point

matrix MCMCEnv::mu1;
matrix MCMCEnv::mu0;
matrix MCMCEnv::sigma0;
matrix MCMCEnv::sigma1;
double MCMCEnv::beta0;
double MCMCEnv::beta1;
double MCMCEnv::deltar;
double MCMCEnv::alpha;
double MCMCEnv::shape;
double MCMCEnv::bk; // For split_merge move and birth_death move, we give equal probability.
double MCMCEnv::ranc; // Tail birth is given much more porbability.
double MCMCEnv::k1;
double MCMCEnv::k0;
double MCMCEnv::constSigma0;
double MCMCEnv::constSigma1;
bool MCMCEnv::buildForest;

MCMCEnv::SortedULContainer MCMCEnv::TPStartIndex,
	MCMCEnv::ReverseTPStartIndex;	// the starting index for each time point

using namespace std;

// matrix calculation
// In this function, cluster indicator z is assumbed to be reflaged using 'reflag' function.
// tao_NULL = 1, delete tao=0, otherwise delete tao=1.
// the return matrix is the shrinked version of the data matrix
// taoIndices: the indices for tao == 1
matrix MCMCEnv::reData(unsigned long long ID, bool tao_NULL, const TaoSet &newTao, vector<unsigned long> &taoIndices) const {
			// the indices of genes that will go to the result matrix
	unsigned long index = 0;
	taoIndices.clear();
	for(TaoSet::const_iterator itor = newTao.begin(); itor != newTao.end(); itor++, index++) {
		if(tao_NULL == *itor) {
			taoIndices.push_back(index);
		}
	}
	matrix result(taoIndices.size(), true, NumOfData);
	//std::vector<std::vector<long> > mclust = Flags.onezerocluster(ID);
	unsigned long length = 0;
	for(unsigned long i = 0; i < Flags.rows(); i++){
		for(unsigned long j = 0; j < Flags.cols(i); j++) {
			if(Flags.getFlag(i, j) == ID || !tao_NULL) {
				dvector resultvec;
				for(vector<unsigned long>::const_iterator itor = taoIndices.begin();
					itor != taoIndices.end(); itor++) {
						resultvec.push_back(Data.getValue(length + j, *itor));
				}
				result.addRow(resultvec);
			}
			//for(k=0; k<one.cols(); k++) {
			//	one.setValue(length + j, k, Data.getValue(length + j, k) 
			//		* ((tao_NULL == newTao[k])? (newTao[k]? mclust.at(i).at(j): 1): 0));
			//}
		} 

		length += Flags.cols(i);
	}
	//cout << "reData" << endl << one;
	return result;

}

matrix MCMCEnv::reData(unsigned long long ID, bool tao_NULL, const TaoSet &newTao, ostream &os) const {
	vector<unsigned long> taoIndices;
	matrix one = reData(ID, tao_NULL, newTao, taoIndices);
	for(TaoSet::const_iterator itor = newTao.begin();
		itor != newTao.end(); itor++) {
			os << *itor << ' ';
	}
	os << endl << "reData" << endl << one;
	return one;
}

matrix MCMCEnv::datacov(unsigned long long flag, bool tao_NULL, const TaoSet &newTao, ostream &os) const {
	matrix result = datacov(flag, tao_NULL, newTao);
	os << "Datacov" << endl << result;
	return result;
}

// Notice that when tao_NULL = true, flag is not used in the function
// so any data will be considered
matrix MCMCEnv::datacov(unsigned long long flag, bool tao_NULL, const TaoSet &newTao) const {
	long x = (tao_NULL? Flags.flagnum(flag): NumOfData);
	vector<unsigned long> taoIndices;
	matrix subData = reData(flag, tao_NULL, newTao, taoIndices);
	//matrix two = (tao_NULL? mu1: mu0);
	unsigned long cols = taoIndices.size();
	
	matrix two(1, cols);
	for(unsigned long k = 0; k < cols; k++){
		two(0, k) = mu1(0, taoIndices[k]);
	}

	matrix one_mean = (subData.datasum() * (1.0 / (double) x));
	//TimeCoordinate xx;

	matrix cov_NULL(cols, cols);
	int i,j;

	for(i = 0; i < subData.rows(); i++)
	{
		//xx = IndexToCoor(i);
		matrix mid(1, cols, subData.getRowVec(i));
		//if((Flags.getFlag(xx.first, xx.second) == flag) || !tao_NULL){
		mid -= one_mean;
		//} else {
		//	mid = one.getRow(i);
		//}

		//std::cout<<"Matrix transpose: "<<std::endl;
		//std::cout<<mid<<std::endl;
		cov_NULL += mid.transTimesSelf();
	}

	two -= one_mean;
	/*std::cout<<"mu0: "<<std::endl;
	std::cout<<two<<std::endl;
	std::cout<<"mean: "<<std::endl;
	std::cout<<one_mean<<std::endl;*/
	return cov_NULL.rbind(two.transTimesSelf());

}

// Caution: if the index is out of bound, will return NumOfTimePoints
MCMCEnv::TimeCoordinate MCMCEnv::IndexToCoor(unsigned long index) const {
	SortedULContainer::const_iterator itor = TPStartIndex.upper_bound(index);
	itor--;
	return TimeCoordinate(itor->second, index - itor->first);
}

unsigned long MCMCEnv::CoorToIndex(const TimeCoordinate &coor) const {
	return ReverseTPStartIndex[coor.first] + coor.second;
}


MCMCEnv::MCMCEnv(const TaoSet &tao, const ClusterFlags &flags, const matrix &data)
	: ID_max(1), rootIDs(RootIDContainer()), Flags(flags), Tao(tao), UpToDate(false), Data(data)
{
	// actually head node need to be initialized here
	ClusterTree &ctree = createTree();
	//cout << "Tree" << flush;
	ctree.fillWeights(1.0);
	//cout << "Tree" << flush;
	ctree.samples = DataNumEachTP;
	
}

MCMCEnv MCMCEnv::initMCMCEnv(unsigned long numOfTimePoints, unsigned long numOfGenes,
	unsigned long numOfData, const std::vector<unsigned long> &numOfDataInEachTimePoint,
	const TaoSet &tao, const ClusterFlags &flags, const matrix &p_data) {
		NumOfTP = numOfTimePoints;
		NumOfGenes = numOfGenes;
		NumOfData = numOfData;
		DataNumEachTP.resize(NumOfTP);
		DataNumEachTP.assign(numOfDataInEachTimePoint.begin(), numOfDataInEachTimePoint.end());

		unsigned long index = 0, time = 0;
		TPStartIndex.clear();
		ReverseTPStartIndex.clear();

		for(ClusterTree::NumOfSampleContainer::const_iterator itor = DataNumEachTP.begin();
			itor != DataNumEachTP.end(); itor++) {
				TPStartIndex.insert(SortedULContainer::value_type(index, time));
				ReverseTPStartIndex.insert(SortedULContainer::value_type(time, index));
				index += *itor;
				time++;
		}
		TPStartIndex.insert(SortedULContainer::value_type(index, time));

		//cout << NumOfTP << flush;
		return MCMCEnv(tao, flags, p_data);
}

MCMCEnv::MCMCEnv(const MCMCEnv &oldenv): MapClusterTrees(oldenv.MapClusterTrees),
	ID_max(oldenv.ID_max), rootIDs(oldenv.rootIDs), Flags(oldenv.Flags), Tao(oldenv.Tao),
	UpToDate(oldenv.UpToDate), LogDensityValue(oldenv.LogDensityValue), Data(oldenv.Data) {
		// NOTICE: I'm not exactly sure if shallow copy of MapClusterTrees works, 
		//         so need verifying
}

MCMCEnv &MCMCEnv::operator=(const MCMCEnv &oldenv) {
	MapClusterTrees = oldenv.MapClusterTrees;
	ID_max = oldenv.ID_max;
	rootIDs = oldenv.rootIDs;
	Flags = oldenv.Flags;
	Tao = oldenv.Tao;
	UpToDate = oldenv.UpToDate;
	LogDensityValue = oldenv.LogDensityValue;
	// doesn't need to change data
	return *this;
}

MCMCEnv::~MCMCEnv(void)
{
	//// when environment is destoryed, destory all the trees
	//for(TreeContainer::const_iterator itor = MapClusterTrees.begin();
	//	itor != MapClusterTrees.end(); itor++) {
	//		if(itor->second) {
	//			delete itor->second;
	//			//itor->second = NULL;
	//		}
	//}
}

unsigned long MCMCEnv::getDataNum() {
	return NumOfData;
}

void MCMCEnv::initializeParameters(const matrix &p_mu1, const matrix &p_mu0, const matrix &p_sigma1,
	const matrix &p_sigma0, double p_beta1, double p_beta0, double p_deltar, double p_alpha,
	double p_shape, double p_bk, double p_ranc, bool bForest) {
		mu1 = p_mu1;
		mu0 = p_mu0;
		sigma1 = p_sigma1;
		sigma0 = p_sigma0;
		beta1 = p_beta1;
		beta0 = p_beta0;
		deltar = p_deltar;
		alpha = p_alpha;
		shape = p_shape;
		bk = p_bk;
		ranc = p_ranc;
		constSigma0 = sigma0(0, 0);
		constSigma1 = sigma1(0, 0);
		buildForest = bForest;
}


//ClusterTree &MCMCEnv::getTreeFromID(unsigned long long ID) {
//	if(MapClusterTrees.find(ID) == MapClusterTrees.end()) {
//		ostringstream ostr;
//		ostr << "ID " << ID << " not in the list!";
//		throw(std::logic_error(ostr.str()));
//	}
//	return MapClusterTrees.find(ID)->second;
//}
//
//const ClusterTree &MCMCEnv::getTreeFromID(unsigned long long ID) const {
//	if(MapClusterTrees.find(ID) == MapClusterTrees.end()) {
//		ostringstream ostr;
//		ostr << "ID " << ID << " not in the list!";
//		throw(std::logic_error(ostr.str()));
//	}
//	return MapClusterTrees.find(ID)->second;
//}
//
unsigned long MCMCEnv::getTreeNumbers() const {
	return rootIDs.size();
}

ClusterTree &MCMCEnv::createTree(unsigned long long parentID, unsigned long time) {
	unsigned long long newID = ID_max++;
	ClusterTree &parent = getTreeFromID(parentID);
	//if(parent.getChildID(time)) {
	//	// there is child here, so throw error
	//	throw(std::logic_error("Already have a child at the time!"));
	//}
	parent.setChild(newID, time);
	MapClusterTrees.insert(TreeContainer::value_type(newID, ClusterTree(parent, newID, time)));
	UpToDate = false;
	return getTreeFromID(newID);
}

ClusterTree &MCMCEnv::createTree() {
	// creating mother
	//if(rootIDs.size()) {
	//	throw(std::logic_error("Mother node already created!"));
	//}
	unsigned long long currRootID = ID_max++;
	rootIDs.insert(RootIDContainer::value_type(currRootID, 0));
	MapClusterTrees.insert(TreeContainer::value_type(currRootID, ClusterTree(currRootID)));
	UpToDate = false;
	return getTreeFromID(currRootID);
}

void MCMCEnv::removeTree(unsigned long long TreeID) {
	if(getTreeFromID(TreeID).getSampleNum()) {
		writeTree(cerr);
		cerr << "Tree to remove: " << TreeID << endl;
		for(unsigned long t = 0; t < NumOfTP; t++) {
			cerr << getTreeFromID(TreeID).weights[t] << ' ';
		}
		for(unsigned long t = 0; t < NumOfTP; t++) {
			cerr << getTreeFromID(TreeID).samples[t] << ' ';
		}
		cerr << endl;
		cerr << "Flags" << endl << Flags << endl;
		throw(std::logic_error("Tree has data points, cannot be removed."));
	}
	// check whether the tree is a root
	if(!getTreeFromID(TreeID).getBornTime()) {
		// is a root
		if(rootIDs.size() <= 1) {
			throw(std::logic_error("Last tree here, cannot be removed."));
		}
		rootIDs.erase(TreeID);
	}
	MapClusterTrees.erase(TreeID);
	UpToDate = false;
}

//ClusterTree &MCMCEnv::getTreeHead() {
//	return MapClusterTrees.find(root_ID)->second;
//}
//
//const ClusterTree &MCMCEnv::getTreeHead() const {
//	return MapClusterTrees.find(root_ID)->second;
//}

MCMCEnv::TreeSet MCMCEnv::getSplitSet(bool hasForest, bool forestOnly) const {

	MCMCEnv::TreeSet result;
	result.reserve(NumOfTP * 2);

	if(buildForest && hasForest) {
	// get the new tree split set
		for(RootIDContainer::const_iterator itor = rootIDs.begin();
			itor != rootIDs.end(); itor++) {
				if(getTreeFromID(itor->first).samples[0]) {
					result.push_back(TreeSet::value_type(0, itor->first));
				}
		}
	}
	
	//cout << "Generating Split Set ... " << endl;
	if(!forestOnly) {
		for(TreeContainer::const_iterator itor = MapClusterTrees.begin();
			itor != MapClusterTrees.end(); itor++) {
				//cout << "ID: " << itor->second.ID << " Parent ID: " 
				//	<< itor->second.parentID 
				//	<< " Born time: " << itor->second.getBornTime() << endl;
				for(unsigned long t = itor->second.getBornTime() + 1; t < NumOfTP; t++) {
					// the born time of its child must be at least 1 more than this
					//cout << "Samples(" << t << "): " << itor->second.samples[t] << endl;
					if(itor->second.samples[t] && !itor->second.children[t]) {
						// has sample and no child at time t
						result.push_back(MCMCEnv::TreeSet::value_type(t, itor->first));
					}
				}
		}
	}
	return result;
	//long i,j;
	//lvector splitone;
	////Collect nonempty cluster vector.
	//std::vector<long> flags = newz.getTtoTindicator(n_timepoints);

	//long midf,midt,mids,sr = 0;
	//long fl = flags.size();
	//for(i=0; i<fl; i++){
	//	midf = flags.at(i);
	//	midt = newz.getEarliestTime(midf);
	//	

	//	for(j=1; j<n_timepoints; j++){ 

	//		mids = (j==midt-1)?(newchome(0,midf-1)-1):(midf-1);

	//		if((newtree(j-1,mids)==0)&&(newz.flaginT(midf,j+1))){
	//		   splitone.push_back(j+1);
	//		   splitone.push_back(midf);
	//		   sr++;
	//		}
	//	}
	//	
	//}


	//return matrixl(sr,2,splitone);

}

MCMCEnv::TreeMergeSet MCMCEnv::getMergeSet(bool hasForest, bool forestOnly) const {
	//std::vector<long> flags = newz.getTtoTindicator(n_timepoints);
	//std::vector<long> freeflags;
	//std::vector<long> nonfreeflags;
	//lvector mergeone;
	//long i,j,k,s=0;
	//long time1, time2, time0;
	//long sr=0;
	//long n,m,l;

	unsigned long long childID, parentID;

	TreeMergeSet result;
	result.reserve(NumOfTP * 2);

	if(buildForest && hasForest) {
		for(RootIDContainer::const_iterator itor = rootIDs.begin();
			itor != rootIDs.end(); itor++) {
				if(!getTreeFromID(itor->first).getNumOfChildren()) {
					// root does not have any children
					for(RootIDContainer::const_iterator itor_parent = rootIDs.begin();
						itor_parent != rootIDs.end(); itor_parent++) {
							if(itor_parent->first != itor->first) {
								result.push_back(TreeMergeSet::value_type(itor_parent->first,
									itor->first));
							}
					}

				}
		}
	}

	if(!forestOnly) {
		for(TreeContainer::const_iterator itor = MapClusterTrees.begin();
			itor != MapClusterTrees.end(); itor++) {
				// first check if the tree doesn't have any children
				if(!itor->second.getNumOfChildren()) {
					// this tree does not have any children
					childID = itor->second.ID;
					for(TreeContainer::const_iterator itor_parent = MapClusterTrees.begin();
						itor_parent != MapClusterTrees.end(); itor_parent++) {
							if(itor_parent->second.ID != itor->second.ID 
								&& itor_parent->second.getBornTime() < itor->second.getBornTime()) {
									// this parent is OK
									result.push_back(TreeMergeSet::value_type(itor_parent->second.ID,
										itor->second.ID));
							}
					}
				}
		}
	}

	return result;

}

MCMCEnv::TreeMergeSet MCMCEnv::getDeathSet(bool hasForest, bool forestOnly) const {
	TreeMergeSet result;
	result.reserve(NumOfTP * 2);
	
	unsigned long long childID, parentID;

	if(buildForest && hasForest) {
		for(RootIDContainer::const_iterator itor = rootIDs.begin();
			itor != rootIDs.end(); itor++) {
				if(!getTreeFromID(itor->first).getNumOfChildren() && !getTreeFromID(itor->first).getSampleNum()) {
					// root does not have any children and empty
					for(RootIDContainer::const_iterator itor_parent = rootIDs.begin();
						itor_parent != rootIDs.end(); itor_parent++) {
							if(itor_parent->first != itor->first) {
								result.push_back(TreeMergeSet::value_type(itor_parent->first,
									itor->first));
							}
					}

				}
		}
	}

	if(!forestOnly) {
		for(TreeContainer::const_iterator itor = MapClusterTrees.begin();
			itor != MapClusterTrees.end(); itor++) {
				// first check if the tree doesn't have any children
				if(!itor->second.getNumOfChildren()) {
					// this tree does not have any children
					if(!itor->second.getSampleNum()) {
						// this tree is empty
						childID = itor->second.ID;
						for(TreeContainer::const_iterator itor_parent = MapClusterTrees.begin();
							itor_parent != MapClusterTrees.end(); itor_parent++) {
								if(itor_parent->second.ID != itor->second.ID 
									&& itor_parent->second.getBornTime() < itor->second.getBornTime()) {
										// this parent is OK
										result.push_back(TreeMergeSet::value_type(itor_parent->second.ID,
											itor->second.ID));
								}
						}
					}
				}
		}
	}

	return result;
}

MCMCEnv::TreeSet MCMCEnv::getTailBirthSet(bool hasForest, bool forestOnly) const {
	// This actually has nothing to do with forests so both flags will not be used at all.

	MCMCEnv::TreeSet result;
	result.reserve(NumOfTP * 2);

	//cout << "Generating Split Set ... " << endl;
	for(TreeContainer::const_iterator itor = MapClusterTrees.begin();
		itor != MapClusterTrees.end(); itor++) {
			//cout << "ID: " << itor->second.ID << " Parent ID: " 
			//	<< itor->second.parentID 
			//	<< " Born time: " << itor->second.getBornTime() << endl;
			if(itor->second.getDeathTime() < NumOfTP) {
				// can birth
				result.push_back(MCMCEnv::TreeSet::value_type(itor->second.getDeathTime(), itor->first));
			}
	}
	return result;
}

MCMCEnv::TreeSet MCMCEnv::getTailDeathSet(bool hasForest, bool forestOnly) const {
	// This actually has nothing to do with forests so both flags will not be used at all.

	MCMCEnv::TreeSet result;
	result.reserve(NumOfTP * 2);

	//cout << "Generating Split Set ... " << endl;
	for(TreeContainer::const_iterator itor = MapClusterTrees.begin();
		itor != MapClusterTrees.end(); itor++) {
			//cout << "ID: " << itor->second.ID << " Parent ID: " 
			//	<< itor->second.parentID 
			//	<< " Born time: " << itor->second.getBornTime() << endl;
			if(itor->second.getDeathTime() >= NumOfTP && itor->second.getEarliestDeathableTime() < NumOfTP) {
				// can death
				for(unsigned long t = itor->second.getEarliestDeathableTime(); t < itor->second.getDeathTime(); t++) {
					result.push_back(MCMCEnv::TreeSet::value_type(t, itor->first));
				}
			}
	}
	return result;
}

unsigned long long MCMCEnv::flagSplit(double &P_alloc, bool &NULLset, double &f_ui,  
	const MCMCEnv::TreeSet::value_type &split_pair) {
		ClusterTree &parent = getTreeFromID(split_pair.second);
		double ut;
		f_ui = 1.0;
		long s = 0;
		unsigned long long childID;
		if(parent.samples[split_pair.first] && split_pair.first < NumOfTP) {
			ClusterTree &child = (split_pair.first)? 
				createTree(split_pair.second, split_pair.first): createTree();
			childID = child.ID;
			for(unsigned long i = split_pair.first; i < NumOfTP; i++) {

				/*weight*/
				ut = ran_beta(2.0,2.0);
				double old_weight = parent.weights[i];
				f_ui *= (old_weight > SMALLNUM)? betad(ut, 2.0, 2.0): 1.0;
				/* This is a good choice for resetting the cluster indicators and weight parameter */
				parent.weights[i] = ut * old_weight;
				child.weights[i] = (1 - ut) * old_weight;

				// flags
				if(parent.samples[i]) {
					// there are samples in this time point for the parent
					for(long k = 0; k < Flags.cols(i); k++){
						if(Flags(i,k) == parent.ID) {
							Flags(i,k) = ran_nber(parent.ID, child.ID, ut);
							if(Flags(i,k) == parent.ID){
								//one(i-1,k) = one(i-1,k);
								P_alloc *= ut;
								// nothing is changed in terms of sample number
							} else {
								P_alloc *= (1-ut);
								(parent.samples[i])--;
								(child.samples[i])++;
								s++;
							}			
						}

						// cluster indicator reloaded.
					}
				}
			}

			NULLset = (s == 0);
			/*std::cout<<"new"<<std::endl;
			std::cout<<one<<std::endl;
			std::cout<<two<<std::endl;
			std::cout<<newtree<<std::endl;
			std::cout<<newchome<<std::endl;
			std::cout<<P_alloc<<std::endl;
			std::cout<<NULLset<<std::endl;*/
			UpToDate = false;

		} else {
			std::cout<<"This is an empty cluster!"<<std::endl;
			std::cout<<"This cluster can not be splitted!"<<std::endl;
		}
		return childID;
}

void MCMCEnv::flagSplitTest(double &P_alloc, bool &NULLset, double &f_ui,  
	const MCMCEnv::TreeSet::value_type &split_pair) {
		ClusterTree &parent = getTreeFromID(split_pair.second);
		double ut;
		f_ui = 1.0;
		long s = 0;
	
		if((split_pair.first) && (split_pair.first < NumOfTP) && 
			(parent.samples[split_pair.first])) {
				// first create a new tree if the split_pair
				ClusterTree &child = createTree(split_pair.second, split_pair.first);
				for(unsigned long i = split_pair.first; i < NumOfTP; i++) {

					/*weight*/
					ut = ran_beta(2.0,2.0);
					double old_weight = parent.weights[i];
					f_ui *= 1.0;
					/* This is a good choice for resetting the cluster indicators and weight parameter */

					// flags
					parent.samples[i] = 0;
					child.samples[i] = 0;
					//if(parent.samples[i]) {
					//	// there are samples in this time point for the parent
						for(long k = 0; k < Flags.cols(i); k++){
							if(Flags(i,k) == parent.ID) {
								////Flags(i,k) = ran_nber(parent.ID, child.ID, ut);
								//if(Flags(i,k) == parent.ID){
								//	//one(i-1,k) = one(i-1,k);
								//	P_alloc *= ut;
								//	// nothing is changed in terms of sample number
								//} else {
									//P_alloc *= (1-ut);
									(parent.samples[i])++;
							} else if(Flags(i, k) == child.ID) {
									(child.samples[i])++;
									//s++;
								//}			
							}

							// cluster indicator reloaded.
						}
					//}
					parent.weights[i] = (double) parent.samples[i] 
						/ (double) (parent.samples[i] + child.samples[i]);
					child.weights[i] = (double) child.samples[i] 
						/ (double) (parent.samples[i] + child.samples[i]);
				}

	
				NULLset = false;
				/*std::cout<<"new"<<std::endl;
				std::cout<<one<<std::endl;
				std::cout<<two<<std::endl;
				std::cout<<newtree<<std::endl;
				std::cout<<newchome<<std::endl;
				std::cout<<P_alloc<<std::endl;
				std::cout<<NULLset<<std::endl;*/
				UpToDate = false;

		} else {
			std::cout<<"This is an empty cluster!"<<std::endl;
			std::cout<<"This cluster can not be splitted!"<<std::endl;
		}
}


void MCMCEnv::flagMerge(double &P_alloc, double &f_ui, const TreeMergeSet::value_type &merge_pair) {

	unsigned long parentNum, childNum;
	double ut;

	P_alloc = 1.0;
	f_ui = 1.0;
	ClusterTree &parent = getTreeFromID(merge_pair.first),
		&child = getTreeFromID(merge_pair.second);

	unsigned long time = child.getBornTime();
	if(parent.getDeathTime() < child.getDeathTime()) {
		parent.setDeathTime(child.getDeathTime());
	}
	for(unsigned long t = time; t < parent.getDeathTime(); t++) {

		if(child.weights[t] > SMALLNUM && parent.weights[t] > SMALLNUM) {
			ut = parent.weights[t] / (parent.weights[t] + child.weights[t]);
			parentNum = Flags.tflagnum(parent.ID, t);
			childNum = Flags.tflagnum(child.ID, t);
			P_alloc *= pow(ut, (int) parentNum) * pow(1 - ut, (int) childNum);
			f_ui *= betad(ut,2.0,2.0);
			parent.weights[t] += child.weights[t];
			// calc P_alloc and change weights and sample count
		}

		for(unsigned long k = 0; k < Flags.cols(t); k++){
			if(Flags(t, k) == child.ID) {
				Flags(t, k) = parent.ID;
			}
			// change cluster indicator
		}
		parent.samples[t] += child.samples[t];
		child.samples[t] = 0;
	}
	if(child.getBornTime() > 0) {
		// not a root
		getTreeFromID(child.parentID).removeChild(time);
	}
	// remove child from its original parent
	//child.parentID = 0;			// reset the parentID of the child, mark it for deletion
	try {
		removeTree(child.ID);
	} catch(std::logic_error &e) {
		cerr << "Merge at time=" << time << endl;
		throw e;
	}
	UpToDate = false;

}

void MCMCEnv::flagBirth(const TreeSet::value_type &splitPair) {
	double weightratio, jacobi, f_wstar;
	return flagBirth(splitPair, weightratio, jacobi, f_wstar);
}

void MCMCEnv::flagBirth(const TreeSet::value_type &splitPair,
	double &weightratio, double &jacobi, double &f_wstar) {
		long j,k,n,timec;
		double w_star;

		jacobi = 1.0;
		f_wstar = 1.0;
		weightratio = 1.0;
		unsigned long time = splitPair.first;

		const ClusterTree &parent = getTreeFromID(splitPair.second);

		if((time < NumOfTP) && parent.getSampleNum(time)) {
			/*Justify whether clusterindex is an empty cluster.*/

			ClusterTree &child = (time? createTree(parent.ID, time): createTree());

			for(unsigned long t = time; t < child.getDeathTime(); t++){
				/*weight*/
				timec = getNoZeroWeights(t);
				w_star = ran_beta(1.0, (double) timec);
				jacobi *= pow(1 - w_star, timec - 1);
				f_wstar *= betad(w_star, 1.0, (double) timec);
				weightratio *= (pow(w_star, alpha - 1) * pow(1 - w_star, 
					(double)timec * (alpha - 1)) / betaf(alpha, (double) timec * alpha));

				/* This is a good choice for resetting the cluster indicators and weight parameter */
				//two(i-1,clusterindex-1) = (1.0 - w_star) * weight(i-1,clusterindex-1);
				child.weights[t] = w_star;

				for(TreeContainer::iterator itor = MapClusterTrees.begin();
					itor != MapClusterTrees.end(); itor++) {
						if(itor->first != child.ID){
							itor->second.weights[t] = (1.0 - w_star) * itor->second.weights[t];
						}
				}
				/*weight relocated.*/
			}


			// Tree reset
			// This step is implemented according to the weight setting.

			/*std::cout<<"Birth"<<std::endl;
			std::cout<<"weightratio"<<std::endl;
			std::cout<<weightratio<<std::endl;
			std::cout<<"Jacobi"<<std::endl;
			std::cout<<jacobi<<std::endl;
			std::cout<<"f_wstar"<<std::endl;
			std::cout<<f_wstar<<std::endl;*/
			UpToDate = false;

		} else {
			std::cout<<"This is an empty cluster!"<<std::endl;
			std::cout<<"This cluster can not be splitted!"<<std::endl;
		}

}

void MCMCEnv::flagDeath(const TreeMergeSet::value_type &merge_pair, 
	double &jacobi, double &weightratio, double &f_wstar) {

		long mergeflag, leftflag;
		long mergenum, leftnum;
		long i,j,k;
		double midw,midn;

		jacobi = 1.0;
		weightratio = 1.0;
		f_wstar = 1.0;

		ClusterTree &parent = getTreeFromID(merge_pair.first),
			&child = getTreeFromID(merge_pair.second);

		unsigned long time = child.getBornTime();
		//// Mother
		//leftflag = clusterindex1;
		//// Son
		//mergeflag = clusterindex2;
		for(unsigned long t = time; t < child.getDeathTime(); t++){
			midw = child.weights[t];
			if(midw > SMALLNUM){
				midn = getNoZeroWeights(t);
				f_wstar *= betad(midw, 1.0, (double) midn);
				weightratio *= (betaf(alpha, (double) (midn - 1) * alpha)
					/ pow(midw, (alpha - 1.0)) 
					/ pow(1.0 - midw, (double) (midn - 1) * (alpha - 1.0)));
				jacobi *= pow(1.0 - midw, -1.0 * (double) (midn - 2));

				for(TreeContainer::iterator itor = MapClusterTrees.begin();
					itor != MapClusterTrees.end(); itor++) {
						if(itor->first != child.ID){
							itor->second.weights[t] = itor->second.weights[t] / (1.0 - midw);
						}
				}
			}
			// cluster indicator reloaded.
		}

		if(child.getBornTime()){
			getTreeFromID(child.parentID).removeChild(time);
		}
		// remove child from its original parent
		//child.parentID = 0;			// reset the parentID of the child, mark it for deletion
		try {
			removeTree(child.ID);
		} catch(std::logic_error &e) {
			cerr << "Death" << endl;
			throw e;
		}
		UpToDate = false;

		/*std::cout<<"Death"<<std::endl;
		std::cout<<"weightratio"<<std::endl;
		std::cout<<weightratio<<std::endl;
		std::cout<<"Jacobi"<<std::endl;
		std::cout<<jacobi<<std::endl;
		std::cout<<"f_wstar"<<std::endl;
		std::cout<<f_wstar<<std::endl;*/

}

void MCMCEnv::flagTailBirth(const TreeSet::value_type &splitPair,
	double &weightratio, double &jacobi, double &f_wstar) {
		long j,k,n,timec;
		double w_star;

		jacobi = 1.0;
		f_wstar = 1.0;
		weightratio = 1.0;

		// technically time should just be parent.getDeathTime();

		ClusterTree &parent = getTreeFromID(splitPair.second);
		unsigned long time = parent.getDeathTime();

		parent.tailBirth();

		if(time < NumOfTP) {
			/*Justify whether clusterindex is an empty cluster.*/

			//ClusterTree &child = (time? createTree(parent.ID, time): createTree());

			for(unsigned long t = time; t < NumOfTP; t++){
				/*weight*/
				timec = getNoZeroWeights(t);
				w_star = ran_beta(1.0, (double) timec);
				jacobi *= pow(1 - w_star, timec - 1);
				f_wstar *= betad(w_star, 1.0, (double) timec);
				weightratio *= (pow(w_star, alpha - 1) * pow(1 - w_star, 
					(double)timec * (alpha - 1)) / betaf(alpha, (double) timec * alpha));

				/* This is a good choice for resetting the cluster indicators and weight parameter */
				//two(i-1,clusterindex-1) = (1.0 - w_star) * weight(i-1,clusterindex-1);
				parent.weights[t] = w_star;

				for(TreeContainer::iterator itor = MapClusterTrees.begin();
					itor != MapClusterTrees.end(); itor++) {
						if(itor->first != parent.ID){
							itor->second.weights[t] = (1.0 - w_star) * itor->second.weights[t];
						}
				}
				/*weight relocated.*/
			}


			// Tree reset
			// This step is implemented according to the weight setting.

			/*std::cout<<"Birth"<<std::endl;
			std::cout<<"weightratio"<<std::endl;
			std::cout<<weightratio<<std::endl;
			std::cout<<"Jacobi"<<std::endl;
			std::cout<<jacobi<<std::endl;
			std::cout<<"f_wstar"<<std::endl;
			std::cout<<f_wstar<<std::endl;*/
			UpToDate = false;

		} else {
			std::cout<<"This is an empty cluster!"<<std::endl;
			std::cout<<"This cluster can not be splitted!"<<std::endl;
		}

}

void MCMCEnv::flagTailDeath(const TreeSet::value_type &death_set, 
	double &jacobi, double &weightratio, double &f_wstar) {

		long mergeflag, leftflag;
		long mergenum, leftnum;
		long i,j,k;
		double midw,midn;

		jacobi = 1.0;
		weightratio = 1.0;
		f_wstar = 1.0;

		ClusterTree &parent = getTreeFromID(death_set.second);

		unsigned long time = death_set.first;
		//// Mother
		//leftflag = clusterindex1;
		//// Son
		//mergeflag = clusterindex2;
		for(unsigned long t = time; t < parent.getDeathTime(); t++){
			midw = parent.weights[t];
			if(midw > SMALLNUM){
				midn = getNoZeroWeights(t);
				f_wstar *= betad(midw, 1.0, (double) midn);
				weightratio *= (betaf(alpha, (double) (midn - 1) * alpha)
					/ pow(midw, (alpha - 1.0)) 
					/ pow(1.0 - midw, (double) (midn - 1) * (alpha - 1.0)));
				jacobi *= pow(1.0 - midw, -1.0 * (double) (midn - 2));

				for(TreeContainer::iterator itor = MapClusterTrees.begin();
					itor != MapClusterTrees.end(); itor++) {
						if(itor->first != parent.ID){
							itor->second.weights[t] = itor->second.weights[t] / (1.0 - midw);
						}
				}
			}
			// cluster indicator reloaded.
		}

		// remove child from its original parent
		//child.parentID = 0;			// reset the parentID of the child, mark it for deletion
		try {
			parent.tailDeath(time);
		} catch(std::logic_error &e) {
			cerr << "Death" << endl;
			throw e;
		}
		UpToDate = false;

		/*std::cout<<"Death"<<std::endl;
		std::cout<<"weightratio"<<std::endl;
		std::cout<<weightratio<<std::endl;
		std::cout<<"Jacobi"<<std::endl;
		std::cout<<jacobi<<std::endl;
		std::cout<<"f_wstar"<<std::endl;
		std::cout<<f_wstar<<std::endl;*/

}


double MCMCEnv::calcLogDensity() const {
	if(!UpToDate) {
		long taonum = taoCount(Tao);
		LogDensityValue = 0.0;
		//LogDensityValueNonFlag = 0.0;
		//LogDensityValueFlag = 0.0;
		/*std::cout<<"taonum"<<std::endl;
		std::cout<<taonum<<std::endl;*/
		if((taonum > 0) && (taonum < NumOfGenes)){
			//cerr << "datacov" << taonum << endl;
			LogDensityValue += calcLogDensityNonFlagPart(taonum) + calcLogDensityFlagPart(taonum);
			//for(flag=1; flag<=clusternum; flag++){
		} else if(taonum==0){
			LogDensityValue += calcLogDensityNonFlagPart(taonum);
		} else {
			LogDensityValue += calcLogDensityFlagPart(taonum);
		}
		UpToDate = true;
	} 
	return LogDensityValue;

}

long MCMCEnv::taoCount(const TaoSet &newTao) const {
	long taonum = 0;
	for(unsigned long k = 0; k < NumOfGenes; k++){
		taonum += (newTao[k])? 1: 0;
	}
	//cerr << "acceptedp" << endl;
	return taonum;
}

long MCMCEnv::taoCount() const {
	return taoCount(Tao);
}

//void MCMCEnv::changeSigmaToNewTao(const TaoSet &newTao) const {
//
//	//long flag;
//	//long clusternum = Flags.getMaxFlag();
//
//	for(unsigned long k = 0; k < NumOfGenes; k++){
//		sigma1(k, k) = (newTao[k]? constSigma1: 1.0);
//		sigma0(k, k) = ((!newTao[k])? constSigma0: 1.0);
//	}
//	//cerr << "acceptedp" << endl;
//}

double MCMCEnv::calcLogDensityNonFlagPart(long taonum) const {
	if(!UpToDate) {
		LogDensityValueNonFlag = calcLogDensityNonFlagPart(Tao, taonum);
	}
	return LogDensityValueNonFlag;
}

double MCMCEnv::calcLogDensityFlagPart(long taonum) const {
	if(!UpToDate) {
		LogDensityValueFlag = calcLogDensityFlagPart(Tao, taonum);
	}
	return LogDensityValueFlag;
}

double MCMCEnv::calcLogDensityNonFlagPart(const TaoSet &newTao, long taonum) const {
	double gammaprod = 0.0, den = 0.0, K0 = 0.0;
	for(unsigned long j = 1; j <= NumOfGenes - taonum; j++) {
		gammaprod += gammaln(0.5*(NumOfData+deltar+NumOfGenes-taonum-j))
			- gammaln(0.5*(deltar+NumOfGenes-taonum-j));
	}

	K0 = log(beta0*NumOfData+1) * (-0.5*(NumOfGenes-taonum)) + gammaprod;
	//std::cout<<"K0"<<std::endl;
	//std::cout<<K0<<std::endl;

	matrix mid(datacov(0l, false, newTao));
	//std::cout<<"mid"<<std::endl;
	//std::cout<<mid<<std::endl;
	matrix S0(mid.getBlock(1,NumOfGenes - taonum,1,NumOfGenes - taonum));
	S0 += mid.getBlock(NumOfGenes - taonum + 1, 2 * (NumOfGenes - taonum),
		1, NumOfGenes - taonum) * (NumOfData / (beta0 * NumOfData + 1));
	for(unsigned long i = 0; i < NumOfGenes - taonum; i++) {
		S0(i, i) += constSigma0;
	}
	//std::cout<<"S0:"<<std::endl;
	//std::cout<<S0<<std::endl;
	//std::cout<<S0.determinant()<<std::endl;
	try {
		den = K0 + log(fabs(constSigma0)) * (double) (NumOfGenes - taonum) * (0.5*(deltar+NumOfGenes-taonum-1))
			+ log(fabs(S0.determinant())) * (-0.5*(NumOfData+deltar+NumOfGenes-taonum-1));
	} catch(std::logic_error &e) {
		cerr << "Mid" << endl << mid
			<< "Sigma 0:" << endl << sigma0 
			<< "S0" << endl << S0;
		throw e;
	}
	return den;
}

double MCMCEnv::calcLogDensityFlagPart(const TaoSet &newTao, long taonum) const {
	double gammaprod = 0.0, den = 0.0, K1 = 0.0;
	long cclusternum;
	if(taonum) {
		for(TreeContainer::const_iterator itor = MapClusterTrees.begin();
			itor != MapClusterTrees.end(); itor++) {
				cclusternum = itor->second.getSampleNum();
				if(cclusternum){
					for(long j=1; j<=taonum; j++) {
						gammaprod += gammaln(0.5 * (cclusternum + deltar + taonum - j))
							- gammaln(0.5 * (deltar + taonum - j));
					}

					K1 = log(itor->second.calcKcWeight()) 
						+ log(beta1 * cclusternum + 1) * (-(double) taonum / 2.0) + gammaprod;
					//std::cout<<"K1"<<std::endl;
					//std::cout<<K1<<std::endl;
					gammaprod = 0.0;

					matrix mid(datacov(itor->first, true, newTao));
					//cerr << "datacov1" << endl;

					//std::cout<<"S1mid:"<<std::endl;
					//std::cout<<mid<<std::endl;
					matrix S1(mid.getBlock(1, taonum, 1, taonum));
					S1 += mid.getBlock(taonum + 1, 2 * taonum, 1, taonum) 
						* (cclusternum / (beta1 * cclusternum + 1));
					for(unsigned long i = 0; i < taonum; i++) {
						S1(i, i) += constSigma1;
					}
					//std::cout<<"S1:"<<std::endl;
					//std::cout<<S1<<std::endl;
					//std::cout<<S1.determinant()<<std::endl;

					try {
						den += K1 + log(fabs(constSigma1)) * (double) taonum * (0.5 * (deltar+taonum-1))
							+ log(fabs(S1.determinant())) * (-0.5*(cclusternum+deltar+taonum-1));
					} catch(std::logic_error &e) {
						cerr << "Sigma 1:" << endl << sigma1 << "S1" << endl << S1;
						throw e;
					}
				}

		}
	}
	return den;
}

//void MCMCEnv::returnSigmaFromTao() const {
//	for(unsigned long k = 0; k < NumOfGenes; k++){
//		sigma1(k, k) = (Tao[k]? constSigma1: 1.0);
//		sigma0(k, k) = ((!Tao[k])? constSigma0: 1.0);
//	}
//}
//
double MCMCEnv::calcLogDensity(const TaoSet &newTao) const {
	double sigmaConst0, sigmaConst1;
	long taonum = taoCount(newTao);
	//changeSigmaToNewTao(newTao);
	double den = 0.0;
	/*std::cout<<"taonum"<<std::endl;
	std::cout<<taonum<<std::endl;*/
	if((taonum > 0) && (taonum < NumOfGenes)){
		//cerr << "datacov" << taonum << endl;
		den += calcLogDensityNonFlagPart(newTao, taonum) + calcLogDensityFlagPart(newTao, taonum);
		//for(flag=1; flag<=clusternum; flag++){
	} else if(taonum==0){
		den += calcLogDensityNonFlagPart(newTao, taonum);
	} else {
		den += calcLogDensityFlagPart(newTao, taonum);
	}
	//returnSigmaFromTao();
	return den;
}

long MCMCEnv::getNoZeroWeights(unsigned long time) const {
	// get the number of clusters that have non-zero weights at time
	long result = 0;
	//cout << time << "/";
	for(TreeContainer::const_iterator itor = MapClusterTrees.begin(); 
		itor != MapClusterTrees.end(); itor++) {
			//cout << itor->second.weights[time] << ' ';
			if(itor->second.getDeathTime() > time && itor->second.weights[time] > SMALLNUM) {
				result++;
			}
	}
	//cout << "|" << result << endl;
	return result;
}


double MCMCEnv::calcPChos(const TreeMergeSet::value_type &parentChildPair, const bool NULLset,
	const TreeSet &SplitSet) const {
		const ClusterTree &parent = getTreeFromID(parentChildPair.first),
			&child = getTreeFromID(parentChildPair.second);
		unsigned long time = getTreeFromID(parentChildPair.second).getBornTime();
		long G_star = getClusterNumber();
		long G1_star = getNonEmptyClusterNumber();
		long G0_star = G_star - G1_star;
		long n;
		double s = 0.0;
		/*std::cout<<"Three chosen probabilities."<<std::endl;
		std::cout<<1.0/((double)(G_star*G1_star))<<std::endl;
		std::cout<<1.0/(double)G_star<<std::endl;
		std::cout<<2.0/(double)G_star<<std::endl;*/
		if(NULLset){
			return 1.0 / ((double)(G_star * G1_star));
		}else{
			n = calcSimilarity(time, parent, child, NULLset, SplitSet);

			s = ((n == 2)? (2.0 / (double) G_star): (1.0 / (double) G_star));

			return s;
			//return s/(double)time;
			//return pow(s,time-1);
			//return pow(s,time);
		}
}

// Return an integar indicating the relationship of the two splitted clusters
// 0: a cluster is empty.
// 1: The relationship is not closest.
// 2: The relationship is closest.
// NULLset tells clusterson is empty or not.
double MCMCEnv::calcSimilarity(unsigned long time, const ClusterTree &parent, const ClusterTree &child,
	const bool NULLset, const TreeSet &SplitSet) const {
		long n = SplitSet.size();
		long i, m = 2;
		double dist_ms, dist_mi;
		if(NULLset){
			return 0;
		}else{
			dist_ms = this->calcClusterDist(time, parent, time, child);
			for(TreeSet::const_iterator itor = SplitSet.begin(); itor != SplitSet.end(); itor++){
				dist_mi = this->calcClusterDist(time, parent, itor->first, getTreeFromID(itor->second));
				if(dist_ms > dist_mi){
					m = 1;
					break;
				}
			}
		}
		return m;
}

double MCMCEnv::calcClusterDist(const unsigned long time1, const ClusterTree &cluster1, 
	const long time2, const ClusterTree &cluster2) const {
		double c = 0.0;
		matrix a,b;
		long i;
		a = calcClusterMean(time1, cluster1);
		b = calcClusterMean(time2, cluster2);
		/*std::cout<<a<<std::endl;
		std::cout<<b<<std::endl;*/
		for(i = 0; i < NumOfGenes; i++){
			c += pow(a(0,i)-b(0,i),2);
		}
		return c;

}

matrix MCMCEnv::calcClusterMean(const unsigned long time, const ClusterTree &cluster) const {

	matrix gdata(1, NumOfGenes);
	long i, j, s = 0, n = 0;
	for(i = time; i < NumOfTP; i++){
		for(j = 0; j < Flags.cols(i); j++){
			if(Flags(i, j) == cluster.ID){
				gdata += Data.getRow(CoorToIndex(TimeCoordinate(i, j)));
				n++;
			}
		}
	}

	gdata *= (1/(double)n);

	return gdata;
}

// Accept probability for split move, the input cluster index must be nonempty at least in the beginning.
// newz is new cluster indicator after split move
// newweight is new weight matrix after split move
// This function follows a random sample from 'splitset' result.
double MCMCEnv::apSplit(const MCMCEnv::TreeSet::value_type &split_pair, const bool NULLset, 
	const MCMCEnv::TreeSet &splitsets, const long splitnum, const double P_alloc, const double f_ui,
	const MCMCEnv &oldenv, unsigned long long newBornTreeID) const {
		double densityratio, modelratio;
		double weightratio = 1.0;
		double P_chos = 1.0;
		double AP = 1.0;
		double dbratio = 1.0;
		double jacobi = 1.0;
		long i;
		if(splitnum == 0){
			return 0.0;
		}else{
			long taonum = taoCount(Tao);
			densityratio = exp(calcLogDensityFlagPart(Tao, taonum) - oldenv.calcLogDensityFlagPart(Tao, taonum));
			const ClusterTree &oldTree = oldenv.getTreeFromID(split_pair.second),
				&newTree = getTreeFromID(split_pair.second),
				&newChildTree = buildForest? getTreeFromID(newBornTreeID)
					: getTreeFromID(newTree.getChildID(split_pair.first));
			unsigned long time = split_pair.first;
			for(i = split_pair.first; i < NumOfTP; i++){
				// If the first cluster is an empty cluster, directly return 0.0.
				if(oldTree.weights[i] > SMALLNUM){
					weightratio *= pow(newTree.weights[i], alpha - 1) * pow(newChildTree.weights[i], alpha - 1)
						/ pow(oldTree.weights[i], alpha - 1) 
						/ betaf(alpha, ((double)oldenv.getNoZeroWeights(i))*alpha);
					jacobi *= oldTree.weights[i];
				}
			}
			// The ratio of f(Tree_new) and f(Tree_old).
			// If we set the same unchange and change probability, then it is equal to 1.0; 
			modelratio = 1.0;
			// The ratio of d_new and b_old, here dbratio = 1.0, because d_new = b_old = 0.5.
			dbratio = (1.0-bk)/bk;
			// Parameter 'splitnum': The number of cluster vector that can be used to do split move.
			//splitnum = 1.0;//
			// A function must be written to summarize the number of individuals in new and old clusters
			//P_alloc = 1.0;
			// Probabilities must be designed to describe the merge clusters.
			P_chos = calcPChos(TreeMergeSet::value_type(newTree.ID, newChildTree.ID), NULLset, splitsets);


			AP = densityratio * weightratio * modelratio * dbratio * P_chos 
				* ((double)splitnum) / f_ui / P_alloc * jacobi;

			//std::cout<<"Split"<<std::endl;
			//std::cout<<"densityratio"<<std::endl;
			//std::cout<<densityratio<<std::endl;
			//std::cout<<"weightratio"<<std::endl;
			//std::cout<<weightratio<<std::endl;
			//std::cout<<"P_chos"<<std::endl;
			//std::cout<<P_chos<<std::endl;
			//std::cout<<"splitnum"<<std::endl;
			//std::cout<<splitnum<<std::endl;
			//std::cout<<"f_ui"<<std::endl;
			//std::cout<<f_ui<<std::endl;
			//std::cout<<"P_alloc"<<std::endl;
			//std::cout<<P_alloc<<std::endl;
			//std::cout<<"Jacobi"<<std::endl;
			//std::cout<<jacobi<<std::endl;
			////std::cout<<"Split accepted probability: "<<std::endl;
			//cerr << "AP: " << AP << endl;
			return f_fmin(AP,1.0);
		}
}

double MCMCEnv::apMerge(const TreeMergeSet::value_type &mergePair, const TreeSet &splitSetsAfterMerge,
	const long mergesetsrows, const long splitnumaftermerge, const double P_alloc, const double f_ui,
	const MCMCEnv &oldenv) const {

		double densityratio, modelratio, bdratio; 
		double weightratio = 1.0;
		double jacobi = 1.0;
		double P_chos;
		double AP;
		bool NULLset = false;
		//cout << "Mergepair: first = " << mergePair.first << "; second = "
		//	<< mergePair.second << endl;
		const ClusterTree &oldparent = oldenv.getTreeFromID(mergePair.first), 
			&oldchild = oldenv.getTreeFromID(mergePair.second),
			&newparent = getTreeFromID(mergePair.first);
		unsigned long time = oldchild.getBornTime();

		if((!mergesetsrows) || (!splitnumaftermerge)){
			return 0.0;
		}else{
			long taonum = taoCount(Tao);
			densityratio = exp(calcLogDensityFlagPart(Tao, taonum) - oldenv.calcLogDensityFlagPart(Tao, taonum));
			//cerr << "New density: " << calcLogDensity() << endl;
			//cerr << "Old density: " << oldenv.calcLogDensity() << endl;
			/*std::cout<<"newz"<<std::endl;
			std::cout<<newz<<std::endl;
			std::cout<<"newweight"<<std::endl;
			std::cout<<newweight<<std::endl;
			std::cout<<"NEW DENSITY"<<std::endl;
			std::cout<<zw_density(newz,newweight,mu1,mu0,sigma0,sigma1,beta0,beta1,deltar)<<std::endl;
			std::cout<<"OLD DENSITY"<<std::endl;
			std::cout<<density(mu1,mu0,sigma0,sigma1,beta0,beta1,deltar)<<std::endl;*/
            
			for(unsigned long i = time; i < NumOfTP; i++){
				if(oldparent.weights[i] > SMALLNUM && oldchild.weights[i] > SMALLNUM ){
					weightratio *= pow(newparent.weights[i], alpha - 1) *
					    betaf(alpha, ((double) getNoZeroWeights(time)) * alpha) 
						/ pow(oldparent.weights[i], alpha - 1) 
						/ pow(oldchild.weights[i], alpha - 1);
					jacobi /= newparent.weights[i];
			
				}

			}
			modelratio = 1.0;
			bdratio = bk / (1.0 - bk);
			try {
				P_chos = oldenv.calcPChos(mergePair, NULLset, splitSetsAfterMerge);
			} catch(std::logic_error &e) {
				cout << "Logic error catched calculating P_chos." << endl;
				throw e;
			}

		}


		AP = densityratio * weightratio * modelratio * bdratio * P_alloc
			/ ((double)splitnumaftermerge) * f_ui / P_chos * jacobi;

		//std::cout<<"Merge"<<std::endl;
		//std::cout<<"densityratio"<<std::endl;
		//std::cout<<densityratio<<std::endl;
		//std::cout<<"weightratio"<<std::endl;
		//std::cout<<weightratio<<std::endl;
		//std::cout<<"P_chos"<<std::endl;
		//std::cout<<P_chos<<std::endl;
		//std::cout<<"splitnumaftermerge"<<std::endl;
		//std::cout<<splitnumaftermerge<<std::endl;
		//std::cout<<"f_ui"<<std::endl;
		//std::cout<<f_ui<<std::endl;
		//std::cout<<"P_alloc"<<std::endl;
		//std::cout<<P_alloc<<std::endl;
		//std::cout<<"Jacobi"<<std::endl;
		//std::cout<<jacobi<<std::endl;
		//std::cout<<"Merge accepted probability: "<<std::endl;
		//cerr << AP << endl;

		return f_fmin(AP,1.0);
}

double MCMCEnv::apBirth(const TreeSet::value_type &splitPair, const long splitnumafterbirth,
	const double weightratio, const double jacobi, const double f_wstar,
	const MCMCEnv &oldenv) const {
		double densityratio, modelratio;
		double AP = 1.0;
		double dbratio = 1.0;
		long i;
		if(splitnumafterbirth == 0){
			return 0.0;
		}else{
			long taonum = taoCount(Tao);
			densityratio = exp(calcLogDensityFlagPart(Tao, taonum) - oldenv.calcLogDensityFlagPart(Tao, taonum));

			// The ratio of f(Tree_new) and f(Tree_old).
			// If we set the same unchange and change probability, then it is equal to 1.0; 
			modelratio = 1.0;
			// The ratio of d_new and b_old, here dbratio = 1.0, because d_new = b_old = 0.5.
			dbratio = (1.0 - bk) / bk;

			AP = densityratio * weightratio * modelratio * dbratio * jacobi
				/ (double) splitnumafterbirth / f_wstar;

			/*std::cout<<"densityratio"<<std::endl;
			std::cout<<densityratio<<std::endl;
			std::cout<<"weightratio"<<std::endl;
			std::cout<<weightratio<<std::endl;
			std::cout<<"splitnumafterbirth"<<std::endl;
			std::cout<<splitnumafterbirth<<std::endl;
			std::cout<<"f_wstar"<<std::endl;
			std::cout<<f_wstar<<std::endl;
			std::cout<<"Jacobi"<<std::endl;
			std::cout<<jacobi<<std::endl;*/
			return f_fmin(AP,1.0);
		}
}

double MCMCEnv::apDeath(const TreeMergeSet::value_type &mergePair, const long emptynumbeforedeath,
	const double weightratio, const double jacobi, const double f_wstar, 
	const MCMCEnv &oldenv) const {
		double densityratio,modelratio,bdratio;
		double AP;

		if(emptynumbeforedeath == 0){
			return 0.0;
		}else{
			long taonum = taoCount(Tao);
			densityratio = exp(calcLogDensityFlagPart(Tao, taonum) - oldenv.calcLogDensityFlagPart(Tao, taonum));
			modelratio = 1.0;
			bdratio = bk / (1.0 - bk);
		}

		AP = densityratio * weightratio * modelratio * bdratio 
			* f_wstar * emptynumbeforedeath * jacobi;

		return f_fmin(AP,1.0);
}

double MCMCEnv::apTailBirth(const TreeSet::value_type &splitPair, const long splitnumafterbirth,
	const double weightratio, const double jacobi, const double f_wstar,
	const MCMCEnv &oldenv) const {
		double densityratio, modelratio;
		double AP = 1.0;
		double dbratio = 1.0;
		long i;
		if(splitnumafterbirth == 0){
			return 0.0;
		}else{
			long taonum = taoCount(Tao);
			densityratio = exp(calcLogDensityFlagPart(Tao, taonum) - oldenv.calcLogDensityFlagPart(Tao, taonum));

			// The ratio of f(Tree_new) and f(Tree_old).
			// If we set the same unchange and change probability, then it is equal to 1.0; 
			modelratio = 1.0;
			// The ratio of d_new and b_old, here dbratio = 1.0, because d_new = b_old = 0.5.
			dbratio = (1.0 - bk) / bk;

			AP = densityratio * weightratio * modelratio * dbratio * jacobi
				/ (double) splitnumafterbirth / f_wstar;

			/*std::cout<<"densityratio"<<std::endl;
			std::cout<<densityratio<<std::endl;
			std::cout<<"weightratio"<<std::endl;
			std::cout<<weightratio<<std::endl;
			std::cout<<"splitnumafterbirth"<<std::endl;
			std::cout<<splitnumafterbirth<<std::endl;
			std::cout<<"f_wstar"<<std::endl;
			std::cout<<f_wstar<<std::endl;
			std::cout<<"Jacobi"<<std::endl;
			std::cout<<jacobi<<std::endl;*/
			return f_fmin(AP,1.0);
		}
}

double MCMCEnv::apTailDeath(const TreeSet::value_type &mergePair, const long emptynumbeforedeath,
	const double weightratio, const double jacobi, const double f_wstar, 
	const MCMCEnv &oldenv) const {
		double densityratio,modelratio,bdratio;
		double AP;

		if(emptynumbeforedeath == 0){
			return 0.0;
		}else{
			long taonum = taoCount(Tao);
			densityratio = exp(calcLogDensityFlagPart(Tao, taonum) - oldenv.calcLogDensityFlagPart(Tao, taonum));
			modelratio = 1.0;
			bdratio = bk / (1.0 - bk);
		}

		AP = densityratio * weightratio * modelratio * bdratio 
			* f_wstar * emptynumbeforedeath * jacobi;

		return f_fmin(AP,1.0);
}

unsigned long MCMCEnv::getClusterNumber() const {
	// get the number of trees
	return MapClusterTrees.size();
}

unsigned long MCMCEnv::getNonEmptyClusterNumber() const {
	// get the number of trees that are not empty throughout the time
	unsigned long count = 0;
	for(TreeContainer::const_iterator itor = MapClusterTrees.begin();
		itor != MapClusterTrees.end(); itor++) {
			if(itor->second.getSampleNum()) {
				count++;
			}
	}
	return count;
}

unsigned long MCMCEnv::getClusterNumber(unsigned long time) const {
	// get the number of trees at time t
	unsigned long count = 0;
	for(TreeContainer::const_iterator itor = MapClusterTrees.begin();
		itor != MapClusterTrees.end(); itor++) {
			if(itor->second.getBornTime() <= time) {
				count++;
			}
	}
	return count;
}

unsigned long MCMCEnv::getNonEmptyClusterNumber(unsigned long time) const {
	// get the number of trees that are not empty throughout the time
	unsigned long count = 0;
	for(TreeContainer::const_iterator itor = MapClusterTrees.begin();
		itor != MapClusterTrees.end(); itor++) {
			if(itor->second.getSampleNum(time)) {
				count++;
			}
	}
	return count;
}

double MCMCEnv::getWeightFromSample(unsigned long time, unsigned long sampleNum) const {
	return getTreeFromID(Flags(time, sampleNum)).weights[time];
}

void MCMCEnv::sampleWeight() {

	double wnum;
	std::vector<long> midv;
	for(unsigned long t = 0; t < NumOfTP; t++) {
		
		wnum = 0.0;
		for(TreeContainer::iterator itor = MapClusterTrees.begin();
			itor != MapClusterTrees.end(); itor++) {
				
				if(itor->second.getBornTime() <= t && itor->second.getDeathTime() > t) {
					itor->second.weights[t] = ran_gamma(alpha 
						+ itor->second.getSampleNum(t), shape);
					//cerr << '[' << itor->second.ID << ',' << t << ']' 
					//	<< itor->second.getSampleNum(t)
					//	<< ' ' << itor->second.weights[t] << endl; 
					wnum += itor->second.weights[t];
					//fnum = z.tflagnum(flag,t+1);
					//weight.setValue(t,flag-1,);
					//wnum += weight.getValue(t,flag-1);
				}
		}

		// normalize all weights
		for(TreeContainer::iterator itor = MapClusterTrees.begin();
			itor != MapClusterTrees.end(); itor++) {
				if(itor->second.getBornTime() <= t && itor->second.getDeathTime() > t) {
					if(wnum > SMALLNUM) {
						itor->second.weights[t] = itor->second.weights[t] / wnum;
					} else {
						itor->second.weights[t] = 1.0;
					}
				}
		}
		//midv.clear();
	}
	UpToDate = false;

}

void MCMCEnv::sampleTao() {

	long coin2, coin3;

	TaoSet newTao = Tao;

	if(ran_ber(0.5) || newTao.size() <= 1) {
		coin2 = ran_iunif(0, newTao.size() - 1);
		newTao[coin2] = !newTao[coin2]; //1. Change value.
	} else {
		coin2 = ran_iunif(0, newTao.size() - 1);
		coin3 = ran_iunif(0, newTao.size() - 2);
		if(coin3 >= coin2){
			coin3++;
		}
		bool tao_mid = newTao[coin2];
		newTao[coin2] = newTao[coin3];
		newTao[coin3] = tao_mid;          //or 2. Swap;
	}

	double acceptedp = f_fmin(1.0, exp(calcLogDensity(newTao) - calcLogDensity()));
	//std::cout << "current tao: " << taoCount() << "; new tao: " << taoCount(newTao) << "; tao p: " << acceptedp << std::endl;

	if(ran_ber(acceptedp)){
		Tao = newTao;
		//returnSigmaFromTao();
		UpToDate = false;
	} 

	/*std::cout<<"Tao"<<std::endl;
	for(i=0; i<tao_length; i++)
	std::cout<<tao.at(i)<<' ';
	std::cout<<std::endl;*/

}

void MCMCEnv::sampleZ() {
	long k,l,midnum;
	std::vector<unsigned long long> flagl;
	std::vector<double> prob_z;
	//ClusterFlags z_mid = z;
	// reset the first timepoint to root_ID
	//for(unsigned long j = 0; j < Flags.cols(0); j++) {
	//	Flags(0, j) = root_ID;
	//}

	long taonum = taoCount(Tao);
	//double denNonFlag = calcLogDensityNonFlagPart(taonum);

	if(taonum) {
		for(unsigned long t = 0; t < Flags.rows(); t++){
		
			flagl.clear();
			for(TreeContainer::iterator itor = MapClusterTrees.begin();
				itor != MapClusterTrees.end(); itor++) {
					if(itor->second.weights[t] > SMALLNUM) {
						flagl.push_back(itor->second.ID);
					}
			}

			for(unsigned long j = 0; j < Flags.cols(t); j++){
				(getTreeFromID(Flags(t, j)).samples[t])--;
				if(flagl.size() > 1) {
					prob_z.clear();
					for(std::vector<unsigned long long>::const_iterator itor = flagl.begin();
						itor != flagl.end(); itor++) {
							Flags(t, j) = *itor;
							(getTreeFromID(*itor).samples[t])++;
							UpToDate = false;
							prob_z.push_back(calcLogDensityFlagPart(taonum));
							(getTreeFromID(*itor).samples[t])--;
					}
					Flags(t, j) = ran_num_log(prob_z, flagl);
				} else {
					Flags(t, j) = flagl[0];
				}
				(getTreeFromID(Flags(t, j)).samples[t])++;


				//for(k=0; k<flagl.size(); k++){
				//	z(i,j) = flagl.at(k);
				//	prob_z.push_back(z_density(z, mu1, mu0, sigma0, sigma1, beta0, beta1, deltar));
				//	// optimize?
				//		
				//}
				//z(i,j) = ran_num(prob_z,flagl);
				//z_mid(i,j) = z(i,j);
				/*std::cout<<"Data"<<"("<<i<<","<<j<<")"<<std::endl;
					for(l=0; l<flagl.size(); l++)
				std::cout<<prob_z.at(l)<<' ';
				std::cout<<std::endl;*/
				/*if(i==3&&j==0){
					for(l=0; l<flagl.size(); l++)
						std::cout<<prob_z.at(l)<<' ';
					std::cout<<std::endl;
				}*/
				
			}	 
		}
		UpToDate = false;
	}
	// Do we need to reflag the cluster indicators?
	//z.reflag();
	/*std::cout<<"Sampled z"<<std::endl;
	std::cout<<z<<std::endl;*/
}

double MCMCEnv::treeSummary() const {
	// Tree Coding method
	double value = 0;
	for(TreeContainer::const_iterator itor = MapClusterTrees.begin();
		itor != MapClusterTrees.end(); itor++) {
			if(itor->second.getSampleNum()) {
				value += pow(4, (double) (NumOfTP - 1 - itor->second.getBornTime()));
			}
	}
	return value;
}

ostream &MCMCEnv::writeTree(ostream &os) const {
	os << "[] Tree was born;" << endl 
		<< "() A new Tree with this ID is born from the current tree." << endl;
	deque<unsigned long long> IdToShow;
	for(RootIDContainer::const_iterator itor = rootIDs.begin();
		itor != rootIDs.end(); itor++) {
			IdToShow.push_back(itor->first);
			try {
				while(!IdToShow.empty()) {
					const ClusterTree &currTree = getTreeFromID(IdToShow.front());
					IdToShow.pop_front();
					for(unsigned long t = 0; t < NumOfTP; t++) {
						if(t < currTree.getBornTime()) {
							os << "   " << setfill(' ') << setw(TREE_ID_SIZE + 2) << ' ';
						} else if(t == currTree.getBornTime()) {
							os << "   [" << setfill(' ') << setw(TREE_ID_SIZE) << currTree.ID << "]";
						} else {
							if(currTree.getChildID(t)) {
								os << " - (" << setfill(' ') << setw(TREE_ID_SIZE) << currTree.getChildID(t) << ")";
								IdToShow.push_back(currTree.getChildID(t));
							} else {
								os << " - " << setfill('-') << setw(TREE_ID_SIZE + 2) << '-';
							}
						}
					}
					os << setfill(' ') << endl;
				}
			} catch(std::logic_error &e) {
				cout << "Logic error is caught in writeTree." << endl;
				throw e;
			}
	}
	return os;
}

ostream &MCMCEnv::writeTao(ostream &os) const {
	for(unsigned long i = 0; i < NumOfGenes; i++) {
		os << (Tao[i]? '1': '0') << ' ';
	}
	return os << std::endl;
}

ostream &MCMCEnv::writeFlags(ostream &os) const {
	return os << Flags;
}

ostream &MCMCEnv::writeWeights(ostream &os) const {
	for(MCMCEnv::TreeContainer::const_iterator itor = MapClusterTrees.begin();
		itor != MapClusterTrees.end(); itor++) {
			itor->second.writeWeights(os);
	}
	return os;
}

ostream &MCMCEnv::writeStatus(ostream &os) const {
	os<<"Cluster proportions in each cluster: "<<std::endl;
	writeWeights(os);
	os<<std::endl;
	os<<"Featured genes: "<<std::endl;
	writeTao(os);
	os << endl;
	os<<"Cluster index: "<<std::endl;
	writeFlags(os);
	os<<std::endl;
	os<<"Tree: "<<std::endl;
	writeTree(os);
	os<<std::endl;
	return os;
}

void MCMCEnv::testNumberSet() {

	ClusterTree &child = createTree(rootIDs[0], 1), &parent = getTreeFromID(rootIDs[0]);

	cout << "Test number set." << endl;

	parent.samples[0] = 3;
	parent.samples[1] = 6;
	parent.samples[2] = 9;
	parent.samples[3] = 9;

	child.samples[0] = 0;
	child.samples[1] = 0;
	child.samples[2] = 0;
	child.samples[3] = 3;

	parent.weights[0] = 1;
	parent.weights[1] = 1;
	parent.weights[2] = 1;
	parent.weights[3] = 0.7171;

	child.weights[0] = 0;
	child.weights[1] = 0;
	child.weights[2] = 0;
	child.weights[3] = 0.2829;

	Flags.setFlag(0,0,1);
	Flags.setFlag(0,1,1);
	Flags.setFlag(0,2,1);
	


	Flags.setFlag(1,0,1);
	Flags.setFlag(1,1,1);
	Flags.setFlag(1,2,1);
	Flags.setFlag(1,3,1);
	Flags.setFlag(1,4,1);
	Flags.setFlag(1,5,1);
	



	Flags.setFlag(2,0,1);
	Flags.setFlag(2,1,1);
	Flags.setFlag(2,2,1);
	Flags.setFlag(2,3,1);
	Flags.setFlag(2,4,1);
	Flags.setFlag(2,5,2);
	Flags.setFlag(2,6,1);
	Flags.setFlag(2,7,1);
	Flags.setFlag(2,8,1);
	


	Flags.setFlag(3,0,1);
	Flags.setFlag(3,1,1);
	Flags.setFlag(3,2,1);
	Flags.setFlag(3,3,1);
	Flags.setFlag(3,4,1);
	Flags.setFlag(3,5,1);
	Flags.setFlag(3,6,1);
	Flags.setFlag(3,7,1);
	Flags.setFlag(3,8,1);
	Flags.setFlag(3,9,2);
	Flags.setFlag(3,10,2);
	Flags.setFlag(3,11,2);
}

ostream &operator<<(ostream &os, const MCMCEnv &env) {
	long i_t = MCMCEnv::NumOfTP;
	long i_g = MCMCEnv::NumOfGenes;
	int i;
	os<<"Time length:"<<std::endl;
	os<<i_t;
	os<<std::endl;
	os<<"Gene number:"<<std::endl;
	os<<i_g;
	os<<std::endl;
	os<<"Data number:"<<std::endl;
	os<<MCMCEnv::getDataNum();
	os<<std::endl;
	os<<"Data number in each cluster:"<<std::endl;
	for(i=0; i<i_t; i++)
		os<<MCMCEnv::DataNumEachTP[i]<<' ';
	os<<std::endl;
	os<<"Cluster number in each cluster:"<<std::endl;
	for(i=0; i<i_t; i++)
		os<<env.getClusterNumber(i)<<' ';
	os<<std::endl;
	os<<"The original data: "<<std::endl;
	os << env.Data;
	os<<std::endl;

	return os;
}

double & MCMCEnv::weight(unsigned long long ID, unsigned long time) {
	return getTreeFromID(ID).weights[time];
}

std::ostream &MCMCEnv::writeSet(std::ostream &os, const TreeSet &set) const {
	os << "Split set:" << endl;
	for(TreeSet::const_iterator itor = set.begin(); itor != set.end(); itor++) {
		os << "Time: " << itor->first << " Split from: " << itor->second << endl;
	}
	return os;
}

std::ostream &MCMCEnv::writeSet(std::ostream &os, const TreeMergeSet &set) const {
	os << "Merge set:" << endl;
	for(TreeMergeSet::const_iterator itor = set.begin(); itor != set.end(); itor++) {
		os << "Parent: " << itor->first << " Child: " << itor->second << endl;
	}
	return os;
}


#undef SMALLNUM
#undef TREE_ID_SIZE