#include "ClusterFlags.h"
#include "Matrix.h"
#include <stdexcept>
using namespace std;
using namespace std::tr1;

// Initialize the clusterflags
// n_timepoints: number of time points in this study
// samples: a vector indicating the number of samples at each time point
//			will throw exceptions if samples does not contain enough items
ClusterFlags::ClusterFlags()
{

}
ClusterFlags::ClusterFlags(long n_timepoints, const std::vector<unsigned long> &samples)
	: data(vector<vector<long> >(n_timepoints))
{
	for(long i=0; i < n_timepoints; i++) {
		data[i].resize(samples.at(i));
	}
}


ClusterFlags::ClusterFlags(const std::vector<std::vector<long> > &one)
{
  data = one;
}

ClusterFlags::~ClusterFlags(void)
{
}


long ClusterFlags::operator()(long row, long col) const {
	return getFlag(row, col);
}

long & ClusterFlags::operator()(long row, long col) {
	return data.at(row).at(col);
}

// Reload = 
ClusterFlags &ClusterFlags::operator=(const ClusterFlags &clu)
{
 data = clu.data;
	return *this;
}


long ClusterFlags::rows() const {
	return data.size();
}

long ClusterFlags::cols(int i) const {
	if(data.size() == 0) {
		return 0;
	}
	return data.at(i).size();
}


// Fill all the flags with the value
// Note that cluster index begins with 0
// Therefore the default value of -1 means flag not set
void ClusterFlags::fillFlags(long flag) {
	for(vector<vector<long> >::iterator itor = data.begin(); itor != data.end(); itor++) {
		for(vector<long>::iterator itor_time = itor->begin(); itor_time != itor->end(); itor_time++) {
			*itor_time = flag;
		}
	}
}

// Get the flag for a specific time/sample
// notice that both time and sample is 0-based
// will throw out_of_bounds if so
long ClusterFlags::getFlag(long time, long sample) const {
	return data.at(time).at(sample);
}

// Set the flag for a specific time/sample
// will throw out_of_bounds if so
void ClusterFlags::setFlag(long time, long sample, long flag) {
	data.at(time).at(sample) = flag;
}


// Get the earliest time point for a specific flag (Real time point).
// If not found, return -1
long ClusterFlags::getEarliestTime(long flag) const {
	for(long i = 0; i < data.size(); i++) {
		for(long j = 0; j < data.at(i).size(); j++) {
			if(data.at(i).at(j) == flag) {
				return i+1;
			}
		}
	}
	return -1;
}

// Get the number of cluster flags in the set
// note that this is 1 larger than the largest number
long ClusterFlags::getMaxFlag() const {
	long max = -1;
	for(vector<vector<long> >::const_iterator itor = data.begin(); itor != data.end(); itor++) {
		for(vector<long>::const_iterator itor_time = itor->begin(); itor_time != itor->end(); itor_time++) {
			if(*itor_time > max) {
				max = *itor_time;
			}
		}
	}
	return max;
}

std::vector<std::vector<long> > ClusterFlags::getTimeSampleByCluster(long flag) const {
	vector<vector<long> > result;
	for(vector<vector<long> >::const_iterator itor = data.begin(); itor != data.end(); itor++) {
		vector<long> current_vec;
		for(long i = 0; i < itor->size(); i++) {
			if(itor->at(i) == flag) {
				current_vec.push_back(i);
			}
		}
		result.push_back(current_vec);
	}
	return result;
}

// If this 'flag' is not exist in clusterflags, a 0 matrix will return.
vector<vector<long> > ClusterFlags::onezerocluster(long flag) const
{
 vector<vector<long> > result;
	for(vector<vector<long> >::const_iterator itor = data.begin(); itor != data.end(); itor++) {
		vector<long> current_vec;
		for(long i = 0; i < itor->size(); i++) {
			if(itor->at(i) == flag) {
				current_vec.push_back(1);
			}
			else {
				current_vec.push_back(0);
			}
		}
		result.push_back(current_vec);
	}
	return result;
}



ClusterFlags ClusterFlags::sortz(const ClusterFlags &one, long n_timepoints, vector<long> NT){
	ClusterFlags two;
	return two;

}


void ClusterFlags::reflag()
{
   vector<long> x;
   vector<long> one;
   vector<long> two;
   long i,j, l=0, k=1;
   long length = 0;
   for(i=0; i<data.size(); i++)
	   for(j=0; j<data.at(i).size(); j++)
		   one.push_back(data.at(i).at(j));

   x.push_back(one.at(0));
   two.push_back(k);
 for(i=1; i<one.size();i++)
 {
	 if(one.at(i)==one.at(i-1)){
	   
		 two.push_back(two.at(i-1));  
	 
	 } else {
	    
		 while(l < k)
		   {
		    if(one.at(i)==x.at(l))
				break;
			l++;
		   }
		 if(l < k)
			 two.push_back(l+1);
		 else {
			 k++;
		     x.push_back(one.at(i));
			 two.push_back(k);
		     
		 }
		 l = 0;
		  
	 }
 }

   for(i=0; i<data.size(); i++){
	  
	   for(j=0; j<data.at(i).size(); j++)
		   data.at(i).at(j) = two.at(length + j);
       length += data.at(i).size();

   }
}



//void reflag(std::vector<long> &one)
//{
// std::vector<long> x;
// std::vector<long> two;
// int i,j=0,k=1;
// x.push_back(one.at(0));
// two.push_back(k);
// for(i=1; i<one.size();i++)
// {
//	 if(one.at(i)==one.at(i-1)){
//	   
//		 two.push_back(two.at(i-1));  
//	 
//	 } else {
//	    
//		 while(j < k)
//		   {
//		    if(one.at(i)==x.at(j))
//				break;
//			j++;
//		   }
//		 if(j < k)
//			 two.push_back(j+1);
//		 else {
//			 k++;
//		     x.push_back(one.at(i));
//			 two.push_back(k);
//		     
//		 }
//		 j = 0;
//		  
//	 }
// }
// one = two;
//}



// Summarize the data number for some cluster 'flag'.
// With a large flag, k should just be equal to 0, it is not need to output an error information.
long ClusterFlags::flagnum(long flag) const
{
	long k = 0;
	long i,j;
	if((flag <= (this->getMaxFlag()))&&(flag>0))
	{
		for(i=0; i< data.size(); i++){
			for(j=0; j<data.at(i).size(); j++){
			   if(data.at(i).at(j)==flag)
			   	   k++;
				
			}


		}

	} 
	//else {
	//
	////cout<<"Flag is out of range!!!"<<endl;
	//throw(exception("Flag is out of range!!!"));
	//}
	return k;
}

// Return the number of cluster indicator for some clusternumber(flag) in time point ntime.
// ntime is 0-based time point!
long ClusterFlags::tflagnum(long flag, long ntime) const{

 long k=0;
 long i;
 if((flag <= (this->getMaxFlag()))&&(flag>0)&&(ntime<data.size()))
	{
		for(i=0; i< data.at(ntime).size(); i++){
			 if(data.at(ntime).at(i)==flag)
				   k++;
			}

	} 
    /*else {
	
	throw(exception("Flag is out of range!!!"));
	}*/
    //cout<<data.size()<<endl;
	return k;

}


// Useless!!!!!!
// Get the clusternumber in time stage 'ntime'.
// Will return the maximum cluster number even if some small cluster number missed.
std::vector<long> ClusterFlags::getTflag(long ntime)const{
	//long clusternum = this->getMaxFlag();
	std::vector<long> one(ntime);
	one.at(0) = 1;
	long i,j;
	long mid = 1;
	if(ntime<data.size()){
	for(i=1; i<ntime; i++){
		mid = data.at(i).at(0);
		for(j=0; j<(data.at(i).size()-1); j++){
		   if(data.at(i).at(j+1)>data.at(i).at(j))
			  mid = data.at(i).at(j+1);
		}

		one.at(i) =  mid;
		if(one.at(i)<one.at(i-1)){
		  one.at(i) = one.at(i-1);
		}
	}
	} else{
	  std::cout<<"ntime is out of range."<<std::endl;
	}
	return one;

}



// Get different cluster indicators in time "ntime"
std::vector<long> ClusterFlags::getTindicator(long ntime)const{
	long n, i, j;
	std::vector<long> one;
	std::vector<long> two;
	bool init = false;
	if(ntime<data.size()){
	  n = data.at(ntime).size();
	  one.push_back(data.at(ntime).at(0));
	  for(i=1; i<n; i++){
		  for(j=0;j<one.size(); j++){
			  if(data.at(ntime).at(i)==one.at(j)){
			    init = true;
				break;
			  }
		  }
		  if(!init){
		    one.push_back(data.at(ntime).at(i));
		  }
		  init = false;
	  }
	}else
	{
	 std::cout<<"ntime is out of range."<<std::endl;
	}

	return one;

}


// Justify whether flagx appears in the set of ntime cluster indicators.
bool ClusterFlags::flaginT(long flagx, long ntime) const{
	long i;
	bool init = false;
	long n = data.at(ntime).size();
	for(i=0; i<n; i++){
		if(flagx == data.at(ntime).at(i)){
		  init = true;
		  break;
		}

	}
	return init;

}








// Get different cluster indicators from time 1 to time "ntime"
std::vector<long> ClusterFlags::getTtoTindicator(long ntime)const{
    long n, i, j, k;
	std::vector<long> one;
	bool init = false;
	if(ntime<data.size()){
      one.push_back(data.at(0).at(0));
	  for(k=0; k<=ntime; k++){
	     n = data.at(k).size();
	  
	     for(i=1; i<n; i++){
		     for(j=0;j<one.size(); j++){
			     if(data.at(k).at(i)==one.at(j)){
			       init = true;
				    break;
			     }
		     }
		     if(!init){
		       one.push_back(data.at(k).at(i));
		     }
		     init = false;
	     }
	  }
	}else
	{
	 std::cout<<"ntime is out of range."<<std::endl;
	}

	return one;
}


// Return nonempty cluster number
long ClusterFlags::nonemptynum()const{
	long n, i, j, k, s=1;
	long m = data.size();
	std::vector<long> one;
	bool init = false;
	
      one.push_back(data.at(0).at(0));
	  for(k=0; k<m; k++){
	     n = data.at(k).size();
	  
	     for(i=0; i<n; i++){
		     for(j=0;j<one.size(); j++){
			     if(data.at(k).at(i)==one.at(j)){
			       init = true;
				    break;
			     }
		     }
		     if(!init){
		       one.push_back(data.at(k).at(i));
			   s++;
		     }
		     init = false;
	     }
	  }
	

	return s;
	
}



std::ostream & operator<<(std::ostream &os, const ClusterFlags &clust) {
	for(long i = 0; i < clust.rows(); i++) {
		for(long j = 0; j < clust.cols(i); j++) {
			os << clust(i,j) << '\t';
		}
		os << endl;
	}
	return os;
}