#ifndef CLUSTERFLAGS_H
#define CLUSTERFLAGS_H
//#include <array>
#include <vector>
#include <iostream>

#ifdef LINUX
#include <tr1/array>
#else
#include <array>
#endif


// ClusterFlags:
// Store the cluster number of the sample
// row = time of the sample
// column = sample number
// value = cluster number 
class ClusterFlags
{
	std::vector<std::vector<long> > data;
public:
	ClusterFlags();
	ClusterFlags(long n_timepoints, const std::vector<unsigned long> &samples);
	ClusterFlags(const std::vector<std::vector<long> > &one);
	~ClusterFlags();
	ClusterFlags &operator=(const ClusterFlags &clu);
	long operator()(long row, long col) const;
	long & operator()(long row, long col);
	
	long rows() const;
	long cols(int i) const;

	void fillFlags(long flag = -1);

	long getFlag(long time, long sample) const;
	void setFlag(long time, long sample, long flag);

	long getEarliestTime(long flag) const;

	long getMaxFlag() const;

	std::vector<std::vector<long> > getTimeSampleByCluster(long flag) const;

	ClusterFlags sortz(const ClusterFlags &one, long n_timepoints, std::vector<long> NT);

	std::vector<std::vector<long> > onezerocluster(long flag) const;

	void reflag();

	long flagnum(long flag) const;
	long tflagnum(long flag, long ntime) const;
	std::vector<long> getTflag(long ntime)const;
	std::vector<long> getTindicator(long ntime)const;
	bool flaginT(long flagx, long ntime) const;
	std::vector<long> getTtoTindicator(long ntime)const;
	long nonemptynum()const;

};


std::ostream & operator<<(std::ostream &os, const ClusterFlags &clust);




#endif