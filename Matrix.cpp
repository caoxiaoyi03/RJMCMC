#include <iostream>
#include <cmath>
#include "Matrix.h"
#include <stdexcept>
#include "Random.h"
using namespace std;
//using namespace std::tr1;


#define SMALLNUM 1e-200

////////////////////////// Class basic functions   //////////////////////
//// construct
matrix::matrix(): colNum(0)
{
}

matrix::matrix(long cols, bool isEmpty, long row_reserves): colNum(cols) {
	if(!isEmpty) {
		data = vector<dvector>(row_reserves, dvector(cols));
	} else {
		data.reserve(row_reserves);
	}
}

matrix::matrix(long rows, long cols, double value)
	: data(vector<dvector>(rows)), colNum(cols) {
	for(vector<dvector>::iterator itor = data.begin(); itor != data.end(); itor++) {
		itor->resize(cols, value);
	}
}

matrix::matrix(long rows, long cols, const double* value)
	: data(vector<dvector>(rows)), colNum(cols) {
	for(long i = 0; i < rows; i++) {
		data.at(i).resize(cols);
		for(long j = 0; j < cols; j++) {
			data.at(i).at(j) = value[i * cols + j];
		}
	}
}

matrix::matrix(long rows, long cols, const dvector &value)
	: data(vector<dvector>(rows)), colNum(cols) {
		for(long i=0; i<rows; i++){
		   data.at(i).resize(cols);
		   for(long j=0; j<cols; j++){
		      data.at(i).at(j) = value.at(i * cols + j);
		   
		   }
		}
}


matrix::matrix(long rows, long cols, const double** value) 
: data(vector<dvector>(rows)), colNum(cols) {
	for(long i = 0; i < rows; i++) {
		data.at(i).resize(cols);
		for(long j = 0; j < cols; j++) {
			data.at(i).at(j) = value[i][j];
		}
	}
}


// Identical matrix
matrix::matrix(long rows, double value, bool x)
:data(vector<dvector>(rows)), colNum(rows) {
	for(long i=0; i<rows; i++){
	    data.at(i).resize(rows);
		for(long j=0; j<rows; j++){
		   data.at(i).at(j) = ((i==j)? value:0.0);
		}
	}
}


matrix::matrix(const vector<dvector> &value): data(value), colNum(0) {
	if(value.size()) {
		colNum = value[0].size();
	}
}



matrix::matrix(const matrix &oldmatrix): data(oldmatrix.data), colNum(oldmatrix.colNum) {
}

matrix::matrix(const dvector &vec, bool isRow)
	: data(vector<dvector>(isRow? 1: vec.size())), colNum(isRow? vec.size(): 1) {
		if(isRow) {
			data.front() = vec;
		} else {
			for(long i = 0; i < (signed) vec.size(); i++) {
				data.at(i).resize(1);
				data.at(i).front() = vec.at(i);
			}
		}
}

matrix::~matrix() {
}

double matrix::getValue(long row, long col) const {
	return data.at(row).at(col);
}

void matrix::setValue(long row, long col, double value) {
	data.at(row).at(col) = value;
}

vector<dvector> matrix::getData()const{
    return data;
}

double matrix::operator()(long row, long col) const {
	return getValue(row, col);
}

double & matrix::operator()(long row, long col) {
	return data.at(row).at(col);
}


bool matrix::assertEqualDim(const matrix &mat, bool throwexcept) const {
	if(rows() != mat.rows() || cols() != mat.cols()) {
		if(throwexcept) {
			throw(logic_error("Matrices dimension not equal!"));
		}
		return false;
	}
	return true;
}


matrix &matrix::operator=(const matrix &mat) {
	data = mat.data;
	colNum = mat.colNum;
	return *this;
}



//  decide equal
bool matrix::operator==(const matrix &mat) const {
	if(!assertEqualDim(mat, false)) {
		// different in dimension
		return false;
	}
	for(long i = 0; i < rows(); i++) {
		for(long j = 0; j < cols(); j++) {
			if((*this)(i, j) != mat(i, j)) {
				return false;
			}
		}
	}
	return true;
}

// operator add
matrix matrix::operator+(const matrix &mat) const {
	return (matrix(*this) += mat);
}

matrix matrix::operator+(double value) const {
	return (matrix(*this) += value);
}

// operator substract
matrix matrix::operator-(const matrix &mat) const {
	return (matrix(*this) -= mat);
}

matrix matrix::operator-(double value) const {
	return (matrix(*this) -= value);
}

// operator multiply
matrix matrix::operator*(const matrix &mat) const {
	matrix result(rows(), mat.cols(), 0.0L);
	if(cols() != mat.rows()) {
		throw(logic_error("Matrices dimension not multiplicable!"));
	}
	for(long i = 0; i < rows(); i++) {
		for(long j = 0; j < mat.cols(); j++) {
			for(long traverse = 0; traverse < cols(); traverse++) {
				result(i, j) += (*this)(i, traverse) * mat(traverse, j);
			}
		}
	}
	return result;
}

matrix matrix::operator*(double value) const {
	return (matrix(*this) *= value);
}

// negative
matrix matrix::operator-() const {
	return (*this) * -1L;
}

// operator addto
matrix &matrix::operator+=(const matrix &mat) {
	assertEqualDim(mat);
	for(long i = 0; i < rows(); i++) {
		for(long j = 0; j < cols(); j++) {
			(*this)(i, j) += mat(i, j);
		}
	}
	return (*this);
}

matrix &matrix::operator+=(double value) {
	for(long i = 0; i < rows(); i++) {
		for(long j = 0; j < cols(); j++) {
			(*this)(i, j) += value;
		}
	}
	return (*this);
}

// operator substractto
matrix &matrix::operator-=(const matrix &mat) {
	assertEqualDim(mat);
	for(long i = 0; i < rows(); i++) {
		for(long j = 0; j < cols(); j++) {
			(*this)(i, j) -= mat(i, j);
		}
	}
	return (*this);
}

matrix &matrix::operator-=(double value) {
	for(long i = 0; i < rows(); i++) {
		for(long j = 0; j < cols(); j++) {
			(*this)(i, j) -= value;
		}
	}
	return (*this);
}

// operator multiplyto
matrix &matrix::operator*=(const matrix &mat) {
	data = ((*this) * mat).data;
	return *this;
}

matrix &matrix::operator*=(double value) {
	for(long i = 0; i < rows(); i++) {
		for(long j = 0; j < cols(); j++) {
			(*this)(i, j) *= value;
		}
	}
	return (*this);
}

matrix matrix::transTimesSelf() const {
	matrix result(colNum, colNum, 0.0L);
	for(long i = 0; i < colNum; i++) {
		for(long j = 0; j < colNum; j++) {
			for(long traverse = 0; traverse < rows(); traverse++) {
				result(i, j) += (*this)(traverse, i) * (*this)(traverse, j);
			}
		}
	}
	return result;
}

matrix matrix:: rbind(const matrix &mat) const{
 
	long i,j;
	matrix one(data.size()+mat.rows(),mat.rows());
	if(mat.cols()==data.at(0).size()){
	   for(i=0; i<data.size(); i++){
		  for(j=0; j<mat.cols(); j++){
			  //one.getData().at(i).push_back(data.at(i).at(j));
			  one.setValue(i,j,data.at(i).at(j));
		  }
	  }
	    for(i=data.size(); i<one.getData().size(); i++){
		  for(j=0; j<mat.cols(); j++){
			  one.setValue(i,j,mat.getValue(i-data.size(),j));
			  //one.getData().at(i).push_back(mat.getValue(i-data.size(),j));
		  }
	  }
	} else {
	  cout<<"These two matrix can not be merged!!!"<<endl;
	}


	return one;
}


// transpose
matrix matrix::transpose() const {
	matrix newmatrix(cols(), rows());
	for(long i = 0; i < rows(); i++) {
		for(long j = 0; j < cols(); j++) {
			newmatrix(j, i) = (*this)(i, j);
		}
	}
	return newmatrix;
}

double matrix::determinant() const {
	const long n = rows();
	const long m = cols();
	matrix newmatrix(*this);
	const double TINY=1.0e-20;
	long i, imax, j, k;
	double big, dum, sum, temp;
	double *vv = new double[n];
	double d=1.0;
	if(n == m)
	{
		for (i=0;i<n;i++) {
			big=0.0;
			for (j=0;j<n;j++) {
				if ((temp=fabs(newmatrix(i, j))) > big) {
					big=temp;
				}
			}
			if (big == 0.0) {
				cerr << *this;
				throw(std::logic_error("Singular matrix in routine ludcmp."));
			}
			vv[i]=1.0/big;
		}
		for (j=0;j<n;j++) {
			for (i=0;i<j;i++) {
				sum=newmatrix(i, j);
				for (k=0;k<i;k++) sum -= newmatrix(i, k) * newmatrix(k, j);
				newmatrix(i, j)=sum;
			}
			big=0.0;
			for (i=j;i<n;i++) {
				sum=newmatrix(i, j);
				for (k=0;k<j;k++) sum -= newmatrix(i, k) * newmatrix(k, j);
				newmatrix(i, j)=sum;
				if ((dum=vv[i]*fabs(sum)) >= big) {
					big=dum;
					imax=i;
				}
			}
			if (j != imax) {
				for (k=0;k<n;k++) {
					dum=newmatrix(imax, k);
					newmatrix(imax, k)=newmatrix(j, k);
					newmatrix(j, k)=dum;
				}
				d = -d;
				vv[imax]=vv[j];
			}
			//indx[j]=imax;
			if (newmatrix(j, j) == 0.0) newmatrix(j, j)=TINY;
			if (j != n-1) {
				dum=1.0/(newmatrix(j, j));
				for (i=j+1;i<n;i++) newmatrix(i, j) *= dum;
			}
		}
		delete []vv;
		vv = NULL;
		for(i=0; i<n; i++)
		{
			d *= newmatrix(i, i);
		}
		return d;
	} else {
		throw(logic_error("Matrix is not a square matrix!"));
	}
	return 0;
}


matrix matrix::matrix_mean()const{
	long i,j;
	long n = data.size();
	long m = data.at(0).size();
	matrix average(1,m);
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
		    average(0,i) += data.at(j).at(i);
		}
		average(0,i) = average(0,i)/((double)n);
	}
	return average;
   
}


long matrix::rows() const {
	return data.size();
}

long matrix::cols() const {
	return colNum;
}

matrix matrix::getRow(const long index) const {
	matrix one(1,data.at(index).size(),data.at(index));
	return one;
}

const dvector &matrix::getRowVec(const long index) const {
	return data.at(index);
}

long matrix::addRow(const dvector& row) {
	data.push_back(row);
	return data.size();
}

matrix matrix::getCol(const long index) const {
	dvector newvec(rows());
	for(long i = 0; i < rows(); i++) {
		newvec.at(i) = (*this)(i, index);
	}
	matrix one(newvec.size(),1,newvec);
	return one;
}

// Note that row_begin relates to the row_begin row in the original matrix.
matrix matrix::getBlock(long row_begin, long row_end, long col_begin, long col_end){
	row_begin = row_begin - 1;
	row_end = row_end - 1;
	col_begin = col_begin - 1;
	col_end = col_end - 1;
	if((row_end<=data.size()-1)&&(col_end<=data.at(0).size()-1)&&(row_begin>=0)&&(col_begin>=0)&&(row_begin<=row_end)&&(col_begin<=col_end)){
	    long n = row_end-row_begin+1;
		long m = col_end-col_begin+1;
		matrix one(n, m);
		for(long i=0; i<n; i++)
			for(long j=0; j<m; j++)
				one.setValue(i,j,data.at(i+row_begin).at(j+col_begin));
		return one;
	
	
	
	} else {
		throw(std::logic_error("Can not get this matrix block!!!"));
      //return NULL;	
	}
}

matrix matrix::datasum()
{
 matrix one(1,data.at(0).size());
 //dvector one(data.at(0).size());
 int i;
 for(i=0; i<data.size(); i++)
	 one += data.at(i);
 return one;
}

std::vector<long> matrix::dim() const {
	vector<long> dimvec(2);
	dimvec.at(0) = rows();
	dimvec.at(1) = cols();
	return dimvec;
}



std::ostream & operator<<(std::ostream &os, const matrix &mat) {
	for(long i = 0; i < mat.rows(); i++) {
		for(long j = 0; j < mat.cols(); j++) {
			os << mat(i, j) << ' ';
		}
		os << endl;
	}
	return os;

}


// Return the number of nonzeor elements in each matrix row.
std::vector<long> matrix::l_nzero() const{
	std::vector<long> one(data.size());
    long i,j;
	long n,m;
	n = data.size();
	m = data.at(0).size();
	for(i=0; i<n; i++){
	   one.at(i)=0;
	   for(j=0; j<m; j++){
	      //one.at(i)=(data.at(i).at(j)==0)?one.at(i):(one.at(i)+1);
		   one.at(i)=(data.at(i).at(j)<(SMALLNUM))?one.at(i):(one.at(i)+1);
	   }
	}
	return one;

}




// Return all cluster index of 'time' point???????????
std::vector<long> matrix::tclusterindex(const long ntime) const{
	std::vector<long> one;
	long i,n;
	
	if((ntime>=1)&&(ntime<=data.size())){
	  n = data.at(ntime-1).size();
	  for(i=0; i<n; i++){
		  if(data.at(ntime-1).at(i)>SMALLNUM){
			  one.push_back(i+1);
		  }
	  }
	}else
	{
	 std::cout<<"ntime is out of range."<<std::endl;
	}

	return one;
}





// Get the earliest time point for a specific flag (Real time point).
// If not found, return -1
// Note that we should not use clusterflags to determine the earliest appearing time for some flag.
long matrix::getEarliestTime(long flag) const {
	long n = data.size();
	long m = data.at(0).size();
	for(long i = 0; i < n; i++) {
		if(data.at(i).at(flag-1) > SMALLNUM)
		  return i+1;


		/*for(long j = 0; j < m; j++) {
			if(data.at(i).at(j) == flag) {
				return i+1;
			}
		}*/
	}
	return -1;
}

// Find a position in weight matrix to set an empty cluster.
// The return value is "realposition - 1".
long matrix::getEndGap(long time) const {
	long n = data.at(time-1).size();
	long i;
	for(i = n-1; i>=0; i--){
	    if(data.at(time-1).at(i) > SMALLNUM)
		   return i+1;
	}
	return -1;

}


long matrix::tclusternumber(long time)const{
   long n = data.at(0).size();
   long i,s=0;
   for(i=0; i<n; i++){
	   if(data.at(time-1).at(i)>SMALLNUM){
	      s++;
	   }
   }
   return s;

}

long matrix::clusternumber() const{
	long n = data.at(0).size();
	long m = data.size();
	long i,j, s=0;
	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			if(data.at(j).at(i)>SMALLNUM){
			   s++;
			}
		}
	}
    return s;
}






////////////////////////// Class matrixl basic functions   //////////////////////
//// construct
matrixl::matrixl()
{

}

matrixl::matrixl(long rows, long cols, long value)
	: data(vector<lvector>(rows)) {
	for(vector<lvector>::iterator itor = data.begin(); itor != data.end(); itor++) {
		itor->resize(cols, value);
	}
}

matrixl::matrixl(long rows, long cols, const long* value)
	: data(vector<lvector>(rows)) {
	for(long i = 0; i < rows; i++) {
		data.at(i).resize(cols);
		for(long j = 0; j < cols; j++) {
			data.at(i).at(j) = value[i * cols + j];
		}
	}
}

matrixl::matrixl(long rows, long cols, const lvector &value)
	:data(vector<lvector>(rows)){
		for(long i=0; i<rows; i++){
		   data.at(i).resize(cols);
		   for(long j=0; j<cols; j++){
		      data.at(i).at(j) = value.at(i * cols + j);
		   
		   }
		}
}


matrixl::matrixl(long rows, long cols, const long** value) 
:data(vector<lvector>(rows)){
	for(long i = 0; i < rows; i++) {
		data.at(i).resize(cols);
		for(long j = 0; j < cols; j++) {
			data.at(i).at(j) = value[i][j];
		}
	}
}


// Identical matrix
matrixl::matrixl(long rows, long value, bool x)
:data(vector<lvector>(rows)){
	for(long i=0; i<rows; i++){
	    data.at(i).resize(rows);
		for(long j=0; j<rows; j++){
		   data.at(i).at(j) = ((i==j)? value:0.0);
		}
	}
}


matrixl::matrixl(const vector<lvector> &value): data(value){
}



matrixl::matrixl(const matrixl &oldmatrix): data(oldmatrix.data) {
}

matrixl::matrixl(const lvector &vec, bool isRow)
	: data(vector<lvector>(isRow? 1: vec.size())) {
		if(isRow) {
			data.front() = vec;
		} else {
			for(long i = 0; i < (signed) vec.size(); i++) {
				data.at(i).resize(1);
				data.at(i).front() = vec.at(i);
			}
		}
}

matrixl::~matrixl() {
}

long matrixl::getValue(long row, long col) const {
	return data.at(row).at(col);
}

void matrixl::setValue(long row, long col, long value) {
	data.at(row).at(col) = value;
}

vector<lvector> matrixl::getData()const{
    return data;
}

long matrixl::operator()(long row, long col) const {
	return getValue(row, col);
}

long & matrixl::operator()(long row, long col) {
	return data.at(row).at(col);
}


bool matrixl::assertEqualDim(const matrixl &mat, bool throwexcept) const {
	if(rows() != mat.rows() || cols() != mat.cols()) {
		if(throwexcept) {
			throw(std::logic_error("Matrices dimension not equal!"));
		}
		return false;
	}
	return true;
}


matrixl &matrixl::operator=(const matrixl &mat) {
	data = mat.data;
	return *this;
}



//  decide equal
bool matrixl::operator==(const matrixl &mat) const {
	if(!assertEqualDim(mat, false)) {
		// different in dimension
		return false;
	}
	for(long i = 0; i < rows(); i++) {
		for(long j = 0; j < cols(); j++) {
			if((*this)(i, j) != mat(i, j)) {
				return false;
			}
		}
	}
	return true;
}

// operator add
matrixl matrixl::operator+(const matrixl &mat) const {
	return (matrixl(*this) += mat);
}

matrixl matrixl::operator+(long value) const {
	return (matrixl(*this) += value);
}

// operator substract
matrixl matrixl::operator-(const matrixl &mat) const {
	return (matrixl(*this) -= mat);
}

matrixl matrixl::operator-(long value) const {
	return (matrixl(*this) -= value);
}

// operator multiply
matrixl matrixl::operator*(const matrixl &mat) const {
	matrixl result(rows(), mat.cols(), 0L);
	if(cols() != mat.rows()) {
		throw(logic_error("Matrices dimension not multiplicable!"));
	}
	for(long i = 0; i < rows(); i++) {
		for(long j = 0; j < mat.cols(); j++) {
			for(long traverse = 0; traverse < cols(); traverse++) {
				result(i, j) += (*this)(i, traverse) * mat(traverse, j);
			}
		}
	}
	return result;
}

matrixl matrixl::operator*(long value) const {
	return (matrixl(*this) *= value);
}

// negative
matrixl matrixl::operator-() const {
	return (*this) * -1L;
}

// operator addto
matrixl &matrixl::operator+=(const matrixl &mat) {
	assertEqualDim(mat);
	for(long i = 0; i < rows(); i++) {
		for(long j = 0; j < cols(); j++) {
			(*this)(i, j) += mat(i, j);
		}
	}
	return (*this);
}

matrixl &matrixl::operator+=(long value) {
	for(long i = 0; i < rows(); i++) {
		for(long j = 0; j < cols(); j++) {
			(*this)(i, j) += value;
		}
	}
	return (*this);
}

// operator substractto
matrixl &matrixl::operator-=(const matrixl &mat) {
	assertEqualDim(mat);
	for(long i = 0; i < rows(); i++) {
		for(long j = 0; j < cols(); j++) {
			(*this)(i, j) -= mat(i, j);
		}
	}
	return (*this);
}

matrixl &matrixl::operator-=(long value) {
	for(long i = 0; i < rows(); i++) {
		for(long j = 0; j < cols(); j++) {
			(*this)(i, j) -= value;
		}
	}
	return (*this);
}

// operator multiplyto
matrixl &matrixl::operator*=(const matrixl &mat) {
	data = ((*this) * mat).data;
	return *this;
}

matrixl &matrixl::operator*=(long value) {
	for(long i = 0; i < rows(); i++) {
		for(long j = 0; j < cols(); j++) {
			(*this)(i, j) *= value;
		}
	}
	return (*this);
}

matrixl matrixl:: rbind(const matrixl &mat) const{
 
	long i,j;
	matrixl one(data.size()+mat.rows(),mat.rows());
	if(mat.cols()==data.at(0).size()){
	   for(i=0; i<data.size(); i++){
		  for(j=0; j<mat.cols(); j++){
			  //one.getData().at(i).push_back(data.at(i).at(j));
			  one.setValue(i,j,data.at(i).at(j));
		  }
	  }
	    for(i=data.size(); i<one.getData().size(); i++){
		  for(j=0; j<mat.cols(); j++){
			  one.setValue(i,j,mat.getValue(i-data.size(),j));
			  //one.getData().at(i).push_back(mat.getValue(i-data.size(),j));
		  }
	  }
	} else {
	  cout<<"These two matrix can not be merged!!!"<<endl;
	}


	return one;
}


// transpose
matrixl matrixl::transpose() const {
	matrixl newmatrix(cols(), rows());
	for(long i = 0; i < rows(); i++) {
		for(long j = 0; j < cols(); j++) {
			newmatrix(j, i) = (*this)(i, j);
		}
	}
	return newmatrix;
}



//matrix matrix::inverse() const {
//}
//
long matrixl::rows() const {
	return data.size();
}

long matrixl::cols() const {
	if(data.size() == 0) {
		return 0;
	}
	return data.at(0).size();
}

matrixl matrixl::getRow(const long index) const {
	matrixl one(1,data.at(index).size(),data.at(index));
	return one;
}

matrixl matrixl::ran_getRow()const{
	return this->getRow(ran_iunif(0, this->rows() - 1));
}


matrixl matrixl::getCol(const long index) const {
	lvector newvec(rows());
	for(long i = 0; i < rows(); i++) {
		newvec.at(i) = (*this)(i, index);
	}
	matrixl one(newvec.size(),1,newvec);
	return one;
}

// Note that row_begin relates to the row_begin row in the original matrix.
matrixl matrixl::getBlock(long row_begin, long row_end, long col_begin, long col_end){
	row_begin = row_begin - 1;
	row_end = row_end - 1;
	col_begin = col_begin - 1;
	col_end = col_end - 1;
	if((row_end<=data.size()-1)&&(col_end<=data.at(0).size()-1)&&(row_begin>=0)&&(col_begin>=0)&&(row_begin<=row_end)&&(col_begin<=col_end)){
	    long n = row_end-row_begin+1;
		long m = col_end-col_begin+1;
		matrixl one(n, m);
		for(long i=0; i<n; i++)
			for(long j=0; j<m; j++)
				one.setValue(i,j,data.at(i+row_begin).at(j+col_begin));
		return one;
	
	
	
	} else {
		throw(std::logic_error("Can not get this matrix block!!!"));
      //return NULL;	
	}
}

matrixl matrixl::datasum()
{
 matrixl one(1,data.at(0).size());
 //dvector one(data.at(0).size());
 int i;
 for(i=0; i<data.size(); i++)
	 one += data.at(i);
 return one;
}

std::vector<long> matrixl::dim() const {
	vector<long> dimvec(2);
	dimvec.at(0) = rows();
	dimvec.at(1) = cols();
	return dimvec;
}



std::ostream & operator<<(std::ostream &os, const matrixl &mat) {
	for(long i = 0; i < mat.rows(); i++) {
		for(long j = 0; j < mat.cols(); j++) {
			os << mat(i, j) << ' ';
		}
		os << endl;
	}
	return os;

}


// Return the number of nonzeor elements in each matrix row.
std::vector<long> matrixl::l_nzero() const{
	std::vector<long> one(data.size());
    long i,j;
	long n,m;
	n = data.size();
	m = data.at(0).size();
	for(i=0; i<n; i++){
	   one.at(i)=0;
	   for(j=0; j<m; j++){
	      //one.at(i)=(data.at(i).at(j)==0)?one.at(i):(one.at(i)+1);
		   one.at(i)=(data.at(i).at(j)==0)?one.at(i):(one.at(i)+1);
	   }
	}
	return one;

}




// Get the earliest time point for a specific flag (Real time point).
// If not found, return -1
// Note that we should not use clusterflags to determine the earliest appearing time for some flag.
long matrixl::getEarliestTime(long flag) const {
	long n = data.size();
	long m = data.at(0).size();
	for(long i = 0; i < n; i++) {
		if(data.at(i).at(flag-1) > SMALLNUM)
		  return i+1;


		/*for(long j = 0; j < m; j++) {
			if(data.at(i).at(j) == flag) {
				return i+1;
			}
		}*/
	}
	return -1;
}

// Find a position in weight matrix to set an empty cluster.
// The return value is "realposition - 1".
long matrixl::getEndGap(long time) const {
	long n = data.at(time-1).size();
	long i;
	for(i = n-1; i>=0; i--){
	    if(data.at(time-1).at(i) > SMALLNUM)
		   return i+1;
	}
	return -1;

}




double matrixl::treesummary()const{
   long n = data.size();
   long m = data.at(0).size();
   long i,j;
   double s = 0.0;
   for(j=0; j<m; j++){
	   for(i=0; i<n; i++){
		   if(data.at(i).at(j) != 0){
		      s += pow(10.0,j)*(1.0+0.5*(double)i);
		   }
	   }
   }
   return s;
}






#undef SMALLNUM