#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <stdexcept>

#ifdef LINUX
#include <tr1/array>
#else
#include <array>
#endif

class dvector;
class lvector;

class matrix
{
	// use vector to implement at first
	// may use array later for speed concern
	std::vector<dvector> data;
	long colNum;
public:
	matrix();
	matrix(long cols, bool isEmpty, long row_reserves);
	matrix(long rows, long cols, double value=0L);
	matrix(long rows, long cols, const double* value);
	matrix(long rows, long cols, const dvector &value);
	matrix(long rows, long cols, const double** value);
	// Identical matrix
	matrix(long rows, double value, bool x);
	matrix(const std::vector<dvector> &value);
	matrix(const matrix &oldmatrix);
	matrix(const dvector &vec, bool isRow = true);
	~matrix();

	double getValue(long row, long col) const;
	void setValue(long row, long col, double value);
	std::vector<dvector> getData()const;

	// try to see if the two matrices have equal dimension, if not return false or throw an exception
	bool assertEqualDim(const matrix &mat, bool throwexcept = true) const throw(std::logic_error);

	// can use matrix(row, col) to get a reference of the value
	double operator()(long row, long col) const;
	double & operator()(long row, long col);

	matrix &operator=(const matrix &mat);

	bool operator==(const matrix &mat) const;

	matrix operator*(const matrix &mat) const;
	matrix operator*(double value) const;
	matrix operator+(const matrix &mat) const;
	matrix operator+(double value) const;
	matrix operator-(const matrix &mat) const;
	matrix operator-(double value) const;

	matrix operator-() const;

	matrix &operator*=(const matrix &mat);
	matrix &operator*=(double value);
	matrix &operator+=(const matrix &mat);
	matrix &operator+=(double value);
	matrix &operator-=(const matrix &mat);
	matrix &operator-=(double value);

	matrix transTimesSelf() const;
	// M.transTimesSelf() = (M)^T * M

	matrix rbind(const matrix &mat) const;
	//matrix rbind(const dvector &vec) const;

	//matrix cbind(const matrix &mat) const;
	//matrix cbind(const dvector &vec) const;

	matrix transpose() const;
	double determinant() const;
	matrix matrix_mean()const;
	//matrix inverse() const;

	long rows() const;
	long cols() const;

	matrix getRow(const long index)const;
	const dvector &getRowVec(const long index) const;
	matrix getCol(const long index)const;
	matrix getBlock(long row_begin, long row_end, long col_begin, long col_end);

	matrix datasum();
	long getEarliestTime(long flag) const;

	std::vector<long> dim() const;
	std::vector<long> l_nzero() const;
	long getEndGap(long time) const;
	long clusternumber() const;
	long tclusternumber(long time)const;
	std::vector<long> tclusterindex(const long time) const;

	long addRow(const dvector& row);
	

	//matrix &fill(double value);
	//matrix &fillRow(long index, double value);
	//matrix &fillCol(long index, double value);

	friend matrix operator+(double, const matrix&);
	friend matrix operator-(double, const matrix&);
	friend matrix operator*(double, const matrix&);

	friend std::ostream & operator<<(std::ostream &os, const matrix &mat);

};

matrix operator+(double, const matrix&);
matrix operator-(double, const matrix&);
matrix operator*(double, const matrix&);

std::ostream & operator<<(std::ostream &os, const matrix &mat);


class dvector: public std::vector<double> {
	//unsigned long size() const;
public:
	dvector(long size): std::vector<double>(size) { }
	dvector() {}
	~dvector() {}
};





class matrixl
{
	// use vector to implement at first
	// may use array later for speed concern
	std::vector<lvector> data;
public:
	matrixl();
	matrixl(long rows, long cols, long value=0L);
	matrixl(long rows, long cols, const long* value);
	matrixl(long rows, long cols, const lvector &value);
	matrixl(long rows, long cols, const long** value);
	// Identical matrix
	matrixl(long rows, long value, bool x);
	matrixl(const std::vector<lvector> &value);
	matrixl(const matrixl &oldmatrix);
	matrixl(const lvector &vec, bool isRow = true);
	~matrixl();

	long getValue(long row, long col) const;
	void setValue(long row, long col, long value);
	std::vector<lvector> getData()const;

	// try to see if the two matrices have equal dimension, if not return false or throw an exception
	bool assertEqualDim(const matrixl &mat, bool throwexcept = true) const throw(std::logic_error);

	// can use matrix(row, col) to get a reference of the value
	long operator()(long row, long col) const;
	long & operator()(long row, long col);

	matrixl &operator=(const matrixl &mat);

	bool operator==(const matrixl &mat) const;

	matrixl operator*(const matrixl &mat) const;
	matrixl operator*(long value) const;
	matrixl operator+(const matrixl &mat) const;
	matrixl operator+(long value) const;
	matrixl operator-(const matrixl &mat) const;
	matrixl operator-(long value) const;

	matrixl operator-() const;

	matrixl &operator*=(const matrixl &mat);
	matrixl &operator*=(long value);
	matrixl &operator+=(const matrixl &mat);
	matrixl &operator+=(long value);
	matrixl &operator-=(const matrixl &mat);
	matrixl &operator-=(long value);

	matrixl rbind(const matrixl &mat) const;
	//matrix rbind(const dvector &vec) const;

	//matrix cbind(const matrix &mat) const;
	//matrix cbind(const dvector &vec) const;

	matrixl transpose() const;

	//matrix inverse() const;

	long rows() const;
	long cols() const;

	matrixl getRow(const long index) const;
	matrixl ran_getRow()const;
	matrixl getCol(const long index) const;
	matrixl getBlock(long row_begin, long row_end, long col_begin, long col_end);

	matrixl datasum();
	long getEarliestTime(long flag) const;

	std::vector<long> dim() const;
	std::vector<long> l_nzero() const;
	
	long getEndGap(long time) const;

	double treesummary()const;

	//matrix &fill(double value);
	//matrix &fillRow(long index, double value);
	//matrix &fillCol(long index, double value);

	friend matrixl operator+(long, const matrixl&);
	friend matrixl operator-(long, const matrixl&);
	friend matrixl operator*(long, const matrixl&);

	friend std::ostream & operator<<(std::ostream &os, const matrixl &mat);

};

matrixl operator+(long, const matrixl&);
matrixl operator-(long, const matrixl&);
matrixl operator*(long, const matrixl&);

std::ostream & operator<<(std::ostream &os, const matrixl &mat);











class lvector: public std::vector<long> {
	//unsigned long size() const;
public:
	lvector(long size): std::vector<long>(size) { }
	lvector() {}
	~lvector() {}
};

#endif