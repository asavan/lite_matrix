#ifndef __MATRIX_H__
#define __MATRIX_H__
#include <vector>
#include <strstream>
#include <string>
#include <fstream>
#include <exception>
namespace matrix {

extern int default_matrix_size;

template <class data_t>
class cmatr
{
	void get_mem(void);
	void set_zero(void);
	int size_stl, size_str; //kol-vo stolbcov, stro4ek
	std::vector<std::vector<data_t> > arr;
public:
	std::vector<data_t> &operator[](int i);
	std::vector<data_t> const & operator[](int i) const;
	int get_size_stl(void) const {return size_stl;}
	int get_size_str(void) const {return size_str;}
	explicit cmatr(int newsize_str = 0, int newsize_stl = 0);	
	explicit cmatr(const char *filename);
	
	cmatr(const cmatr &a);
	~cmatr();
	const cmatr transponate(void) const;

	cmatr& operator+=(const cmatr &b);
	const cmatr operator+(const cmatr &b) const;
	
	cmatr& operator-=(const cmatr &b);
	const cmatr operator-(const cmatr &b) const;


	cmatr& operator/=( data_t n);
	const cmatr operator/( data_t n) const;
	
	cmatr& operator*=(data_t n);
	const cmatr operator*(data_t n) const;
	
	cmatr &operator=(const cmatr &a);
	
	cmatr operator*(const cmatr & a) const;

	template <class data_t>
	friend std::istream& operator>>(std::istream& is, cmatr<data_t>& a);

};
}
#include "matrix.hpp"

#endif
