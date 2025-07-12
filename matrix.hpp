#include <sstream>
namespace matrix {
template <class data_t>
cmatr<data_t> cmatr<data_t>::operator*(const cmatr<data_t>& A) const {
    if (size_stl != A.size_str) {
        throw std::exception("can't multiplay matrices");
    }
    cmatr<data_t> newmatr(size_str, A.size_stl);
    for (int j = 0; j < size_str; ++j) {
        for (int i = 0; i < size_stl; ++i) {
            data_t temp = 0;
            for (int k = 0; k < size_stl; ++k) {
                temp += (*this)[j][k] * A[k][i];
            }
            newmatr[j][i] = temp;
        }
    }
    return newmatr;
}

template <class data_t>
std::ostream& operator<<(std::ostream& os, const cmatr<data_t>& a) {
    os << "Matrix print..." << std::endl;
    for (int i = 0; i < a.get_size_str(); i++) {
        for (int j = 0; j < a.get_size_stl(); j++)
            os << a[i][j] << " ";
        os << std::endl;
    }
    os << "---------------" << std::endl;
    return os;
}
template <class data_t>
std::vector<data_t>& cmatr<data_t>::operator[](int i) {
    return arr[i];
}

template <class data_t>
std::vector<data_t> const& cmatr<data_t>::operator[](int i) const {
    return arr[i];
}

template <class data_t>
cmatr<data_t>& cmatr<data_t>::operator+=(const cmatr& b) {
    int i, j;
    if (size_str < b.size_str || size_stl < b.size_stl) {
        throw std::exception("can't plus= matrices");
    }
    for (i = 0; i < b.size_str; i++)
        for (j = 0; j < b.size_stl; j++)
            (*this)[i][j] += b[i][j];
    return (*this);
}

template <class data_t>
cmatr<data_t>& cmatr<data_t>::operator-=(const cmatr& b) {
    int i, j;
    if (size_str < b.size_str || size_stl < b.size_stl) {
        throw std::exception("can't plus= matrices");
    }
    for (i = 0; i < b.size_str; i++)
        for (j = 0; j < b.size_stl; j++)
            (*this)[i][j] -= b[i][j];
    return (*this);
}

template <class data_t>
cmatr<data_t>& cmatr<data_t>::operator=(const cmatr& a) {
    if (size_str < a.size_str || size_stl < a.size_stl) {
        throw std::exception("can't use operator=");
    }
    for (int i = 0; i < a.size_str; i++)
        for (int j = 0; j < a.size_stl; j++)
            (*this)[i][j] = a[i][j];
    return *this;
}

template <class data_t>
cmatr<data_t>& cmatr<data_t>::operator*=(data_t n) {
    for (int i = 0; i < size_str; i++)
        for (int j = 0; j < size_stl; j++)
            (*this)[i][j] *= n;
    return (*this);
}

template <class data_t>
cmatr<data_t>& cmatr<data_t>::operator/=(data_t n) {
    for (int i = 0; i < size_str; i++)
        for (int j = 0; j < size_stl; j++)
            (*this)[i][j] /= n;
    return *this;
}

template <class data_t>
const cmatr<data_t> cmatr<data_t>::operator+(const cmatr& b) const {
    return data_t(*this) += b;
}

template <class data_t>
const cmatr<data_t> cmatr<data_t>::operator-(const cmatr& b) const {
    return data_t(*this) -= b;
}

template <class data_t>
const cmatr<data_t> cmatr<data_t>::operator/(data_t b) const {
    return data_t(*this) /= b;
}

template <class data_t>
const cmatr<data_t> cmatr<data_t>::operator*(data_t b) const {
    return data_t(*this) *= b;
}

template <class data_t>
cmatr<data_t>::cmatr(const cmatr& a) {
    size_str = a.size_str;
    size_stl = a.size_stl;
    get_mem();
    for (int i = 0; i < size_str; i++)
        for (int j = 0; j < size_stl; j++)
            (*this)[i][j] = a[i][j];
}

template <class data_t>
cmatr<data_t>::~cmatr() {
}

template <class data_t>
const cmatr<data_t> cmatr<data_t>::transponate(void) const {
    cmatr<data_t> newmatr(size_stl, size_str);
    for (int i = 0; i < size_str; i++)
        for (int j = 0; j < size_stl; j++)
            newmatr[j][i] = (*this)[i][j];
    return newmatr;
}

template <class data_t>
cmatr<data_t>::cmatr(const char* file_name) {
    std::ifstream is(file_name);
    if (!is) {
        std::string s("can't open matrix file ");
        s += file_name;
        throw std::exception(s.c_str());
    }
    size_str = 0;
    is >> (*this);
}

template <class data_t>
std::istream& operator>>(std::istream& is, cmatr<data_t>& a) {
    while (is) {
        std::vector<data_t> v;

        const int str_size = 32 * 1024;
        char str[str_size];
        is.getline(str, str_size - 1);
        if (!is)
            break;
        std::istringstream istr(str);
        a.size_stl = 0;
        while (istr) {
            data_t temp;
            istr >> temp;
            if (!istr)
                break;
            v.push_back(temp);
            ++(a.size_stl);
        }
        a.arr.push_back(v);
        a.size_str++;
    }
    return is;
}

template <class data_t>
cmatr<data_t>::cmatr(int newsize_str, int newsize_stl) {
    if (newsize_str == 0)
        newsize_str = default_matrix_size;
    if (newsize_stl == 0)
        newsize_stl = newsize_str;
    size_str = newsize_str;
    size_stl = newsize_stl;
    get_mem();
    set_zero();
}

template <class data_t>
void cmatr<data_t>::set_zero() {
    for (int i = 0; i < size_str; i++)
        for (int j = 0; j < size_stl; j++)
            (*this)[i][j] = 0;
}
template <class data_t>
void cmatr<data_t>::get_mem(void) {
    arr.resize(size_str);
    for (int i = 0; i < size_str; i++) {
        arr[i].resize(size_stl);
    }
}

} // namespace matrix