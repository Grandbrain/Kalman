#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

namespace Kalman {

    template <typename T>
    class Matrix {
    public:
        inline Matrix() = default;
        inline Matrix(unsigned, unsigned);
        inline Matrix(unsigned, unsigned, const T&);
        inline Matrix(unsigned, unsigned, const T*);

    public:
        inline T& operator()(unsigned, unsigned);
        inline const T& operator()(unsigned, unsigned) const;

    public:
        inline std::size_t Rows() const;
        inline std::size_t Columns() const;
        inline void Assign(unsigned, unsigned, const T&);
        inline void Assign(unsigned, unsigned, const T*);
        inline void Resize(unsigned, unsigned);
        inline void Swap(Matrix&);

    private:
        std::vector<std::vector<T> > data;
    };


    template <typename T>
    inline Matrix<T>::Matrix(unsigned m, unsigned n) {
        Resize(m, n);
    }

    template <typename T>
    inline Matrix<T>::Matrix(unsigned m, unsigned n, const T& a) {
        Assign(m, n, a);
    }

    template <typename T>
    inline Matrix<T>::Matrix(unsigned m, unsigned n, const T* v) {
        Assign(m, n, v);
    }

    template <typename T>
    inline T& Matrix<T>::operator()(unsigned i, unsigned j) {
        return data[i][j];
    }

    template <typename T>
    inline const T& Matrix<T>::operator()(unsigned i, unsigned j) const {
        return data[i][j];
    }

    template <typename T>
    inline std::size_t Matrix<T>::Rows() const {
        return data.size();
    }

    template <typename T>
    inline std::size_t Matrix<T>::Columns() const {
        return data.empty() ? 0 : data.back().size();
    }

    template <typename T>
    inline void Matrix<T>::Assign(unsigned m, unsigned n, const T& a) {
        data.resize(m);
        for (unsigned i = 0; i < data.size(); i++)
            data[i].assign(n, a);
    }

    template <typename T>
    inline void Matrix<T>::Assign(unsigned m, unsigned n, const T* v) {
        data.resize(m);
        for (unsigned i = 0; i < data.size(); i++) {
            data[i].resize(n);
            for (unsigned j = 0; j < data[i].size(); j++)
                data[i][j] = v[i * m + j];
        }
    }

    template <typename T>
    inline void Matrix<T>::Resize(unsigned m, unsigned n) {
        data.resize(m);
        for (unsigned i = 0; i < data.size(); i++)
            data[i].resize(n);
    }
    
    template <typename T>
    inline void Matrix<T>::Swap(Matrix& M) {
        data.swap(M.data);
    }
}

#endif
