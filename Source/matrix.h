#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

namespace Kalman
{
    template<typename T>
    class Matrix
    {
    public:

        inline Matrix();
        inline Matrix(unsigned, unsigned);
        inline Matrix(unsigned, unsigned, const T&);
	inline Matrix(unsigned, unsigned, const T*);

    public:

        inline T& operator()(unsigned, unsigned);
        inline const T& operator()(unsigned, unsigned) const;

    public:

        inline unsigned Rows() const;
        inline unsigned Columns() const;
        inline void Resize(unsigned, unsigned);
        inline void Resize(unsigned, unsigned, const T&);
	inline void Resize(unsigned, unsigned, const T*)

    private:

        std::vector<std::vector<T> > data;
    };


    template<typename T>
    inline Matrix<T>::Matrix()
    {
    }


    template<typename T>
    inline Matrix<T>::Matrix(unsigned m, unsigned n)
    {
        Resize(m, n);
    }


    template<typename T>
    inline Matrix<T>::Matrix(unsigned m, unsigned n, const T& a)
    {
        Resize(m, n, a);
    }
	
    template<typename T>
    inline Matrix<T>::Matrix(unsigned m, unsigned n, const T* v)
    {
        Resize(m, n, v);
    }

    template<typename T>
    inline T& Matrix<T>::operator()(unsigned i, unsigned j)
    {
        return data[i][j];
    }


    template<typename T>
    inline const T& Matrix<T>::operator()(unsigned i, unsigned j) const
    {
        return data[i][j];
    }


    template<typename T>
    inline unsigned Matrix<T>::Rows() const
    {
        return data.size();
    }


    template<typename T>
    inline unsigned Matrix<T>::Columns() const
    {
        return data.empty() ? 0 : data.back().size();
    }


    template<typename T>
    inline void Matrix<T>::Resize(unsigned m, unsigned n)
    {
        data.resize(m);
        for (unsigned i = 0; i < data.size(); i++)
            data[i].assign(n);
    }

    template<typename T>
    inline void Matrix<T>::Resize(unsigned m, unsigned n, const T& a)
    {
        data.resize(m);
        for (unsigned i = 0; i < data.size(); i++)
            data[i].assign(n, a);
    }
	
    template<typename T>
    inline void Matrix<T>::Resize(unsigned m, unsigned n, const T* v)
    {
        data.resize(m);
        for (unsigned i = 0; i < data.size(); i++)
	{
            data[i].assign(n);
	    for (unsigned j = 0; j < data[i].size(); i++)
		data[i][j] = v[i * m + j];
	}
    }
}

#endif
