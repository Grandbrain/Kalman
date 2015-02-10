#ifndef KVECTOR_H
#define VECTOR_H

#include <vector>

namespace Kalman
{
    template<typename T>
    class Vector
    {
    public:

        inline Vector();
        inline Vector(unsigned);
        inline Vector(unsigned, const T&);

    public:

        inline T& operator()(unsigned);
        inline const T& operator() (unsigned) const;

    public:

        inline unsigned Size() const;
        inline void Resize(unsigned);
        inline void Resize(unsigned, const T&);
        inline void Swap(Vector&);

    private:

        std::vector<T> data;
    };


    template<typename T>
    inline Vector<T>::Vector()
    {
    }


    template<typename T>
    inline Vector<T>::Vector(unsigned n) : data(n)
    {
    }


    template<typename T>
    inline Vector<T>::Vector(unsigned n, const T& a) : data(n, a)
    {
    }

    template<typename T>
    inline T& Vector<T>::operator()(unsigned i)
    {
        return data[i];
    }


    template<typename T>
    inline const T& Vector<T>::operator()(unsigned i) const
    {
        return data[i];
    }


    template<typename T>
    inline unsigned Vector<T>::Size() const
    {
        return data.size();
    }


    template<typename T>
    inline void Vector<T>::Resize(unsigned n)
    {
        data.resize(n);
    }


    template<typename T>
    inline void Vector<T>::Resize(unsigned n, const T& a)
    {
        data.resize(n, a);
    }


    template<typename T>
    inline void Vector<T>::Swap(Vector& v)
    {
        data.swap(v.data);
    }
}

#endif