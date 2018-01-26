#ifndef VECTOR_H
#define VECTOR_H

#include <vector>

namespace Kalman {

    template <typename T>
    class Vector {
    public:
        inline Vector() = default;
        inline Vector(const Vector&) = default;
        inline Vector(unsigned);
        inline Vector(unsigned, const T&);
        inline Vector(unsigned, const T*);

    public:
        inline Vector& operator=(const T&);
        inline Vector& operator=(const Vector&) = default;
        inline T& operator()(unsigned);
        inline const T& operator()(unsigned) const;

    public:
        inline std::size_t Size() const;
        inline void Resize(unsigned);
        inline void Resize(unsigned, const T&);
        inline void Resize(unsigned, const T*);
        inline void Swap(Vector&);

    private:
        std::vector<T> data;
    };


    template <typename T>
    inline Vector<T>::Vector(unsigned n) {
        Resize(n);
    }


    template <typename T>
    inline Vector<T>::Vector(unsigned n, const T& a) {
        Resize(n, a);
    }

    template <typename T>
    inline Vector<T>::Vector(unsigned n, const T* v) {
        Resize(n, v);
    }

    template <typename T>
    inline Vector<T>& Vector<T>::operator=(const T& a) {
        std::fill(data.begin(), data.end(), a);
        return *this;
    }

    template <typename T>
    inline T& Vector<T>::operator()(unsigned i) {
        return data[i];
    }


    template <typename T>
    inline const T& Vector<T>::operator()(unsigned i) const {
        return data[i];
    }


    template <typename T>
    inline std::size_t Vector<T>::Size() const {
        return data.size();
    }


    template <typename T>
    inline void Vector<T>::Resize(unsigned n) {
        data.resize(n);
    }


    template <typename T>
    inline void Vector<T>::Resize(unsigned n, const T& a) {
        data.resize(n, a);
    }


    template <typename T>
    inline void Vector<T>::Resize(unsigned n, const T* v) {
        data.resize(n);
        for (unsigned i = 0; i < data.size(); i++)
            data[i] = v[i];
    }


    template <typename T>
    inline void Vector<T>::Swap(Vector& v) {
        data.swap(v.data);
    }
}

#endif
