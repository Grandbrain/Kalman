#ifndef BASIC_H
#define BASIC_H

#include "extended.h"

namespace Kalman {

    template <typename T>
    class Basic : public Extended<T> {
    public:
        virtual ~Basic() = 0;

    protected:
        virtual void makeBaseB();
        virtual void makeB();
        virtual void makeProcess();
        virtual void makeMeasure();
        virtual void sizeUpdate();

    protected:
        Vector<T> px;
        Matrix<T> B;
    };


    template <typename T>
    void Basic<T>::makeBaseB() {
    }

    template <typename T>
    void Basic<T>::makeB() {
    }

    template <typename T>
    void Basic<T>::makeProcess() {
        makeB();
        px.Resize(n);

        for (unsigned i = 0; i < n; i++) {
            px(i) = T(0.0);

            for (unsigned j = 0; j < n; j++)
                px(i) += A(i, j) * x(j);

            for (unsigned j = 0; j < nu; j++)
                px(i) += B(i, j) * u(j);
        }

        x.Swap(px);
    }

    template <typename T>
    void Basic<T>::makeMeasure() {
        z.Resize(m);

        for (unsigned i = 0; i < m; i++) {
            z(i) = T(0.0);

            for (unsigned j = 0; j < n; j++)
                z(i) += H(i, j) * x(j);
        }
    }

    template <typename T>
    void Basic<T>::sizeUpdate() {
        if (flags & (ModifiedN | ModifiedNU)) {
            B.Resize(n, nu);
            makeBaseB();
        }

        Extended<T>::sizeUpdate();
    }
}

#endif
