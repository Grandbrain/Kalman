#ifndef KFILTER_H
#define KFILTER_H

#include <kEFilter.h>

namespace Kalman
{
    template<typename T>
    class Filter : public EFilter<T>
    {
    public:

        virtual ~Filter() = 0;

    protected:

        virtual void makeBaseB();
        virtual void makeB();
        virtual void makeProcess();
        virtual void makeMeasure();
        virtual void sizeUpdate();

        Vector px;
        Matrix B;
    };


    template<typename T>
    Filter<T>::~Filter()
    {
    }

    template<typename T>
    void Filter<T>::makeBaseB()
    {
    }

    template<typename T>
    void Filter<T>::makeB()
    {
    }

    template<typename T>
    void Filter<T>::makeProcess()
    {
        makeB();

        px.resize(n);

        for (unsigned i = 0; i < n; ++i)
        {
          px(i) = T(0.0);

          for (unsigned j = 0; j < n; ++j)
            px(i) += A(i,j) * x(j);

          for (unsigned j = 0; j < nu; ++j)
            px(i) += B(i,j) * u(j);

        }

        x.swap(px);

    }

    template<typename T>
    void Filter<T>::makeMeasure()
    {
        z.resize(m);
        for (unsigned i = 0; i < m; ++i)
        {
          z(i) = T(0.0);

          for (unsigned j = 0; j < n; ++j)
            z(i) += H(i,j) * x(j);
        }
    }

    template<typename T>
    void Filter<T>::sizeUpdate()
    {
        if (flags & ( KALMAN_N_MODIFIED | KALMAN_NU_MODIFIED ) )
        {
            B.resize(n, nu);
            makeBaseB();
        }

        EFilter<T>::sizeUpdate();
    }
}

#endif
