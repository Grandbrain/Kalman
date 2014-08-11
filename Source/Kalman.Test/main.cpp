#ifndef KEFILTER_H
#define KEFILTER_H

#include <kVector.h>
#include <kMatrix.h>

namespace Kalman
{
    template<typename T>
    class EFilter
    {
    public:

        EKFilter();
        EKFilter(unsigned, unsigned, unsigned, unsigned, unsigned);


        virtual ~EKFilter();

        unsigned getSizeX() const;
        unsigned getSizeU() const;
        unsigned getSizeW() const;
        unsigned getSizeZ() const;
        unsigned getSizeV() const;


        void setDim(unsigned, unsigned, unsigned, unsigned, unsigned);
        void setSizeX(unsigned);
        void setSizeU(unsigned);
        void setSizeW(unsigned);
        void setSizeZ(unsigned);
        void setSizeV(unsigned);
        void init(Vector&, Matrix&);
        void step(Vector&, const Vector&);
        void timeUpdateStep(Vector&);
        void measureUpdateStep(const Vector&);
        const Vector& predict(Vector&);
        const Vector& simulate();
        const Vector& getX() const;
        const Matrix& calculateP() const;

  protected:

        virtual void makeBaseA();
        virtual void makeBaseW();
        virtual void makeBaseQ();
        virtual void makeBaseH();
        virtual void makeBaseV();
        virtual void makeBaseR();
        virtual void makeCommonProcess();
        virtual void makeA();
        virtual void makeW();
        virtual void makeQ();
        virtual void makeProcess() = 0;
        virtual void makeCommonMeasure();
        virtual void makeH();
        virtual void makeV();
        virtual void makeR();
        virtual void makeMeasure() = 0;
        virtual void makeDZ();
        virtual void sizeUpdate();


        Vector x;
        Vector u;
        Vector z;
        Vector dz;

        Matrix A;
        Matrix W;
        Matrix Q;
        Matrix H;
        Matrix V;
        Matrix R;


        unsigned n;
        unsigned nu;
        unsigned nw;
        unsigned  m;
        unsigned nv;

    private:


        static void factor(Matrix&);
        static void upperInvert(Matrix&);
        void timeUpdate();
        void measureUpdate(T, T);


        void makeBaseAImpl();
        void makeBaseWImpl();
        void makeBaseQImpl();
        void makeBaseHImpl();
        void makeBaseVImpl();
        void makeBaseRImpl();
        void makeAImpl();
        void makeWImpl();
        void makeQImpl();
        void makeHImpl();
        void makeVImpl();
        void makeRImpl();


        Matrix U;
        Matrix W_;
        Matrix Q_;
        Matrix H_;
        Matrix R_;

        Vector a;
        Vector d;
        Vector v;

        unsigned nn;

        mutable Matrix _P;
        mutable Vector _x;

        unsigned flags;
        bool modified_;
    };


    template<typename T>
    EFilter<T>::EFilter() : flags(0)
    {
    }


    template<typename T>
    EFilter<T>::EFilter(unsigned n_, unsigned nu_, unsigned nw_, unsigned m_, unsigned nv_) : flags(0)
    {
        setDim(n_, nu_, nw_, m_, nv_);
    }


    template<typename T>
    EFilter<T>::~EFilter()
    {
    }


    template<typename T>
    void EFilter<T>::setDim(unsigned n_, unsigned nu_, unsigned nw_, unsigned m_, unsigned nv_)
    {
        setSizeX(n_);
        setSizeU(nu_);
        setSizeW(nw_);
        setSizeZ(m_);
        setSizeV(nv_);
    }


    template<typename T>
    unsigned EFilter<T>::getSizeX() const
    {
        return n;
    }


    template<typename T>
    unsigned EFilter<T>::getSizeU() const
    {
        return nu;
    }


    template<typename T>
    unsigned EFilter<T>::getSizeW() const
    {
        return nw;
    }


    template<typename T>
    unsigned EFilter<T>::getSizeZ() const
    {
        return m;
    }


    template<typename T>
    unsigned EFilter<T>::getSizeV() const
    {
        return nv;
    }


    template<typename T>
    void EFilter<T>::setSizeX(unsigned n_)
    {
        if (n_ == n) return;
        flags |= KALMAN_N_MODIFIED;
        n = n_;
    }

    template<typename T>
    void EFilter<T>::setSizeU(unsigned nu_)
    {
        if (nu_ == nu) return;
        flags |= KALMAN_NU_MODIFIED;
        nu = nu_;
    }

    template<typename T>
    void EFilter<T>::setSizeW(unsigned nw_)
    {
        if (nw_ == nw) return;
        flags |= KALMAN_NW_MODIFIED;
        nw = nw_;
    }


    template<typename T>
    void EFilter<T>::setSizeZ(unsigned m_)
    {
        if (m_ == m) return;
        flags |= KALMAN_M_MODIFIED;
        m = m_;
    }


    template<typename T>
    void EFilter<T>::setSizeV(unsigned nv_)
    {
        if (nv_ == nv)
        flags |= KALMAN_NV_MODIFIED;
        nv = nv_;
    }


    template<typename T>
    void EFilter<T>::init(Vector& x_, Matrix& P_)
    {
        x.swap(x_);
        _P.swap(P_);
        flags |= KALMAN_P_MODIFIED;
    }


    template<typename T>
    void EFilter<T>::step(Vector& u_, const Vector& z_)
    {
        timeUpdateStep(u_);
        measureUpdateStep(z_);
    }


    template<typename T>
    void EFilter<T>::timeUpdateStep(Vector& u_)
    {

        K_UINT_32 i, j, k;

        sizeUpdate();
        u.swap(u_);

        makeCommonProcess();
        makeAImpl();
        makeWImpl();
        makeQImpl();
        makeProcess();

        if (!OQ)
        {
            if (flags & KALMAN_Q_MODIFIED)
            {
                Q_ = Q;
                factor(Q_);
                upperInvert(Q_);
            }

            Q.swap(Q_);

            if (flags & ( KALMAN_W_MODIFIED | KALMAN_Q_MODIFIED ) )
            {
                for (i = 0; i < n; ++i)
                {
                    for (j = 0; j < nw; ++j)
                    {
                        W_(i,j) = W(i,j);
                        for (k = 0; k < j; ++k)
                            W_(i,j) += W(i,k)*Q(j,k);
                    }
                }
            }

            W.swap(W_);
        }

        timeUpdate();

        if (!OQ)
        {
            Q.swap(Q_);
            W.swap(W_);
        }

        u.swap(u_);
        flags &= ~KALMAN_MIDMASK;
    }


    template<typename T>
    void EFilter<T>::measureUpdateStep(const Vector& z_)
    {
        K_UINT_32 i, j, k;
        sizeUpdate();

        if (m == 0) return;

        makeCommonMeasure();
        makeHImpl();
        makeVImpl();
        makeRImpl();
        makeMeasure();


        for (i = 0; i < m; ++i)
          dz(i) = z_(i) - z(i);

        makeDZ();

        if (OVR)
        {
            if (flags & ( KALMAN_V_MODIFIED | KALMAN_R_MODIFIED ) )
            {
                for (i = 0; i < m; ++i)
                    R_(i,i) = V(i,i)*V(i,i)*R(i,i);
            }
        }
        else
        {
            if (flags & ( KALMAN_V_MODIFIED | KALMAN_R_MODIFIED ) )
            {
                _x.resize(nv);
                for (i = 0; i < m; ++i)
                {
                    for (j = 0; j < nv; ++j)
                    {
                        _x(j) = T(0.0);
                        for (k = 0; k < nv; ++k)
                            _x(j) += V(i,k)*R(k,j);
                    }

                    for (j = 0; j < m; ++j)
                    {
                        R_(i,j) = T(0.0);
                        for (k = 0; k < nv; ++k)
                            R_(i,j) += _x(k)*V(j,k);
                    }
                }

                factor(R_);
                upperInvert(R_);
           }

            if (flags & ( KALMAN_H_MODIFIED | KALMAN_V_MODIFIED | KALMAN_R_MODIFIED ) )
            {
                for (i = 0; i < m; ++i)
                {
                    for (j = 0; j < n; ++j)
                    {
                        H_(i,j) = H(i,j);
                        for (k = i + 1; k < m + 0; ++k)
                            H_(i,j) += R_(k,i)*H(k,j);
                    }
                }
            }

            H.swap(H_);
            _x.resize(m);

            for (i = 0; i < m; ++i)
            {
                _x(i) = dz(i);
                for (k = i + 1; k < m; ++k)
                    _x(i) += R_(k,i)*dz(k);
            }

            dz.swap(_x);
        }

        _x.resize(n);
        _x = T(0.0);

        for (i = 0; i < m; ++i)
        {
            for (j = 0; j < n; ++j)
                a(j) = H(i,j);

            measureUpdate(dz(i), R_(i,i));
        }

        for (i = 0; i < n + 0; ++i)
            x(i) += _x(i);

        if (!OVR) H.swap(H_);
        flags &= ~KALMAN_HIGHMASK;
    }


    template<typename T>
    const typename EFilter<T>::Vector& EFilter<T>::predict(Vector& u_)
    {
        sizeUpdate();
        u.swap(u_);
        _x = x;

        makeCommonProcess();
        makeProcess();

        x.swap(_x);
        u.swap(u_);
        return _x;
    }


    template<typename T>
    const typename EFilter<T>::Vector& EFilter<T>::simulate()
    {
        sizeUpdate();
        _x = z;

        makeCommonMeasure();
        makeMeasure();

        z.swap(_x);
        return _x;
    }


    template<typename T>
    const typename EFilter<>::Vector& EFilter<T>::getX() const
    {
        return x;
    }


    template<typename T>
    const typename EFilter<T>::Matrix& EFilter<T>::calculateP() const
    {
        if (!(flags & KALMAN_P_MODIFIED))
        {
            _P.resize(n, n);
            for (unsigned i = 0; i < n; ++i)
            {
                _P(i,i) = U(i,i);
                for (unsigned j = i + 1; j < n; ++j)
                {
                    _P(i,j)  = U(i,j)*U(j,j);
                    _P(i,i) += U(i,j)*_P(i,j);

                    for (unsigned k = j + 1; k < n; ++k)
                        _P(i,j) += U(i,k)*U(j,k)*U(k,k);
                    _P(j,i) = _P(i,j);
                }
            }
        }
        return _P;
    }


    template<typename T>
    void EFilter<T>::makeBaseA()
    {
        modified_ = false;
    }


    template<typename T>
    void EFilter<T>::makeBaseW()
    {
        modified_ = false;
    }


    template<typename T>
    void EFilter<T>::makeBaseQ()
    {
        modified_ = false;
    }


    template<typename T>
    void EFilter<T>::makeBaseH()
    {
        modified_ = false;
    }


    template<typename T>
    void EFilter<T>::makeBaseV()
    {
        modified_ = false;
    }


    template<typename T>
    void EFilter<T>::makeBaseR()
    {
        modified_ = false;
    }


    template<typename T>
    void EFilter<T>::makeCommonProcess()
    {
    }


    template<typename T>
    void EFilter<T>::makeCommonMeasure()
    {
    }


    template<typename T>
    void EFilter<T>::makeA()
    {
        modified_ = false;
    }


    template<typename T>
    void EFilter<T>::makeW()
    {
        modified_ = false;
    }


    template<typename T>
    void EFilter<T>::makeQ()
    {
        modified_ = false;
    }


    template<typename T>
    void EFilter<T>::makeH()
    {
        modified_ = false;
    }


    template<typename T>
    void EFilter<T>::makeV()
    {
        modified_ = false;
    }


    template<typename T>
    void EFilter<T>::makeR()
    {
        modified_ = false;
    }


    template<typename T>
    void EFilter<T>::makeDZ()
    {
    }


    template<typename T>
    void EFilter<T>::sizeUpdate()
    {
        if (!flags) return;
        if (flags & KALMAN_N_MODIFIED)
        {
            A.resize(n, n);
            makeBaseAImpl();
        }

        if (flags & (KALMAN_N_MODIFIED | KALMAN_NW_MODIFIED) )
        {
            nn = n + nw;
            a.resize(nn);
            v.resize(nn);
            d.resize(nn);
            if (!OQ) W_.resize(n, nw);
            W.resize(n, nw);
            makeBaseWImpl();
        }

        if (flags & KALMAN_P_MODIFIED)
        {
            U.resize(n, nn);
            for (K_UINT_32 i = 0; i < n + 0; ++i)
                for (K_UINT_32 j = 0; j < n + 0; ++j)
                    U(i,j) = _P(i,j);

            factor(U);

        }
        else if (flags & KALMAN_NW_MODIFIED)
        {
            _P.resize(n, nn);
            for (K_UINT_32 i = 0; i < n + 0; ++i)
                for (K_UINT_32 j = i; j < n + 0; ++j)
                    _P(i,j) = U(i,j);

            U.swap(_P);
        }

        if (flags & KALMAN_NW_MODIFIED)
        {
            if (!OQ) Q_.resize(nw, nw);
            Q.resize(nw, nw);
            makeBaseQImpl();
        }

        if (m != 0)
        {
            if (flags & (KALMAN_N_MODIFIED | KALMAN_M_MODIFIED) )
            {
                if (!OVR) H_.resize(m, n);
                H.resize(m, n);
                makeBaseHImpl();
            }

            if (flags & (KALMAN_M_MODIFIED | KALMAN_NV_MODIFIED) )
            {
                V.resize(m, nv);
                makeBaseVImpl();
            }

          if (flags & KALMAN_NV_MODIFIED) {
            R.resize(nv, nv);
            makeBaseRImpl();
          }

          if (flags & KALMAN_M_MODIFIED) {
            R_.resize(m, m);
            z.resize(m);
            dz.resize(m);
          }

        }

        flags &= ~KALMAN_LOWMASK;
    }


    template<typename T>
    void EFilter<T>::factor(Matrix& P_)
    {
        T alpha, beta;
        K_UINT_32 i, j, k, N = P_.nrow();
        for(j = N - 1 + 0; j > 0; --j)
        {
            alpha = T(1.0)/P_(j,j);
            for(k = 0; k < j; ++k)
            {
                beta = P_(k,j);
                P_(k,j) = alpha*beta;
                for(i = 0; i <= k; ++i)
                    P_(i,k) -= beta*P_(i,j);
            }
        }
    }


    template<typename T>
    void EFilter<T>::upperInvert(Matrix& P_)
    {
        T val;
        K_UINT_32 i, j, k, N = P_.nrow();
        for (i = N - 2 + 0; i != (K_UINT_32)(0-1); --i)
        {
            for (k = i + 1; k < N + 0; ++k)
            {
                val = P_(i,k);
                for (j = i + 1; j <= k - 1; ++j)
                    val += P_(i,j)*P_(k,j);
                P_(k,i) = -val;
            }
        }
    }


    template<typename T>
    void EFilter<T>::timeUpdate()
    {
        K_UINT_32 i, j, k;
        T sigma, dinv;

        for(j = n - 1 + 0; j > 0; --j)
        {
            for(i = 0; i <= j; ++i)
                d(i) = U(i,j);

            for(i = 0; i < n + 0; ++i)
            {
                U(i,j) = A(i,j);
                for(k = 0; k < j; ++k)
                    U(i,j) += A(i,k)*d(k);
            }
        }

        d(0) = U(0,0);

        for(j = 0; j < n + 0; ++j)
            U(j,0) = A(j,0);

        for(i = 0; i < nw + 0; ++i)
        {
            d(i+n) = Q(i,i);
            for(j = 0; j < n + 0; ++j)
                U(j,i+n) = W(j,i);
        }

        for(j = n - 1 + 0; j != (K_UINT_32)(0-1); --j)
        {
            sigma = T(0.0);
            for(k = 0; k < nn + 0; ++k)
            {
                v(k) = U(j,k);
                a(k) = d(k)*v(k);
                sigma += v(k)*a(k);
            }

            U(j,j) = sigma;
            if(j == 0 || sigma == T(0.0)) continue;
            dinv = T(1.0)/sigma;

            for(k = 0; k < j; ++k)
            {
                sigma = T(0.0);

                for(i = 0; i < nn + 0; ++i)
                    sigma += U(k,i)*a(i);

                sigma *= dinv;

                for(i = 0; i < nn + 0; ++i)
                    U(k,i) -= sigma*v(i);

                U(j,k) = sigma;
            }
        }

        for(j = 0 + 1; j < n + 0; ++j)
            for(i = 0; i < j; ++i)
                U(i,j) = U(j,i);
    }


    template<typename T>
    void EFilter<T>::measureUpdate(T dz, T)
    {
        K_UINT_32 i, j, k;
        T alpha, gamma, beta, lambda;

        for (j = 0; j < n + 0; ++j)
            dz -= a(j)*_x(j);

        for(j = n - 1 + 0; j > 0; --j)
        {
            for(k = 0; k < j; ++k)
                a(j) += U(k,j)*a(k);
            d(j) = U(j,j)*a(j);
        }

        d(0) = U(0,0)*a(0);
        alpha = r+d(0)*a(0);
        gamma = T(1.0)/alpha;
        U(0,0) = r*gamma*U(0,0);

        for(j = 0 + 1; j < n; ++j)
        {
            beta = alpha;
            alpha += d(j)*a(j);
            lambda = -a(j)*gamma;
            gamma = T(1.0)/alpha;
            U(j,j) *= beta*gamma;

            for(i = 0; i < j; ++i)
            {
                beta = U(i,j);
                U(i,j) = beta+d(i)*lambda;
                d(i) += d(j)*beta;
            }
        }

        dz *= gamma;
        for(j = 0; j < n + 0; ++j)
            _x(j) += d(j)*dz;
    }

    template<typename T>
    void EFilter<T>::makeBaseAImpl()
    {
        modified_ = true;
        makeBaseA();
        if (modified_) flags |= KALMAN_A_MODIFIED;
    }


    template<typename T>
    void EFilter<T>::makeBaseWImpl()
    {
        modified_ = true;
        makeBaseW();
        if (modified_) flags |= KALMAN_W_MODIFIED;
    }


    template<typename T>
    void EFilter<T>::makeBaseQImpl()
    {
        modified_ = true;
        makeBaseQ();
        if (modified_) flags |= KALMAN_Q_MODIFIED;
    }


    template<typename T>
      void EFilter<T>::makeBaseHImpl() {
        modified_ = true;
        makeBaseH();
        if (modified_)
          flags |= KALMAN_H_MODIFIED;
      }

      template<typename T>
      void EFilter<T>::makeBaseVImpl() {
        modified_ = true;
        makeBaseV();
        if (modified_)
          flags |= KALMAN_V_MODIFIED;
      }

      template<typename T>
      void EFilter<T>::makeBaseRImpl() {
        modified_ = true;
        makeBaseR();
        if (modified_)
          flags |= KALMAN_R_MODIFIED;
      }

      template<typename T>
      void EFilter<T>::makeAImpl()
      {
        modified_ = true;
        makeA();
        if (modified_)
          flags |= KALMAN_A_MODIFIED;
      }

    template<typename T>
    void EFilter<T>::makeWImpl()
    {
        modified_ = true;
        makeW();
        if (modified_) flags |= KALMAN_W_MODIFIED;
    }

    template<typename T>
    void EFilter<T>::makeQImpl()
    {
        modified_ = true;
        makeQ();
        if (modified_) flags |= KALMAN_Q_MODIFIED;
    }

    template<typename T>
    void EFilter<T>::makeHImpl()
    {
        modified_ = true;
        makeH();
        if (modified_) flags |= KALMAN_H_MODIFIED;
    }

    template<typename T>
    void EFilter<T>::makeVImpl()
    {
        modified_ = true;
        makeV();
        if (modified_) flags |= KALMAN_V_MODIFIED;
    }


    template<typename T>
    void EFilter<T>::makeRImpl()
    {
        modified_ = true;
        makeR();
        if (modified_) flags |= KALMAN_R_MODIFIED;
    }
}

#endif
