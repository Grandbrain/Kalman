#ifndef EXTENDED_H
#define EXTENDED_H

#include "vector.h"
#include "matrix.h"

namespace Kalman
{
    enum State
    {
        ModifiedN   = 1,
        ModifiedNU  = (1 << 1),
        ModifiedNW  = (1 << 2),
        ModifiedM   = (1 << 3),
        ModifiedNV  = (1 << 4),
        ModifiedP   = (1 << 5),
        Lowmask     = ((1 << 8) - 1),
        ModifiedA   = (1 << 8),
        ModifiedW   = (1 << 9),
        ModifiedQ   = (1 << 10),
        Midmask     = (((1 << 4) - 1) << 8),
        ModifiedH   = (1 << 12),
        ModifiedV   = (1 << 13),
        ModifiedR   = (1 << 14),
        Highmask    = (((1 << 4) - 1) << 12)
    };

    template<typename T>
    class Extended
    {
    public:

        explicit Extended();
        explicit Extended(unsigned, unsigned, unsigned, unsigned, unsigned);
        virtual ~Extended();

        unsigned getSizeX() const;
        unsigned getSizeU() const;
        unsigned getSizeW() const;
        unsigned getSizeZ() const;
        unsigned getSizeV() const;

        void setSizeX(unsigned);
        void setSizeU(unsigned);
        void setSizeW(unsigned);
        void setSizeZ(unsigned);
        void setSizeV(unsigned);

        void init(Vector<T>&, Matrix<T>&);
        void step(Vector<T>&, const Vector<T>&);
        void timeUpdateStep(Vector<T>&);
        void measureUpdateStep(const Vector<T>&);
        const Vector<T>& predict(Vector<T>&);
        const Vector<T>& simulate();
        const Vector<T>& getX() const;
        const Matrix<T>& calculateP() const;

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


        Vector<T> x;
        Vector<T> u;
        Vector<T> z;
        Vector<T> dz;

        Matrix<T> A;
        Matrix<T> W;
        Matrix<T> Q;
        Matrix<T> H;
        Matrix<T> V;
        Matrix<T> R;


        unsigned n;
        unsigned nu;
        unsigned nw;
        unsigned  m;
        unsigned nv;

    private:


        static void factor(Matrix<T>&);
        static void upperInvert(Matrix<T>&);
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


        Matrix<T> U;
        Matrix<T> W_;
        Matrix<T> Q_;
        Matrix<T> H_;
        Matrix<T> R_;

        Vector<T> a;
        Vector<T> d;
        Vector<T> v;

        unsigned nn;

        mutable Matrix<T> _P;
        mutable Vector<T> _x;

        unsigned flags;
        bool modified_;
    };


    template<typename T>
    Extended<T>::Extended() : flags(0)
    {
    }


    template<typename T>
    Extended<T>::Extended(unsigned n_, unsigned nu_, unsigned nw_, unsigned m_, unsigned nv_) : flags(0)
    {
        setSizeX(n_);
        setSizeU(nu_);
        setSizeW(nw_);
        setSizeZ(m_);
        setSizeV(nv_);
    }


    template<typename T>
    Extended<T>::~Extended()
    {
    }


    template<typename T>
    unsigned Extended<T>::getSizeX() const
    {
        return n;
    }


    template<typename T>
    unsigned Extended<T>::getSizeU() const
    {
        return nu;
    }


    template<typename T>
    unsigned Extended<T>::getSizeW() const
    {
        return nw;
    }


    template<typename T>
    unsigned Extended<T>::getSizeZ() const
    {
        return m;
    }


    template<typename T>
    unsigned Extended<T>::getSizeV() const
    {
        return nv;
    }


    template<typename T>
    void Extended<T>::setSizeX(unsigned n_)
    {
        if (n_ == n) return;
        flags |= ModifiedN;
        n = n_;
    }


    template<typename T>
    void Extended<T>::setSizeU(unsigned nu_)
    {
        if (nu_ == nu) return;
        flags |= ModifiedNU;
        nu = nu_;
    }


    template<typename T>
    void Extended<T>::setSizeW(unsigned nw_)
    {
        if (nw_ == nw) return;
        flags |= ModifiedNW;
        nw = nw_;
    }


    template<typename T>
    void Extended<T>::setSizeZ(unsigned m_)
    {
        if (m_ == m) return;
        flags |= ModifiedM;
        m = m_;
    }


    template<typename T>
    void Extended<T>::setSizeV(unsigned nv_)
    {
        if (nv_ == nv)
            flags |= ModifiedNV;
        nv = nv_;
    }


    template<typename T>
    void Extended<T>::init(Vector<T>& x_, Matrix<T>& P_)
    {
        x.swap(x_);
        _P.swap(P_);
        flags |= ModifiedP;
    }


    template<typename T>
    void Extended<T>::step(Vector<T>& u_, const Vector<T>& z_)
    {
        timeUpdateStep(u_);
        measureUpdateStep(z_);
    }


    template<typename T>
    void Extended<T>::timeUpdateStep(Vector<T>& u_)
    {
        unsigned i, j, k;

        sizeUpdate();
        u.swap(u_);

        makeCommonProcess();
        makeAImpl();
        makeWImpl();
        makeQImpl();
        makeProcess();

        if (flags & ModifiedQ)
        {
            Q_ = Q;
            factor(Q_);
            upperInvert(Q_);
        }

        Q.swap(Q_);

        if (flags & (ModifiedW | ModifiedQ))
        {
            for (i = 0; i < n; ++i)
            {
                for (j = 0; j < nw; ++j)
                {
                    W_(i, j) = W(i, j);
                    for (k = 0; k < j; ++k)
                        W_(i, j) += W(i, k)*Q(j, k);
                }
            }
        }

        W.swap(W_);
        timeUpdate();

        Q.swap(Q_);
        W.swap(W_);

        u.swap(u_);
        flags &= ~Midmask;
    }


    template<typename T>
    void Extended<T>::measureUpdateStep(const Vector<T>& z_)
    {
        unsigned i, j, k;
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

        if (flags & (KALMAN_V_MODIFIED | KALMAN_R_MODIFIED))
        {
            _x.resize(nv);
            for (i = 0; i < m; ++i)
            {
                for (j = 0; j < nv; ++j)
                {
                    _x(j) = T(0.0);
                    for (k = 0; k < nv; ++k)
                        _x(j) += V(i, k)*R(k, j);
                }

                for (j = 0; j < m; ++j)
                {
                    R_(i, j) = T(0.0);
                    for (k = 0; k < nv; ++k)
                        R_(i, j) += _x(k)*V(j, k);
                }
            }

            factor(R_);
            upperInvert(R_);
        }

        if (flags & (KALMAN_H_MODIFIED | KALMAN_V_MODIFIED | KALMAN_R_MODIFIED))
        {
            for (i = 0; i < m; ++i)
            {
                for (j = 0; j < n; ++j)
                {
                    H_(i, j) = H(i, j);
                    for (k = i + 1; k < m + 0; ++k)
                        H_(i, j) += R_(k, i)*H(k, j);
                }
            }
        }

        H.swap(H_);
        _x.resize(m);

        for (i = 0; i < m; ++i)
        {
            _x(i) = dz(i);
            for (k = i + 1; k < m; ++k)
                _x(i) += R_(k, i)*dz(k);
        }

        dz.swap(_x);

        _x.resize(n);
        _x = T(0.0);

        for (i = 0; i < m; ++i)
        {
            for (j = 0; j < n; ++j)
                a(j) = H(i, j);

            measureUpdate(dz(i), R_(i, i));
        }

        for (i = 0; i < n + 0; ++i)
            x(i) += _x(i);

        if (!OVR) H.swap(H_);
        flags &= ~KALMAN_HIGHMASK;
    }


    template<typename T>
    const Vector<T>& Extended<T>::predict(Vector<T>& u_)
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
    const Vector<T>& Extended<T>::simulate()
    {
        sizeUpdate();
        _x = z;

        makeCommonMeasure();
        makeMeasure();

        z.swap(_x);
        return _x;
    }


    template<typename T>
    const Vector<T>& Extended<T>::getX() const
    {
        return x;
    }


    template<typename T>
    const Matrix<T>& Extended<T>::calculateP() const
    {
        if (!(flags & ModifiedP))
        {
            _P.resize(n, n);
            for (unsigned i = 0; i < n; ++i)
            {
                _P(i, i) = U(i, i);
                for (unsigned j = i + 1; j < n; ++j)
                {
                    _P(i, j) = U(i, j)*U(j, j);
                    _P(i, i) += U(i, j)*_P(i, j);

                    for (unsigned k = j + 1; k < n; ++k)
                        _P(i, j) += U(i, k)*U(j, k)*U(k, k);
                    _P(j, i) = _P(i, j);
                }
            }
        }
        return _P;
    }


    template<typename T>
    void Extended<T>::makeBaseA()
    {
        modified_ = false;
    }


    template<typename T>
    void Extended<T>::makeBaseW()
    {
        modified_ = false;
    }


    template<typename T>
    void Extended<T>::makeBaseQ()
    {
        modified_ = false;
    }


    template<typename T>
    void Extended<T>::makeBaseH()
    {
        modified_ = false;
    }


    template<typename T>
    void Extended<T>::makeBaseV()
    {
        modified_ = false;
    }


    template<typename T>
    void Extended<T>::makeBaseR()
    {
        modified_ = false;
    }


    template<typename T>
    void Extended<T>::makeCommonProcess()
    {
    }


    template<typename T>
    void Extended<T>::makeCommonMeasure()
    {
    }


    template<typename T>
    void Extended<T>::makeA()
    {
        modified_ = false;
    }


    template<typename T>
    void Extended<T>::makeW()
    {
        modified_ = false;
    }


    template<typename T>
    void Extended<T>::makeQ()
    {
        modified_ = false;
    }


    template<typename T>
    void Extended<T>::makeH()
    {
        modified_ = false;
    }


    template<typename T>
    void Extended<T>::makeV()
    {
        modified_ = false;
    }


    template<typename T>
    void Extended<T>::makeR()
    {
        modified_ = false;
    }


    template<typename T>
    void Extended<T>::makeDZ()
    {
    }


    template<typename T>
    void Extended<T>::sizeUpdate()
    {
        if (!flags) return;
        if (flags & ModifiedN)
        {
            A.resize(n, n);
            makeBaseAImpl();
        }

        if (flags & (ModifiedN | ModifiedNW))
        {
            nn = n + nw;
            a.resize(nn);
            v.resize(nn);
            d.resize(nn);
            W_.resize(n, nw);
            W.resize(n, nw);
            makeBaseWImpl();
        }

        if (flags & ModifiedP)
        {
            U.resize(n, nn);
            for (unsigned i = 0; i < n; ++i)
                for (unsigned j = 0; j < n; ++j)
                    U(i, j) = _P(i, j);

            factor(U);

        }
        else if (flags & ModifiedNW)
        {
            _P.Resize(n, nn);
            for (unsigned i = 0; i < n; ++i)
                for (unsigned j = i; j < n; ++j)
                    _P(i, j) = U(i, j);

            U.swap(_P);
        }

        if (flags & ModifiedNW)
        {
            Q_.Resize(nw, nw);
            Q.Resize(nw, nw);
            makeBaseQImpl();
        }

        if (m != 0)
        {
            if (flags & (ModifiedN | ModifiedM))
            {
                H_.Resize(m, n);
                H.Resize(m, n);
                makeBaseHImpl();
            }

            if (flags & (ModifiedM | ModifiedNV))
            {
                V.Resize(m, nv);
                makeBaseVImpl();
            }

            if (flags & ModifiedNV) {
                R.Resize(nv, nv);
                makeBaseRImpl();
            }

            if (flags & ModifiedM) {
                R_.Resize(m, m);
                z.Resize(m);
                dz.Resize(m);
            }

        }

        flags &= ~Lowmask;
    }


    template<typename T>
    void Extended<T>::factor(Matrix<T>& P_)
    {
        T alpha, beta;
        unsigned i, j, k, N = P_.Rows();
        for (j = N - 1; j > 0; --j)
        {
            alpha = T(1.0) / P_(j, j);
            for (k = 0; k < j; ++k)
            {
                beta = P_(k, j);
                P_(k, j) = alpha*beta;
                for (i = 0; i <= k; ++i)
                    P_(i, k) -= beta*P_(i, j);
            }
        }
    }


    template<typename T>
    void Extended<T>::upperInvert(Matrix<T>& P_)
    {
        T val;
        unsigned i, j, k, N = P_.Rows();
        for (i = N - 2; i != -1; --i)
        {
            for (k = i + 1; k < N; ++k)
            {
                val = P_(i, k);
                for (j = i + 1; j <= k - 1; ++j)
                    val += P_(i, j)*P_(k, j);
                P_(k, i) = -val;
            }
        }
    }


    template<typename T>
    void Extended<T>::timeUpdate()
    {
        unsigned i, j, k;
        T sigma, dinv;

        for (j = n - 1; j > 0; --j)
        {
            for (i = 0; i <= j; ++i)
                d(i) = U(i, j);

            for (i = 0; i < n; ++i)
            {
                U(i, j) = A(i, j);
                for (k = 0; k < j; ++k)
                    U(i, j) += A(i, k)*d(k);
            }
        }

        d(0) = U(0, 0);

        for (j = 0; j < n; ++j)
            U(j, 0) = A(j, 0);

        for (i = 0; i < nw; ++i)
        {
            d(i + n) = Q(i, i);
            for (j = 0; j < n; ++j)
                U(j, i + n) = W(j, i);
        }

        for (j = n - 1; j != -1; --j)
        {
            sigma = T(0.0);
            for (k = 0; k < nn + 0; ++k)
            {
                v(k) = U(j, k);
                a(k) = d(k)*v(k);
                sigma += v(k)*a(k);
            }

            U(j, j) = sigma;
            if (j == 0 || sigma == T(0.0)) continue;
            dinv = T(1.0) / sigma;

            for (k = 0; k < j; ++k)
            {
                sigma = T(0.0);

                for (i = 0; i < nn; ++i)
                    sigma += U(k, i)*a(i);

                sigma *= dinv;

                for (i = 0; i < nn; ++i)
                    U(k, i) -= sigma*v(i);

                U(j, k) = sigma;
            }
        }

        for (j = 1; j < n; ++j)
            for (i = 0; i < j; ++i)
                U(i, j) = U(j, i);
    }


    template<typename T>
    void Extended<T>::measureUpdate(T dz, T)
    {
        unsigned i, j, k;
        T alpha, gamma, beta, lambda;

        for (j = 0; j < n; ++j)
            dz -= a(j)*_x(j);

        for (j = n - 1; j > 0; --j)
        {
            for (k = 0; k < j; ++k)
                a(j) += U(k, j)*a(k);
            d(j) = U(j, j)*a(j);
        }

        d(0) = U(0, 0)*a(0);
        alpha = r + d(0)*a(0);
        gamma = T(1.0) / alpha;
        U(0, 0) = r*gamma*U(0, 0);

        for (j = 1; j < n; ++j)
        {
            beta = alpha;
            alpha += d(j)*a(j);
            lambda = -a(j)*gamma;
            gamma = T(1.0) / alpha;
            U(j, j) *= beta*gamma;

            for (i = 0; i < j; ++i)
            {
                beta = U(i, j);
                U(i, j) = beta + d(i)*lambda;
                d(i) += d(j)*beta;
            }
        }

        dz *= gamma;
        for (j = 0; j < n; ++j)
            _x(j) += d(j)*dz;
    }


    template<typename T>
    void Extended<T>::makeBaseAImpl()
    {
        modified_ = true;
        makeBaseA();
        if (modified_) flags |= ModifiedA;
    }


    template<typename T>
    void Extended<T>::makeBaseWImpl()
    {
        modified_ = true;
        makeBaseW();
        if (modified_) flags |= ModifiedW;
    }


    template<typename T>
    void Extended<T>::makeBaseQImpl()
    {
        modified_ = true;
        makeBaseQ();
        if (modified_) flags |= ModifiedQ;
    }


    template<typename T>
    void Extended<T>::makeBaseHImpl()
    {
        modified_ = true;
        makeBaseH();
        if (modified_) flags |= ModifiedH;
    }


    template<typename T>
    void Extended<T>::makeBaseVImpl()
    {
        modified_ = true;
        makeBaseV();
        if (modified_) flags |= ModifiedV;
    }


    template<typename T>
    void Extended<T>::makeBaseRImpl()
    {
        modified_ = true;
        makeBaseR();
        if (modified_) flags |= ModifiedR;
    }


    template<typename T>
    void Extended<T>::makeAImpl()
    {
        modified_ = true;
        makeA();
        if (modified_) flags |= ModifiedA;
    }


    template<typename T>
    void Extended<T>::makeWImpl()
    {
        modified_ = true;
        makeW();
        if (modified_) flags |= ModifiedW;
    }


    template<typename T>
    void Extended<T>::makeQImpl()
    {
        modified_ = true;
        makeQ();
        if (modified_) flags |= ModifiedQ;
    }


    template<typename T>
    void Extended<T>::makeHImpl()
    {
        modified_ = true;
        makeH();
        if (modified_) flags |= ModifiedH;
    }


    template<typename T>
    void Extended<T>::makeVImpl()
    {
        modified_ = true;
        makeV();
        if (modified_) flags |= ModifiedV;
    }


    template<typename T>
    void Extended<T>::makeRImpl()
    {
        modified_ = true;
        makeR();
        if (modified_) flags |= ModifiedR;
    }
}

#endif