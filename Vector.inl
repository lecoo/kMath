/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Vector.inl
/// \brief   N维向量的实现
/// \author  DaiMingZhuang
/// \version 
/// \date    2009-07
/////////////////////////////////////////////////////////////////////////////////
template <int N, class Real>
inline Vector<N,Real>::Vector()
{
    //留空以加速
}
//----------------------------------------------------------------------------
template <int N, class Real>
inline Vector<N,Real>::Vector (const Vector& rkV)
{
    memcpy(this, &rkV, N*sizeof(Real));
}
//----------------------------------------------------------------------------
template <int N, class Real>
inline Vector<N,Real>::Vector (const Real* p)
{
    memcpy(this, p, N*sizeof(Real));
}
//----------------------------------------------------------------------------
template <int N, class Real>
inline Real Vector<N,Real>::operator[] (int i) const
{
    return ((Real*)this)[i];
}
//----------------------------------------------------------------------------
template <int N, class Real>
inline Real& Vector<N,Real>::operator[] (int i)
{
    return ((Real*)this)[i];
}
//----------------------------------------------------------------------------
template <int N, class Real>
inline Vector<N,Real>& Vector<N,Real>::operator= (const Vector& rkV)
{
    memcpy(this,&rkV,N*sizeof(Real));
    return *this;
}
//----------------------------------------------------------------------------
template <int N, class Real>
inline bool Vector<N,Real>::operator== (const Vector& rkV) const
{
    return memcmp(this, &rkV, N*sizeof(Real)) == 0;
}
//----------------------------------------------------------------------------
template <int N, class Real>
inline bool Vector<N,Real>::operator!= (const Vector& rkV) const
{
    return memcmp(this, &rkV, N*sizeof(Real)) != 0;
}
//----------------------------------------------------------------------------
template <int N, class Real>
inline bool Vector<N,Real>::operator< (const Vector& rkV) const
{
    return memcmp(this, &rkV, N*sizeof(Real)) < 0;
}
//----------------------------------------------------------------------------
template <int N, class Real>
Vector<N,Real>& Vector<N,Real>::operator+= (const Vector& rkV)
{
    for (int i = 0; i < N; i++)
        (*this)[i] += rkV[i];
    return *this;
}
//----------------------------------------------------------------------------
template <int N, class Real>
Vector<N,Real>& Vector<N,Real>::operator-= (const Vector& rkV)
{
    for (int i = 0; i < N; i++)
        (*this)[i] -= rkV[i];
    return *this;
}
//----------------------------------------------------------------------------
template <int N, class Real>
Vector<N,Real>& Vector<N,Real>::operator*= (Real fScalar)
{
    for (int i = 0; i < N; i++)
        (*this)[i] *= fScalar;
    return *this;
}
//----------------------------------------------------------------------------
template <int N, class Real>
Vector<N,Real>& Vector<N,Real>::operator/= (Real fScalar)
{
    Real fInvScalar = ((Real)1.0)/fScalar;
    for (int i = 0; i < N; i++)
        (*this)[i] *= fInvScalar;

    return *this;
}
//----------------------------------------------------------------------------
template <int N, class Real>
Real Vector<N,Real>::Length () const
{
    Real fSqrLen = (Real)0.0;
    for (int i = 0; i < N; i++)
        fSqrLen += MathUtil::Square((*this)[i]);
    return MathUtil::Sqrt(fSqrLen);
}
//----------------------------------------------------------------------------
template <int N, class Real>
Real Vector<N,Real>::LengthSquared () const
{
    Real fSqrLen = (Real)0.0;
    for (int i = 0; i < N; i++)
        fSqrLen += MathUtil::Square((*this)[i]);
    return fSqrLen;
}
//----------------------------------------------------------------------------
template <int N, class Real>
Real Vector<N,Real>::Dot (const Vector& rkV) const
{
    Real fDot = (Real)0.0;
    for (int i = 0; i < N; i++)
        fDot += (*this)[i]*rkV[i];
    return fDot;
}
//----------------------------------------------------------------------------
template <int N, class Real>
void Vector<N,Real>::Normalize ()
{
    Real fLength = Length();

    if ( fLength > (Real)0.0 )
    {
        Real fInvLength = ((Real)1.0)/fLength;
        for (int i = 0; i < N; i++)
            (*this)[i] *= fInvLength;
    }
    else
    {
        K_MATH_LOG("ZERO Vector cannot be normalized.");
    }
}
//----------------------------------------------------------------------------
template <int N, class Real>
inline Real Dot(const Vector<N,Real>& v1, const Vector<N,Real>& v2)
{
    return v1.Dot(v2);
}
//----------------------------------------------------------------------------
template <int N, class Real>
Real Length(const Vector<N,Real>& v)
{
	return v.Length();
}
//----------------------------------------------------------------------------
