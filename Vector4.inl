/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Vector4.inl
/// \brief   4维向量的实现
/// \author  DaiMingZhuang
/// \version 
/// \date    2009-07
/////////////////////////////////////////////////////////////////////////////////
/// 零向量，全部数值为0
template <typename Real>
const Vector4<Real> Vector4<Real>::ZERO(0.0f);
/// 一向量，全部数值为1
template <typename Real>
const Vector4<Real> Vector4<Real>::ONE(1.0f);
/// X轴坐标为1，Y、Z、W轴坐标为0
template <typename Real>
const Vector4<Real> Vector4<Real>::UNIT_X(1.0f, 0.0f, 0.0f, 0.0f);
/// Y轴坐标为1，X、Z、W轴坐标为0
template <typename Real>
const Vector4<Real> Vector4<Real>::UNIT_Y(0.0f, 1.0f, 0.0f, 0.0f);
/// Z轴坐标为1，X、Y、W轴坐标为0
template <typename Real>
const Vector4<Real> Vector4<Real>::UNIT_Z(0.0f, 0.0f, 1.0f, 0.0f);
/// W轴坐标为1，X、Y、Z轴坐标为0
template <typename Real>
const Vector4<Real> Vector4<Real>::UNIT_W(0.0f, 0.0f, 0.0f, 1.0f);

//----------------------------------------------------------------------------
template <typename Real>
inline Vector4<Real>::Vector4(void)
{
    //memcpy(this, ZERO.m_afTuple, sizeof(THISP(m_afTuple)) );
    //留空以加速
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector4<Real>::Vector4(const Vector4& v)
{
    memcpy(this, &v, 4 * sizeof(Real));
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector4<Real>::Vector4(Real xyzw)
{
    x = y = z = w = xyzw;
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector4<Real>::Vector4(Real _x, Real _y, Real _z, Real _w)
{
    x = _x;
    y = _y;
    z = _z;
    w = _w;
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector4<Real>::Vector4(const Vector<4,Real>& v)
{
    memcpy(this, &v, 4 * sizeof(Real));
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector4<Real>& Vector4<Real>::operator= (const Vector<4,Real>& v)
{
    memcpy(this, &v, 4 * sizeof(Real));
    return(*this);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector4<Real>& Vector4<Real>::operator= (const Vector4& v)
{
    memcpy(this, &v, 4 * sizeof(Real));
    return(*this);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real> Vector4<Real>::ToVector3() const
{
    return Vector3<Real>(x, y, z);
}
//----------------------------------------------------------------------------
template <typename Real>
inline void Vector4<Real>::Homogenize()
{
    if(w != (Real)0.0)
    {
        Real inv = (Real)1.0 /w;
        x *= inv;
        y *= inv;
        z *= inv;
        w = (Real)1.0;
    }
}
//----------------------------------------------------------------------------
template <typename Real>
Vector4<Real> Normalize (const Vector4<Real>& rkV)
{
    Real fLength = rkV.Length();
    Vector4<Real> kV;

    if ( fLength > (Real)0.0 )
    {
        Real fInvLength = ((Real)1.0)/fLength;
        for (int i = 0; i < 4; i++)
            kV[i] = rkV[i] * fInvLength;
    }
    else
    {
        K_MATH_LOG("ZERO Vector cannot be normalized.");
        for (int i = 0; i < 4; i++)
            kV[i] = (Real)0.0;
    }

    return kV;
}
//----------------------------------------------------------------------------
template <typename Real>
Vector4<Real> Homogenize( const Vector4<Real>& rkV )
{
    Vector4<Real> res(rkV);
    res.Homogenize();
    return res;
}
//----------------------------------------------------------------------------

//算数运算
//----------------------------------------------------------------------------
template <typename Real>
Vector4<Real> Vector4<Real>::operator+ (const Vector4<Real>& rkV) const
{
    Vector4<Real> kSum;
    for (int i = 0; i < 4; i++)
        kSum[i] = (*this)[i] + rkV[i];
    return kSum;
}
//----------------------------------------------------------------------------
template <typename Real>
Vector4<Real> Vector4<Real>::operator- (const Vector4<Real>& rkV) const
{
    Vector4<Real> kDiff;
    for (int i = 0; i < 4; i++)
        kDiff[i] = (*this)[i] - rkV[i];
    return kDiff;
}
//----------------------------------------------------------------------------
template <typename Real>
Vector4<Real> Vector4<Real>::operator* (Real fScalar) const
{
    Vector4<Real> kProd;
    for (int i = 0; i < 4; i++)
        kProd[i] = fScalar * (*this)[i];
    return kProd;
}
//----------------------------------------------------------------------------
template <typename Real>
Vector4<Real> operator* (Real fScalar, const Vector4<Real>& rkV)
{
    Vector4<Real> kProd;
    for (int i = 0; i < 4; i++)
        kProd[i] = fScalar * rkV[i];
    return kProd;
}
//----------------------------------------------------------------------------
template <typename Real>
Vector4<Real> Vector4<Real>::operator/ (Real fScalar) const
{
    Vector4<Real> kQuot;
    Real fInvScalar = ((Real)1.0)/fScalar;
    for (int i = 0; i < 4; i++)
        kQuot[i] = fInvScalar* (*this)[i];
    return kQuot;
}
//----------------------------------------------------------------------------
template <typename Real>
Vector4<Real> Vector4<Real>::operator+ () const
{
	return *this;	
}
//----------------------------------------------------------------------------
template <typename Real>
Vector4<Real> Vector4<Real>::operator- () const
{
    Vector4<Real> kNeg;
    for (int i = 0; i < 4; i++)
        kNeg[i] = -(*this)[i];
    return kNeg;
}
//----------------------------------------------------------------------------
