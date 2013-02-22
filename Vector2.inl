/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Vector2.inl
/// \brief   Nά����,�����Ļ���
/// \author  DaiMingZhuang
/// \version 0.9
/// \date    2009-07
/// \history 
/////////////////////////////////////////////////////////////////////////////////
/// ��������ȫ����ֵΪ0
template <typename Real>
const Vector2<Real> Vector2<Real>::ZERO(0.0f);
/// һ������ȫ����ֵΪ1
template <typename Real>
const Vector2<Real> Vector2<Real>::ONE(1.0f);
/// X������Ϊ1��Y������Ϊ0
template <typename Real>
const Vector2<Real> Vector2<Real>::UNIT_X(1.0f, 0.0f);
/// X������Ϊ0��Y������Ϊ1
template <typename Real>
const Vector2<Real> Vector2<Real>::UNIT_Y(0.0f, 1.0f);

//----------------------------------------------------------------------------
template <typename Real>
inline Vector2<Real>::Vector2()
{
    //�����Լ���
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector2<Real>::Vector2(const Vector2& v)
{
    memcpy(this, &v, 2*sizeof(Real));
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector2<Real>::Vector2(Real xy)
{
    x = y = xy;
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector2<Real>::Vector2(Real _x, Real _y)
{
    x = _x;
    y = _y;
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector2<Real>::Vector2(const Vector<2,Real>& v)
{
    memcpy(this, &v , 2 * sizeof(Real));
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector2<Real>& Vector2<Real>::operator= (const Vector<2,Real>& v)
{
    memcpy(this, &v , 2 * sizeof(Real));
    return(*this);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector2<Real>& Vector2<Real>::operator= (const Vector2& v)
{
    memcpy(this, &v , 2 * sizeof(Real));
    return(*this);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Real Vector2<Real>::Cross(const Vector2& v) const
{
    return (x * v.y - y * v.x);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real> Vector2<Real>::ToVector3(Real z )const
{
    return Vector3<Real>(x, y, z);
}
//----------------------------------------------------------------------------

template <typename Real>
Vector2<Real> Normalize (const Vector2<Real>& rkV)
{
    Real fLength = rkV.Length();
    Vector2<Real> kV;

    if ( fLength > (Real)0.0 )
    {
        Real fInvLength = ((Real)1.0)/fLength;
        
        kV.x = rkV.x * fInvLength;
        kV.y = rkV.y * fInvLength;
    }
    else
    {
        K_MATH_LOG("ZERO Vector cannot be normalized.");
        kV.x = kV.y = (Real)0.0;
    }

    return kV;
}
//----------------------------------------------------------------------------

//��������
//----------------------------------------------------------------------------
template <typename Real>
inline Vector2<Real> Vector2<Real>::operator+ (const Vector2<Real>& rkV) const
{
    return Vector2<Real>(x + rkV.x, y + rkV.y);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector2<Real> Vector2<Real>::operator- (const Vector2<Real>& rkV) const
{
    return Vector2<Real>(x - rkV.x, y - rkV.y);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector2<Real> Vector2<Real>::operator* (Real fScalar) const
{
    return Vector2<Real>(fScalar * x, fScalar * y);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector2<Real> operator* (Real fScalar, const Vector2<Real>& rkV)
{
    return Vector2<Real>(fScalar * rkV.x, fScalar * rkV.y);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector2<Real> Vector2<Real>::operator/ (Real fScalar) const
{
    Real fInvScalar = ((Real)1.0)/fScalar;
    return Vector2<Real>(fInvScalar* x, fInvScalar* y);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector2<Real> Vector2<Real>::operator+ () const
{
	return *this;
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector2<Real> Vector2<Real>::operator- () const
{
    return Vector2<Real>(-x, -y);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Real PseudoCross(const Vector2<Real>& v1, const Vector2<Real>& v2)
{
    // (-v1.y, v1.x)��(v2.x, v2.y)
    return v1.x * v2.y - v1.y * v2.x;
}
//----------------------------------------------------------------------------
