/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Vector3.inl
/// \brief   3维向量的实现
/// \author  DaiMingZhuang
/// \version 
/// \date    2009-07
/////////////////////////////////////////////////////////////////////////////////
/// 零向量，全部数值为0
template <typename Real>
const Vector3<Real> Vector3<Real>::ZERO(0.0f);
/// 一向量，全部数值为1
template <typename Real>
const Vector3<Real> Vector3<Real>::ONE(1.0f);
/// X轴坐标为1，Y、Z轴坐标为0
template <typename Real>
const Vector3<Real> Vector3<Real>::UNIT_X(1.0f, 0.0f, 0.0f);
/// Y轴坐标为1，X、Z轴坐标为0
template <typename Real>
const Vector3<Real> Vector3<Real>::UNIT_Y(0.0f, 1.0f, 0.0f);
/// Z轴坐标为1，X、Y轴坐标为0
template <typename Real>
const Vector3<Real> Vector3<Real>::UNIT_Z(0.0f, 0.0f, 1.0f);

//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real>::Vector3(void)
{
    //memcpy(THISP(m_afTuple), ZERO.m_afTuple, sizeof(THISP(m_afTuple)) );
    //留空以加速
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real>::Vector3(const Vector3& v)
{
    memcpy(this, &v, 3 * sizeof(Real));
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real>::Vector3(Real xyz)
{
    x = y = z = xyz;
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real>::Vector3(Real _x, Real _y, Real _z)
{
    x = _x;
    y = _y;
    z = _z;
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real>::Vector3(const Vector<3,Real>& v)
{
    memcpy(this, &v, 3 * sizeof(Real));
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real>& Vector3<Real>::operator= (const Vector<3,Real>& v)
{
    memcpy(this, &v, 3 * sizeof(Real));
    return(*this);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real>& Vector3<Real>::operator= (const Vector3<Real>& v)
{
    memcpy(this, &v, 3 * sizeof(Real));
    return(*this);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Real Vector3<Real>::Dot(const Vector3<Real>& v) const
{
    return x*v.x + y*v.y + z*v.z;
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real> Vector3<Real>::Cross(const Vector3<Real>& v) const
{
    return Vector3(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Matrix3<Real> Vector3<Real>::Tensor( const Vector3<Real>& v ) const
{
	return Matrix3f(
		x * v.x, x * v.y, x * v.z,
		y * v.x, y * v.y, y * v.z,
		z * v.x, z * v.y, z * v.z
		);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real> Vector3<Real>::TransformCoord(const Matrix4<Real>& mat)const
{
    return (*this) * mat;
}
//----------------------------------------------------------------------------
template <typename Real>
Vector3<Real> Vector3<Real>::TransformNormal(const Matrix4<Real>& mat)const
{
    Vector3<Real> v3;
    v3.x = x * mat.m[0][0] + y * mat.m[1][0] + z * mat.m[2][0];
    v3.y = x * mat.m[0][1] + y * mat.m[1][1] + z * mat.m[2][1];
    v3.z = x * mat.m[0][2] + y * mat.m[1][2] + z * mat.m[2][2];
    return v3;
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector4<Real> Vector3<Real>::ToVector4(Real w )const
{
    return Vector4<Real>(x, y, z, w);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector2<Real> Vector3<Real>::ToVector2() const
{
    return Vector2<Real>(x, y);
}
//----------------------------------------------------------------------------
template <typename Real>
Vector3<Real> Normalize (const Vector3<Real>& rkV)
{
    Real fLength = rkV.Length();
    Vector3<Real> kV;

    if ( fLength > (Real)0.0 )
    {
        Real fInvLength = ((Real)1.0)/fLength;
        
        kV.x = rkV.x * fInvLength;
        kV.y = rkV.y * fInvLength;
        kV.z = rkV.z * fInvLength;
    }
    else
    {
        K_MATH_LOG("ZERO Vector cannot be normalized.");
        kV.x = kV.y = kV.z = (Real)0.0;
    }

    return kV;
}
//----------------------------------------------------------------------------

//算数运算
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real> Vector3<Real>::operator+ (const Vector3<Real>& rkV) const
{
    return Vector3<Real>(x + rkV.x, y + rkV.y, z + rkV.z);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real> Vector3<Real>::operator- (const Vector3<Real>& rkV) const
{
    return Vector3<Real>(x - rkV.x, y - rkV.y, z - rkV.z);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real> Vector3<Real>::operator* (Real fScalar) const
{
    return Vector3<Real>(fScalar * x, fScalar * y, fScalar * z);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real> operator* (Real fScalar, const Vector3<Real>& rkV)
{
    return Vector3<Real>(fScalar * rkV.x, fScalar * rkV.y, fScalar * rkV.z);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real> Vector3<Real>::operator/ (Real fScalar) const
{
    Real fInvScalar = ((Real)1.0)/fScalar;
    return Vector3<Real>(fInvScalar* x, fInvScalar* y, fInvScalar* z);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real> Vector3<Real>::operator+ () const
{
	return *this;
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real> Vector3<Real>::operator- () const
{
    return Vector3<Real>(-x, -y, -z);
}
//----------------------------------------------------------------------------
template <typename Real>
inline Vector3<Real> Cross(const Vector3<Real>& v1, const Vector3<Real>& v2)
{
    return v1.Cross(v2);
}
//----------------------------------------------------------------------------
template <typename Real>
Matrix3<Real> Tensor(const Vector3<Real>& u, const Vector3<Real>& v)
{
	return u.Tensor(v);
}
//----------------------------------------------------------------------------
