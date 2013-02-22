/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Quaternion.h
/// \brief   四元数
/// \author  DaiMingzhuang
/// \version 
/// \date    2009-03
/////////////////////////////////////////////////////////////////////////////////
#ifndef __QUATERNION_H__
#define __QUATERNION_H__

#include "MathDef.h"

NS_KMATH_BEGIN

/////////////////////////////////////////////////////////////////////////////////
/// \class Quaternion
/// \brief 四元数
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
class Quaternion
{
public:
	/// \brief 构造函数，将所有数据初始化为0。
	Quaternion(void) {}
	/// \brief 构造函数，设置每个数据的值。
	Quaternion(Real _x, Real _y, Real _z, Real _w) { x = _x; y = _y; z = _z; w = _w; }
	/// \brief 构造函数，用p[4]初始化数据。
	explicit Quaternion(const Real* p) { x = p[0]; y = p[1]; z = p[2]; w = p[3]; }

	/// \brief 设置每个数据的值。
	Quaternion& Set(Real _x, Real _y, Real _z, Real _w) { x = _x; y = _y; z = _z; w = _w; return *this; }
	/// \brief 用p[4]初始化数据。
	Quaternion& Set(const Real* p) { x = p[0]; y = p[1]; z = p[2]; w = p[3]; return *this; }

	/// \brief 拷贝函数。
	Quaternion& operator = (const Quaternion& qtn) { x = qtn.x; y = qtn.y; z = qtn.z; w = qtn.w; return *this; }
	/// \brief 比较两个四元数是否相等。
	bool operator == (const Quaternion& qtn) const { return (w == qtn.w && x == qtn.x && y == qtn.y && z == qtn.z); }
	/// \brief 比较两个四元数是否不等。
	bool operator != (const Quaternion& qtn) const { return !(*this==qtn); }

	/// \brief 四元数加法，按对应位置相加。
	Quaternion operator + (const Quaternion& qtn) const { return Quaternion(x+qtn.x, y+qtn.y, z+qtn.z, w+qtn.w); }
	/// \brief 四元数减法，按对应位置相减。
	Quaternion operator - (const Quaternion& qtn) const { return Quaternion(x-qtn.x, y-qtn.y, z-qtn.z, w-qtn.w); }
	/// \brief 四元数乘法。
	Quaternion operator * (const Quaternion& qtn) const;
	/// \brief 四元数乘以数值，所有数据乘以s。
	Quaternion operator * (Real s) const { return Quaternion(x*s, y*s, z*s, w*s); }
	/// \brief 四元数除以数值，所有数据除以s。
	Quaternion operator / (Real s) const { s = 1.0f / s; return Quaternion(x*s, y*s, z*s, w*s); }
	/// \brief 四元数求反，所有数据取相反值。
	Quaternion operator - (void) const { return Quaternion(-x, -y, -z, -w); }

	/// \brief 四元数加法，按对应位置相加，结果赋给自己。
	Quaternion& operator += (const Quaternion& qtn) { x += qtn.x; y += qtn.y; z += qtn.z; w += qtn.w; return *this; }
	/// \brief 四元数减法，按对应位置相减，结果赋给自己。
	Quaternion& operator -= (const Quaternion& qtn) { x -= qtn.x; y -= qtn.y; z -= qtn.z; w -= qtn.w; return *this; }
	/// \brief 四元数乘法，结果赋给自己。
	Quaternion& operator *= (const Quaternion& qtn) { return *this = Mul(qtn); }
	/// \brief 四元数乘以数值，所有数据乘以s，结果赋给自己。
	Quaternion& operator *= (Real s) { x *= s; y *= s; z *= s; w *= s; return *this; }
	/// \brief 四元数除以数值，所有数据除以s，结果赋给自己。
	Quaternion& operator /= (Real s) { return (*this) *= (1.0f / s); }

	/// \brief 设置单位四元数。
	Quaternion& SetIdentity(void);

    /// \brief 求四元数的共轭四元数。
    Quaternion Conjugate() const;
	/// \brief 求四元数的逆。
	Quaternion Inverse(void) const;

	Quaternion<Real> Slerp(const Quaternion<Real>& qtn2, Real lambda);

	/// \brief 将旋转矩阵转换成四元数。
	void FromRotationMatrix(const Matrix3<Real>& mtx);
	/// \brief 将四元数转换成旋转矩阵。
	void ToRotationMatrix(Matrix3<Real>& mtx) const;
	/// \brief 将旋转轴和转角转换成四元数。
	void FromAxisAngle(const Vector3<Real>& axis, Real ang);
	/// \brief 将四元数转换成旋转轴和转角。
	bool ToAxisAngle(Vector3<Real>& axis, Real& ang) const;

	Vector3<Real> GetAxis() const;
	Real GetAngle() const;

	/// \brief 将欧拉角转换成四元数，按ZXY的次序。
	void FromEulerAnglesZXY(const Vector3<Real>& angles);

	/// \brief 将四元数转换成欧拉角，按ZXY的次序。
	Vector3<Real> GetEulerZXY() const;

	Real Dot(const Quaternion& qtn) const;

	Real Length() const;

public:
	// 相当于axis.x * sin(ang/2)
	Real x;
	// 相当于axis.y * sin(ang/2)
	Real y;
	// 相当于axis.z * sin(ang/2)
	Real z;
	// 相当于cos(ang/2)
	Real w;

public:
	/// 单位四元数，w为1，x、y、z为0
	static const Quaternion IDENTITY;
};
/// float类型的四元数
typedef Quaternion<float> Quaternionf;

/// double类型的四元数
typedef Quaternion<double> Quaterniond;

/// \brief 球面线性插值
template<typename Real>
Quaternion<Real> Slerp(const Quaternion<Real>& qtn1, const Quaternion<Real>& qtn2, Real lambda);

template<typename Real>
Real Dot(const Quaternion<Real>& q1, const Quaternion<Real>& q2);

template<typename Real>
Real Length(const Quaternion<Real>& q);

template<typename Real>
Quaternion<Real> Normalize(const Quaternion<Real>& q);

#include "Quaternion.inl"

NS_KMATH_END

#endif // __QUATERNION_H__
