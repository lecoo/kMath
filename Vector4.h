/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Vector4.h
/// \brief   4维向量
/// \author  DaiMingzhuang
/// \version 
/// \date    2009-02
/////////////////////////////////////////////////////////////////////////////////
#ifndef __VECTOR4_H__
#define __VECTOR4_H__

#include "Vector.h"

NS_KMATH_BEGIN

/////////////////////////////////////////////////////////////////////////////////
/// \class Vector4
/// \brief 四维向量(顶点)
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
class Vector4 : public Vector<4,Real>
{
public:
	/// \brief 构造函数，将x，y，z，w初始化为0。
    Vector4(void);
    /// \brief 拷贝构造函数
    Vector4(const Vector4& v);
	/// \brief 构造函数，设置x，y，z，w为xyzw。
	explicit Vector4(Real xyzw);
	/// \brief 构造函数，分别设置x，y，z，w。
	Vector4(Real _x, Real _y, Real _z, Real _w);
    /// \brief 构造函数，用基类对象进行构造
    Vector4(const Vector<4,Real>& v);
    /// \brief 赋值运算，用基类对象进行赋值
    Vector4& operator= (const Vector<4,Real>& v);
    /// \brief 赋值运算
    Vector4& operator= (const Vector4& v);

    // arithmetic operations
    /// \brief 向量加法运算
    Vector4 operator+ (const Vector4& rkV) const;
    /// \brief 向量减法运算
    Vector4 operator- (const Vector4& rkV) const;
    /// \brief 向量数乘运算
    Vector4 operator* (Real fScalar) const;
    /// \brief 向量除以数(与1/fScalar的数乘运算)
    Vector4 operator/ (Real fScalar) const;
	/// \brief 向量取正运算
	Vector4 operator+ () const;
    /// \brief 向量取负运算
    Vector4 operator- () const;

    /// \brief 转换为Vector3类型数据，丢弃w分量
    Vector3<Real> ToVector3() const;

    /// \brief 将当前向量齐次化
    void Homogenize();

public:
    /// X轴上的坐标
    Real x;
    /// Y轴上的坐标
    Real y;
    /// Z轴上的坐标
    Real z;
    /// W轴上的坐标
    Real w;

public:
	/// 零向量，全部数值为0
	static const Vector4 ZERO;
	/// 一向量，全部数值为1
	static const Vector4 ONE;
	/// X轴坐标为1，Y、Z、W轴坐标为0
	static const Vector4 UNIT_X;
	/// Y轴坐标为1，X、Z、W轴坐标为0
	static const Vector4 UNIT_Y;
	/// Z轴坐标为1，X、Y、W轴坐标为0
	static const Vector4 UNIT_Z;
	/// W轴坐标为1，X、Y、Z轴坐标为0
	static const Vector4 UNIT_W;
};

/// float类型的四维向量
typedef Vector4<float> Vector4f;

/// double类型的四维向量
typedef Vector4<double> Vector4d;

/// \brief 返回参数向量的标准化向量，不改变参数向量
template <typename Real>
Vector4<Real> Normalize (const Vector4<Real>& rkV);

/// \brief 返回参数向量的齐次化向量，不改变参数向量
template <typename Real>
Vector4<Real> Homogenize (const Vector4<Real>& rkV);

/// \brief 向量数乘运算
template <typename Real>
Vector4<Real> operator* (float fScalar, const Vector4<Real>& rkV);

#include "Vector4.inl"

NS_KMATH_END

#endif // __VECTOR4_H__
