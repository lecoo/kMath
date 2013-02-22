/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Vector2.h
/// \brief   2维向量
/// \author  DaiMingzhuang
/// \version 
/// \date    2009-02
/////////////////////////////////////////////////////////////////////////////////
#ifndef __VECTOR2_H__
#define __VECTOR2_H__

#include "Vector.h"

NS_KMATH_BEGIN

/////////////////////////////////////////////////////////////////////////////////
/// \class Vector2
/// \brief 二维向量(顶点)
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
class Vector2 : public Vector<2,Real>
{
public:
	/// \brief 构造函数，将x，y初始化为0。
	Vector2();
    /// \brief 拷贝构造函数
    Vector2(const Vector2& v);
	/// \brief 构造函数，设置x，y为xy。
	explicit Vector2(Real xy);
	/// \brief 构造函数，分别设置x，y。
	Vector2(Real _x, Real _y);
    /// \brief 构造函数，用基类对象进行构造
    Vector2(const Vector<2,Real>& v);
    /// \brief 赋值运算，用基类对象进行赋值
    Vector2& operator= (const Vector<2,Real>& v);
    /// \brief 赋值运算
    Vector2& operator= (const Vector2& v);

    // arithmetic operations
    /// \brief 向量加法运算
    Vector2 operator+ (const Vector2& rkV) const;
    /// \brief 向量减法运算
    Vector2 operator- (const Vector2& rkV) const;
    /// \brief 向量数乘运算
    Vector2 operator* (Real fScalar) const;
    /// \brief 向量除以数(与1/fScalar的数乘运算)
    Vector2 operator/ (Real fScalar) const;
	/// \brief 向量取正运算
	Vector2 operator+ () const;
    /// \brief 向量取负运算
    Vector2 operator- () const;

	/// \brief 求二维向量的叉乘(外积)。
	Real Cross(const Vector2& v) const;

    /// \brief 转换为3维向量
    Vector3<Real> ToVector3(Real z = 1.0f)const;
    
public:
    /// X轴上的坐标
    Real x;
    /// Y轴上的坐标
    Real y;

public:
	/// 零向量，全部数值为0
	static const Vector2 ZERO;
	/// 一向量，全部数值为1
	static const Vector2 ONE;
	/// X轴坐标为1，Y轴坐标为0
	static const Vector2 UNIT_X;
	/// X轴坐标为0，Y轴坐标为1
	static const Vector2 UNIT_Y;
};

/// float类型的二维向量
typedef Vector2<float> Vector2f;

/// double类型的二维向量
typedef Vector2<double> Vector2d;

/// \brief 返回参数向量的标准化向量，不改变参数向量
template <typename Real>
Vector2<Real> Normalize (const Vector2<Real>& rkV);

/// \brief 向量数乘运算
template <typename Real>
Vector2<Real> operator* (Real fScalar, const Vector2<Real>& rkV);

/// \brief 二维向量拟叉积运算
template <typename Real>
Real PseudoCross(const Vector2<Real>& v1, const Vector2<Real>& v2);

#include "Vector2.inl"

NS_KMATH_END

#endif // __VECTOR2_H__
