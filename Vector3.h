/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Vector3.h
/// \brief   3维向量
/// \author  DaiMingzhuang
/// \version 
/// \date    2009-02
/////////////////////////////////////////////////////////////////////////////////
#ifndef __VECTOR3_H__
#define __VECTOR3_H__

#include "Vector.h"

NS_KMATH_BEGIN

/////////////////////////////////////////////////////////////////////////////////
/// \class Vector3
/// \brief 三维向量(顶点)
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
class Vector3 : public Vector<3,Real>
{
public:
	/// \brief 构造函数，将x，y，z初始化为0。
    Vector3(void);
    /// \brief 拷贝构造函数
    Vector3(const Vector3& v);
	/// \brief 构造函数，设置x，y，z为xyz。
	explicit Vector3(Real xyz);
	/// \brief 构造函数，分别设置x，y，z。
	Vector3(Real _x, Real _y, Real _z);
    /// \brief 构造函数，用基类对象进行构造
    Vector3(const Vector<3,Real>& v);
    /// \brief 赋值运算，用基类对象进行赋值
    Vector3& operator= (const Vector<3,Real>& v);
    /// \brief 赋值运算
    Vector3& operator= (const Vector3& v);

    // arithmetic operations
    /// \brief 向量加法运算
    Vector3 operator+ (const Vector3& rkV) const;
    /// \brief 向量减法运算
    Vector3 operator- (const Vector3& rkV) const;
    /// \brief 向量数乘运算
    Vector3 operator* (Real fScalar) const;
    /// \brief 向量除以数(与1/fScalar的数乘运算)
    Vector3 operator/ (Real fScalar) const;
	/// \brief 向量取正运算
	Vector3 operator+ () const;
    /// \brief 向量取负运算
    Vector3 operator- () const;

    /// \brief 数量积
    Real Dot (const Vector3<Real>& v) const;
	/// \brief 向量积
	Vector3<Real> Cross(const Vector3<Real>& v) const;
	/// \brief 张量积
	Matrix3<Real> Tensor(const Vector3<Real>& v) const;

    /// \brief 视向量为(x,y,z,1)用给定矩阵进行变换，并对结果齐次化
    Vector3<Real> TransformCoord(const Matrix4<Real>& mat)const;

    /// \brief 视向量为(x,y,z,0)用给定矩阵进行变换
    /// 如果想做法线变换，给定的矩阵应该是对顶点进行变换的矩阵的逆矩阵的转置矩阵
    Vector3<Real> TransformNormal(const Matrix4<Real>& mat)const;

    /// \brief 转换为4维向量
    Vector4<Real> ToVector4(Real w = 1.0f)const;

    /// \brief 转换为2维向量，丢弃z分量
    Vector2<Real> ToVector2() const;

public:
    /// X轴上的坐标
    Real x;
    /// Y轴上的坐标
    Real y;
    /// Z轴上的坐标
    Real z;

public:
	/// 零向量，全部数值为0
	static const Vector3 ZERO;
	/// 一向量，全部数值为1
	static const Vector3 ONE;
	/// X轴坐标为1，Y、Z轴坐标为0
	static const Vector3 UNIT_X;
	/// Y轴坐标为1，X、Z轴坐标为0
	static const Vector3 UNIT_Y;
	/// Z轴坐标为1，X、Y轴坐标为0
	static const Vector3 UNIT_Z;
};

/// float类型的三维向量
typedef Vector3<float> Vector3f;

/// double类型的三维向量
typedef Vector3<double> Vector3d;

/// \brief 返回参数向量的标准化向量，不改变参数向量
template <typename Real>
Vector3<Real> Normalize (const Vector3<Real>& rkV);

/// \brief 向量数乘运算
template <typename Real>
Vector3<Real> operator* (Real fScalar, const Vector3<Real>& rkV);

/// \brief 3维向量叉乘运算
template <typename Real>
Vector3<Real> Cross(const Vector3<Real>& v1, const Vector3<Real>& v2);

/// \brief 3维向量张量积
template <typename Real>
Matrix3<Real> Tensor(const Vector3<Real>& u, const Vector3<Real>& v);

#include "Vector3.inl"

#ifdef USE_DX_MATH
#include "Vector3f_DX.inl"  //用DX特化Vector3<float>的实现
#endif

NS_KMATH_END

#endif // __VECTOR3_H__
