/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Matrix2.h
/// \brief   2×2的矩阵(右乘矩阵)
/// \author  DaiMingzhuang
/// \version 
/// \date    2009-02
/////////////////////////////////////////////////////////////////////////////////
#ifndef __MATRIX2_H__
#define __MATRIX2_H__

#include "MathDef.h"

NS_KMATH_BEGIN

#if _MSC_VER >= 1200
#pragma warning(push)
#endif
#pragma warning(disable:4201) // anonymous unions warning
/////////////////////////////////////////////////////////////////////////////////
/// \class Matrix2
/// \brief 2×2的矩阵(右乘矩阵)
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
class Matrix2
{
public:
	/// \brief 构造函数，未初始化。
	Matrix2(void);
	/// \brief 构造函数，设置每个数据的值。
	Matrix2(Real m00, Real m01, Real m10, Real m11);
	/// \brief 拷贝构造函数。
	Matrix2(const Matrix2& mtx);

	/// \brief 设置每个数据的值。
	Matrix2& Set(Real m00, Real m01, Real m10, Real m11);

	/// \brief 获取第row行，第col列的数据。
	Real operator() (int row, int col) const;
	/// \brief 获取第row行，第col列的数据，可更改。
	Real& operator() (int row, int col);

	/// \brief 获取第row行的行向量
	const Vector2<Real>& operator[] (int row) const;

	/// \brief 获取第row行的行向量
	Vector2<Real>& operator[] (int row);

	/// \brief 拷贝函数。
	Matrix2& operator = (const Matrix2& mtx);
	/// \brief 比较两个矩阵是否相等。
	bool operator == (const Matrix2& mtx) const;
	/// \brief 比较两个矩阵是否不等。
	bool operator != (const Matrix2& mtx) const;

	/// \brief 矩阵加法，按对应位置相加。
	Matrix2 operator + (const Matrix2& mtx) const;
	/// \brief 矩阵减法，按对应位置相减。
	Matrix2 operator - (const Matrix2& mtx) const;
	/// \brief 矩阵乘法。
	Matrix2 operator * (const Matrix2& mtx) const;
	/// \brief 矩阵乘以数值，所有数据乘以s。
	Matrix2 operator * (Real s) const;
	/// \brief 矩阵除以数值，所有数据除以s。
	Matrix2 operator / (Real s) const;
	/// \brief 矩阵求反，所有数据取相反值。
	Matrix2 operator - (void) const;

	/// \brief 矩阵加法，按对应位置相加，结果赋给自己。
	Matrix2& operator += (const Matrix2& mtx);
	/// \brief 矩阵减法，按对应位置相减，结果赋给自己。
	Matrix2& operator -= (const Matrix2& mtx);
	/// \brief 矩阵乘法，结果赋给自己。
	Matrix2& operator *= (const Matrix2& mtx);
	/// \brief 矩阵乘以数值，所有数据乘以s，结果赋给自己。
	Matrix2& operator *= (Real s);
	/// \brief 矩阵除以数值，所有数据除以s，结果赋给自己。
	Matrix2& operator /= (Real s);

	/// \brief 设置单位矩阵。
	Matrix2& SetIdentity(void);
	/// \brief 设置沿z轴的旋转矩阵。
	Matrix2& Set2DRotation(Real ang);

	/// \brief 求二维矩阵的转置矩阵。
	Matrix2 Transpose(void) const;
	/// \brief 求二维矩阵的逆矩阵。
	Matrix2 Inverse(void) const;
	/// \brief 求二维矩阵的行列式。
	Real Determinant(void) const;
public:
	union
	{
        struct
        {
            Real _11, _12;
            Real _21, _22;
        };
		/// 以[2][2]二维数组形式的数据
		Real m[2][2];
		/// 以[4]一维数组形式的数据
		Real m2[4];
	};

public:
	/// 零矩阵，全部数值为0
	static const Matrix2 ZERO;
	/// 单位矩阵，对角线上数值为1，其余为0
	static const Matrix2 IDENTITY;
};
#if _MSC_VER >= 1200
#pragma warning(pop)
#endif

/// \brief 二维向量和二维矩阵相乘，得到一个新向量。
template <typename Real>
Vector2<Real> operator * (const Vector2<Real>& v, const Matrix2<Real>& mtx);

/// \brief 矩阵数乘运算
template <typename Real>
Matrix2<Real> operator* (Real fScalar, const Matrix2<Real>& mtx);

/// float类型的二维矩阵
typedef Matrix2<float> Matrix2f;

/// double类型的二维矩阵
typedef Matrix2<double> Matrix2d;

#include "Matrix2.inl"

NS_KMATH_END

#endif // __MATRIX2_H__
