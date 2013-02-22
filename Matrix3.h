/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Matrix3.h
/// \brief   3×3矩阵(右乘矩阵)，没有提供缩放控制，并在所有与旋转有关的方法中假定矩阵正交
/// \author  DaiMingzhuang
/// \version 
/// \date    2009-02
/////////////////////////////////////////////////////////////////////////////////
#ifndef __MATRIX3_H__
#define __MATRIX3_H__

#include "MathDef.h"

NS_KMATH_BEGIN

#if _MSC_VER >= 1200
#pragma warning(push)
#pragma warning(disable:4201) // anonymous unions warning
#endif
/////////////////////////////////////////////////////////////////////////////////
/// \class Matrix3
/// \brief 3×3的矩阵(左手系，右乘矩阵)
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
class Matrix3
{
public:
	/// \brief 构造函数，未初始化。
	Matrix3();
	/// \brief 构造函数，设置每个数据的值。
	Matrix3(Real m00, Real m01, Real m02,
			Real m10, Real m11, Real m12,
			Real m20, Real m21, Real m22);
	/// \brief 拷贝构造函数。
	Matrix3(const Matrix3& mtx);

	/// \brief 设置每个数据的值。
	Matrix3& Set(Real m00, Real m01, Real m02,
						Real m10, Real m11, Real m12,
						Real m20, Real m21, Real m22);

	/// \brief 获取第row行，第col列的数据。
	Real operator() (int row, int col) const;
	/// \brief 获取第row行，第col列的数据，可更改。
	Real& operator() (int row, int col);

	/// \brief 获取第row行的行向量
	const Vector3<Real>& operator[] (int row) const;

	/// \brief 获取第row行的行向量
	Vector3<Real>& operator[] (int row);

	/// \brief 拷贝函数。
	Matrix3& operator = (const Matrix3& mtx);
	/// \brief 比较两个矩阵是否相等。
	bool operator == (const Matrix3& mtx) const;
	/// \brief 比较两个矩阵是否不等。
	bool operator != (const Matrix3& mtx) const;

	/// \brief 矩阵加法，按对应位置相加。
	Matrix3 operator + (const Matrix3& mtx) const;
	/// \brief 矩阵减法，按对应位置相减。
	Matrix3 operator - (const Matrix3& mtx) const;
	/// \brief 矩阵乘法。
	Matrix3 operator * (const Matrix3& mtx) const;
	/// \brief 矩阵乘以数值，所有数据乘以s。
	Matrix3 operator * (Real s) const;
	/// \brief 矩阵除以数值，所有数据除以s。
	Matrix3 operator / (Real s) const;
	/// \brief 矩阵求反，所有数据取相反值。
	Matrix3 operator - (void) const;

	/// \brief 矩阵加法，按对应位置相加，结果赋给自己。
	Matrix3& operator += (const Matrix3& mtx);
	/// \brief 矩阵减法，按对应位置相减，结果赋给自己。
	Matrix3& operator -= (const Matrix3& mtx);
	/// \brief 矩阵乘法，结果赋给自己。
	Matrix3& operator *= (const Matrix3& mtx);
	/// \brief 矩阵乘以数值，所有数据乘以s，结果赋给自己。
	Matrix3& operator *= (Real s);
	/// \brief 矩阵除以数值，所有数据除以s，结果赋给自己。
	Matrix3& operator /= (Real s);

public:
	/// \brief 设置为单位矩阵。
	Matrix3& SetIdentity();
	
	/// \brief 设置沿X轴的旋转矩阵。
    /// \param[in] ang 绕轴的转角，单位使用弧度制
	Matrix3& SetRotationX(Real ang);
	/// \brief 设置沿Y轴的旋转矩阵。
    /// \param[in] ang 绕轴的转角，单位使用弧度制
	Matrix3& SetRotationY(Real ang);
	/// \brief 设置沿Z轴的旋转矩阵。
    /// \param[in] ang 绕轴的转角，单位使用弧度制
	Matrix3& SetRotationZ(Real ang);

	/// \brief 设置沿空间某根轴的旋转矩阵。
    /// \param[in] axis 指定轴正方向的向量
    /// \param[in] ang 绕轴的转角，单位使用弧度制
	Matrix3& SetRotationAxis(const Vector3<Real>& axis, Real ang);

    /// \brief 设置能够从dir0方向转到dir1方向的旋转矩阵。
    /// \param[in] dir0 初始方向
    /// \param[in] dir1 目标方向
    Matrix3& SetRotationBy2Vector(const Vector3<Real>& dir0, const Vector3<Real>& dir1);

	/// \brief 设置矩阵旋转参数根据四元数
	Matrix3& SetRotationQuaternion(const Quaternion<Real>& qtn);

	/// \brief 取得矩阵旋转参数根据四元数
	Quaternion<Real> GetRotationQuaternion() const;

	/// \brief 设置旋转矩阵根据欧拉角，按XYZ的次序。
	Matrix3& SetEulerXYZ(const Vector3<Real>& angles);
	/// \brief 设置旋转矩阵根据欧拉角，按XZY的次序。
	Matrix3& SetEulerXZY(const Vector3<Real>& angles);
	/// \brief 设置旋转矩阵根据欧拉角，按YXZ的次序。
	Matrix3& SetEulerYXZ(const Vector3<Real>& angles);
	/// \brief 设置旋转矩阵根据欧拉角，按YZX的次序。
	Matrix3& SetEulerYZX(const Vector3<Real>& angles);
	/// \brief 设置旋转矩阵根据欧拉角，按ZXY的次序。
	Matrix3& SetEulerZXY(const Vector3<Real>& angles);
	/// \brief 设置旋转矩阵根据欧拉角，按ZYX的次序。
	Matrix3& SetEulerZYX(const Vector3<Real>& angles);

	/// \brief 析出旋转分量的欧拉角表示，按XYZ的次序。
	Vector3<Real> GetEulerXYZ() const;
	/// \brief 析出旋转分量的欧拉角表示，按XZY的次序。
	Vector3<Real> GetEulerXZY() const;
	/// \brief 析出旋转分量的欧拉角表示，按YXZ的次序。
	Vector3<Real> GetEulerYXZ() const;
	/// \brief 析出旋转分量的欧拉角表示，按YZX的次序。
	Vector3<Real> GetEulerYZX() const;
	/// \brief 析出旋转分量的欧拉角表示，按ZXY的次序。
	Vector3<Real> GetEulerZXY() const;
	/// \brief 析出旋转分量的欧拉角表示，按ZYX的次序。
	Vector3<Real> GetEulerZYX() const;

	/// \brief 求三维矩阵的转置矩阵。
	Matrix3 Transpose() const;
	/// \brief 求三维矩阵的逆矩阵。
	Matrix3 Inverse() const;
	/// \brief 求三维矩阵的行列式。
	Real Determinant() const;

    /// \brief 将三维矩阵扩展为四维矩阵。
    Matrix4<Real> ToMatrix4() const;

public:
    /// \brief 创建单位矩阵。
    static Matrix3 Identity(void) { return Matrix3<Real>::IDENTITY; }
    /// \brief 创建已知矩阵的转置矩阵。
    static Matrix3 Transpose(const Matrix3<Real>& mtx) { return mtx.Transpose(); }
    /// \brief 创建已知矩阵的逆矩阵。
    static Matrix3 Inverse(const Matrix3<Real>& mtx) { return mtx.Inverse(); }

public:
    /// \brief 创建沿X轴的旋转矩阵。
    /// \param[in] ang 绕轴的转角，单位使用弧度制
    static Matrix3 MakeRotationX(Real ang) { return Matrix3<Real>().SetRotationX(ang); }
    /// \brief 创建沿Y轴的旋转矩阵。
    /// \param[in] ang 绕轴的转角，单位使用弧度制
    static Matrix3 MakeRotationY(Real ang) { return Matrix3<Real>().SetRotationY(ang); }
    /// \brief 创建沿Z轴的旋转矩阵。
    /// \param[in] ang 绕轴的转角，单位使用弧度制
    static Matrix3 MakeRotationZ(Real ang) { return Matrix3<Real>().SetRotationZ(ang); }

    /// \brief 创建沿空间中某根轴的旋转矩阵。
    /// \param[in] axis 指定轴正方向的向量
    /// \param[in] ang 绕轴的转角，单位使用弧度制
    static Matrix3 MakeRotationAxis(const Vector3<Real>& axis, Real ang) { return Matrix3<Real>().SetRotationAxis(axis, ang); }

    /// \brief 创建能够从dir0方向转到dir1方向的旋转矩阵。
    /// \param[in] dir0 初始方向
    /// \param[in] dir1 目标方向
    static Matrix3 MakeRotationBy2Vector(const Vector3<Real>& dir0, const Vector3<Real>& dir1) { return Matrix3<Real>().SetRotationBy2Vector(dir0, dir1); }

	/// \brief 创建旋转矩阵根据欧拉角，按XYZ的次序。
	static Matrix3 MakeRotationEulerXYZ(const Vector3<Real>& angles) { return Matrix3<Real>().SetEulerXYZ(angles); }
	/// \brief 创建旋转矩阵根据欧拉角，按XZY的次序。
	static Matrix3 MakeRotationEulerXZY(const Vector3<Real>& angles) { return Matrix3<Real>().SetEulerXZY(angles); }
	/// \brief 创建旋转矩阵根据欧拉角，按YXZ的次序。
	static Matrix3 MakeRotationEulerYXZ(const Vector3<Real>& angles) { return Matrix3<Real>().SetEulerYXZ(angles); }
	/// \brief 创建旋转矩阵根据欧拉角，按YZX的次序。
	static Matrix3 MakeRotationEulerYZX(const Vector3<Real>& angles) { return Matrix3<Real>().SetEulerYZX(angles); }
	/// \brief 创建旋转矩阵根据欧拉角，按ZXY的次序。
	static Matrix3 MakeRotationEulerZXY(const Vector3<Real>& angles) { return Matrix3<Real>().SetEulerZXY(angles); }
	/// \brief 创建旋转矩阵根据欧拉角，按ZYX的次序。
	static Matrix3 MakeRotationEulerZYX(const Vector3<Real>& angles) { return Matrix3<Real>().SetEulerZYX(angles); }

public:
	union
	{
        struct
        {
            Real _11, _12, _13;
            Real _21, _22, _23;
            Real _31, _32, _33;
        };
		/// 以[3][3]二维数组形式的数据
		Real m[3][3];
		/// 以[9]一维数组形式的数据
		Real m2[9];
	};

public:
	/// 零矩阵，全部数值为0
	static const Matrix3 ZERO;
	/// 单位矩阵，对角线上数值为1，其余为0
	static const Matrix3 IDENTITY;
};

#if _MSC_VER >= 1200
#pragma warning(pop)
#endif

/// \brief 三维向量和三维矩阵相乘，得到一个新向量。
template <typename Real>
Vector3<Real> operator * (const Vector3<Real>& v, const Matrix3<Real>& mtx);

/// \brief 矩阵数乘运算
template <typename Real>
Matrix3<Real> operator* (Real fScalar, const Matrix3<Real>& mtx);

/// float类型的三维矩阵
typedef Matrix3<float> Matrix3f;

/// double类型的三维矩阵
typedef Matrix3<double> Matrix3d;

#include "Matrix3.inl"

NS_KMATH_END

#endif // __MATRIX3_H__
