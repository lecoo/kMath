/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Matrix4.h
/// \brief   4×4矩阵(右乘矩阵)
/// \author  DaiMingzhuang
/// \version 
/// \date    2009-02
/////////////////////////////////////////////////////////////////////////////////
#ifndef __MATRIX4_H__
#define __MATRIX4_H__

#include "MathDef.h"

NS_KMATH_BEGIN

#if _MSC_VER >= 1200
#pragma warning(push)
#pragma warning(disable:4201) // anonymous unions warning
#endif
/////////////////////////////////////////////////////////////////////////////////
/// \class Matrix4
/// \brief 4×4的矩阵(右乘矩阵)
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
class Matrix4
{
public:
	/// \brief 构造函数，未初始化。
	Matrix4(void);

	/// \brief 构造函数，设置每个数据的值。
	Matrix4(Real m00, Real m01, Real m02, Real m03,
			Real m10, Real m11, Real m12, Real m13,
			Real m20, Real m21, Real m22, Real m23,
			Real m30, Real m31, Real m32, Real m33);
	/// \brief 拷贝构造函数。
	Matrix4(const Matrix4& mtx);

	/// \brief 设置每个数据的值。
	Matrix4& Set(Real m00, Real m01, Real m02, Real m03,
						Real m10, Real m11, Real m12, Real m13,
						Real m20, Real m21, Real m22, Real m23,
						Real m30, Real m31, Real m32, Real m33);

	/// \brief 获取第row行，第col列的数据。
	Real operator() (int row, int col) const;
	/// \brief 获取第row行，第col列的数据，可更改。
	Real& operator() (int row, int col);

	/// \brief 获取第row行的行向量
	const Vector4<Real>& operator[] (int row) const;

	/// \brief 获取第row行的行向量
	Vector4<Real>& operator[] (int row);

	/// \brief 拷贝函数。
	Matrix4& operator = (const Matrix4& mtx);
	/// \brief 比较两个矩阵是否相等。
	bool operator == (const Matrix4& mtx) const;
	/// \brief 比较两个矩阵是否不等。
	bool operator != (const Matrix4& mtx) const;

	/// \brief 将三维矩阵转换为四维矩阵。
	Matrix4& FromMatrix3(const Matrix3<Real>& mtx);
	/// \brief 将四维矩阵转换为三维矩阵。
	Matrix3<Real> ToMatrix3(void) const;

	/// \brief 矩阵加法，按对应位置相加。
	Matrix4 operator + (const Matrix4& mtx) const;
	/// \brief 矩阵减法，按对应位置相减。
	Matrix4 operator - (const Matrix4& mtx) const;
	/// \brief 矩阵乘法。
    Matrix4 operator * (const Matrix4& mtx) const;

	/// \brief 矩阵乘以数值，所有数据乘以s。
	Matrix4 operator * (Real s) const;
	/// \brief 矩阵除以数值，所有数据除以s。
	Matrix4 operator / (Real s) const;
	/// \brief 矩阵求反，所有数据取相反值。
	Matrix4 operator - (void) const;

	/// \brief 矩阵加法，按对应位置相加，结果赋给自己。
	Matrix4& operator += (const Matrix4& mtx);
	/// \brief 矩阵减法，按对应位置相减，结果赋给自己。
	Matrix4& operator -= (const Matrix4& mtx);
	/// \brief 矩阵乘法，结果赋给自己。
	Matrix4& operator *= (const Matrix4& mtx);
	/// \brief 矩阵乘以数值，所有数据乘以s，结果赋给自己。
	Matrix4& operator *= (Real s);
	/// \brief 矩阵除以数值，所有数据除以s，结果赋给自己。
	Matrix4& operator /= (Real s);

	/// \brief 设置为单位矩阵。
	Matrix4& SetIdentity();

	/// \brief 设置矩阵的平移分量。
	Matrix4& SetTranslation(const Vector3<Real>& v);
	/// \brief 设置矩阵的平移分量。
	Matrix4& SetTranslation(Real x, Real y, Real z);

	/// \brief 取得矩阵的平移分量。
	Vector3<Real> GetTranslation() const;

	/// \brief 设置矩阵旋转参数为绕X轴旋转ang弧度
    /// \param[in] ang 绕轴的转角，单位使用弧度制
	Matrix4& SetRotationX(Real ang);
	/// \brief 设置矩阵旋转参数为绕Y轴旋转ang弧度
    /// \param[in] ang 绕轴的转角，单位使用弧度制
	Matrix4& SetRotationY(Real ang);
	/// \brief 设置矩阵旋转参数为绕Z轴旋转ang弧度
    /// \param[in] ang 绕轴的转角，单位使用弧度制
	Matrix4& SetRotationZ(Real ang);
	/// \brief 设置矩阵旋转参数为绕axis轴旋转ang弧度
    /// \param[in] axis 指定轴正方向的向量
    /// \param[in] ang 绕轴的转角，单位使用弧度制
	Matrix4& SetRotationAxis(const Vector3<Real>& axis, Real ang);
    /// \brief 设置矩阵旋转参数为从dir0方向转到dir1方向
    /// \param[in] dir0 初始方向
    /// \param[in] dir1 目标方向
    Matrix4& SetRotationBy2Vector(const Vector3<Real>& dir0, const Vector3<Real>& dir1);

	/// \brief 设置矩阵旋转参数根据四元数
	Matrix4& SetRotationQuaternion(const Quaternion<Real>& qtn);

	/// \brief 取得矩阵旋转参数根据四元数
	Quaternion<Real> GetRotationQuaternion() const;

	/// \brief 设置矩阵旋转参数根据欧拉角，按XYZ的次序。
	Matrix4& SetEulerXYZ(const Vector3<Real>& angles);
	/// \brief 设置矩阵旋转参数根据欧拉角，按XZY的次序。
	Matrix4& SetEulerXZY(const Vector3<Real>& angles);
	/// \brief 设置矩阵旋转参数根据欧拉角，按YXZ的次序。
	Matrix4& SetEulerYXZ(const Vector3<Real>& angles);
	/// \brief 设置矩阵旋转参数根据欧拉角，按YZX的次序。
	Matrix4& SetEulerYZX(const Vector3<Real>& angles);
	/// \brief 设置矩阵旋转参数根据欧拉角，按ZXY的次序。
	Matrix4& SetEulerZXY(const Vector3<Real>& angles);
	/// \brief 设置矩阵旋转参数根据欧拉角，按ZYX的次序。
	Matrix4& SetEulerZYX(const Vector3<Real>& angles);

	/// \brief 析出矩阵旋转参数的欧拉角表示，按XYZ的次序。
	Vector3<Real> GetEulerXYZ() const;
	/// \brief 析出矩阵旋转参数的欧拉角表示，按XZY的次序。
	Vector3<Real> GetEulerXZY() const;
	/// \brief 析出矩阵旋转参数的欧拉角表示，按YXZ的次序。
	Vector3<Real> GetEulerYXZ() const;
	/// \brief 析出矩阵旋转参数的欧拉角表示，按YZX的次序。
	Vector3<Real> GetEulerYZX() const;
	/// \brief 析出矩阵旋转参数的欧拉角表示，按ZXY的次序。
	Vector3<Real> GetEulerZXY() const;
	/// \brief 析出矩阵旋转参数的欧拉角表示，按ZYX的次序。
	Vector3<Real> GetEulerZYX() const;

	/// \brief 设置矩阵的缩放参数
	Matrix4& SetScaling(const Vector3<Real>& scaling);
	/// \brief 析出矩阵的缩放参数
	Vector3<Real> GetScaling() const;

	/// \brief 求四维矩阵的转置矩阵。
	Matrix4 Transpose(void) const;
	/// \brief 求四维矩阵的逆矩阵。
	Matrix4 Inverse(void) const;
	/// \brief 求四维矩阵的行列式。
	Real Determinant(void) const;

	/// \brief 清除矩阵缩放参数。
	Matrix4& ClearScaling(void);

	/// \brief 清除矩阵旋转参数。
	Matrix4& ClearRotation(void);

	/// \brief 清除矩阵平移分量。
	Matrix4& ClearTranslation(void);

	/// \brief 在当前矩阵的基础上做preScale变换
	Matrix4& Scale(const Vector3<Real>& v);

	/// \brief 在当前矩阵的基础上做postTranslate变换
	Matrix4& Translate(const Vector3<Real>& v);

public:
    /// \brief 返回单位矩阵。
    static Matrix4 Identity(void) { return Matrix4<Real>::IDENTITY; }
    /// \brief 返回已知矩阵的转置矩阵。
    static Matrix4 Transpose(const Matrix4<Real>& mtx) { return mtx.Transpose(); }
    /// \brief 返回已知矩阵的逆矩阵。
    static Matrix4 Inverse(const Matrix4<Real>& mtx) { return mtx.Inverse(); }

public:
    /// \brief 创建平移矩阵。
    static Matrix4 MakeTranslation(const Vector3<Real>& v) { return Matrix4<Real>().SetIdentity().SetTranslation(v); }
    /// \brief 创建平移矩阵。
    static Matrix4 MakeTranslation(Real x, Real y, Real z) { return Matrix4<Real>().SetIdentity().SetTranslation(x, y, z); }

    /// \brief 创建沿X轴的旋转矩阵。
    /// \param[in] ang 绕轴的转角，单位使用弧度制
    static Matrix4 MakeRotationX(Real ang) { return Matrix4<Real>().SetIdentity().SetRotationX(ang); }
    /// \brief 创建沿Y轴的旋转矩阵。
    /// \param[in] ang 绕轴的转角，单位使用弧度制
    static Matrix4 MakeRotationY(Real ang) { return Matrix4<Real>().SetIdentity().SetRotationY(ang); }
    /// \brief 创建沿Z轴的旋转矩阵。
    /// \param[in] ang 绕轴的转角，单位使用弧度制
    static Matrix4 MakeRotationZ(Real ang) { return Matrix4<Real>().SetIdentity().SetRotationZ(ang); }

    /// \brief 创建沿空间中某根轴的旋转矩阵。
    /// \param[in] axis 指定轴正方向的向量
    /// \param[in] ang 绕轴的转角，单位使用弧度制
    static Matrix4 MakeRotationAxis(const Vector3<Real>& axis, Real ang) { return Matrix4<Real>().SetIdentity().SetRotationAxis(axis, ang); }

    /// \brief 创建能够从dir0方向转到dir1方向的旋转矩阵。
    /// \param[in] dir0 初始方向
    /// \param[in] dir1 目标方向
    static Matrix4 MakeRotationBy2Vector(const Vector3<Real>& dir0, const Vector3<Real>& dir1) { return Matrix4<Real>().SetIdentity().SetRotationBy2Vector(dir0, dir1); }

	/// \brief 创建旋转矩阵根据四元数
	static Matrix4 MakeRotationQuaternion(const Quaternion<Real>& qtn) { return Matrix4<Real>().SetIdentity().SetRotationQuaternion(qtn); }

	/// \brief 创建旋转矩阵根据欧拉角，按XYZ的次序。
	static Matrix4 MakeRotationEulerXYZ(const Vector3<Real>& angles) { return Matrix4<Real>().SetIdentity().SetEulerXYZ(angles); }
	/// \brief 创建旋转矩阵根据欧拉角，按XZY的次序。
	static Matrix4 MakeRotationEulerXZY(const Vector3<Real>& angles) { return Matrix4<Real>().SetIdentity().SetEulerXZY(angles); }
	/// \brief 创建旋转矩阵根据欧拉角，按YXZ的次序。
	static Matrix4 MakeRotationEulerYXZ(const Vector3<Real>& angles) { return Matrix4<Real>().SetIdentity().SetEulerYXZ(angles); }
	/// \brief 创建旋转矩阵根据欧拉角，按YZX的次序。
	static Matrix4 MakeRotationEulerYZX(const Vector3<Real>& angles) { return Matrix4<Real>().SetIdentity().SetEulerYZX(angles); }
	/// \brief 创建旋转矩阵根据欧拉角，按ZXY的次序。
	static Matrix4 MakeRotationEulerZXY(const Vector3<Real>& angles) { return Matrix4<Real>().SetIdentity().SetEulerZXY(angles); }
	/// \brief 创建旋转矩阵根据欧拉角，按ZYX的次序。
	static Matrix4 MakeRotationEulerZYX(const Vector3<Real>& angles) { return Matrix4<Real>().SetIdentity().SetEulerZYX(angles); }

    /// \brief 创建缩放矩阵，缩放矩阵转置对称。
    static Matrix4 MakeScaling(const Vector3<Real>& v);
    /// \brief 创建缩放矩阵，缩放矩阵转置对称。
    static Matrix4 MakeScaling(Real x, Real y, Real z);

public:
    /// \brief 创建相对于平面的对称矩阵，根据Real[4]。
    /// \param[in] p    p[0]，p[1]，p[2]表示平面方向，p[3]表示平面位置。
    static Matrix4 Reflect(const Real* p);

    /// \brief    Builds a left-handed, look-at matrix.
    /// \param[in] eye  视角的位置。
    /// \param[in] at   目标点的位置。
    /// \param[in] up   辅助向量，一般为Y轴正向。
    static Matrix4 LookAtLH(const Vector3<Real>& eye, const Vector3<Real>& at, const Vector3<Real>& up);

    /// \brief    Builds a right-handed, look-at matrix.
    /// \param[in] eye  视角的位置。
    /// \param[in] at   目标点的位置。
    /// \param[in] up   辅助向量，一般为Y轴正向。
    static Matrix4 LookAtRH(const Vector3<Real>& eye, const Vector3<Real>& at, const Vector3<Real>& up);

    /// \brief    Builds a left-handed perspective projection matrix based on a field of view.
    /// \param[in] w      Width of the view volume at the near view-plane.
    /// \param[in] h      Height of the view volume at the near view-plane.
    /// \param[in] zn     Z-value of the near view-plane.
    /// \param[in] zf     Z-value of the far view-plane.
    static Matrix4 PerspectiveLH(Real w, Real h, Real zn, Real zf);

    /// \brief    Builds a right-handed perspective projection matrix based on a field of view.
    /// \param[in] w      Width of the view volume at the near view-plane.
    /// \param[in] h      Height of the view volume at the near view-plane.
    /// \param[in] zn     Z-value of the near view-plane.
    /// \param[in] zf     Z-value of the far view-plane.
    static Matrix4 PerspectiveRH(Real w, Real h, Real zn, Real zf);

    /// \brief    Builds a left-handed perspective projection matrix based on a field of view.
    /// \param[in] fovY   Y方向张角的弧度。
    /// \param[in] aspect 视空间的宽度和长度的比例。
    /// \param[in] zn     裁剪近平面。
    /// \param[in] zf     裁剪远平面。
    static Matrix4 PerspectiveFovLH(Real fovY, Real aspect, Real zn, Real zf);

    /// \brief    Builds a right-handed perspective projection matrix based on a field of view.
    /// \param[in] fovY   Y方向张角的弧度。
    /// \param[in] aspect 视空间的宽度和长度的比例。
    /// \param[in] zn     裁剪近平面。
    /// \param[in] zf     裁剪远平面。
    static Matrix4 PerspectiveFovRH(Real fovY, Real aspect, Real zn, Real zf);

    /// \brief    Builds a left-handed orthogonal projection matrix.
    /// \param[in] w      Width of the view volume at the near view-plane.
    /// \param[in] h      Height of the view volume at the near view-plane.
    /// \param[in] zn     Minimum z-value of the view volume which is referred to as z-near.
    /// \param[in] zf     Maximum z-value of the view volume which is referred to as z-far.
    static Matrix4 OrthoLH(Real w, Real h, Real zn, Real zf);

    /// \brief    Builds a right-handed orthogonal projection matrix.
    /// \param[in] w      Width of the view volume at the near view-plane.
    /// \param[in] h      Height of the view volume at the near view-plane.
    /// \param[in] zn     Minimum z-value of the view volume which is referred to as z-near.
    /// \param[in] zf     Maximum z-value of the view volume which is referred to as z-far.
    static Matrix4 OrthoRH(Real w, Real h, Real zn, Real zf);

public:
	union
	{
        struct
        {
            Real _11, _12, _13, _14;
            Real _21, _22, _23, _24;
            Real _31, _32, _33, _34;
            Real _41, _42, _43, _44;
        };
		/// 以[4][4]二维数组形式的数据
		Real m[4][4];
		/// 以[16]一维数组形式的数据
		Real m2[16];
	};

public:
	/// 零矩阵，全部数值为0
	static const Matrix4 ZERO;
	/// 单位矩阵，对角线上数值为1，其余为0
	static const Matrix4 IDENTITY;
};

#if _MSC_VER >= 1200
#pragma warning(pop)
#endif

/// \brief 三维向量(扩展为四维点 w=1)和四维矩阵相乘，得到一个新三维向量，并对其齐次化
template <typename Real>
Vector3<Real> operator * (const Vector3<Real>& v, const Matrix4<Real>& mtx);

/// \brief 四维向量和四维矩阵相乘，得到一个新四维向量。
template <typename Real>
Vector4<Real> operator * (const Vector4<Real>& v, const Matrix4<Real>& mtx);

/// \brief 矩阵数乘运算
template <typename Real>
Matrix4<Real> operator* (Real fScalar, const Matrix4<Real>& mtx);

/// float类型的四维矩阵
typedef Matrix4<float> Matrix4f;

/// double类型的四维矩阵
typedef Matrix4<double> Matrix4d;

#include "Matrix4.inl"

#ifdef USE_DX_MATH
#include "Matrix4f_DX.inl"  //用DX特化Matrix4<float>的实现
#endif

NS_KMATH_END

#endif // __MATRIX4_H__
