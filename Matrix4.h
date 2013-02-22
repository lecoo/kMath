/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Matrix4.h
/// \brief   4��4����(�ҳ˾���)
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
/// \brief 4��4�ľ���(�ҳ˾���)
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
class Matrix4
{
public:
	/// \brief ���캯����δ��ʼ����
	Matrix4(void);

	/// \brief ���캯��������ÿ�����ݵ�ֵ��
	Matrix4(Real m00, Real m01, Real m02, Real m03,
			Real m10, Real m11, Real m12, Real m13,
			Real m20, Real m21, Real m22, Real m23,
			Real m30, Real m31, Real m32, Real m33);
	/// \brief �������캯����
	Matrix4(const Matrix4& mtx);

	/// \brief ����ÿ�����ݵ�ֵ��
	Matrix4& Set(Real m00, Real m01, Real m02, Real m03,
						Real m10, Real m11, Real m12, Real m13,
						Real m20, Real m21, Real m22, Real m23,
						Real m30, Real m31, Real m32, Real m33);

	/// \brief ��ȡ��row�У���col�е����ݡ�
	Real operator() (int row, int col) const;
	/// \brief ��ȡ��row�У���col�е����ݣ��ɸ��ġ�
	Real& operator() (int row, int col);

	/// \brief ��ȡ��row�е�������
	const Vector4<Real>& operator[] (int row) const;

	/// \brief ��ȡ��row�е�������
	Vector4<Real>& operator[] (int row);

	/// \brief ����������
	Matrix4& operator = (const Matrix4& mtx);
	/// \brief �Ƚ����������Ƿ���ȡ�
	bool operator == (const Matrix4& mtx) const;
	/// \brief �Ƚ����������Ƿ񲻵ȡ�
	bool operator != (const Matrix4& mtx) const;

	/// \brief ����ά����ת��Ϊ��ά����
	Matrix4& FromMatrix3(const Matrix3<Real>& mtx);
	/// \brief ����ά����ת��Ϊ��ά����
	Matrix3<Real> ToMatrix3(void) const;

	/// \brief ����ӷ�������Ӧλ����ӡ�
	Matrix4 operator + (const Matrix4& mtx) const;
	/// \brief �������������Ӧλ�������
	Matrix4 operator - (const Matrix4& mtx) const;
	/// \brief ����˷���
    Matrix4 operator * (const Matrix4& mtx) const;

	/// \brief ���������ֵ���������ݳ���s��
	Matrix4 operator * (Real s) const;
	/// \brief ���������ֵ���������ݳ���s��
	Matrix4 operator / (Real s) const;
	/// \brief �����󷴣���������ȡ�෴ֵ��
	Matrix4 operator - (void) const;

	/// \brief ����ӷ�������Ӧλ����ӣ���������Լ���
	Matrix4& operator += (const Matrix4& mtx);
	/// \brief �������������Ӧλ���������������Լ���
	Matrix4& operator -= (const Matrix4& mtx);
	/// \brief ����˷�����������Լ���
	Matrix4& operator *= (const Matrix4& mtx);
	/// \brief ���������ֵ���������ݳ���s����������Լ���
	Matrix4& operator *= (Real s);
	/// \brief ���������ֵ���������ݳ���s����������Լ���
	Matrix4& operator /= (Real s);

	/// \brief ����Ϊ��λ����
	Matrix4& SetIdentity();

	/// \brief ���þ����ƽ�Ʒ�����
	Matrix4& SetTranslation(const Vector3<Real>& v);
	/// \brief ���þ����ƽ�Ʒ�����
	Matrix4& SetTranslation(Real x, Real y, Real z);

	/// \brief ȡ�þ����ƽ�Ʒ�����
	Vector3<Real> GetTranslation() const;

	/// \brief ���þ�����ת����Ϊ��X����תang����
    /// \param[in] ang �����ת�ǣ���λʹ�û�����
	Matrix4& SetRotationX(Real ang);
	/// \brief ���þ�����ת����Ϊ��Y����תang����
    /// \param[in] ang �����ת�ǣ���λʹ�û�����
	Matrix4& SetRotationY(Real ang);
	/// \brief ���þ�����ת����Ϊ��Z����תang����
    /// \param[in] ang �����ת�ǣ���λʹ�û�����
	Matrix4& SetRotationZ(Real ang);
	/// \brief ���þ�����ת����Ϊ��axis����תang����
    /// \param[in] axis ָ���������������
    /// \param[in] ang �����ת�ǣ���λʹ�û�����
	Matrix4& SetRotationAxis(const Vector3<Real>& axis, Real ang);
    /// \brief ���þ�����ת����Ϊ��dir0����ת��dir1����
    /// \param[in] dir0 ��ʼ����
    /// \param[in] dir1 Ŀ�귽��
    Matrix4& SetRotationBy2Vector(const Vector3<Real>& dir0, const Vector3<Real>& dir1);

	/// \brief ���þ�����ת����������Ԫ��
	Matrix4& SetRotationQuaternion(const Quaternion<Real>& qtn);

	/// \brief ȡ�þ�����ת����������Ԫ��
	Quaternion<Real> GetRotationQuaternion() const;

	/// \brief ���þ�����ת��������ŷ���ǣ���XYZ�Ĵ���
	Matrix4& SetEulerXYZ(const Vector3<Real>& angles);
	/// \brief ���þ�����ת��������ŷ���ǣ���XZY�Ĵ���
	Matrix4& SetEulerXZY(const Vector3<Real>& angles);
	/// \brief ���þ�����ת��������ŷ���ǣ���YXZ�Ĵ���
	Matrix4& SetEulerYXZ(const Vector3<Real>& angles);
	/// \brief ���þ�����ת��������ŷ���ǣ���YZX�Ĵ���
	Matrix4& SetEulerYZX(const Vector3<Real>& angles);
	/// \brief ���þ�����ת��������ŷ���ǣ���ZXY�Ĵ���
	Matrix4& SetEulerZXY(const Vector3<Real>& angles);
	/// \brief ���þ�����ת��������ŷ���ǣ���ZYX�Ĵ���
	Matrix4& SetEulerZYX(const Vector3<Real>& angles);

	/// \brief ����������ת������ŷ���Ǳ�ʾ����XYZ�Ĵ���
	Vector3<Real> GetEulerXYZ() const;
	/// \brief ����������ת������ŷ���Ǳ�ʾ����XZY�Ĵ���
	Vector3<Real> GetEulerXZY() const;
	/// \brief ����������ת������ŷ���Ǳ�ʾ����YXZ�Ĵ���
	Vector3<Real> GetEulerYXZ() const;
	/// \brief ����������ת������ŷ���Ǳ�ʾ����YZX�Ĵ���
	Vector3<Real> GetEulerYZX() const;
	/// \brief ����������ת������ŷ���Ǳ�ʾ����ZXY�Ĵ���
	Vector3<Real> GetEulerZXY() const;
	/// \brief ����������ת������ŷ���Ǳ�ʾ����ZYX�Ĵ���
	Vector3<Real> GetEulerZYX() const;

	/// \brief ���þ�������Ų���
	Matrix4& SetScaling(const Vector3<Real>& scaling);
	/// \brief ������������Ų���
	Vector3<Real> GetScaling() const;

	/// \brief ����ά�����ת�þ���
	Matrix4 Transpose(void) const;
	/// \brief ����ά����������
	Matrix4 Inverse(void) const;
	/// \brief ����ά���������ʽ��
	Real Determinant(void) const;

	/// \brief ����������Ų�����
	Matrix4& ClearScaling(void);

	/// \brief ���������ת������
	Matrix4& ClearRotation(void);

	/// \brief �������ƽ�Ʒ�����
	Matrix4& ClearTranslation(void);

	/// \brief �ڵ�ǰ����Ļ�������preScale�任
	Matrix4& Scale(const Vector3<Real>& v);

	/// \brief �ڵ�ǰ����Ļ�������postTranslate�任
	Matrix4& Translate(const Vector3<Real>& v);

public:
    /// \brief ���ص�λ����
    static Matrix4 Identity(void) { return Matrix4<Real>::IDENTITY; }
    /// \brief ������֪�����ת�þ���
    static Matrix4 Transpose(const Matrix4<Real>& mtx) { return mtx.Transpose(); }
    /// \brief ������֪����������
    static Matrix4 Inverse(const Matrix4<Real>& mtx) { return mtx.Inverse(); }

public:
    /// \brief ����ƽ�ƾ���
    static Matrix4 MakeTranslation(const Vector3<Real>& v) { return Matrix4<Real>().SetIdentity().SetTranslation(v); }
    /// \brief ����ƽ�ƾ���
    static Matrix4 MakeTranslation(Real x, Real y, Real z) { return Matrix4<Real>().SetIdentity().SetTranslation(x, y, z); }

    /// \brief ������X�����ת����
    /// \param[in] ang �����ת�ǣ���λʹ�û�����
    static Matrix4 MakeRotationX(Real ang) { return Matrix4<Real>().SetIdentity().SetRotationX(ang); }
    /// \brief ������Y�����ת����
    /// \param[in] ang �����ת�ǣ���λʹ�û�����
    static Matrix4 MakeRotationY(Real ang) { return Matrix4<Real>().SetIdentity().SetRotationY(ang); }
    /// \brief ������Z�����ת����
    /// \param[in] ang �����ת�ǣ���λʹ�û�����
    static Matrix4 MakeRotationZ(Real ang) { return Matrix4<Real>().SetIdentity().SetRotationZ(ang); }

    /// \brief �����ؿռ���ĳ�������ת����
    /// \param[in] axis ָ���������������
    /// \param[in] ang �����ת�ǣ���λʹ�û�����
    static Matrix4 MakeRotationAxis(const Vector3<Real>& axis, Real ang) { return Matrix4<Real>().SetIdentity().SetRotationAxis(axis, ang); }

    /// \brief �����ܹ���dir0����ת��dir1�������ת����
    /// \param[in] dir0 ��ʼ����
    /// \param[in] dir1 Ŀ�귽��
    static Matrix4 MakeRotationBy2Vector(const Vector3<Real>& dir0, const Vector3<Real>& dir1) { return Matrix4<Real>().SetIdentity().SetRotationBy2Vector(dir0, dir1); }

	/// \brief ������ת���������Ԫ��
	static Matrix4 MakeRotationQuaternion(const Quaternion<Real>& qtn) { return Matrix4<Real>().SetIdentity().SetRotationQuaternion(qtn); }

	/// \brief ������ת�������ŷ���ǣ���XYZ�Ĵ���
	static Matrix4 MakeRotationEulerXYZ(const Vector3<Real>& angles) { return Matrix4<Real>().SetIdentity().SetEulerXYZ(angles); }
	/// \brief ������ת�������ŷ���ǣ���XZY�Ĵ���
	static Matrix4 MakeRotationEulerXZY(const Vector3<Real>& angles) { return Matrix4<Real>().SetIdentity().SetEulerXZY(angles); }
	/// \brief ������ת�������ŷ���ǣ���YXZ�Ĵ���
	static Matrix4 MakeRotationEulerYXZ(const Vector3<Real>& angles) { return Matrix4<Real>().SetIdentity().SetEulerYXZ(angles); }
	/// \brief ������ת�������ŷ���ǣ���YZX�Ĵ���
	static Matrix4 MakeRotationEulerYZX(const Vector3<Real>& angles) { return Matrix4<Real>().SetIdentity().SetEulerYZX(angles); }
	/// \brief ������ת�������ŷ���ǣ���ZXY�Ĵ���
	static Matrix4 MakeRotationEulerZXY(const Vector3<Real>& angles) { return Matrix4<Real>().SetIdentity().SetEulerZXY(angles); }
	/// \brief ������ת�������ŷ���ǣ���ZYX�Ĵ���
	static Matrix4 MakeRotationEulerZYX(const Vector3<Real>& angles) { return Matrix4<Real>().SetIdentity().SetEulerZYX(angles); }

    /// \brief �������ž������ž���ת�öԳơ�
    static Matrix4 MakeScaling(const Vector3<Real>& v);
    /// \brief �������ž������ž���ת�öԳơ�
    static Matrix4 MakeScaling(Real x, Real y, Real z);

public:
    /// \brief ���������ƽ��ĶԳƾ��󣬸���Real[4]��
    /// \param[in] p    p[0]��p[1]��p[2]��ʾƽ�淽��p[3]��ʾƽ��λ�á�
    static Matrix4 Reflect(const Real* p);

    /// \brief    Builds a left-handed, look-at matrix.
    /// \param[in] eye  �ӽǵ�λ�á�
    /// \param[in] at   Ŀ����λ�á�
    /// \param[in] up   ����������һ��ΪY������
    static Matrix4 LookAtLH(const Vector3<Real>& eye, const Vector3<Real>& at, const Vector3<Real>& up);

    /// \brief    Builds a right-handed, look-at matrix.
    /// \param[in] eye  �ӽǵ�λ�á�
    /// \param[in] at   Ŀ����λ�á�
    /// \param[in] up   ����������һ��ΪY������
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
    /// \param[in] fovY   Y�����ŽǵĻ��ȡ�
    /// \param[in] aspect �ӿռ�Ŀ�Ⱥͳ��ȵı�����
    /// \param[in] zn     �ü���ƽ�档
    /// \param[in] zf     �ü�Զƽ�档
    static Matrix4 PerspectiveFovLH(Real fovY, Real aspect, Real zn, Real zf);

    /// \brief    Builds a right-handed perspective projection matrix based on a field of view.
    /// \param[in] fovY   Y�����ŽǵĻ��ȡ�
    /// \param[in] aspect �ӿռ�Ŀ�Ⱥͳ��ȵı�����
    /// \param[in] zn     �ü���ƽ�档
    /// \param[in] zf     �ü�Զƽ�档
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
		/// ��[4][4]��ά������ʽ������
		Real m[4][4];
		/// ��[16]һά������ʽ������
		Real m2[16];
	};

public:
	/// �����ȫ����ֵΪ0
	static const Matrix4 ZERO;
	/// ��λ���󣬶Խ�������ֵΪ1������Ϊ0
	static const Matrix4 IDENTITY;
};

#if _MSC_VER >= 1200
#pragma warning(pop)
#endif

/// \brief ��ά����(��չΪ��ά�� w=1)����ά������ˣ��õ�һ������ά��������������λ�
template <typename Real>
Vector3<Real> operator * (const Vector3<Real>& v, const Matrix4<Real>& mtx);

/// \brief ��ά��������ά������ˣ��õ�һ������ά������
template <typename Real>
Vector4<Real> operator * (const Vector4<Real>& v, const Matrix4<Real>& mtx);

/// \brief ������������
template <typename Real>
Matrix4<Real> operator* (Real fScalar, const Matrix4<Real>& mtx);

/// float���͵���ά����
typedef Matrix4<float> Matrix4f;

/// double���͵���ά����
typedef Matrix4<double> Matrix4d;

#include "Matrix4.inl"

#ifdef USE_DX_MATH
#include "Matrix4f_DX.inl"  //��DX�ػ�Matrix4<float>��ʵ��
#endif

NS_KMATH_END

#endif // __MATRIX4_H__
