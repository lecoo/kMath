/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Matrix3.h
/// \brief   3��3����(�ҳ˾���)��û���ṩ���ſ��ƣ�������������ת�йصķ����мٶ���������
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
/// \brief 3��3�ľ���(����ϵ���ҳ˾���)
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
class Matrix3
{
public:
	/// \brief ���캯����δ��ʼ����
	Matrix3();
	/// \brief ���캯��������ÿ�����ݵ�ֵ��
	Matrix3(Real m00, Real m01, Real m02,
			Real m10, Real m11, Real m12,
			Real m20, Real m21, Real m22);
	/// \brief �������캯����
	Matrix3(const Matrix3& mtx);

	/// \brief ����ÿ�����ݵ�ֵ��
	Matrix3& Set(Real m00, Real m01, Real m02,
						Real m10, Real m11, Real m12,
						Real m20, Real m21, Real m22);

	/// \brief ��ȡ��row�У���col�е����ݡ�
	Real operator() (int row, int col) const;
	/// \brief ��ȡ��row�У���col�е����ݣ��ɸ��ġ�
	Real& operator() (int row, int col);

	/// \brief ��ȡ��row�е�������
	const Vector3<Real>& operator[] (int row) const;

	/// \brief ��ȡ��row�е�������
	Vector3<Real>& operator[] (int row);

	/// \brief ����������
	Matrix3& operator = (const Matrix3& mtx);
	/// \brief �Ƚ����������Ƿ���ȡ�
	bool operator == (const Matrix3& mtx) const;
	/// \brief �Ƚ����������Ƿ񲻵ȡ�
	bool operator != (const Matrix3& mtx) const;

	/// \brief ����ӷ�������Ӧλ����ӡ�
	Matrix3 operator + (const Matrix3& mtx) const;
	/// \brief �������������Ӧλ�������
	Matrix3 operator - (const Matrix3& mtx) const;
	/// \brief ����˷���
	Matrix3 operator * (const Matrix3& mtx) const;
	/// \brief ���������ֵ���������ݳ���s��
	Matrix3 operator * (Real s) const;
	/// \brief ���������ֵ���������ݳ���s��
	Matrix3 operator / (Real s) const;
	/// \brief �����󷴣���������ȡ�෴ֵ��
	Matrix3 operator - (void) const;

	/// \brief ����ӷ�������Ӧλ����ӣ���������Լ���
	Matrix3& operator += (const Matrix3& mtx);
	/// \brief �������������Ӧλ���������������Լ���
	Matrix3& operator -= (const Matrix3& mtx);
	/// \brief ����˷�����������Լ���
	Matrix3& operator *= (const Matrix3& mtx);
	/// \brief ���������ֵ���������ݳ���s����������Լ���
	Matrix3& operator *= (Real s);
	/// \brief ���������ֵ���������ݳ���s����������Լ���
	Matrix3& operator /= (Real s);

public:
	/// \brief ����Ϊ��λ����
	Matrix3& SetIdentity();
	
	/// \brief ������X�����ת����
    /// \param[in] ang �����ת�ǣ���λʹ�û�����
	Matrix3& SetRotationX(Real ang);
	/// \brief ������Y�����ת����
    /// \param[in] ang �����ת�ǣ���λʹ�û�����
	Matrix3& SetRotationY(Real ang);
	/// \brief ������Z�����ת����
    /// \param[in] ang �����ת�ǣ���λʹ�û�����
	Matrix3& SetRotationZ(Real ang);

	/// \brief �����ؿռ�ĳ�������ת����
    /// \param[in] axis ָ���������������
    /// \param[in] ang �����ת�ǣ���λʹ�û�����
	Matrix3& SetRotationAxis(const Vector3<Real>& axis, Real ang);

    /// \brief �����ܹ���dir0����ת��dir1�������ת����
    /// \param[in] dir0 ��ʼ����
    /// \param[in] dir1 Ŀ�귽��
    Matrix3& SetRotationBy2Vector(const Vector3<Real>& dir0, const Vector3<Real>& dir1);

	/// \brief ���þ�����ת����������Ԫ��
	Matrix3& SetRotationQuaternion(const Quaternion<Real>& qtn);

	/// \brief ȡ�þ�����ת����������Ԫ��
	Quaternion<Real> GetRotationQuaternion() const;

	/// \brief ������ת�������ŷ���ǣ���XYZ�Ĵ���
	Matrix3& SetEulerXYZ(const Vector3<Real>& angles);
	/// \brief ������ת�������ŷ���ǣ���XZY�Ĵ���
	Matrix3& SetEulerXZY(const Vector3<Real>& angles);
	/// \brief ������ת�������ŷ���ǣ���YXZ�Ĵ���
	Matrix3& SetEulerYXZ(const Vector3<Real>& angles);
	/// \brief ������ת�������ŷ���ǣ���YZX�Ĵ���
	Matrix3& SetEulerYZX(const Vector3<Real>& angles);
	/// \brief ������ת�������ŷ���ǣ���ZXY�Ĵ���
	Matrix3& SetEulerZXY(const Vector3<Real>& angles);
	/// \brief ������ת�������ŷ���ǣ���ZYX�Ĵ���
	Matrix3& SetEulerZYX(const Vector3<Real>& angles);

	/// \brief ������ת������ŷ���Ǳ�ʾ����XYZ�Ĵ���
	Vector3<Real> GetEulerXYZ() const;
	/// \brief ������ת������ŷ���Ǳ�ʾ����XZY�Ĵ���
	Vector3<Real> GetEulerXZY() const;
	/// \brief ������ת������ŷ���Ǳ�ʾ����YXZ�Ĵ���
	Vector3<Real> GetEulerYXZ() const;
	/// \brief ������ת������ŷ���Ǳ�ʾ����YZX�Ĵ���
	Vector3<Real> GetEulerYZX() const;
	/// \brief ������ת������ŷ���Ǳ�ʾ����ZXY�Ĵ���
	Vector3<Real> GetEulerZXY() const;
	/// \brief ������ת������ŷ���Ǳ�ʾ����ZYX�Ĵ���
	Vector3<Real> GetEulerZYX() const;

	/// \brief ����ά�����ת�þ���
	Matrix3 Transpose() const;
	/// \brief ����ά����������
	Matrix3 Inverse() const;
	/// \brief ����ά���������ʽ��
	Real Determinant() const;

    /// \brief ����ά������չΪ��ά����
    Matrix4<Real> ToMatrix4() const;

public:
    /// \brief ������λ����
    static Matrix3 Identity(void) { return Matrix3<Real>::IDENTITY; }
    /// \brief ������֪�����ת�þ���
    static Matrix3 Transpose(const Matrix3<Real>& mtx) { return mtx.Transpose(); }
    /// \brief ������֪����������
    static Matrix3 Inverse(const Matrix3<Real>& mtx) { return mtx.Inverse(); }

public:
    /// \brief ������X�����ת����
    /// \param[in] ang �����ת�ǣ���λʹ�û�����
    static Matrix3 MakeRotationX(Real ang) { return Matrix3<Real>().SetRotationX(ang); }
    /// \brief ������Y�����ת����
    /// \param[in] ang �����ת�ǣ���λʹ�û�����
    static Matrix3 MakeRotationY(Real ang) { return Matrix3<Real>().SetRotationY(ang); }
    /// \brief ������Z�����ת����
    /// \param[in] ang �����ת�ǣ���λʹ�û�����
    static Matrix3 MakeRotationZ(Real ang) { return Matrix3<Real>().SetRotationZ(ang); }

    /// \brief �����ؿռ���ĳ�������ת����
    /// \param[in] axis ָ���������������
    /// \param[in] ang �����ת�ǣ���λʹ�û�����
    static Matrix3 MakeRotationAxis(const Vector3<Real>& axis, Real ang) { return Matrix3<Real>().SetRotationAxis(axis, ang); }

    /// \brief �����ܹ���dir0����ת��dir1�������ת����
    /// \param[in] dir0 ��ʼ����
    /// \param[in] dir1 Ŀ�귽��
    static Matrix3 MakeRotationBy2Vector(const Vector3<Real>& dir0, const Vector3<Real>& dir1) { return Matrix3<Real>().SetRotationBy2Vector(dir0, dir1); }

	/// \brief ������ת�������ŷ���ǣ���XYZ�Ĵ���
	static Matrix3 MakeRotationEulerXYZ(const Vector3<Real>& angles) { return Matrix3<Real>().SetEulerXYZ(angles); }
	/// \brief ������ת�������ŷ���ǣ���XZY�Ĵ���
	static Matrix3 MakeRotationEulerXZY(const Vector3<Real>& angles) { return Matrix3<Real>().SetEulerXZY(angles); }
	/// \brief ������ת�������ŷ���ǣ���YXZ�Ĵ���
	static Matrix3 MakeRotationEulerYXZ(const Vector3<Real>& angles) { return Matrix3<Real>().SetEulerYXZ(angles); }
	/// \brief ������ת�������ŷ���ǣ���YZX�Ĵ���
	static Matrix3 MakeRotationEulerYZX(const Vector3<Real>& angles) { return Matrix3<Real>().SetEulerYZX(angles); }
	/// \brief ������ת�������ŷ���ǣ���ZXY�Ĵ���
	static Matrix3 MakeRotationEulerZXY(const Vector3<Real>& angles) { return Matrix3<Real>().SetEulerZXY(angles); }
	/// \brief ������ת�������ŷ���ǣ���ZYX�Ĵ���
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
		/// ��[3][3]��ά������ʽ������
		Real m[3][3];
		/// ��[9]һά������ʽ������
		Real m2[9];
	};

public:
	/// �����ȫ����ֵΪ0
	static const Matrix3 ZERO;
	/// ��λ���󣬶Խ�������ֵΪ1������Ϊ0
	static const Matrix3 IDENTITY;
};

#if _MSC_VER >= 1200
#pragma warning(pop)
#endif

/// \brief ��ά��������ά������ˣ��õ�һ����������
template <typename Real>
Vector3<Real> operator * (const Vector3<Real>& v, const Matrix3<Real>& mtx);

/// \brief ������������
template <typename Real>
Matrix3<Real> operator* (Real fScalar, const Matrix3<Real>& mtx);

/// float���͵���ά����
typedef Matrix3<float> Matrix3f;

/// double���͵���ά����
typedef Matrix3<double> Matrix3d;

#include "Matrix3.inl"

NS_KMATH_END

#endif // __MATRIX3_H__
