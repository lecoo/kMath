/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Matrix2.h
/// \brief   2��2�ľ���(�ҳ˾���)
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
/// \brief 2��2�ľ���(�ҳ˾���)
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
class Matrix2
{
public:
	/// \brief ���캯����δ��ʼ����
	Matrix2(void);
	/// \brief ���캯��������ÿ�����ݵ�ֵ��
	Matrix2(Real m00, Real m01, Real m10, Real m11);
	/// \brief �������캯����
	Matrix2(const Matrix2& mtx);

	/// \brief ����ÿ�����ݵ�ֵ��
	Matrix2& Set(Real m00, Real m01, Real m10, Real m11);

	/// \brief ��ȡ��row�У���col�е����ݡ�
	Real operator() (int row, int col) const;
	/// \brief ��ȡ��row�У���col�е����ݣ��ɸ��ġ�
	Real& operator() (int row, int col);

	/// \brief ��ȡ��row�е�������
	const Vector2<Real>& operator[] (int row) const;

	/// \brief ��ȡ��row�е�������
	Vector2<Real>& operator[] (int row);

	/// \brief ����������
	Matrix2& operator = (const Matrix2& mtx);
	/// \brief �Ƚ����������Ƿ���ȡ�
	bool operator == (const Matrix2& mtx) const;
	/// \brief �Ƚ����������Ƿ񲻵ȡ�
	bool operator != (const Matrix2& mtx) const;

	/// \brief ����ӷ�������Ӧλ����ӡ�
	Matrix2 operator + (const Matrix2& mtx) const;
	/// \brief �������������Ӧλ�������
	Matrix2 operator - (const Matrix2& mtx) const;
	/// \brief ����˷���
	Matrix2 operator * (const Matrix2& mtx) const;
	/// \brief ���������ֵ���������ݳ���s��
	Matrix2 operator * (Real s) const;
	/// \brief ���������ֵ���������ݳ���s��
	Matrix2 operator / (Real s) const;
	/// \brief �����󷴣���������ȡ�෴ֵ��
	Matrix2 operator - (void) const;

	/// \brief ����ӷ�������Ӧλ����ӣ���������Լ���
	Matrix2& operator += (const Matrix2& mtx);
	/// \brief �������������Ӧλ���������������Լ���
	Matrix2& operator -= (const Matrix2& mtx);
	/// \brief ����˷�����������Լ���
	Matrix2& operator *= (const Matrix2& mtx);
	/// \brief ���������ֵ���������ݳ���s����������Լ���
	Matrix2& operator *= (Real s);
	/// \brief ���������ֵ���������ݳ���s����������Լ���
	Matrix2& operator /= (Real s);

	/// \brief ���õ�λ����
	Matrix2& SetIdentity(void);
	/// \brief ������z�����ת����
	Matrix2& Set2DRotation(Real ang);

	/// \brief ���ά�����ת�þ���
	Matrix2 Transpose(void) const;
	/// \brief ���ά����������
	Matrix2 Inverse(void) const;
	/// \brief ���ά���������ʽ��
	Real Determinant(void) const;
public:
	union
	{
        struct
        {
            Real _11, _12;
            Real _21, _22;
        };
		/// ��[2][2]��ά������ʽ������
		Real m[2][2];
		/// ��[4]һά������ʽ������
		Real m2[4];
	};

public:
	/// �����ȫ����ֵΪ0
	static const Matrix2 ZERO;
	/// ��λ���󣬶Խ�������ֵΪ1������Ϊ0
	static const Matrix2 IDENTITY;
};
#if _MSC_VER >= 1200
#pragma warning(pop)
#endif

/// \brief ��ά�����Ͷ�ά������ˣ��õ�һ����������
template <typename Real>
Vector2<Real> operator * (const Vector2<Real>& v, const Matrix2<Real>& mtx);

/// \brief ������������
template <typename Real>
Matrix2<Real> operator* (Real fScalar, const Matrix2<Real>& mtx);

/// float���͵Ķ�ά����
typedef Matrix2<float> Matrix2f;

/// double���͵Ķ�ά����
typedef Matrix2<double> Matrix2d;

#include "Matrix2.inl"

NS_KMATH_END

#endif // __MATRIX2_H__
