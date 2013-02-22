/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Vector2.h
/// \brief   2ά����
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
/// \brief ��ά����(����)
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
class Vector2 : public Vector<2,Real>
{
public:
	/// \brief ���캯������x��y��ʼ��Ϊ0��
	Vector2();
    /// \brief �������캯��
    Vector2(const Vector2& v);
	/// \brief ���캯��������x��yΪxy��
	explicit Vector2(Real xy);
	/// \brief ���캯�����ֱ�����x��y��
	Vector2(Real _x, Real _y);
    /// \brief ���캯�����û��������й���
    Vector2(const Vector<2,Real>& v);
    /// \brief ��ֵ���㣬�û��������и�ֵ
    Vector2& operator= (const Vector<2,Real>& v);
    /// \brief ��ֵ����
    Vector2& operator= (const Vector2& v);

    // arithmetic operations
    /// \brief �����ӷ�����
    Vector2 operator+ (const Vector2& rkV) const;
    /// \brief ������������
    Vector2 operator- (const Vector2& rkV) const;
    /// \brief ������������
    Vector2 operator* (Real fScalar) const;
    /// \brief ����������(��1/fScalar����������)
    Vector2 operator/ (Real fScalar) const;
	/// \brief ����ȡ������
	Vector2 operator+ () const;
    /// \brief ����ȡ������
    Vector2 operator- () const;

	/// \brief ���ά�����Ĳ��(���)��
	Real Cross(const Vector2& v) const;

    /// \brief ת��Ϊ3ά����
    Vector3<Real> ToVector3(Real z = 1.0f)const;
    
public:
    /// X���ϵ�����
    Real x;
    /// Y���ϵ�����
    Real y;

public:
	/// ��������ȫ����ֵΪ0
	static const Vector2 ZERO;
	/// һ������ȫ����ֵΪ1
	static const Vector2 ONE;
	/// X������Ϊ1��Y������Ϊ0
	static const Vector2 UNIT_X;
	/// X������Ϊ0��Y������Ϊ1
	static const Vector2 UNIT_Y;
};

/// float���͵Ķ�ά����
typedef Vector2<float> Vector2f;

/// double���͵Ķ�ά����
typedef Vector2<double> Vector2d;

/// \brief ���ز��������ı�׼�����������ı��������
template <typename Real>
Vector2<Real> Normalize (const Vector2<Real>& rkV);

/// \brief ������������
template <typename Real>
Vector2<Real> operator* (Real fScalar, const Vector2<Real>& rkV);

/// \brief ��ά������������
template <typename Real>
Real PseudoCross(const Vector2<Real>& v1, const Vector2<Real>& v2);

#include "Vector2.inl"

NS_KMATH_END

#endif // __VECTOR2_H__
