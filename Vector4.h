/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Vector4.h
/// \brief   4ά����
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
/// \brief ��ά����(����)
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
class Vector4 : public Vector<4,Real>
{
public:
	/// \brief ���캯������x��y��z��w��ʼ��Ϊ0��
    Vector4(void);
    /// \brief �������캯��
    Vector4(const Vector4& v);
	/// \brief ���캯��������x��y��z��wΪxyzw��
	explicit Vector4(Real xyzw);
	/// \brief ���캯�����ֱ�����x��y��z��w��
	Vector4(Real _x, Real _y, Real _z, Real _w);
    /// \brief ���캯�����û��������й���
    Vector4(const Vector<4,Real>& v);
    /// \brief ��ֵ���㣬�û��������и�ֵ
    Vector4& operator= (const Vector<4,Real>& v);
    /// \brief ��ֵ����
    Vector4& operator= (const Vector4& v);

    // arithmetic operations
    /// \brief �����ӷ�����
    Vector4 operator+ (const Vector4& rkV) const;
    /// \brief ������������
    Vector4 operator- (const Vector4& rkV) const;
    /// \brief ������������
    Vector4 operator* (Real fScalar) const;
    /// \brief ����������(��1/fScalar����������)
    Vector4 operator/ (Real fScalar) const;
	/// \brief ����ȡ������
	Vector4 operator+ () const;
    /// \brief ����ȡ������
    Vector4 operator- () const;

    /// \brief ת��ΪVector3�������ݣ�����w����
    Vector3<Real> ToVector3() const;

    /// \brief ����ǰ������λ�
    void Homogenize();

public:
    /// X���ϵ�����
    Real x;
    /// Y���ϵ�����
    Real y;
    /// Z���ϵ�����
    Real z;
    /// W���ϵ�����
    Real w;

public:
	/// ��������ȫ����ֵΪ0
	static const Vector4 ZERO;
	/// һ������ȫ����ֵΪ1
	static const Vector4 ONE;
	/// X������Ϊ1��Y��Z��W������Ϊ0
	static const Vector4 UNIT_X;
	/// Y������Ϊ1��X��Z��W������Ϊ0
	static const Vector4 UNIT_Y;
	/// Z������Ϊ1��X��Y��W������Ϊ0
	static const Vector4 UNIT_Z;
	/// W������Ϊ1��X��Y��Z������Ϊ0
	static const Vector4 UNIT_W;
};

/// float���͵���ά����
typedef Vector4<float> Vector4f;

/// double���͵���ά����
typedef Vector4<double> Vector4d;

/// \brief ���ز��������ı�׼�����������ı��������
template <typename Real>
Vector4<Real> Normalize (const Vector4<Real>& rkV);

/// \brief ���ز�����������λ����������ı��������
template <typename Real>
Vector4<Real> Homogenize (const Vector4<Real>& rkV);

/// \brief ������������
template <typename Real>
Vector4<Real> operator* (float fScalar, const Vector4<Real>& rkV);

#include "Vector4.inl"

NS_KMATH_END

#endif // __VECTOR4_H__
