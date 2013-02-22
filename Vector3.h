/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Vector3.h
/// \brief   3ά����
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
/// \brief ��ά����(����)
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
class Vector3 : public Vector<3,Real>
{
public:
	/// \brief ���캯������x��y��z��ʼ��Ϊ0��
    Vector3(void);
    /// \brief �������캯��
    Vector3(const Vector3& v);
	/// \brief ���캯��������x��y��zΪxyz��
	explicit Vector3(Real xyz);
	/// \brief ���캯�����ֱ�����x��y��z��
	Vector3(Real _x, Real _y, Real _z);
    /// \brief ���캯�����û��������й���
    Vector3(const Vector<3,Real>& v);
    /// \brief ��ֵ���㣬�û��������и�ֵ
    Vector3& operator= (const Vector<3,Real>& v);
    /// \brief ��ֵ����
    Vector3& operator= (const Vector3& v);

    // arithmetic operations
    /// \brief �����ӷ�����
    Vector3 operator+ (const Vector3& rkV) const;
    /// \brief ������������
    Vector3 operator- (const Vector3& rkV) const;
    /// \brief ������������
    Vector3 operator* (Real fScalar) const;
    /// \brief ����������(��1/fScalar����������)
    Vector3 operator/ (Real fScalar) const;
	/// \brief ����ȡ������
	Vector3 operator+ () const;
    /// \brief ����ȡ������
    Vector3 operator- () const;

    /// \brief ������
    Real Dot (const Vector3<Real>& v) const;
	/// \brief ������
	Vector3<Real> Cross(const Vector3<Real>& v) const;
	/// \brief ������
	Matrix3<Real> Tensor(const Vector3<Real>& v) const;

    /// \brief ������Ϊ(x,y,z,1)�ø���������б任�����Խ����λ�
    Vector3<Real> TransformCoord(const Matrix4<Real>& mat)const;

    /// \brief ������Ϊ(x,y,z,0)�ø���������б任
    /// ����������߱任�������ľ���Ӧ���ǶԶ�����б任�ľ����������ת�þ���
    Vector3<Real> TransformNormal(const Matrix4<Real>& mat)const;

    /// \brief ת��Ϊ4ά����
    Vector4<Real> ToVector4(Real w = 1.0f)const;

    /// \brief ת��Ϊ2ά����������z����
    Vector2<Real> ToVector2() const;

public:
    /// X���ϵ�����
    Real x;
    /// Y���ϵ�����
    Real y;
    /// Z���ϵ�����
    Real z;

public:
	/// ��������ȫ����ֵΪ0
	static const Vector3 ZERO;
	/// һ������ȫ����ֵΪ1
	static const Vector3 ONE;
	/// X������Ϊ1��Y��Z������Ϊ0
	static const Vector3 UNIT_X;
	/// Y������Ϊ1��X��Z������Ϊ0
	static const Vector3 UNIT_Y;
	/// Z������Ϊ1��X��Y������Ϊ0
	static const Vector3 UNIT_Z;
};

/// float���͵���ά����
typedef Vector3<float> Vector3f;

/// double���͵���ά����
typedef Vector3<double> Vector3d;

/// \brief ���ز��������ı�׼�����������ı��������
template <typename Real>
Vector3<Real> Normalize (const Vector3<Real>& rkV);

/// \brief ������������
template <typename Real>
Vector3<Real> operator* (Real fScalar, const Vector3<Real>& rkV);

/// \brief 3ά�����������
template <typename Real>
Vector3<Real> Cross(const Vector3<Real>& v1, const Vector3<Real>& v2);

/// \brief 3ά����������
template <typename Real>
Matrix3<Real> Tensor(const Vector3<Real>& u, const Vector3<Real>& v);

#include "Vector3.inl"

#ifdef USE_DX_MATH
#include "Vector3f_DX.inl"  //��DX�ػ�Vector3<float>��ʵ��
#endif

NS_KMATH_END

#endif // __VECTOR3_H__
