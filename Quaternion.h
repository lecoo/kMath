/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Quaternion.h
/// \brief   ��Ԫ��
/// \author  DaiMingzhuang
/// \version 
/// \date    2009-03
/////////////////////////////////////////////////////////////////////////////////
#ifndef __QUATERNION_H__
#define __QUATERNION_H__

#include "MathDef.h"

NS_KMATH_BEGIN

/////////////////////////////////////////////////////////////////////////////////
/// \class Quaternion
/// \brief ��Ԫ��
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
class Quaternion
{
public:
	/// \brief ���캯�������������ݳ�ʼ��Ϊ0��
	Quaternion(void) {}
	/// \brief ���캯��������ÿ�����ݵ�ֵ��
	Quaternion(Real _x, Real _y, Real _z, Real _w) { x = _x; y = _y; z = _z; w = _w; }
	/// \brief ���캯������p[4]��ʼ�����ݡ�
	explicit Quaternion(const Real* p) { x = p[0]; y = p[1]; z = p[2]; w = p[3]; }

	/// \brief ����ÿ�����ݵ�ֵ��
	Quaternion& Set(Real _x, Real _y, Real _z, Real _w) { x = _x; y = _y; z = _z; w = _w; return *this; }
	/// \brief ��p[4]��ʼ�����ݡ�
	Quaternion& Set(const Real* p) { x = p[0]; y = p[1]; z = p[2]; w = p[3]; return *this; }

	/// \brief ����������
	Quaternion& operator = (const Quaternion& qtn) { x = qtn.x; y = qtn.y; z = qtn.z; w = qtn.w; return *this; }
	/// \brief �Ƚ�������Ԫ���Ƿ���ȡ�
	bool operator == (const Quaternion& qtn) const { return (w == qtn.w && x == qtn.x && y == qtn.y && z == qtn.z); }
	/// \brief �Ƚ�������Ԫ���Ƿ񲻵ȡ�
	bool operator != (const Quaternion& qtn) const { return !(*this==qtn); }

	/// \brief ��Ԫ���ӷ�������Ӧλ����ӡ�
	Quaternion operator + (const Quaternion& qtn) const { return Quaternion(x+qtn.x, y+qtn.y, z+qtn.z, w+qtn.w); }
	/// \brief ��Ԫ������������Ӧλ�������
	Quaternion operator - (const Quaternion& qtn) const { return Quaternion(x-qtn.x, y-qtn.y, z-qtn.z, w-qtn.w); }
	/// \brief ��Ԫ���˷���
	Quaternion operator * (const Quaternion& qtn) const;
	/// \brief ��Ԫ��������ֵ���������ݳ���s��
	Quaternion operator * (Real s) const { return Quaternion(x*s, y*s, z*s, w*s); }
	/// \brief ��Ԫ��������ֵ���������ݳ���s��
	Quaternion operator / (Real s) const { s = 1.0f / s; return Quaternion(x*s, y*s, z*s, w*s); }
	/// \brief ��Ԫ���󷴣���������ȡ�෴ֵ��
	Quaternion operator - (void) const { return Quaternion(-x, -y, -z, -w); }

	/// \brief ��Ԫ���ӷ�������Ӧλ����ӣ���������Լ���
	Quaternion& operator += (const Quaternion& qtn) { x += qtn.x; y += qtn.y; z += qtn.z; w += qtn.w; return *this; }
	/// \brief ��Ԫ������������Ӧλ���������������Լ���
	Quaternion& operator -= (const Quaternion& qtn) { x -= qtn.x; y -= qtn.y; z -= qtn.z; w -= qtn.w; return *this; }
	/// \brief ��Ԫ���˷�����������Լ���
	Quaternion& operator *= (const Quaternion& qtn) { return *this = Mul(qtn); }
	/// \brief ��Ԫ��������ֵ���������ݳ���s����������Լ���
	Quaternion& operator *= (Real s) { x *= s; y *= s; z *= s; w *= s; return *this; }
	/// \brief ��Ԫ��������ֵ���������ݳ���s����������Լ���
	Quaternion& operator /= (Real s) { return (*this) *= (1.0f / s); }

	/// \brief ���õ�λ��Ԫ����
	Quaternion& SetIdentity(void);

    /// \brief ����Ԫ���Ĺ�����Ԫ����
    Quaternion Conjugate() const;
	/// \brief ����Ԫ�����档
	Quaternion Inverse(void) const;

	Quaternion<Real> Slerp(const Quaternion<Real>& qtn2, Real lambda);

	/// \brief ����ת����ת������Ԫ����
	void FromRotationMatrix(const Matrix3<Real>& mtx);
	/// \brief ����Ԫ��ת������ת����
	void ToRotationMatrix(Matrix3<Real>& mtx) const;
	/// \brief ����ת���ת��ת������Ԫ����
	void FromAxisAngle(const Vector3<Real>& axis, Real ang);
	/// \brief ����Ԫ��ת������ת���ת�ǡ�
	bool ToAxisAngle(Vector3<Real>& axis, Real& ang) const;

	Vector3<Real> GetAxis() const;
	Real GetAngle() const;

	/// \brief ��ŷ����ת������Ԫ������ZXY�Ĵ���
	void FromEulerAnglesZXY(const Vector3<Real>& angles);

	/// \brief ����Ԫ��ת����ŷ���ǣ���ZXY�Ĵ���
	Vector3<Real> GetEulerZXY() const;

	Real Dot(const Quaternion& qtn) const;

	Real Length() const;

public:
	// �൱��axis.x * sin(ang/2)
	Real x;
	// �൱��axis.y * sin(ang/2)
	Real y;
	// �൱��axis.z * sin(ang/2)
	Real z;
	// �൱��cos(ang/2)
	Real w;

public:
	/// ��λ��Ԫ����wΪ1��x��y��zΪ0
	static const Quaternion IDENTITY;
};
/// float���͵���Ԫ��
typedef Quaternion<float> Quaternionf;

/// double���͵���Ԫ��
typedef Quaternion<double> Quaterniond;

/// \brief �������Բ�ֵ
template<typename Real>
Quaternion<Real> Slerp(const Quaternion<Real>& qtn1, const Quaternion<Real>& qtn2, Real lambda);

template<typename Real>
Real Dot(const Quaternion<Real>& q1, const Quaternion<Real>& q2);

template<typename Real>
Real Length(const Quaternion<Real>& q);

template<typename Real>
Quaternion<Real> Normalize(const Quaternion<Real>& q);

#include "Quaternion.inl"

NS_KMATH_END

#endif // __QUATERNION_H__
