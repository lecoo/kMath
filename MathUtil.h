/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    MathUtil.h
/// \brief   ��ѧ���ù���
/// \author  DaiMingzhuang
/// \version 
/// \date    2009-02
/////////////////////////////////////////////////////////////////////////////////
#ifndef	__MATH_UTIL_H__
#define	__MATH_UTIL_H__

#include "MathDef.h"
#include "AngleUtil.h"
#include "Vector2.h"
#include "Vector3.h"

NS_KMATH_BEGIN

namespace MathUtil
{
#ifdef FLT_EPSILON
    static const float FLOAT_EPS = FLT_EPSILON;
#else
    static const float FLOAT_EPS = 1.192092896e-07F; //��IEEE754��׼��1.0f�ͱ�1.0f�����С��float��������ֵ
#endif

#ifdef FLT_MAX
    static const float FLOAT_INFINITE = FLT_MAX;
#else
    static const float FLOAT_INFINITE = 3.402823466e+38F; //��IEEE754��׼��float�ܱ�ʾ�����ֵ
#endif

    /// \brief ����a����b��������Ľ��
    inline int DivRounded(int a, int b)
    {
        return (a + (b>>1)) / b;
    }

    /// \brief ����a����b����ȡ���Ľ��
    inline int DivRoundedUp(int a, int b)
    {
        return (a + b - 1) / b;
    }

    /// \brief v mod m, m!=0
    /// \returns ��m>0, ����v��[0,m)�����ϵĵȼ������Ԫ
    ///          ��m<0, ����v��(m,0]�����ϵĵȼ������Ԫ
    inline int Mod(int v, int m)
    {
        return (v%m + m)%m;
    }

    /// \brief ����vģm��v,mӦΪʵ����m>0
    template <typename T>
    inline T FMod(T v, T m)
    {
		return std::fmod(v, m);
    }

	/// \brief ȡһ��С�ڵ���v���������
	template <typename T>
	inline T Floor(T v)
	{
		return std::floor(v);
	}

	/// \brief ȡһ�����ڵ���v����С����
	template <typename T>
	inline T Ceil(T v)
	{
		return std::ceil(v);
	}

	/// \brief ��������v���������
	template <typename T>
	inline T Round(T v)
	{
		return Floor(v+T(0.5f));
	}

	/// \brief ȡһ����������С������
	///        ���� v = Floor(v) + Frac(v)
	template <typename T>
	inline T Frac(T v)
	{
		return v - Floor(v);
	}

    /// \brief ��������֮������ֵ
    template <typename T>
    inline T Max(T t1, T t2)
    {
        return (t1 > t2)? t1 : t2;
    }

    /// \brief ��������֮������ֵ
    template <typename T>
    inline T Max(T v1, T v2, T v3)
    {
        T tmp = (v1 > v2)? v1 : v2;
        tmp = (tmp > v3)? tmp : v3;
        return tmp;
    }

    /// \brief ���ĸ���֮������ֵ
    template <typename T>
    T Max(T v1, T v2, T v3, T v4)
	{
        T tmp = (v1 > v2)? v1 : v2;
        tmp = (tmp > v3)? tmp : v3;
        tmp = (tmp > v4)? tmp : v4;
        return tmp;
    }

    /// \brief ��������֮�����Сֵ
    template <typename T>
    inline T Min(T t1, T t2)
    {
        return (t1 < t2)? t1 : t2;
    }

    /// \brief ��������֮�����Сֵ
    template <typename T>
    inline T Min(T v1, T v2, T v3)
    {
        T tmp = (v1 < v2)? v1 : v2;
        tmp = (tmp < v3)? tmp : v3;
        return tmp;
    }

    /// \brief ���ĸ���֮�����Сֵ
    template <typename T>
    T Min(T v1, T v2, T v3, T v4)
    {
        T tmp = (v1 < v2)? v1 : v2;
        tmp = (tmp < v3)? tmp : v3;
        tmp = (tmp < v4)? tmp : v4;
        return tmp;
    }

    /// \brief ���ź���������0ʱ����1��С��0�Ƿ���-1������0�Ƿ���0��
    template <typename T>
    inline int Sign(const T& t)
    {
        return (t > 0)? 1 : (t < 0)? -1 : 0;
    }

    /// \brief �����ľ���ֵ��֧��ģ�塣
    template <typename T>
    inline T Abs(const T& t)
    {
        return (t < 0)? -t : t;
    }

    /// \brief ��ƽ��
    template <typename T>
    inline T Square(T t)
    {
        return t * t;
    }

    /// \brief ������
    template <typename T>
    inline T Cube(T t)
    {
        return t * t * t;
    }

    /// \brief ��ƽ����
    template <typename T>
    inline T Sqrt(T t)
    {
        return std::sqrt(t);
    }

	/// \brief ��a��b����
	template <typename T>
	inline T Pow(T a, T b)
	{
		return std::pow(a, b);
	}

    /// \brief ����
    template <typename T>
    inline T Clamp(T n, T a, T b)
    {
        if (a < b)
            return (n < a)? a : (n > b)? b : n;
        else
            return (n < b)? b : (n > a)? a : n;
    }

    /// \brief ��x������[0,1]
    template<typename T>
    inline T Saturate(const T& x)
    {
        return (x < T(1.0f) ? (x > T(0.0f) ? x : T(0.0f)) : T(1.0f));
    }

    /// \brief �����еȣ�epsΪ���������
    template <typename T>
    inline bool Eq(T x, T y, T eps = (T)1e-5)
    {
        return (x-y<=eps && y-x<=eps);
    }

    /// \brief �����еȣ�epsΪ��������(�ο�x,y�еĽϴ��ߣ����Լ�Сֵ��������)
    /// ��eps = 1e-5 ��x��y�Ĳ�ֵС�ڵ���x,y�нϴ��ߵľ���ֵ��10���֮һ(1e-5)ʱ����true
    template <typename T>
    inline bool Eqr(T x, T y, T eps = (T)1e-5)
    {
        return (Abs(x-y) <= eps * Max(Abs(x),Abs(y)));
        //return (Abs(x-y) <= eps * (Abs(x)+Abs(y)+(T)1.0));
    }

	/// \brief ����
	template <typename T>
	inline T Sin(T x)
	{
		return sin(x);
	}

	/// \brief ����
	template <typename T>
	inline T Cos(T x)
	{
		return cos(x);
	}

	/// \brief ����
	template <typename T>
	inline T Tan(T x)
	{
		return tan(x);
	}

	/// \brief ������
	template <typename T>
	inline T Asin(T x)
	{
		x = Clamp(x, T(-1.0f), T(1.0f));
		return asin(x);
	}

	/// \brief ������
	template <typename T>
	inline T Acos(T x)
	{
		x = Clamp(x, T(-1.0f), T(1.0f));
		return acos(x);
	}

	/// \brief ������
	template <typename T>
	inline T Atan(T x)
	{
		return atan(x);
	}

	/// \brief (y/x) �ķ�����
	template <typename T>
	inline T Atan2(T x, T y)
	{
		return atan2(x, y);
	}

    /// \brief �Ѵ�����ȡ�������ݽ�ѹ����ʵ���ǰ�[0,1]ӳ�䵽[-1,1]
    template<typename T>
    inline T Expand(const T& v)
    {
        return T(2.0f)*(v-T(0.5f));
    }

    /// \brief ���Բ�ֵ
    template<typename T>
    inline T Lerp(const T& x, const T& y, float lambda)
    {
        return x*(1.0f-lambda) + y*lambda;
    }

    /// \brief ƽ����ֵ
    template<typename T>
    T Smorp(const T& x, const T& y, float lambda)
    {
        lambda = Saturate(lambda);
        float lambda2 = lambda * lambda;
        float lambda3 = lambda2 * lambda;
        float t = -2.0f * lambda3 + 3.0f * lambda2;
        return x * (1-t) + y * t;
    }

	/// \brief ͨ�������������Բ�ֵ
	template<typename V>
	V SlerpV(const V& x, const V& y, float lambda)
	{
		float s = Length(x) * Length(y);
		if(s < FLOAT_EPS)
			return x;

		float fCos = Dot(x, y) / s;
		if (Abs(fCos) >= 1.0f)
			return x;

		float fAngle = Acos(fCos);
		float fInvSin = 1.0f / Sin(fAngle);
		float fCoeff0 = Sin((1.0f - lambda) * fAngle);
		float fCoeff1 = Sin(lambda * fAngle);
		return (x * fCoeff0 + y * fCoeff1) * fInvSin;
	}

	/// \brief �������Բ�ֵ(����)
	template<typename Real>
	Vector3<Real> Slerp(const Vector3<Real>& x, const Vector3<Real>& y, Real lambda)
	{
		Real s = Length(x) * Length(y);
		if(s < FLOAT_EPS)
			return x;

		Real fCos = Dot(x, y) / s;
		if (Abs(fCos) >= 1.0f)
			return x;

		Real fAngle = Acos(fCos);
		Real fInvSin = 1.0f / Sin(fAngle);
		Real fCoeff0 = Sin((1.0f - lambda) * fAngle);
		Real fCoeff1 = Sin(lambda * fAngle);
		return (x * fCoeff0 + y * fCoeff1) * fInvSin;
	}

	/// \brief ͨ����Ԫ���������Բ�ֵ
	template<typename Q>
	Q SlerpQ(const Q& x, const Q& y, float lambda)
	{
		float s = Length(x) * Length(y);
		if(s < FLOAT_EPS)
			return x;

		float fCos = Dot(x, y) / s;
		if (Abs(fCos) >= 1.0f)
			return x;

		float sign = (fCos < 0) ? -1.0f : 1.0f;
		float fAngle = Acos(sign * fCos);
		float fInvSin = 1.0f / Sin(fAngle);
		float fCoeff0 = Sin((1.0f - lambda) * fAngle);
		float fCoeff1 = Sin(sign * lambda * fAngle);
		return (x * fCoeff0 + y * fCoeff1) * fInvSin;
	}

	/// \brief �������Բ�ֵ(��Ԫ��)
	template<typename Real>
	Quaternion<Real> Slerp(const Quaternion<Real>& x, const Quaternion<Real>& y, Real lambda)
	{
		Real s = Length(x) * Length(y);
		if(s < FLOAT_EPS)
			return x;

		Real fCos = Dot(x, y) / s;
		if (Abs(fCos) >= 1.0f)
			return x;

		Real sign = (fCos < 0) ? -1.0f : 1.0f;
		Real fAngle = Acos(sign * fCos);
		Real fInvSin = 1.0f / Sin(fAngle);
		Real fCoeff0 = Sin((1.0f - lambda) * fAngle);
		Real fCoeff1 = Sin(sign * lambda * fAngle);
		return (x * fCoeff0 + y * fCoeff1) * fInvSin;
	}

    /// \brief �����β�ֵ
    template<typename T>
    inline T Terp(const T& v1, const T& v2, const T& v3, float s, float t)
    {
        return v1 * s + v2 * t + v3 * (1.0f-s-t);
    }

	/// \brief ����Hermite��ֵ
	/// \param[in] pa ���λ��
	/// \param[in] va ����ٶ�
	/// \param[in] pb �յ�λ��
	/// \param[in] vb �յ��ٶ�
	/// \param[in] lambda 
	template<typename T>
	inline T Herp(const T& pa, const T& va, const T& pb, const T& vb, float lambda)
	{
		const float lambda2 = lambda * lambda;
		const float lambda3 = lambda2 * lambda;
		return pa * (lambda3 * 2.0f - lambda2 * 3.0f + 1.0f) + va * (lambda3 - lambda2 * 2.0f + lambda) + vb * (lambda3 - lambda2) + pb * (lambda3 * -2.0f + lambda2 * 3.0f);
	}

    /// \brief ���㷴������
    /// \param[in] I eyeָ��objpos��������
    /// \param[in] N objpos���ķ�����
    template<typename T>
    inline T Reflect(const T& I, const T& N)
    {
        return I - N * (2.0f * Dot(I,N));
    }

    /// \brief ������������
    /// \param[in] I eyeָ��objpos��������
    /// \param[in] N objpos���ķ�����
    /// \param[in] etaRatio ������������ʺͳ�����������ʵı�ֵ
    template<typename T>
    inline T Refract(const T& I, const T& N, float etaRatio)
    {
        float cosI = Dot(-I,N);
        float cosT2 = 1- (1 - cosI*cosI) * etaRatio*etaRatio;
        float cosT = Sqrt(Abs(cosT2));
        return cosT2 > 0.0f ? -N*cosT + (I + N*cosI) * etaRatio : T(0.0f);
    }

	/// \brief �����ж�2D
    // ���Orient2D(A,B,P)=0�������㹲��
    // ���Orient2D(A,B,P)>0����������ϵ��A-B-P�����(������ϵ��A-B-P���ҹ�)
    // ���Orient2D(A,B,P)<0����������ϵ��A-B-P���ҹ�(������ϵ��A-B-P�����)
    template <typename V2>
    inline float Orient2D(const V2& A, const V2& B, const V2& P)
    {
        return PseudoCross(A-P, B-P);
    }

	/// \brief �����ж�3D
    // ���Orient3D(A,B,C,P)=0�����ĵ㹲��
    // ���Orient3D(A,B,C,P)>0����������ϵ�´�P�㿴A-B-CΪ���
    // ���Orient3D(A,B,C,P)<0����������ϵ�´�P�㿴A-B-CΪ�ҹ�
    template <typename V3>
    inline float Orient3D(const V3& A, const V3& B, const V3& C, const V3& P)
    {
        return Dot(A-P, Cross(B-P, C-P));
    }

	/// \brief �������
	template <typename T>
	void Shuffle(T* b, T* e, int (*RAND)(void) = &std::rand )
	{
		int n = e-b;
		for(int i=0;i<n-1;i++)
		{
			int x = RAND()%(n-i);
			std::swap(b[x], b[n-1-i]);
		}
	}
}

NS_KMATH_END

#endif // __MATH_H__
