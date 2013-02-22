/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    MathUtil.h
/// \brief   数学常用工具
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
    static const float FLOAT_EPS = 1.192092896e-07F; //在IEEE754标准下1.0f和比1.0f大的最小的float型数相差的值
#endif

#ifdef FLT_MAX
    static const float FLOAT_INFINITE = FLT_MAX;
#else
    static const float FLOAT_INFINITE = 3.402823466e+38F; //在IEEE754标准下float能表示的最大值
#endif

    /// \brief 返回a除以b四舍五入的结果
    inline int DivRounded(int a, int b)
    {
        return (a + (b>>1)) / b;
    }

    /// \brief 返回a除以b向上取整的结果
    inline int DivRoundedUp(int a, int b)
    {
        return (a + b - 1) / b;
    }

    /// \brief v mod m, m!=0
    /// \returns 若m>0, 返回v在[0,m)区间上的等价类代表元
    ///          若m<0, 返回v在(m,0]区间上的等价类代表元
    inline int Mod(int v, int m)
    {
        return (v%m + m)%m;
    }

    /// \brief 返回v模m，v,m应为实数，m>0
    template <typename T>
    inline T FMod(T v, T m)
    {
		return std::fmod(v, m);
    }

	/// \brief 取一个小于等于v的最大整数
	template <typename T>
	inline T Floor(T v)
	{
		return std::floor(v);
	}

	/// \brief 取一个大于等于v的最小整数
	template <typename T>
	inline T Ceil(T v)
	{
		return std::ceil(v);
	}

	/// \brief 四舍五入v到最近整数
	template <typename T>
	inline T Round(T v)
	{
		return Floor(v+T(0.5f));
	}

	/// \brief 取一个浮点数的小数部分
	///        恒有 v = Floor(v) + Frac(v)
	template <typename T>
	inline T Frac(T v)
	{
		return v - Floor(v);
	}

    /// \brief 求两个数之间的最大值
    template <typename T>
    inline T Max(T t1, T t2)
    {
        return (t1 > t2)? t1 : t2;
    }

    /// \brief 求三个数之间的最大值
    template <typename T>
    inline T Max(T v1, T v2, T v3)
    {
        T tmp = (v1 > v2)? v1 : v2;
        tmp = (tmp > v3)? tmp : v3;
        return tmp;
    }

    /// \brief 求四个数之间的最大值
    template <typename T>
    T Max(T v1, T v2, T v3, T v4)
	{
        T tmp = (v1 > v2)? v1 : v2;
        tmp = (tmp > v3)? tmp : v3;
        tmp = (tmp > v4)? tmp : v4;
        return tmp;
    }

    /// \brief 求两个数之间的最小值
    template <typename T>
    inline T Min(T t1, T t2)
    {
        return (t1 < t2)? t1 : t2;
    }

    /// \brief 求三个数之间的最小值
    template <typename T>
    inline T Min(T v1, T v2, T v3)
    {
        T tmp = (v1 < v2)? v1 : v2;
        tmp = (tmp < v3)? tmp : v3;
        return tmp;
    }

    /// \brief 求四个数之间的最小值
    template <typename T>
    T Min(T v1, T v2, T v3, T v4)
    {
        T tmp = (v1 < v2)? v1 : v2;
        tmp = (tmp < v3)? tmp : v3;
        tmp = (tmp < v4)? tmp : v4;
        return tmp;
    }

    /// \brief 符号函数，大于0时返回1，小于0是返回-1，等于0是返回0。
    template <typename T>
    inline int Sign(const T& t)
    {
        return (t > 0)? 1 : (t < 0)? -1 : 0;
    }

    /// \brief 求数的绝对值，支持模板。
    template <typename T>
    inline T Abs(const T& t)
    {
        return (t < 0)? -t : t;
    }

    /// \brief 求平方
    template <typename T>
    inline T Square(T t)
    {
        return t * t;
    }

    /// \brief 求立方
    template <typename T>
    inline T Cube(T t)
    {
        return t * t * t;
    }

    /// \brief 求平方根
    template <typename T>
    inline T Sqrt(T t)
    {
        return std::sqrt(t);
    }

	/// \brief 求a的b次幂
	template <typename T>
	inline T Pow(T a, T b)
	{
		return std::pow(a, b);
	}

    /// \brief 夹齐
    template <typename T>
    inline T Clamp(T n, T a, T b)
    {
        if (a < b)
            return (n < a)? a : (n > b)? b : n;
        else
            return (n < b)? b : (n > a)? a : n;
    }

    /// \brief 把x限制在[0,1]
    template<typename T>
    inline T Saturate(const T& x)
    {
        return (x < T(1.0f) ? (x > T(0.0f) ? x : T(0.0f)) : T(1.0f));
    }

    /// \brief 浮点判等，eps为绝对误差限
    template <typename T>
    inline bool Eq(T x, T y, T eps = (T)1e-5)
    {
        return (x-y<=eps && y-x<=eps);
    }

    /// \brief 浮点判等，eps为相对误差限(参考x,y中的较大者，并对极小值做过特判)
    /// 如eps = 1e-5 则当x与y的差值小于等于x,y中较大者的绝对值的10万分之一(1e-5)时返回true
    template <typename T>
    inline bool Eqr(T x, T y, T eps = (T)1e-5)
    {
        return (Abs(x-y) <= eps * Max(Abs(x),Abs(y)));
        //return (Abs(x-y) <= eps * (Abs(x)+Abs(y)+(T)1.0));
    }

	/// \brief 正弦
	template <typename T>
	inline T Sin(T x)
	{
		return sin(x);
	}

	/// \brief 余弦
	template <typename T>
	inline T Cos(T x)
	{
		return cos(x);
	}

	/// \brief 正切
	template <typename T>
	inline T Tan(T x)
	{
		return tan(x);
	}

	/// \brief 反正弦
	template <typename T>
	inline T Asin(T x)
	{
		x = Clamp(x, T(-1.0f), T(1.0f));
		return asin(x);
	}

	/// \brief 反余弦
	template <typename T>
	inline T Acos(T x)
	{
		x = Clamp(x, T(-1.0f), T(1.0f));
		return acos(x);
	}

	/// \brief 反正切
	template <typename T>
	inline T Atan(T x)
	{
		return atan(x);
	}

	/// \brief (y/x) 的反正切
	template <typename T>
	inline T Atan2(T x, T y)
	{
		return atan2(x, y);
	}

    /// \brief 把从纹理取出的数据解压，其实就是把[0,1]映射到[-1,1]
    template<typename T>
    inline T Expand(const T& v)
    {
        return T(2.0f)*(v-T(0.5f));
    }

    /// \brief 线性插值
    template<typename T>
    inline T Lerp(const T& x, const T& y, float lambda)
    {
        return x*(1.0f-lambda) + y*lambda;
    }

    /// \brief 平滑插值
    template<typename T>
    T Smorp(const T& x, const T& y, float lambda)
    {
        lambda = Saturate(lambda);
        float lambda2 = lambda * lambda;
        float lambda3 = lambda2 * lambda;
        float t = -2.0f * lambda3 + 3.0f * lambda2;
        return x * (1-t) + y * t;
    }

	/// \brief 通用向量球面线性插值
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

	/// \brief 球面线性插值(向量)
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

	/// \brief 通用四元数球面线性插值
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

	/// \brief 球面线性插值(四元数)
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

    /// \brief 三角形插值
    template<typename T>
    inline T Terp(const T& v1, const T& v2, const T& v3, float s, float t)
    {
        return v1 * s + v2 * t + v3 * (1.0f-s-t);
    }

	/// \brief 三次Hermite插值
	/// \param[in] pa 起点位置
	/// \param[in] va 起点速度
	/// \param[in] pb 终点位置
	/// \param[in] vb 终点速度
	/// \param[in] lambda 
	template<typename T>
	inline T Herp(const T& pa, const T& va, const T& pb, const T& vb, float lambda)
	{
		const float lambda2 = lambda * lambda;
		const float lambda3 = lambda2 * lambda;
		return pa * (lambda3 * 2.0f - lambda2 * 3.0f + 1.0f) + va * (lambda3 - lambda2 * 2.0f + lambda) + vb * (lambda3 - lambda2) + pb * (lambda3 * -2.0f + lambda2 * 3.0f);
	}

    /// \brief 计算反射向量
    /// \param[in] I eye指向objpos的视向量
    /// \param[in] N objpos处的法向量
    template<typename T>
    inline T Reflect(const T& I, const T& N)
    {
        return I - N * (2.0f * Dot(I,N));
    }

    /// \brief 计算折射向量
    /// \param[in] I eye指向objpos的视向量
    /// \param[in] N objpos处的法向量
    /// \param[in] etaRatio 入射介质折射率和出射介质折射率的比值
    template<typename T>
    inline T Refract(const T& I, const T& N, float etaRatio)
    {
        float cosI = Dot(-I,N);
        float cosT2 = 1- (1 - cosI*cosI) * etaRatio*etaRatio;
        float cosT = Sqrt(Abs(cosT2));
        return cosT2 > 0.0f ? -N*cosT + (I + N*cosI) * etaRatio : T(0.0f);
    }

	/// \brief 拐向判断2D
    // 如果Orient2D(A,B,P)=0，则三点共线
    // 如果Orient2D(A,B,P)>0，则在左手系下A-B-P是左拐(在右手系下A-B-P是右拐)
    // 如果Orient2D(A,B,P)<0，则在左手系下A-B-P是右拐(在右手系下A-B-P是左拐)
    template <typename V2>
    inline float Orient2D(const V2& A, const V2& B, const V2& P)
    {
        return PseudoCross(A-P, B-P);
    }

	/// \brief 拐向判断3D
    // 如果Orient3D(A,B,C,P)=0，则四点共面
    // 如果Orient3D(A,B,C,P)>0，则在左手系下从P点看A-B-C为左拐
    // 如果Orient3D(A,B,C,P)<0，则在左手系下从P点看A-B-C为右拐
    template <typename V3>
    inline float Orient3D(const V3& A, const V3& B, const V3& C, const V3& P)
    {
        return Dot(A-P, Cross(B-P, C-P));
    }

	/// \brief 随机打乱
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
