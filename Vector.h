/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Vector.h
/// \brief   N维向量,向量的基类
/// \author  DaiMingZhuang
/// \version 
/// \date    2009-07
/////////////////////////////////////////////////////////////////////////////////
#ifndef __VECTOR_H__
#define __VECTOR_H__

#include "MathDef.h"

NS_KMATH_BEGIN

/////////////////////////////////////////////////////////////////////////////////
/// \class Vector
/// \brief N维向量,向量的基类
/////////////////////////////////////////////////////////////////////////////////
template <int N, class Real>
class Vector
{
public:
    /// \brief 返回向量第i分量
    Real operator[] (int i) const;
    /// \brief 返回向量第i分量的引用
    Real& operator[] (int i);
    /// \brief 向量赋值运算
    Vector& operator= (const Vector& rkV);
    /// \brief 向量判等运算
    bool operator== (const Vector& rkV) const;
    /// \brief 向量判不等运算
    bool operator!= (const Vector& rkV) const;
    /// \brief 向量小于运算(无数学意义, 以便支持作为map,set的关键字)
    bool operator<  (const Vector& rkV) const;

    // arithmetic updates
    /// \brief 向量的+=运算
    Vector& operator+= (const Vector& rkV);
    /// \brief 向量的-=运算
    Vector& operator-= (const Vector& rkV);
    /// \brief 向量与数的*=运算
    Vector& operator*= (Real fScalar);
    /// \brief 向量与数的/=运算
    Vector& operator/= (Real fScalar);

    // vector operations
    /// \brief 向量取模运算
    Real Length () const;
    /// \brief 取向量模的平方
    Real LengthSquared () const;
    /// \brief 向量点乘运算
    Real Dot (const Vector& rkV) const;
    /// \brief 标准化当前向量
    void Normalize ();

protected:
    /// \brief 无参数构造函数，不做初始化
    Vector ();
    /// \brief 拷贝构造函数
    Vector (const Vector& rkV);
    /// \brief 构造函数，用p[N]初始化m_afTuple[N]
    explicit Vector (const Real* p);
};

#include "Vector.inl"

/// \brief 向量点乘运算
template <int N, class Real>
Real Dot(const Vector<N,Real>& v1, const Vector<N,Real>& v2);

template <int N, class Real>
Real Length(const Vector<N,Real>& v);

NS_KMATH_END

#endif //__VECTOR_H__
