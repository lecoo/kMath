/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Vector.h
/// \brief   Nά����,�����Ļ���
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
/// \brief Nά����,�����Ļ���
/////////////////////////////////////////////////////////////////////////////////
template <int N, class Real>
class Vector
{
public:
    /// \brief ����������i����
    Real operator[] (int i) const;
    /// \brief ����������i����������
    Real& operator[] (int i);
    /// \brief ������ֵ����
    Vector& operator= (const Vector& rkV);
    /// \brief �����е�����
    bool operator== (const Vector& rkV) const;
    /// \brief �����в�������
    bool operator!= (const Vector& rkV) const;
    /// \brief ����С������(����ѧ����, �Ա�֧����Ϊmap,set�Ĺؼ���)
    bool operator<  (const Vector& rkV) const;

    // arithmetic updates
    /// \brief ������+=����
    Vector& operator+= (const Vector& rkV);
    /// \brief ������-=����
    Vector& operator-= (const Vector& rkV);
    /// \brief ����������*=����
    Vector& operator*= (Real fScalar);
    /// \brief ����������/=����
    Vector& operator/= (Real fScalar);

    // vector operations
    /// \brief ����ȡģ����
    Real Length () const;
    /// \brief ȡ����ģ��ƽ��
    Real LengthSquared () const;
    /// \brief �����������
    Real Dot (const Vector& rkV) const;
    /// \brief ��׼����ǰ����
    void Normalize ();

protected:
    /// \brief �޲������캯����������ʼ��
    Vector ();
    /// \brief �������캯��
    Vector (const Vector& rkV);
    /// \brief ���캯������p[N]��ʼ��m_afTuple[N]
    explicit Vector (const Real* p);
};

#include "Vector.inl"

/// \brief �����������
template <int N, class Real>
Real Dot(const Vector<N,Real>& v1, const Vector<N,Real>& v2);

template <int N, class Real>
Real Length(const Vector<N,Real>& v);

NS_KMATH_END

#endif //__VECTOR_H__
