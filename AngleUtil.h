/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    AngleUtil.h
/// \brief   
/// \author  DaiMingZhuang
/// \version 
/// \date    2009-10
/////////////////////////////////////////////////////////////////////////////////
#ifndef __ANGLE_UTIL_H__
#define __ANGLE_UTIL_H__

#include "MathDef.h"

NS_KMATH_BEGIN

/////////////////////////////////////////////////////////////////////////////////
/// \class AngleUtil
/// \brief 
/////////////////////////////////////////////////////////////////////////////////
template <typename T>
class AngleUtil
{
public:
    /// \brief �Ƕ�ת����
    static T DegreeToRadian(T degree);
    /// \brief ����ת�Ƕ�
    static T RadianToDegree(T radian);

public:
    /// Բ���ʦ�
    static const T PI;
    /// Բ���ʦе�����
    static const T DBL_PI;
    /// Բ���ʦе�һ��
    static const T HALF_PI;
    /// Բ���ʦеĵ���
    static const T INV_PI;

    /// �Ƕȵ����ȵ�ϵ��
    static const T DEG_TO_RAD;
    /// ���ȵ��Ƕȵ�ϵ��
    static const T RAD_TO_DEG;
    
    //���ýǶ�Ӧ�Ļ���ֵ
    static const T ANGLE_1;
    static const T ANGLE_5;
    static const T ANGLE_10;
    static const T ANGLE_15;
    static const T ANGLE_30;
    static const T ANGLE_45;
    static const T ANGLE_60;
    static const T ANGLE_90;
    static const T ANGLE_120;
    static const T ANGLE_150;
    static const T ANGLE_180;
    static const T ANGLE_270;
    static const T ANGLE_360;
};

#include "AngleUtil.inl"

NS_KMATH_END

#endif //__ANGLE_UTIL_H__
