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
    /// \brief 角度转弧度
    static T DegreeToRadian(T degree);
    /// \brief 弧度转角度
    static T RadianToDegree(T radian);

public:
    /// 圆周率π
    static const T PI;
    /// 圆周率π的两倍
    static const T DBL_PI;
    /// 圆周率π的一半
    static const T HALF_PI;
    /// 圆周率π的倒数
    static const T INV_PI;

    /// 角度到弧度的系数
    static const T DEG_TO_RAD;
    /// 弧度到角度的系数
    static const T RAD_TO_DEG;
    
    //常用角对应的弧度值
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
