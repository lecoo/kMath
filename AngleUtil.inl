/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    AngleUtil.inl
/// \brief   
/// \author  DaiMingZhuang
/// \version 
/// \date    2009-10
/////////////////////////////////////////////////////////////////////////////////
template <typename T>
inline
T AngleUtil<T>::DegreeToRadian(T degree)
{
    return (T)(degree * DEG_TO_RAD);
}

template <typename T>
inline
T AngleUtil<T>::RadianToDegree(T radian)
{
    return (T)(radian * RAD_TO_DEG);
}

/// 圆周率π
template<typename T> const T AngleUtil<T>::PI         = (T)3.14159265358979323846;
/// 圆周率π的两倍
template<typename T> const T AngleUtil<T>::DBL_PI     = (T)6.28318530717958647693;
/// 圆周率π的一半
template<typename T> const T AngleUtil<T>::HALF_PI    = (T)1.57079632679489661923;
/// 圆周率π的倒数
template<typename T> const T AngleUtil<T>::INV_PI     = (T)0.31830988618379067154;

/// 角度到弧度的系数
template<typename T> const T AngleUtil<T>::DEG_TO_RAD = (T)0.01745329251994329577;
/// 弧度到角度的系数
template<typename T> const T AngleUtil<T>::RAD_TO_DEG = (T)57.2957795130823208768;

//常用角对应的弧度值
template<typename T> const T AngleUtil<T>::ANGLE_1    = (T)0.01745329251994329577;
template<typename T> const T AngleUtil<T>::ANGLE_5    = (T)0.08726646259971647885;
template<typename T> const T AngleUtil<T>::ANGLE_10   = (T)0.17453292519943295769;
template<typename T> const T AngleUtil<T>::ANGLE_15   = (T)0.26179938779914943654;
template<typename T> const T AngleUtil<T>::ANGLE_30   = (T)0.52359877559829887308;
template<typename T> const T AngleUtil<T>::ANGLE_45   = (T)0.78539816339744830962;
template<typename T> const T AngleUtil<T>::ANGLE_60   = (T)1.04719755119659774615;
template<typename T> const T AngleUtil<T>::ANGLE_90   = (T)1.57079632679489661923;
template<typename T> const T AngleUtil<T>::ANGLE_120  = (T)2.09439510239319549231;
template<typename T> const T AngleUtil<T>::ANGLE_150  = (T)2.61799387799149436539;
template<typename T> const T AngleUtil<T>::ANGLE_180  = (T)3.14159265358979323846;
template<typename T> const T AngleUtil<T>::ANGLE_270  = (T)4.71238898038468985769;
template<typename T> const T AngleUtil<T>::ANGLE_360  = (T)6.28318530717958647693;
