/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    MathDef.h
/// \brief   
/// \author  DaiMingzhuang
/// \version 
/// \date    2011-07
/////////////////////////////////////////////////////////////////////////////////
#ifndef __MATH_DEF_H__
#define __MATH_DEF_H__

#define NS_KMATH        K::MATH
#define NS_KMATH_BEGIN  namespace K{ namespace MATH{
#define NS_KMATH_END    }}
#define NS_KMATH_USING  using namespace NS_KMATH;

#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <cstring>

#ifndef K_ASSERT
#include <cassert>
#define K_ASSERT(P) assert(P)
#endif

NS_KMATH_BEGIN
template <typename Real>
class Vector2;
template <typename Real>
class Vector3;
template <typename Real>
class Vector4;
template <typename Real>
class Matrix2;
template <typename Real>
class Matrix3;
template <typename Real>
class Matrix4;
template <typename Real>
class Quaternion;
NS_KMATH_END

// #define USE_DX_MATH
#ifdef USE_DX_MATH
#include <d3dx9math.h>
#endif

#define K_MATH_LOG(format,...)

#endif // __MATH_DEF_H__
