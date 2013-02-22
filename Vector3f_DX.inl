/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Vector3f_DX.inl
/// \brief   Vector3fµÄDirectXÌØ»¯°æ
/// \author  DaiMingZhuang
/// \version 0.9
/// \date    2009-07
/// \history 
/////////////////////////////////////////////////////////////////////////////////
template <> inline
float Vector3<float>::Dot(const Vector3<float>& v) const
{
    return D3DXVec3Dot( (CONST D3DXVECTOR3*)this, (CONST D3DXVECTOR3 *)&v);
}

template <> inline
Vector3<float> Vector3<float>::Cross(const Vector3<float>& v) const
{
    Vector3<float> res;
    D3DXVec3Cross( (D3DXVECTOR3*)&res, (CONST D3DXVECTOR3*)this, (CONST D3DXVECTOR3 *)&v);
    return res;
}
