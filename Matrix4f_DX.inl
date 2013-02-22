/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Matrix4f_DX.inl
/// \brief   Matrix4fµÄDirectXÌØ»¯°æ
/// \author  DaiMingZhuang
/// \version 
/// \date    2009-10
/////////////////////////////////////////////////////////////////////////////////

// V3 * M4
template <> inline
Vector3<float> operator * (const Vector3<float>& v, const Matrix4<float>& mtx)
{
    Vector3<float> res;
    D3DXVec3TransformCoord((D3DXVECTOR3*)&res, (CONST D3DXVECTOR3*)&v, (CONST D3DXMATRIX*)&mtx);
    return res;
}

// V4 * M4
template <> inline
Vector4<float> operator * (const Vector4<float>& v, const Matrix4<float>& mtx)
{
    Vector4<float> res;
    D3DXVec4Transform((D3DXVECTOR4*)&res, (CONST D3DXVECTOR4*)&v, (CONST D3DXMATRIX*)&mtx);
    return res;
}

template<> inline
Matrix4<float> Matrix4<float>::operator * (const Matrix4<float>& mtx) const
{
    Matrix4<float> result;
    D3DXMatrixMultiply((D3DXMATRIX*)&result, (CONST D3DXMATRIX*)this, (CONST D3DXMATRIX*)&mtx);
    return result;
}

template <> inline
Matrix4<float>& Matrix4<float>::operator *= (const Matrix4<float>& mtx)
{
    D3DXMatrixMultiply((D3DXMATRIX*)this, (CONST D3DXMATRIX*)this, (CONST D3DXMATRIX*)&mtx);
    return *this;
}

template <> inline
Matrix4<float> Matrix4<float>::MakeTranslation(const Vector3<float>& v)
{
	Matrix4<float> res;
    D3DXMatrixTranslation((D3DXMATRIX*)&res, v.x, v.y, v.z);
    return res;
}

template <> inline
Matrix4<float> Matrix4<float>::MakeTranslation(float x, float y, float z)
{
	Matrix4<float> res;
    D3DXMatrixTranslation((D3DXMATRIX*)&res, x, y, z);
    return res;
}

template <> inline
Matrix4<float> Matrix4<float>::MakeRotationX(float ang)
{
	Matrix4<float> res;
    D3DXMatrixRotationX((D3DXMATRIX*)&res, ang);
    return res;
}

template <> inline
Matrix4<float> Matrix4<float>::MakeRotationY(float ang)
{
	Matrix4<float> res;
    D3DXMatrixRotationY((D3DXMATRIX*)&res, ang);
    return res;
}

template <> inline
Matrix4<float> Matrix4<float>::MakeRotationZ(float ang)
{
	Matrix4<float> res;
    D3DXMatrixRotationZ((D3DXMATRIX*)&res, ang);
    return res;
}

template <> inline
Matrix4<float> Matrix4<float>::MakeRotationAxis(const Vector3<float>& axis, float ang)
{
	Matrix4<float> res;
    D3DXMatrixRotationAxis((D3DXMATRIX*)&res, (CONST D3DXVECTOR3*)&axis, ang);
    return res;
}

template <> inline
Matrix4<float> Matrix4<float>::MakeScaling(const Vector3<float>& v)
{
	Matrix4<float> res;
    D3DXMatrixScaling((D3DXMATRIX*)&res, v.x, v.y, v.z);
    return res;
}

template <> inline
Matrix4<float> Matrix4<float>::MakeScaling(float x, float y, float z)
{
	Matrix4<float> res;
    D3DXMatrixScaling((D3DXMATRIX*)&res, x, y, z);
    return res;
}

template <> inline
Matrix4<float>& Matrix4<float>::SetIdentity(void)
{
    D3DXMatrixIdentity((D3DXMATRIX*)this);
    return *this;
}

template <> inline
Matrix4<float> Matrix4<float>::Transpose(void) const
{
    Matrix4<float> result;
    D3DXMatrixTranspose((D3DXMATRIX*)&result, (CONST D3DXMATRIX*)this);
    return result;
}

template <> inline
Matrix4<float> Matrix4<float>::Inverse(void) const
{
    Matrix4<float> res;
    D3DXMatrixInverse((D3DXMATRIX*)&res, NULL, (CONST D3DXMATRIX*)this);
    return res;
}

template <> inline
float Matrix4<float>::Determinant(void) const
{
    return (float)D3DXMatrixDeterminant((D3DXMATRIX*)this);
}

template <> inline
Matrix4<float> Matrix4<float>::LookAtLH(const Vector3<float>& eye, const Vector3<float>& at, const Vector3<float>& up)
{
    Matrix4<float> res;
    D3DXMatrixLookAtLH( (D3DXMATRIX*)&res, (CONST D3DXVECTOR3*)&eye, (CONST D3DXVECTOR3*)&at, (CONST D3DXVECTOR3*)&up);
    return res;
}

template <> inline
Matrix4<float> Matrix4<float>::LookAtRH(const Vector3<float>& eye, const Vector3<float>& at, const Vector3<float>& up)
{
    Matrix4<float> res;
    D3DXMatrixLookAtRH( (D3DXMATRIX*)&res, (CONST D3DXVECTOR3*)&eye, (CONST D3DXVECTOR3*)&at, (CONST D3DXVECTOR3*)&up);
    return res;
}

template <> inline
Matrix4<float> Matrix4<float>::PerspectiveLH(float w, float h, float zn, float zf)
{
    Matrix4<float> res;
    D3DXMatrixPerspectiveLH((D3DXMATRIX*)&res, w, h, zn, zf);
    return res;
}

template <> inline
Matrix4<float> Matrix4<float>::PerspectiveRH(float w, float h, float zn, float zf)
{
    Matrix4<float> res;
    D3DXMatrixPerspectiveRH((D3DXMATRIX*)&res, w, h, zn, zf);
    return res;
}

template <> inline
Matrix4<float> Matrix4<float>::PerspectiveFovLH(float fovY, float aspect, float zn, float zf)
{
    Matrix4<float> res;
    D3DXMatrixPerspectiveFovLH((D3DXMATRIX*)&res, fovY, aspect, zn, zf);
    return res;
}

template <> inline
Matrix4<float> Matrix4<float>::PerspectiveFovRH(float fovY, float aspect, float zn, float zf)
{
    Matrix4<float> res;
    D3DXMatrixPerspectiveFovRH((D3DXMATRIX*)&res, fovY, aspect, zn, zf);
    return res;
}

template <> inline
Matrix4<float> Matrix4<float>::OrthoLH(float w, float h, float zn, float zf)
{
    Matrix4<float> res;
    D3DXMatrixOrthoLH((D3DXMATRIX*)&res, w, h, zn, zf);
    return res;
}

template <> inline
Matrix4<float> Matrix4<float>::OrthoRH(float w, float h, float zn, float zf)
{
    Matrix4<float> res;
    D3DXMatrixOrthoRH((D3DXMATRIX*)&res, w, h, zn, zf);
    return res;
}

//////////////////////////////////////////////////////////////////////////

template <> inline
Vector3<float> Vector3<float>::TransformCoord(const Matrix4<float>& mat)const
{
	Vector3<float> res;
	D3DXVec3TransformCoord((D3DXVECTOR3*)&res, (CONST D3DXVECTOR3*)this, (CONST D3DXMATRIX*)&mat);
	return res;
}

template <> inline
Vector3<float> Vector3<float>::TransformNormal(const Matrix4<float>& mat)const
{
	Vector3<float> res;
	D3DXVec3TransformNormal((D3DXVECTOR3*)&res, (CONST D3DXVECTOR3*)this, (CONST D3DXMATRIX*)&mat);
	return res;
}
