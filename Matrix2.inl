/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Matrix2.inl
/// \brief   2°¡2æÿ’Ûµƒ µœ÷
/// \author  DaiMingzhuang
/// \version 
/// \date    2009-02
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
const Matrix2<Real> Matrix2<Real>::ZERO(
    0.0f, 0.0f,
    0.0f, 0.0f);
template <typename Real>
const Matrix2<Real> Matrix2<Real>::IDENTITY(
    1.0f, 0.0f,
    0.0f, 1.0f);

template <typename Real>
inline Matrix2<Real>::Matrix2(void)
{
}

template <typename Real>
inline Matrix2<Real>::Matrix2(const Matrix2<Real>& mtx)
{
	memcpy(m2, mtx.m2, sizeof(m2));
}

template <typename Real>
inline Matrix2<Real>::Matrix2(Real m00, Real m01, Real m10, Real m11)
{
	Set(m00, m01, m10, m11);
}

template <typename Real>
Matrix2<Real>& Matrix2<Real>::Set(Real m00, Real m01, Real m10, Real m11)
{
	m[0][0] = m00;
	m[0][1] = m01;
	m[1][0] = m10;
	m[1][1] = m11;
	return *this;
}

template <typename Real>
inline Real Matrix2<Real>::operator() (int row, int col) const
{
	K_ASSERT(0<=row*2+col && row*2+col<4);
	return m[row][col];
}

template <typename Real>
inline Real& Matrix2<Real>::operator() (int row, int col)
{
	K_ASSERT(0<=row*2+col && row*2+col<4);
	return m[row][col];
}

template <typename Real>
inline const Vector2<Real>& Matrix2<Real>::operator[]( int row ) const
{
	K_ASSERT(0<=row && row<2);
	return *reinterpret_cast<const Vector2<Real>*>(&m[row]);
}

template <typename Real>
inline Vector2<Real>& Matrix2<Real>::operator[]( int row )
{
	K_ASSERT(0<=row && row<2);
	return *reinterpret_cast<Vector2<Real>*>(&m[row]);
}

template <typename Real>
inline Matrix2<Real>& Matrix2<Real>::operator = (const Matrix2& mtx)
{
	memcpy(m2, mtx.m2, sizeof(m2));
	return *this;
}

template <typename Real>
inline bool Matrix2<Real>::operator == (const Matrix2<Real>& mtx) const
{
	return (memcmp(&m2, &mtx.m2, sizeof(m2)) == 0);
}

template <typename Real>
inline bool Matrix2<Real>::operator != (const Matrix2<Real>& mtx) const
{
	return (memcmp(&m2, &mtx.m2, sizeof(m2)) != 0);
}

template <typename Real>
Matrix2<Real> Matrix2<Real>::operator + (const Matrix2<Real>& mtx) const
{
	Matrix2<Real> result(*this);
	for (int i=0; i<4; ++i)
		result.m2[i] += mtx.m2[i];
	return result;
}

template <typename Real>
Matrix2<Real> Matrix2<Real>::operator - (const Matrix2<Real>& mtx) const
{
	Matrix2<Real> result(*this);
	for (int i=0; i<4; ++i)
		result.m2[i] -= mtx.m2[i];
	return result;
}

template <typename Real>
Matrix2<Real> Matrix2<Real>::operator * (const Matrix2<Real>& mtx) const
{
	Matrix2<Real> result;
	result.m[0][0] = m[0][0] * mtx.m[0][0] + m[0][1] * mtx.m[1][0];
	result.m[0][1] = m[0][0] * mtx.m[0][1] + m[0][1] * mtx.m[1][1];
	result.m[1][0] = m[1][0] * mtx.m[0][0] + m[1][1] * mtx.m[1][0];
	result.m[1][1] = m[1][0] * mtx.m[0][1] + m[1][1] * mtx.m[1][1];
	return result;
}

template <typename Real>
Matrix2<Real> Matrix2<Real>::operator * (Real s) const
{
	Matrix2<Real> result(*this);
	for (int i=0; i<4; ++i)
		result.m2[i] *= s;
	return result;
}

template <typename Real>
Matrix2<Real> Matrix2<Real>::operator / (Real s) const
{
	s = 1.0f / s;
	return (*this) * s;
}

template <typename Real>
Matrix2<Real> Matrix2<Real>::operator - (void) const
{
	Matrix2<Real> result;
	for (int i=0; i<4; ++i)
		result.m2[i] = -m2[i];
	return result;
}

template <typename Real>
Matrix2<Real>& Matrix2<Real>::operator += (const Matrix2<Real>& mtx)
{
	for (int i=0; i<4; ++i)
		m2[i] += mtx.m2[i];
	return *this;
}

template <typename Real>
Matrix2<Real>& Matrix2<Real>::operator -= (const Matrix2<Real>& mtx)
{
	for (int i=0; i<4; ++i)
		m2[i] -= mtx.m2[i];
	return *this;
}

template <typename Real>
inline Matrix2<Real>& Matrix2<Real>::operator *= (const Matrix2<Real>& mtx)
{
	return *this = *this * mtx;
}

template <typename Real>
Matrix2<Real>& Matrix2<Real>::operator *= (Real s)
{
	for (int i=0; i<4; ++i)
		m2[i] *= s;
	return *this;
}

template <typename Real>
inline Matrix2<Real>& Matrix2<Real>::operator /= (Real s)
{
	s = 1.0f / s;
	return (*this) *= (s);
}

template <typename Real>
inline Matrix2<Real>& Matrix2<Real>::SetIdentity(void)
{
	return (*this = IDENTITY);
}

template <typename Real>
Matrix2<Real>& Matrix2<Real>::Set2DRotation(Real ang)
{
	Real fCos = cos(ang);
	Real fSin = sin(ang);

	m[0][0] = fCos;
	m[0][1] = fSin;
	m[1][0] = -fSin;
	m[1][1] = fCos;
	return *this;
}

template <typename Real>
Matrix2<Real> Matrix2<Real>::Transpose(void) const
{
	Matrix2<Real> result;
	result.m[0][0] = m[0][0];
	result.m[0][1] = m[1][0];
	result.m[1][0] = m[0][1];
	result.m[1][1] = m[1][1];
	return result;
}

template <typename Real>
Matrix2<Real> Matrix2<Real>::Inverse(void) const
{
	Real fDet = m2[0] *m2[3] - m2[1] * m2[2];
	if (MathUtil::Abs(fDet) <= (Real)0.0)
    {
        K_MATH_LOG("Singular matrix has no inversion.");
		return Matrix2<Real>::ZERO;
    }

	Matrix2<Real> result;
	result.m2[0] =  m2[3];
	result.m2[1] = -m2[1];
	result.m2[2] = -m2[2];
	result.m2[3] =  m2[0];
	return result *= 1.0f / fDet;
}

template <typename Real>
inline Real Matrix2<Real>::Determinant(void) const
{
	return (m2[0] * m2[3] - m2[1] * m2[2]);
}

template <typename Real>
Vector2<Real> operator * (const Vector2<Real>& v, const Matrix2<Real>& mtx)
{
	return Vector2<Real>(v.x * mtx.m[0][0] + v.y * mtx.m[1][0], v.x * mtx.m[0][1] + v.y * mtx.m[1][1]);
}

template <typename Real>
Matrix2<Real> operator* (Real fScalar, const Matrix2<Real>& mtx)
{
	return mtx * fScalar;
}
