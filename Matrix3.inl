/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Matrix3.inl
/// \brief   3×3矩阵的实现
/// \author  DaiMingzhuang
/// \version 
/// \date    2009-02
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
const Matrix3<Real> Matrix3<Real>::ZERO(
    0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f);
template <typename Real>
const Matrix3<Real> Matrix3<Real>::IDENTITY(
    1.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 1.0f);

template <typename Real>
inline Matrix3<Real>::Matrix3(void)
{
}

template <typename Real>
inline Matrix3<Real>::Matrix3(const Matrix3<Real>& mtx)
{
	memcpy(m2, mtx.m2, sizeof(m2));
}

template <typename Real>
inline Matrix3<Real>::Matrix3(Real m00, Real m01, Real m02,
					   Real m10, Real m11, Real m12,
					   Real m20, Real m21, Real m22)
{
	Set(m00, m01, m02,
		m10, m11, m12,
		m20, m21, m22);
}

template <typename Real>
Matrix3<Real>& Matrix3<Real>::Set(Real m00, Real m01, Real m02,
								  Real m10, Real m11, Real m12,
								  Real m20, Real m21, Real m22)
{
	m[0][0] = m00;
	m[0][1] = m01;
	m[0][2] = m02;
	m[1][0] = m10;
	m[1][1] = m11;
	m[1][2] = m12;
	m[2][0] = m20;
	m[2][1] = m21;
	m[2][2] = m22;
	return *this;
}

template <typename Real>
inline Real Matrix3<Real>::operator() (int row, int col) const
{
	K_ASSERT(0<=row*3+col && row*3+col<9);
	return m[row][col];
}

template <typename Real>
inline Real& Matrix3<Real>::operator() (int row, int col)
{
	K_ASSERT(0<=row*3+col && row*3+col<9);
	return m[row][col];
}

template <typename Real>
inline const Vector3<Real>& Matrix3<Real>::operator[]( int row ) const
{
	K_ASSERT(0<=row && row<3);
	return *reinterpret_cast<const Vector3<Real>*>(&m[row]);
}

template <typename Real>
inline Vector3<Real>& Matrix3<Real>::operator[]( int row )
{
	K_ASSERT(0<=row && row<3);
	return *reinterpret_cast<Vector3<Real>*>(&m[row]);
}

template <typename Real>
inline Matrix3<Real>& Matrix3<Real>::operator = (const Matrix3& mtx)
{
	memcpy(m2, mtx.m2, sizeof(m2));
	return *this;
}

template <typename Real>
inline bool Matrix3<Real>::operator == (const Matrix3<Real>& mtx) const
{
	return (memcmp(&m2, &mtx.m2, sizeof(m2)) == 0);
}

template <typename Real>
inline bool Matrix3<Real>::operator != (const Matrix3<Real>& mtx) const
{
	return (memcmp(&m2, &mtx.m2, sizeof(m2)) != 0);
}

template <typename Real>
Matrix3<Real> Matrix3<Real>::operator + (const Matrix3<Real>& mtx) const
{
	Matrix3<Real> result(*this);
	for (int i=0; i<9; ++i)
		result.m2[i] += mtx.m2[i];
	return result;
}

template <typename Real>
Matrix3<Real> Matrix3<Real>::operator - (const Matrix3<Real>& mtx) const
{
	Matrix3<Real> result(*this);
	for (int i=0; i<9; ++i)
		result.m2[i] -= mtx.m2[i];
	return result;
}

template <typename Real>
Matrix3<Real> Matrix3<Real>::operator * (const Matrix3<Real>& mtx) const
{
	Matrix3<Real> result(ZERO);

	for (int i=0; i<3; ++i)
	{
		for (int j=0; j<3; ++j)
		{
			for (int k=0; k<3; ++k)
			{
				result.m[i][j] += m[i][k] * mtx.m[k][j];
			}
		}
	}
	return result;
}

template <typename Real>
Matrix3<Real> Matrix3<Real>::operator * (Real s) const
{
	Matrix3<Real> result(*this);
	for (int i=0; i<9; ++i)
		result.m2[i] *= s;
	return result;
}

template <typename Real>
inline Matrix3<Real> Matrix3<Real>::operator / (Real s) const
{
	s = 1.0f / s;
	return (*this) * s;
}

template <typename Real>
Matrix3<Real> Matrix3<Real>::operator - (void) const
{
	Matrix3<Real> result;
	for (int i=0; i<9; ++i)
		result.m2[i] = -m2[i];
	return result;
}

template <typename Real>
Matrix3<Real>& Matrix3<Real>::operator += (const Matrix3<Real>& mtx)
{
	for (int i=0; i<9; ++i)
		m2[i] += mtx.m2[i];
	return *this;
}

template <typename Real>
Matrix3<Real>& Matrix3<Real>::operator -= (const Matrix3<Real>& mtx)
{
	for (int i=0; i<9; ++i)
		m2[i] -= mtx.m2[i];
	return *this;
}

template <typename Real>
inline Matrix3<Real>& Matrix3<Real>::operator *= (const Matrix3<Real>& mtx)
{
	return *this = *this * mtx;
}

template <typename Real>
Matrix3<Real>& Matrix3<Real>::operator *= (Real s)
{
	for (int i=0; i<9; ++i)
		m2[i] *= s;
	return *this;
}

template <typename Real>
inline Matrix3<Real>& Matrix3<Real>::operator /= (Real s)
{
	s = 1.0f / s;
	return (*this) *= (s);
}

template <typename Real>
inline Matrix3<Real>& Matrix3<Real>::SetIdentity(void)
{
	memcpy(m2, IDENTITY.m2, sizeof(m2));
	return *this;
}

template <typename Real>
Matrix3<Real>& Matrix3<Real>::SetRotationX(Real ang)
{
	SetIdentity();

	Real fCos = cos(ang);
	Real fSin = sin(ang);

	m[1][1] = fCos;
	m[1][2] = fSin;
	m[2][1] = -fSin;
	m[2][2] = fCos;
	return *this;
}

template <typename Real>
Matrix3<Real>& Matrix3<Real>::SetRotationY(Real ang)
{
	SetIdentity();

	Real fCos = cos(ang);
	Real fSin = sin(ang);

	m[0][0] = fCos;
	m[0][2] = -fSin;
	m[2][0] = fSin;
	m[2][2] = fCos;
	return *this;
}

template <typename Real>
Matrix3<Real>& Matrix3<Real>::SetRotationZ(Real ang)
{
	SetIdentity();

	Real fCos = cos(ang);
	Real fSin = sin(ang);

	m[0][0] = fCos;
	m[0][1] = fSin;
	m[1][0] = -fSin;
	m[1][1] = fCos;
	return *this;
}

template <typename Real>
Matrix3<Real>& Matrix3<Real>::SetRotationAxis(const Vector3<Real>& axis, Real ang)
{
	Vector3<Real> axis_one = Normalize(axis);

	Real fCos = cos(ang);
	Real fSin = sin(ang);
	Real fOneMinusCos = 1.0f - fCos;

	Real fX2   = axis_one.x * axis_one.x;
	Real fY2   = axis_one.y * axis_one.y;
	Real fZ2   = axis_one.z * axis_one.z;
	Real fXYM  = axis_one.x * axis_one.y * fOneMinusCos;
	Real fXZM  = axis_one.x * axis_one.z * fOneMinusCos;
	Real fYZM  = axis_one.y * axis_one.z * fOneMinusCos;
	Real fXSin = axis_one.x * fSin;
	Real fYSin = axis_one.y * fSin;
	Real fZSin = axis_one.z * fSin;

	m[0][0] = fX2 * fOneMinusCos + fCos;
    m[0][1] = fXYM + fZSin;
    m[0][2] = fXZM - fYSin;

	m[1][0] = fXYM - fZSin;
    m[1][1] = fY2 * fOneMinusCos + fCos;
    m[1][2] = fYZM + fXSin;

	m[2][0] = fXZM + fYSin;
	m[2][1] = fYZM - fXSin;
	m[2][2] = fZ2 * fOneMinusCos + fCos;

	return *this;
}

template <typename Real>
Matrix3<Real>& Matrix3<Real>::SetRotationBy2Vector(const Vector3<Real>& dir0, const Vector3<Real>& dir1)
{
    if( MathUtil::Eq(dir0.Length(), (Real)0.0) ||  MathUtil::Eq(dir1.Length(), (Real)0.0) )
    { //如果有任一向量是零向量则返回单位矩阵。零向量的方向是任意的，总是可以认为两个向量已经是同向的
        this->SetIdentity();
        return *this;
    }

    Real dotP = Normalize(dir0).Dot(Normalize(dir1)); //得到两向量夹角的余弦值

    Real fAngle;

    if(dotP <= -1.0)
    { //反向
        fAngle = AngleUtil<Real>::PI;
    }
    else if(dotP >= 1.0)
    { //同向
        fAngle = 0;
    }
    else
    {
        fAngle = acos( dotP );
    }

    Vector3<Real> rkAxis = dir0.Cross(dir1); //计算转动轴

    //对转动轴长度为零的情况进行修正
    if( MathUtil::Eq(rkAxis.Length(), (Real)0.0) ) //dir0,dir1平行导致叉积长度为0 ，任意取一个垂直于dir0(dir1)的向量做转动轴即可
    {
        if(MathUtil::Eq(fAngle, (Real)0.0)) //如果同向
        {
            this->SetIdentity();
            return *this;
        }
        else //如果异向，注意不要生成零向量，同时尽量保证 SetRotationBy2Vector(d2,d1)应是SetRotationBy2Vector(d1,d2)的逆
        {
            if(dir0.x>0)
            {
                rkAxis = Vector3<Real>(dir0.y, -dir0.x, 0.0f);
            }
            else if(dir1.x>0)
            {
                rkAxis = Vector3<Real>(dir1.y, -dir1.x, 0.0f);
            }
            else if(dir0.y>0)
            {
                rkAxis = Vector3<Real>(dir0.y, -dir0.x, 0.0f);
            }
            else if(dir1.y>0)
            {
                rkAxis = Vector3<Real>(dir1.y, -dir1.x, 0.0f);
            }
            else if(dir0.z>0)
            {
                rkAxis = Vector3<Real>(0.0f, dir0.z, -dir0.y);
            }
            else if(dir1.z>0)
            {
                rkAxis = Vector3<Real>(0.0f, dir1.z, -dir1.y);
            }
            else
            {
                this->SetIdentity();
                return *this;
            }
        }
    }

    this->SetRotationAxis( Normalize(rkAxis) ,fAngle );
    return *this;
}

template <typename Real>
Matrix3<Real>& K::MATH::Matrix3<Real>::SetRotationQuaternion( const Quaternion<Real>& qtn )
{
	Real fTx  = qtn.x * 2.0f;
	Real fTy  = qtn.y * 2.0f;
	Real fTz  = qtn.z * 2.0f;
	Real fTwx = fTx * qtn.w;
	Real fTwy = fTy * qtn.w;
	Real fTwz = fTz * qtn.w;
	Real fTxx = fTx * qtn.x;
	Real fTxy = fTy * qtn.x;
	Real fTxz = fTz * qtn.x;
	Real fTyy = fTy * qtn.y;
	Real fTyz = fTz * qtn.y;
	Real fTzz = fTz * qtn.z;

	m[0][0] = 1.0f - (fTyy+fTzz);
	m[0][1] = fTxy + fTwz;
	m[0][2] = fTxz - fTwy;

	m[1][0] = fTxy - fTwz;
	m[1][1] = 1.0f - (fTxx+fTzz);
	m[1][2] = fTyz + fTwx;

	m[2][0] = fTxz + fTwy;
	m[2][1] = fTyz - fTwx;
	m[2][2] = 1.0f - (fTxx+fTyy);

	return *this;
}

template <typename Real>
Quaternion<Real> K::MATH::Matrix3<Real>::GetRotationQuaternion() const
{
	Quaternion<Real> q;
	Real m[3][3] = {
		this->m[0][0], this->m[0][1], this->m[0][2],
		this->m[1][0], this->m[1][1], this->m[1][2],
		this->m[2][0], this->m[2][1], this->m[2][2],
	};
	Real fTrace = m[0][0] + m[1][1] + m[2][2];

	if (fTrace > 0.0f)
	{
		Real fRoot = MathUtil::Sqrt(fTrace + 1.0f);
		q.w = fRoot * 0.5f;

		fRoot = 0.5f / fRoot;
		q.x = (m[1][2] - m[2][1]) * fRoot;
		q.y = (m[2][0] - m[0][2]) * fRoot;
		q.z = (m[0][1] - m[1][0]) * fRoot;
	}
	else
	{
		if (m[0][0] >= m[1][1] && m[0][0] >= m[2][2])
		{
			Real fRoot = MathUtil::Sqrt(m[0][0] - m[1][1] - m[2][2] + 1.0f);
			q.x = fRoot * 0.5f;

			fRoot = 0.5f / fRoot;
			q.w = (m[1][2] - m[2][1]) * fRoot;
			q.y = (m[0][1] + m[1][0]) * fRoot;
			q.z = (m[0][2] + m[2][0]) * fRoot;
		}
		else if (m[1][1] >= m[0][0] && m[1][1] >= m[2][2])
		{
			Real fRoot = MathUtil::Sqrt(m[1][1] - m[2][2] - m[0][0] + 1.0f);
			q.y = fRoot * 0.5f;

			fRoot = 0.5f / fRoot;
			q.w = (m[2][0] - m[0][2]) * fRoot;
			q.z = (m[1][2] + m[2][1]) * fRoot;
			q.x = (m[1][0] + m[0][1]) * fRoot;
		}
		else
		{
			Real fRoot = MathUtil::Sqrt(m[2][2] - m[0][0] - m[1][1] + 1.0f);
			q.z = fRoot * 0.5f;

			fRoot = 0.5f / fRoot;
			q.w = (m[0][1] - m[1][0]) * fRoot;
			q.x = (m[2][0] + m[0][2]) * fRoot;
			q.y = (m[2][1] + m[1][2]) * fRoot;
		}
	}
	return q;
}

template <typename Real>
Matrix3<Real>& Matrix3<Real>::SetEulerXYZ(const Vector3<Real>& angles)
{
	return (*this = MakeRotationX(angles.x) * MakeRotationY(angles.y) * MakeRotationZ(angles.z));
}

template <typename Real>
Matrix3<Real>& Matrix3<Real>::SetEulerXZY(const Vector3<Real>& angles)
{
	return (*this = MakeRotationX(angles.x) * MakeRotationZ(angles.z) * MakeRotationY(angles.y));
}

template <typename Real>
Matrix3<Real>& Matrix3<Real>::SetEulerYXZ(const Vector3<Real>& angles)
{
	return (*this = MakeRotationY(angles.y) * MakeRotationX(angles.x) * MakeRotationZ(angles.z));
}

template <typename Real>
Matrix3<Real>& Matrix3<Real>::SetEulerYZX(const Vector3<Real>& angles)
{
	return (*this = MakeRotationY(angles.y) * MakeRotationZ(angles.z) * MakeRotationX(angles.x));
}

template <typename Real>
Matrix3<Real>& Matrix3<Real>::SetEulerZXY(const Vector3<Real>& angles)
{
	return (*this = MakeRotationZ(angles.z) * MakeRotationX(angles.x) * MakeRotationY(angles.y));
}

template <typename Real>
Matrix3<Real>& Matrix3<Real>::SetEulerZYX(const Vector3<Real>& angles)
{
	return (*this = MakeRotationZ(angles.z) * MakeRotationY(angles.y) * MakeRotationX(angles.x));
}

template <typename Real>
Vector3<Real> Matrix3<Real>::GetEulerXYZ() const
{
	// rot =  cy*cz           cy*sz           -sy
	//        sx*sy*cz-cx*sz  sx*sy*sz+cx*cz  sx*cy
	//        cx*sy*cz+sx*sz  cx*sy*sz-sx*cz  cx*cy
	if (m[0][2] >= (Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>(-MathUtil::Atan2(m[1][0], m[1][1]), -AngleUtil<Real>::HALF_PI, (Real)0.0);
	}
	else if (m[0][2] <= -(Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>(MathUtil::Atan2(m[1][0], m[1][1]), AngleUtil<Real>::HALF_PI, (Real)0.0);
	}
	else
	{
		return Vector3<Real>(MathUtil::Atan2(m[1][2], m[2][2]), MathUtil::Asin(-m[0][2]), MathUtil::Atan2(m[0][1], m[0][0]));
	}
}

template <typename Real>
Vector3<Real> Matrix3<Real>::GetEulerXZY() const
{
	// rot =  cy*cz           sz             -sy*cz
	//       -cx*cy*sz+sx*sy  cx*cz           cx*sy*sz+sx*cy
	//        sx*cy*sz+cx*sy -sx*cz          -sx*sy*sz+cx*cy
	if (m2[1] >= (Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>(MathUtil::Atan2(m[2][0], m[2][2]), (Real)0.0, AngleUtil<Real>::HALF_PI);
	}
	else if (m2[1] <= -(Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>(MathUtil::Atan2(-m[2][0], m[2][2]), (Real)0.0, -AngleUtil<Real>::HALF_PI);
	}
	else
	{
		return Vector3<Real>(MathUtil::Atan2(-m[2][1], m[1][1]), MathUtil::Atan2(-m[0][2], m[0][0]), MathUtil::Asin(m[0][1]));
	}
}

template <typename Real>
Vector3<Real> Matrix3<Real>::GetEulerYXZ() const
{
	if (m2[5] >= (Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>(AngleUtil<Real>::HALF_PI, MathUtil::Atan2(m2[1], m2[0]), (Real)0.0);
	}
	else if (m2[5] <= -(Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>(-AngleUtil<Real>::HALF_PI, MathUtil::Atan2(-m2[1], m2[0]), (Real)0.0);
	}
	else
	{
		return Vector3<Real>(MathUtil::Asin(m2[5]), MathUtil::Atan2(-m2[2], m2[8]), MathUtil::Atan2(-m2[3], m2[4]));
	}
}

template <typename Real>
Vector3<Real> Matrix3<Real>::GetEulerYZX() const
{
	if (m2[3] >= (Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>((Real)0.0, -MathUtil::Atan2(m2[7], m2[8]), -AngleUtil<Real>::HALF_PI);
	}
	else if (m2[3] <= -(Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>((Real)0.0, MathUtil::Atan2(m2[7], m2[8]), AngleUtil<Real>::HALF_PI);
	}
	else
	{
		return Vector3<Real>(MathUtil::Atan2(m2[5], m2[4]), MathUtil::Atan2(m2[6], m2[0]), MathUtil::Asin(-m2[3]));
	}
}

template <typename Real>
Vector3<Real> Matrix3<Real>::GetEulerZXY() const
{
	if (m2[7] >= (Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>(-AngleUtil<Real>::HALF_PI, (Real)0.0, -MathUtil::Atan2(m2[2], m2[0]));
	}
	else if (m2[7] <= -(Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>(AngleUtil<Real>::HALF_PI, (Real)0.0, MathUtil::Atan2(m2[2], m2[0]));
	}
	else
	{
		return Vector3<Real>(MathUtil::Asin(-m2[7]), MathUtil::Atan2(m2[6], m2[8]), MathUtil::Atan2(m2[1], m2[4]));
	}
}

template <typename Real>
Vector3<Real> Matrix3<Real>::GetEulerZYX() const
{
	// rot =  cy*cz           cx*sz+sx*sy*cz  sx*sz-cx*cz*sy
	//       -cy*sz           cx*cz-sx*sy*sz  sx*cz+cx*sy*sz
	//        sy             -sx*cy           cx*cy
	if (m2[6] >= (Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>((Real)0.0, AngleUtil<Real>::HALF_PI, MathUtil::MathUtil::Atan2(m[0][1], -m[0][2]));
	}
	else if (m2[6] <= -(Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>((Real)0.0, -AngleUtil<Real>::HALF_PI, MathUtil::Atan2(m[0][1], m[0][2]));
	}
	else
	{
		return Vector3<Real>(MathUtil::MathUtil::Atan2(-m[2][1], m[2][2]), MathUtil::MathUtil::Asin(m[2][0]), MathUtil::MathUtil::Atan2(-m[1][0], m[0][0]));
	}
}

template <typename Real>
Matrix3<Real> Matrix3<Real>::Transpose(void) const
{
	Matrix3<Real> result;

	for (int i=0; i<3; ++i)
	{
		for (int j=0; j<3; ++j)
		{
			result.m[i][j] = m[j][i];
		}
	}

	return result;
}

template <typename Real>
Matrix3<Real> Matrix3<Real>::Inverse(void) const
{
	Real fCo00 = m2[4]*m2[8] - m2[5]*m2[7];
	Real fCo10 = m2[5]*m2[6] - m2[3]*m2[8];
	Real fCo20 = m2[3]*m2[7] - m2[4]*m2[6];

	Real fDet = m2[0]*fCo00 + m2[1]*fCo10 + m2[2]*fCo20;
	if (MathUtil::Abs(fDet) <= (Real)0.0)
	{
		K_MATH_LOG("Singular matrix has no inversion.");
		return Matrix3<Real>::ZERO;
	}

	Matrix3<Real> result;
	result.m[0][0] = fCo00;
	result.m[0][1] = m2[2]*m2[7] - m2[1]*m2[8];
	result.m[0][2] = m2[1]*m2[5] - m2[2]*m2[4];
	result.m[1][0] = fCo10;
	result.m[1][1] = m2[0]*m2[8] - m2[2]*m2[6];
	result.m[1][2] = m2[2]*m2[3] - m2[0]*m2[5];
	result.m[2][0] = fCo20;
	result.m[2][1] = m2[1]*m2[6] - m2[0]*m2[7];
	result.m[2][2] = m2[0]*m2[4] - m2[1]*m2[3];
	return result *= 1.0f / fDet;
}

template <typename Real>
Real Matrix3<Real>::Determinant(void) const
{
	Real fCo00 = m2[4]*m2[8] - m2[5]*m2[7];
	Real fCo10 = m2[5]*m2[6] - m2[3]*m2[8];
	Real fCo20 = m2[3]*m2[7] - m2[4]*m2[6];
	Real fDet = m2[0]*fCo00 + m2[1]*fCo10 + m2[2]*fCo20;
	return fDet;
}

template <typename Real>
Matrix4<Real> Matrix3<Real>::ToMatrix4() const
{
    Matrix4<Real> res;
    res.SetIdentity();
    memcpy(res.m[0], m[0], sizeof(m[0]));
    memcpy(res.m[1], m[1], sizeof(m[1]));
    memcpy(res.m[2], m[2], sizeof(m[2]));
    return res;
}

template <typename Real>
Vector3<Real> operator * (const Vector3<Real>& v, const Matrix3<Real>& mtx)
{
	return Vector3<Real>(
		v.x * mtx.m[0][0] + v.y * mtx.m[1][0] + v.z * mtx.m[2][0],
		v.x * mtx.m[0][1] + v.y * mtx.m[1][1] + v.z * mtx.m[2][1],
		v.x * mtx.m[0][2] + v.y * mtx.m[1][2] + v.z * mtx.m[2][2]
		);
}

template <typename Real>
Matrix3<Real> operator* (Real fScalar, const Matrix3<Real>& mtx)
{
	return mtx * fScalar;
}
