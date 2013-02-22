/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Matrix4.inl
/// \brief   4×4矩阵的实现
/// \author  DaiMingzhuang
/// \version 
/// \date    2009-02
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
const Matrix4<Real> Matrix4<Real>::ZERO(
    0.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 0.0f);
template <typename Real>
const Matrix4<Real> Matrix4<Real>::IDENTITY(
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f);

template <typename Real>
inline Matrix4<Real>::Matrix4(void)
{
}

template <typename Real>
inline Matrix4<Real>::Matrix4(Real m00, Real m01, Real m02, Real m03,
					   Real m10, Real m11, Real m12, Real m13,
					   Real m20, Real m21, Real m22, Real m23,
					   Real m30, Real m31, Real m32, Real m33)
{
	Set(m00, m01, m02, m03,
		m10, m11, m12, m13,
		m20, m21, m22, m23,
		m30, m31, m32, m33);
}

template <typename Real>
inline Matrix4<Real>::Matrix4(const Matrix4<Real>& mtx)
{
	memcpy(m2, mtx.m2, sizeof(m2));
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::Set(Real m00, Real m01, Real m02, Real m03,
								  Real m10, Real m11, Real m12, Real m13,
								  Real m20, Real m21, Real m22, Real m23,
								  Real m30, Real m31, Real m32, Real m33)
{
	m[0][0] = m00;
	m[0][1] = m01;
	m[0][2] = m02;
	m[0][3] = m03;
	m[1][0] = m10;
	m[1][1] = m11;
	m[1][2] = m12;
	m[1][3] = m13;
	m[2][0] = m20;
	m[2][1] = m21;
	m[2][2] = m22;
	m[2][3] = m23;
	m[3][0] = m30;
	m[3][1] = m31;
	m[3][2] = m32;
	m[3][3] = m33;
	return *this;
}

template <typename Real>
inline Matrix4<Real>& Matrix4<Real>::operator = (const Matrix4& mtx)
{
	memcpy(m2, mtx.m2, sizeof(m2));
	return *this;
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::FromMatrix3(const Matrix3<Real>& mtx)
{
	SetIdentity();
	memcpy(m[0], mtx.m[0], sizeof(mtx.m[0]));
	memcpy(m[1], mtx.m[1], sizeof(mtx.m[1]));
	memcpy(m[2], mtx.m[2], sizeof(mtx.m[2]));
	return *this;
}

template <typename Real>
Matrix3<Real> Matrix4<Real>::ToMatrix3(void) const
{
	Matrix3<Real> mtx;
	memcpy(mtx.m[0], m[0], sizeof(mtx.m[0]));
	memcpy(mtx.m[1], m[1], sizeof(mtx.m[1]));
	memcpy(mtx.m[2], m[2], sizeof(mtx.m[2]));
	return mtx;
}

template <typename Real>
inline Real Matrix4<Real>::operator() (int row, int col) const
{
	K_ASSERT(0<=row*4+col && row*4+col<16);
	return m[row][col];
}

template <typename Real>
inline Real& Matrix4<Real>::operator() (int row, int col)
{
	K_ASSERT(0<=row*4+col && row*4+col<16);
	return m[row][col];
}

template <typename Real>
inline const Vector4<Real>& Matrix4<Real>::operator[]( int row ) const
{
	K_ASSERT(0<=row && row<4);
	return *reinterpret_cast<const Vector4<Real>*>(&m[row]);
}

template <typename Real>
inline Vector4<Real>& Matrix4<Real>::operator[]( int row )
{
	K_ASSERT(0<=row && row<4);
	return *reinterpret_cast<Vector4<Real>*>(&m[row]);
}

template <typename Real>
inline bool Matrix4<Real>::operator == (const Matrix4<Real>& mtx) const
{
	return (memcmp(&m2, &mtx.m2, sizeof(m2)) == 0);
}

template <typename Real>
inline bool Matrix4<Real>::operator != (const Matrix4<Real>& mtx) const
{
	return (memcmp(&m2, &mtx.m2, sizeof(m2)) != 0);
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::operator + (const Matrix4<Real>& mtx) const
{
	Matrix4<Real> result(*this);
	for (int i=0; i<16; ++i)
		result.m2[i] += mtx.m2[i];
	return result;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::operator - (const Matrix4<Real>& mtx) const
{
	Matrix4<Real> result(*this);
	for (int i=0; i<16; ++i)
		result.m2[i] -= mtx.m2[i];
	return result;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::operator * (const Matrix4<Real>& mtx) const
{
	Matrix4<Real> result(ZERO);

	for (int i=0; i<4; ++i)
	{
		for (int j=0; j<4; ++j)
		{
			for (int k=0; k<4; ++k)
			{
				result.m[i][j] += m[i][k] * mtx.m[k][j];
			}
		}
	}

	return result;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::operator * (Real s) const
{
	Matrix4<Real> result(*this);
	for (int i=0; i<16; ++i)
		result.m2[i] *= s;
	return result;
}

template <typename Real>
inline Matrix4<Real> Matrix4<Real>::operator / (Real s) const
{
	s = 1.0f / s;
	return (*this) * s;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::operator - (void) const
{
	Matrix4<Real> result;
	for (int i=0; i<16; ++i)
		result.m2[i] = -m2[i];
	return result;
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::operator += (const Matrix4<Real>& mtx)
{
	for (int i=0; i<16; ++i)
		m2[i] += mtx.m2[i];
	return *this;
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::operator -= (const Matrix4<Real>& mtx)
{
	for (int i=0; i<16; ++i)
		m2[i] -= mtx.m2[i];
	return *this;
}

template <typename Real>
inline Matrix4<Real>& Matrix4<Real>::operator *= (const Matrix4<Real>& mtx)
{
	return *this = *this * mtx;
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::operator *= (Real s)
{
	for (int i=0; i<16; ++i)
		m2[i] *= s;
	return *this;
}

template <typename Real>
inline Matrix4<Real>& Matrix4<Real>::operator /= (Real s)
{
	s = 1.0f / s;
	return (*this) *= (s);
}

template <typename Real>
inline Matrix4<Real>& Matrix4<Real>::SetIdentity(void)
{
	memcpy(m2, IDENTITY.m2, sizeof(m2));
	return *this;
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::SetTranslation(const Vector3<Real>& v)
{
	m[3][0] = v.x;
	m[3][1] = v.y;
	m[3][2] = v.z;
	return *this;
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::SetTranslation(Real x, Real y, Real z)
{
	m[3][0] = x;
	m[3][1] = y;
	m[3][2] = z;
	return *this;
}

template <typename Real>
Vector3<Real> Matrix4<Real>::GetTranslation() const
{
	return Vector3<Real>(m[3][0], m[3][1], m[3][2]);
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::SetRotationX(Real ang)
{
	Vector3<Real> scaling = GetScaling();

	Real fCos = MathUtil::Cos(ang);
	Real fSin = MathUtil::Sin(ang);

	m[0][0] = scaling[0];
	m[0][1] = 0;
	m[0][2] = 0;

	m[1][0] = 0;
	m[1][1] = scaling[1] * fCos;
	m[1][2] = scaling[1] * fSin;

	m[2][0] = 0
	m[2][1] = scaling[2] * -fSin;
	m[2][2] = scaling[2] * fCos;

	return *this;
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::SetRotationY(Real ang)
{
	Vector3<Real> scaling = GetScaling();

	Real fCos = MathUtil::Cos(ang);
	Real fSin = MathUtil::Sin(ang);

	m[0][0] = scaling[0] * fCos;
	m[0][1] = 0;
	m[0][2] = scaling[0] * -fSin;

	m[1][0] = 0;
	m[1][1] = scaling[1];
	m[1][2] = 0;

	m[2][0] = scaling[2] * fSin;
	m[2][1] = 0;
	m[2][2] = scaling[2] * fCos;

	return *this;
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::SetRotationZ(Real ang)
{
	Vector3<Real> scaling = GetScaling();

	Real fCos = MathUtil::Cos(ang);
	Real fSin = MathUtil::Sin(ang);

	m[0][0] = scaling[0] * fCos;
	m[0][1] = scaling[0] * fSin;
	m[0][2] = 0;

	m[1][0] = scaling[1] * -fSin;
	m[1][1] = scaling[1] * fCos;
	m[1][2] = 0;

	m[2][0] = 0;
	m[2][1] = 0;
	m[2][2] = scaling[2];

	return (*this);
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::SetRotationAxis(const Vector3<Real>& axis, Real ang)
{
	Vector3<Real> scaling = GetScaling();
	Vector3<Real> axis_one = Normalize(axis);

	Real fCos = cos(ang);
	Real fSin = sin(ang);
	Real fOneMinusCos = (Real)1.0 - fCos;

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
	return Scale(scaling);
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::SetRotationBy2Vector(const Vector3<Real>& dir0, const Vector3<Real>& dir1)
{
    if( MathUtil::Eq(dir0.Length(), (Real)0.0) ||  MathUtil::Eq(dir1.Length(), (Real)0.0) )
    { //如果有任一向量是零向量则返回单位矩阵。零向量的方向是任意的，总是可以认为两个向量已经是同向的
        this->ClearRotation();
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
            this->ClearRotation();
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
                this->ClearRotation();
                return *this;
            }
        }
    }

    this->SetRotationAxis( Normalize(rkAxis) ,fAngle );
    return *this;
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::SetRotationQuaternion( const Quaternion<Real>& qtn )
{
	Vector3<Real> scaling = GetScaling();
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

	return Scale(scaling);
}

template <typename Real>
Quaternion<Real> Matrix4<Real>::GetRotationQuaternion() const
{
	Quaternion<Real> q;
	Vector3<Real> scaling = GetScaling();
	Real m[3][3] = {
		this->m[0][0] / scaling[0], this->m[0][1] / scaling[0], this->m[0][2] / scaling[0],
		this->m[1][0] / scaling[1], this->m[1][1] / scaling[1], this->m[1][2] / scaling[1],
		this->m[2][0] / scaling[2], this->m[2][1] / scaling[2], this->m[2][2] / scaling[2],
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
Matrix4<Real>& Matrix4<Real>::SetEulerXYZ(const Vector3<Real>& angles)
{
	Vector3<Real> scaling = GetScaling();
	Matrix3<Real> rot = Matrix3<Real>::MakeRotationX(angles.x) * Matrix3<Real>::MakeRotationY(angles.y) * Matrix3<Real>::MakeRotationZ(angles.z);
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			m[i][j] = scaling[i] * rot.m[i][j];
		}
	}
	return (*this);
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::SetEulerXZY(const Vector3<Real>& angles)
{
	Vector3<Real> scaling = GetScaling();
	Matrix3<Real> rot = Matrix3<Real>::MakeRotationX(angles.x) * Matrix3<Real>::MakeRotationZ(angles.z) * Matrix3<Real>::MakeRotationY(angles.y);
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			m[i][j] = scaling[i] * rot.m[i][j];
		}
	}
	return (*this);
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::SetEulerYXZ(const Vector3<Real>& angles)
{
	Vector3<Real> scaling = GetScaling();
	Matrix3<Real> rot = Matrix3<Real>::MakeRotationY(angles.y) * Matrix3<Real>::MakeRotationX(angles.x) * Matrix3<Real>::MakeRotationZ(angles.z);
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			m[i][j] = scaling[i] * rot.m[i][j];
		}
	}
	return (*this);
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::SetEulerYZX(const Vector3<Real>& angles)
{
	Vector3<Real> scaling = GetScaling();
	Matrix3<Real> rot = Matrix3<Real>::MakeRotationY(angles.y) * Matrix3<Real>::MakeRotationZ(angles.z) * Matrix3<Real>::MakeRotationX(angles.x);
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			m[i][j] = scaling[i] * rot.m[i][j];
		}
	}
	return (*this);
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::SetEulerZXY(const Vector3<Real>& angles)
{
	Vector3<Real> scaling = GetScaling();
	Matrix3<Real> rot = Matrix3<Real>::MakeRotationZ(angles.z) * Matrix3<Real>::MakeRotationX(angles.x) * Matrix3<Real>::MakeRotationY(angles.y);
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			m[i][j] = scaling[i] * rot.m[i][j];
		}
	}
	return (*this);
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::SetEulerZYX(const Vector3<Real>& angles)
{
	Vector3<Real> scaling = GetScaling();
	Matrix3<Real> rot = Matrix3<Real>::MakeRotationZ(angles.z) * Matrix3<Real>::MakeRotationY(angles.y) * Matrix3<Real>::MakeRotationX(angles.x);
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			m[i][j] = scaling[i] * rot.m[i][j];
		}
	}
	return (*this);
}

template <typename Real>
Vector3<Real> Matrix4<Real>::GetEulerXYZ() const
{
	Vector3<Real> scaling = GetScaling();
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
		return Vector3<Real>(MathUtil::Atan2(m[1][2]*scaling[2], m[2][2]*scaling[1]), MathUtil::Asin(-m[0][2]/scaling[0]), MathUtil::Atan2(m[0][1], m[0][0]));
	}
}

template <typename Real>
Vector3<Real> Matrix4<Real>::GetEulerXZY() const
{
	Vector3<Real> scaling = GetScaling();
	// rot =  cy*cz           sz             -sy*cz
	//       -cx*cy*sz+sx*sy  cx*cz           cx*sy*sz+sx*cy
	//        sx*cy*sz+cx*sy -sx*cz          -sx*sy*sz+cx*cy
	if (m[0][1] >= (Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>(MathUtil::Atan2(m[2][0], m[2][2]), (Real)0.0, AngleUtil<Real>::HALF_PI);
	}
	else if (m[0][1] <= -(Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>(MathUtil::Atan2(-m[2][0], m[2][2]), (Real)0.0, -AngleUtil<Real>::HALF_PI);
	}
	else
	{
		return Vector3<Real>(MathUtil::Atan2(-m[2][1]*scaling[1], m[1][1]*scaling[2]), MathUtil::Atan2(-m[0][2], m[0][0]), MathUtil::Asin(m[0][1]/scaling[0]));
	}
}

template <typename Real>
Vector3<Real> Matrix4<Real>::GetEulerYXZ() const
{
	Vector3<Real> scaling = GetScaling();
	if (m[1][2] >= (Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>(AngleUtil<Real>::HALF_PI, MathUtil::Atan2(m[0][1], m[0][0]), (Real)0.0);
	}
	else if (m[1][2] <= -(Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>(-AngleUtil<Real>::HALF_PI, MathUtil::Atan2(-m[0][1], m[0][0]), (Real)0.0);
	}
	else
	{
		return Vector3<Real>(MathUtil::Asin(m[1][2]/scaling[1]), MathUtil::Atan2(-m[0][2]*scaling[2], m[2][2]*scaling[0]), MathUtil::Atan2(-m[1][0], m[1][1]));
	}
}

template <typename Real>
Vector3<Real> Matrix4<Real>::GetEulerYZX() const
{
	Vector3<Real> scaling = GetScaling();
	if (m[1][0] >= (Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>((Real)0.0, -MathUtil::Atan2(m[2][1], m[2][2]), -AngleUtil<Real>::HALF_PI);
	}
	else if (m[1][0] <= -(Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>((Real)0.0, MathUtil::Atan2(m[2][1], m[2][2]), AngleUtil<Real>::HALF_PI);
	}
	else
	{
		return Vector3<Real>(MathUtil::Atan2(m[1][2], m[1][1]), MathUtil::Atan2(m[2][0]*scaling[0], m[0][0]*scaling[2]), MathUtil::Asin(-m[1][0]/scaling[1]));
	}
}

template <typename Real>
Vector3<Real> Matrix4<Real>::GetEulerZXY() const
{
	Vector3<Real> scaling = GetScaling();
	if (m[2][1] >= (Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>(-AngleUtil<Real>::HALF_PI, (Real)0.0, -MathUtil::Atan2(m[0][2], m[0][0]));
	}
	else if (m[2][1] <= -(Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>(AngleUtil<Real>::HALF_PI, (Real)0.0, MathUtil::Atan2(m[0][2], m[0][0]));
	}
	else
	{
		return Vector3<Real>(MathUtil::Asin(-m[2][1]/scaling[2]), MathUtil::Atan2(m[2][0], m[2][2]), MathUtil::Atan2(m[0][1]*scaling[1], m[1][1]*scaling[0]));
	}
}

template <typename Real>
Vector3<Real> Matrix4<Real>::GetEulerZYX() const
{
	Vector3<Real> scaling = GetScaling();
	// rot =  cy*cz           cx*sz+sx*sy*cz  sx*sz-cx*cz*sy
	//       -cy*sz           cx*cz-sx*sy*sz  sx*cz+cx*sy*sz
	//        sy             -sx*cy           cx*cy
	if (m[2][0] >= (Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>((Real)0.0, AngleUtil<Real>::HALF_PI, MathUtil::MathUtil::Atan2(m[0][1], -m[0][2]));
	}
	else if (m[2][0] <= -(Real)1.0)
	{
		// WARNING. Maybe not unique.
		return Vector3<Real>((Real)0.0, -AngleUtil<Real>::HALF_PI, MathUtil::Atan2(m[0][1], m[0][2]));
	}
	else
	{
		return Vector3<Real>(MathUtil::MathUtil::Atan2(-m[2][1], m[2][2]), MathUtil::MathUtil::Asin(m[2][0]/scaling[2]), MathUtil::MathUtil::Atan2(-m[1][0]*scaling[0], m[0][0]*scaling[1]));
	}
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::SetScaling(const Vector3<Real>& scaling)
{
	Vector3<Real> oldScaling = GetScaling();
	return Scale(Vector3<Real>(scaling.x / oldScaling.x, scaling.y / oldScaling.y, scaling.z / oldScaling.z));
}

template <typename Real>
Vector3<Real> Matrix4<Real>::GetScaling() const
{
	Vector3<Real> scaling;
	for(int i=0;i<3;i++)
	{
		scaling[i] = MathUtil::Sqrt((*this)[i].Dot((*this)[i]));
	}
	return scaling;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::Transpose(void) const
{
	Matrix4<Real> result;
	for (int i=0; i<4; ++i)
	{
		for (int j=0; j<4; ++j)
		{
			result.m[i][j] = m[j][i];
		}
	}
	return result;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::Inverse(void) const
{
	Real fA0 = m2[ 0]*m2[ 5] - m2[ 1]*m2[ 4];
	Real fA1 = m2[ 0]*m2[ 6] - m2[ 2]*m2[ 4];
	Real fA2 = m2[ 0]*m2[ 7] - m2[ 3]*m2[ 4];
	Real fA3 = m2[ 1]*m2[ 6] - m2[ 2]*m2[ 5];
	Real fA4 = m2[ 1]*m2[ 7] - m2[ 3]*m2[ 5];
	Real fA5 = m2[ 2]*m2[ 7] - m2[ 3]*m2[ 6];
	Real fB0 = m2[ 8]*m2[13] - m2[ 9]*m2[12];
	Real fB1 = m2[ 8]*m2[14] - m2[10]*m2[12];
	Real fB2 = m2[ 8]*m2[15] - m2[11]*m2[12];
	Real fB3 = m2[ 9]*m2[14] - m2[10]*m2[13];
	Real fB4 = m2[ 9]*m2[15] - m2[11]*m2[13];
	Real fB5 = m2[10]*m2[15] - m2[11]*m2[14];

	Real fDet = fA0*fB5 - fA1*fB4 + fA2*fB3 + fA3*fB2 - fA4*fB1 + fA5*fB0;
	if (MathUtil::Abs(fDet) <= (Real)0.0)
	{
		K_MATH_LOG("Singular matrix has no inversion.");
		return Matrix4<Real>::ZERO;
	}

	Matrix4<Real> result;
	result.m[0][0] = + m2[ 5]*fB5 - m2[ 6]*fB4 + m2[ 7]*fB3;
	result.m[1][0] = - m2[ 4]*fB5 + m2[ 6]*fB2 - m2[ 7]*fB1;
	result.m[2][0] = + m2[ 4]*fB4 - m2[ 5]*fB2 + m2[ 7]*fB0;
	result.m[3][0] = - m2[ 4]*fB3 + m2[ 5]*fB1 - m2[ 6]*fB0;
	result.m[0][1] = - m2[ 1]*fB5 + m2[ 2]*fB4 - m2[ 3]*fB3;
	result.m[1][1] = + m2[ 0]*fB5 - m2[ 2]*fB2 + m2[ 3]*fB1;
	result.m[2][1] = - m2[ 0]*fB4 + m2[ 1]*fB2 - m2[ 3]*fB0;
	result.m[3][1] = + m2[ 0]*fB3 - m2[ 1]*fB1 + m2[ 2]*fB0;
	result.m[0][2] = + m2[13]*fA5 - m2[14]*fA4 + m2[15]*fA3;
	result.m[1][2] = - m2[12]*fA5 + m2[14]*fA2 - m2[15]*fA1;
	result.m[2][2] = + m2[12]*fA4 - m2[13]*fA2 + m2[15]*fA0;
	result.m[3][2] = - m2[12]*fA3 + m2[13]*fA1 - m2[14]*fA0;
	result.m[0][3] = - m2[ 9]*fA5 + m2[10]*fA4 - m2[11]*fA3;
	result.m[1][3] = + m2[ 8]*fA5 - m2[10]*fA2 + m2[11]*fA1;
	result.m[2][3] = - m2[ 8]*fA4 + m2[ 9]*fA2 - m2[11]*fA0;
	result.m[3][3] = + m2[ 8]*fA3 - m2[ 9]*fA1 + m2[10]*fA0;
	return (result *= 1.0f / fDet);
}

template <typename Real>
Real Matrix4<Real>::Determinant(void) const
{
	Real fA0 = m2[ 0]*m2[ 5] - m2[ 1]*m2[ 4];
	Real fA1 = m2[ 0]*m2[ 6] - m2[ 2]*m2[ 4];
	Real fA2 = m2[ 0]*m2[ 7] - m2[ 3]*m2[ 4];
	Real fA3 = m2[ 1]*m2[ 6] - m2[ 2]*m2[ 5];
	Real fA4 = m2[ 1]*m2[ 7] - m2[ 3]*m2[ 5];
	Real fA5 = m2[ 2]*m2[ 7] - m2[ 3]*m2[ 6];
	Real fB0 = m2[ 8]*m2[13] - m2[ 9]*m2[12];
	Real fB1 = m2[ 8]*m2[14] - m2[10]*m2[12];
	Real fB2 = m2[ 8]*m2[15] - m2[11]*m2[12];
	Real fB3 = m2[ 9]*m2[14] - m2[10]*m2[13];
	Real fB4 = m2[ 9]*m2[15] - m2[11]*m2[13];
	Real fB5 = m2[10]*m2[15] - m2[11]*m2[14];
	Real fDet = fA0*fB5 - fA1*fB4 + fA2*fB3 + fA3*fB2 - fA4*fB1 + fA5*fB0;
	return fDet;
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::ClearScaling(void)
{
	Vector3f scaling = GetScaling();
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			m[i][j] /= scaling[i];
		}
	}
	return *this;
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::ClearRotation(void)
{
	Vector3f scaling = GetScaling();
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			m[i][j] = (i==j?scaling[i]:0.0f);
		}
	}
	return *this;
}

template <typename Real>
inline Matrix4<Real>& Matrix4<Real>::ClearTranslation(void)
{
	m[3][0] = m[3][1] = m[3][2] = 0.0f;
	return *this;
}

template <typename Real>
Matrix4<Real>& K::MATH::Matrix4<Real>::Scale( const Vector3<Real>& v )
{
	for(int i=0;i<3;i++)
	{
		(*this)[i] *= v[i];
	}
	return *this;
}

template <typename Real>
Matrix4<Real>& Matrix4<Real>::Translate( const Vector3<Real>& v )
{
	m[3][0] += v.x;
	m[3][1] += v.y;
	m[3][2] += v.z;
	return *this;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::MakeScaling( const Vector3<Real>& v )
{
	Matrix4<Real> res(IDENTITY);
	m[0][0] = v.x;
	m[1][1] = v.y;
	m[2][2] = v.z;
	return res;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::MakeScaling( Real x, Real y, Real z )
{
	Matrix4<Real> res(IDENTITY);
	m[0][0] = x;
	m[1][1] = y;
	m[2][2] = z;
	return res;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::Reflect(const Real* p)
{
    // 求相对于平面的对称矩阵，算法如下。
    // -2 * P.a * P.a + 1  -2 * P.b * P.a      -2 * P.c * P.a        0
    // -2 * P.a * P.b      -2 * P.b * P.b + 1  -2 * P.c * P.b        0
    // -2 * P.a * P.c      -2 * P.b * P.c      -2 * P.c * P.c + 1    0
    // -2 * P.a * P.d      -2 * P.b * P.d      -2 * P.c * P.d        1

    Real f2A2 = p[0] * p[0] * 2.0f;
    Real f2B2 = p[1] * p[1] * 2.0f;
    Real f2C2 = p[2] * p[2] * 2.0f;
    Real f2AB = p[0] * p[1] * 2.0f;
    Real f2BC = p[1] * p[2] * 2.0f;
    Real f2AC = p[0] * p[2] * 2.0f;
    Real f2AD = p[0] * p[3] * 2.0f;
    Real f2BD = p[1] * p[3] * 2.0f;
    Real f2CD = p[2] * p[3] * 2.0f;

    Matrix4<Real> result;
    result.SetIdentity();

    result.m[0][0] -= f2A2;
    result.m[0][1] -= f2AB;
    result.m[0][2] -= f2AC;

    result.m[1][0] -= f2AB;
    result.m[1][1] -= f2B2;
    result.m[1][2] -= f2BC;

    result.m[2][0] -= f2AC;
    result.m[2][1] -= f2BC;
    result.m[2][2] -= f2C2;

    result.m[3][0] -= f2AD;
    result.m[3][1] -= f2BD;
    result.m[3][2] -= f2CD;

    return result;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::LookAtLH(const Vector3<Real>& eye, const Vector3<Real>& at, const Vector3<Real>& up)
{
    //  zaxis = normal(at - eye)
    //	xaxis = normal(cross(up, zaxis))
    //	yaxis = cross(zaxis, xaxis)
    //
    //	xaxis.x           yaxis.x           zaxis.x           0
    //	xaxis.y           yaxis.y           zaxis.y           0
    //	xaxis.z           yaxis.z           zaxis.z           0
    //	-dot(xaxis, eye)  -dot(yaxis, eye)  -dot(zaxis, eye)  1

    Vector3<Real> zAxis = Normalize(at - eye);
    Vector3<Real> xAxis = Normalize(up.Cross(zAxis));
    Vector3<Real> yAxis = zAxis.Cross(xAxis);

    Matrix4<Real> result;
    result.m[0][0] = xAxis.x;
    result.m[0][1] = yAxis.x;
    result.m[0][2] = zAxis.x;
    result.m[0][3] = 0.0f;
    result.m[1][0] = xAxis.y;
    result.m[1][1] = yAxis.y;
    result.m[1][2] = zAxis.y;
    result.m[1][3] = 0.0f;
    result.m[2][0] = xAxis.z;
    result.m[2][1] = yAxis.z;
    result.m[2][2] = zAxis.z;
    result.m[2][3] = 0.0f;
    result.m[3][0] = -xAxis.Dot(eye);
    result.m[3][1] = -yAxis.Dot(eye);
    result.m[3][2] = -zAxis.Dot(eye);
    result.m[3][3] = 1.0f;
    return result;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::LookAtRH(const Vector3<Real>& eye, const Vector3<Real>& at, const Vector3<Real>& up)
{
    //  zaxis = normal(Eye - At)
    //	xaxis = normal(cross(up, zaxis))
    //	yaxis = cross(zaxis, xaxis)
    //
    //	xaxis.x           yaxis.x           zaxis.x           0
    //	xaxis.y           yaxis.y           zaxis.y           0
    //	xaxis.z           yaxis.z           zaxis.z           0
    //	-dot(xaxis, eye)  -dot(yaxis, eye)  -dot(zaxis, eye)  1

    Vector3<Real> zAxis = Normalize(eye - at);
    Vector3<Real> xAxis = Normalize(up.Cross(zAxis));
    Vector3<Real> yAxis = zAxis.Cross(xAxis);

    Matrix4<Real> result;
    result.m[0][0] = xAxis.x;
    result.m[0][1] = yAxis.x;
    result.m[0][2] = zAxis.x;
    result.m[0][3] = 0.0f;
    result.m[1][0] = xAxis.y;
    result.m[1][1] = yAxis.y;
    result.m[1][2] = zAxis.y;
    result.m[1][3] = 0.0f;
    result.m[2][0] = xAxis.z;
    result.m[2][1] = yAxis.z;
    result.m[2][2] = zAxis.z;
    result.m[2][3] = 0.0f;
    result.m[3][0] = -xAxis.Dot(eye);
    result.m[3][1] = -yAxis.Dot(eye);
    result.m[3][2] = -zAxis.Dot(eye);
    result.m[3][3] = 1.0f;
    return result;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::PerspectiveLH(Real w, Real h, Real zn, Real zf)
{
    // 2*zn/w  0       0              0
    // 0       2*zn/h  0              0
    // 0       0       zf/(zf-zn)     1
    // 0       0       zn*zf/(zn-zf)  0

    Matrix4<Real> result;
    result.m[0][0] = 2 * zn / w;
    result.m[0][1] = 0.0;
    result.m[0][2] = 0.0;
    result.m[0][3] = 0.0;

    result.m[1][0] = 0.0;
    result.m[1][1] = 2 * zn / h;
    result.m[1][2] = 0.0;
    result.m[1][3] = 0.0;

    result.m[2][0] = 0.0;
    result.m[2][1] = 0.0;
    result.m[2][2] = zf / (zf - zn);
    result.m[2][3] = Real(1.0);

    result.m[3][0] = 0.0;
    result.m[3][1] = 0.0;
    result.m[3][2] = zn * zf / (zn - zf);
    result.m[3][3] = 0.0;
    return result;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::PerspectiveRH(Real w, Real h, Real zn, Real zf)
{
    // 2*zn/w  0       0              0
    // 0       2*zn/h  0              0
    // 0       0       zf/(zn-zf)    -1
    // 0       0       zn*zf/(zn-zf)  0

    Matrix4<Real> result;
    result.m[0][0] = 2 * zn / w;
    result.m[0][1] = 0.0;
    result.m[0][2] = 0.0;
    result.m[0][3] = 0.0;

    result.m[1][0] = 0.0;
    result.m[1][1] = 2 * zn / h;
    result.m[1][2] = 0.0;
    result.m[1][3] = 0.0;

    result.m[2][0] = 0.0;
    result.m[2][1] = 0.0;
    result.m[2][2] = zf / (zn - zf);
    result.m[2][3] = Real(-1.0);

    result.m[3][0] = 0.0;
    result.m[3][1] = 0.0;
    result.m[3][2] = zn * zf / (zn - zf);
    result.m[3][3] = 0.0;
    return result;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::PerspectiveFovLH(Real fovY, Real aspect, Real zn, Real zf)
{
    //  xScale     0          0              0
    //	0          yScale     0              0
    //	0          0       zf/(zf-zn)        1
    //	0          0       -zn*zf/(zf-zn)    0
    //where:
    // yScale = cot(fovY/2)
    // xScale = yScale / aspect

    Real yScale = 1.0f / tan(fovY * 0.5f);
    Real xScale = yScale / aspect;

    Matrix4<Real> result;
    result.SetIdentity();
    result.m[0][0] = xScale;
    result.m[1][1] = yScale;
    result.m[2][2] = zf/(zf-zn);
    result.m[2][3] = 1.0f;
    result.m[3][2] = zn*zf/(zn-zf);
    result.m[3][3] = 0.0f;
    return result;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::PerspectiveFovRH(Real fovY, Real aspect, Real zn, Real zf)
{
    //  xScale     0          0              0
    //	0          yScale     0              0
    //	0          0        zf/(zn-zf)      -1
    //	0          0        zn*zf/(zn-zf)    0
    //where:
    // yScale = cot(fovY/2)
    // xScale = yScale / aspect

    Real yScale = 1.0f / tan(fovY * 0.5f);
    Real xScale = yScale / aspect;

    Matrix4<Real> result;
    result.SetIdentity();
    result.m[0][0] = xScale;
    result.m[1][1] = yScale;
    result.m[2][2] = zf/(zn-zf);
    result.m[2][3] = -1.0f;
    result.m[3][2] = zn*zf/(zn-zf);
    result.m[3][3] = 0.0f;
    return result;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::OrthoLH(Real w, Real h, Real zn, Real zf)
{
    // 2/w  0    0           0
    // 0    2/h  0           0
    // 0    0    1/(zf-zn)   0
    // 0    0    zn/(zn-zf)  1

    Matrix4<Real> result;
    result.SetIdentity();
    result.m[0][0] = 2 / w;
    result.m[1][1] = 2 / h;
    result.m[2][2] = 1 / (zf - zn);
    result.m[3][2] = zn / (zn - zf);
    return result;
}

template <typename Real>
Matrix4<Real> Matrix4<Real>::OrthoRH(Real w, Real h, Real zn, Real zf)
{
    // 2/w  0    0           0
    // 0    2/h  0           0
    // 0    0    1/(zn-zf)   0
    // 0    0    zn/(zn-zf)  l

    Matrix4<Real> result;
    result.SetIdentity();
    result.m[0][0] = 2 / w;
    result.m[1][1] = 2 / h;
    result.m[2][2] = 1 / (zn - zf);
    result.m[3][2] = zn / (zn - zf);
    return result;
}

/// \brief 三维向量(扩展为四维点 w=1)和四维矩阵相乘，得到一个新三维向量，并对其齐次化
template <typename Real>
Vector3<Real> operator * (const Vector3<Real>& v, const Matrix4<Real>& mtx)
{
	Vector3<Real> result;
	result.x = v.x * mtx.m[0][0] + v.y * mtx.m[1][0] + v.z * mtx.m[2][0] + mtx.m[3][0];
	result.y = v.x * mtx.m[0][1] + v.y * mtx.m[1][1] + v.z * mtx.m[2][1] + mtx.m[3][1];
	result.z = v.x * mtx.m[0][2] + v.y * mtx.m[1][2] + v.z * mtx.m[2][2] + mtx.m[3][2];
	Real w   = v.x * mtx.m[0][3] + v.y * mtx.m[1][3] + v.z * mtx.m[2][3] + mtx.m[3][3];
	if(w != Real(0.0))
	{
		Real invW= Real(1.0f) / w;
		result.x *= invW;
		result.y *= invW;
		result.z *= invW;
	}
	return result;
}

/// \brief 四维向量和四维矩阵相乘，得到一个新四维向量。
template <typename Real>
Vector4<Real> operator * (const Vector4<Real>& v, const Matrix4<Real>& mtx)
{
	return Vector4<Real>(
		v.x * mtx.m[0][0] + v.y * mtx.m[1][0] + v.z * mtx.m[2][0] + v.w * mtx.m[3][0],
		v.x * mtx.m[0][1] + v.y * mtx.m[1][1] + v.z * mtx.m[2][1] + v.w * mtx.m[3][1],
		v.x * mtx.m[0][2] + v.y * mtx.m[1][2] + v.z * mtx.m[2][2] + v.w * mtx.m[3][2],
		v.x * mtx.m[0][3] + v.y * mtx.m[1][3] + v.z * mtx.m[2][3] + v.w * mtx.m[3][3]
		);
}

/// \brief 矩阵数乘运算
template <typename Real>
Matrix4<Real> operator* (Real fScalar, const Matrix4<Real>& mtx)
{
	return mtx * fScalar;
}
