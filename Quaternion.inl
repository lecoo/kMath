/////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) 2009-2012 DaiMingzhuang. All rights reserved.
/// \file    Quaternion.inl
/// \brief   模板Quaternion的inline函数的实现
/// \author  DaiMingzhuang
/// \version 
/// \date    2009-03
/////////////////////////////////////////////////////////////////////////////////
template <typename Real>
const Quaternion<Real> Quaternion<Real>::IDENTITY(0.0f, 0.0f, 0.0f,1.0f);

template <typename Real>
Quaternion<Real> Quaternion<Real>::operator* (const Quaternion<Real>& qtn) const
{
	Quaternion<Real> result;
	result.x = w*qtn.x + x*qtn.w + y*qtn.z - z*qtn.y;
	result.y = w*qtn.y + y*qtn.w + z*qtn.x - x*qtn.z;
	result.z = w*qtn.z + z*qtn.w + x*qtn.y - y*qtn.x;
	result.w = w*qtn.w - x*qtn.x - y*qtn.y - z*qtn.z;
	return result;
}

template <typename Real>
Quaternion<Real>& Quaternion<Real>::SetIdentity(void)
{
	memcpy(this, &IDENTITY, sizeof(IDENTITY));
	return *this;
}

template <typename Real>
Quaternion<Real> Quaternion<Real>::Conjugate(void) const
{
    return Quaternion<Real>(-x, -y, -z, w);
}

template <typename Real>
Quaternion<Real> Quaternion<Real>::Inverse(void) const
{
	Quaternion<Real> result;

	Real length = x*x + y*y + z*z + w*w;

	if (length > 0.0f)
	{
		Real invLength = 1.0f / length;
		result.x = -x * invLength;
		result.y = -y * invLength;
		result.z = -z * invLength;
		result.w =  w * invLength;
	}
	else
	{
		result.x = 0;
		result.y = 0;
		result.z = 0;
		result.w = 0;
	}

	return result;
}

template <typename Real>
Quaternion<Real> Quaternion<Real>::Slerp(const Quaternion<Real>& qtn2, Real lambda)
{
	return MathUtil::Slerp(*this, qtn2, lambda);
}

template <typename Real>
void Quaternion<Real>::FromRotationMatrix(const Matrix3<Real>& mtx)
{
	Real fTrace = mtx.m[0][0] + mtx.m[1][1] + mtx.m[2][2];

	if (fTrace > 0.0f)
	{
		Real fRoot = MathUtil::Sqrt(fTrace + 1.0f);
		w = fRoot * 0.5f;

		fRoot = 0.5f / fRoot;
        x = (mtx.m[1][2] - mtx.m[2][1]) * fRoot;
        y = (mtx.m[2][0] - mtx.m[0][2]) * fRoot;
        z = (mtx.m[0][1] - mtx.m[1][0]) * fRoot;
	}
	else
	{
		if (mtx.m[0][0] >= mtx.m[1][1] && mtx.m[0][0] >= mtx.m[2][2])
		{
			Real fRoot = MathUtil::Sqrt(mtx.m[0][0] - mtx.m[1][1] - mtx.m[2][2] + 1.0f);
			x = fRoot * 0.5f;

			fRoot = 0.5f / fRoot;
            w = (mtx.m[1][2] - mtx.m[2][1]) * fRoot;
            y = (mtx.m[0][1] + mtx.m[1][0]) * fRoot;
            z = (mtx.m[0][2] + mtx.m[2][0]) * fRoot;
		}
		else if (mtx.m[1][1] >= mtx.m[0][0] && mtx.m[1][1] >= mtx.m[2][2])
		{
			Real fRoot = MathUtil::Sqrt(mtx.m[1][1] - mtx.m[2][2] - mtx.m[0][0] + 1.0f);
			y = fRoot * 0.5f;

			fRoot = 0.5f / fRoot;
            w = (mtx.m[2][0] - mtx.m[0][2]) * fRoot;
            z = (mtx.m[1][2] + mtx.m[2][1]) * fRoot;
            x = (mtx.m[1][0] + mtx.m[0][1]) * fRoot;
		}
		else
		{
			Real fRoot = MathUtil::Sqrt(mtx.m[2][2] - mtx.m[0][0] - mtx.m[1][1] + 1.0f);
			z = fRoot * 0.5f;

			fRoot = 0.5f / fRoot;
            w = (mtx.m[0][1] - mtx.m[1][0]) * fRoot;
            x = (mtx.m[2][0] + mtx.m[0][2]) * fRoot;
            y = (mtx.m[2][1] + mtx.m[1][2]) * fRoot;
		}
	}
}

template <typename Real>
void Quaternion<Real>::ToRotationMatrix(Matrix3<Real>& mtx) const
{
	Real fTx  = x * 2.0f;
	Real fTy  = y * 2.0f;
	Real fTz  = z * 2.0f;
	Real fTwx = fTx * w;
	Real fTwy = fTy * w;
	Real fTwz = fTz * w;
	Real fTxx = fTx * x;
	Real fTxy = fTy * x;
	Real fTxz = fTz * x;
	Real fTyy = fTy * y;
	Real fTyz = fTz * y;
	Real fTzz = fTz * z;

	mtx.m[0][0] = 1.0f - (fTyy+fTzz);
	mtx.m[1][0] = fTxy - fTwz;
	mtx.m[2][0] = fTxz + fTwy;
	mtx.m[0][1] = fTxy + fTwz;
	mtx.m[1][1] = 1.0f - (fTxx+fTzz);
	mtx.m[2][1] = fTyz - fTwx;
	mtx.m[0][2] = fTxz - fTwy;
	mtx.m[1][2] = fTyz + fTwx;
	mtx.m[2][2] = 1.0f - (fTxx+fTyy);
}

template <typename Real>
void Quaternion<Real>::FromAxisAngle(const Vector3<Real>& axis, Real ang)
{
	Vector3<Real> axis_one = Normalize(axis);

	Real fHalfAngle = ang * 0.5f;

	Real fSin = MathUtil::Sin(fHalfAngle);
	Real fCos = MathUtil::Cos(fHalfAngle);

	w = fCos;
	x = axis_one.x * fSin;
	y = axis_one.y * fSin;
	z = axis_one.z * fSin;
}

template <typename Real>
bool Quaternion<Real>::ToAxisAngle(Vector3<Real>& axis, Real& ang) const
{
	Real fSqrLength = x*x + y*y + z*z;
	if (fSqrLength > MathUtil::FLOAT_EPS)
	{
		ang = MathUtil::Acos(w) * 2.0f;

		axis = Vector3<Real>(x, y, z);
		axis *= MathUtil::Sqrt(1.0f / fSqrLength);
		return true;
	}
	else
	{
		// angle is 0 (mod 2*pi), so any axis will do
		ang = (Real)0.0;
		axis = Vector3<Real>(1.0f, 0.0f, 0.0f);
		return false;
	}
}

template <typename Real>
Vector3<Real> Quaternion<Real>::GetAxis() const
{
	return Vector3<Real>(x, y, z);
}

template <typename Real>
Real Quaternion<Real>::GetAngle() const
{
	return MathUtil::Acos(w) * 2.0f;
}

template <typename Real>
void Quaternion<Real>::FromEulerAnglesZXY(const Vector3<Real>& angles)
{
	Real halfYaw = angles.y * Real(0.5);
	Real halfPitch = angles.x * Real(0.5);
	Real halfRoll = angles.z * Real(0.5);
	Real cosYaw = MathUtil::Cos(halfYaw);
	Real sinYaw = MathUtil::Sin(halfYaw);
	Real cosPitch = MathUtil::Cos(halfPitch);
	Real sinPitch = MathUtil::Sin(halfPitch);
	Real cosRoll = MathUtil::Cos(halfRoll);
	Real sinRoll = MathUtil::Sin(halfRoll);
	Set(cosRoll * sinPitch * cosYaw + sinRoll * cosPitch * sinYaw,
		cosRoll * cosPitch * sinYaw - sinRoll * sinPitch * cosYaw,
		sinRoll * cosPitch * cosYaw - cosRoll * sinPitch * sinYaw,
		cosRoll * cosPitch * cosYaw + sinRoll * sinPitch * sinYaw);
}

template <typename Real>
Vector3<Real> Quaternion<Real>::GetEulerZXY() const
{
	Matrix3<Real> mtx;
	ToRotationMatrix(mtx);
	return mtx.GetEulerZXY();
}

template <typename Real>
Real Quaternion<Real>::Dot(const Quaternion& q) const
{
	return x*q.x + y*q.y + z*q.z + w*q.w;
}

template <typename Real>
Real Quaternion<Real>::Length() const
{
	return MathUtil::Sqrt(x*x+y*y+z*z+w*w);
}

template <typename Real>
Quaternion<Real> operator * (const Quaternion<Real>& q, const Vector3<Real>& v)
{
	return Quaternion<Real>( q.w * v.x + q.y * v.z - q.z * v.y,
		q.w * v.y + q.z * v.x - q.x * v.z,
		q.w * v.z + q.x * v.y - q.y * v.x,
		-q.x * v.x - q.y * v.y - q.z * v.z); 
}

template <typename Real>
Quaternion<Real> operator * (const Vector3<Real>& v, const Quaternion<Real>& q)
{
	return Quaternion<Real>( v.x * q.w + v.y * q.z - v.z * q.y,
		v.y * q.w + v.z * q.x - v.x * q.z,
		v.z * q.w + v.x * q.y - v.y * q.x,
		-v.x * q.x - v.y * q.y - v.z * q.z); 
}

template<typename Real>
Real Dot(const Quaternion<Real>& q1, const Quaternion<Real>& q2)
{
	return q1.Dot(q2);
}

template<typename Real>
Real Length(const Quaternion<Real>& q)
{
	return q.Length();
}

template<typename Real>
Quaternion<Real> Normalize(const Quaternion<Real>& q)
{
	return q/q.Length();
}
