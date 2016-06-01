/* Start Header =========================================================================
File Name:		Vector.h
Purpose:		Implements a 3-Vector.
Language:		C++, Visual Studio 2013 compiler
Author:			Hew Jun-Wei
== End Header =========================================================================*/

// Header Guard:
//---------------------------------------------------------------------------------------
#ifndef VECTOR_H
#define VECTOR_H


/****************************************************************************************
Class: Vector3

Implements a 3d vector.
****************************************************************************************/
struct Vector3
{
	Vector3();
	~Vector3();
	Vector3(const Vector3& v);
	Vector3(float x, float y, float z);

	Vector3		Normalize()						const;
	Vector3&	NormalizeThis();
	
	float		Length()						const;
	float		SquaredLength()					const;

	Vector3		operator-()						const;
	Vector3&	operator=(const Vector3& rhs);

	Vector3		operator+(const Vector3& rhs)	const;
	Vector3		operator-(const Vector3& rhs)	const;
	Vector3&	operator+=(const Vector3& rhs);
	Vector3&	operator-=(const Vector3& rhs);

	float		operator*(const Vector3& rhs)	const;	// dot product
	Vector3		operator%(const Vector3& rhs)	const;	// cross product

	bool		operator==(const Vector3& rhs)	const;

	Vector3		operator/(const float rhs)		const;
	Vector3&	operator*=(const float rhs);
	Vector3&	operator/=(const float rhs);

	union
	{
		struct
		{
			float x, y, z;
		};
		float data[3];
	};
};

Vector3 operator*(const float lhs, const Vector3& rhs);

typedef		Vector3		Vec3;

#endif