/* Start Header =========================================================================
File Name:		Matrix.h
Purpose:		Implements 3x3 matrix.
Language:		C++, Visual Studio 2013 compiler
Author:			Hew Jun-Wei
== End Header =========================================================================*/

// Header Guard:
//---------------------------------------------------------------------------------------
#ifndef MATRIX_H
#define MATRIX_H

// Includes:
//---------------------------------------------------------------------------------------
#include "Vector.h"

/****************************************************************************************
Class: Matrix33

Implements a 3x3 matrix. Vectors and entries are stored COLUMN-MAJOR.
****************************************************************************************/
struct Matrix33
{
	// Constructors/Destructors:
	//-------------------------------------------------
	Matrix33();
	~Matrix33();

	Matrix33(float m00, float m01, float m02,
			 float m10, float m11, float m12,
			 float m20, float m21, float m22);
	Matrix33(Vector3 a1, Vector3 a2, Vector3 a3);
	Matrix33(const Matrix33& rhs);

	float& ColRow(unsigned column, unsigned row);
	float  ColRow(unsigned column, unsigned row) const;
	float& operator[](unsigned index);
	float  operator[](unsigned index) const;
	float& operator()(unsigned column, unsigned row);
	float  operator()(unsigned column, unsigned row) const;

	Matrix33 Transpose() const;
	Matrix33& TransposeThis();

	Matrix33 Inverse() const;
	Matrix33& InverseThis();

	float Determinant() const;

	Matrix33 operator+(const Matrix33& rhs) const;
	Matrix33 operator-(const Matrix33& rhs) const;
	Matrix33 operator*(const Matrix33& rhs) const;

	Matrix33& operator+=(const Matrix33& rhs);
	Matrix33& operator-=(const Matrix33& rhs);

	Vector3 operator*(const Vector3& rhs) const;

	Matrix33 operator*(const float rhs) const;
	Matrix33 operator/(const float rhs) const;
	Matrix33& operator*=(const float rhs);
	Matrix33& operator/=(const float rhs);

	static Matrix33 Identity();
	static Matrix33 Zero();
	static Matrix33 Scale(float x, float y, float z = 1.0f);
	static Matrix33 Rotate(float x, float y, float z, float angle);
	static Matrix33 Rotate(Vector3 axis, float angle);

	union
	{
		struct
		{
			// matrix data is stored COLUMN-MAJOR
			// naming convention is m[col][row]
			float	m00, m01, m02,
					m10, m11, m12,
					m20, m21, m22;
		};
		struct
		{
			Vector3 a1, a2, a3;
		};
		float data[9];
	};
};

Matrix33 operator*(float lhs, const Matrix33& rhs);

#endif