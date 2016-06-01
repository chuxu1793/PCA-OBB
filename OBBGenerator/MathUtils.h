/* Start Header =========================================================================
File Name:		MathUtils.h
Purpose:		Helpful math functions and constants.
Language:		C++, Visual Studio 2013 compiler
Author:			Hew Jun-Wei
== End Header =========================================================================*/

// Header Guard:
//---------------------------------------------------------------------------------------
#ifndef MATHUTILS_H
#define MATHUTILS_H

const float INV_SQRT_TWO	= 0.707106781f;

template <typename T>
inline T Max(const T x, const T y)
{
	return (x < y) ? y : x;
}

template <typename T>
inline T Min(const T x, const T y)
{
	return (x > y) ? y : x;
}

#endif