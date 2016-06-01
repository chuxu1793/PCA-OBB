/* Start Header =========================================================================
File Name:		OBB.h
Purpose:		Implements an oriented bounding box.
Language:		C++, Visual Studio 2013 compiler
Author:			Hew Jun-Wei
== End Header =========================================================================*/

// Header Guard:
//---------------------------------------------------------------------------------------
#ifndef OBB_H
#define OBB_H

// Includes:
//---------------------------------------------------------------------------------------
#include "Matrix.h"

/****************************************************************************************
****************************************************************************************/
class OBB
{
public:
	OBB::OBB() : position(0.0f, 0.0f, 0.0f), halfExtents(1.0f, 1.0f, 1.0f)
	{
		orientation = Matrix33::Identity();
	}

	~OBB() {}

	Vector3		position;
	Matrix33	orientation;
	Vector3		halfExtents;	// half-extents along local x, y, z axes
};

#endif