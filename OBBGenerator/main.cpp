/* Start Header =========================================================================
File Name:		main.cpp
Purpose:		Generates an oriented bounding box (OBB) from a provided set of vertices
				via principal-component analysis (PCA).
Language:		C++, Visual Studio 2013 compiler
Author:			Hew Jun-Wei
== End Header =========================================================================*/

// Includes:
//---------------------------------------------------------------------------------------
#include "OBB.h"
#include "MathUtils.h"
#include <vector>
#include <iostream>
#include <sstream>

// Function prototypes:
//---------------------------------------------------------------------------------------
OBB				GenerateOBB(const Vector3* pVertices, const int numVertices);
Matrix33		ComputeCovarianceMatrix(const Vector3* pVertices, const int numVertices);
void			GrahmSchmidt(Vector3& v0, Vector3& v1, Vector3& v2);
void			JacobiSolver(Matrix33 m, float eValues[3], Vector3 eVectors[3]);
std::ostream&	operator<<(std::ostream& os, const Vector3& v);


/****************************************************************************************
Function: main

Entry point.
****************************************************************************************/
int main(int argc, char* argv[])
{
	std::vector<Vector3> vertices;
	std::string userInput;

	std::cout << "Enter a list of Vec3s, one triple per line, as\n\tx y z\n";
	std::cout << "Press ENTER on a blank vector to finish.\r\n";

	do 
	{
		std::cout << "p[" << vertices.size() << "] = ";

		getline(std::cin, userInput);
		std::stringstream ss(userInput);
		Vector3 newVec;
		ss >> newVec.x >> newVec.y >> newVec.z;

		if (!userInput.empty())
		{
			if (!ss)
			{
				std::cout << "Invalid Vec3. Try again.\r\n";
			}
			else
			{
				vertices.push_back(newVec);
			}
		}
	} while (!userInput.empty());

	if (vertices.empty())
	{
		std::cout << "\r\nUsing demo set:\r\n";

		vertices.push_back(Vector3(3.7f, 1.7f, 0.0f));
		vertices.push_back(Vector3(4.1f, 3.8f, 0.0f));
		vertices.push_back(Vector3(4.7f, 2.9f, 0.0f));
		vertices.push_back(Vector3(5.2f, 2.8f, 0.0f));
		vertices.push_back(Vector3(6.0f, 4.0f, 0.0f));
		vertices.push_back(Vector3(6.3f, 3.6f, 0.0f));
		vertices.push_back(Vector3(9.7f, 6.3f, 0.0f));
		vertices.push_back(Vector3(10.0f, 4.9f, 0.0f));
		vertices.push_back(Vector3(11.0f, 3.6f, 0.0f));
		vertices.push_back(Vector3(12.5f, 6.4f, 0.0f));

		for (size_t i = 0; i < vertices.size(); ++i)
		{
			std::cout << "p[" << i << "] = " << vertices[i] << ")\r\n";
		}
	}

	OBB obb = GenerateOBB(vertices.data(), vertices.size());

	std::cout << "------------------------------------------------------------\r\n";
	std::cout << "Generated OBB:\r\n";
	std::cout << "Center: " << obb.position << "\r\n";
	std::cout << "Axis 1: " << obb.halfExtents.x << " * " << obb.orientation.a1 << "\r\n";
	std::cout << "Axis 2: " << obb.halfExtents.y << " * " << obb.orientation.a2 << "\r\n";
	std::cout << "Axis 3: " << obb.halfExtents.z << " * " << obb.orientation.a3 << "\r\n" << std::endl;

	system("pause");

	return 0;
}


/****************************************************************************************
Function: GenerateOBB

Computes a bounding oriented bounding box (OBB) given the supplied vertices.

pVertex		: A pointer to the first Vector3 in the series.
numVertices : The number of vertices to compute the covariance matrix over.
****************************************************************************************/
OBB GenerateOBB(const Vector3* pVertices, const int numVertices)
{
	OBB			result;
	float		eValue[3];
	Vector3		eVector[3];
	Matrix33	covariance;
	Vector3		axis;

	for (int i = 0; i < numVertices; ++i)
	{
		result.position += pVertices[i];
	}

	result.position /= static_cast<float>(numVertices);

	covariance = ComputeCovarianceMatrix(pVertices, numVertices);

	JacobiSolver(covariance, eValue, eVector);

	float temp;
	Vector3 tempVec;

	// sort to obtain eValue[0] <= eValue[1] <= eValue[2]:
	if (eValue[0] <= eValue[1])
	{
		if (eValue[1] > eValue[2])
		{
			if (eValue[0] < eValue[2])
			{
				temp = eValue[0];
				tempVec = eVector[0];

				eValue[0] = eValue[2];
				eValue[2] = eValue[1];
				eValue[1] = temp;

				eVector[0] = eVector[2];
				eVector[2] = eVector[1];
				eVector[1] = tempVec;
			}
			else
			{
				temp = eValue[1];
				tempVec = eVector[1];

				eValue[1] = eValue[2];
				eValue[2] = temp;

				eVector[1] = eVector[2];
				eVector[2] = tempVec;
			}
		}
		// else do nothing
	}
	else
	{
		if (eValue[0] > eValue[2])
		{
			if (eValue[1] < eValue[2])
			{
				temp = eValue[0];
				tempVec = eVector[0];

				eValue[0] = eValue[1];
				eValue[1] = eValue[2];
				eValue[2] = temp;

				eVector[0] = eVector[1];
				eVector[1] = eVector[2];
				eVector[2] = tempVec;
			}
			else
			{
				temp = eValue[0];
				tempVec = eVector[0];

				eValue[0] = eValue[2];
				eValue[2] = temp;

				eVector[0] = eVector[2];
				eVector[2] = tempVec;
			}
		}
		else
		{
			temp = eValue[0];
			tempVec = eVector[0];

			eValue[0] = eValue[1];
			eValue[1] = temp;

			eVector[0] = eVector[1];
			eVector[1] = tempVec;
		}
	}

	result.orientation.a1 = eVector[2];
	result.orientation.a2 = eVector[1];
	result.orientation.a3 = eVector[0];

	// perform Grahm-Schmidt orthogonalization using the eigenvector corresponding to the
	// largest eigenvalue as the base vector
	GrahmSchmidt(result.orientation.a1, result.orientation.a2, result.orientation.a3);

	// eigenbasis set- now center the OBB in the middle
	float infinity = std::numeric_limits<float>::infinity();

	Vector3 minExtents(infinity, infinity, infinity);
	Vector3 maxExtents(-infinity, -infinity, -infinity);

	for (int index = 0; index < numVertices; ++index)
	{
		Vector3 vec = pVertices[index];
		Vector3 displacement = vec - result.position;

		minExtents.x = Min(minExtents.x, displacement * result.orientation.a1);
		minExtents.y = Min(minExtents.y, displacement * result.orientation.a2);
		minExtents.z = Min(minExtents.z, displacement * result.orientation.a3);

		maxExtents.x = Max(maxExtents.x, displacement * result.orientation.a1);
		maxExtents.y = Max(maxExtents.y, displacement * result.orientation.a2);
		maxExtents.z = Max(maxExtents.z, displacement * result.orientation.a3);
	}

	Vector3 offset = (maxExtents - minExtents) / 2.0f + minExtents;

	result.position += (offset.x * result.orientation.a1) +
					   (offset.y * result.orientation.a2) +
					   (offset.z * result.orientation.a3);

	for (int i = 0; i < 3; ++i)
	{
		result.halfExtents.data[i] = (maxExtents.data[i] - minExtents.data[i]) / 2.0f;
	}

	return result;
}


/****************************************************************************************
Function: ComputeCovarianceMatrix

Computes a covariance matrix given a series of Vector3s.

pVec		: A pointer to the first Vector3 in the series.
numVertices : The number of vertices to compute the covariance matrix over.
****************************************************************************************/
Matrix33 ComputeCovarianceMatrix(const Vector3* pVertices, const int numVertices)
{
	Matrix33 covariance;

	// duplicate the vector array
	Vector3* pVectors = new Vector3[numVertices];

	// compute the average x, y, z values
	Vector3 avg(0.0f, 0.0f, 0.0f);
	for (int i = 0; i < numVertices; ++i)
	{
		pVectors[i] = pVertices[i];
		avg += pVertices[i];
	}
	avg /= static_cast<float>(numVertices);

	for (int i = 0; i < numVertices; ++i)
	{
		pVectors[i] -= avg;
	}

	// compute the covariance (we are computing the lower-triangular entries then using
	// the symmetric property):
	for (int c = 0; c < 3; ++c)
	{
		for (int r = c; r < 3; ++r)
		{
			float& acc = covariance.ColRow(c, r);
			acc = 0.0f;

			// cov(X, Y) = E[(X - x)(Y - y)]
			for (int i = 0; i < numVertices; ++i)
			{
				acc += pVectors[i].data[c] * pVectors[i].data[r];
			}

			acc /= static_cast<float>(numVertices);

			covariance.ColRow(r, c) = acc;	// covariance matrix is symmetric
		}
	}

	delete[] pVectors;
	return covariance;
}


/****************************************************************************************
Function: GrahmSchmidt

Orthonormalizes three given vectors in the order provided.

v0 :	The major axis vector. This vector's orientation will be unchanged.
v1 :	The secondary axis vector.
v2 :	The tertiary axis vector.
****************************************************************************************/
void GrahmSchmidt(Vector3& v0, Vector3& v1, Vector3& v2)
{
	v0.NormalizeThis();

	v1 -= (v1 * v0) * v0;
	v1.NormalizeThis();

	v2 = v0 % v1;	// no need to normalize
}


/****************************************************************************************
Function: JacobiSolver

Computes the eigenvalues and corresponding eigenvectors of a given 3x3 matrix using
Jacobi's method.

Re: Numerical Methods with Computer Programs in C++ by Pallab Ghosh, Section 5.7.1

m			: The 3x3 matrix to compute the eigensystem for.
eValues		: An array of floats to hold the resulting eigenvalues.
eVectors	: An array of Vector3s to hold the resulting eigenvectors.
****************************************************************************************/
void JacobiSolver(Matrix33 m, float eValues[3], Vector3 eVectors[3])
{
	const float eps1 = static_cast<float>(1.e-5);	// Error tolerances
	const float eps2 = static_cast<float>(1.e-5);
	const float eps3 = static_cast<float>(1.e-5);

	float p, q, spq;
	float cosa, sina;					// holds cos(alpha) and sin(alpha)
	float temp;							// used for temporary storage
	float s1 = 0.0f;					// sums of squares of diagonal
	float s2;							// elements

	bool flag = true;					// determines whether to iterate again.
	int iteration = 0;					// iteration counter

	Vector3 mik;						// used for temporary storage of m[i][k]

	Matrix33 t = Matrix33::Identity();	// stores the product of the rotation matrices.
										// Its columns ultimately hold the eigenvectors

	do
	{
		iteration++;

		for (int i = 0; i < 2; ++i)
		{
			for (int j = i + 1; j < 3; ++j)
			{
				if ((fabs(m.ColRow(j, i)) < eps1))
				{
					m.ColRow(j, i) = 0.0f;
				}
				else
				{
					q = fabs(m.ColRow(i, i) - m.ColRow(j, j));

					if (q > eps2)
					{
						p = 2.0f * m.ColRow(j, i) * q / (m.ColRow(i, i) - m.ColRow(j, j));
						spq = sqrt(p * p + q * q);
						cosa = sqrt((1.0f + q / spq) / 2.0f);
						sina = p / (2.0f * cosa * spq);
					}
					else
					{
						sina = cosa = INV_SQRT_TWO;
					}

					for (int k = 0; k < 3; ++k)
					{
						temp = t.ColRow(i, k);
						t.ColRow(i, k) = temp * cosa + t.ColRow(j, k) * sina;
						t.ColRow(j, k) = temp * sina - t.ColRow(j, k) * cosa;
					}

					for (int k = i; k < 3; ++k)
					{
						if (k > j)
						{
							temp = m.ColRow(k, i);
							m.ColRow(k, i) = cosa * temp + sina * m.ColRow(k, j);
							m.ColRow(k, j) = sina * temp - cosa * m.ColRow(k, j);
						}
						else
						{
							mik.data[k] = m.ColRow(k, i);
							m.ColRow(k, i) = cosa * mik.data[k] + sina * m.ColRow(j, k);

							if (k == j)
							{
								m.ColRow(k, j) = sina * mik.data[k] - cosa * m.ColRow(k, j);
							}
						}
					}

					mik.data[j] = sina * mik.data[i] - cosa * mik.data[j];

					for (int k = 0; k <= j; ++k)
					{
						if (k <= i)
						{
							temp = m.ColRow(i, k);
							m.ColRow(i, k) = cosa * temp + sina * m.ColRow(j, k);
							m.ColRow(j, k) = sina * temp - cosa * m.ColRow(j, k);
						}
						else
						{
							m.ColRow(j, k) = sina * mik.data[k] - cosa * m.ColRow(j, k);
						}
					}
				}
			}
		}

		s2 = 0.0f;

		for (int i = 0; i < 3; ++i)
		{
			eValues[i] = m.ColRow(i, i);
			s2 += eValues[i] * eValues[i];
		}

		if (fabs(s2) < static_cast<float>(1.e-5) || fabs(1 - s1 / s2) < eps3)
		{
			flag = false;
		}
		else
		{
			s1 = s2;
		}
	} while (flag);

	eVectors[0] = t.a1;
	eVectors[1] = t.a2;
	eVectors[2] = t.a3;

	// preserve righthanded-ness:
	if ((eVectors[0] % eVectors[1]) * eVectors[2] < 0.0f)
	{
		eVectors[2] = -eVectors[2];
	}

	std::cout << "\nEigenvectors converged in " << iteration << " iteration(s)" << std::endl;
}


/****************************************************************************************
Function: operator<<

Overloading the << operator to insert a Vector3 into a std::ostream for convenience.

os			: The std::ostream object to insert the Vector3 to.
v			: The Vector3 object to insert.
****************************************************************************************/
std::ostream& operator<<(std::ostream& os, const Vector3& v)
{
	os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
	return os;
}