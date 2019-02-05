#ifndef LIQUID_FUN_H_DEFINE
#define LIQUID_FUN_H_DEFINE

#include <cmath>
#include <math.h>
#include <limits.h>
#include <memory.h>
#include <string.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <stdint.h>
#include <new>
#include <algorithm>
#include <cfloat>
#include <stddef.h>
#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#define B2_NOT_USED(x) ((void)(x))
#if DEBUG && !defined(NDEBUG)
#define b2Assert(A) assert(A)
#define B2_ASSERT_ENABLED 1
#else
#define b2Assert(A)
#define B2_ASSERT_ENABLED 0
#endif

// Statement which is compiled out when DEBUG isn't defined.
#if DEBUG
#define B2_DEBUG_STATEMENT(A) A
#else
#define B2_DEBUG_STATEMENT(A)
#endif  // DEBUG

// Calculate the size of a static array.
#define B2_ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

typedef signed char	int8;
typedef signed short int16;
typedef signed int int32;
typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;
typedef float float32;
typedef double float64;

#ifdef WIN32
typedef __int64   int64;
typedef unsigned __int64   uint64;
#else // !WIN32
typedef long long int64;
typedef unsigned long long uint64;
#endif

#define	b2_maxFloat		FLT_MAX
#define	b2_epsilon		FLT_EPSILON
#define b2_pi			3.14159265359f

#if !defined(b2Inline)
#if defined(__GNUC__)
#define b2Inline __attribute__((always_inline))
#else
#define b2Inline inline
#endif // defined(__GNUC__)
#endif // !defined(b2Inline)

// We expand the API so that other languages (e.g. Java) can call into
// our C++ more easily. Only set if when the flag is not externally defined.
#if !defined(LIQUIDFUN_EXTERNAL_LANGUAGE_API)
#if defined(SWIG) || defined(LIQUIDFUN_UNIT_TESTS)
#define LIQUIDFUN_EXTERNAL_LANGUAGE_API 1
#else
#define LIQUIDFUN_EXTERNAL_LANGUAGE_API 0
#endif
#endif

/// @file
/// Global tuning constants based on meters-kilograms-seconds (MKS) units.
///

// Collision

/// The maximum number of contact points between two convex shapes. Do
/// not change this value.
#define b2_maxManifoldPoints	2

/// The maximum number of vertices on a convex polygon. You cannot increase
/// this too much because b2BlockAllocator has a maximum object size.
#define b2_maxPolygonVertices	8

/// This is used to fatten AABBs in the dynamic tree. This allows proxies
/// to move by a small amount without triggering a tree adjustment.
/// This is in meters.
#define b2_aabbExtension		0.1f

/// This is used to fatten AABBs in the dynamic tree. This is used to predict
/// the future position based on the current displacement.
/// This is a dimensionless multiplier.
#define b2_aabbMultiplier		2.0f

/// A small length used as a collision and constraint tolerance. Usually it is
/// chosen to be numerically significant, but visually insignificant.
#define b2_linearSlop			0.005f

/// A small angle used as a collision and constraint tolerance. Usually it is
/// chosen to be numerically significant, but visually insignificant.
#define b2_angularSlop			(2.0f / 180.0f * b2_pi)

/// The radius of the polygon/edge shape skin. This should not be modified. Making
/// this smaller means polygons will have an insufficient buffer for continuous collision.
/// Making it larger may create artifacts for vertex collision.
#define b2_polygonRadius		(2.0f * b2_linearSlop)

/// Maximum number of sub-steps per contact in continuous physics simulation.
#define b2_maxSubSteps			8


// Dynamics

/// Maximum number of contacts to be handled to solve a TOI impact.
#define b2_maxTOIContacts			32

/// A velocity threshold for elastic collisions. Any collision with a relative linear
/// velocity below this threshold will be treated as inelastic.
#define b2_velocityThreshold		1.0f

/// The maximum linear position correction used when solving constraints. This helps to
/// prevent overshoot.
#define b2_maxLinearCorrection		0.2f

/// The maximum angular position correction used when solving constraints. This helps to
/// prevent overshoot.
#define b2_maxAngularCorrection		(8.0f / 180.0f * b2_pi)

/// The maximum linear velocity of a body. This limit is very large and is used
/// to prevent numerical problems. You shouldn't need to adjust this.
#define b2_maxTranslation			2.0f
#define b2_maxTranslationSquared	(b2_maxTranslation * b2_maxTranslation)

/// The maximum angular velocity of a body. This limit is very large and is used
/// to prevent numerical problems. You shouldn't need to adjust this.
#define b2_maxRotation				(0.5f * b2_pi)
#define b2_maxRotationSquared		(b2_maxRotation * b2_maxRotation)

/// This scale factor controls how fast overlap is resolved. Ideally this would be 1 so
/// that overlap is removed in one time step. However using values close to 1 often lead
/// to overshoot.
#define b2_baumgarte				0.2f
#define b2_toiBaugarte				0.75f


// Particle

/// NEON SIMD requires 16-bit particle indices
#if !defined(B2_USE_16_BIT_PARTICLE_INDICES) && defined(LIQUIDFUN_SIMD_NEON)
#define B2_USE_16_BIT_PARTICLE_INDICES
#endif

/// A symbolic constant that stands for particle allocation error.
#define b2_invalidParticleIndex		(-1)

#ifdef B2_USE_16_BIT_PARTICLE_INDICES
#define b2_maxParticleIndex			0x7FFF
#else
#define b2_maxParticleIndex			0x7FFFFFFF
#endif

/// The default distance between particles, multiplied by the particle diameter.
#define b2_particleStride			0.75f

/// The minimum particle weight that produces pressure.
#define b2_minParticleWeight			1.0f

/// The upper limit for particle pressure.
#define b2_maxParticlePressure		0.25f

/// The upper limit for force between particles.
#define b2_maxParticleForce		0.5f

/// The maximum distance between particles in a triad, multiplied by the
/// particle diameter.
#define b2_maxTriadDistance			2
#define b2_maxTriadDistanceSquared		(b2_maxTriadDistance * b2_maxTriadDistance)

/// The initial size of particle data buffers.
#define b2_minParticleSystemBufferCapacity	256

/// The time into the future that collisions against barrier particles will be detected.
#define b2_barrierCollisionTime 2.5f

// Sleep

/// The time that a body must be still before it will go to sleep.
#define b2_timeToSleep				0.5f

/// A body cannot sleep if its linear velocity is above this tolerance.
#define b2_linearSleepTolerance		0.01f

/// A body cannot sleep if its angular velocity is above this tolerance.
#define b2_angularSleepTolerance	(2.0f / 180.0f * b2_pi)

// Memory Allocation

/// Implement this function to use your own memory allocator.
void* b2Alloc(int32 size);

/// If you implement b2Alloc, you should also implement this function.
void b2Free(void* mem);

/// Use this function to override b2Alloc() without recompiling this library.
typedef void* (*b2AllocFunction)(int32 size, void* callbackData);
/// Use this function to override b2Free() without recompiling this library.
typedef void (*b2FreeFunction)(void* mem, void* callbackData);

/// Set alloc and free callbacks to override the default behavior of using
/// malloc() and free() for dynamic memory allocation.
/// Set allocCallback and freeCallback to NULL to restore the default
/// allocator (malloc / free).
void b2SetAllocFreeCallbacks(b2AllocFunction allocCallback,
							 b2FreeFunction freeCallback,
							 void* callbackData);

/// Set the number of calls to b2Alloc minus the number of calls to b2Free.
/// This can be used to disable the empty heap check in
/// b2SetAllocFreeCallbacks() which can be useful for testing.
void b2SetNumAllocs(const int32 numAllocs);

/// Get number of calls to b2Alloc minus number of calls to b2Free.
int32 b2GetNumAllocs();

/// Logging function.
void b2Log(const char* string, ...);

/// Version numbering scheme.
/// See http://en.wikipedia.org/wiki/Software_versioning
struct b2Version
{
	int32 major;		///< significant changes
	int32 minor;		///< incremental changes
	int32 revision;		///< bug fixes
};

/// Current version.
/// Version of Box2D, LiquidFun is based upon.
extern b2Version b2_version;

/// Global variable is used to identify the version of LiquidFun.
extern const b2Version b2_liquidFunVersion;
/// String which identifies the current version of LiquidFun.
/// b2_liquidFunVersionString is used by Google developers to identify which
/// applications uploaded to Google Play are using this library.  This allows
/// the development team at Google to determine the popularity of the library.
/// How it works: Applications that are uploaded to the Google Play Store are
/// scanned for this version string.  We track which applications are using it
/// to measure popularity.  You are free to remove it (of course) but we would
/// appreciate if you left it in.
extern const char *b2_liquidFunVersionString;

// end of Settings.h

b2Version b2_version = {2, 3, 0};

#define LIQUIDFUN_VERSION_MAJOR 1
#define LIQUIDFUN_VERSION_MINOR 1
#define LIQUIDFUN_VERSION_REVISION 0
#define LIQUIDFUN_STRING_EXPAND(X) #X
#define LIQUIDFUN_STRING(X) LIQUIDFUN_STRING_EXPAND(X)

static void* b2AllocDefault(int32 size, void* callbackData);
static void b2FreeDefault(void* mem, void* callbackData);

const b2Version b2_liquidFunVersion = {
	LIQUIDFUN_VERSION_MAJOR, LIQUIDFUN_VERSION_MINOR,
	LIQUIDFUN_VERSION_REVISION,
};

const char *b2_liquidFunVersionString =
	"LiquidFun "
	LIQUIDFUN_STRING(LIQUIDFUN_VERSION_MAJOR) "."
	LIQUIDFUN_STRING(LIQUIDFUN_VERSION_MINOR) "."
	LIQUIDFUN_STRING(LIQUIDFUN_VERSION_REVISION);

static int32 b2_numAllocs = 0;

// Initialize default allocator.
static b2AllocFunction b2_allocCallback = b2AllocDefault;
static b2FreeFunction b2_freeCallback = b2FreeDefault;
static void *b2_callbackData = NULL;

// Default implementation of b2AllocFunction.
static void* b2AllocDefault(int32 size, void* callbackData)
{
	B2_NOT_USED(callbackData);
	return malloc(size);
}

// Default implementation of b2FreeFunction.
static void b2FreeDefault(void* mem, void* callbackData)
{
	B2_NOT_USED(callbackData);
	free(mem);
}

/// Set alloc and free callbacks to override the default behavior of using
/// malloc() and free() for dynamic memory allocation.
/// Set allocCallback and freeCallback to NULL to restore the default
/// allocator (malloc / free).
void b2SetAllocFreeCallbacks(b2AllocFunction allocCallback,
							 b2FreeFunction freeCallback, void* callbackData)
{
	b2Assert((allocCallback && freeCallback) ||
			 (!allocCallback && !freeCallback));
	b2Assert(0 == b2GetNumAllocs());
	if (allocCallback && freeCallback)
	{
		b2_allocCallback = allocCallback;
		b2_freeCallback = freeCallback;
		b2_callbackData = callbackData;
	}
	else
	{
		b2_allocCallback = b2AllocDefault;
		b2_freeCallback = b2FreeDefault;
		b2_callbackData = NULL;
	}
}

// Memory allocators. Modify these to use your own allocator.
void* b2Alloc(int32 size)
{
	b2_numAllocs++;
	return b2_allocCallback(size, b2_callbackData);
}

void b2Free(void* mem)
{
	b2_numAllocs--;
	b2_freeCallback(mem, b2_callbackData);
}

void b2SetNumAllocs(const int32 numAllocs)
{
	b2_numAllocs = numAllocs;
}

int32 b2GetNumAllocs()
{
	return b2_numAllocs;
}

// You can modify this to use your logging facility.
void b2Log(const char* string, ...)
{
#if DEBUG
	va_list args;
	va_start(args, string);
	vprintf(string, args);
	va_end(args);
#else
	B2_NOT_USED(string);
#endif
}

class Validator
{
public:
	Validator()
	{
		b2Assert(sizeof(uint64)==8);
		b2Assert(sizeof(int64)==8);
	}
} validate;

// end of Settings.cpp

/// This function is used to ensure that a floating point number is not a NaN or infinity.
inline bool b2IsValid(float32 x)
{
	union {
		float32 f;
		int32 i;
	} v = { x };
	return (v.i & 0x7f800000) != 0x7f800000;
}

/// This is a approximate yet fast inverse square-root.
inline float32 b2InvSqrt(float32 x)
{
	union
	{
		float32 x;
		int32 i;
	} convert;

	convert.x = x;
	float32 xhalf = 0.5f * x;
	convert.i = 0x5f3759df - (convert.i >> 1);
	x = convert.x;
	x = x * (1.5f - xhalf * x * x);
	return x;
}

#define	b2Sqrt(x)	sqrtf(x)
#define	b2Atan2(y, x)	atan2f(y, x)

/// A 2D column vector.
struct b2Vec2
{
	/// Default constructor does nothing (for performance).
	b2Vec2() {}

	/// Construct using coordinates.
	b2Vec2(float32 x, float32 y) : x(x), y(y) {}

	/// Set this vector to all zeros.
	void SetZero() { x = 0.0f; y = 0.0f; }

	/// Set this vector to some specified coordinates.
	void Set(float32 x_, float32 y_) { x = x_; y = y_; }

	/// Negate this vector.
	b2Vec2 operator -() const { b2Vec2 v; v.Set(-x, -y); return v; }

	/// Read from and indexed element.
	float32 operator () (int32 i) const
	{
		return (&x)[i];
	}

	/// Write to an indexed element.
	float32& operator () (int32 i)
	{
		return (&x)[i];
	}

	/// Add a vector to this vector.
	void operator += (const b2Vec2& v)
	{
		x += v.x; y += v.y;
	}

	/// Subtract a vector from this vector.
	void operator -= (const b2Vec2& v)
	{
		x -= v.x; y -= v.y;
	}

	/// Multiply this vector by a scalar.
	void operator *= (float32 a)
	{
		x *= a; y *= a;
	}

	/// Get the length of this vector (the norm).
	float32 Length() const
	{
		return b2Sqrt(x * x + y * y);
	}

	/// Get the length squared. For performance, use this instead of
	/// b2Vec2::Length (if possible).
	float32 LengthSquared() const
	{
		return x * x + y * y;
	}

	/// Convert this vector into a unit vector. Returns the length.
	float32 Normalize()
	{
		float32 length = Length();
		if (length < b2_epsilon)
		{
			return 0.0f;
		}
		float32 invLength = 1.0f / length;
		x *= invLength;
		y *= invLength;

		return length;
	}

	/// Does this vector contain finite coordinates?
	bool IsValid() const
	{
		return b2IsValid(x) && b2IsValid(y);
	}

	/// Get the skew vector such that dot(skew_vec, other) == cross(vec, other)
	b2Vec2 Skew() const
	{
		return b2Vec2(-y, x);
	}

	float32 x, y;
};

/// Add a float to a vector.
inline b2Vec2 operator + (const b2Vec2& v, float f)
{
	return b2Vec2(v.x + f, v.y + f);
}

/// Substract a float from a vector.
inline b2Vec2 operator - (const b2Vec2& v, float f)
{
	return b2Vec2(v.x - f, v.y - f);
}

/// Multiply a float with a vector.
inline b2Vec2 operator * (const b2Vec2& v, float f)
{
	return b2Vec2(v.x * f, v.y * f);
}

/// Divide a vector by a float.
inline b2Vec2 operator / (const b2Vec2& v, float f)
{
	return b2Vec2(v.x / f, v.y / f);
}

/// A 3D column vector with 3 elements.
struct b2Vec3
{
	/// Default constructor does nothing (for performance).
	b2Vec3() {}

	/// Construct using coordinates.
	b2Vec3(float32 x, float32 y, float32 z) : x(x), y(y), z(z) {}

	/// Set this vector to all zeros.
	void SetZero() { x = 0.0f; y = 0.0f; z = 0.0f; }

	/// Set this vector to some specified coordinates.
	void Set(float32 x_, float32 y_, float32 z_) { x = x_; y = y_; z = z_; }

	/// Negate this vector.
	b2Vec3 operator -() const { b2Vec3 v; v.Set(-x, -y, -z); return v; }

	/// Add a vector to this vector.
	void operator += (const b2Vec3& v)
	{
		x += v.x; y += v.y; z += v.z;
	}

	/// Subtract a vector from this vector.
	void operator -= (const b2Vec3& v)
	{
		x -= v.x; y -= v.y; z -= v.z;
	}

	/// Multiply this vector by a scalar.
	void operator *= (float32 s)
	{
		x *= s; y *= s; z *= s;
	}

		/// Get the length of this vector (the norm).
	float32 Length() const
	{
		return b2Sqrt(x * x + y * y + z * z);
	}

	/// Convert this vector into a unit vector. Returns the length.
	float32 Normalize()
	{
		float32 length = Length();
		if (length < b2_epsilon)
		{
			return 0.0f;
		}
		float32 invLength = 1.0f / length;
		x *= invLength;
		y *= invLength;
		z *= invLength;

		return length;
	}

	float32 x, y, z;
};

/// A 4D column vector with 4 elements.
struct b2Vec4
{
	/// Default constructor does nothing (for performance).
	b2Vec4() {}

	/// Construct using coordinates.
	b2Vec4(float32 x, float32 y, float32 z, float32 w) : x(x), y(y), z(z), w(w) {}

	float32 x, y, z, w;
};

/// A 2-by-2 matrix. Stored in column-major order.
struct b2Mat22
{
	/// The default constructor does nothing (for performance).
	b2Mat22() {}

	/// Construct this matrix using columns.
	b2Mat22(const b2Vec2& c1, const b2Vec2& c2)
	{
		ex = c1;
		ey = c2;
	}

	/// Construct this matrix using scalars.
	b2Mat22(float32 a11, float32 a12, float32 a21, float32 a22)
	{
		ex.x = a11; ex.y = a21;
		ey.x = a12; ey.y = a22;
	}

	/// Initialize this matrix using columns.
	void Set(const b2Vec2& c1, const b2Vec2& c2)
	{
		ex = c1;
		ey = c2;
	}

	/// Set this to the identity matrix.
	void SetIdentity()
	{
		ex.x = 1.0f; ey.x = 0.0f;
		ex.y = 0.0f; ey.y = 1.0f;
	}

	/// Set this matrix to all zeros.
	void SetZero()
	{
		ex.x = 0.0f; ey.x = 0.0f;
		ex.y = 0.0f; ey.y = 0.0f;
	}

	b2Mat22 GetInverse() const
	{
		float32 a = ex.x, b = ey.x, c = ex.y, d = ey.y;
		b2Mat22 B;
		float32 det = a * d - b * c;
		if (det != 0.0f)
		{
			det = 1.0f / det;
		}
		B.ex.x =  det * d;	B.ey.x = -det * b;
		B.ex.y = -det * c;	B.ey.y =  det * a;
		return B;
	}

	/// Solve A * x = b, where b is a column vector. This is more efficient
	/// than computing the inverse in one-shot cases.
	b2Vec2 Solve(const b2Vec2& b) const
	{
		float32 a11 = ex.x, a12 = ey.x, a21 = ex.y, a22 = ey.y;
		float32 det = a11 * a22 - a12 * a21;
		if (det != 0.0f)
		{
			det = 1.0f / det;
		}
		b2Vec2 x;
		x.x = det * (a22 * b.x - a12 * b.y);
		x.y = det * (a11 * b.y - a21 * b.x);
		return x;
	}

	b2Vec2 ex, ey;
};

/// A 3-by-3 matrix. Stored in column-major order.
struct b2Mat33
{
	/// The default constructor does nothing (for performance).
	b2Mat33() {}

	/// Construct this matrix using columns.
	b2Mat33(const b2Vec3& c1, const b2Vec3& c2, const b2Vec3& c3)
	{
		ex = c1;
		ey = c2;
		ez = c3;
	}

	/// Set this matrix to all zeros.
	void SetZero()
	{
		ex.SetZero();
		ey.SetZero();
		ez.SetZero();
	}

	/// Solve A * x = b, where b is a column vector. This is more efficient
	/// than computing the inverse in one-shot cases.
	b2Vec3 Solve33(const b2Vec3& b) const;

	/// Solve A * x = b, where b is a column vector. This is more efficient
	/// than computing the inverse in one-shot cases. Solve only the upper
	/// 2-by-2 matrix equation.
	b2Vec2 Solve22(const b2Vec2& b) const;

	/// Get the inverse of this matrix as a 2-by-2.
	/// Returns the zero matrix if singular.
	void GetInverse22(b2Mat33* M) const;

	/// Get the symmetric inverse of this matrix as a 3-by-3.
	/// Returns the zero matrix if singular.
	void GetSymInverse33(b2Mat33* M) const;

	b2Vec3 ex, ey, ez;
};

/// Rotation
struct b2Rot
{
	b2Rot() {}

	/// Initialize from an angle in radians
	explicit b2Rot(float32 angle)
	{
		/// TODO_ERIN optimize
		s = sinf(angle);
		c = cosf(angle);
	}

	/// Set using an angle in radians.
	void Set(float32 angle)
	{
		/// TODO_ERIN optimize
		s = sinf(angle);
		c = cosf(angle);
	}

	/// Set to the identity rotation
	void SetIdentity()
	{
		s = 0.0f;
		c = 1.0f;
	}

	/// Get the angle in radians
	float32 GetAngle() const
	{
		return b2Atan2(s, c);
	}

	/// Get the x-axis
	b2Vec2 GetXAxis() const
	{
		return b2Vec2(c, s);
	}

	/// Get the u-axis
	b2Vec2 GetYAxis() const
	{
		return b2Vec2(-s, c);
	}

	/// Sine and cosine
	float32 s, c;
};

/// A transform contains translation and rotation. It is used to represent
/// the position and orientation of rigid frames.
struct b2Transform
{
	/// The default constructor does nothing.
	b2Transform() {}

	/// Initialize using a position vector and a rotation.
	b2Transform(const b2Vec2& position, const b2Rot& rotation) : p(position), q(rotation) {}

	/// Set this to the identity transform.
	void SetIdentity()
	{
		p.SetZero();
		q.SetIdentity();
	}

	/// Set this based on the position and angle.
	void Set(const b2Vec2& position, float32 angle)
	{
		p = position;
		q.Set(angle);
	}

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
	/// Get x-coordinate of p.
	float32 GetPositionX() const { return p.x; }

	/// Get y-coordinate of p.
	float32 GetPositionY() const { return p.y; }

	/// Get sine-component of q.
	float32 GetRotationSin() const { return q.s; }

	/// Get cosine-component of q.
	float32 GetRotationCos() const { return q.c; }
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

	b2Vec2 p;
	b2Rot q;
};

/// This describes the motion of a body/shape for TOI computation.
/// Shapes are defined with respect to the body origin, which may
/// no coincide with the center of mass. However, to support dynamics
/// we must interpolate the center of mass position.
struct b2Sweep
{
	/// Get the interpolated transform at a specific time.
	/// @param beta is a factor in [0,1], where 0 indicates alpha0.
	void GetTransform(b2Transform* xfb, float32 beta) const;

	/// Advance the sweep forward, yielding a new initial state.
	/// @param alpha the new initial time.
	void Advance(float32 alpha);

	/// Normalize the angles.
	void Normalize();

	b2Vec2 localCenter;	///< local center of mass position
	b2Vec2 c0, c;		///< center world positions
	float32 a0, a;		///< world angles

	/// Fraction of the current time step in the range [0,1]
	/// c0 and a0 are the positions at alpha0.
	float32 alpha0;
};

/// Useful constant
extern const b2Vec2 b2Vec2_zero;

/// Perform the dot product on two vectors.
inline float32 b2Dot(const b2Vec2& a, const b2Vec2& b)
{
	return a.x * b.x + a.y * b.y;
}

/// Perform the cross product on two vectors. In 2D this produces a scalar.
inline float32 b2Cross(const b2Vec2& a, const b2Vec2& b)
{
	return a.x * b.y - a.y * b.x;
}

/// Perform the cross product on a vector and a scalar. In 2D this produces
/// a vector.
inline b2Vec2 b2Cross(const b2Vec2& a, float32 s)
{
	return b2Vec2(s * a.y, -s * a.x);
}

/// Perform the cross product on a scalar and a vector. In 2D this produces
/// a vector.
inline b2Vec2 b2Cross(float32 s, const b2Vec2& a)
{
	return b2Vec2(-s * a.y, s * a.x);
}

/// Multiply a matrix times a vector. If a rotation matrix is provided,
/// then this transforms the vector from one frame to another.
inline b2Vec2 b2Mul(const b2Mat22& A, const b2Vec2& v)
{
	return b2Vec2(A.ex.x * v.x + A.ey.x * v.y, A.ex.y * v.x + A.ey.y * v.y);
}

/// Multiply a matrix transpose times a vector. If a rotation matrix is provided,
/// then this transforms the vector from one frame to another (inverse transform).
inline b2Vec2 b2MulT(const b2Mat22& A, const b2Vec2& v)
{
	return b2Vec2(b2Dot(v, A.ex), b2Dot(v, A.ey));
}

/// Add two vectors component-wise.
inline b2Vec2 operator + (const b2Vec2& a, const b2Vec2& b)
{
	return b2Vec2(a.x + b.x, a.y + b.y);
}

/// Subtract two vectors component-wise.
inline b2Vec2 operator - (const b2Vec2& a, const b2Vec2& b)
{
	return b2Vec2(a.x - b.x, a.y - b.y);
}

inline b2Vec2 operator * (float32 s, const b2Vec2& a)
{
	return b2Vec2(s * a.x, s * a.y);
}

inline bool operator == (const b2Vec2& a, const b2Vec2& b)
{
	return a.x == b.x && a.y == b.y;
}

inline bool operator != (const b2Vec2& a, const b2Vec2& b)
{
	return !operator==(a, b);
}

inline float32 b2Distance(const b2Vec2& a, const b2Vec2& b)
{
	b2Vec2 c = a - b;
	return c.Length();
}

inline float32 b2DistanceSquared(const b2Vec2& a, const b2Vec2& b)
{
	b2Vec2 c = a - b;
	return b2Dot(c, c);
}

inline b2Vec3 operator * (float32 s, const b2Vec3& a)
{
	return b2Vec3(s * a.x, s * a.y, s * a.z);
}

/// Add two vectors component-wise.
inline b2Vec3 operator + (const b2Vec3& a, const b2Vec3& b)
{
	return b2Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

/// Subtract two vectors component-wise.
inline b2Vec3 operator - (const b2Vec3& a, const b2Vec3& b)
{
	return b2Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

/// Perform the dot product on two vectors.
inline float32 b2Dot(const b2Vec3& a, const b2Vec3& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

/// Perform the cross product on two vectors.
inline b2Vec3 b2Cross(const b2Vec3& a, const b2Vec3& b)
{
	return b2Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

inline b2Mat22 operator + (const b2Mat22& A, const b2Mat22& B)
{
	return b2Mat22(A.ex + B.ex, A.ey + B.ey);
}

// A * B
inline b2Mat22 b2Mul(const b2Mat22& A, const b2Mat22& B)
{
	return b2Mat22(b2Mul(A, B.ex), b2Mul(A, B.ey));
}

// A^T * B
inline b2Mat22 b2MulT(const b2Mat22& A, const b2Mat22& B)
{
	b2Vec2 c1(b2Dot(A.ex, B.ex), b2Dot(A.ey, B.ex));
	b2Vec2 c2(b2Dot(A.ex, B.ey), b2Dot(A.ey, B.ey));
	return b2Mat22(c1, c2);
}

/// Multiply a matrix times a vector.
inline b2Vec3 b2Mul(const b2Mat33& A, const b2Vec3& v)
{
	return v.x * A.ex + v.y * A.ey + v.z * A.ez;
}

/// Multiply a matrix times a vector.
inline b2Vec2 b2Mul22(const b2Mat33& A, const b2Vec2& v)
{
	return b2Vec2(A.ex.x * v.x + A.ey.x * v.y, A.ex.y * v.x + A.ey.y * v.y);
}

/// Multiply two rotations: q * r
inline b2Rot b2Mul(const b2Rot& q, const b2Rot& r)
{
	// [qc -qs] * [rc -rs] = [qc*rc-qs*rs -qc*rs-qs*rc]
	// [qs  qc]   [rs  rc]   [qs*rc+qc*rs -qs*rs+qc*rc]
	// s = qs * rc + qc * rs
	// c = qc * rc - qs * rs
	b2Rot qr;
	qr.s = q.s * r.c + q.c * r.s;
	qr.c = q.c * r.c - q.s * r.s;
	return qr;
}

/// Transpose multiply two rotations: qT * r
inline b2Rot b2MulT(const b2Rot& q, const b2Rot& r)
{
	// [ qc qs] * [rc -rs] = [qc*rc+qs*rs -qc*rs+qs*rc]
	// [-qs qc]   [rs  rc]   [-qs*rc+qc*rs qs*rs+qc*rc]
	// s = qc * rs - qs * rc
	// c = qc * rc + qs * rs
	b2Rot qr;
	qr.s = q.c * r.s - q.s * r.c;
	qr.c = q.c * r.c + q.s * r.s;
	return qr;
}

/// Rotate a vector
inline b2Vec2 b2Mul(const b2Rot& q, const b2Vec2& v)
{
	return b2Vec2(q.c * v.x - q.s * v.y, q.s * v.x + q.c * v.y);
}

/// Inverse rotate a vector
inline b2Vec2 b2MulT(const b2Rot& q, const b2Vec2& v)
{
	return b2Vec2(q.c * v.x + q.s * v.y, -q.s * v.x + q.c * v.y);
}

inline b2Vec2 b2Mul(const b2Transform& T, const b2Vec2& v)
{
	float32 x = (T.q.c * v.x - T.q.s * v.y) + T.p.x;
	float32 y = (T.q.s * v.x + T.q.c * v.y) + T.p.y;

	return b2Vec2(x, y);
}

inline b2Vec2 b2MulT(const b2Transform& T, const b2Vec2& v)
{
	float32 px = v.x - T.p.x;
	float32 py = v.y - T.p.y;
	float32 x = (T.q.c * px + T.q.s * py);
	float32 y = (-T.q.s * px + T.q.c * py);

	return b2Vec2(x, y);
}

// v2 = A.q.Rot(B.q.Rot(v1) + B.p) + A.p
//    = (A.q * B.q).Rot(v1) + A.q.Rot(B.p) + A.p
inline b2Transform b2Mul(const b2Transform& A, const b2Transform& B)
{
	b2Transform C;
	C.q = b2Mul(A.q, B.q);
	C.p = b2Mul(A.q, B.p) + A.p;
	return C;
}

// v2 = A.q' * (B.q * v1 + B.p - A.p)
//    = A.q' * B.q * v1 + A.q' * (B.p - A.p)
inline b2Transform b2MulT(const b2Transform& A, const b2Transform& B)
{
	b2Transform C;
	C.q = b2MulT(A.q, B.q);
	C.p = b2MulT(A.q, B.p - A.p);
	return C;
}

template <typename T>
inline T b2Abs(T a)
{
	return a > T(0) ? a : -a;
}

inline b2Vec2 b2Abs(const b2Vec2& a)
{
	return b2Vec2(b2Abs(a.x), b2Abs(a.y));
}

inline b2Mat22 b2Abs(const b2Mat22& A)
{
	return b2Mat22(b2Abs(A.ex), b2Abs(A.ey));
}

template <typename T>
inline T b2Min(T a, T b)
{
	return a < b ? a : b;
}

inline b2Vec2 b2Min(const b2Vec2& a, const b2Vec2& b)
{
	return b2Vec2(b2Min(a.x, b.x), b2Min(a.y, b.y));
}

template <typename T>
inline T b2Max(T a, T b)
{
	return a > b ? a : b;
}

inline b2Vec2 b2Max(const b2Vec2& a, const b2Vec2& b)
{
	return b2Vec2(b2Max(a.x, b.x), b2Max(a.y, b.y));
}

template <typename T>
inline T b2Clamp(T a, T low, T high)
{
	return b2Max(low, b2Min(a, high));
}

inline b2Vec2 b2Clamp(const b2Vec2& a, const b2Vec2& low, const b2Vec2& high)
{
	return b2Max(low, b2Min(a, high));
}

template<typename T> inline void b2Swap(T& a, T& b)
{
	T tmp = a;
	a = b;
	b = tmp;
}

/// "Next Largest Power of 2
/// Given a binary integer value x, the next largest power of 2 can be computed by a SWAR algorithm
/// that recursively "folds" the upper bits into the lower bits. This process yields a bit vector with
/// the same most significant 1 as x, but all 1's below it. Adding 1 to that value yields the next
/// largest power of 2. For a 32-bit value:"
inline uint32 b2NextPowerOfTwo(uint32 x)
{
	x |= (x >> 1);
	x |= (x >> 2);
	x |= (x >> 4);
	x |= (x >> 8);
	x |= (x >> 16);
	return x + 1;
}

inline bool b2IsPowerOfTwo(uint32 x)
{
	bool result = x > 0 && (x & (x - 1)) == 0;
	return result;
}

inline void b2Sweep::GetTransform(b2Transform* xf, float32 beta) const
{
	xf->p = (1.0f - beta) * c0 + beta * c;
	float32 angle = (1.0f - beta) * a0 + beta * a;
	xf->q.Set(angle);

	// Shift to origin
	xf->p -= b2Mul(xf->q, localCenter);
}

inline void b2Sweep::Advance(float32 alpha)
{
	b2Assert(alpha0 < 1.0f);
	float32 beta = (alpha - alpha0) / (1.0f - alpha0);
	c0 += beta * (c - c0);
	a0 += beta * (a - a0);
	alpha0 = alpha;
}

/// Normalize an angle in radians to be between -pi and pi
inline void b2Sweep::Normalize()
{
	float32 twoPi = 2.0f * b2_pi;
	float32 d =  twoPi * floorf(a0 / twoPi);
	a0 -= d;
	a -= d;
}

// end of Math.h

/// Calculates min/max/mean of a set of samples
class b2Stat
{
public:
	b2Stat();

	/// Record a sample
	void Record( float32 t );

	/// Returns the number of recorded samples
	int GetCount() const;

	/// Returns the mean of all recorded samples,
	/// Returns 0 if there are no recorded samples
	float32 GetMean() const;

	/// Returns the min of all recorded samples,
	/// FLT_MAX if there are no recorded samples
	float32 GetMin() const;

	/// Returns the max of all recorded samples,
	/// -FLT_MAX if there are no recorded samples
	float32 GetMax() const;

	/// Erase all recorded samples
	void Clear();
private:

	int m_count;
	float64 m_total;
	float32 m_min;
	float32 m_max;
};

// end of Stat.h

b2Stat::b2Stat()
{
	Clear();
}

void b2Stat::Record( float32 t )
{
	m_total += t;
	m_min = std::min(m_min,t);
	m_max = std::max(m_max,t);
	m_count++;
}

int b2Stat::GetCount() const
{
	return m_count;
}

float32 b2Stat::GetMean() const
{
	if (m_count == 0)
	{
		return 0.0f;
	}
	return (float32)(m_total / m_count);
}

float32 b2Stat::GetMin() const
{
	return m_min;
}

float32 b2Stat::GetMax() const
{
	return m_max;
}

void b2Stat::Clear()
{
	m_count = 0;
	m_total = 0;
	m_min = FLT_MAX;
	m_max = -FLT_MAX;
}

// end of Stat.cpp

/// Timer for profiling. This has platform specific code and may
/// not work on every platform.
class b2Timer
{
public:

	/// Constructor
	b2Timer();

	/// Reset the timer.
	void Reset();

	/// Get the time since construction or the last reset.
	float32 GetMilliseconds() const;

private:
	/// Get platform specific tick count
	static int64 GetTicks();

#if defined(_WIN32)
	static float64 s_invFrequency;
#endif
	int64 m_start;
};


// end of Timer.h

#if defined(_WIN32)

float64 b2Timer::s_invFrequency = 0.0f;

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

typedef BOOL (WINAPI *SystemGetTimeFunc)(_Out_ LARGE_INTEGER *lpFrequency);
SystemGetTimeFunc systemGetTimeFunc = ::QueryPerformanceCounter;
SystemGetTimeFunc systemGetFreqFunc = ::QueryPerformanceFrequency;

int64 b2Timer::GetTicks()
{
	LARGE_INTEGER largeInteger;
	systemGetTimeFunc(&largeInteger);
	return largeInteger.QuadPart;
}

b2Timer::b2Timer()
{
	LARGE_INTEGER largeInteger;

	if (s_invFrequency == 0.0f)
	{
		systemGetFreqFunc(&largeInteger);
		s_invFrequency = float64(largeInteger.QuadPart);
		if (s_invFrequency > 0.0f)
		{
			s_invFrequency = 1000.0f / s_invFrequency;
		}
	}

	m_start = GetTicks();
}

void b2Timer::Reset()
{
	m_start = GetTicks();
}

float32 b2Timer::GetMilliseconds() const
{
	int64 elapsed = GetTicks() - m_start;
	return (float32)(s_invFrequency * elapsed);
}

#elif defined(__linux__) || defined (__APPLE__)

#include <sys/time.h>
#include <time.h>

// systemGetTimeFunc is defined with external linkage to allow unit
// test to mock out the system time function

#if defined(__linux__)

typedef int (*SystemGetTimeFunc)(clockid_t clk_id, struct timespec *tp);
SystemGetTimeFunc systemGetTimeFunc = ::clock_gettime;

#elif defined(__APPLE__)

typedef int (*SystemGetTimeFunc)(struct timeval * tp, void * tzp);
SystemGetTimeFunc systemGetTimeFunc = ::gettimeofday;

#endif

int64 b2Timer::GetTicks()
{
	static const int NSEC_PER_SEC = 1000000000;

#ifdef __linux__
	timespec ts;
	systemGetTimeFunc(CLOCK_MONOTONIC,&ts);
	return ((int64)ts.tv_sec) * NSEC_PER_SEC + ts.tv_nsec;
#else
	timeval t;
	systemGetTimeFunc(&t, 0);
	return ((int64)t.tv_sec) * NSEC_PER_SEC + t.tv_usec * 1000;
#endif
}

b2Timer::b2Timer()
{
	Reset();
}

void b2Timer::Reset()
{
	m_start = GetTicks();
}

float32 b2Timer::GetMilliseconds() const
{
	static const float32 kTicksToMs = 0.000001f;
	return kTicksToMs * (float32)(GetTicks() - m_start);
}

#else

b2Timer::b2Timer()
{
}

void b2Timer::Reset()
{
}

float32 b2Timer::GetMilliseconds() const
{
	return 0.0f;
}

#endif

// end of Timer.cpp

// Whether to enable b2IntrusiveList::ValidateList().
// Be careful when enabling this since this changes the size of
// b2IntrusiveListNode so make sure *all* projects that include Box2D.h
// also define this value in the same way to avoid data corruption.
#ifndef B2_INTRUSIVE_LIST_VALIDATE
#define B2_INTRUSIVE_LIST_VALIDATE 0
#endif  // B2_INTRUSIVE_LIST_VALIDATE

/// b2IntrusiveListNode is used to implement an intrusive doubly-linked
/// list.
///
/// For example:
///
/// class MyClass {
/// public:
/// 	MyClass(const char *msg) : m_msg(msg) {}
/// 	const char* GetMessage() const { return m_msg; }
/// 	B2_INTRUSIVE_LIST_GET_NODE(m_node);
/// 	B2_INTRUSIVE_LIST_NODE_GET_CLASS(MyClass, m_node);
/// private:
/// 	b2IntrusiveListNode m_node;
/// 	const char *m_msg;
/// };
///
/// int main(int argc, char *argv[]) {
/// 	b2IntrusiveListNode list; // NOTE: type is NOT MyClass
/// 	MyClass a("this");
/// 	MyClass b("is");
/// 	MyClass c("a");
/// 	MyClass d("test");
/// 	list.InsertBefore(a.GetListNode());
/// 	list.InsertBefore(b.GetListNode());
/// 	list.InsertBefore(c.GetListNode());
/// 	list.InsertBefore(d.GetListNode());
/// 	for (b2IntrusiveListNode* node = list.GetNext();
/// 		 node != list.GetTerminator(); node = node->GetNext()) {
/// 		MyClass *cls = MyClass::GetInstanceFromListNode(node);
/// 		printf("%s\n", cls->GetMessage());
/// 	}
/// 	return 0;
/// }
class b2IntrusiveListNode
{
public:
	/// Initialize the node.
	b2IntrusiveListNode()
	{
		Initialize();
#if B2_INTRUSIVE_LIST_VALIDATE
		m_magic = k_magic;
#endif // B2_INTRUSIVE_LIST_VALIDATE
	}

	/// If the node is in a list, remove it from the list.
	~b2IntrusiveListNode()
	{
		Remove();
#if B2_INTRUSIVE_LIST_VALIDATE
		m_magic = 0;
#endif // B2_INTRUSIVE_LIST_VALIDATE
	}

	/// Insert this node after the specified node.
	void InsertAfter(b2IntrusiveListNode* const node)
	{
		b2Assert(!node->InList());
		node->m_next = m_next;
		node->m_prev = this;
		m_next->m_prev = node;
		m_next = node;
	}

	/// Insert this node before the specified node.
	void InsertBefore(b2IntrusiveListNode* const node)
	{
		b2Assert(!node->InList());
		node->m_next = this;
		node->m_prev = m_prev;
		m_prev->m_next = node;
		m_prev = node;
	}

	/// Get the terminator of the list.
	const b2IntrusiveListNode* GetTerminator() const
	{
		return this;
	}

	/// Remove this node from the list it's currently in.
	b2IntrusiveListNode* Remove()
	{
		m_prev->m_next = m_next;
		m_next->m_prev = m_prev;
		Initialize();
		return this;
	}

	/// Determine whether this list is empty or the node isn't in a list.
	bool IsEmpty() const
	{
	  return GetNext() == this;
	}

	/// Determine whether this node is in a list or the list contains nodes.
	bool InList() const
	{
	  return !IsEmpty();
	}

	/// Calculate the length of the list.
	uint32 GetLength() const
	{
		uint32 length = 0;
		const b2IntrusiveListNode * const terminator = GetTerminator();
		for (const b2IntrusiveListNode* node = GetNext();
			 node != terminator; node = node->GetNext())
		{
			length++;
		}
		return length;
	}

	/// Get the next node in the list.
	b2IntrusiveListNode* GetNext() const
	{
		return m_next;
	}

	/// Get the previous node in the list.
	b2IntrusiveListNode* GetPrevious() const
	{
		return m_prev;
	}

	/// If B2_INTRUSIVE_LIST_VALIDATE is 1 perform a very rough validation
	/// of all nodes in the list.
	bool ValidateList() const
	{
#if B2_INTRUSIVE_LIST_VALIDATE
	  if (m_magic != k_magic) return false;
	  const b2IntrusiveListNode * const terminator = GetTerminator();
	  for (b2IntrusiveListNode *node = GetNext(); node != terminator;
		   node = node->GetNext()) {
		if (node->m_magic != k_magic) return false;
	  }
#endif  // B2_INTRUSIVE_LIST_VALIDATE
	  return true;
	}

	/// Determine whether the specified node is present in this list.
	bool FindNodeInList(b2IntrusiveListNode* const nodeToFind) const
	{
		const b2IntrusiveListNode * const terminator = GetTerminator();
		for (b2IntrusiveListNode *node = GetNext(); node != terminator;
			 node = node->GetNext())
		{
			if (nodeToFind == node) return true;
		}
		return false;
	}

private:
	/// Initialize the list node.
	void Initialize()
	{
		m_next = this;
		m_prev = this;
	}

private:
#if B2_INTRUSIVE_LIST_VALIDATE
	uint32 m_magic;
#endif  // B2_INTRUSIVE_LIST_VALIDATE
	/// The next node in the list.
	b2IntrusiveListNode *m_prev;
	/// The previous node in the list.
	b2IntrusiveListNode *m_next;

private:
#if B2_INTRUSIVE_LIST_VALIDATE
	static const uint32 k_magic = 0x7157ac01;
#endif  // B2_INTRUSIVE_LIST_VALIDATE
};

/// Declares the member function GetListNode() of Class to retrieve a pointer
/// to NodeMemberName.
/// See #B2_INTRUSIVE_LIST_NODE_GET_CLASS_ACCESSOR()
#define B2_INTRUSIVE_LIST_GET_NODE(NodeMemberName) \
	b2IntrusiveListNode* GetListNode() { return &NodeMemberName; } \
	const b2IntrusiveListNode* GetListNode() const { return &NodeMemberName; }

/// Declares the member function FunctionName of Class to retrieve a pointer
/// to a Class instance from a list node pointer.   NodeMemberName references
/// the name of the b2IntrusiveListNode member of Class.
#define B2_INTRUSIVE_LIST_NODE_GET_CLASS_ACCESSOR( \
	Class, NodeMemberName, FunctionName) \
	static Class* FunctionName(b2IntrusiveListNode *node) \
	{ \
		Class *cls = NULL; \
		/* This effectively performs offsetof(Class, NodeMemberName) */ \
		/* which ends up in the undefined behavior realm of C++ but in */ \
		/* practice this works with most compilers. */ \
		return reinterpret_cast<Class*>((uint8*)(node) - \
										(uint8*)(&cls->NodeMemberName)); \
	} \
	\
	static const Class* FunctionName(const b2IntrusiveListNode *node) \
	{ \
		return FunctionName(const_cast<b2IntrusiveListNode*>(node)); \
	}

/// Declares the member function GetInstanceFromListNode() of Class to retrieve
/// a pointer to a Class instance from a list node pointer.  NodeMemberName
/// reference the name of the b2IntrusiveListNode member of Class.
#define B2_INTRUSIVE_LIST_NODE_GET_CLASS(Class, NodeMemberName) \
	B2_INTRUSIVE_LIST_NODE_GET_CLASS_ACCESSOR(Class, NodeMemberName, \
											  GetInstanceFromListNode)

/// b2TypedIntrusiveListNode which supports inserting an object into a single
/// doubly linked list.  For objects that need to be inserted in multiple
/// doubly linked lists, use b2IntrusiveListNode.
///
/// For example:
///
/// class IntegerItem : public b2TypedIntrusiveListNode<IntegerItem>
/// {
/// public:
/// 	IntegerItem(int32 value) : m_value(value) { }
/// 	~IntegerItem() { }
/// 	int32 GetValue() const { return m_value; }
/// private:
/// 	int32 m_value;
/// };
///
/// int main(int argc, const char *arvg[]) {
/// 	b2TypedIntrusiveListNode<IntegerItem> list;
/// 	IntegerItem a(1);
/// 	IntegerItem b(2);
/// 	IntegerItem c(3);
/// 	list.InsertBefore(&a);
/// 	list.InsertBefore(&b);
/// 	list.InsertBefore(&c);
/// 	for (IntegerItem* item = list.GetNext();
/// 		 item != list.GetTerminator(); item = item->GetNext())
/// 	{
/// 		printf("%d\n", item->GetValue());
/// 	}
/// }
template<typename T>
class b2TypedIntrusiveListNode
{
public:
	b2TypedIntrusiveListNode() { }
	~b2TypedIntrusiveListNode() { }

	/// Insert this object after the specified object.
	void InsertAfter(T* const obj)
	{
		b2Assert(obj);
		GetListNode()->InsertAfter(obj->GetListNode());
	}

	/// Insert this object before the specified object.
	void InsertBefore(T* const obj)
	{
		b2Assert(obj);
		GetListNode()->InsertBefore(obj->GetListNode());
	}

	/// Get the next object in the list.
	/// Check against GetTerminator() before deferencing the object.
	T* GetNext() const
	{
		return GetInstanceFromListNode(GetListNode()->GetNext());
	}

	/// Get the previous object in the list.
	/// Check against GetTerminator() before deferencing the object.
	T* GetPrevious() const
	{
		return GetInstanceFromListNode(GetListNode()->GetPrevious());
	}

	/// Get the terminator of the list.
	/// This should not be dereferenced as it is a pointer to
	/// b2TypedIntrusiveListNode<T> *not* T.
	T* GetTerminator() const
	{
		return (T*)GetListNode();
	}

	/// Remove this object from the list it's currently in.
	T* Remove()
	{
		GetListNode()->Remove();
		return GetInstanceFromListNode(GetListNode());
	}

	/// Determine whether this object is in a list.
	bool InList() const
	{
		return GetListNode()->InList();
	}

	// Determine whether this list is empty.
	bool IsEmpty() const
	{
		return GetListNode()->IsEmpty();
	}

	/// Calculate the length of the list.
	uint32 GetLength() const
	{
		return GetListNode()->GetLength();
	}

	B2_INTRUSIVE_LIST_GET_NODE(m_node);

private:
	// Node within an intrusive list.
	b2IntrusiveListNode m_node;

public:
	/// Get a pointer to the instance of T that contains "node".
	static T* GetInstanceFromListNode(b2IntrusiveListNode* const node)
	{
		b2Assert(node);
		// Calculate the pointer to T from the offset.
		return (T*)((uint8*)node - GetNodeOffset(node));
	}

private:
	// Get the offset of m_node within this class.
	static int32 GetNodeOffset(b2IntrusiveListNode* const node)
	{
		b2Assert(node);
		// Perform some type punning to calculate the offset of m_node in T.
		// WARNING: This could result in undefined behavior with some C++
		// compilers.
		T* obj = (T*)node;
		int32 nodeOffset = (int32)((uint8*)&obj->m_node - (uint8*)obj);
		return nodeOffset;
	}
};

// end of IntrusiveList.h

/// Alignment (in bytes) of user memory associated with b2TrackedBlock.
const int32 b2_mallocAlignment = 32;

/// Allocated block of memory that can be tracked in a b2IntrusiveList.
class b2TrackedBlock : public b2TypedIntrusiveListNode<b2TrackedBlock>
{
private:
	// Initialize this block with a reference to "this".
	b2TrackedBlock();
	// Remove the block from the list.
	~b2TrackedBlock() { }

public:
	/// Get the allocated memory associated with this block.
	void* GetMemory() const;

private:
	// Padding required to align the pointer to user memory in the block
	// to b2_mallocAlignment.
	uint8 m_padding[b2_mallocAlignment + sizeof(b2TrackedBlock**)];

public:
	/// Allocate a b2TrackedBlock returning a pointer to memory of size
	/// bytes that can be used by the caller.
	static void* Allocate(uint32 size);

	/// Get a b2TrackedBlock from a pointer to memory returned by
	/// b2TrackedBlock::Allocate().
	static b2TrackedBlock* GetFromMemory(void *memory);

	/// Free a block of memory returned by b2TrackedBlock::Allocate()
	static void Free(void *memory);

	/// Free a b2TrackedBlock.
	static void Free(b2TrackedBlock *block);
};

/// Allocator of blocks which are tracked in a list.
class b2TrackedBlockAllocator
{
public:
	/// Initialize.
	b2TrackedBlockAllocator() {}
	/// Free all allocated blocks.
	~b2TrackedBlockAllocator()
	{
		FreeAll();
	}

	/// Allocate a block of size bytes using b2TrackedBlock::Allocate().
	void* Allocate(uint32 size);

	/// Free a block returned by Allocate().
	void Free(void *memory);

	/// Free all allocated blocks.
	void FreeAll();

	// Get the list of allocated blocks.
	const b2TypedIntrusiveListNode<b2TrackedBlock>& GetList() const
	{
		return m_blocks;
	}

private:
	b2TypedIntrusiveListNode<b2TrackedBlock> m_blocks;
};

// end of TrackedBlock.h


// Initialize this block with a reference to "this".
b2TrackedBlock::b2TrackedBlock()
{
	b2TrackedBlock** pointerToThis =
		(b2TrackedBlock**)((uint8*)GetMemory() - sizeof(b2TrackedBlock**));
	*pointerToThis = this;
}

/// Get the allocated memory associated with this block.
void* b2TrackedBlock::GetMemory() const
{
	// The size of data in this without padding.
	static const uint32 kSizeOfThisWithNoPadding =
		sizeof(*this) - sizeof(m_padding) + sizeof(b2TrackedBlock**);

	// Make sure b2_mallocAlignment is base2.
	b2Assert(((b2_mallocAlignment - 1) & b2_mallocAlignment) == 0);

	// Round the pointer following data in this to b2_mallocAlignment.
	uint8* const aligned = (uint8*)(
		((uintptr_t)this + kSizeOfThisWithNoPadding + b2_mallocAlignment - 1) &
		~((uintptr_t)b2_mallocAlignment - 1));
	// Verify offset doesn't overlap data in this.
	b2Assert((uintptr_t)aligned - (uintptr_t)this >= kSizeOfThisWithNoPadding);
	return aligned;
}

/// Allocate a b2TrackedBlock returning a pointer to memory of size
/// bytes that can be used by the caller.
void* b2TrackedBlock::Allocate(uint32 size)
{
	void* memory = (b2TrackedBlock*)b2Alloc(sizeof(b2TrackedBlock) +
											size);
	if (!memory)
	{
		return NULL;
	}
	return (new(memory) b2TrackedBlock)->GetMemory();
}

/// Get a b2TrackedBlock from a pointer to memory returned by
/// b2TrackedBlock::Allocate().
b2TrackedBlock* b2TrackedBlock::GetFromMemory(void *memory)
{
	uint8* const aligned = (uint8*)memory;
	b2Assert(memory);
	b2TrackedBlock **blockPtr = (b2TrackedBlock**)(aligned -
												   sizeof(b2TrackedBlock**));
	b2Assert(*blockPtr);
	return *blockPtr;
}

/// Free a block of memory returned by b2TrackedBlock::Allocate()
void b2TrackedBlock::Free(void *memory)
{
	Free(GetFromMemory(memory));
}

/// Free a b2TrackedBlock.
void b2TrackedBlock::Free(b2TrackedBlock *block)
{
	b2Assert(block);
	block->~b2TrackedBlock();
	b2Free(block);
}

/// Allocate a block of size bytes using b2TrackedBlock::Allocate().
void* b2TrackedBlockAllocator::Allocate(uint32 size)
{
	void *memory = b2TrackedBlock::Allocate(size);
	m_blocks.InsertBefore(b2TrackedBlock::GetFromMemory(memory));
	return memory;
}

/// Free a block returned by Allocate().
void b2TrackedBlockAllocator::Free(void *memory)
{
	b2TrackedBlock::Free(memory);
}

/// Free all allocated blocks.
void b2TrackedBlockAllocator::FreeAll()
{
	while (!m_blocks.IsEmpty())
	{
		b2TrackedBlock::Free(m_blocks.GetNext());
	}
}

// end of TrackedBlock.cpp

const int32 b2_stackSize = 100 * 1024;	// 100k
const int32 b2_maxStackEntries = 32;

struct b2StackEntry
{
	char* data;
	int32 size;
	bool usedMalloc;
};

// This is a stack allocator used for fast per step allocations.
// You must nest allocate/free pairs. The code will assert
// if you try to interleave multiple allocate/free pairs.
class b2StackAllocator
{
public:
	enum { MIN_ALIGNMENT = sizeof(void*) }; // Must be a power of 2
	enum { ALIGN_MASK = MIN_ALIGNMENT - 1 };

	b2StackAllocator();
	~b2StackAllocator();

	void* Allocate(int32 size);
	void* Reallocate(void* p, int32 size);
	void Free(void* p);

	int32 GetMaxAllocation() const;

private:

	char m_data[b2_stackSize];
	int32 m_index;

	int32 m_allocation;
	int32 m_maxAllocation;

	b2StackEntry m_entries[b2_maxStackEntries];
	int32 m_entryCount;
};

// end of StackAllocator.h

b2StackAllocator::b2StackAllocator()
{
	m_index = 0;
	m_allocation = 0;
	m_maxAllocation = 0;
	m_entryCount = 0;
}

b2StackAllocator::~b2StackAllocator()
{
	b2Assert(m_index == 0);
	b2Assert(m_entryCount == 0);
}

void* b2StackAllocator::Allocate(int32 size)
{
	b2Assert(m_entryCount < b2_maxStackEntries);
	const int32 roundedSize = (size + ALIGN_MASK) & ~ALIGN_MASK;
	b2StackEntry* entry = m_entries + m_entryCount;
	entry->size = roundedSize;
	if (m_index + roundedSize > b2_stackSize)
	{
		entry->data = (char*)b2Alloc(roundedSize);
		entry->usedMalloc = true;
	}
	else
	{
		entry->data = m_data + m_index;
		entry->usedMalloc = false;
		m_index += roundedSize;
	}

	m_allocation += roundedSize;
	m_maxAllocation = b2Max(m_maxAllocation, m_allocation);
	++m_entryCount;

	return entry->data;
}

void* b2StackAllocator::Reallocate(void* p, int32 size)
{
	b2Assert(m_entryCount > 0);
	b2StackEntry* entry = m_entries + m_entryCount - 1;
	b2Assert(p == entry->data);
	B2_NOT_USED(p);
	int32 incrementSize = size - entry->size;
	if (incrementSize > 0)
	{
		if (entry->usedMalloc)
		{
			void* data = b2Alloc(size);
			memcpy(data, entry->data, entry->size);
			b2Free(entry->data);
			entry->data = (char*)data;
		}
		else if (m_index + incrementSize > b2_stackSize)
		{
			void* data = b2Alloc(size);
			memcpy(data, entry->data, entry->size);
			m_index -= entry->size;
			entry->data = (char*)data;
			entry->usedMalloc = true;
		}
		else
		{
			m_index += incrementSize;
			m_allocation += incrementSize;
			m_maxAllocation = b2Max(m_maxAllocation, m_allocation);
		}
		entry->size = size;
	}

	return entry->data;
}

void b2StackAllocator::Free(void* p)
{
	b2Assert(m_entryCount > 0);
	b2StackEntry* entry = m_entries + m_entryCount - 1;
	b2Assert(p == entry->data);
	if (entry->usedMalloc)
	{
		b2Free(p);
	}
	else
	{
		m_index -= entry->size;
	}
	m_allocation -= entry->size;
	--m_entryCount;

	p = NULL;
}

int32 b2StackAllocator::GetMaxAllocation() const
{
	return m_maxAllocation;
}

// end of StackAllocater.cpp

/// When B2_FREE_LIST_CHECK_ALLOCATED_ON_FREE is 1, b2FreeList::Free() will
/// check that the deallocated node was allocated from the freelist.
#ifndef B2_FREE_LIST_CHECK_ALLOCATED_ON_FREE
#define B2_FREE_LIST_CHECK_ALLOCATED_ON_FREE 0
#endif // B2_FREE_LIST_CHECK_ALLOCATED_ON_FREE


/// Fast - O(1) - list based allocator for items that can be inserted into
/// b2IntrusiveListNode lists.
class b2FreeList
{
public:
	/// Construct the free list.
	b2FreeList() { }

	/// Destroy the free list.
	~b2FreeList() { }

	/// Allocate an item from the freelist.
	b2IntrusiveListNode* Allocate();

	/// Free an item from the freelist.
	void Free(b2IntrusiveListNode* node);

	/// Add an item to the freelist so that it can be allocated using
	/// b2FreeList::Allocate().
	void AddToFreeList(b2IntrusiveListNode* node);

	/// Remove all items (allocated and free) from the freelist.
	void RemoveAll();

	/// Get the list which tracks allocated items.
	const b2IntrusiveListNode& GetAllocatedList() const {
		return m_allocated;
	}

	/// Get the list which tracks free items.
	const b2IntrusiveListNode& GetFreeList() const {
		return m_free;
	}

protected:
	/// List of allocated items.
	b2IntrusiveListNode m_allocated;
	/// List of free items.
	b2IntrusiveListNode m_free;
};


/// Typed b2FreeList which manages items of type T assuming T implements
/// the GetInstanceFromListNode() and GetListNode() methods.
template<typename T>
class b2TypedFreeList {
public:
	/// Construct the free list.
	b2TypedFreeList() { }

	/// Destroy the free list.
	~b2TypedFreeList() { }

	/// Allocate an item from the free list.
	T* Allocate() {
		b2IntrusiveListNode* const node = m_freeList.Allocate();
		if (!node) return NULL;
		return T::GetInstanceFromListNode(node);
	}

	/// Free an item.
	void Free(T* instance) {
		b2Assert(instance);
		m_freeList.Free(instance->GetListNode());
	}

	/// Add an item to the freelist so that it can be allocated with
	/// b2TypedFreeList::Allocate().
	void AddToFreeList(T* instance)
	{
		b2Assert(instance);
		m_freeList.AddToFreeList(instance->GetListNode());
	}

	// Get the underlying b2FreeList.
	b2FreeList* GetFreeList() { return &m_freeList; }
	const b2FreeList* GetFreeList() const { return &m_freeList; }

protected:
	b2FreeList m_freeList;
};

// end of FreeList.h

/// Allocate an item from the freelist.
b2IntrusiveListNode* b2FreeList::Allocate()
{
	if (m_free.IsEmpty()) return NULL;
	b2IntrusiveListNode * const node = m_free.GetNext();
	node->Remove();
	m_allocated.InsertBefore(node);
	return node;
}

void b2FreeList::Free(b2IntrusiveListNode* node)
{
	b2Assert(node);
#if B2_FREE_LIST_CHECK_ALLOCATED_ON_FREE
	b2Assert(m_allocated.FindNodeInList(node));
#endif // B2_FREE_LIST_CHECK_ALLOCATED_ON_FREE
	node->Remove();
	m_free.InsertAfter(node);
}

void b2FreeList::AddToFreeList(b2IntrusiveListNode* node)
{
	b2Assert(node);
	b2Assert(!node->InList());
	m_free.InsertBefore(node);
}

void b2FreeList::RemoveAll()
{
	while (!m_allocated.IsEmpty()) {
		m_allocated.GetNext()->Remove();
	}
	while (!m_free.IsEmpty()) {
		m_free.GetNext()->Remove();
	}
}

// end of FreeList.cpp

/// Freelist based allocator for fixed sized items from slabs (memory
/// preallocated from the heap).
/// T should be a class which has a default constructor and implements the
/// member function "b2IntrusiveList* GetListNode()".
/// All objects in a slab are constructed when a slab is created and destructed
/// when a slab is freed.
template<typename T>
class b2SlabAllocator
{
private:
	// Information about a slab.
	class Slab
	{
	public:
		/// Initialize a slab with the number of items it contains.
		Slab(uint32 numberOfItems) :
			m_numberOfItems(numberOfItems)
		{
			B2_NOT_USED(m_padding);
			// This assumes that this class is packed on at least a 4-byte
			// boundary with no padding.  Verify the assumption.
			b2Assert(sizeof(*this) == b2_mallocAlignment);
		}

		/// Empty destructor.
		~Slab() { }

		/// Get the number of items in this slab.
		uint32 GetNumberOfItems() const { return m_numberOfItems; }

		/// Get a pointer to the first item in the slab.
		T* GetFirstItem() const
		{
			return (T*)((uint8*)(this + 1));
		}

		/// Get a pointer to the end of the slab.
		/// NOTE: This is a pointer after the last byte of the slab not the
		/// last item in the slab.
		T* GetItemEnd() const { return GetFirstItem() + GetNumberOfItems(); }

	private:
		/// Number of items in the slab.
		uint32 m_numberOfItems;
		/// Padding to align the first item in the slab to b2_mallocAlignment.
		uint8 m_padding[b2_mallocAlignment - sizeof(uint32)];
	};

public:
	/// Initialize the allocator to allocate itemsPerSlab of type T for each
	/// slab that is allocated.
	b2SlabAllocator(const uint32 itemsPerSlab) :
		m_itemsPerSlab(itemsPerSlab)
	{
	}

	/// Free all allocated slabs.
	~b2SlabAllocator()
	{
		FreeAllSlabs();
	}

	/// Set size of the next allocated slab using the number of items per
	/// slab.  Setting this value to zero disables further slab allocation.
	void SetItemsPerSlab(uint32 itemsPerSlab)
	{
		m_itemsPerSlab = itemsPerSlab;
	}

	// Get the size of the next allocated slab.
	uint32 GetItemsPerSlab() const
	{
		return m_itemsPerSlab;
	}

	/// Allocate a item from the slab.
	T* Allocate()
	{
		// Allocate a slab if needed here.
		if (m_freeList.GetFreeList()->GetFreeList().IsEmpty() &&
			!AllocateSlab())
			return NULL;
		return m_freeList.Allocate();
	}

	/// Free an item from the slab.
	void Free(T *object)
	{
		m_freeList.Free(object);
	}

	/// Allocate a slab, construct instances of T and add them to the free
	/// pool.
	bool AllocateSlab()
	{
		if (!m_itemsPerSlab) return false;
		const uint32 slabSize = sizeof(Slab) + (sizeof(T) * m_itemsPerSlab);
		void* const memory = m_slabs.Allocate(slabSize);
		if (!memory) return false;

		Slab* const slab = new (BlockGetSlab(memory)) Slab(m_itemsPerSlab);
		T* item = slab->GetFirstItem();
		for (uint32 i = 0; i < m_itemsPerSlab; ++i, ++item)
		{
			m_freeList.AddToFreeList(new (item) T);
		}
		return true;
	}

	/// Free all slabs.
	void FreeAllSlabs()
	{
		const b2TypedIntrusiveListNode<b2TrackedBlock>& slabList =
			m_slabs.GetList();
		while (!slabList.IsEmpty())
		{
			FreeSlab(BlockGetSlab(slabList.GetNext()->GetMemory()));
		}
	}

	/// Free all empty slabs.
	/// This method is slow - O(M^N) - since this class doesn't track
	/// the association between each item and slab.
	void FreeEmptySlabs()
	{
		const b2IntrusiveListNode& freeItemList =
			m_freeList.GetFreeList()->GetFreeList();
		const b2IntrusiveListNode* freeItemListTerminator =
			freeItemList.GetTerminator();
		const b2TypedIntrusiveListNode<b2TrackedBlock>& slabList =
			m_slabs.GetList();
		const b2TypedIntrusiveListNode<b2TrackedBlock>* slabListTerminator =
			slabList.GetTerminator();
		b2TrackedBlock* block = slabList.GetNext();
		while (block != slabListTerminator)
		{
			// Get the Slab from the memory associated with the block.
			Slab* const slab = BlockGetSlab(block->GetMemory());
			block = block->GetNext();

			// Determine the range of memory the Slab owns.
			const uint8* const slabItemStart = (uint8*)slab->GetFirstItem();
			const uint8* const slabItemEnd = (uint8*)slab->GetItemEnd();

			// Count all free items that are owned by the current slab.
			uint8 freeItems = 0;
			bool empty = false;
			for (b2IntrusiveListNode* itemNode = freeItemList.GetNext();
				 itemNode != freeItemListTerminator;
				 itemNode = itemNode->GetNext())
			{
				const uint8* itemNodeAddress = (uint8*)itemNode;
				if (itemNodeAddress >= slabItemStart &&
					itemNodeAddress <= slabItemEnd)
				{
					++freeItems;
					if (slab->GetNumberOfItems() == freeItems)
					{
						empty = true;
						break;
					}
				}
			}
			// If a slab is empty, free it.
			if (empty)
			{
				FreeSlab(slab);
			}
		}
	}

	/// Get the item allocator freelist.
	const b2TypedFreeList<T>& GetFreeList() const
	{
		return m_freeList;
	}

private:
	/// Destroy all objects in a slab and free the slab.
	void FreeSlab(Slab * const slab)
	{
		b2Assert(slab);
		const uint32 numberOfItems = slab->GetNumberOfItems();
		T* item = slab->GetFirstItem();
		for (uint32 i = 0; i < numberOfItems; ++i, ++item)
		{
			item->~T();
		}
		slab->~Slab();
		m_slabs.Free(slab);
	}

	/// Get a pointer to a Slab from a block of memory in m_slabs.
	Slab* BlockGetSlab(void *memory)
	{
		return (Slab*)memory;
	}

	/// Get a pointer to the first item in the array of items referenced by a
	/// Slab.
	T* SlabGetFirstItem(Slab* slab)
	{
		return (T*)(slab + 1);
	}

private:
	/// Contains a list of b2TrackedBlock instances where each b2TrackedBlock's
	/// associated user memory contains a Slab followed by instances of T.
	b2TrackedBlockAllocator m_slabs;
	/// Number of items to allocate in the next allocated slab.
	uint32 m_itemsPerSlab;
	/// Freelist which contains instances of T.
	b2TypedFreeList<T> m_freeList;
};

// end of SlabAllocater.h

const int32 b2_chunkSize = 16 * 1024;
const int32 b2_maxBlockSize = 640;
const int32 b2_blockSizes = 14;
const int32 b2_chunkArrayIncrement = 128;

struct b2Block;
struct b2Chunk;

/// This is a small object allocator used for allocating small
/// objects that persist for more than one time step.
/// See: http://www.codeproject.com/useritems/Small_Block_Allocator.asp
class b2BlockAllocator
{
public:
	b2BlockAllocator();
	~b2BlockAllocator();

	/// Allocate memory. This uses b2Alloc if the size is larger than b2_maxBlockSize.
	void* Allocate(int32 size);

	/// Free memory. This uses b2Free if the size is larger than b2_maxBlockSize.
	void Free(void* p, int32 size);

	void Clear();

	/// Returns the number of allocations larger than the max block size.
	uint32 GetNumGiantAllocations() const;

private:
	b2Chunk* m_chunks;
	int32 m_chunkCount;
	int32 m_chunkSpace;

	b2Block* m_freeLists[b2_blockSizes];

	// Record giant allocations--ones bigger than the max block size
	b2TrackedBlockAllocator m_giants;

	static int32 s_blockSizes[b2_blockSizes];
	static uint8 s_blockSizeLookup[b2_maxBlockSize + 1];
	static bool s_blockSizeLookupInitialized;
};

// end of BlockAllocater.h

int32 b2BlockAllocator::s_blockSizes[b2_blockSizes] =
{
	16,		// 0
	32,		// 1
	64,		// 2
	96,		// 3
	128,	// 4
	160,	// 5
	192,	// 6
	224,	// 7
	256,	// 8
	320,	// 9
	384,	// 10
	448,	// 11
	512,	// 12
	640,	// 13
};
uint8 b2BlockAllocator::s_blockSizeLookup[b2_maxBlockSize + 1];
bool b2BlockAllocator::s_blockSizeLookupInitialized;

struct b2Chunk
{
	int32 blockSize;
	b2Block* blocks;
};

struct b2Block
{
	b2Block* next;
};

b2BlockAllocator::b2BlockAllocator()
{
	b2Assert((uint32)b2_blockSizes < UCHAR_MAX);

	m_chunkSpace = b2_chunkArrayIncrement;
	m_chunkCount = 0;
	m_chunks = (b2Chunk*)b2Alloc(m_chunkSpace * sizeof(b2Chunk));

	memset(m_chunks, 0, m_chunkSpace * sizeof(b2Chunk));
	memset(m_freeLists, 0, sizeof(m_freeLists));

	if (s_blockSizeLookupInitialized == false)
	{
		int32 j = 0;
		for (int32 i = 1; i <= b2_maxBlockSize; ++i)
		{
			b2Assert(j < b2_blockSizes);
			if (i <= s_blockSizes[j])
			{
				s_blockSizeLookup[i] = (uint8)j;
			}
			else
			{
				++j;
				s_blockSizeLookup[i] = (uint8)j;
			}
		}

		s_blockSizeLookupInitialized = true;
	}
}

b2BlockAllocator::~b2BlockAllocator()
{
	for (int32 i = 0; i < m_chunkCount; ++i)
	{
		b2Free(m_chunks[i].blocks);
	}

	b2Free(m_chunks);
}

uint32 b2BlockAllocator::GetNumGiantAllocations() const
{
	return m_giants.GetList().GetLength();
}

void* b2BlockAllocator::Allocate(int32 size)
{
	if (size == 0)
		return NULL;

	b2Assert(0 < size);

	if (size > b2_maxBlockSize)
	{
		return m_giants.Allocate(size);
	}

	int32 index = s_blockSizeLookup[size];
	b2Assert(0 <= index && index < b2_blockSizes);

	if (m_freeLists[index])
	{
		b2Block* block = m_freeLists[index];
		m_freeLists[index] = block->next;
		return block;
	}
	else
	{
		if (m_chunkCount == m_chunkSpace)
		{
			b2Chunk* oldChunks = m_chunks;
			m_chunkSpace += b2_chunkArrayIncrement;
			m_chunks = (b2Chunk*)b2Alloc(m_chunkSpace * sizeof(b2Chunk));
			memcpy(m_chunks, oldChunks, m_chunkCount * sizeof(b2Chunk));
			memset(m_chunks + m_chunkCount, 0, b2_chunkArrayIncrement * sizeof(b2Chunk));
			b2Free(oldChunks);
		}

		b2Chunk* chunk = m_chunks + m_chunkCount;
		chunk->blocks = (b2Block*)b2Alloc(b2_chunkSize);
#if DEBUG
		memset(chunk->blocks, 0xcd, b2_chunkSize);
#endif
		int32 blockSize = s_blockSizes[index];
		chunk->blockSize = blockSize;
		int32 blockCount = b2_chunkSize / blockSize;
		b2Assert(blockCount * blockSize <= b2_chunkSize);
		for (int32 i = 0; i < blockCount - 1; ++i)
		{
			b2Block* block = (b2Block*)((int8*)chunk->blocks + blockSize * i);
			b2Block* next = (b2Block*)((int8*)chunk->blocks + blockSize * (i + 1));
			block->next = next;
		}
		b2Block* last = (b2Block*)((int8*)chunk->blocks + blockSize * (blockCount - 1));
		last->next = NULL;

		m_freeLists[index] = chunk->blocks->next;
		++m_chunkCount;

		return chunk->blocks;
	}
}

void b2BlockAllocator::Free(void* p, int32 size)
{
	if (size == 0)
	{
		return;
	}

	b2Assert(0 < size);

	if (size > b2_maxBlockSize)
	{
		m_giants.Free(p);
		return;
	}

	int32 index = s_blockSizeLookup[size];
	b2Assert(0 <= index && index < b2_blockSizes);

#if B2_ASSERT_ENABLED
	// Verify the memory address and size is valid.
	int32 blockSize = s_blockSizes[index];
	bool found = false;
	for (int32 i = 0; i < m_chunkCount; ++i)
	{
		b2Chunk* chunk = m_chunks + i;
		if (chunk->blockSize != blockSize)
		{
			b2Assert(	(int8*)p + blockSize <= (int8*)chunk->blocks ||
						(int8*)chunk->blocks + b2_chunkSize <= (int8*)p);
		}
		else
		{
			if ((int8*)chunk->blocks <= (int8*)p && (int8*)p + blockSize <= (int8*)chunk->blocks + b2_chunkSize)
			{
				found = true;
			}
		}
	}

	b2Assert(found);
#endif // B2_ASSERT_ENABLED

#if DEBUG
	memset(p, 0xfd, s_blockSizes[index]);
#endif

	b2Block* block = (b2Block*)p;
	block->next = m_freeLists[index];
	m_freeLists[index] = block;
}

void b2BlockAllocator::Clear()
{
	for (int32 i = 0; i < m_chunkCount; ++i)
	{
		b2Free(m_chunks[i].blocks);
	}

	m_chunkCount = 0;
	memset(m_chunks, 0, m_chunkSpace * sizeof(b2Chunk));

	memset(m_freeLists, 0, sizeof(m_freeLists));
}

// end of BlockAllocater.cpp

/// A simple array-like container, similar to std::vector.
/// If we ever start using stl, we should replace this with std::vector.
template <typename T>
class b2GrowableBuffer
{
public:
	b2GrowableBuffer(b2BlockAllocator& allocator) :
		data(NULL),
		count(0),
		capacity(0),
		allocator(&allocator)
	{
	#if defined(LIQUIDFUN_SIMD_NEON)
		// b2ParticleAssembly.neon.s assumes these values are at fixed offsets.
        // If this assert fails, be sure to update the assembly offsets!
		// ldr r3, [r9, #0] @ r3 = out = contacts.data
        // ldr r6, [r9, #8] @ r6 = contacts.capacity
		b2Assert((intptr_t)&data - (intptr_t)(this) == 0
			  && (intptr_t)&capacity - (intptr_t)(this) == 8);
	#endif // defined(LIQUIDFUN_SIMD_NEON)
	}

	b2GrowableBuffer(const b2GrowableBuffer<T>& rhs) :
		data(NULL),
		count(rhs.count),
		capacity(rhs.capacity),
		allocator(rhs.allocator)
	{
		if (rhs.data != NULL)
		{
			data = (T*) allocator->Allocate(sizeof(T) * capacity);
			memcpy(data, rhs.data, sizeof(T) * count);
		}
	}

	~b2GrowableBuffer()
	{
		Free();
	}

	T& Append()
	{
		if (count >= capacity)
		{
			Grow();
		}
		return data[count++];
	}

	void Reserve(int32 newCapacity)
	{
		if (capacity >= newCapacity)
			return;

		// Reallocate and copy.
		T* newData = (T*) allocator->Allocate(sizeof(T) * newCapacity);
		if (data)
		{
			memcpy(newData, data, sizeof(T) * count);
			allocator->Free(data, sizeof(T) * capacity);
		}

		// Update pointer and capacity.
		capacity = newCapacity;
		data = newData;
	}

	void Grow()
	{
		// Double the capacity.
		int32 newCapacity = capacity ? 2 * capacity
						  : b2_minParticleSystemBufferCapacity;
		b2Assert(newCapacity > capacity);
		Reserve(newCapacity);
	}

	void Free()
	{
		if (data == NULL)
			return;

		allocator->Free(data, sizeof(data[0]) * capacity);
		data = NULL;
		capacity = 0;
		count = 0;
	}

	void Shorten(const T* newEnd)
	{
		b2Assert(newEnd >= data);
		count = (int32) (newEnd - data);
	}

	T& operator[](int i)
	{
		return data[i];
	}

	const T& operator[](int i) const
	{
		return data[i];
	}

	T* Data()
	{
		return data;
	}

	const T* Data() const
	{
		return data;
	}

	T* Begin()
	{
		return data;
	}

	const T* Begin() const
	{
		return data;
	}

	T* End()
	{
		return &data[count];
	}

	const T* End() const
	{
		return &data[count];
	}

	int32 GetCount() const
	{
		return count;
	}

	void SetCount(int32 newCount)
	{
		b2Assert(0 <= newCount && newCount <= capacity);
		count = newCount;
	}

	int32 GetCapacity() const
	{
		return capacity;
	}

	template<class UnaryPredicate>
	T* RemoveIf(UnaryPredicate pred)
	{
		T* newEnd = std::remove_if(data, data + count, pred);
		Shorten(newEnd);
		return newEnd;
	}

	template<class BinaryPredicate>
	T* Unique(BinaryPredicate pred)
	{
		T* newEnd = std::unique(data, data + count, pred);
		Shorten(newEnd);
		return newEnd;
	}

private:
	T* data;
	int32 count;
	int32 capacity;
	b2BlockAllocator* allocator;
};

// end of GrowableBuffer.h

/// This is a growable LIFO stack with an initial capacity of N.
/// If the stack size exceeds the initial capacity, the heap is used
/// to increase the size of the stack.
template <typename T, int32 N>
class b2GrowableStack
{
public:
	b2GrowableStack()
	{
		m_stack = m_array;
		m_count = 0;
		m_capacity = N;
	}

	~b2GrowableStack()
	{
		if (m_stack != m_array)
		{
			b2Free(m_stack);
			m_stack = NULL;
		}
	}

	void Push(const T& element)
	{
		if (m_count == m_capacity)
		{
			T* old = m_stack;
			m_capacity *= 2;
			m_stack = (T*)b2Alloc(m_capacity * sizeof(T));
			memcpy(m_stack, old, m_count * sizeof(T));
			if (old != m_array)
			{
				b2Free(old);
			}
		}

		m_stack[m_count] = element;
		++m_count;
	}

	T Pop()
	{
		b2Assert(m_count > 0);
		--m_count;
		return m_stack[m_count];
	}

	int32 GetCount()
	{
		return m_count;
	}

private:
	T* m_stack;
	T m_array[N];
	int32 m_count;
	int32 m_capacity;
};

// end of GrowableStack

const b2Vec2 b2Vec2_zero(0.0f, 0.0f);

/// Solve A * x = b, where b is a column vector. This is more efficient
/// than computing the inverse in one-shot cases.
b2Vec3 b2Mat33::Solve33(const b2Vec3& b) const
{
	float32 det = b2Dot(ex, b2Cross(ey, ez));
	if (det != 0.0f)
	{
		det = 1.0f / det;
	}
	b2Vec3 x;
	x.x = det * b2Dot(b, b2Cross(ey, ez));
	x.y = det * b2Dot(ex, b2Cross(b, ez));
	x.z = det * b2Dot(ex, b2Cross(ey, b));
	return x;
}

/// Solve A * x = b, where b is a column vector. This is more efficient
/// than computing the inverse in one-shot cases.
b2Vec2 b2Mat33::Solve22(const b2Vec2& b) const
{
	float32 a11 = ex.x, a12 = ey.x, a21 = ex.y, a22 = ey.y;
	float32 det = a11 * a22 - a12 * a21;
	if (det != 0.0f)
	{
		det = 1.0f / det;
	}
	b2Vec2 x;
	x.x = det * (a22 * b.x - a12 * b.y);
	x.y = det * (a11 * b.y - a21 * b.x);
	return x;
}

///
void b2Mat33::GetInverse22(b2Mat33* M) const
{
	float32 a = ex.x, b = ey.x, c = ex.y, d = ey.y;
	float32 det = a * d - b * c;
	if (det != 0.0f)
	{
		det = 1.0f / det;
	}

	M->ex.x =  det * d;	M->ey.x = -det * b; M->ex.z = 0.0f;
	M->ex.y = -det * c;	M->ey.y =  det * a; M->ey.z = 0.0f;
	M->ez.x = 0.0f; M->ez.y = 0.0f; M->ez.z = 0.0f;
}

/// Returns the zero matrix if singular.
void b2Mat33::GetSymInverse33(b2Mat33* M) const
{
	float32 det = b2Dot(ex, b2Cross(ey, ez));
	if (det != 0.0f)
	{
		det = 1.0f / det;
	}

	float32 a11 = ex.x, a12 = ey.x, a13 = ez.x;
	float32 a22 = ey.y, a23 = ez.y;
	float32 a33 = ez.z;

	M->ex.x = det * (a22 * a33 - a23 * a23);
	M->ex.y = det * (a13 * a23 - a12 * a33);
	M->ex.z = det * (a12 * a23 - a13 * a22);

	M->ey.x = M->ex.y;
	M->ey.y = det * (a11 * a33 - a13 * a13);
	M->ey.z = det * (a13 * a12 - a11 * a23);

	M->ez.x = M->ex.z;
	M->ez.y = M->ey.z;
	M->ez.z = det * (a11 * a22 - a12 * a12);
}

// end of Math.cpp

struct b2Color;
class b2ParticleGroup;

/// @file

/// The particle type. Can be combined with the | operator.
enum b2ParticleFlag
{
	/// Water particle.
	b2_waterParticle = 0,
	/// Removed after next simulation step.
	b2_zombieParticle = 1 << 1,
	/// Zero velocity.
	b2_wallParticle = 1 << 2,
	/// With restitution from stretching.
	b2_springParticle = 1 << 3,
	/// With restitution from deformation.
	b2_elasticParticle = 1 << 4,
	/// With viscosity.
	b2_viscousParticle = 1 << 5,
	/// Without isotropic pressure.
	b2_powderParticle = 1 << 6,
	/// With surface tension.
	b2_tensileParticle = 1 << 7,
	/// Mix color between contacting particles.
	b2_colorMixingParticle = 1 << 8,
	/// Call b2DestructionListener on destruction.
	b2_destructionListenerParticle = 1 << 9,
	/// Prevents other particles from leaking.
	b2_barrierParticle = 1 << 10,
	/// Less compressibility.
	b2_staticPressureParticle = 1 << 11,
	/// Makes pairs or triads with other particles.
	b2_reactiveParticle = 1 << 12,
	/// With high repulsive force.
	b2_repulsiveParticle = 1 << 13,
	/// Call b2ContactListener when this particle is about to interact with
	/// a rigid body or stops interacting with a rigid body.
	/// This results in an expensive operation compared to using
	/// b2_fixtureContactFilterParticle to detect collisions between
	/// particles.
	b2_fixtureContactListenerParticle = 1 << 14,
	/// Call b2ContactListener when this particle is about to interact with
	/// another particle or stops interacting with another particle.
	/// This results in an expensive operation compared to using
	/// b2_particleContactFilterParticle to detect collisions between
	/// particles.
	b2_particleContactListenerParticle = 1 << 15,
	/// Call b2ContactFilter when this particle interacts with rigid bodies.
	b2_fixtureContactFilterParticle = 1 << 16,
	/// Call b2ContactFilter when this particle interacts with other
	/// particles.
	b2_particleContactFilterParticle = 1 << 17,
};

/// Small color object for each particle
class b2ParticleColor
{
public:
	b2ParticleColor() {}
	/// Constructor with four elements: r (red), g (green), b (blue), and a
	/// (opacity).
	/// Each element can be specified 0 to 255.
	b2Inline b2ParticleColor(uint8 r, uint8 g, uint8 b, uint8 a)
	{
		Set(r, g, b, a);
	}

	/// Constructor that initializes the above four elements with the value of
	/// the b2Color object.
	b2ParticleColor(const b2Color& color);

	/// True when all four color elements equal 0. When true, a particle color
	/// buffer isn't allocated by CreateParticle().
	///
	bool IsZero() const
	{
		return !r && !g && !b && !a;
	}

	/// Used internally to convert the value of b2Color.
	///
	b2Color GetColor() const;

	/// Sets color for current object using the four elements described above.
	///
	b2Inline void Set(uint8 r_, uint8 g_, uint8 b_, uint8 a_)
	{
		r = r_;
		g = g_;
		b = b_;
		a = a_;
	}

	/// Initializes the object with the value of the b2Color.
	///
	void Set(const b2Color& color);

	/// Assign a b2ParticleColor to this instance.
	b2ParticleColor& operator = (const b2ParticleColor &color)
	{
		Set(color.r, color.g, color.b, color.a);
		return *this;
	}

	/// Multiplies r, g, b, a members by s where s is a value between 0.0
	/// and 1.0.
	b2ParticleColor& operator *= (float32 s)
	{
		Set((uint8)(r * s), (uint8)(g * s), (uint8)(b * s), (uint8)(a * s));
		return *this;
	}

	/// Scales r, g, b, a members by s where s is a value between 0 and 255.
	b2ParticleColor& operator *= (uint8 s)
	{
		// 1..256 to maintain the complete dynamic range.
		const int32 scale = (int32)s + 1;
		Set((uint8)(((int32)r * scale) >> k_bitsPerComponent),
			(uint8)(((int32)g * scale) >> k_bitsPerComponent),
			(uint8)(((int32)b * scale) >> k_bitsPerComponent),
			(uint8)(((int32)a * scale) >> k_bitsPerComponent));
		return *this;
	}

	/// Scales r, g, b, a members by s returning the modified b2ParticleColor.
	b2ParticleColor operator * (float32 s) const
	{
		return MultiplyByScalar(s);
	}

	/// Scales r, g, b, a members by s returning the modified b2ParticleColor.
	b2ParticleColor operator * (uint8 s) const
	{
		return MultiplyByScalar(s);
	}

	/// Add two colors.  This is a non-saturating addition so values
	/// overflows will wrap.
	b2Inline b2ParticleColor& operator += (const b2ParticleColor &color)
	{
		r += color.r;
		g += color.g;
		b += color.b;
		a += color.a;
		return *this;
	}

	/// Add two colors.  This is a non-saturating addition so values
	/// overflows will wrap.
	b2ParticleColor operator + (const b2ParticleColor &color) const
	{
		b2ParticleColor newColor(*this);
		newColor += color;
		return newColor;
	}

	/// Subtract a color from this color.  This is a subtraction without
	/// saturation so underflows will wrap.
	b2Inline b2ParticleColor& operator -= (const b2ParticleColor &color)
	{
		r -= color.r;
		g -= color.g;
		b -= color.b;
		a -= color.a;
		return *this;
	}

	/// Subtract a color from this color returning the result.  This is a
	/// subtraction without saturation so underflows will wrap.
	b2ParticleColor operator - (const b2ParticleColor &color) const
	{
		b2ParticleColor newColor(*this);
		newColor -= color;
		return newColor;
	}

	/// Compare this color with the specified color.
	bool operator == (const b2ParticleColor &color) const
	{
		return r == color.r && g == color.g && b == color.b && a == color.a;
	}

	/// Mix mixColor with this color using strength to control how much of
	/// mixColor is mixed with this color and vice versa.  The range of
	/// strength is 0..128 where 0 results in no color mixing and 128 results
	/// in an equal mix of both colors.  strength 0..128 is analogous to an
	/// alpha channel value between 0.0f..0.5f.
	b2Inline void Mix(b2ParticleColor * const mixColor, const int32 strength)
	{
		MixColors(this, mixColor, strength);
	}

	/// Mix colorA with colorB using strength to control how much of
	/// colorA is mixed with colorB and vice versa.  The range of
	/// strength is 0..128 where 0 results in no color mixing and 128 results
	/// in an equal mix of both colors.  strength 0..128 is analogous to an
	/// alpha channel value between 0.0f..0.5f.
	static b2Inline void MixColors(b2ParticleColor * const colorA,
							 b2ParticleColor * const colorB,
							 const int32 strength)
	{
		const uint8 dr = (uint8)((strength * (colorB->r - colorA->r)) >>
								 k_bitsPerComponent);
		const uint8 dg = (uint8)((strength * (colorB->g - colorA->g)) >>
								 k_bitsPerComponent);
		const uint8 db = (uint8)((strength * (colorB->b - colorA->b)) >>
								 k_bitsPerComponent);
		const uint8 da = (uint8)((strength * (colorB->a - colorA->a)) >>
								 k_bitsPerComponent);
		colorA->r += dr;
		colorA->g += dg;
		colorA->b += db;
		colorA->a += da;
		colorB->r -= dr;
		colorB->g -= dg;
		colorB->b -= db;
		colorB->a -= da;
	}

private:
	/// Generalization of the multiply operator using a scalar in-place
	/// multiplication.
	template <typename T>
	b2ParticleColor MultiplyByScalar(T s) const
	{
		b2ParticleColor color(*this);
		color *= s;
		return color;
	}

public:
	uint8 r, g, b, a;

protected:
	/// Maximum value of a b2ParticleColor component.
	static const float32 k_maxValue;
	/// 1.0 / k_maxValue.
	static const float32 k_inverseMaxValue;
	/// Number of bits used to store each b2ParticleColor component.
	static const uint8 k_bitsPerComponent;
};

extern b2ParticleColor b2ParticleColor_zero;

/// A particle definition holds all the data needed to construct a particle.
/// You can safely re-use these definitions.
struct b2ParticleDef
{
	b2ParticleDef()
	{
		flags = 0;
		position = b2Vec2_zero;
		velocity = b2Vec2_zero;
		color = b2ParticleColor_zero;
		lifetime = 0.0f;
		userData = NULL;
		group = NULL;
	}

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
	/// Set position with direct floats
	void SetPosition(float32 x, float32 y);

	/// Set color with direct ints.
	void SetColor(int32 r, int32 g, int32 b, int32 a);
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

	/// \brief Specifies the type of particle (see #b2ParticleFlag).
	///
	/// A particle may be more than one type.
	/// Multiple types are chained by logical sums, for example:
	/// pd.flags = b2_elasticParticle | b2_viscousParticle
	uint32 flags;

	/// The world position of the particle.
	b2Vec2 position;

	/// The linear velocity of the particle in world co-ordinates.
	b2Vec2 velocity;

	/// The color of the particle.
	b2ParticleColor color;

	/// Lifetime of the particle in seconds.  A value <= 0.0f indicates a
	/// particle with infinite lifetime.
	float32 lifetime;

	/// Use this to store application-specific body data.
	void* userData;

	/// An existing particle group to which the particle will be added.
	b2ParticleGroup* group;

};

/// A helper function to calculate the optimal number of iterations.
int32 b2CalculateParticleIterations(
	float32 gravity, float32 radius, float32 timeStep);

/// Handle to a particle. Particle indices are ephemeral: the same index might
/// refer to a different particle, from frame-to-frame. If you need to keep a
/// reference to a particular particle across frames, you should acquire a
/// b2ParticleHandle. Use #b2ParticleSystem::GetParticleHandleFromIndex() to
/// retrieve the b2ParticleHandle of a particle from the particle system.
class b2ParticleHandle : public b2TypedIntrusiveListNode<b2ParticleHandle>
{
	// Allow b2ParticleSystem to use SetIndex() to associate particle handles
	// with particle indices.
	friend class b2ParticleSystem;

public:
	/// Initialize the index associated with the handle to an invalid index.
	b2ParticleHandle() : m_index(b2_invalidParticleIndex) { }
	/// Empty destructor.
	~b2ParticleHandle() { }

	/// Get the index of the particle associated with this handle.
	int32 GetIndex() const { return m_index; }

private:
	/// Set the index of the particle associated with this handle.
	void SetIndex(int32 index) { m_index = index; }

private:
	// Index of the particle within the particle system.
	int32 m_index;
};

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
inline void b2ParticleDef::SetPosition(float32 x, float32 y)
{
	position.Set(x, y);
}

inline void b2ParticleDef::SetColor(int32 r, int32 g, int32 b, int32 a)
{
	color.Set((uint8)r, (uint8)g, (uint8)b, (uint8)a);
}
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

// end of Particle.h

/// Color for debug drawing. Each value has the range [0,1].
struct b2Color
{
	b2Color() {}
	b2Color(float32 r, float32 g, float32 b) : r(r), g(g), b(b) {}
	void Set(float32 ri, float32 gi, float32 bi) { r = ri; g = gi; b = bi; }
	float32 r, g, b;
};

/// Implement and register this class with a b2World to provide debug drawing of physics
/// entities in your game.
class b2Draw
{
public:
	b2Draw();

	virtual ~b2Draw() {}

	enum
	{
		e_shapeBit				= 0x0001,	///< draw shapes
		e_jointBit				= 0x0002,	///< draw joint connections
		e_aabbBit				= 0x0004,	///< draw axis aligned bounding boxes
		e_pairBit				= 0x0008,	///< draw broad-phase pairs
		e_centerOfMassBit			= 0x0010,	///< draw center of mass frame
		e_particleBit				= 0x0020  ///< draw particles
	};

	/// Set the drawing flags.
	void SetFlags(uint32 flags);

	/// Get the drawing flags.
	uint32 GetFlags() const;

	/// Append flags to the current flags.
	void AppendFlags(uint32 flags);

	/// Clear flags from the current flags.
	void ClearFlags(uint32 flags);

	/// Draw a closed polygon provided in CCW order.
	virtual void DrawPolygon(const b2Vec2* vertices, int32 vertexCount, const b2Color& color) = 0;

	/// Draw a solid closed polygon provided in CCW order.
	virtual void DrawSolidPolygon(const b2Vec2* vertices, int32 vertexCount, const b2Color& color) = 0;

	/// Draw a circle.
	virtual void DrawCircle(const b2Vec2& center, float32 radius, const b2Color& color) = 0;

	/// Draw a solid circle.
	virtual void DrawSolidCircle(const b2Vec2& center, float32 radius, const b2Vec2& axis, const b2Color& color) = 0;

	/// Draw a particle array
	virtual void DrawParticles(const b2Vec2 *centers, float32 radius, const b2ParticleColor *colors, int32 count) = 0;

	/// Draw a line segment.
	virtual void DrawSegment(const b2Vec2& p1, const b2Vec2& p2, const b2Color& color) = 0;

	/// Draw a transform. Choose your own length scale.
	/// @param xf a transform.
	virtual void DrawTransform(const b2Transform& xf) = 0;

protected:
	uint32 m_drawFlags;
};

// end of Draw.h

b2Draw::b2Draw()
{
	m_drawFlags = 0;
}

void b2Draw::SetFlags(uint32 flags)
{
	m_drawFlags = flags;
}

uint32 b2Draw::GetFlags() const
{
	return m_drawFlags;
}

void b2Draw::AppendFlags(uint32 flags)
{
	m_drawFlags |= flags;
}

void b2Draw::ClearFlags(uint32 flags)
{
	m_drawFlags &= ~flags;
}

// end of Draw.cpp

class b2Shape;

/// A distance proxy is used by the GJK algorithm.
/// It encapsulates any shape.
struct b2DistanceProxy
{
	b2DistanceProxy() : m_vertices(NULL), m_count(0), m_radius(0.0f) {}

	/// Initialize the proxy using the given shape. The shape
	/// must remain in scope while the proxy is in use.
	void Set(const b2Shape* shape, int32 index);

	/// Get the supporting vertex index in the given direction.
	int32 GetSupport(const b2Vec2& d) const;

	/// Get the supporting vertex in the given direction.
	const b2Vec2& GetSupportVertex(const b2Vec2& d) const;

	/// Get the vertex count.
	int32 GetVertexCount() const;

	/// Get a vertex by index. Used by b2Distance.
	const b2Vec2& GetVertex(int32 index) const;

	b2Vec2 m_buffer[2];
	const b2Vec2* m_vertices;
	int32 m_count;
	float32 m_radius;
};

/// Used to warm start b2Distance.
/// Set count to zero on first call.
struct b2SimplexCache
{
	float32 metric;		///< length or area
	uint16 count;
	uint8 indexA[3];	///< vertices on shape A
	uint8 indexB[3];	///< vertices on shape B
};

/// Input for b2Distance.
/// You have to option to use the shape radii
/// in the computation. Even 
struct b2DistanceInput
{
	b2DistanceProxy proxyA;
	b2DistanceProxy proxyB;
	b2Transform transformA;
	b2Transform transformB;
	bool useRadii;
};

/// Output for b2Distance.
struct b2DistanceOutput
{
	b2Vec2 pointA;		///< closest point on shapeA
	b2Vec2 pointB;		///< closest point on shapeB
	float32 distance;
	int32 iterations;	///< number of GJK iterations used
};

/// Compute the closest points between two shapes. Supports any combination of:
/// b2CircleShape, b2PolygonShape, b2EdgeShape. The simplex cache is input/output.
/// On the first call set b2SimplexCache.count to zero.
void b2Distance(b2DistanceOutput* output,
				b2SimplexCache* cache, 
				const b2DistanceInput* input);


//////////////////////////////////////////////////////////////////////////

inline int32 b2DistanceProxy::GetVertexCount() const
{
	return m_count;
}

inline const b2Vec2& b2DistanceProxy::GetVertex(int32 index) const
{
	b2Assert(0 <= index && index < m_count);
	return m_vertices[index];
}

inline int32 b2DistanceProxy::GetSupport(const b2Vec2& d) const
{
	int32 bestIndex = 0;
	float32 bestValue = b2Dot(m_vertices[0], d);
	for (int32 i = 1; i < m_count; ++i)
	{
		float32 value = b2Dot(m_vertices[i], d);
		if (value > bestValue)
		{
			bestIndex = i;
			bestValue = value;
		}
	}

	return bestIndex;
}

inline const b2Vec2& b2DistanceProxy::GetSupportVertex(const b2Vec2& d) const
{
	int32 bestIndex = 0;
	float32 bestValue = b2Dot(m_vertices[0], d);
	for (int32 i = 1; i < m_count; ++i)
	{
		float32 value = b2Dot(m_vertices[i], d);
		if (value > bestValue)
		{
			bestIndex = i;
			bestValue = value;
		}
	}

	return m_vertices[bestIndex];
}

// end of Distance.h

/// @file
/// Structures and functions used for computing contact points, distance
/// queries, and TOI queries.

class b2Shape;
class b2CircleShape;
class b2EdgeShape;
class b2PolygonShape;

const uint8 b2_nullFeature = UCHAR_MAX;

/// The features that intersect to form the contact point
/// This must be 4 bytes or less.
struct b2ContactFeature
{
	enum Type
	{
		e_vertex = 0,
		e_face = 1
	};

	uint8 indexA;		///< Feature index on shapeA
	uint8 indexB;		///< Feature index on shapeB
	uint8 typeA;		///< The feature type on shapeA
	uint8 typeB;		///< The feature type on shapeB
};

/// Contact ids to facilitate warm starting.
union b2ContactID
{
	b2ContactFeature cf;
	uint32 key;					///< Used to quickly compare contact ids.
};

/// A manifold point is a contact point belonging to a contact
/// manifold. It holds details related to the geometry and dynamics
/// of the contact points.
/// The local point usage depends on the manifold type:
/// -e_circles: the local center of circleB
/// -e_faceA: the local center of cirlceB or the clip point of polygonB
/// -e_faceB: the clip point of polygonA
/// This structure is stored across time steps, so we keep it small.
/// Note: the impulses are used for internal caching and may not
/// provide reliable contact forces, especially for high speed collisions.
struct b2ManifoldPoint
{
	b2Vec2 localPoint;		///< usage depends on manifold type
	float32 normalImpulse;	///< the non-penetration impulse
	float32 tangentImpulse;	///< the friction impulse
	b2ContactID id;			///< uniquely identifies a contact point between two shapes
};

/// A manifold for two touching convex shapes.
/// Box2D supports multiple types of contact:
/// - clip point versus plane with radius
/// - point versus point with radius (circles)
/// The local point usage depends on the manifold type:
/// -e_circles: the local center of circleA
/// -e_faceA: the center of faceA
/// -e_faceB: the center of faceB
/// Similarly the local normal usage:
/// -e_circles: not used
/// -e_faceA: the normal on polygonA
/// -e_faceB: the normal on polygonB
/// We store contacts in this way so that position correction can
/// account for movement, which is critical for continuous physics.
/// All contact scenarios must be expressed in one of these types.
/// This structure is stored across time steps, so we keep it small.
struct b2Manifold
{
	enum Type
	{
		e_circles,
		e_faceA,
		e_faceB
	};

	b2ManifoldPoint points[b2_maxManifoldPoints];	///< the points of contact
	b2Vec2 localNormal;								///< not use for Type::e_points
	b2Vec2 localPoint;								///< usage depends on manifold type
	Type type;
	int32 pointCount;								///< the number of manifold points
};

/// This is used to compute the current state of a contact manifold.
struct b2WorldManifold
{
	/// Evaluate the manifold with supplied transforms. This assumes
	/// modest motion from the original state. This does not change the
	/// point count, impulses, etc. The radii must come from the shapes
	/// that generated the manifold.
	void Initialize(const b2Manifold* manifold,
					const b2Transform& xfA, float32 radiusA,
					const b2Transform& xfB, float32 radiusB);

	b2Vec2 normal;								///< world vector pointing from A to B
	b2Vec2 points[b2_maxManifoldPoints];		///< world contact point (point of intersection)
	float32 separations[b2_maxManifoldPoints];	///< a negative value indicates overlap, in meters
};

/// This is used for determining the state of contact points.
enum b2PointState
{
	b2_nullState,		///< point does not exist
	b2_addState,		///< point was added in the update
	b2_persistState,	///< point persisted across the update
	b2_removeState		///< point was removed in the update
};

/// Compute the point states given two manifolds. The states pertain to the transition from manifold1
/// to manifold2. So state1 is either persist or remove while state2 is either add or persist.
void b2GetPointStates(b2PointState state1[b2_maxManifoldPoints], b2PointState state2[b2_maxManifoldPoints],
					  const b2Manifold* manifold1, const b2Manifold* manifold2);

/// Used for computing contact manifolds.
struct b2ClipVertex
{
	b2Vec2 v;
	b2ContactID id;
};

/// Ray-cast input data. The ray extends from p1 to p1 + maxFraction * (p2 - p1).
struct b2RayCastInput
{
	b2Vec2 p1, p2;
	float32 maxFraction;
};

/// Ray-cast output data. The ray hits at p1 + fraction * (p2 - p1), where p1 and p2
/// come from b2RayCastInput.
struct b2RayCastOutput
{
	b2Vec2 normal;
	float32 fraction;
};

/// An axis aligned bounding box.
struct b2AABB
{
	/// Verify that the bounds are sorted.
	bool IsValid() const;

	/// Get the center of the AABB.
	b2Vec2 GetCenter() const
	{
		return 0.5f * (lowerBound + upperBound);
	}

	/// Get the extents of the AABB (half-widths).
	b2Vec2 GetExtents() const
	{
		return 0.5f * (upperBound - lowerBound);
	}

	/// Get the perimeter length
	float32 GetPerimeter() const
	{
		float32 wx = upperBound.x - lowerBound.x;
		float32 wy = upperBound.y - lowerBound.y;
		return 2.0f * (wx + wy);
	}

	/// Combine an AABB into this one.
	void Combine(const b2AABB& aabb)
	{
		lowerBound = b2Min(lowerBound, aabb.lowerBound);
		upperBound = b2Max(upperBound, aabb.upperBound);
	}

	/// Combine two AABBs into this one.
	void Combine(const b2AABB& aabb1, const b2AABB& aabb2)
	{
		lowerBound = b2Min(aabb1.lowerBound, aabb2.lowerBound);
		upperBound = b2Max(aabb1.upperBound, aabb2.upperBound);
	}

	/// Does this aabb contain the provided AABB.
	bool Contains(const b2AABB& aabb) const
	{
		bool result = true;
		result = result && lowerBound.x <= aabb.lowerBound.x;
		result = result && lowerBound.y <= aabb.lowerBound.y;
		result = result && aabb.upperBound.x <= upperBound.x;
		result = result && aabb.upperBound.y <= upperBound.y;
		return result;
	}

	bool RayCast(b2RayCastOutput* output, const b2RayCastInput& input) const;

	b2Vec2 lowerBound;	///< the lower vertex
	b2Vec2 upperBound;	///< the upper vertex
};

/// Compute the collision manifold between two circles.
void b2CollideCircles(b2Manifold* manifold,
					  const b2CircleShape* circleA, const b2Transform& xfA,
					  const b2CircleShape* circleB, const b2Transform& xfB);

/// Compute the collision manifold between a polygon and a circle.
void b2CollidePolygonAndCircle(b2Manifold* manifold,
							   const b2PolygonShape* polygonA, const b2Transform& xfA,
							   const b2CircleShape* circleB, const b2Transform& xfB);

/// Compute the collision manifold between two polygons.
void b2CollidePolygons(b2Manifold* manifold,
					   const b2PolygonShape* polygonA, const b2Transform& xfA,
					   const b2PolygonShape* polygonB, const b2Transform& xfB);

/// Compute the collision manifold between an edge and a circle.
void b2CollideEdgeAndCircle(b2Manifold* manifold,
							   const b2EdgeShape* polygonA, const b2Transform& xfA,
							   const b2CircleShape* circleB, const b2Transform& xfB);

/// Compute the collision manifold between an edge and a circle.
void b2CollideEdgeAndPolygon(b2Manifold* manifold,
							   const b2EdgeShape* edgeA, const b2Transform& xfA,
							   const b2PolygonShape* circleB, const b2Transform& xfB);

/// Clipping for contact manifolds.
int32 b2ClipSegmentToLine(b2ClipVertex vOut[2], const b2ClipVertex vIn[2],
							const b2Vec2& normal, float32 offset, int32 vertexIndexA);

/// Determine if two generic shapes overlap.
bool b2TestOverlap(	const b2Shape* shapeA, int32 indexA,
					const b2Shape* shapeB, int32 indexB,
					const b2Transform& xfA, const b2Transform& xfB);

// ---------------- Inline Functions ------------------------------------------

inline bool b2AABB::IsValid() const
{
	b2Vec2 d = upperBound - lowerBound;
	bool valid = d.x >= 0.0f && d.y >= 0.0f;
	valid = valid && lowerBound.IsValid() && upperBound.IsValid();
	return valid;
}

inline bool b2TestOverlap(const b2AABB& a, const b2AABB& b)
{
	b2Vec2 d1, d2;
	d1 = b.lowerBound - a.upperBound;
	d2 = a.lowerBound - b.upperBound;

	if (d1.x > 0.0f || d1.y > 0.0f)
		return false;

	if (d2.x > 0.0f || d2.y > 0.0f)
		return false;

	return true;
}

// end of Collision.h

/// This holds the mass data computed for a shape.
struct b2MassData
{
	/// The mass of the shape, usually in kilograms.
	float32 mass;

	/// The position of the shape's centroid relative to the shape's origin.
	b2Vec2 center;

	/// The rotational inertia of the shape about the local origin.
	float32 I;
};

/// A shape is used for collision detection. You can create a shape however you like.
/// Shapes used for simulation in b2World are created automatically when a b2Fixture
/// is created. Shapes may encapsulate a one or more child shapes.
class b2Shape
{
public:
	
	enum Type
	{
		e_circle = 0,
		e_edge = 1,
		e_polygon = 2,
		e_chain = 3,
		e_typeCount = 4
	};

	virtual ~b2Shape() {}

	/// Clone the concrete shape using the provided allocator.
	virtual b2Shape* Clone(b2BlockAllocator* allocator) const = 0;

	/// Get the type of this shape. You can use this to down cast to the concrete shape.
	/// @return the shape type.
	Type GetType() const;

	/// Get the number of child primitives.
	virtual int32 GetChildCount() const = 0;

	/// Test a point for containment in this shape. This only works for convex shapes.
	/// @param xf the shape world transform.
	/// @param p a point in world coordinates.
	virtual bool TestPoint(const b2Transform& xf, const b2Vec2& p) const = 0;

	/// Compute the distance from the current shape to the specified point. This only works for convex shapes.
	/// @param xf the shape world transform.
	/// @param p a point in world coordinates.
	/// @param distance returns the distance from the current shape.
	/// @param normal returns the direction in which the distance increases.
	virtual void ComputeDistance(const b2Transform& xf, const b2Vec2& p, float32* distance, b2Vec2* normal, int32 childIndex) const= 0;

	/// Cast a ray against a child shape.
	/// @param output the ray-cast results.
	/// @param input the ray-cast input parameters.
	/// @param transform the transform to be applied to the shape.
	/// @param childIndex the child shape index
	virtual bool RayCast(b2RayCastOutput* output, const b2RayCastInput& input,
						const b2Transform& transform, int32 childIndex) const = 0;

	/// Given a transform, compute the associated axis aligned bounding box for a child shape.
	/// @param aabb returns the axis aligned box.
	/// @param xf the world transform of the shape.
	/// @param childIndex the child shape
	virtual void ComputeAABB(b2AABB* aabb, const b2Transform& xf, int32 childIndex) const = 0;

	/// Compute the mass properties of this shape using its dimensions and density.
	/// The inertia tensor is computed about the local origin.
	/// @param massData returns the mass data for this shape.
	/// @param density the density in kilograms per meter squared.
	virtual void ComputeMass(b2MassData* massData, float32 density) const = 0;

	Type m_type;
	float32 m_radius;
};

inline b2Shape::Type b2Shape::GetType() const
{
	return m_type;
}

// end of Shape.h

/// A circle shape.
class b2CircleShape : public b2Shape
{
public:
	b2CircleShape();

	/// Implement b2Shape.
	b2Shape* Clone(b2BlockAllocator* allocator) const;

	/// @see b2Shape::GetChildCount
	int32 GetChildCount() const;

	/// Implement b2Shape.
	bool TestPoint(const b2Transform& transform, const b2Vec2& p) const;

	// @see b2Shape::ComputeDistance
	void ComputeDistance(const b2Transform& xf, const b2Vec2& p, float32* distance, b2Vec2* normal, int32 childIndex) const;

	/// Implement b2Shape.
	bool RayCast(b2RayCastOutput* output, const b2RayCastInput& input,
				const b2Transform& transform, int32 childIndex) const;

	/// @see b2Shape::ComputeAABB
	void ComputeAABB(b2AABB* aabb, const b2Transform& transform, int32 childIndex) const;

	/// @see b2Shape::ComputeMass
	void ComputeMass(b2MassData* massData, float32 density) const;

	/// Get the supporting vertex index in the given direction.
	int32 GetSupport(const b2Vec2& d) const;

	/// Get the supporting vertex in the given direction.
	const b2Vec2& GetSupportVertex(const b2Vec2& d) const;

	/// Get the vertex count.
	int32 GetVertexCount() const { return 1; }

	/// Get a vertex by index. Used by b2Distance.
	const b2Vec2& GetVertex(int32 index) const;

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
public:
	/// Set position with direct floats.
	void SetPosition(float32 x, float32 y) { m_p.Set(x, y); }

	/// Get x-coordinate of position.
	float32 GetPositionX() const { return m_p.x; }

	/// Get y-coordinate of position.
	float32 GetPositionY() const { return m_p.y; }
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

	/// Position
	b2Vec2 m_p;
};

inline b2CircleShape::b2CircleShape()
{
	m_type = e_circle;
	m_radius = 0.0f;
	m_p.SetZero();
}

inline int32 b2CircleShape::GetSupport(const b2Vec2 &d) const
{
	B2_NOT_USED(d);
	return 0;
}

inline const b2Vec2& b2CircleShape::GetSupportVertex(const b2Vec2 &d) const
{
	B2_NOT_USED(d);
	return m_p;
}

inline const b2Vec2& b2CircleShape::GetVertex(int32 index) const
{
	B2_NOT_USED(index);
	b2Assert(index == 0);
	return m_p;
}

// end of CircleShape.h

b2Shape* b2CircleShape::Clone(b2BlockAllocator* allocator) const
{
	void* mem = allocator->Allocate(sizeof(b2CircleShape));
	b2CircleShape* clone = new (mem) b2CircleShape;
	*clone = *this;
	return clone;
}

int32 b2CircleShape::GetChildCount() const
{
	return 1;
}

bool b2CircleShape::TestPoint(const b2Transform& transform, const b2Vec2& p) const
{
	b2Vec2 center = transform.p + b2Mul(transform.q, m_p);
	b2Vec2 d = p - center;
	return b2Dot(d, d) <= m_radius * m_radius;
}

void b2CircleShape::ComputeDistance(const b2Transform& transform, const b2Vec2& p, float32* distance, b2Vec2* normal, int32 childIndex) const
{
	B2_NOT_USED(childIndex);

	b2Vec2 center = transform.p + b2Mul(transform.q, m_p);
	b2Vec2 d = p - center;
	float32 d1 = d.Length();
	*distance = d1 - m_radius;
	*normal = 1 / d1 * d;
}

// Collision Detection in Interactive 3D Environments by Gino van den Bergen
// From Section 3.1.2
// x = s + a * r
// norm(x) = radius
bool b2CircleShape::RayCast(b2RayCastOutput* output, const b2RayCastInput& input,
							const b2Transform& transform, int32 childIndex) const
{
	B2_NOT_USED(childIndex);

	b2Vec2 position = transform.p + b2Mul(transform.q, m_p);
	b2Vec2 s = input.p1 - position;
	float32 b = b2Dot(s, s) - m_radius * m_radius;

	// Solve quadratic equation.
	b2Vec2 r = input.p2 - input.p1;
	float32 c =  b2Dot(s, r);
	float32 rr = b2Dot(r, r);
	float32 sigma = c * c - rr * b;

	// Check for negative discriminant and short segment.
	if (sigma < 0.0f || rr < b2_epsilon)
	{
		return false;
	}

	// Find the point of intersection of the line with the circle.
	float32 a = -(c + b2Sqrt(sigma));

	// Is the intersection point on the segment?
	if (0.0f <= a && a <= input.maxFraction * rr)
	{
		a /= rr;
		output->fraction = a;
		output->normal = s + a * r;
		output->normal.Normalize();
		return true;
	}

	return false;
}

void b2CircleShape::ComputeAABB(b2AABB* aabb, const b2Transform& transform, int32 childIndex) const
{
	B2_NOT_USED(childIndex);

	b2Vec2 p = transform.p + b2Mul(transform.q, m_p);
	aabb->lowerBound.Set(p.x - m_radius, p.y - m_radius);
	aabb->upperBound.Set(p.x + m_radius, p.y + m_radius);
}

void b2CircleShape::ComputeMass(b2MassData* massData, float32 density) const
{
	massData->mass = density * b2_pi * m_radius * m_radius;
	massData->center = m_p;

	// inertia about the local origin
	massData->I = massData->mass * (0.5f * m_radius * m_radius + b2Dot(m_p, m_p));
}

// end of CircleShape.cpp

/// A line segment (edge) shape. These can be connected in chains or loops
/// to other edge shapes. The connectivity information is used to ensure
/// correct contact normals.
class b2EdgeShape : public b2Shape
{
public:
	b2EdgeShape();

	/// Set this as an isolated edge.
	void Set(const b2Vec2& v1, const b2Vec2& v2);

	/// Implement b2Shape.
	b2Shape* Clone(b2BlockAllocator* allocator) const;

	/// @see b2Shape::GetChildCount
	int32 GetChildCount() const;

	/// @see b2Shape::TestPoint
	bool TestPoint(const b2Transform& transform, const b2Vec2& p) const;

	// @see b2Shape::ComputeDistance
	void ComputeDistance(const b2Transform& xf, const b2Vec2& p, float32* distance, b2Vec2* normal, int32 childIndex) const;

	/// Implement b2Shape.
	bool RayCast(b2RayCastOutput* output, const b2RayCastInput& input,
				const b2Transform& transform, int32 childIndex) const;

	/// @see b2Shape::ComputeAABB
	void ComputeAABB(b2AABB* aabb, const b2Transform& transform, int32 childIndex) const;

	/// @see b2Shape::ComputeMass
	void ComputeMass(b2MassData* massData, float32 density) const;

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
public:
	/// Set this as an isolated edge, with direct floats.
	void Set(float32 vx1, float32 vy1, float32 vx2, float32 vy2);
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

	/// These are the edge vertices
	b2Vec2 m_vertex1, m_vertex2;

	/// Optional adjacent vertices. These are used for smooth collision.
	b2Vec2 m_vertex0, m_vertex3;
	bool m_hasVertex0, m_hasVertex3;
};

inline b2EdgeShape::b2EdgeShape()
{
	m_type = e_edge;
	m_radius = b2_polygonRadius;
	m_vertex0.x = 0.0f;
	m_vertex0.y = 0.0f;
	m_vertex3.x = 0.0f;
	m_vertex3.y = 0.0f;
	m_hasVertex0 = false;
	m_hasVertex3 = false;
}

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
inline void b2EdgeShape::Set(float32 vx1,
														 float32 vy1,
														 float32 vx2,
														 float32 vy2) {
	Set(b2Vec2(vx1, vy1), b2Vec2(vx2, vy2));
}
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

// end of EdgeShape.h

void b2EdgeShape::Set(const b2Vec2& v1, const b2Vec2& v2)
{
	m_vertex1 = v1;
	m_vertex2 = v2;
	m_hasVertex0 = false;
	m_hasVertex3 = false;
}

b2Shape* b2EdgeShape::Clone(b2BlockAllocator* allocator) const
{
	void* mem = allocator->Allocate(sizeof(b2EdgeShape));
	b2EdgeShape* clone = new (mem) b2EdgeShape;
	*clone = *this;
	return clone;
}

int32 b2EdgeShape::GetChildCount() const
{
	return 1;
}

bool b2EdgeShape::TestPoint(const b2Transform& xf, const b2Vec2& p) const
{
	B2_NOT_USED(xf);
	B2_NOT_USED(p);
	return false;
}

void b2EdgeShape::ComputeDistance(const b2Transform& xf, const b2Vec2& p, float32* distance, b2Vec2* normal, int32 childIndex) const
{
	B2_NOT_USED(childIndex);

	b2Vec2 v1 = b2Mul(xf, m_vertex1);
	b2Vec2 v2 = b2Mul(xf, m_vertex2);

	b2Vec2 d = p - v1;
	b2Vec2 s = v2 - v1;
	float32 ds = b2Dot(d, s);
	if (ds > 0)
	{
		float32 s2 = b2Dot(s, s);
		if (ds > s2)
		{
			d = p - v2;
		}
		else
		{
			d -= ds / s2 * s;
		}
	}

	float32 d1 = d.Length();
	*distance = d1;
	*normal = d1 > 0 ? 1 / d1 * d : b2Vec2_zero;

}

// p = p1 + t * d
// v = v1 + s * e
// p1 + t * d = v1 + s * e
// s * e - t * d = p1 - v1
bool b2EdgeShape::RayCast(b2RayCastOutput* output, const b2RayCastInput& input,
							const b2Transform& xf, int32 childIndex) const
{
	B2_NOT_USED(childIndex);

	// Put the ray into the edge's frame of reference.
	b2Vec2 p1 = b2MulT(xf.q, input.p1 - xf.p);
	b2Vec2 p2 = b2MulT(xf.q, input.p2 - xf.p);
	b2Vec2 d = p2 - p1;

	b2Vec2 v1 = m_vertex1;
	b2Vec2 v2 = m_vertex2;
	b2Vec2 e = v2 - v1;
	b2Vec2 normal(e.y, -e.x);
	normal.Normalize();

	// q = p1 + t * d
	// dot(normal, q - v1) = 0
	// dot(normal, p1 - v1) + t * dot(normal, d) = 0
	float32 numerator = b2Dot(normal, v1 - p1);
	float32 denominator = b2Dot(normal, d);

	if (denominator == 0.0f)
	{
		return false;
	}

	float32 t = numerator / denominator;
	if (t < 0.0f || input.maxFraction < t)
	{
		return false;
	}

	b2Vec2 q = p1 + t * d;

	// q = v1 + s * r
	// s = dot(q - v1, r) / dot(r, r)
	b2Vec2 r = v2 - v1;
	float32 rr = b2Dot(r, r);
	if (rr == 0.0f)
	{
		return false;
	}

	float32 s = b2Dot(q - v1, r) / rr;
	if (s < 0.0f || 1.0f < s)
	{
		return false;
	}

	output->fraction = t;
	if (numerator > 0.0f)
	{
		output->normal = -b2Mul(xf.q, normal);
	}
	else
	{
		output->normal = b2Mul(xf.q, normal);
	}
	return true;
}

void b2EdgeShape::ComputeAABB(b2AABB* aabb, const b2Transform& xf, int32 childIndex) const
{
	B2_NOT_USED(childIndex);

	b2Vec2 v1 = b2Mul(xf, m_vertex1);
	b2Vec2 v2 = b2Mul(xf, m_vertex2);

	b2Vec2 lower = b2Min(v1, v2);
	b2Vec2 upper = b2Max(v1, v2);

	b2Vec2 r(m_radius, m_radius);
	aabb->lowerBound = lower - r;
	aabb->upperBound = upper + r;
}

void b2EdgeShape::ComputeMass(b2MassData* massData, float32 density) const
{
	B2_NOT_USED(density);

	massData->mass = 0.0f;
	massData->center = 0.5f * (m_vertex1 + m_vertex2);
	massData->I = 0.0f;
}

// end of EdgeShape.cpp

class b2EdgeShape;

/// A chain shape is a free form sequence of line segments.
/// The chain has two-sided collision, so you can use inside and outside collision.
/// Therefore, you may use any winding order.
/// Since there may be many vertices, they are allocated using b2Alloc.
/// Connectivity information is used to create smooth collisions.
/// WARNING: The chain will not collide properly if there are self-intersections.
class b2ChainShape : public b2Shape
{
public:
	b2ChainShape();

	/// The destructor frees the vertices using b2Free.
	~b2ChainShape();

	/// Create a loop. This automatically adjusts connectivity.
	/// @param vertices an array of vertices, these are copied
	/// @param count the vertex count
	void CreateLoop(const b2Vec2* vertices, int32 count);

	/// Create a chain with isolated end vertices.
	/// @param vertices an array of vertices, these are copied
	/// @param count the vertex count
	void CreateChain(const b2Vec2* vertices, int32 count);

	/// Establish connectivity to a vertex that precedes the first vertex.
	/// Don't call this for loops.
	void SetPrevVertex(const b2Vec2& prevVertex);

	/// Establish connectivity to a vertex that follows the last vertex.
	/// Don't call this for loops.
	void SetNextVertex(const b2Vec2& nextVertex);

	/// Implement b2Shape. Vertices are cloned using b2Alloc.
	b2Shape* Clone(b2BlockAllocator* allocator) const;

	/// @see b2Shape::GetChildCount
	int32 GetChildCount() const;

	/// Get a child edge.
	void GetChildEdge(b2EdgeShape* edge, int32 index) const;

	/// This always return false.
	/// @see b2Shape::TestPoint
	bool TestPoint(const b2Transform& transform, const b2Vec2& p) const;

	// @see b2Shape::ComputeDistance
	void ComputeDistance(const b2Transform& xf, const b2Vec2& p, float32* distance, b2Vec2* normal, int32 childIndex) const;

	/// Implement b2Shape.
	bool RayCast(b2RayCastOutput* output, const b2RayCastInput& input,
					const b2Transform& transform, int32 childIndex) const;

	/// @see b2Shape::ComputeAABB
	void ComputeAABB(b2AABB* aabb, const b2Transform& transform, int32 childIndex) const;

	/// Chains have zero mass.
	/// @see b2Shape::ComputeMass
	void ComputeMass(b2MassData* massData, float32 density) const;

	/// The vertices. Owned by this class.
	b2Vec2* m_vertices;

	/// The vertex count.
	int32 m_count;

	b2Vec2 m_prevVertex, m_nextVertex;
	bool m_hasPrevVertex, m_hasNextVertex;
};

inline b2ChainShape::b2ChainShape()
{
	m_type = e_chain;
	m_radius = b2_polygonRadius;
	m_vertices = NULL;
	m_count = 0;
	m_hasPrevVertex = false;
	m_hasNextVertex = false;
}

// end of ChainShape.h

b2ChainShape::~b2ChainShape()
{
	b2Free(m_vertices);
	m_vertices = NULL;
	m_count = 0;
}

void b2ChainShape::CreateLoop(const b2Vec2* vertices, int32 count)
{
	b2Assert(m_vertices == NULL && m_count == 0);
	b2Assert(count >= 3);
	for (int32 i = 1; i < count; ++i)
	{
#if B2_ASSERT_ENABLED
		b2Vec2 v1 = vertices[i-1];
		b2Vec2 v2 = vertices[i];
		// If the code crashes here, it means your vertices are too close together.
		b2Assert(b2DistanceSquared(v1, v2) > b2_linearSlop * b2_linearSlop);
#endif // B2_ASSERT_ENABLED
	}

	m_count = count + 1;
	m_vertices = (b2Vec2*)b2Alloc(m_count * sizeof(b2Vec2));
	memcpy(m_vertices, vertices, count * sizeof(b2Vec2));
	m_vertices[count] = m_vertices[0];
	m_prevVertex = m_vertices[m_count - 2];
	m_nextVertex = m_vertices[1];
	m_hasPrevVertex = true;
	m_hasNextVertex = true;
}

void b2ChainShape::CreateChain(const b2Vec2* vertices, int32 count)
{
	b2Assert(m_vertices == NULL && m_count == 0);
	b2Assert(count >= 2);
	for (int32 i = 1; i < count; ++i)
	{
#if B2_ASSERT_ENABLED
		b2Vec2 v1 = vertices[i-1];
		b2Vec2 v2 = vertices[i];
		// If the code crashes here, it means your vertices are too close together.
		b2Assert(b2DistanceSquared(v1, v2) > b2_linearSlop * b2_linearSlop);
#endif // B2_ASSERT_ENABLED
	}

	m_count = count;
	m_vertices = (b2Vec2*)b2Alloc(count * sizeof(b2Vec2));
	memcpy(m_vertices, vertices, m_count * sizeof(b2Vec2));

	m_hasPrevVertex = false;
	m_hasNextVertex = false;

	m_prevVertex.SetZero();
	m_nextVertex.SetZero();
}

void b2ChainShape::SetPrevVertex(const b2Vec2& prevVertex)
{
	m_prevVertex = prevVertex;
	m_hasPrevVertex = true;
}

void b2ChainShape::SetNextVertex(const b2Vec2& nextVertex)
{
	m_nextVertex = nextVertex;
	m_hasNextVertex = true;
}

b2Shape* b2ChainShape::Clone(b2BlockAllocator* allocator) const
{
	void* mem = allocator->Allocate(sizeof(b2ChainShape));
	b2ChainShape* clone = new (mem) b2ChainShape;
	clone->CreateChain(m_vertices, m_count);
	clone->m_prevVertex = m_prevVertex;
	clone->m_nextVertex = m_nextVertex;
	clone->m_hasPrevVertex = m_hasPrevVertex;
	clone->m_hasNextVertex = m_hasNextVertex;
	return clone;
}

int32 b2ChainShape::GetChildCount() const
{
	// edge count = vertex count - 1
	return m_count - 1;
}

void b2ChainShape::GetChildEdge(b2EdgeShape* edge, int32 index) const
{
	b2Assert(0 <= index && index < m_count - 1);
	edge->m_type = b2Shape::e_edge;
	edge->m_radius = m_radius;

	edge->m_vertex1 = m_vertices[index + 0];
	edge->m_vertex2 = m_vertices[index + 1];

	if (index > 0)
	{
		edge->m_vertex0 = m_vertices[index - 1];
		edge->m_hasVertex0 = true;
	}
	else
	{
		edge->m_vertex0 = m_prevVertex;
		edge->m_hasVertex0 = m_hasPrevVertex;
	}

	if (index < m_count - 2)
	{
		edge->m_vertex3 = m_vertices[index + 2];
		edge->m_hasVertex3 = true;
	}
	else
	{
		edge->m_vertex3 = m_nextVertex;
		edge->m_hasVertex3 = m_hasNextVertex;
	}
}

void b2ChainShape::ComputeDistance(const b2Transform& xf, const b2Vec2& p, float32* distance, b2Vec2* normal, int32 childIndex) const
{
	b2EdgeShape edge;
	GetChildEdge(&edge, childIndex);
	edge.ComputeDistance(xf, p, distance, normal, 0);
}

bool b2ChainShape::TestPoint(const b2Transform& xf, const b2Vec2& p) const
{
	B2_NOT_USED(xf);
	B2_NOT_USED(p);
	return false;
}

bool b2ChainShape::RayCast(b2RayCastOutput* output, const b2RayCastInput& input,
							const b2Transform& xf, int32 childIndex) const
{
	b2Assert(childIndex < m_count);

	b2EdgeShape edgeShape;

	int32 i1 = childIndex;
	int32 i2 = childIndex + 1;
	if (i2 == m_count)
	{
		i2 = 0;
	}

	edgeShape.m_vertex1 = m_vertices[i1];
	edgeShape.m_vertex2 = m_vertices[i2];

	return edgeShape.RayCast(output, input, xf, 0);
}

void b2ChainShape::ComputeAABB(b2AABB* aabb, const b2Transform& xf, int32 childIndex) const
{
	b2Assert(childIndex < m_count);

	int32 i1 = childIndex;
	int32 i2 = childIndex + 1;
	if (i2 == m_count)
	{
		i2 = 0;
	}

	b2Vec2 v1 = b2Mul(xf, m_vertices[i1]);
	b2Vec2 v2 = b2Mul(xf, m_vertices[i2]);

	aabb->lowerBound = b2Min(v1, v2);
	aabb->upperBound = b2Max(v1, v2);
}

void b2ChainShape::ComputeMass(b2MassData* massData, float32 density) const
{
	B2_NOT_USED(density);

	massData->mass = 0.0f;
	massData->center.SetZero();
	massData->I = 0.0f;
}

// end of ChainShape.cpp

/// A convex polygon. It is assumed that the interior of the polygon is to
/// the left of each edge.
/// Polygons have a maximum number of vertices equal to b2_maxPolygonVertices.
/// In most cases you should not need many vertices for a convex polygon.
class b2PolygonShape : public b2Shape
{
public:
	b2PolygonShape();

	/// Implement b2Shape.
	b2Shape* Clone(b2BlockAllocator* allocator) const;

	/// @see b2Shape::GetChildCount
	int32 GetChildCount() const;

	/// Create a convex hull from the given array of local points.
	/// The count must be in the range [3, b2_maxPolygonVertices].
	/// @warning the points may be re-ordered, even if they form a convex polygon
	/// @warning collinear points are handled but not removed. Collinear points
	/// may lead to poor stacking behavior.
	void Set(const b2Vec2* points, int32 count);

	/// Build vertices to represent an axis-aligned box centered on the local origin.
	/// @param hx the half-width.
	/// @param hy the half-height.
	void SetAsBox(float32 hx, float32 hy);

	/// Build vertices to represent an oriented box.
	/// @param hx the half-width.
	/// @param hy the half-height.
	/// @param center the center of the box in local coordinates.
	/// @param angle the rotation of the box in local coordinates.
	void SetAsBox(float32 hx, float32 hy, const b2Vec2& center, float32 angle);

	/// @see b2Shape::TestPoint
	bool TestPoint(const b2Transform& transform, const b2Vec2& p) const;

	// @see b2Shape::ComputeDistance
	void ComputeDistance(const b2Transform& xf, const b2Vec2& p, float32* distance, b2Vec2* normal, int32 childIndex) const;

	/// Implement b2Shape.
	bool RayCast(b2RayCastOutput* output, const b2RayCastInput& input,
					const b2Transform& transform, int32 childIndex) const;

	/// @see b2Shape::ComputeAABB
	void ComputeAABB(b2AABB* aabb, const b2Transform& transform, int32 childIndex) const;

	/// @see b2Shape::ComputeMass
	void ComputeMass(b2MassData* massData, float32 density) const;

	/// Get the vertex count.
	int32 GetVertexCount() const { return m_count; }

	/// Get a vertex by index.
	const b2Vec2& GetVertex(int32 index) const;

	/// Validate convexity. This is a very time consuming operation.
	/// @returns true if valid
	bool Validate() const;

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
public:
	/// Set centroid with direct floats.
	void SetCentroid(float32 x, float32 y);

	/// SetAsBox with direct floats for center.
	/// @see b2Shape::SetAsBox
	void SetAsBox(float32 hx,
								float32 hy,
								float32 centerX,
								float32 centerY,
								float32 angle);
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

	b2Vec2 m_centroid;
	b2Vec2 m_vertices[b2_maxPolygonVertices];
	b2Vec2 m_normals[b2_maxPolygonVertices];
	int32 m_count;
};

inline b2PolygonShape::b2PolygonShape()
{
	m_type = e_polygon;
	m_radius = b2_polygonRadius;
	m_count = 0;
	m_centroid.SetZero();
}

inline const b2Vec2& b2PolygonShape::GetVertex(int32 index) const
{
	b2Assert(0 <= index && index < m_count);
	return m_vertices[index];
}

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
inline void b2PolygonShape::SetCentroid(float32 x, float32 y)
{
	m_centroid.Set(x, y);
}

inline void b2PolygonShape::SetAsBox(float32 hx,
										 float32 hy,
										 float32 centerX,
										 float32 centerY,
										 float32 angle) {
	SetAsBox(hx, hy, b2Vec2(centerX, centerY), angle);
}
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

// end of PolygonShape.h

b2Shape* b2PolygonShape::Clone(b2BlockAllocator* allocator) const
{
	void* mem = allocator->Allocate(sizeof(b2PolygonShape));
	b2PolygonShape* clone = new (mem) b2PolygonShape;
	*clone = *this;
	return clone;
}

void b2PolygonShape::SetAsBox(float32 hx, float32 hy)
{
	m_count = 4;
	m_vertices[0].Set(-hx, -hy);
	m_vertices[1].Set( hx, -hy);
	m_vertices[2].Set( hx,  hy);
	m_vertices[3].Set(-hx,  hy);
	m_normals[0].Set(0.0f, -1.0f);
	m_normals[1].Set(1.0f, 0.0f);
	m_normals[2].Set(0.0f, 1.0f);
	m_normals[3].Set(-1.0f, 0.0f);
	m_centroid.SetZero();
}

void b2PolygonShape::SetAsBox(float32 hx, float32 hy, const b2Vec2& center, float32 angle)
{
	m_count = 4;
	m_vertices[0].Set(-hx, -hy);
	m_vertices[1].Set( hx, -hy);
	m_vertices[2].Set( hx,  hy);
	m_vertices[3].Set(-hx,  hy);
	m_normals[0].Set(0.0f, -1.0f);
	m_normals[1].Set(1.0f, 0.0f);
	m_normals[2].Set(0.0f, 1.0f);
	m_normals[3].Set(-1.0f, 0.0f);
	m_centroid = center;

	b2Transform xf;
	xf.p = center;
	xf.q.Set(angle);

	// Transform vertices and normals.
	for (int32 i = 0; i < m_count; ++i)
	{
		m_vertices[i] = b2Mul(xf, m_vertices[i]);
		m_normals[i] = b2Mul(xf.q, m_normals[i]);
	}
}

int32 b2PolygonShape::GetChildCount() const
{
	return 1;
}

static b2Vec2 ComputeCentroid(const b2Vec2* vs, int32 count)
{
	b2Assert(count >= 3);

	b2Vec2 c; c.Set(0.0f, 0.0f);
	float32 area = 0.0f;

	// pRef is the reference point for forming triangles.
	// It's location doesn't change the result (except for rounding error).
	b2Vec2 pRef(0.0f, 0.0f);
#if 0
	// This code would put the reference point inside the polygon.
	for (int32 i = 0; i < count; ++i)
	{
		pRef += vs[i];
	}
	pRef *= 1.0f / count;
#endif

	const float32 inv3 = 1.0f / 3.0f;

	for (int32 i = 0; i < count; ++i)
	{
		// Triangle vertices.
		b2Vec2 p1 = pRef;
		b2Vec2 p2 = vs[i];
		b2Vec2 p3 = i + 1 < count ? vs[i+1] : vs[0];

		b2Vec2 e1 = p2 - p1;
		b2Vec2 e2 = p3 - p1;

		float32 D = b2Cross(e1, e2);

		float32 triangleArea = 0.5f * D;
		area += triangleArea;

		// Area weighted centroid
		c += triangleArea * inv3 * (p1 + p2 + p3);
	}

	// Centroid
	b2Assert(area > b2_epsilon);
	c *= 1.0f / area;
	return c;
}

void b2PolygonShape::Set(const b2Vec2* vertices, int32 count)
{
	b2Assert(3 <= count && count <= b2_maxPolygonVertices);
	if (count < 3)
	{
		SetAsBox(1.0f, 1.0f);
		return;
	}
	
	int32 n = b2Min(count, b2_maxPolygonVertices);

	// Perform welding and copy vertices into local buffer.
	b2Vec2 ps[b2_maxPolygonVertices];
	int32 tempCount = 0;
	for (int32 i = 0; i < n; ++i)
	{
		b2Vec2 v = vertices[i];

		bool unique = true;
		for (int32 j = 0; j < tempCount; ++j)
		{
			if (b2DistanceSquared(v, ps[j]) < 0.5f * b2_linearSlop)
			{
				unique = false;
				break;
			}
		}

		if (unique)
		{
			ps[tempCount++] = v;
		}
	}

	n = tempCount;
	if (n < 3)
	{
		// Polygon is degenerate.
		b2Assert(false);
		SetAsBox(1.0f, 1.0f);
		return;
	}

	// Create the convex hull using the Gift wrapping algorithm
	// http://en.wikipedia.org/wiki/Gift_wrapping_algorithm

	// Find the right most point on the hull
	int32 i0 = 0;
	float32 x0 = ps[0].x;
	for (int32 i = 1; i < n; ++i)
	{
		float32 x = ps[i].x;
		if (x > x0 || (x == x0 && ps[i].y < ps[i0].y))
		{
			i0 = i;
			x0 = x;
		}
	}

	int32 hull[b2_maxPolygonVertices];
	int32 m = 0;
	int32 ih = i0;

	for (;;)
	{
		hull[m] = ih;

		int32 ie = 0;
		for (int32 j = 1; j < n; ++j)
		{
			if (ie == ih)
			{
				ie = j;
				continue;
			}

			b2Vec2 r = ps[ie] - ps[hull[m]];
			b2Vec2 v = ps[j] - ps[hull[m]];
			float32 c = b2Cross(r, v);
			if (c < 0.0f)
			{
				ie = j;
			}

			// Collinearity check
			if (c == 0.0f && v.LengthSquared() > r.LengthSquared())
			{
				ie = j;
			}
		}

		++m;
		ih = ie;

		if (ie == i0)
		{
			break;
		}
	}
	
	m_count = m;

	// Copy vertices.
	for (int32 i = 0; i < m; ++i)
	{
		m_vertices[i] = ps[hull[i]];
	}

	// Compute normals. Ensure the edges have non-zero length.
	for (int32 i = 0; i < m; ++i)
	{
		int32 i1 = i;
		int32 i2 = i + 1 < m ? i + 1 : 0;
		b2Vec2 edge = m_vertices[i2] - m_vertices[i1];
		b2Assert(edge.LengthSquared() > b2_epsilon * b2_epsilon);
		m_normals[i] = b2Cross(edge, 1.0f);
		m_normals[i].Normalize();
	}

	// Compute the polygon centroid.
	m_centroid = ComputeCentroid(m_vertices, m);
}

bool b2PolygonShape::TestPoint(const b2Transform& xf, const b2Vec2& p) const
{
	b2Vec2 pLocal = b2MulT(xf.q, p - xf.p);

	for (int32 i = 0; i < m_count; ++i)
	{
		float32 dot = b2Dot(m_normals[i], pLocal - m_vertices[i]);
		if (dot > 0.0f)
		{
			return false;
		}
	}

	return true;
}

void b2PolygonShape::ComputeDistance(const b2Transform& xf, const b2Vec2& p, float32* distance, b2Vec2* normal, int32 childIndex) const
{
	B2_NOT_USED(childIndex);

	b2Vec2 pLocal = b2MulT(xf.q, p - xf.p);
	float32 maxDistance = -FLT_MAX;
	b2Vec2 normalForMaxDistance = pLocal;

	for (int32 i = 0; i < m_count; ++i)
	{
		float32 dot = b2Dot(m_normals[i], pLocal - m_vertices[i]);
		if (dot > maxDistance)
		{
			maxDistance = dot;
			normalForMaxDistance = m_normals[i];
		}
	}

	if (maxDistance > 0)
	{
		b2Vec2 minDistance = normalForMaxDistance;
		float32 minDistance2 = maxDistance * maxDistance;
		for (int32 i = 0; i < m_count; ++i)
		{
			b2Vec2 distance = pLocal - m_vertices[i];
			float32 distance2 = distance.LengthSquared();
			if (minDistance2 > distance2)
			{
				minDistance = distance;
				minDistance2 = distance2;
			}
		}

		*distance = b2Sqrt(minDistance2);
		*normal = b2Mul(xf.q, minDistance);
		normal->Normalize();
	}
	else
	{
		*distance = maxDistance;
		*normal = b2Mul(xf.q, normalForMaxDistance);
	}
}

bool b2PolygonShape::RayCast(b2RayCastOutput* output, const b2RayCastInput& input,
								const b2Transform& xf, int32 childIndex) const
{
	B2_NOT_USED(childIndex);

	// Put the ray into the polygon's frame of reference.
	b2Vec2 p1 = b2MulT(xf.q, input.p1 - xf.p);
	b2Vec2 p2 = b2MulT(xf.q, input.p2 - xf.p);
	b2Vec2 d = p2 - p1;

	float32 lower = 0.0f, upper = input.maxFraction;

	int32 index = -1;

	for (int32 i = 0; i < m_count; ++i)
	{
		// p = p1 + a * d
		// dot(normal, p - v) = 0
		// dot(normal, p1 - v) + a * dot(normal, d) = 0
		float32 numerator = b2Dot(m_normals[i], m_vertices[i] - p1);
		float32 denominator = b2Dot(m_normals[i], d);

		if (denominator == 0.0f)
		{	
			if (numerator < 0.0f)
			{
				return false;
			}
		}
		else
		{
			// Note: we want this predicate without division:
			// lower < numerator / denominator, where denominator < 0
			// Since denominator < 0, we have to flip the inequality:
			// lower < numerator / denominator <==> denominator * lower > numerator.
			if (denominator < 0.0f && numerator < lower * denominator)
			{
				// Increase lower.
				// The segment enters this half-space.
				lower = numerator / denominator;
				index = i;
			}
			else if (denominator > 0.0f && numerator < upper * denominator)
			{
				// Decrease upper.
				// The segment exits this half-space.
				upper = numerator / denominator;
			}
		}

		// The use of epsilon here causes the assert on lower to trip
		// in some cases. Apparently the use of epsilon was to make edge
		// shapes work, but now those are handled separately.
		//if (upper < lower - b2_epsilon)
		if (upper < lower)
		{
			return false;
		}
	}

	b2Assert(0.0f <= lower && lower <= input.maxFraction);

	if (index >= 0)
	{
		output->fraction = lower;
		output->normal = b2Mul(xf.q, m_normals[index]);
		return true;
	}

	return false;
}

void b2PolygonShape::ComputeAABB(b2AABB* aabb, const b2Transform& xf, int32 childIndex) const
{
	B2_NOT_USED(childIndex);

	b2Vec2 lower = b2Mul(xf, m_vertices[0]);
	b2Vec2 upper = lower;

	for (int32 i = 1; i < m_count; ++i)
	{
		b2Vec2 v = b2Mul(xf, m_vertices[i]);
		lower = b2Min(lower, v);
		upper = b2Max(upper, v);
	}

	b2Vec2 r(m_radius, m_radius);
	aabb->lowerBound = lower - r;
	aabb->upperBound = upper + r;
}

void b2PolygonShape::ComputeMass(b2MassData* massData, float32 density) const
{
	// Polygon mass, centroid, and inertia.
	// Let rho be the polygon density in mass per unit area.
	// Then:
	// mass = rho * int(dA)
	// centroid.x = (1/mass) * rho * int(x * dA)
	// centroid.y = (1/mass) * rho * int(y * dA)
	// I = rho * int((x*x + y*y) * dA)
	//
	// We can compute these integrals by summing all the integrals
	// for each triangle of the polygon. To evaluate the integral
	// for a single triangle, we make a change of variables to
	// the (u,v) coordinates of the triangle:
	// x = x0 + e1x * u + e2x * v
	// y = y0 + e1y * u + e2y * v
	// where 0 <= u && 0 <= v && u + v <= 1.
	//
	// We integrate u from [0,1-v] and then v from [0,1].
	// We also need to use the Jacobian of the transformation:
	// D = cross(e1, e2)
	//
	// Simplification: triangle centroid = (1/3) * (p1 + p2 + p3)
	//
	// The rest of the derivation is handled by computer algebra.

	b2Assert(m_count >= 3);

	b2Vec2 center; center.Set(0.0f, 0.0f);
	float32 area = 0.0f;
	float32 I = 0.0f;

	// s is the reference point for forming triangles.
	// It's location doesn't change the result (except for rounding error).
	b2Vec2 s(0.0f, 0.0f);

	// This code would put the reference point inside the polygon.
	for (int32 i = 0; i < m_count; ++i)
	{
		s += m_vertices[i];
	}
	s *= 1.0f / m_count;

	const float32 k_inv3 = 1.0f / 3.0f;

	for (int32 i = 0; i < m_count; ++i)
	{
		// Triangle vertices.
		b2Vec2 e1 = m_vertices[i] - s;
		b2Vec2 e2 = i + 1 < m_count ? m_vertices[i+1] - s : m_vertices[0] - s;

		float32 D = b2Cross(e1, e2);

		float32 triangleArea = 0.5f * D;
		area += triangleArea;

		// Area weighted centroid
		center += triangleArea * k_inv3 * (e1 + e2);

		float32 ex1 = e1.x, ey1 = e1.y;
		float32 ex2 = e2.x, ey2 = e2.y;

		float32 intx2 = ex1*ex1 + ex2*ex1 + ex2*ex2;
		float32 inty2 = ey1*ey1 + ey2*ey1 + ey2*ey2;

		I += (0.25f * k_inv3 * D) * (intx2 + inty2);
	}

	// Total mass
	massData->mass = density * area;

	// Center of mass
	b2Assert(area > b2_epsilon);
	center *= 1.0f / area;
	massData->center = center + s;

	// Inertia tensor relative to the local origin (point s).
	massData->I = density * I;
	
	// Shift to center of mass then to original body origin.
	massData->I += massData->mass * (b2Dot(massData->center, massData->center) - b2Dot(center, center));
}

bool b2PolygonShape::Validate() const
{
	for (int32 i = 0; i < m_count; ++i)
	{
		int32 i1 = i;
		int32 i2 = i < m_count - 1 ? i1 + 1 : 0;
		b2Vec2 p = m_vertices[i1];
		b2Vec2 e = m_vertices[i2] - p;

		for (int32 j = 0; j < m_count; ++j)
		{
			if (j == i1 || j == i2)
			{
				continue;
			}

			b2Vec2 v = m_vertices[j] - p;
			float32 c = b2Cross(e, v);
			if (c < 0.0f)
			{
				return false;
			}
		}
	}

	return true;
}

// end of PolygonShape.cpp

// GJK using Voronoi regions (Christer Ericson) and Barycentric coordinates.
int32 b2_gjkCalls, b2_gjkIters, b2_gjkMaxIters;

void b2DistanceProxy::Set(const b2Shape* shape, int32 index)
{
	switch (shape->GetType())
	{
	case b2Shape::e_circle:
		{
			const b2CircleShape* circle = static_cast<const b2CircleShape*>(shape);
			m_vertices = &circle->m_p;
			m_count = 1;
			m_radius = circle->m_radius;
		}
		break;

	case b2Shape::e_polygon:
		{
			const b2PolygonShape* polygon = static_cast<const b2PolygonShape*>(shape);
			m_vertices = polygon->m_vertices;
			m_count = polygon->m_count;
			m_radius = polygon->m_radius;
		}
		break;

	case b2Shape::e_chain:
		{
			const b2ChainShape* chain = static_cast<const b2ChainShape*>(shape);
			b2Assert(0 <= index && index < chain->m_count);

			m_buffer[0] = chain->m_vertices[index];
			if (index + 1 < chain->m_count)
			{
				m_buffer[1] = chain->m_vertices[index + 1];
			}
			else
			{
				m_buffer[1] = chain->m_vertices[0];
			}

			m_vertices = m_buffer;
			m_count = 2;
			m_radius = chain->m_radius;
		}
		break;

	case b2Shape::e_edge:
		{
			const b2EdgeShape* edge = static_cast<const b2EdgeShape*>(shape);
			m_vertices = &edge->m_vertex1;
			m_count = 2;
			m_radius = edge->m_radius;
		}
		break;

	default:
		b2Assert(false);
	}
}


struct b2SimplexVertex
{
	b2Vec2 wA;		// support point in proxyA
	b2Vec2 wB;		// support point in proxyB
	b2Vec2 w;		// wB - wA
	float32 a;		// barycentric coordinate for closest point
	int32 indexA;	// wA index
	int32 indexB;	// wB index
};

struct b2Simplex
{
	void ReadCache(	const b2SimplexCache* cache,
					const b2DistanceProxy* proxyA, const b2Transform& transformA,
					const b2DistanceProxy* proxyB, const b2Transform& transformB)
	{
		b2Assert(cache->count <= 3);
		
		// Copy data from cache.
		m_count = cache->count;
		b2SimplexVertex* vertices = &m_v1;
		for (int32 i = 0; i < m_count; ++i)
		{
			b2SimplexVertex* v = vertices + i;
			v->indexA = cache->indexA[i];
			v->indexB = cache->indexB[i];
			b2Vec2 wALocal = proxyA->GetVertex(v->indexA);
			b2Vec2 wBLocal = proxyB->GetVertex(v->indexB);
			v->wA = b2Mul(transformA, wALocal);
			v->wB = b2Mul(transformB, wBLocal);
			v->w = v->wB - v->wA;
			v->a = 0.0f;
		}

		// Compute the new simplex metric, if it is substantially different than
		// old metric then flush the simplex.
		if (m_count > 1)
		{
			float32 metric1 = cache->metric;
			float32 metric2 = GetMetric();
			if (metric2 < 0.5f * metric1 || 2.0f * metric1 < metric2 || metric2 < b2_epsilon)
			{
				// Reset the simplex.
				m_count = 0;
			}
		}

		// If the cache is empty or invalid ...
		if (m_count == 0)
		{
			b2SimplexVertex* v = vertices + 0;
			v->indexA = 0;
			v->indexB = 0;
			b2Vec2 wALocal = proxyA->GetVertex(0);
			b2Vec2 wBLocal = proxyB->GetVertex(0);
			v->wA = b2Mul(transformA, wALocal);
			v->wB = b2Mul(transformB, wBLocal);
			v->w = v->wB - v->wA;
			v->a = 1.0f;
			m_count = 1;
		}
	}

	void WriteCache(b2SimplexCache* cache) const
	{
		cache->metric = GetMetric();
		cache->count = uint16(m_count);
		const b2SimplexVertex* vertices = &m_v1;
		for (int32 i = 0; i < m_count; ++i)
		{
			cache->indexA[i] = uint8(vertices[i].indexA);
			cache->indexB[i] = uint8(vertices[i].indexB);
		}
	}

	b2Vec2 GetSearchDirection() const
	{
		switch (m_count)
		{
		case 1:
			return -m_v1.w;

		case 2:
			{
				b2Vec2 e12 = m_v2.w - m_v1.w;
				float32 sgn = b2Cross(e12, -m_v1.w);
				if (sgn > 0.0f)
				{
					// Origin is left of e12.
					return b2Cross(1.0f, e12);
				}
				else
				{
					// Origin is right of e12.
					return b2Cross(e12, 1.0f);
				}
			}

		default:
			b2Assert(false);
			return b2Vec2_zero;
		}
	}

	b2Vec2 GetClosestPoint() const
	{
		switch (m_count)
		{
		case 0:
			b2Assert(false);
			return b2Vec2_zero;

		case 1:
			return m_v1.w;

		case 2:
			return m_v1.a * m_v1.w + m_v2.a * m_v2.w;

		case 3:
			return b2Vec2_zero;

		default:
			b2Assert(false);
			return b2Vec2_zero;
		}
	}

	void GetWitnessPoints(b2Vec2* pA, b2Vec2* pB) const
	{
		switch (m_count)
		{
		case 0:
			b2Assert(false);
			break;

		case 1:
			*pA = m_v1.wA;
			*pB = m_v1.wB;
			break;

		case 2:
			*pA = m_v1.a * m_v1.wA + m_v2.a * m_v2.wA;
			*pB = m_v1.a * m_v1.wB + m_v2.a * m_v2.wB;
			break;

		case 3:
			*pA = m_v1.a * m_v1.wA + m_v2.a * m_v2.wA + m_v3.a * m_v3.wA;
			*pB = *pA;
			break;

		default:
			b2Assert(false);
			break;
		}
	}

	float32 GetMetric() const
	{
		switch (m_count)
		{
		case 0:
			b2Assert(false);
			return 0.0f;

		case 1:
			return 0.0f;

		case 2:
			return b2Distance(m_v1.w, m_v2.w);

		case 3:
			return b2Cross(m_v2.w - m_v1.w, m_v3.w - m_v1.w);

		default:
			b2Assert(false);
			return 0.0f;
		}
	}

	void Solve2();
	void Solve3();

	b2SimplexVertex m_v1, m_v2, m_v3;
	int32 m_count;
};


// Solve a line segment using barycentric coordinates.
//
// p = a1 * w1 + a2 * w2
// a1 + a2 = 1
//
// The vector from the origin to the closest point on the line is
// perpendicular to the line.
// e12 = w2 - w1
// dot(p, e) = 0
// a1 * dot(w1, e) + a2 * dot(w2, e) = 0
//
// 2-by-2 linear system
// [1      1     ][a1] = [1]
// [w1.e12 w2.e12][a2] = [0]
//
// Define
// d12_1 =  dot(w2, e12)
// d12_2 = -dot(w1, e12)
// d12 = d12_1 + d12_2
//
// Solution
// a1 = d12_1 / d12
// a2 = d12_2 / d12
void b2Simplex::Solve2()
{
	b2Vec2 w1 = m_v1.w;
	b2Vec2 w2 = m_v2.w;
	b2Vec2 e12 = w2 - w1;

	// w1 region
	float32 d12_2 = -b2Dot(w1, e12);
	if (d12_2 <= 0.0f)
	{
		// a2 <= 0, so we clamp it to 0
		m_v1.a = 1.0f;
		m_count = 1;
		return;
	}

	// w2 region
	float32 d12_1 = b2Dot(w2, e12);
	if (d12_1 <= 0.0f)
	{
		// a1 <= 0, so we clamp it to 0
		m_v2.a = 1.0f;
		m_count = 1;
		m_v1 = m_v2;
		return;
	}

	// Must be in e12 region.
	float32 inv_d12 = 1.0f / (d12_1 + d12_2);
	m_v1.a = d12_1 * inv_d12;
	m_v2.a = d12_2 * inv_d12;
	m_count = 2;
}

// Possible regions:
// - points[2]
// - edge points[0]-points[2]
// - edge points[1]-points[2]
// - inside the triangle
void b2Simplex::Solve3()
{
	b2Vec2 w1 = m_v1.w;
	b2Vec2 w2 = m_v2.w;
	b2Vec2 w3 = m_v3.w;

	// Edge12
	// [1      1     ][a1] = [1]
	// [w1.e12 w2.e12][a2] = [0]
	// a3 = 0
	b2Vec2 e12 = w2 - w1;
	float32 w1e12 = b2Dot(w1, e12);
	float32 w2e12 = b2Dot(w2, e12);
	float32 d12_1 = w2e12;
	float32 d12_2 = -w1e12;

	// Edge13
	// [1      1     ][a1] = [1]
	// [w1.e13 w3.e13][a3] = [0]
	// a2 = 0
	b2Vec2 e13 = w3 - w1;
	float32 w1e13 = b2Dot(w1, e13);
	float32 w3e13 = b2Dot(w3, e13);
	float32 d13_1 = w3e13;
	float32 d13_2 = -w1e13;

	// Edge23
	// [1      1     ][a2] = [1]
	// [w2.e23 w3.e23][a3] = [0]
	// a1 = 0
	b2Vec2 e23 = w3 - w2;
	float32 w2e23 = b2Dot(w2, e23);
	float32 w3e23 = b2Dot(w3, e23);
	float32 d23_1 = w3e23;
	float32 d23_2 = -w2e23;
	
	// Triangle123
	float32 n123 = b2Cross(e12, e13);

	float32 d123_1 = n123 * b2Cross(w2, w3);
	float32 d123_2 = n123 * b2Cross(w3, w1);
	float32 d123_3 = n123 * b2Cross(w1, w2);

	// w1 region
	if (d12_2 <= 0.0f && d13_2 <= 0.0f)
	{
		m_v1.a = 1.0f;
		m_count = 1;
		return;
	}

	// e12
	if (d12_1 > 0.0f && d12_2 > 0.0f && d123_3 <= 0.0f)
	{
		float32 inv_d12 = 1.0f / (d12_1 + d12_2);
		m_v1.a = d12_1 * inv_d12;
		m_v2.a = d12_2 * inv_d12;
		m_count = 2;
		return;
	}

	// e13
	if (d13_1 > 0.0f && d13_2 > 0.0f && d123_2 <= 0.0f)
	{
		float32 inv_d13 = 1.0f / (d13_1 + d13_2);
		m_v1.a = d13_1 * inv_d13;
		m_v3.a = d13_2 * inv_d13;
		m_count = 2;
		m_v2 = m_v3;
		return;
	}

	// w2 region
	if (d12_1 <= 0.0f && d23_2 <= 0.0f)
	{
		m_v2.a = 1.0f;
		m_count = 1;
		m_v1 = m_v2;
		return;
	}

	// w3 region
	if (d13_1 <= 0.0f && d23_1 <= 0.0f)
	{
		m_v3.a = 1.0f;
		m_count = 1;
		m_v1 = m_v3;
		return;
	}

	// e23
	if (d23_1 > 0.0f && d23_2 > 0.0f && d123_1 <= 0.0f)
	{
		float32 inv_d23 = 1.0f / (d23_1 + d23_2);
		m_v2.a = d23_1 * inv_d23;
		m_v3.a = d23_2 * inv_d23;
		m_count = 2;
		m_v1 = m_v3;
		return;
	}

	// Must be in triangle123
	float32 inv_d123 = 1.0f / (d123_1 + d123_2 + d123_3);
	m_v1.a = d123_1 * inv_d123;
	m_v2.a = d123_2 * inv_d123;
	m_v3.a = d123_3 * inv_d123;
	m_count = 3;
}

void b2Distance(b2DistanceOutput* output,
				b2SimplexCache* cache,
				const b2DistanceInput* input)
{
	++b2_gjkCalls;

	const b2DistanceProxy* proxyA = &input->proxyA;
	const b2DistanceProxy* proxyB = &input->proxyB;

	b2Transform transformA = input->transformA;
	b2Transform transformB = input->transformB;

	// Initialize the simplex.
	b2Simplex simplex;
	simplex.ReadCache(cache, proxyA, transformA, proxyB, transformB);

	// Get simplex vertices as an array.
	b2SimplexVertex* vertices = &simplex.m_v1;
	const int32 k_maxIters = 20;

	// These store the vertices of the last simplex so that we
	// can check for duplicates and prevent cycling.
	int32 saveA[3], saveB[3];
	int32 saveCount = 0;

	// Work around spurious gcc-4.8.2 warnings when -Wmaybe-uninitialized is
	// enabled by initializing saveA / saveB arrays when they're referenced
	// at the end of the main iteration loop below even though saveCount
	// entries of each array are initialized at the start of the main
	// iteration loop.
	memset(saveA, 0, sizeof(saveA));
	memset(saveB, 0, sizeof(saveB));

	float32 distanceSqr1 = b2_maxFloat;
	float32 distanceSqr2;

	// Main iteration loop.
	int32 iter = 0;
	while (iter < k_maxIters)
	{
		// Copy simplex so we can identify duplicates.
		saveCount = simplex.m_count;
		for (int32 i = 0; i < saveCount; ++i)
		{
			saveA[i] = vertices[i].indexA;
			saveB[i] = vertices[i].indexB;
		}

		switch (simplex.m_count)
		{
		case 1:
			break;

		case 2:
			simplex.Solve2();
			break;

		case 3:
			simplex.Solve3();
			break;

		default:
			b2Assert(false);
		}

		// If we have 3 points, then the origin is in the corresponding triangle.
		if (simplex.m_count == 3)
		{
			break;
		}

		// Compute closest point.
		b2Vec2 p = simplex.GetClosestPoint();
		distanceSqr2 = p.LengthSquared();

		// Ensure progress
		if (distanceSqr2 >= distanceSqr1)
		{
			//break;
		}
		distanceSqr1 = distanceSqr2;

		// Get search direction.
		b2Vec2 d = simplex.GetSearchDirection();

		// Ensure the search direction is numerically fit.
		if (d.LengthSquared() < b2_epsilon * b2_epsilon)
		{
			// The origin is probably contained by a line segment
			// or triangle. Thus the shapes are overlapped.

			// We can't return zero here even though there may be overlap.
			// In case the simplex is a point, segment, or triangle it is difficult
			// to determine if the origin is contained in the CSO or very close to it.
			break;
		}

		// Compute a tentative new simplex vertex using support points.
		b2SimplexVertex* vertex = vertices + simplex.m_count;
		vertex->indexA = proxyA->GetSupport(b2MulT(transformA.q, -d));
		vertex->wA = b2Mul(transformA, proxyA->GetVertex(vertex->indexA));
		b2Vec2 wBLocal;
		vertex->indexB = proxyB->GetSupport(b2MulT(transformB.q, d));
		vertex->wB = b2Mul(transformB, proxyB->GetVertex(vertex->indexB));
		vertex->w = vertex->wB - vertex->wA;

		// Iteration count is equated to the number of support point calls.
		++iter;
		++b2_gjkIters;

		// Check for duplicate support points. This is the main termination criteria.
		bool duplicate = false;
		for (int32 i = 0; i < saveCount; ++i)
		{
			if (vertex->indexA == saveA[i] && vertex->indexB == saveB[i])
			{
				duplicate = true;
				break;
			}
		}

		// If we found a duplicate support point we must exit to avoid cycling.
		if (duplicate)
		{
			break;
		}

		// New vertex is ok and needed.
		++simplex.m_count;
	}

	b2_gjkMaxIters = b2Max(b2_gjkMaxIters, iter);

	// Prepare output.
	simplex.GetWitnessPoints(&output->pointA, &output->pointB);
	output->distance = b2Distance(output->pointA, output->pointB);
	output->iterations = iter;

	// Cache the simplex.
	simplex.WriteCache(cache);

	// Apply radii if requested.
	if (input->useRadii)
	{
		float32 rA = proxyA->m_radius;
		float32 rB = proxyB->m_radius;

		if (output->distance > rA + rB && output->distance > b2_epsilon)
		{
			// Shapes are still no overlapped.
			// Move the witness points to the outer surface.
			output->distance -= rA + rB;
			b2Vec2 normal = output->pointB - output->pointA;
			normal.Normalize();
			output->pointA += rA * normal;
			output->pointB -= rB * normal;
		}
		else
		{
			// Shapes are overlapped when radii are considered.
			// Move the witness points to the middle.
			b2Vec2 p = 0.5f * (output->pointA + output->pointB);
			output->pointA = p;
			output->pointB = p;
			output->distance = 0.0f;
		}
	}
}

// end of Distance.cpp

#define b2_nullNode (-1)

/// A node in the dynamic tree. The client does not interact with this directly.
struct b2TreeNode
{
	bool IsLeaf() const
	{
		return child1 == b2_nullNode;
	}

	/// Enlarged AABB
	b2AABB aabb;

	void* userData;

	union
	{
		int32 parent;
		int32 next;
	};

	int32 child1;
	int32 child2;

	// leaf = 0, free node = -1
	int32 height;
};

/// A dynamic AABB tree broad-phase, inspired by Nathanael Presson's btDbvt.
/// A dynamic tree arranges data in a binary tree to accelerate
/// queries such as volume queries and ray casts. Leafs are proxies
/// with an AABB. In the tree we expand the proxy AABB by b2_fatAABBFactor
/// so that the proxy AABB is bigger than the client object. This allows the client
/// object to move by small amounts without triggering a tree update.
///
/// Nodes are pooled and relocatable, so we use node indices rather than pointers.
class b2DynamicTree
{
public:
	/// Constructing the tree initializes the node pool.
	b2DynamicTree();

	/// Destroy the tree, freeing the node pool.
	~b2DynamicTree();

	/// Create a proxy. Provide a tight fitting AABB and a userData pointer.
	int32 CreateProxy(const b2AABB& aabb, void* userData);

	/// Destroy a proxy. This asserts if the id is invalid.
	void DestroyProxy(int32 proxyId);

	/// Move a proxy with a swepted AABB. If the proxy has moved outside of its fattened AABB,
	/// then the proxy is removed from the tree and re-inserted. Otherwise
	/// the function returns immediately.
	/// @return true if the proxy was re-inserted.
	bool MoveProxy(int32 proxyId, const b2AABB& aabb1, const b2Vec2& displacement);

	/// Get proxy user data.
	/// @return the proxy user data or 0 if the id is invalid.
	void* GetUserData(int32 proxyId) const;

	/// Get the fat AABB for a proxy.
	const b2AABB& GetFatAABB(int32 proxyId) const;

	/// Query an AABB for overlapping proxies. The callback class
	/// is called for each proxy that overlaps the supplied AABB.
	template <typename T>
	void Query(T* callback, const b2AABB& aabb) const;

	/// Ray-cast against the proxies in the tree. This relies on the callback
	/// to perform a exact ray-cast in the case were the proxy contains a shape.
	/// The callback also performs the any collision filtering. This has performance
	/// roughly equal to k * log(n), where k is the number of collisions and n is the
	/// number of proxies in the tree.
	/// @param input the ray-cast input data. The ray extends from p1 to p1 + maxFraction * (p2 - p1).
	/// @param callback a callback class that is called for each proxy that is hit by the ray.
	template <typename T>
	void RayCast(T* callback, const b2RayCastInput& input) const;

	/// Validate this tree. For testing.
	void Validate() const;

	/// Compute the height of the binary tree in O(N) time. Should not be
	/// called often.
	int32 GetHeight() const;

	/// Get the maximum balance of an node in the tree. The balance is the difference
	/// in height of the two children of a node.
	int32 GetMaxBalance() const;

	/// Get the ratio of the sum of the node areas to the root area.
	float32 GetAreaRatio() const;

	/// Build an optimal tree. Very expensive. For testing.
	void RebuildBottomUp();

	/// Shift the world origin. Useful for large worlds.
	/// The shift formula is: position -= newOrigin
	/// @param newOrigin the new origin with respect to the old origin
	void ShiftOrigin(const b2Vec2& newOrigin);

private:

	int32 AllocateNode();
	void FreeNode(int32 node);

	void InsertLeaf(int32 node);
	void RemoveLeaf(int32 node);

	int32 Balance(int32 index);

	int32 ComputeHeight() const;
	int32 ComputeHeight(int32 nodeId) const;

	void ValidateStructure(int32 index) const;
	void ValidateMetrics(int32 index) const;

	int32 m_root;

	b2TreeNode* m_nodes;
	int32 m_nodeCount;
	int32 m_nodeCapacity;

	int32 m_freeList;

	/// This is used to incrementally traverse the tree for re-balancing.
	uint32 m_path;

	int32 m_insertionCount;
};

inline void* b2DynamicTree::GetUserData(int32 proxyId) const
{
	b2Assert(0 <= proxyId && proxyId < m_nodeCapacity);
	return m_nodes[proxyId].userData;
}

inline const b2AABB& b2DynamicTree::GetFatAABB(int32 proxyId) const
{
	b2Assert(0 <= proxyId && proxyId < m_nodeCapacity);
	return m_nodes[proxyId].aabb;
}

template <typename T>
inline void b2DynamicTree::Query(T* callback, const b2AABB& aabb) const
{
	b2GrowableStack<int32, 256> stack;
	stack.Push(m_root);

	while (stack.GetCount() > 0)
	{
		int32 nodeId = stack.Pop();
		if (nodeId == b2_nullNode)
		{
			continue;
		}

		const b2TreeNode* node = m_nodes + nodeId;

		if (b2TestOverlap(node->aabb, aabb))
		{
			if (node->IsLeaf())
			{
				bool proceed = callback->QueryCallback(nodeId);
				if (proceed == false)
				{
					return;
				}
			}
			else
			{
				stack.Push(node->child1);
				stack.Push(node->child2);
			}
		}
	}
}

template <typename T>
inline void b2DynamicTree::RayCast(T* callback, const b2RayCastInput& input) const
{
	b2Vec2 p1 = input.p1;
	b2Vec2 p2 = input.p2;
	b2Vec2 r = p2 - p1;
	b2Assert(r.LengthSquared() > 0.0f);
	r.Normalize();

	// v is perpendicular to the segment.
	b2Vec2 v = b2Cross(1.0f, r);
	b2Vec2 abs_v = b2Abs(v);

	// Separating axis for segment (Gino, p80).
	// |dot(v, p1 - c)| > dot(|v|, h)

	float32 maxFraction = input.maxFraction;

	// Build a bounding box for the segment.
	b2AABB segmentAABB;
	{
		b2Vec2 t = p1 + maxFraction * (p2 - p1);
		segmentAABB.lowerBound = b2Min(p1, t);
		segmentAABB.upperBound = b2Max(p1, t);
	}

	b2GrowableStack<int32, 256> stack;
	stack.Push(m_root);

	while (stack.GetCount() > 0)
	{
		int32 nodeId = stack.Pop();
		if (nodeId == b2_nullNode)
		{
			continue;
		}

		const b2TreeNode* node = m_nodes + nodeId;

		if (b2TestOverlap(node->aabb, segmentAABB) == false)
		{
			continue;
		}

		// Separating axis for segment (Gino, p80).
		// |dot(v, p1 - c)| > dot(|v|, h)
		b2Vec2 c = node->aabb.GetCenter();
		b2Vec2 h = node->aabb.GetExtents();
		float32 separation = b2Abs(b2Dot(v, p1 - c)) - b2Dot(abs_v, h);
		if (separation > 0.0f)
		{
			continue;
		}

		if (node->IsLeaf())
		{
			b2RayCastInput subInput;
			subInput.p1 = input.p1;
			subInput.p2 = input.p2;
			subInput.maxFraction = maxFraction;

			float32 value = callback->RayCastCallback(subInput, nodeId);

			if (value == 0.0f)
			{
				// The client has terminated the ray cast.
				return;
			}

			if (value > 0.0f)
			{
				// Update segment bounding box.
				maxFraction = value;
				b2Vec2 t = p1 + maxFraction * (p2 - p1);
				segmentAABB.lowerBound = b2Min(p1, t);
				segmentAABB.upperBound = b2Max(p1, t);
			}
		}
		else
		{
			stack.Push(node->child1);
			stack.Push(node->child2);
		}
	}
}

// end of DynamicTree.h

b2DynamicTree::b2DynamicTree()
{
	m_root = b2_nullNode;

	m_nodeCapacity = 16;
	m_nodeCount = 0;
	m_nodes = (b2TreeNode*)b2Alloc(m_nodeCapacity * sizeof(b2TreeNode));
	memset(m_nodes, 0, m_nodeCapacity * sizeof(b2TreeNode));

	// Build a linked list for the free list.
	for (int32 i = 0; i < m_nodeCapacity - 1; ++i)
	{
		m_nodes[i].next = i + 1;
		m_nodes[i].height = -1;
	}
	m_nodes[m_nodeCapacity-1].next = b2_nullNode;
	m_nodes[m_nodeCapacity-1].height = -1;
	m_freeList = 0;

	m_path = 0;

	m_insertionCount = 0;
}

b2DynamicTree::~b2DynamicTree()
{
	// This frees the entire tree in one shot.
	b2Free(m_nodes);
}

// Allocate a node from the pool. Grow the pool if necessary.
int32 b2DynamicTree::AllocateNode()
{
	// Expand the node pool as needed.
	if (m_freeList == b2_nullNode)
	{
		b2Assert(m_nodeCount == m_nodeCapacity);

		// The free list is empty. Rebuild a bigger pool.
		b2TreeNode* oldNodes = m_nodes;
		m_nodeCapacity *= 2;
		m_nodes = (b2TreeNode*)b2Alloc(m_nodeCapacity * sizeof(b2TreeNode));
		memcpy(m_nodes, oldNodes, m_nodeCount * sizeof(b2TreeNode));
		b2Free(oldNodes);

		// Build a linked list for the free list. The parent
		// pointer becomes the "next" pointer.
		for (int32 i = m_nodeCount; i < m_nodeCapacity - 1; ++i)
		{
			m_nodes[i].next = i + 1;
			m_nodes[i].height = -1;
		}
		m_nodes[m_nodeCapacity-1].next = b2_nullNode;
		m_nodes[m_nodeCapacity-1].height = -1;
		m_freeList = m_nodeCount;
	}

	// Peel a node off the free list.
	int32 nodeId = m_freeList;
	m_freeList = m_nodes[nodeId].next;
	m_nodes[nodeId].parent = b2_nullNode;
	m_nodes[nodeId].child1 = b2_nullNode;
	m_nodes[nodeId].child2 = b2_nullNode;
	m_nodes[nodeId].height = 0;
	m_nodes[nodeId].userData = NULL;
	++m_nodeCount;
	return nodeId;
}

// Return a node to the pool.
void b2DynamicTree::FreeNode(int32 nodeId)
{
	b2Assert(0 <= nodeId && nodeId < m_nodeCapacity);
	b2Assert(0 < m_nodeCount);
	m_nodes[nodeId].next = m_freeList;
	m_nodes[nodeId].height = -1;
	m_freeList = nodeId;
	--m_nodeCount;
}

// Create a proxy in the tree as a leaf node. We return the index
// of the node instead of a pointer so that we can grow
// the node pool.
int32 b2DynamicTree::CreateProxy(const b2AABB& aabb, void* userData)
{
	int32 proxyId = AllocateNode();

	// Fatten the aabb.
	b2Vec2 r(b2_aabbExtension, b2_aabbExtension);
	m_nodes[proxyId].aabb.lowerBound = aabb.lowerBound - r;
	m_nodes[proxyId].aabb.upperBound = aabb.upperBound + r;
	m_nodes[proxyId].userData = userData;
	m_nodes[proxyId].height = 0;

	InsertLeaf(proxyId);

	return proxyId;
}

void b2DynamicTree::DestroyProxy(int32 proxyId)
{
	b2Assert(0 <= proxyId && proxyId < m_nodeCapacity);
	b2Assert(m_nodes[proxyId].IsLeaf());

	RemoveLeaf(proxyId);
	FreeNode(proxyId);
}

bool b2DynamicTree::MoveProxy(int32 proxyId, const b2AABB& aabb, const b2Vec2& displacement)
{
	b2Assert(0 <= proxyId && proxyId < m_nodeCapacity);

	b2Assert(m_nodes[proxyId].IsLeaf());

	if (m_nodes[proxyId].aabb.Contains(aabb))
	{
		return false;
	}

	RemoveLeaf(proxyId);

	// Extend AABB.
	b2AABB b = aabb;
	b2Vec2 r(b2_aabbExtension, b2_aabbExtension);
	b.lowerBound = b.lowerBound - r;
	b.upperBound = b.upperBound + r;

	// Predict AABB displacement.
	b2Vec2 d = b2_aabbMultiplier * displacement;

	if (d.x < 0.0f)
	{
		b.lowerBound.x += d.x;
	}
	else
	{
		b.upperBound.x += d.x;
	}

	if (d.y < 0.0f)
	{
		b.lowerBound.y += d.y;
	}
	else
	{
		b.upperBound.y += d.y;
	}

	m_nodes[proxyId].aabb = b;

	InsertLeaf(proxyId);
	return true;
}

void b2DynamicTree::InsertLeaf(int32 leaf)
{
	++m_insertionCount;

	if (m_root == b2_nullNode)
	{
		m_root = leaf;
		m_nodes[m_root].parent = b2_nullNode;
		return;
	}

	// Find the best sibling for this node
	b2AABB leafAABB = m_nodes[leaf].aabb;
	int32 index = m_root;
	while (m_nodes[index].IsLeaf() == false)
	{
		int32 child1 = m_nodes[index].child1;
		int32 child2 = m_nodes[index].child2;

		float32 area = m_nodes[index].aabb.GetPerimeter();

		b2AABB combinedAABB;
		combinedAABB.Combine(m_nodes[index].aabb, leafAABB);
		float32 combinedArea = combinedAABB.GetPerimeter();

		// Cost of creating a new parent for this node and the new leaf
		float32 cost = 2.0f * combinedArea;

		// Minimum cost of pushing the leaf further down the tree
		float32 inheritanceCost = 2.0f * (combinedArea - area);

		// Cost of descending into child1
		float32 cost1;
		if (m_nodes[child1].IsLeaf())
		{
			b2AABB aabb;
			aabb.Combine(leafAABB, m_nodes[child1].aabb);
			cost1 = aabb.GetPerimeter() + inheritanceCost;
		}
		else
		{
			b2AABB aabb;
			aabb.Combine(leafAABB, m_nodes[child1].aabb);
			float32 oldArea = m_nodes[child1].aabb.GetPerimeter();
			float32 newArea = aabb.GetPerimeter();
			cost1 = (newArea - oldArea) + inheritanceCost;
		}

		// Cost of descending into child2
		float32 cost2;
		if (m_nodes[child2].IsLeaf())
		{
			b2AABB aabb;
			aabb.Combine(leafAABB, m_nodes[child2].aabb);
			cost2 = aabb.GetPerimeter() + inheritanceCost;
		}
		else
		{
			b2AABB aabb;
			aabb.Combine(leafAABB, m_nodes[child2].aabb);
			float32 oldArea = m_nodes[child2].aabb.GetPerimeter();
			float32 newArea = aabb.GetPerimeter();
			cost2 = newArea - oldArea + inheritanceCost;
		}

		// Descend according to the minimum cost.
		if (cost < cost1 && cost < cost2)
		{
			break;
		}

		// Descend
		if (cost1 < cost2)
		{
			index = child1;
		}
		else
		{
			index = child2;
		}
	}

	int32 sibling = index;

	// Create a new parent.
	int32 oldParent = m_nodes[sibling].parent;
	int32 newParent = AllocateNode();
	m_nodes[newParent].parent = oldParent;
	m_nodes[newParent].userData = NULL;
	m_nodes[newParent].aabb.Combine(leafAABB, m_nodes[sibling].aabb);
	m_nodes[newParent].height = m_nodes[sibling].height + 1;

	if (oldParent != b2_nullNode)
	{
		// The sibling was not the root.
		if (m_nodes[oldParent].child1 == sibling)
		{
			m_nodes[oldParent].child1 = newParent;
		}
		else
		{
			m_nodes[oldParent].child2 = newParent;
		}

		m_nodes[newParent].child1 = sibling;
		m_nodes[newParent].child2 = leaf;
		m_nodes[sibling].parent = newParent;
		m_nodes[leaf].parent = newParent;
	}
	else
	{
		// The sibling was the root.
		m_nodes[newParent].child1 = sibling;
		m_nodes[newParent].child2 = leaf;
		m_nodes[sibling].parent = newParent;
		m_nodes[leaf].parent = newParent;
		m_root = newParent;
	}

	// Walk back up the tree fixing heights and AABBs
	index = m_nodes[leaf].parent;
	while (index != b2_nullNode)
	{
		index = Balance(index);

		int32 child1 = m_nodes[index].child1;
		int32 child2 = m_nodes[index].child2;

		b2Assert(child1 != b2_nullNode);
		b2Assert(child2 != b2_nullNode);

		m_nodes[index].height = 1 + b2Max(m_nodes[child1].height, m_nodes[child2].height);
		m_nodes[index].aabb.Combine(m_nodes[child1].aabb, m_nodes[child2].aabb);

		index = m_nodes[index].parent;
	}

	//Validate();
}

void b2DynamicTree::RemoveLeaf(int32 leaf)
{
	if (leaf == m_root)
	{
		m_root = b2_nullNode;
		return;
	}

	int32 parent = m_nodes[leaf].parent;
	int32 grandParent = m_nodes[parent].parent;
	int32 sibling;
	if (m_nodes[parent].child1 == leaf)
	{
		sibling = m_nodes[parent].child2;
	}
	else
	{
		sibling = m_nodes[parent].child1;
	}

	if (grandParent != b2_nullNode)
	{
		// Destroy parent and connect sibling to grandParent.
		if (m_nodes[grandParent].child1 == parent)
		{
			m_nodes[grandParent].child1 = sibling;
		}
		else
		{
			m_nodes[grandParent].child2 = sibling;
		}
		m_nodes[sibling].parent = grandParent;
		FreeNode(parent);

		// Adjust ancestor bounds.
		int32 index = grandParent;
		while (index != b2_nullNode)
		{
			index = Balance(index);

			int32 child1 = m_nodes[index].child1;
			int32 child2 = m_nodes[index].child2;

			m_nodes[index].aabb.Combine(m_nodes[child1].aabb, m_nodes[child2].aabb);
			m_nodes[index].height = 1 + b2Max(m_nodes[child1].height, m_nodes[child2].height);

			index = m_nodes[index].parent;
		}
	}
	else
	{
		m_root = sibling;
		m_nodes[sibling].parent = b2_nullNode;
		FreeNode(parent);
	}

	//Validate();
}

// Perform a left or right rotation if node A is imbalanced.
// Returns the new root index.
int32 b2DynamicTree::Balance(int32 iA)
{
	b2Assert(iA != b2_nullNode);

	b2TreeNode* A = m_nodes + iA;
	if (A->IsLeaf() || A->height < 2)
	{
		return iA;
	}

	int32 iB = A->child1;
	int32 iC = A->child2;
	b2Assert(0 <= iB && iB < m_nodeCapacity);
	b2Assert(0 <= iC && iC < m_nodeCapacity);

	b2TreeNode* B = m_nodes + iB;
	b2TreeNode* C = m_nodes + iC;

	int32 balance = C->height - B->height;

	// Rotate C up
	if (balance > 1)
	{
		int32 iF = C->child1;
		int32 iG = C->child2;
		b2TreeNode* F = m_nodes + iF;
		b2TreeNode* G = m_nodes + iG;
		b2Assert(0 <= iF && iF < m_nodeCapacity);
		b2Assert(0 <= iG && iG < m_nodeCapacity);

		// Swap A and C
		C->child1 = iA;
		C->parent = A->parent;
		A->parent = iC;

		// A's old parent should point to C
		if (C->parent != b2_nullNode)
		{
			if (m_nodes[C->parent].child1 == iA)
			{
				m_nodes[C->parent].child1 = iC;
			}
			else
			{
				b2Assert(m_nodes[C->parent].child2 == iA);
				m_nodes[C->parent].child2 = iC;
			}
		}
		else
		{
			m_root = iC;
		}

		// Rotate
		if (F->height > G->height)
		{
			C->child2 = iF;
			A->child2 = iG;
			G->parent = iA;
			A->aabb.Combine(B->aabb, G->aabb);
			C->aabb.Combine(A->aabb, F->aabb);

			A->height = 1 + b2Max(B->height, G->height);
			C->height = 1 + b2Max(A->height, F->height);
		}
		else
		{
			C->child2 = iG;
			A->child2 = iF;
			F->parent = iA;
			A->aabb.Combine(B->aabb, F->aabb);
			C->aabb.Combine(A->aabb, G->aabb);

			A->height = 1 + b2Max(B->height, F->height);
			C->height = 1 + b2Max(A->height, G->height);
		}

		return iC;
	}
	
	// Rotate B up
	if (balance < -1)
	{
		int32 iD = B->child1;
		int32 iE = B->child2;
		b2TreeNode* D = m_nodes + iD;
		b2TreeNode* E = m_nodes + iE;
		b2Assert(0 <= iD && iD < m_nodeCapacity);
		b2Assert(0 <= iE && iE < m_nodeCapacity);

		// Swap A and B
		B->child1 = iA;
		B->parent = A->parent;
		A->parent = iB;

		// A's old parent should point to B
		if (B->parent != b2_nullNode)
		{
			if (m_nodes[B->parent].child1 == iA)
			{
				m_nodes[B->parent].child1 = iB;
			}
			else
			{
				b2Assert(m_nodes[B->parent].child2 == iA);
				m_nodes[B->parent].child2 = iB;
			}
		}
		else
		{
			m_root = iB;
		}

		// Rotate
		if (D->height > E->height)
		{
			B->child2 = iD;
			A->child1 = iE;
			E->parent = iA;
			A->aabb.Combine(C->aabb, E->aabb);
			B->aabb.Combine(A->aabb, D->aabb);

			A->height = 1 + b2Max(C->height, E->height);
			B->height = 1 + b2Max(A->height, D->height);
		}
		else
		{
			B->child2 = iE;
			A->child1 = iD;
			D->parent = iA;
			A->aabb.Combine(C->aabb, D->aabb);
			B->aabb.Combine(A->aabb, E->aabb);

			A->height = 1 + b2Max(C->height, D->height);
			B->height = 1 + b2Max(A->height, E->height);
		}

		return iB;
	}

	return iA;
}

int32 b2DynamicTree::GetHeight() const
{
	if (m_root == b2_nullNode)
	{
		return 0;
	}

	return m_nodes[m_root].height;
}

//
float32 b2DynamicTree::GetAreaRatio() const
{
	if (m_root == b2_nullNode)
	{
		return 0.0f;
	}

	const b2TreeNode* root = m_nodes + m_root;
	float32 rootArea = root->aabb.GetPerimeter();

	float32 totalArea = 0.0f;
	for (int32 i = 0; i < m_nodeCapacity; ++i)
	{
		const b2TreeNode* node = m_nodes + i;
		if (node->height < 0)
		{
			// Free node in pool
			continue;
		}

		totalArea += node->aabb.GetPerimeter();
	}

	return totalArea / rootArea;
}

// Compute the height of a sub-tree.
int32 b2DynamicTree::ComputeHeight(int32 nodeId) const
{
	b2Assert(0 <= nodeId && nodeId < m_nodeCapacity);
	b2TreeNode* node = m_nodes + nodeId;

	if (node->IsLeaf())
	{
		return 0;
	}

	int32 height1 = ComputeHeight(node->child1);
	int32 height2 = ComputeHeight(node->child2);
	return 1 + b2Max(height1, height2);
}

int32 b2DynamicTree::ComputeHeight() const
{
	int32 height = ComputeHeight(m_root);
	return height;
}

void b2DynamicTree::ValidateStructure(int32 index) const
{
	if (index == b2_nullNode)
	{
		return;
	}

	if (index == m_root)
	{
		b2Assert(m_nodes[index].parent == b2_nullNode);
	}

	const b2TreeNode* node = m_nodes + index;

#if B2_ASSERT_ENABLED || DEBUG
	int32 child1 = node->child1;
	int32 child2 = node->child2;
#endif  // B2_ASSERT_ENABLED || DEBUG

	if (node->IsLeaf())
	{
		b2Assert(child1 == b2_nullNode);
		b2Assert(child2 == b2_nullNode);
		b2Assert(node->height == 0);
		return;
	}

	b2Assert(0 <= child1 && child1 < m_nodeCapacity);
	b2Assert(0 <= child2 && child2 < m_nodeCapacity);

	b2Assert(m_nodes[child1].parent == index);
	b2Assert(m_nodes[child2].parent == index);

	B2_DEBUG_STATEMENT(ValidateStructure(child1));
	B2_DEBUG_STATEMENT(ValidateStructure(child2));
}

void b2DynamicTree::ValidateMetrics(int32 index) const
{
	if (index == b2_nullNode)
	{
		return;
	}

	const b2TreeNode* node = m_nodes + index;

	int32 child1 = node->child1;
	int32 child2 = node->child2;

	if (node->IsLeaf())
	{
		b2Assert(child1 == b2_nullNode);
		b2Assert(child2 == b2_nullNode);
		b2Assert(node->height == 0);
		return;
	}

	b2Assert(0 <= child1 && child1 < m_nodeCapacity);
	b2Assert(0 <= child2 && child2 < m_nodeCapacity);

#if B2_ASSERT_ENABLED
	int32 height1 = m_nodes[child1].height;
	int32 height2 = m_nodes[child2].height;
	int32 height;
	height = 1 + b2Max(height1, height2);
#endif // B2_ASSERT_ENABLED
	b2Assert(node->height == height);

	b2AABB aabb;
	aabb.Combine(m_nodes[child1].aabb, m_nodes[child2].aabb);

	b2Assert(aabb.lowerBound == node->aabb.lowerBound);
	b2Assert(aabb.upperBound == node->aabb.upperBound);

	ValidateMetrics(child1);
	ValidateMetrics(child2);
}

void b2DynamicTree::Validate() const
{
	B2_DEBUG_STATEMENT(ValidateStructure(m_root));
	B2_DEBUG_STATEMENT(ValidateMetrics(m_root));

	int32 freeCount = 0;
	int32 freeIndex = m_freeList;
	while (freeIndex != b2_nullNode)
	{
		b2Assert(0 <= freeIndex && freeIndex < m_nodeCapacity);
		freeIndex = m_nodes[freeIndex].next;
		++freeCount;
	}

	b2Assert(GetHeight() == ComputeHeight());

	b2Assert(m_nodeCount + freeCount == m_nodeCapacity);
}

int32 b2DynamicTree::GetMaxBalance() const
{
	int32 maxBalance = 0;
	for (int32 i = 0; i < m_nodeCapacity; ++i)
	{
		const b2TreeNode* node = m_nodes + i;
		if (node->height <= 1)
		{
			continue;
		}

		b2Assert(node->IsLeaf() == false);

		int32 child1 = node->child1;
		int32 child2 = node->child2;
		int32 balance = b2Abs(m_nodes[child2].height - m_nodes[child1].height);
		maxBalance = b2Max(maxBalance, balance);
	}

	return maxBalance;
}

void b2DynamicTree::RebuildBottomUp()
{
	int32* nodes = (int32*)b2Alloc(m_nodeCount * sizeof(int32));
	int32 count = 0;

	// Build array of leaves. Free the rest.
	for (int32 i = 0; i < m_nodeCapacity; ++i)
	{
		if (m_nodes[i].height < 0)
		{
			// free node in pool
			continue;
		}

		if (m_nodes[i].IsLeaf())
		{
			m_nodes[i].parent = b2_nullNode;
			nodes[count] = i;
			++count;
		}
		else
		{
			FreeNode(i);
		}
	}

	while (count > 1)
	{
		float32 minCost = b2_maxFloat;
		int32 iMin = -1, jMin = -1;
		for (int32 i = 0; i < count; ++i)
		{
			b2AABB aabbi = m_nodes[nodes[i]].aabb;

			for (int32 j = i + 1; j < count; ++j)
			{
				b2AABB aabbj = m_nodes[nodes[j]].aabb;
				b2AABB b;
				b.Combine(aabbi, aabbj);
				float32 cost = b.GetPerimeter();
				if (cost < minCost)
				{
					iMin = i;
					jMin = j;
					minCost = cost;
				}
			}
		}

		int32 index1 = nodes[iMin];
		int32 index2 = nodes[jMin];
		b2TreeNode* child1 = m_nodes + index1;
		b2TreeNode* child2 = m_nodes + index2;

		int32 parentIndex = AllocateNode();
		b2TreeNode* parent = m_nodes + parentIndex;
		parent->child1 = index1;
		parent->child2 = index2;
		parent->height = 1 + b2Max(child1->height, child2->height);
		parent->aabb.Combine(child1->aabb, child2->aabb);
		parent->parent = b2_nullNode;

		child1->parent = parentIndex;
		child2->parent = parentIndex;

		nodes[jMin] = nodes[count-1];
		nodes[iMin] = parentIndex;
		--count;
	}

	m_root = nodes[0];
	b2Free(nodes);

	B2_DEBUG_STATEMENT(Validate());
}

void b2DynamicTree::ShiftOrigin(const b2Vec2& newOrigin)
{
	// Build array of leaves. Free the rest.
	for (int32 i = 0; i < m_nodeCapacity; ++i)
	{
		m_nodes[i].aabb.lowerBound -= newOrigin;
		m_nodes[i].aabb.upperBound -= newOrigin;
	}
}

// end of DynamicTree.cpp

void b2WorldManifold::Initialize(const b2Manifold* manifold,
						  const b2Transform& xfA, float32 radiusA,
						  const b2Transform& xfB, float32 radiusB)
{
	if (manifold->pointCount == 0)
	{
		return;
	}

	switch (manifold->type)
	{
	case b2Manifold::e_circles:
		{
			normal.Set(1.0f, 0.0f);
			b2Vec2 pointA = b2Mul(xfA, manifold->localPoint);
			b2Vec2 pointB = b2Mul(xfB, manifold->points[0].localPoint);
			if (b2DistanceSquared(pointA, pointB) > b2_epsilon * b2_epsilon)
			{
				normal = pointB - pointA;
				normal.Normalize();
			}

			b2Vec2 cA = pointA + radiusA * normal;
			b2Vec2 cB = pointB - radiusB * normal;
			points[0] = 0.5f * (cA + cB);
			separations[0] = b2Dot(cB - cA, normal);
		}
		break;

	case b2Manifold::e_faceA:
		{
			normal = b2Mul(xfA.q, manifold->localNormal);
			b2Vec2 planePoint = b2Mul(xfA, manifold->localPoint);
			
			for (int32 i = 0; i < manifold->pointCount; ++i)
			{
				b2Vec2 clipPoint = b2Mul(xfB, manifold->points[i].localPoint);
				b2Vec2 cA = clipPoint + (radiusA - b2Dot(clipPoint - planePoint, normal)) * normal;
				b2Vec2 cB = clipPoint - radiusB * normal;
				points[i] = 0.5f * (cA + cB);
				separations[i] = b2Dot(cB - cA, normal);
			}
		}
		break;

	case b2Manifold::e_faceB:
		{
			normal = b2Mul(xfB.q, manifold->localNormal);
			b2Vec2 planePoint = b2Mul(xfB, manifold->localPoint);

			for (int32 i = 0; i < manifold->pointCount; ++i)
			{
				b2Vec2 clipPoint = b2Mul(xfA, manifold->points[i].localPoint);
				b2Vec2 cB = clipPoint + (radiusB - b2Dot(clipPoint - planePoint, normal)) * normal;
				b2Vec2 cA = clipPoint - radiusA * normal;
				points[i] = 0.5f * (cA + cB);
				separations[i] = b2Dot(cA - cB, normal);
			}

			// Ensure normal points from A to B.
			normal = -normal;
		}
		break;
	}
}

void b2GetPointStates(b2PointState state1[b2_maxManifoldPoints], b2PointState state2[b2_maxManifoldPoints],
					  const b2Manifold* manifold1, const b2Manifold* manifold2)
{
	for (int32 i = 0; i < b2_maxManifoldPoints; ++i)
	{
		state1[i] = b2_nullState;
		state2[i] = b2_nullState;
	}

	// Detect persists and removes.
	for (int32 i = 0; i < manifold1->pointCount; ++i)
	{
		b2ContactID id = manifold1->points[i].id;

		state1[i] = b2_removeState;

		for (int32 j = 0; j < manifold2->pointCount; ++j)
		{
			if (manifold2->points[j].id.key == id.key)
			{
				state1[i] = b2_persistState;
				break;
			}
		}
	}

	// Detect persists and adds.
	for (int32 i = 0; i < manifold2->pointCount; ++i)
	{
		b2ContactID id = manifold2->points[i].id;

		state2[i] = b2_addState;

		for (int32 j = 0; j < manifold1->pointCount; ++j)
		{
			if (manifold1->points[j].id.key == id.key)
			{
				state2[i] = b2_persistState;
				break;
			}
		}
	}
}

// From Real-time Collision Detection, p179.
bool b2AABB::RayCast(b2RayCastOutput* output, const b2RayCastInput& input) const
{
	float32 tmin = -b2_maxFloat;
	float32 tmax = b2_maxFloat;

	b2Vec2 p = input.p1;
	b2Vec2 d = input.p2 - input.p1;
	b2Vec2 absD = b2Abs(d);

	b2Vec2 normal;

	for (int32 i = 0; i < 2; ++i)
	{
		if (absD(i) < b2_epsilon)
		{
			// Parallel.
			if (p(i) < lowerBound(i) || upperBound(i) < p(i))
			{
				return false;
			}
		}
		else
		{
			float32 inv_d = 1.0f / d(i);
			float32 t1 = (lowerBound(i) - p(i)) * inv_d;
			float32 t2 = (upperBound(i) - p(i)) * inv_d;

			// Sign of the normal vector.
			float32 s = -1.0f;

			if (t1 > t2)
			{
				b2Swap(t1, t2);
				s = 1.0f;
			}

			// Push the min up
			if (t1 > tmin)
			{
				normal.SetZero();
				normal(i) = s;
				tmin = t1;
			}

			// Pull the max down
			tmax = b2Min(tmax, t2);

			if (tmin > tmax)
			{
				return false;
			}
		}
	}

	// Does the ray start inside the box?
	// Does the ray intersect beyond the max fraction?
	if (tmin < 0.0f || input.maxFraction < tmin)
	{
		return false;
	}

	// Intersection.
	output->fraction = tmin;
	output->normal = normal;
	return true;
}

// Sutherland-Hodgman clipping.
int32 b2ClipSegmentToLine(b2ClipVertex vOut[2], const b2ClipVertex vIn[2],
						const b2Vec2& normal, float32 offset, int32 vertexIndexA)
{
	// Start with no output points
	int32 numOut = 0;

	// Calculate the distance of end points to the line
	float32 distance0 = b2Dot(normal, vIn[0].v) - offset;
	float32 distance1 = b2Dot(normal, vIn[1].v) - offset;

	// If the points are behind the plane
	if (distance0 <= 0.0f) vOut[numOut++] = vIn[0];
	if (distance1 <= 0.0f) vOut[numOut++] = vIn[1];

	// If the points are on different sides of the plane
	if (distance0 * distance1 < 0.0f)
	{
		// Find intersection point of edge and plane
		float32 interp = distance0 / (distance0 - distance1);
		vOut[numOut].v = vIn[0].v + interp * (vIn[1].v - vIn[0].v);

		// VertexA is hitting edgeB.
		vOut[numOut].id.cf.indexA = static_cast<uint8>(vertexIndexA);
		vOut[numOut].id.cf.indexB = vIn[0].id.cf.indexB;
		vOut[numOut].id.cf.typeA = b2ContactFeature::e_vertex;
		vOut[numOut].id.cf.typeB = b2ContactFeature::e_face;
		++numOut;
	}

	return numOut;
}

bool b2TestOverlap(	const b2Shape* shapeA, int32 indexA,
					const b2Shape* shapeB, int32 indexB,
					const b2Transform& xfA, const b2Transform& xfB)
{
	b2DistanceInput input;
	input.proxyA.Set(shapeA, indexA);
	input.proxyB.Set(shapeB, indexB);
	input.transformA = xfA;
	input.transformB = xfB;
	input.useRadii = true;

	b2SimplexCache cache;
	cache.count = 0;

	b2DistanceOutput output;

	b2Distance(&output, &cache, &input);

	return output.distance < 10.0f * b2_epsilon;
}

// end of Collision.cpp

struct b2Pair
{
	int32 proxyIdA;
	int32 proxyIdB;
};

/// The broad-phase is used for computing pairs and performing volume queries and ray casts.
/// This broad-phase does not persist pairs. Instead, this reports potentially new pairs.
/// It is up to the client to consume the new pairs and to track subsequent overlap.
class b2BroadPhase
{
public:

	enum
	{
		e_nullProxy = -1
	};

	b2BroadPhase();
	~b2BroadPhase();

	/// Create a proxy with an initial AABB. Pairs are not reported until
	/// UpdatePairs is called.
	int32 CreateProxy(const b2AABB& aabb, void* userData);

	/// Destroy a proxy. It is up to the client to remove any pairs.
	void DestroyProxy(int32 proxyId);

	/// Call MoveProxy as many times as you like, then when you are done
	/// call UpdatePairs to finalized the proxy pairs (for your time step).
	void MoveProxy(int32 proxyId, const b2AABB& aabb, const b2Vec2& displacement);

	/// Call to trigger a re-processing of it's pairs on the next call to UpdatePairs.
	void TouchProxy(int32 proxyId);

	/// Get the fat AABB for a proxy.
	const b2AABB& GetFatAABB(int32 proxyId) const;

	/// Get user data from a proxy. Returns NULL if the id is invalid.
	void* GetUserData(int32 proxyId) const;

	/// Test overlap of fat AABBs.
	bool TestOverlap(int32 proxyIdA, int32 proxyIdB) const;

	/// Get the number of proxies.
	int32 GetProxyCount() const;

	/// Update the pairs. This results in pair callbacks. This can only add pairs.
	template <typename T>
	void UpdatePairs(T* callback);

	/// Query an AABB for overlapping proxies. The callback class
	/// is called for each proxy that overlaps the supplied AABB.
	template <typename T>
	void Query(T* callback, const b2AABB& aabb) const;

	/// Ray-cast against the proxies in the tree. This relies on the callback
	/// to perform a exact ray-cast in the case were the proxy contains a shape.
	/// The callback also performs the any collision filtering. This has performance
	/// roughly equal to k * log(n), where k is the number of collisions and n is the
	/// number of proxies in the tree.
	/// @param input the ray-cast input data. The ray extends from p1 to p1 + maxFraction * (p2 - p1).
	/// @param callback a callback class that is called for each proxy that is hit by the ray.
	template <typename T>
	void RayCast(T* callback, const b2RayCastInput& input) const;

	/// Get the height of the embedded tree.
	int32 GetTreeHeight() const;

	/// Get the balance of the embedded tree.
	int32 GetTreeBalance() const;

	/// Get the quality metric of the embedded tree.
	float32 GetTreeQuality() const;

	/// Shift the world origin. Useful for large worlds.
	/// The shift formula is: position -= newOrigin
	/// @param newOrigin the new origin with respect to the old origin
	void ShiftOrigin(const b2Vec2& newOrigin);

private:

	friend class b2DynamicTree;

	void BufferMove(int32 proxyId);
	void UnBufferMove(int32 proxyId);

	bool QueryCallback(int32 proxyId);

	b2DynamicTree m_tree;

	int32 m_proxyCount;

	int32* m_moveBuffer;
	int32 m_moveCapacity;
	int32 m_moveCount;

	b2Pair* m_pairBuffer;
	int32 m_pairCapacity;
	int32 m_pairCount;

	int32 m_queryProxyId;
};

/// This is used to sort pairs.
inline bool b2PairLessThan(const b2Pair& pair1, const b2Pair& pair2)
{
	if (pair1.proxyIdA < pair2.proxyIdA)
	{
		return true;
	}

	if (pair1.proxyIdA == pair2.proxyIdA)
	{
		return pair1.proxyIdB < pair2.proxyIdB;
	}

	return false;
}

inline void* b2BroadPhase::GetUserData(int32 proxyId) const
{
	return m_tree.GetUserData(proxyId);
}

inline bool b2BroadPhase::TestOverlap(int32 proxyIdA, int32 proxyIdB) const
{
	const b2AABB& aabbA = m_tree.GetFatAABB(proxyIdA);
	const b2AABB& aabbB = m_tree.GetFatAABB(proxyIdB);
	return b2TestOverlap(aabbA, aabbB);
}

inline const b2AABB& b2BroadPhase::GetFatAABB(int32 proxyId) const
{
	return m_tree.GetFatAABB(proxyId);
}

inline int32 b2BroadPhase::GetProxyCount() const
{
	return m_proxyCount;
}

inline int32 b2BroadPhase::GetTreeHeight() const
{
	return m_tree.GetHeight();
}

inline int32 b2BroadPhase::GetTreeBalance() const
{
	return m_tree.GetMaxBalance();
}

inline float32 b2BroadPhase::GetTreeQuality() const
{
	return m_tree.GetAreaRatio();
}

template <typename T>
void b2BroadPhase::UpdatePairs(T* callback)
{
	// Reset pair buffer
	m_pairCount = 0;

	// Perform tree queries for all moving proxies.
	for (int32 i = 0; i < m_moveCount; ++i)
	{
		m_queryProxyId = m_moveBuffer[i];
		if (m_queryProxyId == e_nullProxy)
		{
			continue;
		}

		// We have to query the tree with the fat AABB so that
		// we don't fail to create a pair that may touch later.
		const b2AABB& fatAABB = m_tree.GetFatAABB(m_queryProxyId);

		// Query tree, create pairs and add them pair buffer.
		m_tree.Query(this, fatAABB);
	}

	// Reset move buffer
	m_moveCount = 0;

	// Sort the pair buffer to expose duplicates.
	std::sort(m_pairBuffer, m_pairBuffer + m_pairCount, b2PairLessThan);

	// Send the pairs back to the client.
	int32 i = 0;
	while (i < m_pairCount)
	{
		b2Pair* primaryPair = m_pairBuffer + i;
		void* userDataA = m_tree.GetUserData(primaryPair->proxyIdA);
		void* userDataB = m_tree.GetUserData(primaryPair->proxyIdB);

		callback->AddPair(userDataA, userDataB);
		++i;

		// Skip any duplicate pairs.
		while (i < m_pairCount)
		{
			b2Pair* pair = m_pairBuffer + i;
			if (pair->proxyIdA != primaryPair->proxyIdA || pair->proxyIdB != primaryPair->proxyIdB)
			{
				break;
			}
			++i;
		}
	}

	// Try to keep the tree balanced.
	//m_tree.Rebalance(4);
}

template <typename T>
inline void b2BroadPhase::Query(T* callback, const b2AABB& aabb) const
{
	m_tree.Query(callback, aabb);
}

template <typename T>
inline void b2BroadPhase::RayCast(T* callback, const b2RayCastInput& input) const
{
	m_tree.RayCast(callback, input);
}

inline void b2BroadPhase::ShiftOrigin(const b2Vec2& newOrigin)
{
	m_tree.ShiftOrigin(newOrigin);
}

// end of BroadPhase.h

b2BroadPhase::b2BroadPhase()
{
	m_proxyCount = 0;

	m_pairCapacity = 16;
	m_pairCount = 0;
	m_pairBuffer = (b2Pair*)b2Alloc(m_pairCapacity * sizeof(b2Pair));

	m_moveCapacity = 16;
	m_moveCount = 0;
	m_moveBuffer = (int32*)b2Alloc(m_moveCapacity * sizeof(int32));
}

b2BroadPhase::~b2BroadPhase()
{
	b2Free(m_moveBuffer);
	b2Free(m_pairBuffer);
}

int32 b2BroadPhase::CreateProxy(const b2AABB& aabb, void* userData)
{
	int32 proxyId = m_tree.CreateProxy(aabb, userData);
	++m_proxyCount;
	BufferMove(proxyId);
	return proxyId;
}

void b2BroadPhase::DestroyProxy(int32 proxyId)
{
	UnBufferMove(proxyId);
	--m_proxyCount;
	m_tree.DestroyProxy(proxyId);
}

void b2BroadPhase::MoveProxy(int32 proxyId, const b2AABB& aabb, const b2Vec2& displacement)
{
	bool buffer = m_tree.MoveProxy(proxyId, aabb, displacement);
	if (buffer)
	{
		BufferMove(proxyId);
	}
}

void b2BroadPhase::TouchProxy(int32 proxyId)
{
	BufferMove(proxyId);
}

void b2BroadPhase::BufferMove(int32 proxyId)
{
	if (m_moveCount == m_moveCapacity)
	{
		int32* oldBuffer = m_moveBuffer;
		m_moveCapacity *= 2;
		m_moveBuffer = (int32*)b2Alloc(m_moveCapacity * sizeof(int32));
		memcpy(m_moveBuffer, oldBuffer, m_moveCount * sizeof(int32));
		b2Free(oldBuffer);
	}

	m_moveBuffer[m_moveCount] = proxyId;
	++m_moveCount;
}

void b2BroadPhase::UnBufferMove(int32 proxyId)
{
	for (int32 i = 0; i < m_moveCount; ++i)
	{
		if (m_moveBuffer[i] == proxyId)
		{
			m_moveBuffer[i] = e_nullProxy;
		}
	}
}

// This is called from b2DynamicTree::Query when we are gathering pairs.
bool b2BroadPhase::QueryCallback(int32 proxyId)
{
	// A proxy cannot form a pair with itself.
	if (proxyId == m_queryProxyId)
	{
		return true;
	}

	// Grow the pair buffer as needed.
	if (m_pairCount == m_pairCapacity)
	{
		b2Pair* oldBuffer = m_pairBuffer;
		m_pairCapacity *= 2;
		m_pairBuffer = (b2Pair*)b2Alloc(m_pairCapacity * sizeof(b2Pair));
		memcpy(m_pairBuffer, oldBuffer, m_pairCount * sizeof(b2Pair));
		b2Free(oldBuffer);
	}

	m_pairBuffer[m_pairCount].proxyIdA = b2Min(proxyId, m_queryProxyId);
	m_pairBuffer[m_pairCount].proxyIdB = b2Max(proxyId, m_queryProxyId);
	++m_pairCount;

	return true;
}

// end of BroadPhase.cpp

/// Input parameters for b2TimeOfImpact
struct b2TOIInput
{
	b2DistanceProxy proxyA;
	b2DistanceProxy proxyB;
	b2Sweep sweepA;
	b2Sweep sweepB;
	float32 tMax;		// defines sweep interval [0, tMax]
};

// Output parameters for b2TimeOfImpact.
struct b2TOIOutput
{
	enum State
	{
		e_unknown,
		e_failed,
		e_overlapped,
		e_touching,
		e_separated
	};

	State state;
	float32 t;
};

/// Compute the upper bound on time before two shapes penetrate. Time is represented as
/// a fraction between [0,tMax]. This uses a swept separating axis and may miss some intermediate,
/// non-tunneling collision. If you change the time interval, you should call this function
/// again.
/// Note: use b2Distance to compute the contact point and normal at the time of impact.
void b2TimeOfImpact(b2TOIOutput* output, const b2TOIInput* input);

// end of TimeOfImpact.h

float32 b2_toiTime, b2_toiMaxTime;
int32 b2_toiCalls, b2_toiIters, b2_toiMaxIters;
int32 b2_toiRootIters, b2_toiMaxRootIters;

//
struct b2SeparationFunction
{
	enum Type
	{
		e_points,
		e_faceA,
		e_faceB
	};

	// TODO_ERIN might not need to return the separation

	float32 Initialize(const b2SimplexCache* cache,
		const b2DistanceProxy* proxyA, const b2Sweep& sweepA,
		const b2DistanceProxy* proxyB, const b2Sweep& sweepB,
		float32 t1)
	{
		m_proxyA = proxyA;
		m_proxyB = proxyB;
		int32 count = cache->count;
		b2Assert(0 < count && count < 3);

		m_sweepA = sweepA;
		m_sweepB = sweepB;

		b2Transform xfA, xfB;
		m_sweepA.GetTransform(&xfA, t1);
		m_sweepB.GetTransform(&xfB, t1);

		if (count == 1)
		{
			m_type = e_points;
			b2Vec2 localPointA = m_proxyA->GetVertex(cache->indexA[0]);
			b2Vec2 localPointB = m_proxyB->GetVertex(cache->indexB[0]);
			b2Vec2 pointA = b2Mul(xfA, localPointA);
			b2Vec2 pointB = b2Mul(xfB, localPointB);
			m_axis = pointB - pointA;
			float32 s = m_axis.Normalize();
			m_localPoint = b2Vec2_zero;
			return s;
		}
		else if (cache->indexA[0] == cache->indexA[1])
		{
			// Two points on B and one on A.
			m_type = e_faceB;
			b2Vec2 localPointB1 = proxyB->GetVertex(cache->indexB[0]);
			b2Vec2 localPointB2 = proxyB->GetVertex(cache->indexB[1]);

			m_axis = b2Cross(localPointB2 - localPointB1, 1.0f);
			m_axis.Normalize();
			b2Vec2 normal = b2Mul(xfB.q, m_axis);

			m_localPoint = 0.5f * (localPointB1 + localPointB2);
			b2Vec2 pointB = b2Mul(xfB, m_localPoint);

			b2Vec2 localPointA = proxyA->GetVertex(cache->indexA[0]);
			b2Vec2 pointA = b2Mul(xfA, localPointA);

			float32 s = b2Dot(pointA - pointB, normal);
			if (s < 0.0f)
			{
				m_axis = -m_axis;
				s = -s;
			}
			return s;
		}
		else
		{
			// Two points on A and one or two points on B.
			m_type = e_faceA;
			b2Vec2 localPointA1 = m_proxyA->GetVertex(cache->indexA[0]);
			b2Vec2 localPointA2 = m_proxyA->GetVertex(cache->indexA[1]);
			
			m_axis = b2Cross(localPointA2 - localPointA1, 1.0f);
			m_axis.Normalize();
			b2Vec2 normal = b2Mul(xfA.q, m_axis);

			m_localPoint = 0.5f * (localPointA1 + localPointA2);
			b2Vec2 pointA = b2Mul(xfA, m_localPoint);

			b2Vec2 localPointB = m_proxyB->GetVertex(cache->indexB[0]);
			b2Vec2 pointB = b2Mul(xfB, localPointB);

			float32 s = b2Dot(pointB - pointA, normal);
			if (s < 0.0f)
			{
				m_axis = -m_axis;
				s = -s;
			}
			return s;
		}
	}

	//
	float32 FindMinSeparation(int32* indexA, int32* indexB, float32 t) const
	{
		b2Transform xfA, xfB;
		m_sweepA.GetTransform(&xfA, t);
		m_sweepB.GetTransform(&xfB, t);

		switch (m_type)
		{
		case e_points:
			{
				b2Vec2 axisA = b2MulT(xfA.q,  m_axis);
				b2Vec2 axisB = b2MulT(xfB.q, -m_axis);

				*indexA = m_proxyA->GetSupport(axisA);
				*indexB = m_proxyB->GetSupport(axisB);

				b2Vec2 localPointA = m_proxyA->GetVertex(*indexA);
				b2Vec2 localPointB = m_proxyB->GetVertex(*indexB);
				
				b2Vec2 pointA = b2Mul(xfA, localPointA);
				b2Vec2 pointB = b2Mul(xfB, localPointB);

				float32 separation = b2Dot(pointB - pointA, m_axis);
				return separation;
			}

		case e_faceA:
			{
				b2Vec2 normal = b2Mul(xfA.q, m_axis);
				b2Vec2 pointA = b2Mul(xfA, m_localPoint);

				b2Vec2 axisB = b2MulT(xfB.q, -normal);
				
				*indexA = -1;
				*indexB = m_proxyB->GetSupport(axisB);

				b2Vec2 localPointB = m_proxyB->GetVertex(*indexB);
				b2Vec2 pointB = b2Mul(xfB, localPointB);

				float32 separation = b2Dot(pointB - pointA, normal);
				return separation;
			}

		case e_faceB:
			{
				b2Vec2 normal = b2Mul(xfB.q, m_axis);
				b2Vec2 pointB = b2Mul(xfB, m_localPoint);

				b2Vec2 axisA = b2MulT(xfA.q, -normal);

				*indexB = -1;
				*indexA = m_proxyA->GetSupport(axisA);

				b2Vec2 localPointA = m_proxyA->GetVertex(*indexA);
				b2Vec2 pointA = b2Mul(xfA, localPointA);

				float32 separation = b2Dot(pointA - pointB, normal);
				return separation;
			}

		default:
			b2Assert(false);
			*indexA = -1;
			*indexB = -1;
			return 0.0f;
		}
	}

	//
	float32 Evaluate(int32 indexA, int32 indexB, float32 t) const
	{
		b2Transform xfA, xfB;
		m_sweepA.GetTransform(&xfA, t);
		m_sweepB.GetTransform(&xfB, t);

		switch (m_type)
		{
		case e_points:
			{
				b2Vec2 localPointA = m_proxyA->GetVertex(indexA);
				b2Vec2 localPointB = m_proxyB->GetVertex(indexB);

				b2Vec2 pointA = b2Mul(xfA, localPointA);
				b2Vec2 pointB = b2Mul(xfB, localPointB);
				float32 separation = b2Dot(pointB - pointA, m_axis);

				return separation;
			}

		case e_faceA:
			{
				b2Vec2 normal = b2Mul(xfA.q, m_axis);
				b2Vec2 pointA = b2Mul(xfA, m_localPoint);

				b2Vec2 localPointB = m_proxyB->GetVertex(indexB);
				b2Vec2 pointB = b2Mul(xfB, localPointB);

				float32 separation = b2Dot(pointB - pointA, normal);
				return separation;
			}

		case e_faceB:
			{
				b2Vec2 normal = b2Mul(xfB.q, m_axis);
				b2Vec2 pointB = b2Mul(xfB, m_localPoint);

				b2Vec2 localPointA = m_proxyA->GetVertex(indexA);
				b2Vec2 pointA = b2Mul(xfA, localPointA);

				float32 separation = b2Dot(pointA - pointB, normal);
				return separation;
			}

		default:
			b2Assert(false);
			return 0.0f;
		}
	}

	const b2DistanceProxy* m_proxyA;
	const b2DistanceProxy* m_proxyB;
	b2Sweep m_sweepA, m_sweepB;
	Type m_type;
	b2Vec2 m_localPoint;
	b2Vec2 m_axis;
};

// CCD via the local separating axis method. This seeks progression
// by computing the largest time at which separation is maintained.
void b2TimeOfImpact(b2TOIOutput* output, const b2TOIInput* input)
{
	b2Timer timer;

	++b2_toiCalls;

	output->state = b2TOIOutput::e_unknown;
	output->t = input->tMax;

	const b2DistanceProxy* proxyA = &input->proxyA;
	const b2DistanceProxy* proxyB = &input->proxyB;

	b2Sweep sweepA = input->sweepA;
	b2Sweep sweepB = input->sweepB;

	// Large rotations can make the root finder fail, so we normalize the
	// sweep angles.
	sweepA.Normalize();
	sweepB.Normalize();

	float32 tMax = input->tMax;

	float32 totalRadius = proxyA->m_radius + proxyB->m_radius;
	float32 target = b2Max(b2_linearSlop, totalRadius - 3.0f * b2_linearSlop);
	float32 tolerance = 0.25f * b2_linearSlop;
	b2Assert(target > tolerance);

	float32 t1 = 0.0f;
	const int32 k_maxIterations = 20;	// TODO_ERIN b2Settings
	int32 iter = 0;

	// Prepare input for distance query.
	b2SimplexCache cache;
	cache.count = 0;
	b2DistanceInput distanceInput;
	distanceInput.proxyA = input->proxyA;
	distanceInput.proxyB = input->proxyB;
	distanceInput.useRadii = false;

	// The outer loop progressively attempts to compute new separating axes.
	// This loop terminates when an axis is repeated (no progress is made).
	for(;;)
	{
		b2Transform xfA, xfB;
		sweepA.GetTransform(&xfA, t1);
		sweepB.GetTransform(&xfB, t1);

		// Get the distance between shapes. We can also use the results
		// to get a separating axis.
		distanceInput.transformA = xfA;
		distanceInput.transformB = xfB;
		b2DistanceOutput distanceOutput;
		b2Distance(&distanceOutput, &cache, &distanceInput);

		// If the shapes are overlapped, we give up on continuous collision.
		if (distanceOutput.distance <= 0.0f)
		{
			// Failure!
			output->state = b2TOIOutput::e_overlapped;
			output->t = 0.0f;
			break;
		}

		if (distanceOutput.distance < target + tolerance)
		{
			// Victory!
			output->state = b2TOIOutput::e_touching;
			output->t = t1;
			break;
		}

		// Initialize the separating axis.
		b2SeparationFunction fcn;
		fcn.Initialize(&cache, proxyA, sweepA, proxyB, sweepB, t1);
#if 0
		// Dump the curve seen by the root finder
		{
			const int32 N = 100;
			float32 dx = 1.0f / N;
			float32 xs[N+1];
			float32 fs[N+1];

			float32 x = 0.0f;

			for (int32 i = 0; i <= N; ++i)
			{
				sweepA.GetTransform(&xfA, x);
				sweepB.GetTransform(&xfB, x);
				float32 f = fcn.Evaluate(xfA, xfB) - target;

				printf("%g %g\n", x, f);

				xs[i] = x;
				fs[i] = f;

				x += dx;
			}
		}
#endif

		// Compute the TOI on the separating axis. We do this by successively
		// resolving the deepest point. This loop is bounded by the number of vertices.
		bool done = false;
		float32 t2 = tMax;
		int32 pushBackIter = 0;
		for (;;)
		{
			// Find the deepest point at t2. Store the witness point indices.
			int32 indexA, indexB;
			float32 s2 = fcn.FindMinSeparation(&indexA, &indexB, t2);

			// Is the final configuration separated?
			if (s2 > target + tolerance)
			{
				// Victory!
				output->state = b2TOIOutput::e_separated;
				output->t = tMax;
				done = true;
				break;
			}

			// Has the separation reached tolerance?
			if (s2 > target - tolerance)
			{
				// Advance the sweeps
				t1 = t2;
				break;
			}

			// Compute the initial separation of the witness points.
			float32 s1 = fcn.Evaluate(indexA, indexB, t1);

			// Check for initial overlap. This might happen if the root finder
			// runs out of iterations.
			if (s1 < target - tolerance)
			{
				output->state = b2TOIOutput::e_failed;
				output->t = t1;
				done = true;
				break;
			}

			// Check for touching
			if (s1 <= target + tolerance)
			{
				// Victory! t1 should hold the TOI (could be 0.0).
				output->state = b2TOIOutput::e_touching;
				output->t = t1;
				done = true;
				break;
			}

			// Compute 1D root of: f(x) - target = 0
			int32 rootIterCount = 0;
			float32 a1 = t1, a2 = t2;
			for (;;)
			{
				// Use a mix of the secant rule and bisection.
				float32 t;
				if (rootIterCount & 1)
				{
					// Secant rule to improve convergence.
					t = a1 + (target - s1) * (a2 - a1) / (s2 - s1);
				}
				else
				{
					// Bisection to guarantee progress.
					t = 0.5f * (a1 + a2);
				}

				++rootIterCount;
				++b2_toiRootIters;

				float32 s = fcn.Evaluate(indexA, indexB, t);

				if (b2Abs(s - target) < tolerance)
				{
					// t2 holds a tentative value for t1
					t2 = t;
					break;
				}

				// Ensure we continue to bracket the root.
				if (s > target)
				{
					a1 = t;
					s1 = s;
				}
				else
				{
					a2 = t;
					s2 = s;
				}
				
				if (rootIterCount == 50)
				{
					break;
				}
			}

			b2_toiMaxRootIters = b2Max(b2_toiMaxRootIters, rootIterCount);

			++pushBackIter;

			if (pushBackIter == b2_maxPolygonVertices)
			{
				break;
			}
		}

		++iter;
		++b2_toiIters;

		if (done)
		{
			break;
		}

		if (iter == k_maxIterations)
		{
			// Root finder got stuck. Semi-victory.
			output->state = b2TOIOutput::e_failed;
			output->t = t1;
			break;
		}
	}

	b2_toiMaxIters = b2Max(b2_toiMaxIters, iter);

	float32 time = timer.GetMilliseconds();
	b2_toiMaxTime = b2Max(b2_toiMaxTime, time);
	b2_toiTime += time;
}

// end of TimeOfImpact.cpp

// Find the max separation between poly1 and poly2 using edge normals from poly1.
static float32 b2FindMaxSeparation(int32* edgeIndex,
								 const b2PolygonShape* poly1, const b2Transform& xf1,
								 const b2PolygonShape* poly2, const b2Transform& xf2)
{
	int32 count1 = poly1->m_count;
	int32 count2 = poly2->m_count;
	const b2Vec2* n1s = poly1->m_normals;
	const b2Vec2* v1s = poly1->m_vertices;
	const b2Vec2* v2s = poly2->m_vertices;
	b2Transform xf = b2MulT(xf2, xf1);

	int32 bestIndex = 0;
	float32 maxSeparation = -b2_maxFloat;
	for (int32 i = 0; i < count1; ++i)
	{
		// Get poly1 normal in frame2.
		b2Vec2 n = b2Mul(xf.q, n1s[i]);
		b2Vec2 v1 = b2Mul(xf, v1s[i]);

		// Find deepest point for normal i.
		float32 si = b2_maxFloat;
		for (int32 j = 0; j < count2; ++j)
		{
			float32 sij = b2Dot(n, v2s[j] - v1);
			if (sij < si)
			{
				si = sij;
			}
		}

		if (si > maxSeparation)
		{
			maxSeparation = si;
			bestIndex = i;
		}
	}

	*edgeIndex = bestIndex;
	return maxSeparation;
}

static void b2FindIncidentEdge(b2ClipVertex c[2],
							 const b2PolygonShape* poly1, const b2Transform& xf1, int32 edge1,
							 const b2PolygonShape* poly2, const b2Transform& xf2)
{
	const b2Vec2* normals1 = poly1->m_normals;

	int32 count2 = poly2->m_count;
	const b2Vec2* vertices2 = poly2->m_vertices;
	const b2Vec2* normals2 = poly2->m_normals;

	b2Assert(0 <= edge1 && edge1 < poly1->m_count);

	// Get the normal of the reference edge in poly2's frame.
	b2Vec2 normal1 = b2MulT(xf2.q, b2Mul(xf1.q, normals1[edge1]));

	// Find the incident edge on poly2.
	int32 index = 0;
	float32 minDot = b2_maxFloat;
	for (int32 i = 0; i < count2; ++i)
	{
		float32 dot = b2Dot(normal1, normals2[i]);
		if (dot < minDot)
		{
			minDot = dot;
			index = i;
		}
	}

	// Build the clip vertices for the incident edge.
	int32 i1 = index;
	int32 i2 = i1 + 1 < count2 ? i1 + 1 : 0;

	c[0].v = b2Mul(xf2, vertices2[i1]);
	c[0].id.cf.indexA = (uint8)edge1;
	c[0].id.cf.indexB = (uint8)i1;
	c[0].id.cf.typeA = b2ContactFeature::e_face;
	c[0].id.cf.typeB = b2ContactFeature::e_vertex;

	c[1].v = b2Mul(xf2, vertices2[i2]);
	c[1].id.cf.indexA = (uint8)edge1;
	c[1].id.cf.indexB = (uint8)i2;
	c[1].id.cf.typeA = b2ContactFeature::e_face;
	c[1].id.cf.typeB = b2ContactFeature::e_vertex;
}

// Find edge normal of max separation on A - return if separating axis is found
// Find edge normal of max separation on B - return if separation axis is found
// Choose reference edge as min(minA, minB)
// Find incident edge
// Clip

// The normal points from 1 to 2
void b2CollidePolygons(b2Manifold* manifold,
					  const b2PolygonShape* polyA, const b2Transform& xfA,
					  const b2PolygonShape* polyB, const b2Transform& xfB)
{
	manifold->pointCount = 0;
	float32 totalRadius = polyA->m_radius + polyB->m_radius;

	int32 edgeA = 0;
	float32 separationA = b2FindMaxSeparation(&edgeA, polyA, xfA, polyB, xfB);
	if (separationA > totalRadius)
		return;

	int32 edgeB = 0;
	float32 separationB = b2FindMaxSeparation(&edgeB, polyB, xfB, polyA, xfA);
	if (separationB > totalRadius)
		return;

	const b2PolygonShape* poly1;	// reference polygon
	const b2PolygonShape* poly2;	// incident polygon
	b2Transform xf1, xf2;
	int32 edge1;					// reference edge
	uint8 flip;
	const float32 k_tol = 0.1f * b2_linearSlop;

	if (separationB > separationA + k_tol)
	{
		poly1 = polyB;
		poly2 = polyA;
		xf1 = xfB;
		xf2 = xfA;
		edge1 = edgeB;
		manifold->type = b2Manifold::e_faceB;
		flip = 1;
	}
	else
	{
		poly1 = polyA;
		poly2 = polyB;
		xf1 = xfA;
		xf2 = xfB;
		edge1 = edgeA;
		manifold->type = b2Manifold::e_faceA;
		flip = 0;
	}

	b2ClipVertex incidentEdge[2];
	b2FindIncidentEdge(incidentEdge, poly1, xf1, edge1, poly2, xf2);

	int32 count1 = poly1->m_count;
	const b2Vec2* vertices1 = poly1->m_vertices;

	int32 iv1 = edge1;
	int32 iv2 = edge1 + 1 < count1 ? edge1 + 1 : 0;

	b2Vec2 v11 = vertices1[iv1];
	b2Vec2 v12 = vertices1[iv2];

	b2Vec2 localTangent = v12 - v11;
	localTangent.Normalize();
	
	b2Vec2 localNormal = b2Cross(localTangent, 1.0f);
	b2Vec2 planePoint = 0.5f * (v11 + v12);

	b2Vec2 tangent = b2Mul(xf1.q, localTangent);
	b2Vec2 normal = b2Cross(tangent, 1.0f);
	
	v11 = b2Mul(xf1, v11);
	v12 = b2Mul(xf1, v12);

	// Face offset.
	float32 frontOffset = b2Dot(normal, v11);

	// Side offsets, extended by polytope skin thickness.
	float32 sideOffset1 = -b2Dot(tangent, v11) + totalRadius;
	float32 sideOffset2 = b2Dot(tangent, v12) + totalRadius;

	// Clip incident edge against extruded edge1 side edges.
	b2ClipVertex clipPoints1[2];
	b2ClipVertex clipPoints2[2];
	int np;

	// Clip to box side 1
	np = b2ClipSegmentToLine(clipPoints1, incidentEdge, -tangent, sideOffset1, iv1);

	if (np < 2)
		return;

	// Clip to negative box side 1
	np = b2ClipSegmentToLine(clipPoints2, clipPoints1,  tangent, sideOffset2, iv2);

	if (np < 2)
	{
		return;
	}

	// Now clipPoints2 contains the clipped points.
	manifold->localNormal = localNormal;
	manifold->localPoint = planePoint;

	int32 pointCount = 0;
	for (int32 i = 0; i < b2_maxManifoldPoints; ++i)
	{
		float32 separation = b2Dot(normal, clipPoints2[i].v) - frontOffset;

		if (separation <= totalRadius)
		{
			b2ManifoldPoint* cp = manifold->points + pointCount;
			cp->localPoint = b2MulT(xf2, clipPoints2[i].v);
			cp->id = clipPoints2[i].id;
			if (flip)
			{
				// Swap features
				b2ContactFeature cf = cp->id.cf;
				cp->id.cf.indexA = cf.indexB;
				cp->id.cf.indexB = cf.indexA;
				cp->id.cf.typeA = cf.typeB;
				cp->id.cf.typeB = cf.typeA;
			}
			++pointCount;
		}
	}

	manifold->pointCount = pointCount;
}

// end of CollidePolygon.cpp

class b2Draw;

/// 
struct b2RopeDef
{
	b2RopeDef()
	{
		vertices = NULL;
		count = 0;
		masses = NULL;
		gravity.SetZero();
		damping = 0.1f;
		k2 = 0.9f;
		k3 = 0.1f;
	}

	///
	b2Vec2* vertices;

	///
	int32 count;

	///
	float32* masses;

	///
	b2Vec2 gravity;

	///
	float32 damping;

	/// Stretching stiffness
	float32 k2;

	/// Bending stiffness. Values above 0.5 can make the simulation blow up.
	float32 k3;
};

/// 
class b2Rope
{
public:
	b2Rope();
	~b2Rope();

	///
	void Initialize(const b2RopeDef* def);

	///
	void Step(float32 timeStep, int32 iterations);

	///
	int32 GetVertexCount() const
	{
		return m_count;
	}

	///
	const b2Vec2* GetVertices() const
	{
		return m_ps;
	}

	///
	void Draw(b2Draw* draw) const;

	///
	void SetAngle(float32 angle);

private:

	void SolveC2();
	void SolveC3();

	int32 m_count;
	b2Vec2* m_ps;
	b2Vec2* m_p0s;
	b2Vec2* m_vs;

	float32* m_ims;

	float32* m_Ls;
	float32* m_as;

	b2Vec2 m_gravity;
	float32 m_damping;

	float32 m_k2;
	float32 m_k3;
};

// end of Rope.h

b2Rope::b2Rope()
{
	m_count = 0;
	m_ps = NULL;
	m_p0s = NULL;
	m_vs = NULL;
	m_ims = NULL;
	m_Ls = NULL;
	m_as = NULL;
	m_gravity.SetZero();
	m_k2 = 1.0f;
	m_k3 = 0.1f;
}

b2Rope::~b2Rope()
{
	b2Free(m_ps);
	b2Free(m_p0s);
	b2Free(m_vs);
	b2Free(m_ims);
	b2Free(m_Ls);
	b2Free(m_as);
}

void b2Rope::Initialize(const b2RopeDef* def)
{
	b2Assert(def->count >= 3);
	m_count = def->count;
	m_ps = (b2Vec2*)b2Alloc(m_count * sizeof(b2Vec2));
	m_p0s = (b2Vec2*)b2Alloc(m_count * sizeof(b2Vec2));
	m_vs = (b2Vec2*)b2Alloc(m_count * sizeof(b2Vec2));
	m_ims = (float32*)b2Alloc(m_count * sizeof(float32));

	for (int32 i = 0; i < m_count; ++i)
	{
		m_ps[i] = def->vertices[i];
		m_p0s[i] = def->vertices[i];
		m_vs[i].SetZero();

		float32 m = def->masses[i];
		if (m > 0.0f)
		{
			m_ims[i] = 1.0f / m;
		}
		else
		{
			m_ims[i] = 0.0f;
		}
	}

	int32 count2 = m_count - 1;
	int32 count3 = m_count - 2;
	m_Ls = (float32*)b2Alloc(count2 * sizeof(float32));
	m_as = (float32*)b2Alloc(count3 * sizeof(float32));

	for (int32 i = 0; i < count2; ++i)
	{
		b2Vec2 p1 = m_ps[i];
		b2Vec2 p2 = m_ps[i+1];
		m_Ls[i] = b2Distance(p1, p2);
	}

	for (int32 i = 0; i < count3; ++i)
	{
		b2Vec2 p1 = m_ps[i];
		b2Vec2 p2 = m_ps[i + 1];
		b2Vec2 p3 = m_ps[i + 2];

		b2Vec2 d1 = p2 - p1;
		b2Vec2 d2 = p3 - p2;

		float32 a = b2Cross(d1, d2);
		float32 b = b2Dot(d1, d2);

		m_as[i] = b2Atan2(a, b);
	}

	m_gravity = def->gravity;
	m_damping = def->damping;
	m_k2 = def->k2;
	m_k3 = def->k3;
}

void b2Rope::Step(float32 h, int32 iterations)
{
	if (h == 0.0)
	{
		return;
	}

	float32 d = expf(- h * m_damping);

	for (int32 i = 0; i < m_count; ++i)
	{
		m_p0s[i] = m_ps[i];
		if (m_ims[i] > 0.0f)
		{
			m_vs[i] += h * m_gravity;
		}
		m_vs[i] *= d;
		m_ps[i] += h * m_vs[i];

	}

	for (int32 i = 0; i < iterations; ++i)
	{
		SolveC2();
		SolveC3();
		SolveC2();
	}

	float32 inv_h = 1.0f / h;
	for (int32 i = 0; i < m_count; ++i)
	{
		m_vs[i] = inv_h * (m_ps[i] - m_p0s[i]);
	}
}

void b2Rope::SolveC2()
{
	int32 count2 = m_count - 1;

	for (int32 i = 0; i < count2; ++i)
	{
		b2Vec2 p1 = m_ps[i];
		b2Vec2 p2 = m_ps[i + 1];

		b2Vec2 d = p2 - p1;
		float32 L = d.Normalize();

		float32 im1 = m_ims[i];
		float32 im2 = m_ims[i + 1];

		if (im1 + im2 == 0.0f)
		{
			continue;
		}

		float32 s1 = im1 / (im1 + im2);
		float32 s2 = im2 / (im1 + im2);

		p1 -= m_k2 * s1 * (m_Ls[i] - L) * d;
		p2 += m_k2 * s2 * (m_Ls[i] - L) * d;

		m_ps[i] = p1;
		m_ps[i + 1] = p2;
	}
}

void b2Rope::SetAngle(float32 angle)
{
	int32 count3 = m_count - 2;
	for (int32 i = 0; i < count3; ++i)
	{
		m_as[i] = angle;
	}
}

void b2Rope::SolveC3()
{
	int32 count3 = m_count - 2;

	for (int32 i = 0; i < count3; ++i)
	{
		b2Vec2 p1 = m_ps[i];
		b2Vec2 p2 = m_ps[i + 1];
		b2Vec2 p3 = m_ps[i + 2];

		float32 m1 = m_ims[i];
		float32 m2 = m_ims[i + 1];
		float32 m3 = m_ims[i + 2];

		b2Vec2 d1 = p2 - p1;
		b2Vec2 d2 = p3 - p2;

		float32 L1sqr = d1.LengthSquared();
		float32 L2sqr = d2.LengthSquared();

		if (L1sqr * L2sqr == 0.0f)
		{
			continue;
		}

		float32 a = b2Cross(d1, d2);
		float32 b = b2Dot(d1, d2);

		float32 angle = b2Atan2(a, b);

		b2Vec2 Jd1 = (-1.0f / L1sqr) * d1.Skew();
		b2Vec2 Jd2 = (1.0f / L2sqr) * d2.Skew();

		b2Vec2 J1 = -Jd1;
		b2Vec2 J2 = Jd1 - Jd2;
		b2Vec2 J3 = Jd2;

		float32 mass = m1 * b2Dot(J1, J1) + m2 * b2Dot(J2, J2) + m3 * b2Dot(J3, J3);
		if (mass == 0.0f)
		{
			continue;
		}

		mass = 1.0f / mass;

		float32 C = angle - m_as[i];

		while (C > b2_pi)
		{
			angle -= 2 * b2_pi;
			C = angle - m_as[i];
		}

		while (C < -b2_pi)
		{
			angle += 2.0f * b2_pi;
			C = angle - m_as[i];
		}

		float32 impulse = - m_k3 * mass * C;

		p1 += (m1 * impulse) * J1;
		p2 += (m2 * impulse) * J2;
		p3 += (m3 * impulse) * J3;

		m_ps[i] = p1;
		m_ps[i + 1] = p2;
		m_ps[i + 2] = p3;
	}
}

void b2Rope::Draw(b2Draw* draw) const
{
	b2Color c(0.4f, 0.5f, 0.7f);

	for (int32 i = 0; i < m_count - 1; ++i)
	{
		draw->DrawSegment(m_ps[i], m_ps[i+1], c);
	}
}

// end of Rope.cpp

#endif
