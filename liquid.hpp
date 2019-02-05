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

void b2CollideCircles(
	b2Manifold* manifold,
	const b2CircleShape* circleA, const b2Transform& xfA,
	const b2CircleShape* circleB, const b2Transform& xfB)
{
	manifold->pointCount = 0;

	b2Vec2 pA = b2Mul(xfA, circleA->m_p);
	b2Vec2 pB = b2Mul(xfB, circleB->m_p);

	b2Vec2 d = pB - pA;
	float32 distSqr = b2Dot(d, d);
	float32 rA = circleA->m_radius, rB = circleB->m_radius;
	float32 radius = rA + rB;
	if (distSqr > radius * radius)
	{
		return;
	}

	manifold->type = b2Manifold::e_circles;
	manifold->localPoint = circleA->m_p;
	manifold->localNormal.SetZero();
	manifold->pointCount = 1;

	manifold->points[0].localPoint = circleB->m_p;
	manifold->points[0].id.key = 0;
}

void b2CollidePolygonAndCircle(
	b2Manifold* manifold,
	const b2PolygonShape* polygonA, const b2Transform& xfA,
	const b2CircleShape* circleB, const b2Transform& xfB)
{
	manifold->pointCount = 0;

	// Compute circle position in the frame of the polygon.
	b2Vec2 c = b2Mul(xfB, circleB->m_p);
	b2Vec2 cLocal = b2MulT(xfA, c);

	// Find the min separating edge.
	int32 normalIndex = 0;
	float32 separation = -b2_maxFloat;
	float32 radius = polygonA->m_radius + circleB->m_radius;
	int32 vertexCount = polygonA->m_count;
	const b2Vec2* vertices = polygonA->m_vertices;
	const b2Vec2* normals = polygonA->m_normals;

	for (int32 i = 0; i < vertexCount; ++i)
	{
		float32 s = b2Dot(normals[i], cLocal - vertices[i]);

		if (s > radius)
		{
			// Early out.
			return;
		}

		if (s > separation)
		{
			separation = s;
			normalIndex = i;
		}
	}

	// Vertices that subtend the incident face.
	int32 vertIndex1 = normalIndex;
	int32 vertIndex2 = vertIndex1 + 1 < vertexCount ? vertIndex1 + 1 : 0;
	b2Vec2 v1 = vertices[vertIndex1];
	b2Vec2 v2 = vertices[vertIndex2];

	// If the center is inside the polygon ...
	if (separation < b2_epsilon)
	{
		manifold->pointCount = 1;
		manifold->type = b2Manifold::e_faceA;
		manifold->localNormal = normals[normalIndex];
		manifold->localPoint = 0.5f * (v1 + v2);
		manifold->points[0].localPoint = circleB->m_p;
		manifold->points[0].id.key = 0;
		return;
	}

	// Compute barycentric coordinates
	float32 u1 = b2Dot(cLocal - v1, v2 - v1);
	float32 u2 = b2Dot(cLocal - v2, v1 - v2);
	if (u1 <= 0.0f)
	{
		if (b2DistanceSquared(cLocal, v1) > radius * radius)
		{
			return;
		}

		manifold->pointCount = 1;
		manifold->type = b2Manifold::e_faceA;
		manifold->localNormal = cLocal - v1;
		manifold->localNormal.Normalize();
		manifold->localPoint = v1;
		manifold->points[0].localPoint = circleB->m_p;
		manifold->points[0].id.key = 0;
	}
	else if (u2 <= 0.0f)
	{
		if (b2DistanceSquared(cLocal, v2) > radius * radius)
		{
			return;
		}

		manifold->pointCount = 1;
		manifold->type = b2Manifold::e_faceA;
		manifold->localNormal = cLocal - v2;
		manifold->localNormal.Normalize();
		manifold->localPoint = v2;
		manifold->points[0].localPoint = circleB->m_p;
		manifold->points[0].id.key = 0;
	}
	else
	{
		b2Vec2 faceCenter = 0.5f * (v1 + v2);
		float32 separation = b2Dot(cLocal - faceCenter, normals[vertIndex1]);
		if (separation > radius)
		{
			return;
		}

		manifold->pointCount = 1;
		manifold->type = b2Manifold::e_faceA;
		manifold->localNormal = normals[vertIndex1];
		manifold->localPoint = faceCenter;
		manifold->points[0].localPoint = circleB->m_p;
		manifold->points[0].id.key = 0;
	}
}

// end of CollideCircle.cpp

// Compute contact points for edge versus circle.
// This accounts for edge connectivity.
void b2CollideEdgeAndCircle(b2Manifold* manifold,
							const b2EdgeShape* edgeA, const b2Transform& xfA,
							const b2CircleShape* circleB, const b2Transform& xfB)
{
	manifold->pointCount = 0;
	
	// Compute circle in frame of edge
	b2Vec2 Q = b2MulT(xfA, b2Mul(xfB, circleB->m_p));
	
	b2Vec2 A = edgeA->m_vertex1, B = edgeA->m_vertex2;
	b2Vec2 e = B - A;
	
	// Barycentric coordinates
	float32 u = b2Dot(e, B - Q);
	float32 v = b2Dot(e, Q - A);
	
	float32 radius = edgeA->m_radius + circleB->m_radius;
	
	b2ContactFeature cf;
	cf.indexB = 0;
	cf.typeB = b2ContactFeature::e_vertex;
	
	// Region A
	if (v <= 0.0f)
	{
		b2Vec2 P = A;
		b2Vec2 d = Q - P;
		float32 dd = b2Dot(d, d);
		if (dd > radius * radius)
		{
			return;
		}
		
		// Is there an edge connected to A?
		if (edgeA->m_hasVertex0)
		{
			b2Vec2 A1 = edgeA->m_vertex0;
			b2Vec2 B1 = A;
			b2Vec2 e1 = B1 - A1;
			float32 u1 = b2Dot(e1, B1 - Q);
			
			// Is the circle in Region AB of the previous edge?
			if (u1 > 0.0f)
			{
				return;
			}
		}
		
		cf.indexA = 0;
		cf.typeA = b2ContactFeature::e_vertex;
		manifold->pointCount = 1;
		manifold->type = b2Manifold::e_circles;
		manifold->localNormal.SetZero();
		manifold->localPoint = P;
		manifold->points[0].id.key = 0;
		manifold->points[0].id.cf = cf;
		manifold->points[0].localPoint = circleB->m_p;
		return;
	}
	
	// Region B
	if (u <= 0.0f)
	{
		b2Vec2 P = B;
		b2Vec2 d = Q - P;
		float32 dd = b2Dot(d, d);
		if (dd > radius * radius)
		{
			return;
		}
		
		// Is there an edge connected to B?
		if (edgeA->m_hasVertex3)
		{
			b2Vec2 B2 = edgeA->m_vertex3;
			b2Vec2 A2 = B;
			b2Vec2 e2 = B2 - A2;
			float32 v2 = b2Dot(e2, Q - A2);
			
			// Is the circle in Region AB of the next edge?
			if (v2 > 0.0f)
			{
				return;
			}
		}
		
		cf.indexA = 1;
		cf.typeA = b2ContactFeature::e_vertex;
		manifold->pointCount = 1;
		manifold->type = b2Manifold::e_circles;
		manifold->localNormal.SetZero();
		manifold->localPoint = P;
		manifold->points[0].id.key = 0;
		manifold->points[0].id.cf = cf;
		manifold->points[0].localPoint = circleB->m_p;
		return;
	}
	
	// Region AB
	float32 den = b2Dot(e, e);
	b2Assert(den > 0.0f);
	b2Vec2 P = (1.0f / den) * (u * A + v * B);
	b2Vec2 d = Q - P;
	float32 dd = b2Dot(d, d);
	if (dd > radius * radius)
	{
		return;
	}
	
	b2Vec2 n(-e.y, e.x);
	if (b2Dot(n, Q - A) < 0.0f)
	{
		n.Set(-n.x, -n.y);
	}
	n.Normalize();
	
	cf.indexA = 0;
	cf.typeA = b2ContactFeature::e_face;
	manifold->pointCount = 1;
	manifold->type = b2Manifold::e_faceA;
	manifold->localNormal = n;
	manifold->localPoint = A;
	manifold->points[0].id.key = 0;
	manifold->points[0].id.cf = cf;
	manifold->points[0].localPoint = circleB->m_p;
}

// This structure is used to keep track of the best separating axis.
struct b2EPAxis
{
	enum Type
	{
		e_unknown,
		e_edgeA,
		e_edgeB
	};
	
	Type type;
	int32 index;
	float32 separation;
};

// This holds polygon B expressed in frame A.
struct b2TempPolygon
{
	b2Vec2 vertices[b2_maxPolygonVertices];
	b2Vec2 normals[b2_maxPolygonVertices];
	int32 count;
};

// Reference face used for clipping
struct b2ReferenceFace
{
	int32 i1, i2;
	
	b2Vec2 v1, v2;
	
	b2Vec2 normal;
	
	b2Vec2 sideNormal1;
	float32 sideOffset1;
	
	b2Vec2 sideNormal2;
	float32 sideOffset2;
};

// This class collides and edge and a polygon, taking into account edge adjacency.
struct b2EPCollider
{
	void Collide(b2Manifold* manifold, const b2EdgeShape* edgeA, const b2Transform& xfA,
				 const b2PolygonShape* polygonB, const b2Transform& xfB);
	b2EPAxis ComputeEdgeSeparation();
	b2EPAxis ComputePolygonSeparation();
	
	enum VertexType
	{
		e_isolated,
		e_concave,
		e_convex
	};
	
	b2TempPolygon m_polygonB;
	
	b2Transform m_xf;
	b2Vec2 m_centroidB;
	b2Vec2 m_v0, m_v1, m_v2, m_v3;
	b2Vec2 m_normal0, m_normal1, m_normal2;
	b2Vec2 m_normal;
	VertexType m_type1, m_type2;
	b2Vec2 m_lowerLimit, m_upperLimit;
	float32 m_radius;
	bool m_front;
};

// Algorithm:
// 1. Classify v1 and v2
// 2. Classify polygon centroid as front or back
// 3. Flip normal if necessary
// 4. Initialize normal range to [-pi, pi] about face normal
// 5. Adjust normal range according to adjacent edges
// 6. Visit each separating axes, only accept axes within the range
// 7. Return if _any_ axis indicates separation
// 8. Clip
void b2EPCollider::Collide(b2Manifold* manifold, const b2EdgeShape* edgeA, const b2Transform& xfA,
						   const b2PolygonShape* polygonB, const b2Transform& xfB)
{
	m_xf = b2MulT(xfA, xfB);
	
	m_centroidB = b2Mul(m_xf, polygonB->m_centroid);
	
	m_v0 = edgeA->m_vertex0;
	m_v1 = edgeA->m_vertex1;
	m_v2 = edgeA->m_vertex2;
	m_v3 = edgeA->m_vertex3;
	
	bool hasVertex0 = edgeA->m_hasVertex0;
	bool hasVertex3 = edgeA->m_hasVertex3;
	
	b2Vec2 edge1 = m_v2 - m_v1;
	edge1.Normalize();
	m_normal1.Set(edge1.y, -edge1.x);
	float32 offset1 = b2Dot(m_normal1, m_centroidB - m_v1);
	float32 offset0 = 0.0f, offset2 = 0.0f;
	bool convex1 = false, convex2 = false;
	
	// Is there a preceding edge?
	if (hasVertex0)
	{
		b2Vec2 edge0 = m_v1 - m_v0;
		edge0.Normalize();
		m_normal0.Set(edge0.y, -edge0.x);
		convex1 = b2Cross(edge0, edge1) >= 0.0f;
		offset0 = b2Dot(m_normal0, m_centroidB - m_v0);
	}
	
	// Is there a following edge?
	if (hasVertex3)
	{
		b2Vec2 edge2 = m_v3 - m_v2;
		edge2.Normalize();
		m_normal2.Set(edge2.y, -edge2.x);
		convex2 = b2Cross(edge1, edge2) > 0.0f;
		offset2 = b2Dot(m_normal2, m_centroidB - m_v2);
	}
	
	// Determine front or back collision. Determine collision normal limits.
	if (hasVertex0 && hasVertex3)
	{
		if (convex1 && convex2)
		{
			m_front = offset0 >= 0.0f || offset1 >= 0.0f || offset2 >= 0.0f;
			if (m_front)
			{
				m_normal = m_normal1;
				m_lowerLimit = m_normal0;
				m_upperLimit = m_normal2;
			}
			else
			{
				m_normal = -m_normal1;
				m_lowerLimit = -m_normal1;
				m_upperLimit = -m_normal1;
			}
		}
		else if (convex1)
		{
			m_front = offset0 >= 0.0f || (offset1 >= 0.0f && offset2 >= 0.0f);
			if (m_front)
			{
				m_normal = m_normal1;
				m_lowerLimit = m_normal0;
				m_upperLimit = m_normal1;
			}
			else
			{
				m_normal = -m_normal1;
				m_lowerLimit = -m_normal2;
				m_upperLimit = -m_normal1;
			}
		}
		else if (convex2)
		{
			m_front = offset2 >= 0.0f || (offset0 >= 0.0f && offset1 >= 0.0f);
			if (m_front)
			{
				m_normal = m_normal1;
				m_lowerLimit = m_normal1;
				m_upperLimit = m_normal2;
			}
			else
			{
				m_normal = -m_normal1;
				m_lowerLimit = -m_normal1;
				m_upperLimit = -m_normal0;
			}
		}
		else
		{
			m_front = offset0 >= 0.0f && offset1 >= 0.0f && offset2 >= 0.0f;
			if (m_front)
			{
				m_normal = m_normal1;
				m_lowerLimit = m_normal1;
				m_upperLimit = m_normal1;
			}
			else
			{
				m_normal = -m_normal1;
				m_lowerLimit = -m_normal2;
				m_upperLimit = -m_normal0;
			}
		}
	}
	else if (hasVertex0)
	{
		if (convex1)
		{
			m_front = offset0 >= 0.0f || offset1 >= 0.0f;
			if (m_front)
			{
				m_normal = m_normal1;
				m_lowerLimit = m_normal0;
				m_upperLimit = -m_normal1;
			}
			else
			{
				m_normal = -m_normal1;
				m_lowerLimit = m_normal1;
				m_upperLimit = -m_normal1;
			}
		}
		else
		{
			m_front = offset0 >= 0.0f && offset1 >= 0.0f;
			if (m_front)
			{
				m_normal = m_normal1;
				m_lowerLimit = m_normal1;
				m_upperLimit = -m_normal1;
			}
			else
			{
				m_normal = -m_normal1;
				m_lowerLimit = m_normal1;
				m_upperLimit = -m_normal0;
			}
		}
	}
	else if (hasVertex3)
	{
		if (convex2)
		{
			m_front = offset1 >= 0.0f || offset2 >= 0.0f;
			if (m_front)
			{
				m_normal = m_normal1;
				m_lowerLimit = -m_normal1;
				m_upperLimit = m_normal2;
			}
			else
			{
				m_normal = -m_normal1;
				m_lowerLimit = -m_normal1;
				m_upperLimit = m_normal1;
			}
		}
		else
		{
			m_front = offset1 >= 0.0f && offset2 >= 0.0f;
			if (m_front)
			{
				m_normal = m_normal1;
				m_lowerLimit = -m_normal1;
				m_upperLimit = m_normal1;
			}
			else
			{
				m_normal = -m_normal1;
				m_lowerLimit = -m_normal2;
				m_upperLimit = m_normal1;
			}
		}		
	}
	else
	{
		m_front = offset1 >= 0.0f;
		if (m_front)
		{
			m_normal = m_normal1;
			m_lowerLimit = -m_normal1;
			m_upperLimit = -m_normal1;
		}
		else
		{
			m_normal = -m_normal1;
			m_lowerLimit = m_normal1;
			m_upperLimit = m_normal1;
		}
	}
	
	// Get polygonB in frameA
	m_polygonB.count = polygonB->m_count;
	for (int32 i = 0; i < polygonB->m_count; ++i)
	{
		m_polygonB.vertices[i] = b2Mul(m_xf, polygonB->m_vertices[i]);
		m_polygonB.normals[i] = b2Mul(m_xf.q, polygonB->m_normals[i]);
	}
	
	m_radius = 2.0f * b2_polygonRadius;
	
	manifold->pointCount = 0;
	
	b2EPAxis edgeAxis = ComputeEdgeSeparation();
	
	// If no valid normal can be found than this edge should not collide.
	if (edgeAxis.type == b2EPAxis::e_unknown)
	{
		return;
	}
	
	if (edgeAxis.separation > m_radius)
	{
		return;
	}
	
	b2EPAxis polygonAxis = ComputePolygonSeparation();
	if (polygonAxis.type != b2EPAxis::e_unknown && polygonAxis.separation > m_radius)
	{
		return;
	}
	
	// Use hysteresis for jitter reduction.
	const float32 k_relativeTol = 0.98f;
	const float32 k_absoluteTol = 0.001f;
	
	b2EPAxis primaryAxis;
	if (polygonAxis.type == b2EPAxis::e_unknown)
	{
		primaryAxis = edgeAxis;
	}
	else if (polygonAxis.separation > k_relativeTol * edgeAxis.separation + k_absoluteTol)
	{
		primaryAxis = polygonAxis;
	}
	else
	{
		primaryAxis = edgeAxis;
	}
	
	b2ClipVertex ie[2];
	b2ReferenceFace rf;
	if (primaryAxis.type == b2EPAxis::e_edgeA)
	{
		manifold->type = b2Manifold::e_faceA;
		
		// Search for the polygon normal that is most anti-parallel to the edge normal.
		int32 bestIndex = 0;
		float32 bestValue = b2Dot(m_normal, m_polygonB.normals[0]);
		for (int32 i = 1; i < m_polygonB.count; ++i)
		{
			float32 value = b2Dot(m_normal, m_polygonB.normals[i]);
			if (value < bestValue)
			{
				bestValue = value;
				bestIndex = i;
			}
		}
		
		int32 i1 = bestIndex;
		int32 i2 = i1 + 1 < m_polygonB.count ? i1 + 1 : 0;
		
		ie[0].v = m_polygonB.vertices[i1];
		ie[0].id.cf.indexA = 0;
		ie[0].id.cf.indexB = static_cast<uint8>(i1);
		ie[0].id.cf.typeA = b2ContactFeature::e_face;
		ie[0].id.cf.typeB = b2ContactFeature::e_vertex;
		
		ie[1].v = m_polygonB.vertices[i2];
		ie[1].id.cf.indexA = 0;
		ie[1].id.cf.indexB = static_cast<uint8>(i2);
		ie[1].id.cf.typeA = b2ContactFeature::e_face;
		ie[1].id.cf.typeB = b2ContactFeature::e_vertex;
		
		if (m_front)
		{
			rf.i1 = 0;
			rf.i2 = 1;
			rf.v1 = m_v1;
			rf.v2 = m_v2;
			rf.normal = m_normal1;
		}
		else
		{
			rf.i1 = 1;
			rf.i2 = 0;
			rf.v1 = m_v2;
			rf.v2 = m_v1;
			rf.normal = -m_normal1;
		}		
	}
	else
	{
		manifold->type = b2Manifold::e_faceB;
		
		ie[0].v = m_v1;
		ie[0].id.cf.indexA = 0;
		ie[0].id.cf.indexB = static_cast<uint8>(primaryAxis.index);
		ie[0].id.cf.typeA = b2ContactFeature::e_vertex;
		ie[0].id.cf.typeB = b2ContactFeature::e_face;
		
		ie[1].v = m_v2;
		ie[1].id.cf.indexA = 0;
		ie[1].id.cf.indexB = static_cast<uint8>(primaryAxis.index);		
		ie[1].id.cf.typeA = b2ContactFeature::e_vertex;
		ie[1].id.cf.typeB = b2ContactFeature::e_face;
		
		rf.i1 = primaryAxis.index;
		rf.i2 = rf.i1 + 1 < m_polygonB.count ? rf.i1 + 1 : 0;
		rf.v1 = m_polygonB.vertices[rf.i1];
		rf.v2 = m_polygonB.vertices[rf.i2];
		rf.normal = m_polygonB.normals[rf.i1];
	}
	
	rf.sideNormal1.Set(rf.normal.y, -rf.normal.x);
	rf.sideNormal2 = -rf.sideNormal1;
	rf.sideOffset1 = b2Dot(rf.sideNormal1, rf.v1);
	rf.sideOffset2 = b2Dot(rf.sideNormal2, rf.v2);
	
	// Clip incident edge against extruded edge1 side edges.
	b2ClipVertex clipPoints1[2];
	b2ClipVertex clipPoints2[2];
	int32 np;
	
	// Clip to box side 1
	np = b2ClipSegmentToLine(clipPoints1, ie, rf.sideNormal1, rf.sideOffset1, rf.i1);
	
	if (np < b2_maxManifoldPoints)
	{
		return;
	}
	
	// Clip to negative box side 1
	np = b2ClipSegmentToLine(clipPoints2, clipPoints1, rf.sideNormal2, rf.sideOffset2, rf.i2);
	
	if (np < b2_maxManifoldPoints)
	{
		return;
	}
	
	// Now clipPoints2 contains the clipped points.
	if (primaryAxis.type == b2EPAxis::e_edgeA)
	{
		manifold->localNormal = rf.normal;
		manifold->localPoint = rf.v1;
	}
	else
	{
		manifold->localNormal = polygonB->m_normals[rf.i1];
		manifold->localPoint = polygonB->m_vertices[rf.i1];
	}
	
	int32 pointCount = 0;
	for (int32 i = 0; i < b2_maxManifoldPoints; ++i)
	{
		float32 separation;
		
		separation = b2Dot(rf.normal, clipPoints2[i].v - rf.v1);
		
		if (separation <= m_radius)
		{
			b2ManifoldPoint* cp = manifold->points + pointCount;
			
			if (primaryAxis.type == b2EPAxis::e_edgeA)
			{
				cp->localPoint = b2MulT(m_xf, clipPoints2[i].v);
				cp->id = clipPoints2[i].id;
			}
			else
			{
				cp->localPoint = clipPoints2[i].v;
				cp->id.cf.typeA = clipPoints2[i].id.cf.typeB;
				cp->id.cf.typeB = clipPoints2[i].id.cf.typeA;
				cp->id.cf.indexA = clipPoints2[i].id.cf.indexB;
				cp->id.cf.indexB = clipPoints2[i].id.cf.indexA;
			}
			
			++pointCount;
		}
	}
	
	manifold->pointCount = pointCount;
}

b2EPAxis b2EPCollider::ComputeEdgeSeparation()
{
	b2EPAxis axis;
	axis.type = b2EPAxis::e_edgeA;
	axis.index = m_front ? 0 : 1;
	axis.separation = FLT_MAX;
	
	for (int32 i = 0; i < m_polygonB.count; ++i)
	{
		float32 s = b2Dot(m_normal, m_polygonB.vertices[i] - m_v1);
		if (s < axis.separation)
		{
			axis.separation = s;
		}
	}
	
	return axis;
}

b2EPAxis b2EPCollider::ComputePolygonSeparation()
{
	b2EPAxis axis;
	axis.type = b2EPAxis::e_unknown;
	axis.index = -1;
	axis.separation = -FLT_MAX;

	b2Vec2 perp(-m_normal.y, m_normal.x);

	for (int32 i = 0; i < m_polygonB.count; ++i)
	{
		b2Vec2 n = -m_polygonB.normals[i];
		
		float32 s1 = b2Dot(n, m_polygonB.vertices[i] - m_v1);
		float32 s2 = b2Dot(n, m_polygonB.vertices[i] - m_v2);
		float32 s = b2Min(s1, s2);
		
		if (s > m_radius)
		{
			// No collision
			axis.type = b2EPAxis::e_edgeB;
			axis.index = i;
			axis.separation = s;
			return axis;
		}
		
		// Adjacency
		if (b2Dot(n, perp) >= 0.0f)
		{
			if (b2Dot(n - m_upperLimit, m_normal) < -b2_angularSlop)
			{
				continue;
			}
		}
		else
		{
			if (b2Dot(n - m_lowerLimit, m_normal) < -b2_angularSlop)
			{
				continue;
			}
		}
		
		if (s > axis.separation)
		{
			axis.type = b2EPAxis::e_edgeB;
			axis.index = i;
			axis.separation = s;
		}
	}
	
	return axis;
}

void b2CollideEdgeAndPolygon(	b2Manifold* manifold,
							 const b2EdgeShape* edgeA, const b2Transform& xfA,
							 const b2PolygonShape* polygonB, const b2Transform& xfB)
{
	b2EPCollider collider;
	collider.Collide(manifold, edgeA, xfA, polygonB, xfB);
}

// end of CollideEdge.cpp


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

class b2Fixture;
class b2Joint;
class b2Contact;
class b2Controller;
class b2World;
struct b2FixtureDef;
struct b2JointEdge;
struct b2ContactEdge;

/// The body type.
/// static: zero mass, zero velocity, may be manually moved
/// kinematic: zero mass, non-zero velocity set by user, moved by solver
/// dynamic: positive mass, non-zero velocity determined by forces, moved by solver
enum b2BodyType
{
	b2_staticBody = 0,
	b2_kinematicBody,
	b2_dynamicBody

	// TODO_ERIN
	//b2_bulletBody,
};

/// A body definition holds all the data needed to construct a rigid body.
/// You can safely re-use body definitions. Shapes are added to a body after construction.
struct b2BodyDef
{
	/// This constructor sets the body definition default values.
	b2BodyDef()
	{
		userData = NULL;
		position.Set(0.0f, 0.0f);
		angle = 0.0f;
		linearVelocity.Set(0.0f, 0.0f);
		angularVelocity = 0.0f;
		linearDamping = 0.0f;
		angularDamping = 0.0f;
		allowSleep = true;
		awake = true;
		fixedRotation = false;
		bullet = false;
		type = b2_staticBody;
		active = true;
		gravityScale = 1.0f;
	}

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
	/// Set position with direct floats.
	void SetPosition(float32 positionX, float32 positionY);
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

	/// The body type: static, kinematic, or dynamic.
	/// Note: if a dynamic body would have zero mass, the mass is set to one.
	b2BodyType type;

	/// The world position of the body. Avoid creating bodies at the origin
	/// since this can lead to many overlapping shapes.
	b2Vec2 position;

	/// The world angle of the body in radians.
	float32 angle;

	/// The linear velocity of the body's origin in world co-ordinates.
	b2Vec2 linearVelocity;

	/// The angular velocity of the body.
	float32 angularVelocity;

	/// Linear damping is use to reduce the linear velocity. The damping parameter
	/// can be larger than 1.0f but the damping effect becomes sensitive to the
	/// time step when the damping parameter is large.
	float32 linearDamping;

	/// Angular damping is use to reduce the angular velocity. The damping parameter
	/// can be larger than 1.0f but the damping effect becomes sensitive to the
	/// time step when the damping parameter is large.
	float32 angularDamping;

	/// Set this flag to false if this body should never fall asleep. Note that
	/// this increases CPU usage.
	bool allowSleep;

	/// Is this body initially awake or sleeping?
	bool awake;

	/// Should this body be prevented from rotating? Useful for characters.
	bool fixedRotation;

	/// Is this a fast moving body that should be prevented from tunneling through
	/// other moving bodies? Note that all bodies are prevented from tunneling through
	/// kinematic and static bodies. This setting is only considered on dynamic bodies.
	/// @warning You should use this flag sparingly since it increases processing time.
	bool bullet;

	/// Does this body start out active?
	bool active;

	/// Use this to store application specific body data.
	void* userData;

	/// Scale the gravity applied to this body.
	float32 gravityScale;
};

/// A rigid body. These are created via b2World::CreateBody.
class b2Body
{
public:
	/// Creates a fixture and attach it to this body. Use this function if you need
	/// to set some fixture parameters, like friction. Otherwise you can create the
	/// fixture directly from a shape.
	/// If the density is non-zero, this function automatically updates the mass of the body.
	/// Contacts are not created until the next time step.
	/// @param def the fixture definition.
	/// @warning This function is locked during callbacks.
	b2Fixture* CreateFixture(const b2FixtureDef* def);

	/// Creates a fixture from a shape and attach it to this body.
	/// This is a convenience function. Use b2FixtureDef if you need to set parameters
	/// like friction, restitution, user data, or filtering.
	/// If the density is non-zero, this function automatically updates the mass of the body.
	/// @param shape the shape to be cloned.
	/// @param density the shape density (set to zero for static bodies).
	/// @warning This function is locked during callbacks.
	b2Fixture* CreateFixture(const b2Shape* shape, float32 density);

	/// Destroy a fixture. This removes the fixture from the broad-phase and
	/// destroys all contacts associated with this fixture. This will
	/// automatically adjust the mass of the body if the body is dynamic and the
	/// fixture has positive density.
	/// All fixtures attached to a body are implicitly destroyed when the body is destroyed.
	/// @param fixture the fixture to be removed.
	/// @warning This function is locked during callbacks.
	void DestroyFixture(b2Fixture* fixture);

	/// Set the position of the body's origin and rotation.
	/// Manipulating a body's transform may cause non-physical behavior.
	/// Note: contacts are updated on the next call to b2World::Step.
	/// @param position the world position of the body's local origin.
	/// @param angle the world rotation in radians.
	void SetTransform(const b2Vec2& position, float32 angle);

	/// Get the body transform for the body's origin.
	/// @return the world transform of the body's origin.
	const b2Transform& GetTransform() const;

	/// Get the world body origin position.
	/// @return the world position of the body's origin.
	const b2Vec2& GetPosition() const;

	/// Get the angle in radians.
	/// @return the current world rotation angle in radians.
	float32 GetAngle() const;

	/// Get the world position of the center of mass.
	const b2Vec2& GetWorldCenter() const;

	/// Get the local position of the center of mass.
	const b2Vec2& GetLocalCenter() const;

	/// Set the linear velocity of the center of mass.
	/// @param v the new linear velocity of the center of mass.
	void SetLinearVelocity(const b2Vec2& v);

	/// Get the linear velocity of the center of mass.
	/// @return the linear velocity of the center of mass.
	const b2Vec2& GetLinearVelocity() const;

	/// Set the angular velocity.
	/// @param omega the new angular velocity in radians/second.
	void SetAngularVelocity(float32 omega);

	/// Get the angular velocity.
	/// @return the angular velocity in radians/second.
	float32 GetAngularVelocity() const;

	/// Apply a force at a world point. If the force is not
	/// applied at the center of mass, it will generate a torque and
	/// affect the angular velocity. This wakes up the body.
	/// @param force the world force vector, usually in Newtons (N).
	/// @param point the world position of the point of application.
	/// @param wake also wake up the body
	void ApplyForce(const b2Vec2& force, const b2Vec2& point, bool wake);

	/// Apply a force to the center of mass. This wakes up the body.
	/// @param force the world force vector, usually in Newtons (N).
	/// @param wake also wake up the body
	void ApplyForceToCenter(const b2Vec2& force, bool wake);

	/// Apply a torque. This affects the angular velocity
	/// without affecting the linear velocity of the center of mass.
	/// This wakes up the body.
	/// @param torque about the z-axis (out of the screen), usually in N-m.
	/// @param wake also wake up the body
	void ApplyTorque(float32 torque, bool wake);

	/// Apply an impulse at a point. This immediately modifies the velocity.
	/// It also modifies the angular velocity if the point of application
	/// is not at the center of mass. This wakes up the body.
	/// @param impulse the world impulse vector, usually in N-seconds or kg-m/s.
	/// @param point the world position of the point of application.
	/// @param wake also wake up the body
	void ApplyLinearImpulse(const b2Vec2& impulse, const b2Vec2& point, bool wake);

	/// Apply an angular impulse.
	/// @param impulse the angular impulse in units of kg*m*m/s
	/// @param wake also wake up the body
	void ApplyAngularImpulse(float32 impulse, bool wake);

	/// Get the total mass of the body.
	/// @return the mass, usually in kilograms (kg).
	float32 GetMass() const;

	/// Get the rotational inertia of the body about the local origin.
	/// @return the rotational inertia, usually in kg-m^2.
	float32 GetInertia() const;

	/// Get the mass data of the body.
	/// @return a struct containing the mass, inertia and center of the body.
	void GetMassData(b2MassData* data) const;

	/// Set the mass properties to override the mass properties of the fixtures.
	/// Note that this changes the center of mass position.
	/// Note that creating or destroying fixtures can also alter the mass.
	/// This function has no effect if the body isn't dynamic.
	/// @param massData the mass properties.
	void SetMassData(const b2MassData* data);

	/// This resets the mass properties to the sum of the mass properties of the fixtures.
	/// This normally does not need to be called unless you called SetMassData to override
	/// the mass and you later want to reset the mass.
	void ResetMassData();

	/// Get the world coordinates of a point given the local coordinates.
	/// @param localPoint a point on the body measured relative the the body's origin.
	/// @return the same point expressed in world coordinates.
	b2Vec2 GetWorldPoint(const b2Vec2& localPoint) const;

	/// Get the world coordinates of a vector given the local coordinates.
	/// @param localVector a vector fixed in the body.
	/// @return the same vector expressed in world coordinates.
	b2Vec2 GetWorldVector(const b2Vec2& localVector) const;

	/// Gets a local point relative to the body's origin given a world point.
	/// @param a point in world coordinates.
	/// @return the corresponding local point relative to the body's origin.
	b2Vec2 GetLocalPoint(const b2Vec2& worldPoint) const;

	/// Gets a local vector given a world vector.
	/// @param a vector in world coordinates.
	/// @return the corresponding local vector.
	b2Vec2 GetLocalVector(const b2Vec2& worldVector) const;

	/// Get the world linear velocity of a world point attached to this body.
	/// @param a point in world coordinates.
	/// @return the world velocity of a point.
	b2Vec2 GetLinearVelocityFromWorldPoint(const b2Vec2& worldPoint) const;

	/// Get the world velocity of a local point.
	/// @param a point in local coordinates.
	/// @return the world velocity of a point.
	b2Vec2 GetLinearVelocityFromLocalPoint(const b2Vec2& localPoint) const;

	/// Get the linear damping of the body.
	float32 GetLinearDamping() const;

	/// Set the linear damping of the body.
	void SetLinearDamping(float32 linearDamping);

	/// Get the angular damping of the body.
	float32 GetAngularDamping() const;

	/// Set the angular damping of the body.
	void SetAngularDamping(float32 angularDamping);

	/// Get the gravity scale of the body.
	float32 GetGravityScale() const;

	/// Set the gravity scale of the body.
	void SetGravityScale(float32 scale);

	/// Set the type of this body. This may alter the mass and velocity.
	void SetType(b2BodyType type);

	/// Get the type of this body.
	b2BodyType GetType() const;

	/// Should this body be treated like a bullet for continuous collision detection?
	void SetBullet(bool flag);

	/// Is this body treated like a bullet for continuous collision detection?
	bool IsBullet() const;

	/// You can disable sleeping on this body. If you disable sleeping, the
	/// body will be woken.
	void SetSleepingAllowed(bool flag);

	/// Is this body allowed to sleep
	bool IsSleepingAllowed() const;

	/// Set the sleep state of the body. A sleeping body has very
	/// low CPU cost.
	/// @param flag set to true to wake the body, false to put it to sleep.
	void SetAwake(bool flag);

	/// Get the sleeping state of this body.
	/// @return true if the body is awake.
	bool IsAwake() const;

	/// Set the active state of the body. An inactive body is not
	/// simulated and cannot be collided with or woken up.
	/// If you pass a flag of true, all fixtures will be added to the
	/// broad-phase.
	/// If you pass a flag of false, all fixtures will be removed from
	/// the broad-phase and all contacts will be destroyed.
	/// Fixtures and joints are otherwise unaffected. You may continue
	/// to create/destroy fixtures and joints on inactive bodies.
	/// Fixtures on an inactive body are implicitly inactive and will
	/// not participate in collisions, ray-casts, or queries.
	/// Joints connected to an inactive body are implicitly inactive.
	/// An inactive body is still owned by a b2World object and remains
	/// in the body list.
	void SetActive(bool flag);

	/// Get the active state of the body.
	bool IsActive() const;

	/// Set this body to have fixed rotation. This causes the mass
	/// to be reset.
	void SetFixedRotation(bool flag);

	/// Does this body have fixed rotation?
	bool IsFixedRotation() const;

	/// Get the list of all fixtures attached to this body.
	b2Fixture* GetFixtureList();
	const b2Fixture* GetFixtureList() const;

	/// Get the list of all joints attached to this body.
	b2JointEdge* GetJointList();
	const b2JointEdge* GetJointList() const;

	/// Get the list of all contacts attached to this body.
	/// @warning this list changes during the time step and you may
	/// miss some collisions if you don't use b2ContactListener.
	b2ContactEdge* GetContactList();
	const b2ContactEdge* GetContactList() const;

	/// Get the next body in the world's body list.
	b2Body* GetNext();
	const b2Body* GetNext() const;

	/// Get the user data pointer that was provided in the body definition.
	void* GetUserData() const;

	/// Set the user data. Use this to store your application specific data.
	void SetUserData(void* data);

	/// Get the parent world of this body.
	b2World* GetWorld();
	const b2World* GetWorld() const;

	/// Dump this body to a log file
	void Dump();

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
public:
	/// Get x-coordinate of position.
	float32 GetPositionX() const { return GetPosition().x; }

	/// Get y-coordinate of position.
	float32 GetPositionY() const { return GetPosition().y; }

	/// Set b2Transform using direct floats.
	void SetTransform(float32 positionX, float32 positionY, float32 angle);
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

private:

	friend class b2World;
	friend class b2Island;
	friend class b2ContactManager;
	friend class b2ContactSolver;
	friend class b2Contact;

	friend class b2DistanceJoint;
	friend class b2FrictionJoint;
	friend class b2GearJoint;
	friend class b2MotorJoint;
	friend class b2MouseJoint;
	friend class b2PrismaticJoint;
	friend class b2PulleyJoint;
	friend class b2RevoluteJoint;
	friend class b2RopeJoint;
	friend class b2WeldJoint;
	friend class b2WheelJoint;

	friend class b2ParticleSystem;
	friend class b2ParticleGroup;

	// m_flags
	enum
	{
		e_islandFlag		= 0x0001,
		e_awakeFlag			= 0x0002,
		e_autoSleepFlag		= 0x0004,
		e_bulletFlag		= 0x0008,
		e_fixedRotationFlag	= 0x0010,
		e_activeFlag		= 0x0020,
		e_toiFlag			= 0x0040
	};

	b2Body(const b2BodyDef* bd, b2World* world);
	~b2Body();

	void SynchronizeFixtures();
	void SynchronizeTransform();

	// This is used to prevent connected bodies from colliding.
	// It may lie, depending on the collideConnected flag.
	bool ShouldCollide(const b2Body* other) const;

	void Advance(float32 t);

	b2BodyType m_type;

	uint16 m_flags;

	int32 m_islandIndex;

	b2Transform m_xf;		// the body origin transform
	b2Transform m_xf0;		// the previous transform for particle simulation
	b2Sweep m_sweep;		// the swept motion for CCD

	b2Vec2 m_linearVelocity;
	float32 m_angularVelocity;

	b2Vec2 m_force;
	float32 m_torque;

	b2World* m_world;
	b2Body* m_prev;
	b2Body* m_next;

	b2Fixture* m_fixtureList;
	int32 m_fixtureCount;

	b2JointEdge* m_jointList;
	b2ContactEdge* m_contactList;

	float32 m_mass, m_invMass;

	// Rotational inertia about the center of mass.
	float32 m_I, m_invI;

	float32 m_linearDamping;
	float32 m_angularDamping;
	float32 m_gravityScale;

	float32 m_sleepTime;

	void* m_userData;
};

inline b2BodyType b2Body::GetType() const
{
	return m_type;
}

inline const b2Transform& b2Body::GetTransform() const
{
	return m_xf;
}

inline const b2Vec2& b2Body::GetPosition() const
{
	return m_xf.p;
}

inline float32 b2Body::GetAngle() const
{
	return m_sweep.a;
}

inline const b2Vec2& b2Body::GetWorldCenter() const
{
	return m_sweep.c;
}

inline const b2Vec2& b2Body::GetLocalCenter() const
{
	return m_sweep.localCenter;
}

inline void b2Body::SetLinearVelocity(const b2Vec2& v)
{
	if (m_type == b2_staticBody)
	{
		return;
	}

	if (b2Dot(v,v) > 0.0f)
	{
		SetAwake(true);
	}

	m_linearVelocity = v;
}

inline const b2Vec2& b2Body::GetLinearVelocity() const
{
	return m_linearVelocity;
}

inline void b2Body::SetAngularVelocity(float32 w)
{
	if (m_type == b2_staticBody)
	{
		return;
	}

	if (w * w > 0.0f)
	{
		SetAwake(true);
	}

	m_angularVelocity = w;
}

inline float32 b2Body::GetAngularVelocity() const
{
	return m_angularVelocity;
}

inline float32 b2Body::GetMass() const
{
	return m_mass;
}

inline float32 b2Body::GetInertia() const
{
	return m_I + m_mass * b2Dot(m_sweep.localCenter, m_sweep.localCenter);
}

inline void b2Body::GetMassData(b2MassData* data) const
{
	data->mass = m_mass;
	data->I = m_I + m_mass * b2Dot(m_sweep.localCenter, m_sweep.localCenter);
	data->center = m_sweep.localCenter;
}

inline b2Vec2 b2Body::GetWorldPoint(const b2Vec2& localPoint) const
{
	return b2Mul(m_xf, localPoint);
}

inline b2Vec2 b2Body::GetWorldVector(const b2Vec2& localVector) const
{
	return b2Mul(m_xf.q, localVector);
}

inline b2Vec2 b2Body::GetLocalPoint(const b2Vec2& worldPoint) const
{
	return b2MulT(m_xf, worldPoint);
}

inline b2Vec2 b2Body::GetLocalVector(const b2Vec2& worldVector) const
{
	return b2MulT(m_xf.q, worldVector);
}

inline b2Vec2 b2Body::GetLinearVelocityFromWorldPoint(const b2Vec2& worldPoint) const
{
	return m_linearVelocity + b2Cross(m_angularVelocity, worldPoint - m_sweep.c);
}

inline b2Vec2 b2Body::GetLinearVelocityFromLocalPoint(const b2Vec2& localPoint) const
{
	return GetLinearVelocityFromWorldPoint(GetWorldPoint(localPoint));
}

inline float32 b2Body::GetLinearDamping() const
{
	return m_linearDamping;
}

inline void b2Body::SetLinearDamping(float32 linearDamping)
{
	m_linearDamping = linearDamping;
}

inline float32 b2Body::GetAngularDamping() const
{
	return m_angularDamping;
}

inline void b2Body::SetAngularDamping(float32 angularDamping)
{
	m_angularDamping = angularDamping;
}

inline float32 b2Body::GetGravityScale() const
{
	return m_gravityScale;
}

inline void b2Body::SetGravityScale(float32 scale)
{
	m_gravityScale = scale;
}

inline void b2Body::SetBullet(bool flag)
{
	if (flag)
	{
		m_flags |= e_bulletFlag;
	}
	else
	{
		m_flags &= ~e_bulletFlag;
	}
}

inline bool b2Body::IsBullet() const
{
	return (m_flags & e_bulletFlag) == e_bulletFlag;
}

inline void b2Body::SetAwake(bool flag)
{
	if (flag)
	{
		if ((m_flags & e_awakeFlag) == 0)
		{
			m_flags |= e_awakeFlag;
			m_sleepTime = 0.0f;
		}
	}
	else
	{
		m_flags &= ~e_awakeFlag;
		m_sleepTime = 0.0f;
		m_linearVelocity.SetZero();
		m_angularVelocity = 0.0f;
		m_force.SetZero();
		m_torque = 0.0f;
	}
}

inline bool b2Body::IsAwake() const
{
	return (m_flags & e_awakeFlag) == e_awakeFlag;
}

inline bool b2Body::IsActive() const
{
	return (m_flags & e_activeFlag) == e_activeFlag;
}

inline bool b2Body::IsFixedRotation() const
{
	return (m_flags & e_fixedRotationFlag) == e_fixedRotationFlag;
}

inline void b2Body::SetSleepingAllowed(bool flag)
{
	if (flag)
	{
		m_flags |= e_autoSleepFlag;
	}
	else
	{
		m_flags &= ~e_autoSleepFlag;
		SetAwake(true);
	}
}

inline bool b2Body::IsSleepingAllowed() const
{
	return (m_flags & e_autoSleepFlag) == e_autoSleepFlag;
}

inline b2Fixture* b2Body::GetFixtureList()
{
	return m_fixtureList;
}

inline const b2Fixture* b2Body::GetFixtureList() const
{
	return m_fixtureList;
}

inline b2JointEdge* b2Body::GetJointList()
{
	return m_jointList;
}

inline const b2JointEdge* b2Body::GetJointList() const
{
	return m_jointList;
}

inline b2ContactEdge* b2Body::GetContactList()
{
	return m_contactList;
}

inline const b2ContactEdge* b2Body::GetContactList() const
{
	return m_contactList;
}

inline b2Body* b2Body::GetNext()
{
	return m_next;
}

inline const b2Body* b2Body::GetNext() const
{
	return m_next;
}

inline void b2Body::SetUserData(void* data)
{
	m_userData = data;
}

inline void* b2Body::GetUserData() const
{
	return m_userData;
}

inline void b2Body::ApplyForce(const b2Vec2& force, const b2Vec2& point, bool wake)
{
	if (m_type != b2_dynamicBody)
	{
		return;
	}

	if (wake && (m_flags & e_awakeFlag) == 0)
	{
		SetAwake(true);
	}

	// Don't accumulate a force if the body is sleeping.
	if (m_flags & e_awakeFlag)
	{
		m_force += force;
		m_torque += b2Cross(point - m_sweep.c, force);
	}
}

inline void b2Body::ApplyForceToCenter(const b2Vec2& force, bool wake)
{
	if (m_type != b2_dynamicBody)
	{
		return;
	}

	if (wake && (m_flags & e_awakeFlag) == 0)
	{
		SetAwake(true);
	}

	// Don't accumulate a force if the body is sleeping
	if (m_flags & e_awakeFlag)
	{
		m_force += force;
	}
}

inline void b2Body::ApplyTorque(float32 torque, bool wake)
{
	if (m_type != b2_dynamicBody)
	{
		return;
	}

	if (wake && (m_flags & e_awakeFlag) == 0)
	{
		SetAwake(true);
	}

	// Don't accumulate a force if the body is sleeping
	if (m_flags & e_awakeFlag)
	{
		m_torque += torque;
	}
}

inline void b2Body::ApplyLinearImpulse(const b2Vec2& impulse, const b2Vec2& point, bool wake)
{
	if (m_type != b2_dynamicBody)
	{
		return;
	}

	if (wake && (m_flags & e_awakeFlag) == 0)
	{
		SetAwake(true);
	}

	// Don't accumulate velocity if the body is sleeping
	if (m_flags & e_awakeFlag)
	{
		m_linearVelocity += m_invMass * impulse;
		m_angularVelocity += m_invI * b2Cross(point - m_sweep.c, impulse);
	}
}

inline void b2Body::ApplyAngularImpulse(float32 impulse, bool wake)
{
	if (m_type != b2_dynamicBody)
	{
		return;
	}

	if (wake && (m_flags & e_awakeFlag) == 0)
	{
		SetAwake(true);
	}

	// Don't accumulate velocity if the body is sleeping
	if (m_flags & e_awakeFlag)
	{
		m_angularVelocity += m_invI * impulse;
	}
}

inline void b2Body::SynchronizeTransform()
{
	m_xf.q.Set(m_sweep.a);
	m_xf.p = m_sweep.c - b2Mul(m_xf.q, m_sweep.localCenter);
}

inline void b2Body::Advance(float32 alpha)
{
	// Advance to the new safe time. This doesn't sync the broad-phase.
	m_sweep.Advance(alpha);
	m_sweep.c = m_sweep.c0;
	m_sweep.a = m_sweep.a0;
	m_xf.q.Set(m_sweep.a);
	m_xf.p = m_sweep.c - b2Mul(m_xf.q, m_sweep.localCenter);
}

inline b2World* b2Body::GetWorld()
{
	return m_world;
}

inline const b2World* b2Body::GetWorld() const
{
	return m_world;
}

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
inline void b2BodyDef::SetPosition(float32 positionX, float32 positionY)
{
	position.Set(positionX, positionY);
}

inline void b2Body::SetTransform(float32 positionX, float32 positionY, float32 angle)
{
	SetTransform(b2Vec2(positionX, positionY), angle);
}
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

// end of Body.h

class b2BlockAllocator;
class b2Body;
class b2BroadPhase;
class b2Fixture;

/// This holds contact filtering data.
struct b2Filter
{
	b2Filter()
	{
		categoryBits = 0x0001;
		maskBits = 0xFFFF;
		groupIndex = 0;
	}

	/// The collision category bits. Normally you would just set one bit.
	uint16 categoryBits;

	/// The collision mask bits. This states the categories that this
	/// shape would accept for collision.
	uint16 maskBits;

	/// Collision groups allow a certain group of objects to never collide (negative)
	/// or always collide (positive). Zero means no collision group. Non-zero group
	/// filtering always wins against the mask bits.
	int16 groupIndex;
};

/// A fixture definition is used to create a fixture. This class defines an
/// abstract fixture definition. You can reuse fixture definitions safely.
struct b2FixtureDef
{
	/// The constructor sets the default fixture definition values.
	b2FixtureDef()
	{
		shape = NULL;
		userData = NULL;
		friction = 0.2f;
		restitution = 0.0f;
		density = 0.0f;
		isSensor = false;
	}

	/// The shape, this must be set. The shape will be cloned, so you
	/// can create the shape on the stack.
	const b2Shape* shape;

	/// Use this to store application specific fixture data.
	void* userData;

	/// The friction coefficient, usually in the range [0,1].
	float32 friction;

	/// The restitution (elasticity) usually in the range [0,1].
	float32 restitution;

	/// The density, usually in kg/m^2.
	float32 density;

	/// A sensor shape collects contact information but never generates a collision
	/// response.
	bool isSensor;

	/// Contact filtering data.
	b2Filter filter;
};

/// This proxy is used internally to connect fixtures to the broad-phase.
struct b2FixtureProxy
{
	b2AABB aabb;
	b2Fixture* fixture;
	int32 childIndex;
	int32 proxyId;
};

/// A fixture is used to attach a shape to a body for collision detection. A fixture
/// inherits its transform from its parent. Fixtures hold additional non-geometric data
/// such as friction, collision filters, etc.
/// Fixtures are created via b2Body::CreateFixture.
/// @warning you cannot reuse fixtures.
class b2Fixture
{
public:
	/// Get the type of the child shape. You can use this to down cast to the concrete shape.
	/// @return the shape type.
	b2Shape::Type GetType() const;

	/// Get the child shape. You can modify the child shape, however you should not change the
	/// number of vertices because this will crash some collision caching mechanisms.
	/// Manipulating the shape may lead to non-physical behavior.
	b2Shape* GetShape();
	const b2Shape* GetShape() const;

	/// Set if this fixture is a sensor.
	void SetSensor(bool sensor);

	/// Is this fixture a sensor (non-solid)?
	/// @return the true if the shape is a sensor.
	bool IsSensor() const;

	/// Set the contact filtering data. This will not update contacts until the next time
	/// step when either parent body is active and awake.
	/// This automatically calls Refilter.
	void SetFilterData(const b2Filter& filter);

	/// Get the contact filtering data.
	const b2Filter& GetFilterData() const;

	/// Call this if you want to establish collision that was previously disabled by b2ContactFilter::ShouldCollide.
	void Refilter();

	/// Get the parent body of this fixture. This is NULL if the fixture is not attached.
	/// @return the parent body.
	b2Body* GetBody();
	const b2Body* GetBody() const;

	/// Get the next fixture in the parent body's fixture list.
	/// @return the next shape.
	b2Fixture* GetNext();
	const b2Fixture* GetNext() const;

	/// Get the user data that was assigned in the fixture definition. Use this to
	/// store your application specific data.
	void* GetUserData() const;

	/// Set the user data. Use this to store your application specific data.
	void SetUserData(void* data);

	/// Test a point for containment in this fixture.
	/// @param p a point in world coordinates.
	bool TestPoint(const b2Vec2& p) const;

	/// Compute the distance from this fixture.
	/// @param p a point in world coordinates.
	void ComputeDistance(const b2Vec2& p, float32* distance, b2Vec2* normal, int32 childIndex) const;

	/// Cast a ray against this shape.
	/// @param output the ray-cast results.
	/// @param input the ray-cast input parameters.
	bool RayCast(b2RayCastOutput* output, const b2RayCastInput& input, int32 childIndex) const;

	/// Get the mass data for this fixture. The mass data is based on the density and
	/// the shape. The rotational inertia is about the shape's origin. This operation
	/// may be expensive.
	void GetMassData(b2MassData* massData) const;

	/// Set the density of this fixture. This will _not_ automatically adjust the mass
	/// of the body. You must call b2Body::ResetMassData to update the body's mass.
	void SetDensity(float32 density);

	/// Get the density of this fixture.
	float32 GetDensity() const;

	/// Get the coefficient of friction.
	float32 GetFriction() const;

	/// Set the coefficient of friction. This will _not_ change the friction of
	/// existing contacts.
	void SetFriction(float32 friction);

	/// Get the coefficient of restitution.
	float32 GetRestitution() const;

	/// Set the coefficient of restitution. This will _not_ change the restitution of
	/// existing contacts.
	void SetRestitution(float32 restitution);

	/// Get the fixture's AABB. This AABB may be enlarge and/or stale.
	/// If you need a more accurate AABB, compute it using the shape and
	/// the body transform.
	const b2AABB& GetAABB(int32 childIndex) const;

	/// Dump this fixture to the log file.
	void Dump(int32 bodyIndex);

protected:

	friend class b2Body;
	friend class b2World;
	friend class b2Contact;
	friend class b2ContactManager;

	b2Fixture();

	// We need separation create/destroy functions from the constructor/destructor because
	// the destructor cannot access the allocator (no destructor arguments allowed by C++).
	void Create(b2BlockAllocator* allocator, b2Body* body, const b2FixtureDef* def);
	void Destroy(b2BlockAllocator* allocator);

	// These support body activation/deactivation.
	void CreateProxies(b2BroadPhase* broadPhase, const b2Transform& xf);
	void DestroyProxies(b2BroadPhase* broadPhase);

	void Synchronize(b2BroadPhase* broadPhase, const b2Transform& xf1, const b2Transform& xf2);

	float32 m_density;

	b2Fixture* m_next;
	b2Body* m_body;

	b2Shape* m_shape;

	float32 m_friction;
	float32 m_restitution;

	b2FixtureProxy* m_proxies;
	int32 m_proxyCount;

	b2Filter m_filter;

	bool m_isSensor;

	void* m_userData;
};

inline b2Shape::Type b2Fixture::GetType() const
{
	return m_shape->GetType();
}

inline b2Shape* b2Fixture::GetShape()
{
	return m_shape;
}

inline const b2Shape* b2Fixture::GetShape() const
{
	return m_shape;
}

inline bool b2Fixture::IsSensor() const
{
	return m_isSensor;
}

inline const b2Filter& b2Fixture::GetFilterData() const
{
	return m_filter;
}

inline void* b2Fixture::GetUserData() const
{
	return m_userData;
}

inline void b2Fixture::SetUserData(void* data)
{
	m_userData = data;
}

inline b2Body* b2Fixture::GetBody()
{
	return m_body;
}

inline const b2Body* b2Fixture::GetBody() const
{
	return m_body;
}

inline b2Fixture* b2Fixture::GetNext()
{
	return m_next;
}

inline const b2Fixture* b2Fixture::GetNext() const
{
	return m_next;
}

inline void b2Fixture::SetDensity(float32 density)
{
	b2Assert(b2IsValid(density) && density >= 0.0f);
	m_density = density;
}

inline float32 b2Fixture::GetDensity() const
{
	return m_density;
}

inline float32 b2Fixture::GetFriction() const
{
	return m_friction;
}

inline void b2Fixture::SetFriction(float32 friction)
{
	m_friction = friction;
}

inline float32 b2Fixture::GetRestitution() const
{
	return m_restitution;
}

inline void b2Fixture::SetRestitution(float32 restitution)
{
	m_restitution = restitution;
}

inline bool b2Fixture::TestPoint(const b2Vec2& p) const
{
	return m_shape->TestPoint(m_body->GetTransform(), p);
}

inline void b2Fixture::ComputeDistance(const b2Vec2& p, float32* d, b2Vec2* n, int32 childIndex) const
{
	m_shape->ComputeDistance(m_body->GetTransform(), p, d, n, childIndex);
}

inline bool b2Fixture::RayCast(b2RayCastOutput* output, const b2RayCastInput& input, int32 childIndex) const
{
	return m_shape->RayCast(output, input, m_body->GetTransform(), childIndex);
}

inline void b2Fixture::GetMassData(b2MassData* massData) const
{
	m_shape->ComputeMass(massData, m_density);
}

inline const b2AABB& b2Fixture::GetAABB(int32 childIndex) const
{
	b2Assert(0 <= childIndex && childIndex < m_proxyCount);
	return m_proxies[childIndex].aabb;
}

// end of Fixture.h

/// Profiling data. Times are in milliseconds.
struct b2Profile
{
	float32 step;
	float32 collide;
	float32 solve;
	float32 solveInit;
	float32 solveVelocity;
	float32 solvePosition;
	float32 broadphase;
	float32 solveTOI;
};

/// This is an internal structure.
struct b2TimeStep
{
	float32 dt;			// time step
	float32 inv_dt;		// inverse time step (0 if dt == 0).
	float32 dtRatio;	// dt * inv_dt0
	int32 velocityIterations;
	int32 positionIterations;
	int32 particleIterations;
	bool warmStarting;
};

/// This is an internal structure.
struct b2Position
{
	b2Vec2 c;
	float32 a;
};

/// This is an internal structure.
struct b2Velocity
{
	b2Vec2 v;
	float32 w;
};

/// Solver Data
struct b2SolverData
{
	b2TimeStep step;
	b2Position* positions;
	b2Velocity* velocities;
};

// end of TimeStep.h

struct b2Vec2;
struct b2Transform;
class b2Fixture;
class b2Body;
class b2Joint;
class b2Contact;
class b2ParticleSystem;
struct b2ContactResult;
struct b2Manifold;
class b2ParticleGroup;
struct b2ParticleBodyContact;
struct b2ParticleContact;

/// Joints and fixtures are destroyed when their associated
/// body is destroyed. Implement this listener so that you
/// may nullify references to these joints and shapes.
class b2DestructionListener
{
public:
	virtual ~b2DestructionListener() {}

	/// Called when any joint is about to be destroyed due
	/// to the destruction of one of its attached bodies.
	virtual void SayGoodbye(b2Joint* joint) = 0;

	/// Called when any fixture is about to be destroyed due
	/// to the destruction of its parent body.
	virtual void SayGoodbye(b2Fixture* fixture) = 0;

	/// Called when any particle group is about to be destroyed.
	virtual void SayGoodbye(b2ParticleGroup* group)
	{
		B2_NOT_USED(group);
	}

	/// Called when a particle is about to be destroyed.
	/// The index can be used in conjunction with
	/// b2ParticleSystem::GetUserDataBuffer() or
	/// b2ParticleSystem::GetParticleHandleFromIndex() to determine which
	/// particle has been destroyed.
	virtual void SayGoodbye(b2ParticleSystem* particleSystem, int32 index)
	{
		B2_NOT_USED(particleSystem);
		B2_NOT_USED(index);
	}
};

/// Implement this class to provide collision filtering. In other words, you can implement
/// this class if you want finer control over contact creation.
class b2ContactFilter
{
public:
	virtual ~b2ContactFilter() {}

	/// Return true if contact calculations should be performed between these two shapes.
	/// @warning for performance reasons this is only called when the AABBs begin to overlap.
	virtual bool ShouldCollide(b2Fixture* fixtureA, b2Fixture* fixtureB);

	/// Return true if contact calculations should be performed between a
	/// fixture and particle.  This is only called if the
	/// b2_fixtureContactListenerParticle flag is set on the particle.
	virtual bool ShouldCollide(b2Fixture* fixture,
							   b2ParticleSystem* particleSystem,
							   int32 particleIndex)
	{
		B2_NOT_USED(fixture);
		B2_NOT_USED(particleIndex);
		B2_NOT_USED(particleSystem);
		return true;
	}

	/// Return true if contact calculations should be performed between two
	/// particles.  This is only called if the
	/// b2_particleContactListenerParticle flag is set on the particle.
	virtual bool ShouldCollide(b2ParticleSystem* particleSystem,
							   int32 particleIndexA, int32 particleIndexB)
	{
		B2_NOT_USED(particleSystem);
		B2_NOT_USED(particleIndexA);
		B2_NOT_USED(particleIndexB);
		return true;
	}
};

/// Contact impulses for reporting. Impulses are used instead of forces because
/// sub-step forces may approach infinity for rigid body collisions. These
/// match up one-to-one with the contact points in b2Manifold.
struct b2ContactImpulse
{
	float32 normalImpulses[b2_maxManifoldPoints];
	float32 tangentImpulses[b2_maxManifoldPoints];
	int32 count;
};

/// Implement this class to get contact information. You can use these results for
/// things like sounds and game logic. You can also get contact results by
/// traversing the contact lists after the time step. However, you might miss
/// some contacts because continuous physics leads to sub-stepping.
/// Additionally you may receive multiple callbacks for the same contact in a
/// single time step.
/// You should strive to make your callbacks efficient because there may be
/// many callbacks per time step.
/// @warning You cannot create/destroy Box2D entities inside these callbacks.
class b2ContactListener
{
public:
	virtual ~b2ContactListener() {}

	/// Called when two fixtures begin to touch.
	virtual void BeginContact(b2Contact* contact) { B2_NOT_USED(contact); }

	/// Called when two fixtures cease to touch.
	virtual void EndContact(b2Contact* contact) { B2_NOT_USED(contact); }

	/// Called when a fixture and particle start touching if the
	/// b2_fixtureContactFilterParticle flag is set on the particle.
	virtual void BeginContact(b2ParticleSystem* particleSystem,
							  b2ParticleBodyContact* particleBodyContact)
	{
		B2_NOT_USED(particleSystem);
		B2_NOT_USED(particleBodyContact);
	}

	/// Called when a fixture and particle stop touching if the
	/// b2_fixtureContactFilterParticle flag is set on the particle.
	virtual void EndContact(b2Fixture* fixture,
							b2ParticleSystem* particleSystem, int32 index)
	{
		B2_NOT_USED(fixture);
		B2_NOT_USED(particleSystem);
		B2_NOT_USED(index);
	}

	/// Called when two particles start touching if
	/// b2_particleContactFilterParticle flag is set on either particle.
	virtual void BeginContact(b2ParticleSystem* particleSystem,
							  b2ParticleContact* particleContact)
	{
		B2_NOT_USED(particleSystem);
		B2_NOT_USED(particleContact);
	}

	/// Called when two particles start touching if
	/// b2_particleContactFilterParticle flag is set on either particle.
	virtual void EndContact(b2ParticleSystem* particleSystem,
							int32 indexA, int32 indexB)
	{
		B2_NOT_USED(particleSystem);
		B2_NOT_USED(indexA);
		B2_NOT_USED(indexB);
	}

	/// This is called after a contact is updated. This allows you to inspect a
	/// contact before it goes to the solver. If you are careful, you can modify the
	/// contact manifold (e.g. disable contact).
	/// A copy of the old manifold is provided so that you can detect changes.
	/// Note: this is called only for awake bodies.
	/// Note: this is called even when the number of contact points is zero.
	/// Note: this is not called for sensors.
	/// Note: if you set the number of contact points to zero, you will not
	/// get an EndContact callback. However, you may get a BeginContact callback
	/// the next step.
	virtual void PreSolve(b2Contact* contact, const b2Manifold* oldManifold)
	{
		B2_NOT_USED(contact);
		B2_NOT_USED(oldManifold);
	}

	/// This lets you inspect a contact after the solver is finished. This is useful
	/// for inspecting impulses.
	/// Note: the contact manifold does not include time of impact impulses, which can be
	/// arbitrarily large if the sub-step is small. Hence the impulse is provided explicitly
	/// in a separate data structure.
	/// Note: this is only called for contacts that are touching, solid, and awake.
	virtual void PostSolve(b2Contact* contact, const b2ContactImpulse* impulse)
	{
		B2_NOT_USED(contact);
		B2_NOT_USED(impulse);
	}
};

/// Callback class for AABB queries.
/// See b2World::Query
class b2QueryCallback
{
public:
	virtual ~b2QueryCallback() {}

	/// Called for each fixture found in the query AABB.
	/// @return false to terminate the query.
	virtual bool ReportFixture(b2Fixture* fixture) = 0;

	/// Called for each particle found in the query AABB.
	/// @return false to terminate the query.
	virtual bool ReportParticle(const b2ParticleSystem* particleSystem,
								int32 index)
	{
		B2_NOT_USED(particleSystem);
		B2_NOT_USED(index);
		return false;
	}

	/// Cull an entire particle system from b2World::QueryAABB. Ignored for
	/// b2ParticleSystem::QueryAABB.
	/// @return true if you want to include particleSystem in the AABB query,
	/// or false to cull particleSystem from the AABB query.
	virtual bool ShouldQueryParticleSystem(
		const b2ParticleSystem* particleSystem)
	{
		B2_NOT_USED(particleSystem);
		return true;
	}
};

/// Callback class for ray casts.
/// See b2World::RayCast
class b2RayCastCallback
{
public:
	virtual ~b2RayCastCallback() {}

	/// Called for each fixture found in the query. You control how the ray cast
	/// proceeds by returning a float:
	/// return -1: ignore this fixture and continue
	/// return 0: terminate the ray cast
	/// return fraction: clip the ray to this point
	/// return 1: don't clip the ray and continue
	/// @param fixture the fixture hit by the ray
	/// @param point the point of initial intersection
	/// @param normal the normal vector at the point of intersection
	/// @return -1 to filter, 0 to terminate, fraction to clip the ray for
	/// closest hit, 1 to continue
	virtual float32 ReportFixture(	b2Fixture* fixture, const b2Vec2& point,
									const b2Vec2& normal, float32 fraction) = 0;

	/// Called for each particle found in the query. You control how the ray
	/// cast proceeds by returning a float:
	/// return <=0: ignore the remaining particles in this particle system
	/// return fraction: ignore particles that are 'fraction' percent farther
	///   along the line from 'point1' to 'point2'. Note that 'point1' and
	///   'point2' are parameters to b2World::RayCast.
	/// @param particleSystem the particle system containing the particle
	/// @param index the index of the particle in particleSystem
	/// @param point the point of intersection bt the ray and the particle
	/// @param normal the normal vector at the point of intersection
	/// @param fraction percent (0.0~1.0) from 'point0' to 'point1' along the
	///   ray. Note that 'point1' and 'point2' are parameters to
	///   b2World::RayCast.
	/// @return <=0 to ignore rest of particle system, fraction to ignore
	/// particles that are farther away.
	virtual float32 ReportParticle(const b2ParticleSystem* particleSystem,
								   int32 index, const b2Vec2& point,
								   const b2Vec2& normal, float32 fraction)
	{
		B2_NOT_USED(particleSystem);
		B2_NOT_USED(index);
		B2_NOT_USED(&point);
		B2_NOT_USED(&normal);
		B2_NOT_USED(fraction);
		return 0;
	}

	/// Cull an entire particle system from b2World::RayCast. Ignored in
	/// b2ParticleSystem::RayCast.
	/// @return true if you want to include particleSystem in the RayCast, or
	/// false to cull particleSystem from the RayCast.
	virtual bool ShouldQueryParticleSystem(
		const b2ParticleSystem* particleSystem)
	{
		B2_NOT_USED(particleSystem);
		return true;
	}
};

// end of WorldCallbacks.h

// Return true if contact calculations should be performed between these two shapes.
// If you implement your own collision filter you may want to build from this implementation.
bool b2ContactFilter::ShouldCollide(b2Fixture* fixtureA, b2Fixture* fixtureB)
{
	const b2Filter& filterA = fixtureA->GetFilterData();
	const b2Filter& filterB = fixtureB->GetFilterData();

	if (filterA.groupIndex == filterB.groupIndex && filterA.groupIndex != 0)
	{
		return filterA.groupIndex > 0;
	}

	bool collide = (filterA.maskBits & filterB.categoryBits) != 0 && (filterA.categoryBits & filterB.maskBits) != 0;
	return collide;
}

// end of WorldCallbacks.cpp

class b2Contact;
class b2ContactFilter;
class b2ContactListener;
class b2BlockAllocator;
class b2ParticleSystem;

// Delegate of b2World.
class b2ContactManager
{
public:
	friend class b2ParticleSystem;

	b2ContactManager();

	// Broad-phase callback.
	void AddPair(void* proxyUserDataA, void* proxyUserDataB);

	void FindNewContacts();

	void Destroy(b2Contact* c);

	void Collide();
            
	b2BroadPhase m_broadPhase;
	b2Contact* m_contactList;
	int32 m_contactCount;
	b2ContactFilter* m_contactFilter;
	b2ContactListener* m_contactListener;
	b2BlockAllocator* m_allocator;
};

// end of ContactManager.h

class b2Body;
class b2Contact;
class b2Fixture;
class b2World;
class b2BlockAllocator;
class b2StackAllocator;
class b2ContactListener;

/// Friction mixing law. The idea is to allow either fixture to drive the restitution to zero.
/// For example, anything slides on ice.
inline float32 b2MixFriction(float32 friction1, float32 friction2)
{
	return b2Sqrt(friction1 * friction2);
}

/// Restitution mixing law. The idea is allow for anything to bounce off an inelastic surface.
/// For example, a superball bounces on anything.
inline float32 b2MixRestitution(float32 restitution1, float32 restitution2)
{
	return restitution1 > restitution2 ? restitution1 : restitution2;
}

typedef b2Contact* b2ContactCreateFcn(	b2Fixture* fixtureA, int32 indexA,
										b2Fixture* fixtureB, int32 indexB,
										b2BlockAllocator* allocator);
typedef void b2ContactDestroyFcn(b2Contact* contact, b2BlockAllocator* allocator);

struct b2ContactRegister
{
	b2ContactCreateFcn* createFcn;
	b2ContactDestroyFcn* destroyFcn;
	bool primary;
};

/// A contact edge is used to connect bodies and contacts together
/// in a contact graph where each body is a node and each contact
/// is an edge. A contact edge belongs to a doubly linked list
/// maintained in each attached body. Each contact has two contact
/// nodes, one for each attached body.
struct b2ContactEdge
{
	b2Body* other;			///< provides quick access to the other body attached.
	b2Contact* contact;		///< the contact
	b2ContactEdge* prev;	///< the previous contact edge in the body's contact list
	b2ContactEdge* next;	///< the next contact edge in the body's contact list
};

/// The class manages contact between two shapes. A contact exists for each overlapping
/// AABB in the broad-phase (except if filtered). Therefore a contact object may exist
/// that has no contact points.
class b2Contact
{
public:

	/// Get the contact manifold. Do not modify the manifold unless you understand the
	/// internals of Box2D.
	b2Manifold* GetManifold();
	const b2Manifold* GetManifold() const;

	/// Get the world manifold.
	void GetWorldManifold(b2WorldManifold* worldManifold) const;

	/// Is this contact touching?
	bool IsTouching() const;

	/// Enable/disable this contact. This can be used inside the pre-solve
	/// contact listener. The contact is only disabled for the current
	/// time step (or sub-step in continuous collisions).
	void SetEnabled(bool flag);

	/// Has this contact been disabled?
	bool IsEnabled() const;

	/// Get the next contact in the world's contact list.
	b2Contact* GetNext();
	const b2Contact* GetNext() const;

	/// Get fixture A in this contact.
	b2Fixture* GetFixtureA();
	const b2Fixture* GetFixtureA() const;

	/// Get the child primitive index for fixture A.
	int32 GetChildIndexA() const;

	/// Get fixture B in this contact.
	b2Fixture* GetFixtureB();
	const b2Fixture* GetFixtureB() const;

	/// Get the child primitive index for fixture B.
	int32 GetChildIndexB() const;

	/// Override the default friction mixture. You can call this in b2ContactListener::PreSolve.
	/// This value persists until set or reset.
	void SetFriction(float32 friction);

	/// Get the friction.
	float32 GetFriction() const;

	/// Reset the friction mixture to the default value.
	void ResetFriction();

	/// Override the default restitution mixture. You can call this in b2ContactListener::PreSolve.
	/// The value persists until you set or reset.
	void SetRestitution(float32 restitution);

	/// Get the restitution.
	float32 GetRestitution() const;

	/// Reset the restitution to the default value.
	void ResetRestitution();

	/// Set the desired tangent speed for a conveyor belt behavior. In meters per second.
	void SetTangentSpeed(float32 speed);

	/// Get the desired tangent speed. In meters per second.
	float32 GetTangentSpeed() const;

	/// Evaluate this contact with your own manifold and transforms.
	virtual void Evaluate(b2Manifold* manifold, const b2Transform& xfA, const b2Transform& xfB) = 0;

protected:
	friend class b2ContactManager;
	friend class b2World;
	friend class b2ContactSolver;
	friend class b2Body;
	friend class b2Fixture;

	// Flags stored in m_flags
	enum
	{
		// Used when crawling contact graph when forming islands.
		e_islandFlag		= 0x0001,

        // Set when the shapes are touching.
		e_touchingFlag		= 0x0002,

		// This contact can be disabled (by user)
		e_enabledFlag		= 0x0004,

		// This contact needs filtering because a fixture filter was changed.
		e_filterFlag		= 0x0008,

		// This bullet contact had a TOI event
		e_bulletHitFlag		= 0x0010,

		// This contact has a valid TOI in m_toi
		e_toiFlag			= 0x0020
	};

	/// Flag this contact for filtering. Filtering will occur the next time step.
	void FlagForFiltering();

	static void AddType(b2ContactCreateFcn* createFcn, b2ContactDestroyFcn* destroyFcn,
						b2Shape::Type typeA, b2Shape::Type typeB);
	static void InitializeRegisters();
	static b2Contact* Create(b2Fixture* fixtureA, int32 indexA, b2Fixture* fixtureB, int32 indexB, b2BlockAllocator* allocator);
	static void Destroy(b2Contact* contact, b2Shape::Type typeA, b2Shape::Type typeB, b2BlockAllocator* allocator);
	static void Destroy(b2Contact* contact, b2BlockAllocator* allocator);

	b2Contact() : m_fixtureA(NULL), m_fixtureB(NULL) {}
	b2Contact(b2Fixture* fixtureA, int32 indexA, b2Fixture* fixtureB, int32 indexB);
	virtual ~b2Contact() {}

	void Update(b2ContactListener* listener);

	static b2ContactRegister s_registers[b2Shape::e_typeCount][b2Shape::e_typeCount];
	static bool s_initialized;

	uint32 m_flags;

	// World pool and list pointers.
	b2Contact* m_prev;
	b2Contact* m_next;

	// Nodes for connecting bodies.
	b2ContactEdge m_nodeA;
	b2ContactEdge m_nodeB;

	b2Fixture* m_fixtureA;
	b2Fixture* m_fixtureB;

	int32 m_indexA;
	int32 m_indexB;

	b2Manifold m_manifold;

	int32 m_toiCount;
	float32 m_toi;

	float32 m_friction;
	float32 m_restitution;

	float32 m_tangentSpeed;
};

inline b2Manifold* b2Contact::GetManifold()
{
	return &m_manifold;
}

inline const b2Manifold* b2Contact::GetManifold() const
{
	return &m_manifold;
}

inline void b2Contact::GetWorldManifold(b2WorldManifold* worldManifold) const
{
	const b2Body* bodyA = m_fixtureA->GetBody();
	const b2Body* bodyB = m_fixtureB->GetBody();
	const b2Shape* shapeA = m_fixtureA->GetShape();
	const b2Shape* shapeB = m_fixtureB->GetShape();

	worldManifold->Initialize(&m_manifold, bodyA->GetTransform(), shapeA->m_radius, bodyB->GetTransform(), shapeB->m_radius);
}

inline void b2Contact::SetEnabled(bool flag)
{
	if (flag)
	{
		m_flags |= e_enabledFlag;
	}
	else
	{
		m_flags &= ~e_enabledFlag;
	}
}

inline bool b2Contact::IsEnabled() const
{
	return (m_flags & e_enabledFlag) == e_enabledFlag;
}

inline bool b2Contact::IsTouching() const
{
	return (m_flags & e_touchingFlag) == e_touchingFlag;
}

inline b2Contact* b2Contact::GetNext()
{
	return m_next;
}

inline const b2Contact* b2Contact::GetNext() const
{
	return m_next;
}

inline b2Fixture* b2Contact::GetFixtureA()
{
	return m_fixtureA;
}

inline const b2Fixture* b2Contact::GetFixtureA() const
{
	return m_fixtureA;
}

inline b2Fixture* b2Contact::GetFixtureB()
{
	return m_fixtureB;
}

inline int32 b2Contact::GetChildIndexA() const
{
	return m_indexA;
}

inline const b2Fixture* b2Contact::GetFixtureB() const
{
	return m_fixtureB;
}

inline int32 b2Contact::GetChildIndexB() const
{
	return m_indexB;
}

inline void b2Contact::FlagForFiltering()
{
	m_flags |= e_filterFlag;
}

inline void b2Contact::SetFriction(float32 friction)
{
	m_friction = friction;
}

inline float32 b2Contact::GetFriction() const
{
	return m_friction;
}

inline void b2Contact::ResetFriction()
{
	m_friction = b2MixFriction(m_fixtureA->m_friction, m_fixtureB->m_friction);
}

inline void b2Contact::SetRestitution(float32 restitution)
{
	m_restitution = restitution;
}

inline float32 b2Contact::GetRestitution() const
{
	return m_restitution;
}

inline void b2Contact::ResetRestitution()
{
	m_restitution = b2MixRestitution(m_fixtureA->m_restitution, m_fixtureB->m_restitution);
}

inline void b2Contact::SetTangentSpeed(float32 speed)
{
	m_tangentSpeed = speed;
}

inline float32 b2Contact::GetTangentSpeed() const
{
	return m_tangentSpeed;
}

// Contact.h

b2ContactFilter b2_defaultFilter;
b2ContactListener b2_defaultListener;

b2ContactManager::b2ContactManager()
{
	m_contactList = NULL;
	m_contactCount = 0;
	m_contactFilter = &b2_defaultFilter;
	m_contactListener = &b2_defaultListener;
	m_allocator = NULL;
}

// end of Contact.h

void b2ContactManager::Destroy(b2Contact* c)
{
	b2Fixture* fixtureA = c->GetFixtureA();
	b2Fixture* fixtureB = c->GetFixtureB();
	b2Body* bodyA = fixtureA->GetBody();
	b2Body* bodyB = fixtureB->GetBody();

	if (m_contactListener && c->IsTouching())
	{
		m_contactListener->EndContact(c);
	}

	// Remove from the world.
	if (c->m_prev)
	{
		c->m_prev->m_next = c->m_next;
	}

	if (c->m_next)
	{
		c->m_next->m_prev = c->m_prev;
	}

	if (c == m_contactList)
	{
		m_contactList = c->m_next;
	}

	// Remove from body 1
	if (c->m_nodeA.prev)
	{
		c->m_nodeA.prev->next = c->m_nodeA.next;
	}

	if (c->m_nodeA.next)
	{
		c->m_nodeA.next->prev = c->m_nodeA.prev;
	}

	if (&c->m_nodeA == bodyA->m_contactList)
	{
		bodyA->m_contactList = c->m_nodeA.next;
	}

	// Remove from body 2
	if (c->m_nodeB.prev)
	{
		c->m_nodeB.prev->next = c->m_nodeB.next;
	}

	if (c->m_nodeB.next)
	{
		c->m_nodeB.next->prev = c->m_nodeB.prev;
	}

	if (&c->m_nodeB == bodyB->m_contactList)
	{
		bodyB->m_contactList = c->m_nodeB.next;
	}

	// Call the factory.
	b2Contact::Destroy(c, m_allocator);
	--m_contactCount;
}

// This is the top level collision call for the time step. Here
// all the narrow phase collision is processed for the world
// contact list.
void b2ContactManager::Collide()
{
	// Update awake contacts.
	b2Contact* c = m_contactList;
	while (c)
	{
		b2Fixture* fixtureA = c->GetFixtureA();
		b2Fixture* fixtureB = c->GetFixtureB();
		int32 indexA = c->GetChildIndexA();
		int32 indexB = c->GetChildIndexB();
		b2Body* bodyA = fixtureA->GetBody();
		b2Body* bodyB = fixtureB->GetBody();
		 
		// Is this contact flagged for filtering?
		if (c->m_flags & b2Contact::e_filterFlag)
		{
			// Should these bodies collide?
			if (bodyB->ShouldCollide(bodyA) == false)
			{
				b2Contact* cNuke = c;
				c = cNuke->GetNext();
				Destroy(cNuke);
				continue;
			}

			// Check user filtering.
			if (m_contactFilter && m_contactFilter->ShouldCollide(fixtureA, fixtureB) == false)
			{
				b2Contact* cNuke = c;
				c = cNuke->GetNext();
				Destroy(cNuke);
				continue;
			}

			// Clear the filtering flag.
			c->m_flags &= ~b2Contact::e_filterFlag;
		}

		bool activeA = bodyA->IsAwake() && bodyA->m_type != b2_staticBody;
		bool activeB = bodyB->IsAwake() && bodyB->m_type != b2_staticBody;

		// At least one body must be awake and it must be dynamic or kinematic.
		if (activeA == false && activeB == false)
		{
			c = c->GetNext();
			continue;
		}

		int32 proxyIdA = fixtureA->m_proxies[indexA].proxyId;
		int32 proxyIdB = fixtureB->m_proxies[indexB].proxyId;
		bool overlap = m_broadPhase.TestOverlap(proxyIdA, proxyIdB);

		// Here we destroy contacts that cease to overlap in the broad-phase.
		if (overlap == false)
		{
			b2Contact* cNuke = c;
			c = cNuke->GetNext();
			Destroy(cNuke);
			continue;
		}

		// The contact persists.
		c->Update(m_contactListener);
		c = c->GetNext();
	}
}

void b2ContactManager::FindNewContacts()
{
	m_broadPhase.UpdatePairs(this);
}

void b2ContactManager::AddPair(void* proxyUserDataA, void* proxyUserDataB)
{
	b2FixtureProxy* proxyA = (b2FixtureProxy*)proxyUserDataA;
	b2FixtureProxy* proxyB = (b2FixtureProxy*)proxyUserDataB;

	b2Fixture* fixtureA = proxyA->fixture;
	b2Fixture* fixtureB = proxyB->fixture;

	int32 indexA = proxyA->childIndex;
	int32 indexB = proxyB->childIndex;

	b2Body* bodyA = fixtureA->GetBody();
	b2Body* bodyB = fixtureB->GetBody();

	// Are the fixtures on the same body?
	if (bodyA == bodyB)
	{
		return;
	}

	// TODO_ERIN use a hash table to remove a potential bottleneck when both
	// bodies have a lot of contacts.
	// Does a contact already exist?
	b2ContactEdge* edge = bodyB->GetContactList();
	while (edge)
	{
		if (edge->other == bodyA)
		{
			b2Fixture* fA = edge->contact->GetFixtureA();
			b2Fixture* fB = edge->contact->GetFixtureB();
			int32 iA = edge->contact->GetChildIndexA();
			int32 iB = edge->contact->GetChildIndexB();

			if (fA == fixtureA && fB == fixtureB && iA == indexA && iB == indexB)
			{
				// A contact already exists.
				return;
			}

			if (fA == fixtureB && fB == fixtureA && iA == indexB && iB == indexA)
			{
				// A contact already exists.
				return;
			}
		}

		edge = edge->next;
	}

	// Does a joint override collision? Is at least one body dynamic?
	if (bodyB->ShouldCollide(bodyA) == false)
	{
		return;
	}

	// Check user filtering.
	if (m_contactFilter && m_contactFilter->ShouldCollide(fixtureA, fixtureB) == false)
	{
		return;
	}

	// Call the factory.
	b2Contact* c = b2Contact::Create(fixtureA, indexA, fixtureB, indexB, m_allocator);
	if (c == NULL)
	{
		return;
	}

	// Contact creation may swap fixtures.
	fixtureA = c->GetFixtureA();
	fixtureB = c->GetFixtureB();
	bodyA = fixtureA->GetBody();
	bodyB = fixtureB->GetBody();

	// Insert into the world.
	c->m_prev = NULL;
	c->m_next = m_contactList;
	if (m_contactList != NULL)
	{
		m_contactList->m_prev = c;
	}
	m_contactList = c;

	// Connect to island graph.

	// Connect to body A
	c->m_nodeA.contact = c;
	c->m_nodeA.other = bodyB;

	c->m_nodeA.prev = NULL;
	c->m_nodeA.next = bodyA->m_contactList;
	if (bodyA->m_contactList != NULL)
	{
		bodyA->m_contactList->prev = &c->m_nodeA;
	}
	bodyA->m_contactList = &c->m_nodeA;

	// Connect to body B
	c->m_nodeB.contact = c;
	c->m_nodeB.other = bodyA;

	c->m_nodeB.prev = NULL;
	c->m_nodeB.next = bodyB->m_contactList;
	if (bodyB->m_contactList != NULL)
	{
		bodyB->m_contactList->prev = &c->m_nodeB;
	}
	bodyB->m_contactList = &c->m_nodeB;

	// Wake up the bodies
	if (fixtureA->IsSensor() == false && fixtureB->IsSensor() == false)
	{
		bodyA->SetAwake(true);
		bodyB->SetAwake(true);
	}

	++m_contactCount;
}

// end of ContactManager.cpp

#ifdef LIQUIDFUN_UNIT_TESTS
#include <gtest/gtest.h>
#endif // LIQUIDFUN_UNIT_TESTS

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
#include <cstring>
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

class b2World;
class b2Body;
class b2Shape;
class b2ParticleGroup;
class b2BlockAllocator;
class b2StackAllocator;
class b2QueryCallback;
class b2RayCastCallback;
class b2Fixture;
class b2ContactFilter;
class b2ContactListener;
class b2ParticlePairSet;
class FixtureParticleSet;
struct b2ParticleGroupDef;
struct b2Vec2;
struct b2AABB;
struct FindContactInput;
struct FindContactCheck;

struct b2ParticleContact
{
private:
	// 16-bit particle indices consume less memory and thus improve
	// performance. We iterate through m_contactBuffer many times during
	// b2ParticleSystem::Solve, so reducing the amount of data we churn
	// through speeds things up. Also, FindContactsFromChecks_Simd takes
	// advantage of the reduced size for specific optimizations.
	#ifdef B2_USE_16_BIT_PARTICLE_INDICES
		typedef int16 b2ParticleIndex;
	#else
		typedef int32 b2ParticleIndex;
	#endif

	/// Indices of the respective particles making contact.
	b2ParticleIndex indexA, indexB;

	/// Weight of the contact. A value between 0.0f and 1.0f.
	/// 0.0f ==> particles are just barely touching
	/// 1.0f ==> particles are perfectly on top of each other
	float32 weight;

	/// The normalized direction from A to B.
	b2Vec2 normal;

	/// The logical sum of the particle behaviors that have been set.
	/// See the b2ParticleFlag enum.
	uint32 flags;

public:
	void SetIndices(int32 a, int32 b);
	void SetWeight(float32 w) { weight = w; }
	void SetNormal(const b2Vec2& n) { normal = n; }
	void SetFlags(uint32 f) { flags = f; }

	int32 GetIndexA() const { return indexA; }
	int32 GetIndexB() const { return indexB; }
	float32 GetWeight() const { return weight; }
	const b2Vec2& GetNormal() const { return normal; }
	uint32 GetFlags() const { return flags; }

	bool operator==(const b2ParticleContact& rhs) const;
	bool operator!=(const b2ParticleContact& rhs) const { return !operator==(rhs); }
	bool ApproximatelyEqual(const b2ParticleContact& rhs) const;
};

struct b2ParticleBodyContact
{
	/// Index of the particle making contact.
	int32 index;

	/// The body making contact.
	b2Body* body;

	/// The specific fixture making contact
	b2Fixture* fixture;

	/// Weight of the contact. A value between 0.0f and 1.0f.
	float32 weight;

	/// The normalized direction from the particle to the body.
	b2Vec2 normal;

	/// The effective mass used in calculating force.
	float32 mass;
};

/// Connection between two particles
struct b2ParticlePair
{
	/// Indices of the respective particles making pair.
	int32 indexA, indexB;

	/// The logical sum of the particle flags. See the b2ParticleFlag enum.
	uint32 flags;

	/// The strength of cohesion among the particles.
	float32 strength;

	/// The initial distance of the particles.
	float32 distance;
};

/// Connection between three particles
struct b2ParticleTriad
{
	/// Indices of the respective particles making triad.
	int32 indexA, indexB, indexC;

	/// The logical sum of the particle flags. See the b2ParticleFlag enum.
	uint32 flags;

	/// The strength of cohesion among the particles.
	float32 strength;

	/// Values used for calculation.
	b2Vec2 pa, pb, pc;
	float32 ka, kb, kc, s;
};

struct b2ParticleSystemDef
{
	b2ParticleSystemDef()
	{
		strictContactCheck = false;
		density = 1.0f;
		gravityScale = 1.0f;
		radius = 1.0f;
		maxCount = 0;

		// Initialize physical coefficients to the maximum values that
		// maintain numerical stability.
		pressureStrength = 0.05f;
		dampingStrength = 1.0f;
		elasticStrength = 0.25f;
		springStrength = 0.25f;
		viscousStrength = 0.25f;
		surfaceTensionPressureStrength = 0.2f;
		surfaceTensionNormalStrength = 0.2f;
		repulsiveStrength = 1.0f;
		powderStrength = 0.5f;
		ejectionStrength = 0.5f;
		staticPressureStrength = 0.2f;
		staticPressureRelaxation = 0.2f;
		staticPressureIterations = 8;
		colorMixingStrength = 0.5f;
		destroyByAge = true;
		lifetimeGranularity = 1.0f / 60.0f;
	}

	/// Enable strict Particle/Body contact check.
	/// See SetStrictContactCheck for details.
	bool strictContactCheck;

	/// Set the particle density.
	/// See SetDensity for details.
	float32 density;

	/// Change the particle gravity scale. Adjusts the effect of the global
	/// gravity vector on particles. Default value is 1.0f.
	float32 gravityScale;

	/// Particles behave as circles with this radius. In Box2D units.
	float32 radius;

	/// Set the maximum number of particles.
	/// By default, there is no maximum. The particle buffers can continue to
	/// grow while b2World's block allocator still has memory.
	/// See SetMaxParticleCount for details.
	int32 maxCount;

	/// Increases pressure in response to compression
	/// Smaller values allow more compression
	float32 pressureStrength;

	/// Reduces velocity along the collision normal
	/// Smaller value reduces less
	float32 dampingStrength;

	/// Restores shape of elastic particle groups
	/// Larger values increase elastic particle velocity
	float32 elasticStrength;

	/// Restores length of spring particle groups
	/// Larger values increase spring particle velocity
	float32 springStrength;

	/// Reduces relative velocity of viscous particles
	/// Larger values slow down viscous particles more
	float32 viscousStrength;

	/// Produces pressure on tensile particles
	/// 0~0.2. Larger values increase the amount of surface tension.
	float32 surfaceTensionPressureStrength;

	/// Smoothes outline of tensile particles
	/// 0~0.2. Larger values result in rounder, smoother, water-drop-like
	/// clusters of particles.
	float32 surfaceTensionNormalStrength;

	/// Produces additional pressure on repulsive particles
	/// Larger values repulse more
	/// Negative values mean attraction. The range where particles behave
	/// stably is about -0.2 to 2.0.
	float32 repulsiveStrength;

	/// Produces repulsion between powder particles
	/// Larger values repulse more
	float32 powderStrength;

	/// Pushes particles out of solid particle group
	/// Larger values repulse more
	float32 ejectionStrength;

	/// Produces static pressure
	/// Larger values increase the pressure on neighboring partilces
	/// For a description of static pressure, see
	/// http://en.wikipedia.org/wiki/Static_pressure#Static_pressure_in_fluid_dynamics
	float32 staticPressureStrength;

	/// Reduces instability in static pressure calculation
	/// Larger values make stabilize static pressure with fewer iterations
	float32 staticPressureRelaxation;

	/// Computes static pressure more precisely
	/// See SetStaticPressureIterations for details
	int32 staticPressureIterations;

	/// Determines how fast colors are mixed
	/// 1.0f ==> mixed immediately
	/// 0.5f ==> mixed half way each simulation step (see b2World::Step())
	float32 colorMixingStrength;

	/// Whether to destroy particles by age when no more particles can be
	/// created.  See #b2ParticleSystem::SetDestructionByAge() for
	/// more information.
	bool destroyByAge;

	/// Granularity of particle lifetimes in seconds.  By default this is
	/// set to (1.0f / 60.0f) seconds.  b2ParticleSystem uses a 32-bit signed
	/// value to track particle lifetimes so the maximum lifetime of a
	/// particle is (2^32 - 1) / (1.0f / lifetimeGranularity) seconds.
	/// With the value set to 1/60 the maximum lifetime or age of a particle is
	/// 2.27 years.
	float32 lifetimeGranularity;
};


class b2ParticleSystem
{
public:
	/// Create a particle whose properties have been defined.
	/// No reference to the definition is retained.
	/// A simulation step must occur before it's possible to interact with a
	/// newly created particle.  For example, DestroyParticleInShape() will
	/// not destroy a particle until b2World::Step() has been called.
	/// @warning This function is locked during callbacks.
	/// @return the index of the particle.
	int32 CreateParticle(const b2ParticleDef& def);

	/// Retrieve a handle to the particle at the specified index.
	/// Please see #b2ParticleHandle for why you might want a handle.
	const b2ParticleHandle* GetParticleHandleFromIndex(const int32 index);

	/// Destroy a particle.
	/// The particle is removed after the next simulation step (see
	/// b2World::Step()).
	void DestroyParticle(int32 index)
	{
		DestroyParticle(index, false);
	}

	/// Destroy a particle.
	/// The particle is removed after the next step.
	/// @param Index of the particle to destroy.
	/// @param Whether to call the destruction listener just before the
	/// particle is destroyed.
	void DestroyParticle(int32 index, bool callDestructionListener);

	/// Destroy the Nth oldest particle in the system.
	/// The particle is removed after the next b2World::Step().
	/// @param Index of the Nth oldest particle to destroy, 0 will destroy the
	/// oldest particle in the system, 1 will destroy the next oldest
	/// particle etc.
	/// @param Whether to call the destruction listener just before the
	/// particle is destroyed.
	void DestroyOldestParticle(const int32 index,
							   const bool callDestructionListener);

	/// Destroy particles inside a shape without enabling the destruction
	/// callback for destroyed particles.
	/// This function is locked during callbacks.
	/// For more information see
	/// DestroyParticleInShape(const b2Shape&, const b2Transform&,bool).
	/// @param Shape which encloses particles that should be destroyed.
	/// @param Transform applied to the shape.
	/// @warning This function is locked during callbacks.
	/// @return Number of particles destroyed.
	int32 DestroyParticlesInShape(const b2Shape& shape, const b2Transform& xf)
	{
		return DestroyParticlesInShape(shape, xf, false);
	}

	/// Destroy particles inside a shape.
	/// This function is locked during callbacks.
	/// In addition, this function immediately destroys particles in the shape
	/// in constrast to DestroyParticle() which defers the destruction until
	/// the next simulation step.
	/// @param Shape which encloses particles that should be destroyed.
	/// @param Transform applied to the shape.
	/// @param Whether to call the world b2DestructionListener for each
	/// particle destroyed.
	/// @warning This function is locked during callbacks.
	/// @return Number of particles destroyed.
	int32 DestroyParticlesInShape(const b2Shape& shape, const b2Transform& xf,
	                              bool callDestructionListener);

	/// Create a particle group whose properties have been defined. No
	/// reference to the definition is retained.
	/// @warning This function is locked during callbacks.
	b2ParticleGroup* CreateParticleGroup(const b2ParticleGroupDef& def);

	/// Join two particle groups.
	/// @param the first group. Expands to encompass the second group.
	/// @param the second group. It is destroyed.
	/// @warning This function is locked during callbacks.
	void JoinParticleGroups(b2ParticleGroup* groupA, b2ParticleGroup* groupB);

	/// Split particle group into multiple disconnected groups.
	/// @param the group to be split.
	/// @warning This function is locked during callbacks.
	void SplitParticleGroup(b2ParticleGroup* group);

	/// Get the world particle group list. With the returned group, use
	/// b2ParticleGroup::GetNext to get the next group in the world list.
	/// A NULL group indicates the end of the list.
	/// @return the head of the world particle group list.
	b2ParticleGroup* GetParticleGroupList();
	const b2ParticleGroup* GetParticleGroupList() const;

	/// Get the number of particle groups.
	int32 GetParticleGroupCount() const;

	/// Get the number of particles.
	int32 GetParticleCount() const;

	/// Get the maximum number of particles.
	int32 GetMaxParticleCount() const;

	/// Set the maximum number of particles.
	/// A value of 0 means there is no maximum. The particle buffers can
	/// continue to grow while b2World's block allocator still has memory.
	/// Note: If you try to CreateParticle() with more than this count,
	/// b2_invalidParticleIndex is returned unless
	/// SetDestructionByAge() is used to enable the destruction of the
	/// oldest particles in the system.
	void SetMaxParticleCount(int32 count);

	/// Get all existing particle flags.
	uint32 GetAllParticleFlags() const;

	/// Get all existing particle group flags.
	uint32 GetAllGroupFlags() const;

	/// Pause or unpause the particle system. When paused, b2World::Step()
	/// skips over this particle system. All b2ParticleSystem function calls
	/// still work.
	/// @param paused is true to pause, false to un-pause.
	void SetPaused(bool paused);

	/// @return true if the particle system is being updated in
	/// b2World::Step().
	/// Initially, true, then, the last value passed into SetPaused().
	bool GetPaused() const;

	/// Change the particle density.
	/// Particle density affects the mass of the particles, which in turn
	/// affects how the particles interact with b2Bodies. Note that the density
	/// does not affect how the particles interact with each other.
	void SetDensity(float32 density);

	/// Get the particle density.
	float32 GetDensity() const;

	/// Change the particle gravity scale. Adjusts the effect of the global
	/// gravity vector on particles.
	void SetGravityScale(float32 gravityScale);

	/// Get the particle gravity scale.
	float32 GetGravityScale() const;

	/// Damping is used to reduce the velocity of particles. The damping
	/// parameter can be larger than 1.0f but the damping effect becomes
	/// sensitive to the time step when the damping parameter is large.
	void SetDamping(float32 damping);

	/// Get damping for particles
	float32 GetDamping() const;

	/// Change the number of iterations when calculating the static pressure of
	/// particles. By default, 8 iterations. You can reduce the number of
	/// iterations down to 1 in some situations, but this may cause
	/// instabilities when many particles come together. If you see particles
	/// popping away from each other like popcorn, you may have to increase the
	/// number of iterations.
	/// For a description of static pressure, see
	/// http://en.wikipedia.org/wiki/Static_pressure#Static_pressure_in_fluid_dynamics
	void SetStaticPressureIterations(int32 iterations);

	/// Get the number of iterations for static pressure of particles.
	int32 GetStaticPressureIterations() const;

	/// Change the particle radius.
	/// You should set this only once, on world start.
	/// If you change the radius during execution, existing particles may
	/// explode, shrink, or behave unexpectedly.
	void SetRadius(float32 radius);

	/// Get the particle radius.
	float32 GetRadius() const;

	/// Get the position of each particle
	/// Array is length GetParticleCount()
	/// @return the pointer to the head of the particle positions array.
	b2Vec2* GetPositionBuffer();
	const b2Vec2* GetPositionBuffer() const;

	/// Get the velocity of each particle
	/// Array is length GetParticleCount()
	/// @return the pointer to the head of the particle velocities array.
	b2Vec2* GetVelocityBuffer();
	const b2Vec2* GetVelocityBuffer() const;

	/// Get the color of each particle
	/// Array is length GetParticleCount()
	/// @return the pointer to the head of the particle colors array.
	b2ParticleColor* GetColorBuffer();
	const b2ParticleColor* GetColorBuffer() const;

	/// Get the particle-group of each particle.
	/// Array is length GetParticleCount()
	/// @return the pointer to the head of the particle group array.
	b2ParticleGroup* const* GetGroupBuffer();
	const b2ParticleGroup* const* GetGroupBuffer() const;

	/// Get the weight of each particle
	/// Array is length GetParticleCount()
	/// @return the pointer to the head of the particle positions array.
	float32* GetWeightBuffer();
	const float32* GetWeightBuffer() const;

	/// Get the user-specified data of each particle.
	/// Array is length GetParticleCount()
	/// @return the pointer to the head of the particle user-data array.
	void** GetUserDataBuffer();
	void* const* GetUserDataBuffer() const;

	/// Get the flags for each particle. See the b2ParticleFlag enum.
	/// Array is length GetParticleCount()
	/// @return the pointer to the head of the particle-flags array.
	const uint32* GetFlagsBuffer() const;

	/// Set flags for a particle. See the b2ParticleFlag enum.
	void SetParticleFlags(int32 index, uint32 flags);
	/// Get flags for a particle. See the b2ParticleFlag enum.
	uint32 GetParticleFlags(const int32 index);

	/// Set an external buffer for particle data.
	/// Normally, the b2World's block allocator is used for particle data.
	/// However, sometimes you may have an OpenGL or Java buffer for particle
	/// data. To avoid data duplication, you may supply this external buffer.
	///
	/// Note that, when b2World's block allocator is used, the particle data
	/// buffers can grow as required. However, when external buffers are used,
	/// the maximum number of particles is clamped to the size of the smallest
	/// external buffer.
	///
	/// @param buffer is a pointer to a block of memory.
	/// @param size is the number of values in the block.
	void SetFlagsBuffer(uint32* buffer, int32 capacity);
	void SetPositionBuffer(b2Vec2* buffer, int32 capacity);
	void SetVelocityBuffer(b2Vec2* buffer, int32 capacity);
	void SetColorBuffer(b2ParticleColor* buffer, int32 capacity);
	void SetUserDataBuffer(void** buffer, int32 capacity);

	/// Get contacts between particles
	/// Contact data can be used for many reasons, for example to trigger
	/// rendering or audio effects.
	const b2ParticleContact* GetContacts() const;
	int32 GetContactCount() const;

	/// Get contacts between particles and bodies
	/// Contact data can be used for many reasons, for example to trigger
	/// rendering or audio effects.
	const b2ParticleBodyContact* GetBodyContacts() const;
	int32 GetBodyContactCount() const;

	/// Get array of particle pairs. The particles in a pair:
	///   (1) are contacting,
	///   (2) are in the same particle group,
	///   (3) are part of a rigid particle group, or are spring, elastic,
	///       or wall particles.
	///   (4) have at least one particle that is a spring or barrier
	///       particle (i.e. one of the types in k_pairFlags),
	///   (5) have at least one particle that returns true for
	///       ConnectionFilter::IsNecessary,
	///   (6) are not zombie particles.
	/// Essentially, this is an array of spring or barrier particles that
	/// are interacting. The array is sorted by b2ParticlePair's indexA,
	/// and then indexB. There are no duplicate entries.
	const b2ParticlePair* GetPairs() const;
	int32 GetPairCount() const;

	/// Get array of particle triads. The particles in a triad:
	///   (1) are in the same particle group,
	///   (2) are in a Voronoi triangle together,
	///   (3) are within b2_maxTriadDistance particle diameters of each
	///       other,
	///   (4) return true for ConnectionFilter::ShouldCreateTriad
	///   (5) have at least one particle of type elastic (i.e. one of the
	///       types in k_triadFlags),
	///   (6) are part of a rigid particle group, or are spring, elastic,
	///       or wall particles.
	///   (7) are not zombie particles.
	/// Essentially, this is an array of elastic particles that are
	/// interacting. The array is sorted by b2ParticleTriad's indexA,
	/// then indexB, then indexC. There are no duplicate entries.
	const b2ParticleTriad* GetTriads() const;
	int32 GetTriadCount() const;

	/// Set an optional threshold for the maximum number of
	/// consecutive particle iterations that a particle may contact
	/// multiple bodies before it is considered a candidate for being
	/// "stuck". Setting to zero or less disables.
	void SetStuckThreshold(int32 iterations);

	/// Get potentially stuck particles from the last step; the user must
	/// decide if they are stuck or not, and if so, delete or move them
	const int32* GetStuckCandidates() const;

	/// Get the number of stuck particle candidates from the last step.
	int32 GetStuckCandidateCount() const;

	/// Compute the kinetic energy that can be lost by damping force
	float32 ComputeCollisionEnergy() const;

	/// Set strict Particle/Body contact check.
	/// This is an option that will help ensure correct behavior if there are
	/// corners in the world model where Particle/Body contact is ambiguous.
	/// This option scales at n*log(n) of the number of Particle/Body contacts,
	/// so it is best to only enable if it is necessary for your geometry.
	/// Enable if you see strange particle behavior around b2Body
	/// intersections.
	void SetStrictContactCheck(bool enabled);
	/// Get the status of the strict contact check.
	bool GetStrictContactCheck() const;

	/// Set the lifetime (in seconds) of a particle relative to the current
	/// time.  A lifetime of less than or equal to 0.0f results in the particle
	/// living forever until it's manually destroyed by the application.
	void SetParticleLifetime(const int32 index, const float32 lifetime);
	/// Get the lifetime (in seconds) of a particle relative to the current
	/// time.  A value > 0.0f is returned if the particle is scheduled to be
	/// destroyed in the future, values <= 0.0f indicate the particle has an
	/// infinite lifetime.
	float32 GetParticleLifetime(const int32 index);

	/// Enable / disable destruction of particles in CreateParticle() when
	/// no more particles can be created due to a prior call to
	/// SetMaxParticleCount().  When this is enabled, the oldest particle is
	/// destroyed in CreateParticle() favoring the destruction of particles
	/// with a finite lifetime over particles with infinite lifetimes.
	/// This feature is enabled by default when particle lifetimes are
	/// tracked.  Explicitly enabling this feature using this function enables
	/// particle lifetime tracking.
	void SetDestructionByAge(const bool enable);
	/// Get whether the oldest particle will be destroyed in CreateParticle()
	/// when the maximum number of particles are present in the system.
	bool GetDestructionByAge() const;

	/// Get the array of particle expiration times indexed by particle index.
	/// GetParticleCount() items are in the returned array.
	const int32* GetExpirationTimeBuffer();
	/// Convert a expiration time value in returned by
	/// GetExpirationTimeBuffer() to a time in seconds relative to the
	/// current simulation time.
	float32 ExpirationTimeToLifetime(const int32 expirationTime) const;
	/// Get the array of particle indices ordered by reverse lifetime.
	/// The oldest particle indexes are at the end of the array with the
	/// newest at the start.  Particles with infinite lifetimes
	/// (i.e expiration times less than or equal to 0) are placed at the start
	///  of the array.
	/// ExpirationTimeToLifetime(GetExpirationTimeBuffer()[index])
	/// is equivalent to GetParticleLifetime(index).
	/// GetParticleCount() items are in the returned array.
	const int32* GetIndexByExpirationTimeBuffer();

	/// Apply an impulse to one particle. This immediately modifies the
	/// velocity. Similar to b2Body::ApplyLinearImpulse.
	/// @param index the particle that will be modified.
	/// @param impulse the world impulse vector, usually in N-seconds or
    ///        kg-m/s.
	void ParticleApplyLinearImpulse(int32 index, const b2Vec2& impulse);

	/// Apply an impulse to all particles between 'firstIndex' and 'lastIndex'.
	/// This immediately modifies the velocity. Note that the impulse is
	/// applied to the total mass of all particles. So, calling
	/// ParticleApplyLinearImpulse(0, impulse) and
	/// ParticleApplyLinearImpulse(1, impulse) will impart twice as much
	/// velocity as calling just ApplyLinearImpulse(0, 1, impulse).
	/// @param firstIndex the first particle to be modified.
	/// @param lastIndex the last particle to be modified.
	/// @param impulse the world impulse vector, usually in N-seconds or
    ///        kg-m/s.
	void ApplyLinearImpulse(int32 firstIndex, int32 lastIndex,
							const b2Vec2& impulse);

	/// Apply a force to the center of a particle.
	/// @param index the particle that will be modified.
	/// @param force the world force vector, usually in Newtons (N).
	void ParticleApplyForce(int32 index, const b2Vec2& force);

	/// Distribute a force across several particles. The particles must not be
	/// wall particles. Note that the force is distributed across all the
	/// particles, so calling this function for indices 0..N is not the same as
	/// calling ParticleApplyForce(i, force) for i in 0..N.
	/// @param firstIndex the first particle to be modified.
	/// @param lastIndex the last particle to be modified.
	/// @param force the world force vector, usually in Newtons (N).
	void ApplyForce(int32 firstIndex, int32 lastIndex, const b2Vec2& force);

	/// Get the next particle-system in the world's particle-system list.
	b2ParticleSystem* GetNext();
	const b2ParticleSystem* GetNext() const;

	/// Query the particle system for all particles that potentially overlap
	/// the provided AABB. b2QueryCallback::ShouldQueryParticleSystem is
	/// ignored.
	/// @param callback a user implemented callback class.
	/// @param aabb the query box.
	void QueryAABB(b2QueryCallback* callback, const b2AABB& aabb) const;

	/// Query the particle system for all particles that potentially overlap
	/// the provided shape's AABB. Calls QueryAABB internally.
	/// b2QueryCallback::ShouldQueryParticleSystem is ignored.
	/// @param callback a user implemented callback class.
	/// @param shape the query shape
	/// @param xf the transform of the AABB
	void QueryShapeAABB(b2QueryCallback* callback, const b2Shape& shape,
						const b2Transform& xf) const;

	/// Ray-cast the particle system for all particles in the path of the ray.
	/// Your callback controls whether you get the closest point, any point, or
	/// n-points. The ray-cast ignores particles that contain the starting
	/// point. b2RayCastCallback::ShouldQueryParticleSystem is ignored.
	/// @param callback a user implemented callback class.
	/// @param point1 the ray starting point
	/// @param point2 the ray ending point
	void RayCast(b2RayCastCallback* callback, const b2Vec2& point1,
				 const b2Vec2& point2) const;

	/// Compute the axis-aligned bounding box for all particles contained
	/// within this particle system.
	/// @param aabb Returns the axis-aligned bounding box of the system.
	void ComputeAABB(b2AABB* const aabb) const;

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
public:
	enum b2ExceptionType
	{
		b2_bufferTooSmall,
		b2_particleIndexOutOfBounds,
		b2_numErrors,
		b2_noExceptions,
	};

	/// Set the velocity of particle at index with direct floats.
	void SetParticleVelocity(int32 index, float32 vx, float32 vy);

	/// Get the x-coordinate of particle at index.
	float GetParticlePositionX(int32 index) const;

	/// Get the y-coordinate of particle at index.
	float GetParticlePositionY(int32 index) const;

	/// Copy position buffer into a specified buffer, starting from startIndex.
	int CopyPositionBuffer(int startIndex, int numParticles, void* outBuf,
						   int size) const;

	/// Copy color buffer into a specified buffer, starting from startIndex.
	int CopyColorBuffer(int startIndex, int numParticles, void* outBuf,
					    int size) const;

	/// Copy color buffer into a specified buffer, starting from startIndex.
	int CopyWeightBuffer(int startIndex, int numParticles, void* outBuf,
						 int size) const;

private:
	/// Helper function for buffer copies.
	int CopyBuffer(int startIndex, int numParticles, void* inBufWithOffset,
				   void* outBuf, int outBufSize, int copySize) const;

	/// Check if buffer copy is valid for the Get*Buffer functions that have
	/// a user-supplied output buffer.
	b2ExceptionType IsBufCopyValid(int startIndex, int numParticles,
								   int copySize, int bufSize) const;
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

private:
	friend class b2World;
	friend class b2ParticleGroup;
	friend class b2ParticleBodyContactRemovePredicate;
	friend class b2FixtureParticleQueryCallback;
#ifdef LIQUIDFUN_UNIT_TESTS
	FRIEND_TEST(FunctionTests, GetParticleMass);
	FRIEND_TEST(FunctionTests, AreProxyBuffersTheSame);
#endif // LIQUIDFUN_UNIT_TESTS

	template <typename T>
	struct UserOverridableBuffer
	{
		UserOverridableBuffer()
		{
			data = NULL;
			userSuppliedCapacity = 0;
		}
		T* data;
		int32 userSuppliedCapacity;
	};

	/// Used for detecting particle contacts
	struct Proxy
	{
		int32 index;
		uint32 tag;
		friend inline bool operator<(const Proxy &a, const Proxy &b)
		{
			return a.tag < b.tag;
		}
		friend inline bool operator<(uint32 a, const Proxy &b)
		{
			return a < b.tag;
		}
		friend inline bool operator<(const Proxy &a, uint32 b)
		{
			return a.tag < b;
		}
	};

	/// Class for filtering pairs or triads.
	class ConnectionFilter
	{
	public:
		virtual ~ConnectionFilter() {}
		/// Is the particle necessary for connection?
		/// A pair or a triad should contain at least one 'necessary' particle.
		virtual bool IsNecessary(int32 index) const
		{
			B2_NOT_USED(index);
			return true;
		}
		/// An additional condition for creating a pair.
		virtual bool ShouldCreatePair(int32 a, int32 b) const
		{
			B2_NOT_USED(a);
			B2_NOT_USED(b);
			return true;
		}
		/// An additional condition for creating a triad.
		virtual bool ShouldCreateTriad(int32 a, int32 b, int32 c) const
		{
			B2_NOT_USED(a);
			B2_NOT_USED(b);
			B2_NOT_USED(c);
			return true;
		}
	};

	/// InsideBoundsEnumerator enumerates all particles inside the given bounds.
	class InsideBoundsEnumerator
	{
	public:
		/// Construct an enumerator with bounds of tags and a range of proxies.
		InsideBoundsEnumerator(
			uint32 lower, uint32 upper,
			const Proxy* first, const Proxy* last);

		/// Get index of the next particle. Returns b2_invalidParticleIndex if
		/// there are no more particles.
		int32 GetNext();
	private:
		/// The lower and upper bound of x component in the tag.
		uint32 m_xLower, m_xUpper;
		/// The lower and upper bound of y component in the tag.
		uint32 m_yLower, m_yUpper;
		/// The range of proxies.
		const Proxy* m_first;
		const Proxy* m_last;
	};

	/// Node of linked lists of connected particles
	struct ParticleListNode
	{
		/// The head of the list.
		ParticleListNode* list;
		/// The next node in the list.
		ParticleListNode* next;
		/// Number of entries in the list. Valid only for the node at the head
		/// of the list.
		int32 count;
		/// Particle index.
		int32 index;
	};

	/// All particle types that require creating pairs
	static const int32 k_pairFlags =
		b2_springParticle |
		b2_barrierParticle;
	/// All particle types that require creating triads
	static const int32 k_triadFlags =
		b2_elasticParticle;
	/// All particle types that do not produce dynamic pressure
	static const int32 k_noPressureFlags =
		b2_powderParticle |
		b2_tensileParticle;
	/// All particle types that apply extra damping force with bodies
	static const int32 k_extraDampingFlags =
		b2_staticPressureParticle;

	b2ParticleSystem(const b2ParticleSystemDef* def, b2World* world);
	~b2ParticleSystem();

	template <typename T> void FreeBuffer(T** b, int capacity);
	template <typename T> void FreeUserOverridableBuffer(
		UserOverridableBuffer<T>* b);
	template <typename T> T* ReallocateBuffer(T* buffer, int32 oldCapacity,
											  int32 newCapacity);
	template <typename T> T* ReallocateBuffer(
		T* buffer, int32 userSuppliedCapacity, int32 oldCapacity,
		int32 newCapacity, bool deferred);
	template <typename T> T* ReallocateBuffer(
		UserOverridableBuffer<T>* buffer, int32 oldCapacity, int32 newCapacity,
		bool deferred);
	template <typename T> T* RequestBuffer(T* buffer);

	/// Reallocate the handle / index map and schedule the allocation of a new
	/// pool for handle allocation.
	void ReallocateHandleBuffers(int32 newCapacity);

	void ReallocateInternalAllocatedBuffers(int32 capacity);
	int32 CreateParticleForGroup(
		const b2ParticleGroupDef& groupDef,
		const b2Transform& xf, const b2Vec2& position);
	void CreateParticlesStrokeShapeForGroup(
		const b2Shape* shape,
		const b2ParticleGroupDef& groupDef, const b2Transform& xf);
	void CreateParticlesFillShapeForGroup(
		const b2Shape* shape,
		const b2ParticleGroupDef& groupDef, const b2Transform& xf);
	void CreateParticlesWithShapeForGroup(
		const b2Shape* shape,
		const b2ParticleGroupDef& groupDef, const b2Transform& xf);
	void CreateParticlesWithShapesForGroup(
		const b2Shape* const* shapes, int32 shapeCount,
		const b2ParticleGroupDef& groupDef, const b2Transform& xf);
	int32 CloneParticle(int32 index, b2ParticleGroup* group);
	void DestroyParticleGroup(b2ParticleGroup* group);

	void UpdatePairsAndTriads(
		int32 firstIndex, int32 lastIndex, const ConnectionFilter& filter);
	void UpdatePairsAndTriadsWithReactiveParticles();
	static bool ComparePairIndices(const b2ParticlePair& a, const b2ParticlePair& b);
	static bool MatchPairIndices(const b2ParticlePair& a, const b2ParticlePair& b);
	static bool CompareTriadIndices(const b2ParticleTriad& a, const b2ParticleTriad& b);
	static bool MatchTriadIndices(const b2ParticleTriad& a, const b2ParticleTriad& b);

	static void InitializeParticleLists(
		const b2ParticleGroup* group, ParticleListNode* nodeBuffer);
	void MergeParticleListsInContact(
		const b2ParticleGroup* group, ParticleListNode* nodeBuffer) const;
	static void MergeParticleLists(
		ParticleListNode* listA, ParticleListNode* listB);
	static ParticleListNode* FindLongestParticleList(
		const b2ParticleGroup* group, ParticleListNode* nodeBuffer);
	void MergeZombieParticleListNodes(
		const b2ParticleGroup* group, ParticleListNode* nodeBuffer,
		ParticleListNode* survivingList) const;
	static void MergeParticleListAndNode(
		ParticleListNode* list, ParticleListNode* node);
	void CreateParticleGroupsFromParticleList(
		const b2ParticleGroup* group, ParticleListNode* nodeBuffer,
		const ParticleListNode* survivingList);
	void UpdatePairsAndTriadsWithParticleList(
		const b2ParticleGroup* group, const ParticleListNode* nodeBuffer);

	void ComputeDepth();

	InsideBoundsEnumerator GetInsideBoundsEnumerator(const b2AABB& aabb) const;

	void UpdateAllParticleFlags();
	void UpdateAllGroupFlags();
	void AddContact(int32 a, int32 b,
		b2GrowableBuffer<b2ParticleContact>& contacts) const;
	void FindContacts_Reference(
		b2GrowableBuffer<b2ParticleContact>& contacts) const;
	void ReorderForFindContact(FindContactInput* reordered,
		                       int alignedCount) const;
	void GatherChecksOneParticle(
		const uint32 bound,
		const int startIndex,
		const int particleIndex,
		int* nextUncheckedIndex,
		b2GrowableBuffer<FindContactCheck>& checks) const;
	void GatherChecks(b2GrowableBuffer<FindContactCheck>& checks) const;
	void FindContacts_Simd(
		b2GrowableBuffer<b2ParticleContact>& contacts) const;
	void FindContacts(
		b2GrowableBuffer<b2ParticleContact>& contacts) const;
	static void UpdateProxyTags(
		const uint32* const tags,
		b2GrowableBuffer<Proxy>& proxies);
	static bool ProxyBufferHasIndex(
		int32 index, const Proxy* const a, int count);
	static int NumProxiesWithSameTag(
		const Proxy* const a, const Proxy* const b, int count);
	static bool AreProxyBuffersTheSame(const b2GrowableBuffer<Proxy>& a,
								   	   const b2GrowableBuffer<Proxy>& b);
	void UpdateProxies_Reference(b2GrowableBuffer<Proxy>& proxies) const;
	void UpdateProxies_Simd(b2GrowableBuffer<Proxy>& proxies) const;
	void UpdateProxies(b2GrowableBuffer<Proxy>& proxies) const;
	void SortProxies(b2GrowableBuffer<Proxy>& proxies) const;
	void FilterContacts(b2GrowableBuffer<b2ParticleContact>& contacts);
	void NotifyContactListenerPreContact(
		b2ParticlePairSet* particlePairs) const;
	void NotifyContactListenerPostContact(b2ParticlePairSet& particlePairs);
	void UpdateContacts(bool exceptZombie);
	void NotifyBodyContactListenerPreContact(
		FixtureParticleSet* fixtureSet) const;
	void NotifyBodyContactListenerPostContact(FixtureParticleSet& fixtureSet);
	void UpdateBodyContacts();

	void Solve(const b2TimeStep& step);
	void SolveCollision(const b2TimeStep& step);
	void LimitVelocity(const b2TimeStep& step);
	void SolveGravity(const b2TimeStep& step);
	void SolveBarrier(const b2TimeStep& step);
	void SolveStaticPressure(const b2TimeStep& step);
	void ComputeWeight();
	void SolvePressure(const b2TimeStep& step);
	void SolveDamping(const b2TimeStep& step);
	void SolveRigidDamping();
	void SolveExtraDamping();
	void SolveWall();
	void SolveRigid(const b2TimeStep& step);
	void SolveElastic(const b2TimeStep& step);
	void SolveSpring(const b2TimeStep& step);
	void SolveTensile(const b2TimeStep& step);
	void SolveViscous();
	void SolveRepulsive(const b2TimeStep& step);
	void SolvePowder(const b2TimeStep& step);
	void SolveSolid(const b2TimeStep& step);
	void SolveForce(const b2TimeStep& step);
	void SolveColorMixing();
	void SolveZombie();
	/// Destroy all particles which have outlived their lifetimes set by
	/// SetParticleLifetime().
	void SolveLifetimes(const b2TimeStep& step);
	void RotateBuffer(int32 start, int32 mid, int32 end);

	float32 GetCriticalVelocity(const b2TimeStep& step) const;
	float32 GetCriticalVelocitySquared(const b2TimeStep& step) const;
	float32 GetCriticalPressure(const b2TimeStep& step) const;
	float32 GetParticleStride() const;
	float32 GetParticleMass() const;
	float32 GetParticleInvMass() const;

	// Get the world's contact filter if any particles with the
	// b2_contactFilterParticle flag are present in the system.
	b2ContactFilter* GetFixtureContactFilter() const;

	// Get the world's contact filter if any particles with the
	// b2_particleContactFilterParticle flag are present in the system.
	b2ContactFilter* GetParticleContactFilter() const;

	// Get the world's contact listener if any particles with the
	// b2_fixtureContactListenerParticle flag are present in the system.
	b2ContactListener* GetFixtureContactListener() const;

	// Get the world's contact listener if any particles with the
	// b2_particleContactListenerParticle flag are present in the system.
	b2ContactListener* GetParticleContactListener() const;

	template <typename T> void SetUserOverridableBuffer(
		UserOverridableBuffer<T>* buffer, T* newBufferData, int32 newCapacity);

	void SetGroupFlags(b2ParticleGroup* group, uint32 flags);

	void RemoveSpuriousBodyContacts();
	static bool BodyContactCompare(const b2ParticleBodyContact& lhs,
								   const b2ParticleBodyContact& rhs);

	void DetectStuckParticle(int32 particle);

	/// Determine whether a particle index is valid.
	bool ValidateParticleIndex(const int32 index) const;

	/// Get the time elapsed in b2ParticleSystemDef::lifetimeGranularity.
	int32 GetQuantizedTimeElapsed() const;
	/// Convert a lifetime in seconds to an expiration time.
	int64 LifetimeToExpirationTime(const float32 lifetime) const;

	bool ForceCanBeApplied(uint32 flags) const;
	void PrepareForceBuffer();

	bool IsRigidGroup(b2ParticleGroup *group) const;
	b2Vec2 GetLinearVelocity(
		b2ParticleGroup *group, int32 particleIndex,
		const b2Vec2 &point) const;
	void InitDampingParameter(
		float32* invMass, float32* invInertia, float32* tangentDistance,
		float32 mass, float32 inertia, const b2Vec2& center,
		const b2Vec2& point, const b2Vec2& normal) const;
	void InitDampingParameterWithRigidGroupOrParticle(
		float32* invMass, float32* invInertia, float32* tangentDistance,
		bool isRigidGroup, b2ParticleGroup* group, int32 particleIndex,
		const b2Vec2& point, const b2Vec2& normal) const;
	float32 ComputeDampingImpulse(
		float32 invMassA, float32 invInertiaA, float32 tangentDistanceA,
		float32 invMassB, float32 invInertiaB, float32 tangentDistanceB,
		float32 normalVelocity) const;
	void ApplyDamping(
		float32 invMass, float32 invInertia, float32 tangentDistance,
		bool isRigidGroup, b2ParticleGroup* group, int32 particleIndex,
		float32 impulse, const b2Vec2& normal);

	bool m_paused;
	int32 m_timestamp;
	int32 m_allParticleFlags;
	bool m_needsUpdateAllParticleFlags;
	int32 m_allGroupFlags;
	bool m_needsUpdateAllGroupFlags;
	bool m_hasForce;
	int32 m_iterationIndex;
	float32 m_inverseDensity;
	float32 m_particleDiameter;
	float32 m_inverseDiameter;
	float32 m_squaredDiameter;

	int32 m_count;
	int32 m_internalAllocatedCapacity;
	/// Allocator for b2ParticleHandle instances.
	b2SlabAllocator<b2ParticleHandle> m_handleAllocator;
	/// Maps particle indicies to  handles.
	UserOverridableBuffer<b2ParticleHandle*> m_handleIndexBuffer;
	UserOverridableBuffer<uint32> m_flagsBuffer;
	UserOverridableBuffer<b2Vec2> m_positionBuffer;
	UserOverridableBuffer<b2Vec2> m_velocityBuffer;
	b2Vec2* m_forceBuffer;
	/// m_weightBuffer is populated in ComputeWeight and used in
	/// ComputeDepth(), SolveStaticPressure() and SolvePressure().
	float32* m_weightBuffer;
	/// When any particles have the flag b2_staticPressureParticle,
	/// m_staticPressureBuffer is first allocated and used in
	/// SolveStaticPressure() and SolvePressure().  It will be reallocated on
	/// subsequent CreateParticle() calls.
	float32* m_staticPressureBuffer;
	/// m_accumulationBuffer is used in many functions as a temporary buffer
	/// for scalar values.
	float32* m_accumulationBuffer;
	/// When any particles have the flag b2_tensileParticle,
	/// m_accumulation2Buffer is first allocated and used in SolveTensile()
	/// as a temporary buffer for vector values.  It will be reallocated on
	/// subsequent CreateParticle() calls.
	b2Vec2* m_accumulation2Buffer;
	/// When any particle groups have the flag b2_solidParticleGroup,
	/// m_depthBuffer is first allocated and populated in ComputeDepth() and
	/// used in SolveSolid(). It will be reallocated on subsequent
	/// CreateParticle() calls.
	float32* m_depthBuffer;
	UserOverridableBuffer<b2ParticleColor> m_colorBuffer;
	b2ParticleGroup** m_groupBuffer;
	UserOverridableBuffer<void*> m_userDataBuffer;

	/// Stuck particle detection parameters and record keeping
	int32 m_stuckThreshold;
	UserOverridableBuffer<int32> m_lastBodyContactStepBuffer;
	UserOverridableBuffer<int32> m_bodyContactCountBuffer;
	UserOverridableBuffer<int32> m_consecutiveContactStepsBuffer;
	b2GrowableBuffer<int32> m_stuckParticleBuffer;
	b2GrowableBuffer<Proxy> m_proxyBuffer;
	b2GrowableBuffer<b2ParticleContact> m_contactBuffer;
	b2GrowableBuffer<b2ParticleBodyContact> m_bodyContactBuffer;
	b2GrowableBuffer<b2ParticlePair> m_pairBuffer;
	b2GrowableBuffer<b2ParticleTriad> m_triadBuffer;

	/// Time each particle should be destroyed relative to the last time
	/// m_timeElapsed was initialized.  Each unit of time corresponds to
	/// b2ParticleSystemDef::lifetimeGranularity seconds.
	UserOverridableBuffer<int32> m_expirationTimeBuffer;
	/// List of particle indices sorted by expiration time.
	UserOverridableBuffer<int32> m_indexByExpirationTimeBuffer;
	/// Time elapsed in 32:32 fixed point.  Each non-fractional unit of time
	/// corresponds to b2ParticleSystemDef::lifetimeGranularity seconds.
	int64 m_timeElapsed;
	/// Whether the expiration time buffer has been modified and needs to be
	/// resorted.
	bool m_expirationTimeBufferRequiresSorting;

	int32 m_groupCount;
	b2ParticleGroup* m_groupList;

	b2ParticleSystemDef m_def;

	b2World* m_world;
	b2ParticleSystem* m_prev;
	b2ParticleSystem* m_next;
};

inline void b2ParticleContact::SetIndices(int32 a, int32 b)
{
	b2Assert(a <= b2_maxParticleIndex && b <= b2_maxParticleIndex);
	indexA = (b2ParticleIndex)a;
	indexB = (b2ParticleIndex)b;
}


inline bool b2ParticleContact::operator==(
	const b2ParticleContact& rhs) const
{
	return indexA == rhs.indexA
		&& indexB == rhs.indexB
		&& flags == rhs.flags
		&& weight == rhs.weight
		&& normal == rhs.normal;
}

// The reciprocal sqrt function differs between SIMD and non-SIMD, but they
// should create approximately equal results.
inline bool b2ParticleContact::ApproximatelyEqual(
	const b2ParticleContact& rhs) const
{
	static const float MAX_WEIGHT_DIFF = 0.01f; // Weight 0 ~ 1, so about 1%
	static const float MAX_NORMAL_DIFF = 0.01f; // Normal length = 1, so 1%
	return indexA == rhs.indexA
		&& indexB == rhs.indexB
		&& flags == rhs.flags
		&& b2Abs(weight - rhs.weight) < MAX_WEIGHT_DIFF
		&& (normal - rhs.normal).Length() < MAX_NORMAL_DIFF;
}

inline b2ParticleGroup* b2ParticleSystem::GetParticleGroupList()
{
	return m_groupList;
}

inline const b2ParticleGroup* b2ParticleSystem::GetParticleGroupList() const
{
	return m_groupList;
}

inline int32 b2ParticleSystem::GetParticleGroupCount() const
{
	return m_groupCount;
}

inline int32 b2ParticleSystem::GetParticleCount() const
{
	return m_count;
}

inline void b2ParticleSystem::SetPaused(bool paused)
{
	m_paused = paused;
}

inline bool b2ParticleSystem::GetPaused() const
{
	return m_paused;
}

inline const b2ParticleContact* b2ParticleSystem::GetContacts() const
{
	return m_contactBuffer.Data();
}

inline int32 b2ParticleSystem::GetContactCount() const
{
	return m_contactBuffer.GetCount();
}

inline const b2ParticleBodyContact* b2ParticleSystem::GetBodyContacts() const
{
	return m_bodyContactBuffer.Data();
}

inline int32 b2ParticleSystem::GetBodyContactCount() const
{
	return m_bodyContactBuffer.GetCount();
}

inline const b2ParticlePair* b2ParticleSystem::GetPairs() const
{
	return m_pairBuffer.Data();
}

inline int32 b2ParticleSystem::GetPairCount() const
{
	return m_pairBuffer.GetCount();
}

inline const b2ParticleTriad* b2ParticleSystem::GetTriads() const
{
	return m_triadBuffer.Data();
}

inline int32 b2ParticleSystem::GetTriadCount() const
{
	return m_triadBuffer.GetCount();
}

inline b2ParticleSystem* b2ParticleSystem::GetNext()
{
	return m_next;
}

inline const b2ParticleSystem* b2ParticleSystem::GetNext() const
{
	return m_next;
}

inline const int32* b2ParticleSystem::GetStuckCandidates() const
{
	return m_stuckParticleBuffer.Data();
}

inline int32 b2ParticleSystem::GetStuckCandidateCount() const
{
	return m_stuckParticleBuffer.GetCount();
}

inline void b2ParticleSystem::SetStrictContactCheck(bool enabled)
{
	m_def.strictContactCheck = enabled;
}

inline bool b2ParticleSystem::GetStrictContactCheck() const
{
	return m_def.strictContactCheck;
}

inline void b2ParticleSystem::SetRadius(float32 radius)
{
	m_particleDiameter = 2 * radius;
	m_squaredDiameter = m_particleDiameter * m_particleDiameter;
	m_inverseDiameter = 1 / m_particleDiameter;
}

inline void b2ParticleSystem::SetDensity(float32 density)
{
	m_def.density = density;
	m_inverseDensity =  1 / m_def.density;
}

inline float32 b2ParticleSystem::GetDensity() const
{
	return m_def.density;
}

inline void b2ParticleSystem::SetGravityScale(float32 gravityScale)
{
	m_def.gravityScale = gravityScale;
}

inline float32 b2ParticleSystem::GetGravityScale() const
{
	return m_def.gravityScale;
}

inline void b2ParticleSystem::SetDamping(float32 damping)
{
	m_def.dampingStrength = damping;
}

inline float32 b2ParticleSystem::GetDamping() const
{
	return m_def.dampingStrength;
}

inline void b2ParticleSystem::SetStaticPressureIterations(int32 iterations)
{
	m_def.staticPressureIterations = iterations;
}

inline int32 b2ParticleSystem::GetStaticPressureIterations() const
{
	return m_def.staticPressureIterations;
}

inline float32 b2ParticleSystem::GetRadius() const
{
	return m_particleDiameter / 2;
}

inline float32 b2ParticleSystem::GetCriticalVelocity(const b2TimeStep& step) const
{
	return m_particleDiameter * step.inv_dt;
}

inline float32 b2ParticleSystem::GetCriticalVelocitySquared(
	const b2TimeStep& step) const
{
	float32 velocity = GetCriticalVelocity(step);
	return velocity * velocity;
}

inline float32 b2ParticleSystem::GetCriticalPressure(const b2TimeStep& step) const
{
	return m_def.density * GetCriticalVelocitySquared(step);
}

inline float32 b2ParticleSystem::GetParticleStride() const
{
	return b2_particleStride * m_particleDiameter;
}

inline float32 b2ParticleSystem::GetParticleMass() const
{
	float32 stride = GetParticleStride();
	return m_def.density * stride * stride;
}

inline float32 b2ParticleSystem::GetParticleInvMass() const
{
	// mass = density * stride^2, so we take the inverse of this.
	float32 inverseStride = m_inverseDiameter * (1.0f / b2_particleStride);
	return m_inverseDensity * inverseStride * inverseStride;
}

inline b2Vec2* b2ParticleSystem::GetPositionBuffer()
{
	return m_positionBuffer.data;
}

inline b2Vec2* b2ParticleSystem::GetVelocityBuffer()
{
	return m_velocityBuffer.data;
}

inline float32* b2ParticleSystem::GetWeightBuffer()
{
	return m_weightBuffer;
}

inline int32 b2ParticleSystem::GetMaxParticleCount() const
{
	return m_def.maxCount;
}

inline void b2ParticleSystem::SetMaxParticleCount(int32 count)
{
	b2Assert(m_count <= count);
	m_def.maxCount = count;
}

inline uint32 b2ParticleSystem::GetAllParticleFlags() const
{
	return m_allParticleFlags;
}

inline uint32 b2ParticleSystem::GetAllGroupFlags() const
{
	return m_allGroupFlags;
}

inline const uint32* b2ParticleSystem::GetFlagsBuffer() const
{
	return m_flagsBuffer.data;
}

inline const b2Vec2* b2ParticleSystem::GetPositionBuffer() const
{
	return m_positionBuffer.data;
}

inline const b2Vec2* b2ParticleSystem::GetVelocityBuffer() const
{
	return m_velocityBuffer.data;
}

inline const b2ParticleColor* b2ParticleSystem::GetColorBuffer() const
{
	return ((b2ParticleSystem*) this)->GetColorBuffer();
}

inline const b2ParticleGroup* const* b2ParticleSystem::GetGroupBuffer() const
{
	return m_groupBuffer;
}

inline const float32* b2ParticleSystem::GetWeightBuffer() const
{
	return m_weightBuffer;
}

inline void* const* b2ParticleSystem::GetUserDataBuffer() const
{
	return ((b2ParticleSystem*) this)->GetUserDataBuffer();
}

inline b2ParticleGroup* const* b2ParticleSystem::GetGroupBuffer()
{
	return m_groupBuffer;
}

inline uint32 b2ParticleSystem::GetParticleFlags(int32 index)
{
	return GetFlagsBuffer()[index];
}

inline bool b2ParticleSystem::ValidateParticleIndex(const int32 index) const
{
	return index >= 0 && index < GetParticleCount() &&
		index != b2_invalidParticleIndex;
}

inline bool b2ParticleSystem::GetDestructionByAge() const
{
	return m_def.destroyByAge;
}

inline void b2ParticleSystem::ParticleApplyLinearImpulse(int32 index,
														 const b2Vec2& impulse)
{
	ApplyLinearImpulse(index, index + 1, impulse);
}


// Note: These functions must go in the header so the unit tests will compile
// them. b2ParticleSystem.cpp does not compile with this #define.
#if LIQUIDFUN_EXTERNAL_LANGUAGE_API

inline void b2ParticleSystem::SetParticleVelocity(int32 index,
												  float32 vx,
												  float32 vy)
{
	b2Vec2& v = GetVelocityBuffer()[index];
	v.x = vx;
	v.y = vy;
}

inline float b2ParticleSystem::GetParticlePositionX(int32 index) const
{
	return GetPositionBuffer()[index].x;
}

inline float b2ParticleSystem::GetParticlePositionY(int32 index) const
{
	return GetPositionBuffer()[index].y;
}

inline int b2ParticleSystem::CopyPositionBuffer(int startIndex,
												int numParticles,
												void* outBuf,
												int size) const
{
	int copySize = numParticles * sizeof(b2Vec2);
	void* inBufWithOffset = (void*) (GetPositionBuffer() + startIndex);
	return CopyBuffer(startIndex, numParticles, inBufWithOffset, outBuf, size,
					  copySize);
}

inline int b2ParticleSystem::CopyColorBuffer(int startIndex,
											 int numParticles,
											 void* outBuf,
											 int size) const
{
	int copySize = numParticles * sizeof(b2ParticleColor);
	void* inBufWithOffset = (void*) (GetColorBuffer() + startIndex);
	return CopyBuffer(startIndex, numParticles, inBufWithOffset, outBuf, size,
					  copySize);
}

inline int b2ParticleSystem::CopyWeightBuffer(int startIndex,
											  int numParticles,
											  void* outBuf,
											  int size) const
{
	int copySize = numParticles * sizeof(float32);
	void* inBufWithOffset = (void*) (GetWeightBuffer() + startIndex);
	return CopyBuffer(startIndex, numParticles, inBufWithOffset, outBuf, size,
					  copySize);
}

inline int b2ParticleSystem::CopyBuffer(int startIndex, int numParticles,
										void* inBufWithOffset, void* outBuf,
										int outBufSize, int copySize) const
{
	b2ExceptionType exception = IsBufCopyValid(startIndex, numParticles,
											   copySize, outBufSize);
	if (exception != b2_noExceptions)
	{
		return exception;
	}

	memcpy(outBuf, inBufWithOffset, copySize);
	return b2_noExceptions;
}

#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

// end of ParticleSystem.h

struct b2AABB;
struct b2BodyDef;
struct b2Color;
struct b2JointDef;
class b2Body;
class b2Draw;
class b2Fixture;
class b2Joint;
class b2ParticleGroup;

/// The world class manages all physics entities, dynamic simulation,
/// and asynchronous queries. The world also contains efficient memory
/// management facilities.
class b2World
{
public:
	/// Construct a world object.
	/// @param gravity the world gravity vector.
	b2World(const b2Vec2& gravity);

	/// Destruct the world. All physics entities are destroyed and all heap memory is released.
	~b2World();

	/// Register a destruction listener. The listener is owned by you and must
	/// remain in scope.
	void SetDestructionListener(b2DestructionListener* listener);

	/// Register a contact filter to provide specific control over collision.
	/// Otherwise the default filter is used (b2_defaultFilter). The listener is
	/// owned by you and must remain in scope.
	void SetContactFilter(b2ContactFilter* filter);

	/// Register a contact event listener. The listener is owned by you and must
	/// remain in scope.
	void SetContactListener(b2ContactListener* listener);

	/// Register a routine for debug drawing. The debug draw functions are called
	/// inside with b2World::DrawDebugData method. The debug draw object is owned
	/// by you and must remain in scope.
	void SetDebugDraw(b2Draw* debugDraw);

	/// Create a rigid body given a definition. No reference to the definition
	/// is retained.
	/// @warning This function is locked during callbacks.
	b2Body* CreateBody(const b2BodyDef* def);

	/// Destroy a rigid body.
	/// This function is locked during callbacks.
	/// @warning This automatically deletes all associated shapes and joints.
	/// @warning This function is locked during callbacks.
	void DestroyBody(b2Body* body);

	/// Create a joint to constrain bodies together. No reference to the definition
	/// is retained. This may cause the connected bodies to cease colliding.
	/// @warning This function is locked during callbacks.
	b2Joint* CreateJoint(const b2JointDef* def);

	/// Destroy a joint. This may cause the connected bodies to begin colliding.
	/// @warning This function is locked during callbacks.
	void DestroyJoint(b2Joint* joint);

	/// Create a particle system given a definition. No reference to the
	/// definition is retained.
	/// @warning This function is locked during callbacks.
	b2ParticleSystem* CreateParticleSystem(const b2ParticleSystemDef* def);

	/// Destroy a particle system.
	/// @warning This function is locked during callbacks.
	void DestroyParticleSystem(b2ParticleSystem* p);

	/// Take a time step. This performs collision detection, integration,
	/// and constraint solution.
	/// For the numerical stability of particles, minimize the following
	/// dimensionless gravity acceleration:
	///     gravity / particleRadius * (timeStep / particleIterations)^2
	/// b2CalculateParticleIterations() or
	/// CalculateReasonableParticleIterations() help to determine the optimal
	/// particleIterations.
	/// @param timeStep the amount of time to simulate, this should not vary.
	/// @param velocityIterations for the velocity constraint solver.
	/// @param positionIterations for the position constraint solver.
	/// @param particleIterations for the particle simulation.
	void Step(	float32 timeStep,
				int32 velocityIterations,
				int32 positionIterations,
				int32 particleIterations);

	/// Take a time step. This performs collision detection, integration,
	/// and constraint solution.
	/// @param timeStep the amount of time to simulate, this should not vary.
	/// @param velocityIterations for the velocity constraint solver.
	/// @param positionIterations for the position constraint solver.
	void Step(	float32 timeStep,
				int32 velocityIterations,
				int32 positionIterations)
	{
		Step(timeStep, velocityIterations, positionIterations, 1);
	}

	/// Recommend a value to be used in `Step` for `particleIterations`.
	/// This calculation is necessarily a simplification and should only be
	/// used as a starting point. Please see "Particle Iterations" in the
	/// Programmer's Guide for details.
	/// @param timeStep is the value to be passed into `Step`.
	int CalculateReasonableParticleIterations(float32 timeStep) const;

	/// Manually clear the force buffer on all bodies. By default, forces are cleared automatically
	/// after each call to Step. The default behavior is modified by calling SetAutoClearForces.
	/// The purpose of this function is to support sub-stepping. Sub-stepping is often used to maintain
	/// a fixed sized time step under a variable frame-rate.
	/// When you perform sub-stepping you will disable auto clearing of forces and instead call
	/// ClearForces after all sub-steps are complete in one pass of your game loop.
	/// @see SetAutoClearForces
	void ClearForces();

	/// Call this to draw shapes and other debug draw data. This is intentionally non-const.
	void DrawDebugData();

	/// Query the world for all fixtures that potentially overlap the
	/// provided AABB.
	/// @param callback a user implemented callback class.
	/// @param aabb the query box.
	void QueryAABB(b2QueryCallback* callback, const b2AABB& aabb) const;

	/// Query the world for all fixtures that potentially overlap the
	/// provided shape's AABB. Calls QueryAABB internally.
	/// @param callback a user implemented callback class.
	/// @param shape the query shape
	/// @param xf the transform of the AABB
	void QueryShapeAABB(b2QueryCallback* callback, const b2Shape& shape,
	                    const b2Transform& xf) const;

	/// Ray-cast the world for all fixtures in the path of the ray. Your callback
	/// controls whether you get the closest point, any point, or n-points.
	/// The ray-cast ignores shapes that contain the starting point.
	/// @param callback a user implemented callback class.
	/// @param point1 the ray starting point
	/// @param point2 the ray ending point
	void RayCast(b2RayCastCallback* callback, const b2Vec2& point1, const b2Vec2& point2) const;

	/// Get the world body list. With the returned body, use b2Body::GetNext to get
	/// the next body in the world list. A NULL body indicates the end of the list.
	/// @return the head of the world body list.
	b2Body* GetBodyList();
	const b2Body* GetBodyList() const;

	/// Get the world joint list. With the returned joint, use b2Joint::GetNext to get
	/// the next joint in the world list. A NULL joint indicates the end of the list.
	/// @return the head of the world joint list.
	b2Joint* GetJointList();
	const b2Joint* GetJointList() const;

	/// Get the world particle-system list. With the returned body, use
	/// b2ParticleSystem::GetNext to get the next particle-system in the world
	/// list. A NULL particle-system indicates the end of the list.
	/// @return the head of the world particle-system list.
	b2ParticleSystem* GetParticleSystemList();
	const b2ParticleSystem* GetParticleSystemList() const;

	/// Get the world contact list. With the returned contact, use b2Contact::GetNext to get
	/// the next contact in the world list. A NULL contact indicates the end of the list.
	/// @return the head of the world contact list.
	/// @warning contacts are created and destroyed in the middle of a time step.
	/// Use b2ContactListener to avoid missing contacts.
	b2Contact* GetContactList();
	const b2Contact* GetContactList() const;

	/// Enable/disable sleep.
	void SetAllowSleeping(bool flag);
	bool GetAllowSleeping() const { return m_allowSleep; }

	/// Enable/disable warm starting. For testing.
	void SetWarmStarting(bool flag) { m_warmStarting = flag; }
	bool GetWarmStarting() const { return m_warmStarting; }

	/// Enable/disable continuous physics. For testing.
	void SetContinuousPhysics(bool flag) { m_continuousPhysics = flag; }
	bool GetContinuousPhysics() const { return m_continuousPhysics; }

	/// Enable/disable single stepped continuous physics. For testing.
	void SetSubStepping(bool flag) { m_subStepping = flag; }
	bool GetSubStepping() const { return m_subStepping; }

	/// Get the number of broad-phase proxies.
	int32 GetProxyCount() const;

	/// Get the number of bodies.
	int32 GetBodyCount() const;

	/// Get the number of joints.
	int32 GetJointCount() const;

	/// Get the number of contacts (each may have 0 or more contact points).
	int32 GetContactCount() const;

	/// Get the height of the dynamic tree.
	int32 GetTreeHeight() const;

	/// Get the balance of the dynamic tree.
	int32 GetTreeBalance() const;

	/// Get the quality metric of the dynamic tree. The smaller the better.
	/// The minimum is 1.
	float32 GetTreeQuality() const;

	/// Change the global gravity vector.
	void SetGravity(const b2Vec2& gravity);

	/// Get the global gravity vector.
	b2Vec2 GetGravity() const;

	/// Is the world locked (in the middle of a time step).
	bool IsLocked() const;

	/// Set flag to control automatic clearing of forces after each time step.
	void SetAutoClearForces(bool flag);

	/// Get the flag that controls automatic clearing of forces after each time step.
	bool GetAutoClearForces() const;

	/// Shift the world origin. Useful for large worlds.
	/// The body shift formula is: position -= newOrigin
	/// @param newOrigin the new origin with respect to the old origin
	void ShiftOrigin(const b2Vec2& newOrigin);

	/// Get the contact manager for testing.
	const b2ContactManager& GetContactManager() const;

	/// Get the current profile.
	const b2Profile& GetProfile() const;

	/// Dump the world into the log file.
	/// @warning this should be called outside of a time step.
	void Dump();

	/// Get API version.
	const b2Version* GetVersion() const {
		return m_liquidFunVersion;
	}

	/// Get API version string.
	const char* GetVersionString() const {
		return m_liquidFunVersionString;
	}

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
public:
	/// Constructor which takes direct floats.
	b2World(float32 gravityX, float32 gravityY);

	/// Set gravity with direct floats.
	void SetGravity(float32 gravityX, float32 gravityY);
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

private:

	// m_flags
	enum
	{
		e_newFixture	= 0x0001,
		e_locked		= 0x0002,
		e_clearForces	= 0x0004
	};

	friend class b2Body;
	friend class b2Fixture;
	friend class b2ContactManager;
	friend class b2Controller;
	friend class b2ParticleSystem;

	void Init(const b2Vec2& gravity);

	void Solve(const b2TimeStep& step);
	void SolveTOI(const b2TimeStep& step);

	void DrawJoint(b2Joint* joint);
	void DrawShape(b2Fixture* shape, const b2Transform& xf, const b2Color& color);

	void DrawParticleSystem(const b2ParticleSystem& system);

	b2BlockAllocator m_blockAllocator;
	b2StackAllocator m_stackAllocator;

	int32 m_flags;

	b2ContactManager m_contactManager;

	b2Body* m_bodyList;
	b2Joint* m_jointList;
	b2ParticleSystem* m_particleSystemList;

	int32 m_bodyCount;
	int32 m_jointCount;

	b2Vec2 m_gravity;
	bool m_allowSleep;

	b2DestructionListener* m_destructionListener;
	b2Draw* m_debugDraw;

	// This is used to compute the time step ratio to
	// support a variable time step.
	float32 m_inv_dt0;

	// These are for debugging the solver.
	bool m_warmStarting;
	bool m_continuousPhysics;
	bool m_subStepping;

	bool m_stepComplete;

	b2Profile m_profile;

	/// Used to reference b2_LiquidFunVersion so that it's not stripped from
	/// the static library.
	const b2Version *m_liquidFunVersion;
	const char *m_liquidFunVersionString;
};

inline b2Body* b2World::GetBodyList()
{
	return m_bodyList;
}

inline const b2Body* b2World::GetBodyList() const
{
	return m_bodyList;
}

inline b2Joint* b2World::GetJointList()
{
	return m_jointList;
}

inline const b2Joint* b2World::GetJointList() const
{
	return m_jointList;
}

inline b2ParticleSystem* b2World::GetParticleSystemList()
{
	return m_particleSystemList;
}

inline const b2ParticleSystem* b2World::GetParticleSystemList() const
{
	return m_particleSystemList;
}

inline b2Contact* b2World::GetContactList()
{
	return m_contactManager.m_contactList;
}

inline const b2Contact* b2World::GetContactList() const
{
	return m_contactManager.m_contactList;
}

inline int32 b2World::GetBodyCount() const
{
	return m_bodyCount;
}

inline int32 b2World::GetJointCount() const
{
	return m_jointCount;
}

inline int32 b2World::GetContactCount() const
{
	return m_contactManager.m_contactCount;
}

inline void b2World::SetGravity(const b2Vec2& gravity)
{
	m_gravity = gravity;
}

inline b2Vec2 b2World::GetGravity() const
{
	return m_gravity;
}

inline bool b2World::IsLocked() const
{
	return (m_flags & e_locked) == e_locked;
}

inline void b2World::SetAutoClearForces(bool flag)
{
	if (flag)
	{
		m_flags |= e_clearForces;
	}
	else
	{
		m_flags &= ~e_clearForces;
	}
}

/// Get the flag that controls automatic clearing of forces after each time step.
inline bool b2World::GetAutoClearForces() const
{
	return (m_flags & e_clearForces) == e_clearForces;
}

inline const b2ContactManager& b2World::GetContactManager() const
{
	return m_contactManager;
}

inline const b2Profile& b2World::GetProfile() const
{
	return m_profile;
}

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
inline b2World::b2World(float32 gravityX, float32 gravityY)
{
	Init(b2Vec2(gravityX, gravityY));
}

inline void b2World::SetGravity(float32 gravityX, float32 gravityY)
{
	SetGravity(b2Vec2(gravityX, gravityY));
}
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

// end of World.h

b2Fixture::b2Fixture()
{
	m_userData = NULL;
	m_body = NULL;
	m_next = NULL;
	m_proxies = NULL;
	m_proxyCount = 0;
	m_shape = NULL;
	m_density = 0.0f;
}

void b2Fixture::Create(b2BlockAllocator* allocator, b2Body* body, const b2FixtureDef* def)
{
	m_userData = def->userData;
	m_friction = def->friction;
	m_restitution = def->restitution;

	m_body = body;
	m_next = NULL;

	m_filter = def->filter;

	m_isSensor = def->isSensor;

	m_shape = def->shape->Clone(allocator);

	// Reserve proxy space
	int32 childCount = m_shape->GetChildCount();
	m_proxies = (b2FixtureProxy*)allocator->Allocate(childCount * sizeof(b2FixtureProxy));
	for (int32 i = 0; i < childCount; ++i)
	{
		m_proxies[i].fixture = NULL;
		m_proxies[i].proxyId = b2BroadPhase::e_nullProxy;
	}
	m_proxyCount = 0;

	m_density = def->density;
}

void b2Fixture::Destroy(b2BlockAllocator* allocator)
{
	// The proxies must be destroyed before calling this.
	b2Assert(m_proxyCount == 0);

	// Free the proxy array.
	int32 childCount = m_shape->GetChildCount();
	allocator->Free(m_proxies, childCount * sizeof(b2FixtureProxy));
	m_proxies = NULL;

	// Free the child shape.
	switch (m_shape->m_type)
	{
	case b2Shape::e_circle:
		{
			b2CircleShape* s = (b2CircleShape*)m_shape;
			s->~b2CircleShape();
			allocator->Free(s, sizeof(b2CircleShape));
		}
		break;

	case b2Shape::e_edge:
		{
			b2EdgeShape* s = (b2EdgeShape*)m_shape;
			s->~b2EdgeShape();
			allocator->Free(s, sizeof(b2EdgeShape));
		}
		break;

	case b2Shape::e_polygon:
		{
			b2PolygonShape* s = (b2PolygonShape*)m_shape;
			s->~b2PolygonShape();
			allocator->Free(s, sizeof(b2PolygonShape));
		}
		break;

	case b2Shape::e_chain:
		{
			b2ChainShape* s = (b2ChainShape*)m_shape;
			s->~b2ChainShape();
			allocator->Free(s, sizeof(b2ChainShape));
		}
		break;

	default:
		b2Assert(false);
		break;
	}

	m_shape = NULL;
}

void b2Fixture::CreateProxies(b2BroadPhase* broadPhase, const b2Transform& xf)
{
	b2Assert(m_proxyCount == 0);

	// Create proxies in the broad-phase.
	m_proxyCount = m_shape->GetChildCount();

	for (int32 i = 0; i < m_proxyCount; ++i)
	{
		b2FixtureProxy* proxy = m_proxies + i;
		m_shape->ComputeAABB(&proxy->aabb, xf, i);
		proxy->proxyId = broadPhase->CreateProxy(proxy->aabb, proxy);
		proxy->fixture = this;
		proxy->childIndex = i;
	}
}

void b2Fixture::DestroyProxies(b2BroadPhase* broadPhase)
{
	// Destroy proxies in the broad-phase.
	for (int32 i = 0; i < m_proxyCount; ++i)
	{
		b2FixtureProxy* proxy = m_proxies + i;
		broadPhase->DestroyProxy(proxy->proxyId);
		proxy->proxyId = b2BroadPhase::e_nullProxy;
	}

	m_proxyCount = 0;
}

void b2Fixture::Synchronize(b2BroadPhase* broadPhase, const b2Transform& transform1, const b2Transform& transform2)
{
	if (m_proxyCount == 0)
	{	
		return;
	}

	for (int32 i = 0; i < m_proxyCount; ++i)
	{
		b2FixtureProxy* proxy = m_proxies + i;

		// Compute an AABB that covers the swept shape (may miss some rotation effect).
		b2AABB aabb1, aabb2;
		m_shape->ComputeAABB(&aabb1, transform1, proxy->childIndex);
		m_shape->ComputeAABB(&aabb2, transform2, proxy->childIndex);
	
		proxy->aabb.Combine(aabb1, aabb2);

		b2Vec2 displacement = transform2.p - transform1.p;

		broadPhase->MoveProxy(proxy->proxyId, proxy->aabb, displacement);
	}
}

void b2Fixture::SetFilterData(const b2Filter& filter)
{
	m_filter = filter;

	Refilter();
}

void b2Fixture::Refilter()
{
	if (m_body == NULL)
	{
		return;
	}

	// Flag associated contacts for filtering.
	b2ContactEdge* edge = m_body->GetContactList();
	while (edge)
	{
		b2Contact* contact = edge->contact;
		b2Fixture* fixtureA = contact->GetFixtureA();
		b2Fixture* fixtureB = contact->GetFixtureB();
		if (fixtureA == this || fixtureB == this)
		{
			contact->FlagForFiltering();
		}

		edge = edge->next;
	}

	b2World* world = m_body->GetWorld();

	if (world == NULL)
	{
		return;
	}

	// Touch each proxy so that new pairs may be created
	b2BroadPhase* broadPhase = &world->m_contactManager.m_broadPhase;
	for (int32 i = 0; i < m_proxyCount; ++i)
	{
		broadPhase->TouchProxy(m_proxies[i].proxyId);
	}
}

void b2Fixture::SetSensor(bool sensor)
{
	if (sensor != m_isSensor)
	{
		m_body->SetAwake(true);
		m_isSensor = sensor;
	}
}

void b2Fixture::Dump(int32 bodyIndex)
{
	b2Log("    b2FixtureDef fd;\n");
	b2Log("    fd.friction = %.15lef;\n", m_friction);
	b2Log("    fd.restitution = %.15lef;\n", m_restitution);
	b2Log("    fd.density = %.15lef;\n", m_density);
	b2Log("    fd.isSensor = bool(%d);\n", m_isSensor);
	b2Log("    fd.filter.categoryBits = uint16(%d);\n", m_filter.categoryBits);
	b2Log("    fd.filter.maskBits = uint16(%d);\n", m_filter.maskBits);
	b2Log("    fd.filter.groupIndex = int16(%d);\n", m_filter.groupIndex);

	switch (m_shape->m_type)
	{
	case b2Shape::e_circle:
		{
			b2CircleShape* s = (b2CircleShape*)m_shape;
			b2Log("    b2CircleShape shape;\n");
			b2Log("    shape.m_radius = %.15lef;\n", s->m_radius);
			b2Log("    shape.m_p.Set(%.15lef, %.15lef);\n", s->m_p.x, s->m_p.y);
		}
		break;

	case b2Shape::e_edge:
		{
			b2EdgeShape* s = (b2EdgeShape*)m_shape;
			b2Log("    b2EdgeShape shape;\n");
			b2Log("    shape.m_radius = %.15lef;\n", s->m_radius);
			b2Log("    shape.m_vertex0.Set(%.15lef, %.15lef);\n", s->m_vertex0.x, s->m_vertex0.y);
			b2Log("    shape.m_vertex1.Set(%.15lef, %.15lef);\n", s->m_vertex1.x, s->m_vertex1.y);
			b2Log("    shape.m_vertex2.Set(%.15lef, %.15lef);\n", s->m_vertex2.x, s->m_vertex2.y);
			b2Log("    shape.m_vertex3.Set(%.15lef, %.15lef);\n", s->m_vertex3.x, s->m_vertex3.y);
			b2Log("    shape.m_hasVertex0 = bool(%d);\n", s->m_hasVertex0);
			b2Log("    shape.m_hasVertex3 = bool(%d);\n", s->m_hasVertex3);
		}
		break;

	case b2Shape::e_polygon:
		{
			b2PolygonShape* s = (b2PolygonShape*)m_shape;
			b2Log("    b2PolygonShape shape;\n");
			b2Log("    b2Vec2 vs[%d];\n", b2_maxPolygonVertices);
			for (int32 i = 0; i < s->m_count; ++i)
			{
				b2Log("    vs[%d].Set(%.15lef, %.15lef);\n", i, s->m_vertices[i].x, s->m_vertices[i].y);
			}
			b2Log("    shape.Set(vs, %d);\n", s->m_count);
		}
		break;

	case b2Shape::e_chain:
		{
			b2ChainShape* s = (b2ChainShape*)m_shape;
			b2Log("    b2ChainShape shape;\n");
			b2Log("    b2Vec2 vs[%d];\n", s->m_count);
			for (int32 i = 0; i < s->m_count; ++i)
			{
				b2Log("    vs[%d].Set(%.15lef, %.15lef);\n", i, s->m_vertices[i].x, s->m_vertices[i].y);
			}
			b2Log("    shape.CreateChain(vs, %d);\n", s->m_count);
			b2Log("    shape.m_prevVertex.Set(%.15lef, %.15lef);\n", s->m_prevVertex.x, s->m_prevVertex.y);
			b2Log("    shape.m_nextVertex.Set(%.15lef, %.15lef);\n", s->m_nextVertex.x, s->m_nextVertex.y);
			b2Log("    shape.m_hasPrevVertex = bool(%d);\n", s->m_hasPrevVertex);
			b2Log("    shape.m_hasNextVertex = bool(%d);\n", s->m_hasNextVertex);
		}
		break;

	default:
		return;
	}

	b2Log("\n");
	b2Log("    fd.shape = &shape;\n");
	b2Log("\n");
	b2Log("    bodies[%d]->CreateFixture(&fd);\n", bodyIndex);
}

// end of Fixture.cpp

class b2Body;
class b2Joint;
struct b2SolverData;
class b2BlockAllocator;

enum b2JointType
{
	e_unknownJoint,
	e_revoluteJoint,
	e_prismaticJoint,
	e_distanceJoint,
	e_pulleyJoint,
	e_mouseJoint,
	e_gearJoint,
	e_wheelJoint,
    e_weldJoint,
	e_frictionJoint,
	e_ropeJoint,
	e_motorJoint
};

enum b2LimitState
{
	e_inactiveLimit,
	e_atLowerLimit,
	e_atUpperLimit,
	e_equalLimits
};

struct b2Jacobian
{
	b2Vec2 linear;
	float32 angularA;
	float32 angularB;
};

/// A joint edge is used to connect bodies and joints together
/// in a joint graph where each body is a node and each joint
/// is an edge. A joint edge belongs to a doubly linked list
/// maintained in each attached body. Each joint has two joint
/// nodes, one for each attached body.
struct b2JointEdge
{
	b2Body* other;			///< provides quick access to the other body attached.
	b2Joint* joint;			///< the joint
	b2JointEdge* prev;		///< the previous joint edge in the body's joint list
	b2JointEdge* next;		///< the next joint edge in the body's joint list
};

/// Joint definitions are used to construct joints.
struct b2JointDef
{
	b2JointDef()
	{
		type = e_unknownJoint;
		userData = NULL;
		bodyA = NULL;
		bodyB = NULL;
		collideConnected = false;
	}

	/// The joint type is set automatically for concrete joint types.
	b2JointType type;

	/// Use this to attach application specific data to your joints.
	void* userData;

	/// The first attached body.
	b2Body* bodyA;

	/// The second attached body.
	b2Body* bodyB;

	/// Set this flag to true if the attached bodies should collide.
	bool collideConnected;
};

/// The base joint class. Joints are used to constraint two bodies together in
/// various fashions. Some joints also feature limits and motors.
class b2Joint
{
public:

	/// Get the type of the concrete joint.
	b2JointType GetType() const;

	/// Get the first body attached to this joint.
	b2Body* GetBodyA();

	/// Get the second body attached to this joint.
	b2Body* GetBodyB();

	/// Get the anchor point on bodyA in world coordinates.
	virtual b2Vec2 GetAnchorA() const = 0;

	/// Get the anchor point on bodyB in world coordinates.
	virtual b2Vec2 GetAnchorB() const = 0;

	/// Get the reaction force on bodyB at the joint anchor in Newtons.
	virtual b2Vec2 GetReactionForce(float32 inv_dt) const = 0;

	/// Get the reaction torque on bodyB in N*m.
	virtual float32 GetReactionTorque(float32 inv_dt) const = 0;

	/// Get the next joint the world joint list.
	b2Joint* GetNext();
	const b2Joint* GetNext() const;

	/// Get the user data pointer.
	void* GetUserData() const;

	/// Set the user data pointer.
	void SetUserData(void* data);

	/// Short-cut function to determine if either body is inactive.
	bool IsActive() const;

	/// Get collide connected.
	/// Note: modifying the collide connect flag won't work correctly because
	/// the flag is only checked when fixture AABBs begin to overlap.
	bool GetCollideConnected() const;

	/// Dump this joint to the log file.
	virtual void Dump() { b2Log("// Dump is not supported for this joint type.\n"); }

	/// Shift the origin for any points stored in world coordinates.
	virtual void ShiftOrigin(const b2Vec2& newOrigin) { B2_NOT_USED(newOrigin);  }

protected:
	friend class b2World;
	friend class b2Body;
	friend class b2Island;
	friend class b2GearJoint;

	static b2Joint* Create(const b2JointDef* def, b2BlockAllocator* allocator);
	static void Destroy(b2Joint* joint, b2BlockAllocator* allocator);

	b2Joint(const b2JointDef* def);
	virtual ~b2Joint() {}

	virtual void InitVelocityConstraints(const b2SolverData& data) = 0;
	virtual void SolveVelocityConstraints(const b2SolverData& data) = 0;

	// This returns true if the position errors are within tolerance.
	virtual bool SolvePositionConstraints(const b2SolverData& data) = 0;

	b2JointType m_type;
	b2Joint* m_prev;
	b2Joint* m_next;
	b2JointEdge m_edgeA;
	b2JointEdge m_edgeB;
	b2Body* m_bodyA;
	b2Body* m_bodyB;

	int32 m_index;

	bool m_islandFlag;
	bool m_collideConnected;

	void* m_userData;
};

inline b2JointType b2Joint::GetType() const
{
	return m_type;
}

inline b2Body* b2Joint::GetBodyA()
{
	return m_bodyA;
}

inline b2Body* b2Joint::GetBodyB()
{
	return m_bodyB;
}

inline b2Joint* b2Joint::GetNext()
{
	return m_next;
}

inline const b2Joint* b2Joint::GetNext() const
{
	return m_next;
}

inline void* b2Joint::GetUserData() const
{
	return m_userData;
}

inline void b2Joint::SetUserData(void* data)
{
	m_userData = data;
}

inline bool b2Joint::GetCollideConnected() const
{
	return m_collideConnected;
}

// end of Joint.h

b2Body::b2Body(const b2BodyDef* bd, b2World* world)
{
	b2Assert(bd->position.IsValid());
	b2Assert(bd->linearVelocity.IsValid());
	b2Assert(b2IsValid(bd->angle));
	b2Assert(b2IsValid(bd->angularVelocity));
	b2Assert(b2IsValid(bd->angularDamping) && bd->angularDamping >= 0.0f);
	b2Assert(b2IsValid(bd->linearDamping) && bd->linearDamping >= 0.0f);

	m_flags = 0;

	if (bd->bullet)
	{
		m_flags |= e_bulletFlag;
	}
	if (bd->fixedRotation)
	{
		m_flags |= e_fixedRotationFlag;
	}
	if (bd->allowSleep)
	{
		m_flags |= e_autoSleepFlag;
	}
	if (bd->awake)
	{
		m_flags |= e_awakeFlag;
	}
	if (bd->active)
	{
		m_flags |= e_activeFlag;
	}

	m_world = world;

	m_xf.p = bd->position;
	m_xf.q.Set(bd->angle);
	m_xf0 = m_xf;

	m_sweep.localCenter.SetZero();
	m_sweep.c0 = m_xf.p;
	m_sweep.c = m_xf.p;
	m_sweep.a0 = bd->angle;
	m_sweep.a = bd->angle;
	m_sweep.alpha0 = 0.0f;

	m_jointList = NULL;
	m_contactList = NULL;
	m_prev = NULL;
	m_next = NULL;

	m_linearVelocity = bd->linearVelocity;
	m_angularVelocity = bd->angularVelocity;

	m_linearDamping = bd->linearDamping;
	m_angularDamping = bd->angularDamping;
	m_gravityScale = bd->gravityScale;

	m_force.SetZero();
	m_torque = 0.0f;

	m_sleepTime = 0.0f;

	m_type = bd->type;

	if (m_type == b2_dynamicBody)
	{
		m_mass = 1.0f;
		m_invMass = 1.0f;
	}
	else
	{
		m_mass = 0.0f;
		m_invMass = 0.0f;
	}

	m_I = 0.0f;
	m_invI = 0.0f;

	m_userData = bd->userData;

	m_fixtureList = NULL;
	m_fixtureCount = 0;
}

b2Body::~b2Body()
{
	// shapes and joints are destroyed in b2World::Destroy
}

void b2Body::SetType(b2BodyType type)
{
	b2Assert(m_world->IsLocked() == false);
	if (m_world->IsLocked() == true)
	{
		return;
	}

	if (m_type == type)
	{
		return;
	}

	m_type = type;

	ResetMassData();

	if (m_type == b2_staticBody)
	{
		m_linearVelocity.SetZero();
		m_angularVelocity = 0.0f;
		m_sweep.a0 = m_sweep.a;
		m_sweep.c0 = m_sweep.c;
		SynchronizeFixtures();
	}

	SetAwake(true);

	m_force.SetZero();
	m_torque = 0.0f;

	// Delete the attached contacts.
	b2ContactEdge* ce = m_contactList;
	while (ce)
	{
		b2ContactEdge* ce0 = ce;
		ce = ce->next;
		m_world->m_contactManager.Destroy(ce0->contact);
	}
	m_contactList = NULL;

	// Touch the proxies so that new contacts will be created (when appropriate)
	b2BroadPhase* broadPhase = &m_world->m_contactManager.m_broadPhase;
	for (b2Fixture* f = m_fixtureList; f; f = f->m_next)
	{
		int32 proxyCount = f->m_proxyCount;
		for (int32 i = 0; i < proxyCount; ++i)
		{
			broadPhase->TouchProxy(f->m_proxies[i].proxyId);
		}
	}
}

b2Fixture* b2Body::CreateFixture(const b2FixtureDef* def)
{
	b2Assert(m_world->IsLocked() == false);
	if (m_world->IsLocked() == true)
	{
		return NULL;
	}

	b2BlockAllocator* allocator = &m_world->m_blockAllocator;

	void* memory = allocator->Allocate(sizeof(b2Fixture));
	b2Fixture* fixture = new (memory) b2Fixture;
	fixture->Create(allocator, this, def);

	if (m_flags & e_activeFlag)
	{
		b2BroadPhase* broadPhase = &m_world->m_contactManager.m_broadPhase;
		fixture->CreateProxies(broadPhase, m_xf);
	}

	fixture->m_next = m_fixtureList;
	m_fixtureList = fixture;
	++m_fixtureCount;

	fixture->m_body = this;

	// Adjust mass properties if needed.
	if (fixture->m_density > 0.0f)
	{
		ResetMassData();
	}

	// Let the world know we have a new fixture. This will cause new contacts
	// to be created at the beginning of the next time step.
	m_world->m_flags |= b2World::e_newFixture;

	return fixture;
}

b2Fixture* b2Body::CreateFixture(const b2Shape* shape, float32 density)
{
	b2FixtureDef def;
	def.shape = shape;
	def.density = density;

	return CreateFixture(&def);
}

void b2Body::DestroyFixture(b2Fixture* fixture)
{
	b2Assert(m_world->IsLocked() == false);
	if (m_world->IsLocked() == true)
	{
		return;
	}

	b2Assert(fixture->m_body == this);

	// Remove the fixture from this body's singly linked list.
	b2Assert(m_fixtureCount > 0);
	b2Fixture** node = &m_fixtureList;
#if B2_ASSERT_ENABLED
	bool found = false;
#endif // B2_ASSERT_ENABLED
	while (*node != NULL)
	{
		if (*node == fixture)
		{
			*node = fixture->m_next;
#if B2_ASSERT_ENABLED
			found = true;
#endif // B2_ASSERT_ENABLED
			break;
		}

		node = &(*node)->m_next;
	}

	// You tried to remove a shape that is not attached to this body.
	b2Assert(found);

	// Destroy any contacts associated with the fixture.
	b2ContactEdge* edge = m_contactList;
	while (edge)
	{
		b2Contact* c = edge->contact;
		edge = edge->next;

		b2Fixture* fixtureA = c->GetFixtureA();
		b2Fixture* fixtureB = c->GetFixtureB();

		if (fixture == fixtureA || fixture == fixtureB)
		{
			// This destroys the contact and removes it from
			// this body's contact list.
			m_world->m_contactManager.Destroy(c);
		}
	}

	b2BlockAllocator* allocator = &m_world->m_blockAllocator;

	if (m_flags & e_activeFlag)
	{
		b2BroadPhase* broadPhase = &m_world->m_contactManager.m_broadPhase;
		fixture->DestroyProxies(broadPhase);
	}

	fixture->Destroy(allocator);
	fixture->m_body = NULL;
	fixture->m_next = NULL;
	fixture->~b2Fixture();
	allocator->Free(fixture, sizeof(b2Fixture));

	--m_fixtureCount;

	// Reset the mass data.
	ResetMassData();
}

void b2Body::ResetMassData()
{
	// Compute mass data from shapes. Each shape has its own density.
	m_mass = 0.0f;
	m_invMass = 0.0f;
	m_I = 0.0f;
	m_invI = 0.0f;
	m_sweep.localCenter.SetZero();

	// Static and kinematic bodies have zero mass.
	if (m_type == b2_staticBody || m_type == b2_kinematicBody)
	{
		m_sweep.c0 = m_xf.p;
		m_sweep.c = m_xf.p;
		m_sweep.a0 = m_sweep.a;
		return;
	}

	b2Assert(m_type == b2_dynamicBody);

	// Accumulate mass over all fixtures.
	b2Vec2 localCenter = b2Vec2_zero;
	for (b2Fixture* f = m_fixtureList; f; f = f->m_next)
	{
		if (f->m_density == 0.0f)
		{
			continue;
		}

		b2MassData massData;
		f->GetMassData(&massData);
		m_mass += massData.mass;
		localCenter += massData.mass * massData.center;
		m_I += massData.I;
	}

	// Compute center of mass.
	if (m_mass > 0.0f)
	{
		m_invMass = 1.0f / m_mass;
		localCenter *= m_invMass;
	}
	else
	{
		// Force all dynamic bodies to have a positive mass.
		m_mass = 1.0f;
		m_invMass = 1.0f;
	}

	if (m_I > 0.0f && (m_flags & e_fixedRotationFlag) == 0)
	{
		// Center the inertia about the center of mass.
		m_I -= m_mass * b2Dot(localCenter, localCenter);
		b2Assert(m_I > 0.0f);
		m_invI = 1.0f / m_I;

	}
	else
	{
		m_I = 0.0f;
		m_invI = 0.0f;
	}

	// Move center of mass.
	b2Vec2 oldCenter = m_sweep.c;
	m_sweep.localCenter = localCenter;
	m_sweep.c0 = m_sweep.c = b2Mul(m_xf, m_sweep.localCenter);

	// Update center of mass velocity.
	m_linearVelocity += b2Cross(m_angularVelocity, m_sweep.c - oldCenter);
}

void b2Body::SetMassData(const b2MassData* massData)
{
	b2Assert(m_world->IsLocked() == false);
	if (m_world->IsLocked() == true)
	{
		return;
	}

	if (m_type != b2_dynamicBody)
	{
		return;
	}

	m_invMass = 0.0f;
	m_I = 0.0f;
	m_invI = 0.0f;

	m_mass = massData->mass;
	if (m_mass <= 0.0f)
	{
		m_mass = 1.0f;
	}

	m_invMass = 1.0f / m_mass;

	if (massData->I > 0.0f && (m_flags & b2Body::e_fixedRotationFlag) == 0)
	{
		m_I = massData->I - m_mass * b2Dot(massData->center, massData->center);
		b2Assert(m_I > 0.0f);
		m_invI = 1.0f / m_I;
	}

	// Move center of mass.
	b2Vec2 oldCenter = m_sweep.c;
	m_sweep.localCenter =  massData->center;
	m_sweep.c0 = m_sweep.c = b2Mul(m_xf, m_sweep.localCenter);

	// Update center of mass velocity.
	m_linearVelocity += b2Cross(m_angularVelocity, m_sweep.c - oldCenter);
}

bool b2Body::ShouldCollide(const b2Body* other) const
{
	// At least one body should be dynamic.
	if (m_type != b2_dynamicBody && other->m_type != b2_dynamicBody)
	{
		return false;
	}

	// Does a joint prevent collision?
	for (b2JointEdge* jn = m_jointList; jn; jn = jn->next)
	{
		if (jn->other == other)
		{
			if (jn->joint->m_collideConnected == false)
			{
				return false;
			}
		}
	}

	return true;
}

void b2Body::SetTransform(const b2Vec2& position, float32 angle)
{
	b2Assert(m_world->IsLocked() == false);
	if (m_world->IsLocked() == true)
	{
		return;
	}

	m_xf.q.Set(angle);
	m_xf.p = position;
	m_xf0 = m_xf;

	m_sweep.c = b2Mul(m_xf, m_sweep.localCenter);
	m_sweep.a = angle;

	m_sweep.c0 = m_sweep.c;
	m_sweep.a0 = angle;

	b2BroadPhase* broadPhase = &m_world->m_contactManager.m_broadPhase;
	for (b2Fixture* f = m_fixtureList; f; f = f->m_next)
	{
		f->Synchronize(broadPhase, m_xf, m_xf);
	}
}

void b2Body::SynchronizeFixtures()
{
	b2Transform xf1;
	xf1.q.Set(m_sweep.a0);
	xf1.p = m_sweep.c0 - b2Mul(xf1.q, m_sweep.localCenter);

	b2BroadPhase* broadPhase = &m_world->m_contactManager.m_broadPhase;
	for (b2Fixture* f = m_fixtureList; f; f = f->m_next)
	{
		f->Synchronize(broadPhase, xf1, m_xf);
	}
}

void b2Body::SetActive(bool flag)
{
	b2Assert(m_world->IsLocked() == false);

	if (flag == IsActive())
	{
		return;
	}

	if (flag)
	{
		m_flags |= e_activeFlag;

		// Create all proxies.
		b2BroadPhase* broadPhase = &m_world->m_contactManager.m_broadPhase;
		for (b2Fixture* f = m_fixtureList; f; f = f->m_next)
		{
			f->CreateProxies(broadPhase, m_xf);
		}

		// Contacts are created the next time step.
	}
	else
	{
		m_flags &= ~e_activeFlag;

		// Destroy all proxies.
		b2BroadPhase* broadPhase = &m_world->m_contactManager.m_broadPhase;
		for (b2Fixture* f = m_fixtureList; f; f = f->m_next)
		{
			f->DestroyProxies(broadPhase);
		}

		// Destroy the attached contacts.
		b2ContactEdge* ce = m_contactList;
		while (ce)
		{
			b2ContactEdge* ce0 = ce;
			ce = ce->next;
			m_world->m_contactManager.Destroy(ce0->contact);
		}
		m_contactList = NULL;
	}
}

void b2Body::SetFixedRotation(bool flag)
{
	bool status = (m_flags & e_fixedRotationFlag) == e_fixedRotationFlag;
	if (status == flag)
	{
		return;
	}

	if (flag)
	{
		m_flags |= e_fixedRotationFlag;
	}
	else
	{
		m_flags &= ~e_fixedRotationFlag;
	}

	m_angularVelocity = 0.0f;

	ResetMassData();
}

void b2Body::Dump()
{
	int32 bodyIndex = m_islandIndex;

	b2Log("{\n");
	b2Log("  b2BodyDef bd;\n");
	b2Log("  bd.type = b2BodyType(%d);\n", m_type);
	b2Log("  bd.position.Set(%.15lef, %.15lef);\n", m_xf.p.x, m_xf.p.y);
	b2Log("  bd.angle = %.15lef;\n", m_sweep.a);
	b2Log("  bd.linearVelocity.Set(%.15lef, %.15lef);\n", m_linearVelocity.x, m_linearVelocity.y);
	b2Log("  bd.angularVelocity = %.15lef;\n", m_angularVelocity);
	b2Log("  bd.linearDamping = %.15lef;\n", m_linearDamping);
	b2Log("  bd.angularDamping = %.15lef;\n", m_angularDamping);
	b2Log("  bd.allowSleep = bool(%d);\n", m_flags & e_autoSleepFlag);
	b2Log("  bd.awake = bool(%d);\n", m_flags & e_awakeFlag);
	b2Log("  bd.fixedRotation = bool(%d);\n", m_flags & e_fixedRotationFlag);
	b2Log("  bd.bullet = bool(%d);\n", m_flags & e_bulletFlag);
	b2Log("  bd.active = bool(%d);\n", m_flags & e_activeFlag);
	b2Log("  bd.gravityScale = %.15lef;\n", m_gravityScale);
	b2Log("  bodies[%d] = m_world->CreateBody(&bd);\n", m_islandIndex);
	b2Log("\n");
	for (b2Fixture* f = m_fixtureList; f; f = f->m_next)
	{
		b2Log("  {\n");
		f->Dump(bodyIndex);
		b2Log("  }\n");
	}
	b2Log("}\n");
}

// end of Body.cpp

class b2Contact;
class b2Joint;
class b2StackAllocator;
class b2ContactListener;
struct b2ContactVelocityConstraint;
struct b2Profile;

/// This is an internal class.
class b2Island
{
public:
	b2Island(int32 bodyCapacity, int32 contactCapacity, int32 jointCapacity,
			b2StackAllocator* allocator, b2ContactListener* listener);
	~b2Island();

	void Clear()
	{
		m_bodyCount = 0;
		m_contactCount = 0;
		m_jointCount = 0;
	}

	void Solve(b2Profile* profile, const b2TimeStep& step, const b2Vec2& gravity, bool allowSleep);

	void SolveTOI(const b2TimeStep& subStep, int32 toiIndexA, int32 toiIndexB);

	void Add(b2Body* body)
	{
		b2Assert(m_bodyCount < m_bodyCapacity);
		body->m_islandIndex = m_bodyCount;
		m_bodies[m_bodyCount] = body;
		++m_bodyCount;
	}

	void Add(b2Contact* contact)
	{
		b2Assert(m_contactCount < m_contactCapacity);
		m_contacts[m_contactCount++] = contact;
	}

	void Add(b2Joint* joint)
	{
		b2Assert(m_jointCount < m_jointCapacity);
		m_joints[m_jointCount++] = joint;
	}

	void Report(const b2ContactVelocityConstraint* constraints);

	b2StackAllocator* m_allocator;
	b2ContactListener* m_listener;

	b2Body** m_bodies;
	b2Contact** m_contacts;
	b2Joint** m_joints;

	b2Position* m_positions;
	b2Velocity* m_velocities;

	int32 m_bodyCount;
	int32 m_jointCount;
	int32 m_contactCount;

	int32 m_bodyCapacity;
	int32 m_contactCapacity;
	int32 m_jointCapacity;
};

// end of Island.h

class b2Contact;
class b2Body;
class b2StackAllocator;
struct b2ContactPositionConstraint;

struct b2VelocityConstraintPoint
{
	b2Vec2 rA;
	b2Vec2 rB;
	float32 normalImpulse;
	float32 tangentImpulse;
	float32 normalMass;
	float32 tangentMass;
	float32 velocityBias;
};

struct b2ContactVelocityConstraint
{
	b2VelocityConstraintPoint points[b2_maxManifoldPoints];
	b2Vec2 normal;
	b2Mat22 normalMass;
	b2Mat22 K;
	int32 indexA;
	int32 indexB;
	float32 invMassA, invMassB;
	float32 invIA, invIB;
	float32 friction;
	float32 restitution;
	float32 tangentSpeed;
	int32 pointCount;
	int32 contactIndex;
};

struct b2ContactSolverDef
{
	b2TimeStep step;
	b2Contact** contacts;
	int32 count;
	b2Position* positions;
	b2Velocity* velocities;
	b2StackAllocator* allocator;
};

class b2ContactSolver
{
public:
	b2ContactSolver(b2ContactSolverDef* def);
	~b2ContactSolver();

	void InitializeVelocityConstraints();

	void WarmStart();
	void SolveVelocityConstraints();
	void StoreImpulses();

	bool SolvePositionConstraints();
	bool SolveTOIPositionConstraints(int32 toiIndexA, int32 toiIndexB);

	b2TimeStep m_step;
	b2Position* m_positions;
	b2Velocity* m_velocities;
	b2StackAllocator* m_allocator;
	b2ContactPositionConstraint* m_positionConstraints;
	b2ContactVelocityConstraint* m_velocityConstraints;
	b2Contact** m_contacts;
	int m_count;
};

// end of ContactSolver.h

/*
Position Correction Notes
=========================
I tried the several algorithms for position correction of the 2D revolute joint.
I looked at these systems:
- simple pendulum (1m diameter sphere on massless 5m stick) with initial angular velocity of 100 rad/s.
- suspension bridge with 30 1m long planks of length 1m.
- multi-link chain with 30 1m long links.

Here are the algorithms:

Baumgarte - A fraction of the position error is added to the velocity error. There is no
separate position solver.

Pseudo Velocities - After the velocity solver and position integration,
the position error, Jacobian, and effective mass are recomputed. Then
the velocity constraints are solved with pseudo velocities and a fraction
of the position error is added to the pseudo velocity error. The pseudo
velocities are initialized to zero and there is no warm-starting. After
the position solver, the pseudo velocities are added to the positions.
This is also called the First Order World method or the Position LCP method.

Modified Nonlinear Gauss-Seidel (NGS) - Like Pseudo Velocities except the
position error is re-computed for each constraint and the positions are updated
after the constraint is solved. The radius vectors (aka Jacobians) are
re-computed too (otherwise the algorithm has horrible instability). The pseudo
velocity states are not needed because they are effectively zero at the beginning
of each iteration. Since we have the current position error, we allow the
iterations to terminate early if the error becomes smaller than b2_linearSlop.

Full NGS or just NGS - Like Modified NGS except the effective mass are re-computed
each time a constraint is solved.

Here are the results:
Baumgarte - this is the cheapest algorithm but it has some stability problems,
especially with the bridge. The chain links separate easily close to the root
and they jitter as they struggle to pull together. This is one of the most common
methods in the field. The big drawback is that the position correction artificially
affects the momentum, thus leading to instabilities and false bounce. I used a
bias factor of 0.2. A larger bias factor makes the bridge less stable, a smaller
factor makes joints and contacts more spongy.

Pseudo Velocities - the is more stable than the Baumgarte method. The bridge is
stable. However, joints still separate with large angular velocities. Drag the
simple pendulum in a circle quickly and the joint will separate. The chain separates
easily and does not recover. I used a bias factor of 0.2. A larger value lead to
the bridge collapsing when a heavy cube drops on it.

Modified NGS - this algorithm is better in some ways than Baumgarte and Pseudo
Velocities, but in other ways it is worse. The bridge and chain are much more
stable, but the simple pendulum goes unstable at high angular velocities.

Full NGS - stable in all tests. The joints display good stiffness. The bridge
still sags, but this is better than infinite forces.

Recommendations
Pseudo Velocities are not really worthwhile because the bridge and chain cannot
recover from joint separation. In other cases the benefit over Baumgarte is small.

Modified NGS is not a robust method for the revolute joint due to the violent
instability seen in the simple pendulum. Perhaps it is viable with other constraint
types, especially scalar constraints where the effective mass is a scalar.

This leaves Baumgarte and Full NGS. Baumgarte has small, but manageable instabilities
and is very fast. I don't think we can escape Baumgarte, especially in highly
demanding cases where high constraint fidelity is not needed.

Full NGS is robust and easy on the eyes. I recommend this as an option for
higher fidelity simulation and certainly for suspension bridges and long chains.
Full NGS might be a good choice for ragdolls, especially motorized ragdolls where
joint separation can be problematic. The number of NGS iterations can be reduced
for better performance without harming robustness much.

Each joint in a can be handled differently in the position solver. So I recommend
a system where the user can select the algorithm on a per joint basis. I would
probably default to the slower Full NGS and let the user select the faster
Baumgarte method in performance critical scenarios.
*/

/*
Cache Performance

The Box2D solvers are dominated by cache misses. Data structures are designed
to increase the number of cache hits. Much of misses are due to random access
to body data. The constraint structures are iterated over linearly, which leads
to few cache misses.

The bodies are not accessed during iteration. Instead read only data, such as
the mass values are stored with the constraints. The mutable data are the constraint
impulses and the bodies velocities/positions. The impulses are held inside the
constraint structures. The body velocities/positions are held in compact, temporary
arrays to increase the number of cache hits. Linear and angular velocity are
stored in a single array since multiple arrays lead to multiple misses.
*/

/*
2D Rotation

R = [cos(theta) -sin(theta)]
    [sin(theta) cos(theta) ]

thetaDot = omega

Let q1 = cos(theta), q2 = sin(theta).
R = [q1 -q2]
    [q2  q1]

q1Dot = -thetaDot * q2
q2Dot = thetaDot * q1

q1_new = q1_old - dt * w * q2
q2_new = q2_old + dt * w * q1
then normalize.

This might be faster than computing sin+cos.
However, we can compute sin+cos of the same angle fast.
*/

b2Island::b2Island(
	int32 bodyCapacity,
	int32 contactCapacity,
	int32 jointCapacity,
	b2StackAllocator* allocator,
	b2ContactListener* listener)
{
	m_bodyCapacity = bodyCapacity;
	m_contactCapacity = contactCapacity;
	m_jointCapacity	 = jointCapacity;
	m_bodyCount = 0;
	m_contactCount = 0;
	m_jointCount = 0;

	m_allocator = allocator;
	m_listener = listener;

	m_bodies = (b2Body**)m_allocator->Allocate(bodyCapacity * sizeof(b2Body*));
	m_contacts = (b2Contact**)m_allocator->Allocate(contactCapacity	 * sizeof(b2Contact*));
	m_joints = (b2Joint**)m_allocator->Allocate(jointCapacity * sizeof(b2Joint*));

	m_velocities = (b2Velocity*)m_allocator->Allocate(m_bodyCapacity * sizeof(b2Velocity));
	m_positions = (b2Position*)m_allocator->Allocate(m_bodyCapacity * sizeof(b2Position));
}

b2Island::~b2Island()
{
	// Warning: the order should reverse the constructor order.
	m_allocator->Free(m_positions);
	m_allocator->Free(m_velocities);
	m_allocator->Free(m_joints);
	m_allocator->Free(m_contacts);
	m_allocator->Free(m_bodies);
}

void b2Island::Solve(b2Profile* profile, const b2TimeStep& step, const b2Vec2& gravity, bool allowSleep)
{
	b2Timer timer;

	float32 h = step.dt;

	// Integrate velocities and apply damping. Initialize the body state.
	for (int32 i = 0; i < m_bodyCount; ++i)
	{
		b2Body* b = m_bodies[i];

		b2Vec2 c = b->m_sweep.c;
		float32 a = b->m_sweep.a;
		b2Vec2 v = b->m_linearVelocity;
		float32 w = b->m_angularVelocity;

		// Store positions for continuous collision.
		b->m_sweep.c0 = b->m_sweep.c;
		b->m_sweep.a0 = b->m_sweep.a;

		if (b->m_type == b2_dynamicBody)
		{
			// Integrate velocities.
			v += h * (b->m_gravityScale * gravity + b->m_invMass * b->m_force);
			w += h * b->m_invI * b->m_torque;

			// Apply damping.
			// ODE: dv/dt + c * v = 0
			// Solution: v(t) = v0 * exp(-c * t)
			// Time step: v(t + dt) = v0 * exp(-c * (t + dt)) = v0 * exp(-c * t) * exp(-c * dt) = v * exp(-c * dt)
			// v2 = exp(-c * dt) * v1
			// Pade approximation:
			// v2 = v1 * 1 / (1 + c * dt)
			v *= 1.0f / (1.0f + h * b->m_linearDamping);
			w *= 1.0f / (1.0f + h * b->m_angularDamping);
		}

		m_positions[i].c = c;
		m_positions[i].a = a;
		m_velocities[i].v = v;
		m_velocities[i].w = w;
	}

	timer.Reset();

	// Solver data
	b2SolverData solverData;
	solverData.step = step;
	solverData.positions = m_positions;
	solverData.velocities = m_velocities;

	// Initialize velocity constraints.
	b2ContactSolverDef contactSolverDef;
	contactSolverDef.step = step;
	contactSolverDef.contacts = m_contacts;
	contactSolverDef.count = m_contactCount;
	contactSolverDef.positions = m_positions;
	contactSolverDef.velocities = m_velocities;
	contactSolverDef.allocator = m_allocator;

	b2ContactSolver contactSolver(&contactSolverDef);
	contactSolver.InitializeVelocityConstraints();

	if (step.warmStarting)
	{
		contactSolver.WarmStart();
	}
	
	for (int32 i = 0; i < m_jointCount; ++i)
	{
		m_joints[i]->InitVelocityConstraints(solverData);
	}

	profile->solveInit = timer.GetMilliseconds();

	// Solve velocity constraints
	timer.Reset();
	for (int32 i = 0; i < step.velocityIterations; ++i)
	{
		for (int32 j = 0; j < m_jointCount; ++j)
		{
			m_joints[j]->SolveVelocityConstraints(solverData);
		}

		contactSolver.SolveVelocityConstraints();
	}

	// Store impulses for warm starting
	contactSolver.StoreImpulses();
	profile->solveVelocity = timer.GetMilliseconds();

	// Integrate positions
	for (int32 i = 0; i < m_bodyCount; ++i)
	{
		b2Vec2 c = m_positions[i].c;
		float32 a = m_positions[i].a;
		b2Vec2 v = m_velocities[i].v;
		float32 w = m_velocities[i].w;

		// Check for large velocities
		b2Vec2 translation = h * v;
		if (b2Dot(translation, translation) > b2_maxTranslationSquared)
		{
			float32 ratio = b2_maxTranslation / translation.Length();
			v *= ratio;
		}

		float32 rotation = h * w;
		if (rotation * rotation > b2_maxRotationSquared)
		{
			float32 ratio = b2_maxRotation / b2Abs(rotation);
			w *= ratio;
		}

		// Integrate
		c += h * v;
		a += h * w;

		m_positions[i].c = c;
		m_positions[i].a = a;
		m_velocities[i].v = v;
		m_velocities[i].w = w;
	}

	// Solve position constraints
	timer.Reset();
	bool positionSolved = false;
	for (int32 i = 0; i < step.positionIterations; ++i)
	{
		bool contactsOkay = contactSolver.SolvePositionConstraints();

		bool jointsOkay = true;
		for (int32 i = 0; i < m_jointCount; ++i)
		{
			bool jointOkay = m_joints[i]->SolvePositionConstraints(solverData);
			jointsOkay = jointsOkay && jointOkay;
		}

		if (contactsOkay && jointsOkay)
		{
			// Exit early if the position errors are small.
			positionSolved = true;
			break;
		}
	}

	// Copy state buffers back to the bodies
	for (int32 i = 0; i < m_bodyCount; ++i)
	{
		b2Body* body = m_bodies[i];
		body->m_sweep.c = m_positions[i].c;
		body->m_sweep.a = m_positions[i].a;
		body->m_linearVelocity = m_velocities[i].v;
		body->m_angularVelocity = m_velocities[i].w;
		body->SynchronizeTransform();
	}

	profile->solvePosition = timer.GetMilliseconds();

	Report(contactSolver.m_velocityConstraints);

	if (allowSleep)
	{
		float32 minSleepTime = b2_maxFloat;

		const float32 linTolSqr = b2_linearSleepTolerance * b2_linearSleepTolerance;
		const float32 angTolSqr = b2_angularSleepTolerance * b2_angularSleepTolerance;

		for (int32 i = 0; i < m_bodyCount; ++i)
		{
			b2Body* b = m_bodies[i];
			if (b->GetType() == b2_staticBody)
			{
				continue;
			}

			if ((b->m_flags & b2Body::e_autoSleepFlag) == 0 ||
				b->m_angularVelocity * b->m_angularVelocity > angTolSqr ||
				b2Dot(b->m_linearVelocity, b->m_linearVelocity) > linTolSqr)
			{
				b->m_sleepTime = 0.0f;
				minSleepTime = 0.0f;
			}
			else
			{
				b->m_sleepTime += h;
				minSleepTime = b2Min(minSleepTime, b->m_sleepTime);
			}
		}

		if (minSleepTime >= b2_timeToSleep && positionSolved)
		{
			for (int32 i = 0; i < m_bodyCount; ++i)
			{
				b2Body* b = m_bodies[i];
				b->SetAwake(false);
			}
		}
	}
}

void b2Island::SolveTOI(const b2TimeStep& subStep, int32 toiIndexA, int32 toiIndexB)
{
	b2Assert(toiIndexA < m_bodyCount);
	b2Assert(toiIndexB < m_bodyCount);

	// Initialize the body state.
	for (int32 i = 0; i < m_bodyCount; ++i)
	{
		b2Body* b = m_bodies[i];
		m_positions[i].c = b->m_sweep.c;
		m_positions[i].a = b->m_sweep.a;
		m_velocities[i].v = b->m_linearVelocity;
		m_velocities[i].w = b->m_angularVelocity;
	}

	b2ContactSolverDef contactSolverDef;
	contactSolverDef.contacts = m_contacts;
	contactSolverDef.count = m_contactCount;
	contactSolverDef.allocator = m_allocator;
	contactSolverDef.step = subStep;
	contactSolverDef.positions = m_positions;
	contactSolverDef.velocities = m_velocities;
	b2ContactSolver contactSolver(&contactSolverDef);

	// Solve position constraints.
	for (int32 i = 0; i < subStep.positionIterations; ++i)
	{
		bool contactsOkay = contactSolver.SolveTOIPositionConstraints(toiIndexA, toiIndexB);
		if (contactsOkay)
		{
			break;
		}
	}

#if 0
	// Is the new position really safe?
	for (int32 i = 0; i < m_contactCount; ++i)
	{
		b2Contact* c = m_contacts[i];
		b2Fixture* fA = c->GetFixtureA();
		b2Fixture* fB = c->GetFixtureB();

		b2Body* bA = fA->GetBody();
		b2Body* bB = fB->GetBody();

		int32 indexA = c->GetChildIndexA();
		int32 indexB = c->GetChildIndexB();

		b2DistanceInput input;
		input.proxyA.Set(fA->GetShape(), indexA);
		input.proxyB.Set(fB->GetShape(), indexB);
		input.transformA = bA->GetTransform();
		input.transformB = bB->GetTransform();
		input.useRadii = false;

		b2DistanceOutput output;
		b2SimplexCache cache;
		cache.count = 0;
		b2Distance(&output, &cache, &input);

		if (output.distance == 0 || cache.count == 3)
		{
			cache.count += 0;
		}
	}
#endif

	// Leap of faith to new safe state.
	m_bodies[toiIndexA]->m_sweep.c0 = m_positions[toiIndexA].c;
	m_bodies[toiIndexA]->m_sweep.a0 = m_positions[toiIndexA].a;
	m_bodies[toiIndexB]->m_sweep.c0 = m_positions[toiIndexB].c;
	m_bodies[toiIndexB]->m_sweep.a0 = m_positions[toiIndexB].a;

	// No warm starting is needed for TOI events because warm
	// starting impulses were applied in the discrete solver.
	contactSolver.InitializeVelocityConstraints();

	// Solve velocity constraints.
	for (int32 i = 0; i < subStep.velocityIterations; ++i)
	{
		contactSolver.SolveVelocityConstraints();
	}

	// Don't store the TOI contact forces for warm starting
	// because they can be quite large.

	float32 h = subStep.dt;

	// Integrate positions
	for (int32 i = 0; i < m_bodyCount; ++i)
	{
		b2Vec2 c = m_positions[i].c;
		float32 a = m_positions[i].a;
		b2Vec2 v = m_velocities[i].v;
		float32 w = m_velocities[i].w;

		// Check for large velocities
		b2Vec2 translation = h * v;
		if (b2Dot(translation, translation) > b2_maxTranslationSquared)
		{
			float32 ratio = b2_maxTranslation / translation.Length();
			v *= ratio;
		}

		float32 rotation = h * w;
		if (rotation * rotation > b2_maxRotationSquared)
		{
			float32 ratio = b2_maxRotation / b2Abs(rotation);
			w *= ratio;
		}

		// Integrate
		c += h * v;
		a += h * w;

		m_positions[i].c = c;
		m_positions[i].a = a;
		m_velocities[i].v = v;
		m_velocities[i].w = w;

		// Sync bodies
		b2Body* body = m_bodies[i];
		body->m_sweep.c = c;
		body->m_sweep.a = a;
		body->m_linearVelocity = v;
		body->m_angularVelocity = w;
		body->SynchronizeTransform();
	}

	Report(contactSolver.m_velocityConstraints);
}

void b2Island::Report(const b2ContactVelocityConstraint* constraints)
{
	if (m_listener == NULL)
	{
		return;
	}

	for (int32 i = 0; i < m_contactCount; ++i)
	{
		b2Contact* c = m_contacts[i];

		const b2ContactVelocityConstraint* vc = constraints + i;
		
		b2ContactImpulse impulse;
		impulse.count = vc->pointCount;
		for (int32 j = 0; j < vc->pointCount; ++j)
		{
			impulse.normalImpulses[j] = vc->points[j].normalImpulse;
			impulse.tangentImpulses[j] = vc->points[j].tangentImpulse;
		}

		m_listener->PostSolve(c, &impulse);
	}
}

// end of Island.cpp

const float32 b2_minPulleyLength = 2.0f;

/// Pulley joint definition. This requires two ground anchors,
/// two dynamic body anchor points, and a pulley ratio.
struct b2PulleyJointDef : public b2JointDef
{
	b2PulleyJointDef()
	{
		type = e_pulleyJoint;
		groundAnchorA.Set(-1.0f, 1.0f);
		groundAnchorB.Set(1.0f, 1.0f);
		localAnchorA.Set(-1.0f, 0.0f);
		localAnchorB.Set(1.0f, 0.0f);
		lengthA = 0.0f;
		lengthB = 0.0f;
		ratio = 1.0f;
		collideConnected = true;
	}

	/// Initialize the bodies, anchors, lengths, max lengths, and ratio using the world anchors.
	void Initialize(b2Body* bodyA, b2Body* bodyB,
					const b2Vec2& groundAnchorA, const b2Vec2& groundAnchorB,
					const b2Vec2& anchorA, const b2Vec2& anchorB,
					float32 ratio);

	/// The first ground anchor in world coordinates. This point never moves.
	b2Vec2 groundAnchorA;

	/// The second ground anchor in world coordinates. This point never moves.
	b2Vec2 groundAnchorB;

	/// The local anchor point relative to bodyA's origin.
	b2Vec2 localAnchorA;

	/// The local anchor point relative to bodyB's origin.
	b2Vec2 localAnchorB;

	/// The a reference length for the segment attached to bodyA.
	float32 lengthA;

	/// The a reference length for the segment attached to bodyB.
	float32 lengthB;

	/// The pulley ratio, used to simulate a block-and-tackle.
	float32 ratio;
};

/// The pulley joint is connected to two bodies and two fixed ground points.
/// The pulley supports a ratio such that:
/// length1 + ratio * length2 <= constant
/// Yes, the force transmitted is scaled by the ratio.
/// Warning: the pulley joint can get a bit squirrelly by itself. They often
/// work better when combined with prismatic joints. You should also cover the
/// the anchor points with static shapes to prevent one side from going to
/// zero length.
class b2PulleyJoint : public b2Joint
{
public:
	b2Vec2 GetAnchorA() const;
	b2Vec2 GetAnchorB() const;

	b2Vec2 GetReactionForce(float32 inv_dt) const;
	float32 GetReactionTorque(float32 inv_dt) const;

	/// Get the first ground anchor.
	b2Vec2 GetGroundAnchorA() const;

	/// Get the second ground anchor.
	b2Vec2 GetGroundAnchorB() const;

	/// Get the current length of the segment attached to bodyA.
	float32 GetLengthA() const;

	/// Get the current length of the segment attached to bodyB.
	float32 GetLengthB() const;

	/// Get the pulley ratio.
	float32 GetRatio() const;

	/// Get the current length of the segment attached to bodyA.
	float32 GetCurrentLengthA() const;

	/// Get the current length of the segment attached to bodyB.
	float32 GetCurrentLengthB() const;

	/// Dump joint to dmLog
	void Dump();

	/// Implement b2Joint::ShiftOrigin
	void ShiftOrigin(const b2Vec2& newOrigin);

protected:

	friend class b2Joint;
	b2PulleyJoint(const b2PulleyJointDef* data);

	void InitVelocityConstraints(const b2SolverData& data);
	void SolveVelocityConstraints(const b2SolverData& data);
	bool SolvePositionConstraints(const b2SolverData& data);

	b2Vec2 m_groundAnchorA;
	b2Vec2 m_groundAnchorB;
	float32 m_lengthA;
	float32 m_lengthB;
	
	// Solver shared
	b2Vec2 m_localAnchorA;
	b2Vec2 m_localAnchorB;
	float32 m_constant;
	float32 m_ratio;
	float32 m_impulse;

	// Solver temp
	int32 m_indexA;
	int32 m_indexB;
	b2Vec2 m_uA;
	b2Vec2 m_uB;
	b2Vec2 m_rA;
	b2Vec2 m_rB;
	b2Vec2 m_localCenterA;
	b2Vec2 m_localCenterB;
	float32 m_invMassA;
	float32 m_invMassB;
	float32 m_invIA;
	float32 m_invIB;
	float32 m_mass;
};

// end of PulleyJoint.h

// Pulley:
// length1 = norm(p1 - s1)
// length2 = norm(p2 - s2)
// C0 = (length1 + ratio * length2)_initial
// C = C0 - (length1 + ratio * length2)
// u1 = (p1 - s1) / norm(p1 - s1)
// u2 = (p2 - s2) / norm(p2 - s2)
// Cdot = -dot(u1, v1 + cross(w1, r1)) - ratio * dot(u2, v2 + cross(w2, r2))
// J = -[u1 cross(r1, u1) ratio * u2  ratio * cross(r2, u2)]
// K = J * invM * JT
//   = invMass1 + invI1 * cross(r1, u1)^2 + ratio^2 * (invMass2 + invI2 * cross(r2, u2)^2)

void b2PulleyJointDef::Initialize(b2Body* bA, b2Body* bB,
				const b2Vec2& groundA, const b2Vec2& groundB,
				const b2Vec2& anchorA, const b2Vec2& anchorB,
				float32 r)
{
	bodyA = bA;
	bodyB = bB;
	groundAnchorA = groundA;
	groundAnchorB = groundB;
	localAnchorA = bodyA->GetLocalPoint(anchorA);
	localAnchorB = bodyB->GetLocalPoint(anchorB);
	b2Vec2 dA = anchorA - groundA;
	lengthA = dA.Length();
	b2Vec2 dB = anchorB - groundB;
	lengthB = dB.Length();
	ratio = r;
	b2Assert(ratio > b2_epsilon);
}

b2PulleyJoint::b2PulleyJoint(const b2PulleyJointDef* def)
: b2Joint(def)
{
	m_groundAnchorA = def->groundAnchorA;
	m_groundAnchorB = def->groundAnchorB;
	m_localAnchorA = def->localAnchorA;
	m_localAnchorB = def->localAnchorB;

	m_lengthA = def->lengthA;
	m_lengthB = def->lengthB;

	b2Assert(def->ratio != 0.0f);
	m_ratio = def->ratio;

	m_constant = def->lengthA + m_ratio * def->lengthB;

	m_impulse = 0.0f;
}

void b2PulleyJoint::InitVelocityConstraints(const b2SolverData& data)
{
	m_indexA = m_bodyA->m_islandIndex;
	m_indexB = m_bodyB->m_islandIndex;
	m_localCenterA = m_bodyA->m_sweep.localCenter;
	m_localCenterB = m_bodyB->m_sweep.localCenter;
	m_invMassA = m_bodyA->m_invMass;
	m_invMassB = m_bodyB->m_invMass;
	m_invIA = m_bodyA->m_invI;
	m_invIB = m_bodyB->m_invI;

	b2Vec2 cA = data.positions[m_indexA].c;
	float32 aA = data.positions[m_indexA].a;
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;

	b2Vec2 cB = data.positions[m_indexB].c;
	float32 aB = data.positions[m_indexB].a;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	b2Rot qA(aA), qB(aB);

	m_rA = b2Mul(qA, m_localAnchorA - m_localCenterA);
	m_rB = b2Mul(qB, m_localAnchorB - m_localCenterB);

	// Get the pulley axes.
	m_uA = cA + m_rA - m_groundAnchorA;
	m_uB = cB + m_rB - m_groundAnchorB;

	float32 lengthA = m_uA.Length();
	float32 lengthB = m_uB.Length();

	if (lengthA > 10.0f * b2_linearSlop)
	{
		m_uA *= 1.0f / lengthA;
	}
	else
	{
		m_uA.SetZero();
	}

	if (lengthB > 10.0f * b2_linearSlop)
	{
		m_uB *= 1.0f / lengthB;
	}
	else
	{
		m_uB.SetZero();
	}

	// Compute effective mass.
	float32 ruA = b2Cross(m_rA, m_uA);
	float32 ruB = b2Cross(m_rB, m_uB);

	float32 mA = m_invMassA + m_invIA * ruA * ruA;
	float32 mB = m_invMassB + m_invIB * ruB * ruB;

	m_mass = mA + m_ratio * m_ratio * mB;

	if (m_mass > 0.0f)
	{
		m_mass = 1.0f / m_mass;
	}

	if (data.step.warmStarting)
	{
		// Scale impulses to support variable time steps.
		m_impulse *= data.step.dtRatio;

		// Warm starting.
		b2Vec2 PA = -(m_impulse) * m_uA;
		b2Vec2 PB = (-m_ratio * m_impulse) * m_uB;

		vA += m_invMassA * PA;
		wA += m_invIA * b2Cross(m_rA, PA);
		vB += m_invMassB * PB;
		wB += m_invIB * b2Cross(m_rB, PB);
	}
	else
	{
		m_impulse = 0.0f;
	}

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

void b2PulleyJoint::SolveVelocityConstraints(const b2SolverData& data)
{
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	b2Vec2 vpA = vA + b2Cross(wA, m_rA);
	b2Vec2 vpB = vB + b2Cross(wB, m_rB);

	float32 Cdot = -b2Dot(m_uA, vpA) - m_ratio * b2Dot(m_uB, vpB);
	float32 impulse = -m_mass * Cdot;
	m_impulse += impulse;

	b2Vec2 PA = -impulse * m_uA;
	b2Vec2 PB = -m_ratio * impulse * m_uB;
	vA += m_invMassA * PA;
	wA += m_invIA * b2Cross(m_rA, PA);
	vB += m_invMassB * PB;
	wB += m_invIB * b2Cross(m_rB, PB);

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

bool b2PulleyJoint::SolvePositionConstraints(const b2SolverData& data)
{
	b2Vec2 cA = data.positions[m_indexA].c;
	float32 aA = data.positions[m_indexA].a;
	b2Vec2 cB = data.positions[m_indexB].c;
	float32 aB = data.positions[m_indexB].a;

	b2Rot qA(aA), qB(aB);

	b2Vec2 rA = b2Mul(qA, m_localAnchorA - m_localCenterA);
	b2Vec2 rB = b2Mul(qB, m_localAnchorB - m_localCenterB);

	// Get the pulley axes.
	b2Vec2 uA = cA + rA - m_groundAnchorA;
	b2Vec2 uB = cB + rB - m_groundAnchorB;

	float32 lengthA = uA.Length();
	float32 lengthB = uB.Length();

	if (lengthA > 10.0f * b2_linearSlop)
	{
		uA *= 1.0f / lengthA;
	}
	else
	{
		uA.SetZero();
	}

	if (lengthB > 10.0f * b2_linearSlop)
	{
		uB *= 1.0f / lengthB;
	}
	else
	{
		uB.SetZero();
	}

	// Compute effective mass.
	float32 ruA = b2Cross(rA, uA);
	float32 ruB = b2Cross(rB, uB);

	float32 mA = m_invMassA + m_invIA * ruA * ruA;
	float32 mB = m_invMassB + m_invIB * ruB * ruB;

	float32 mass = mA + m_ratio * m_ratio * mB;

	if (mass > 0.0f)
	{
		mass = 1.0f / mass;
	}

	float32 C = m_constant - lengthA - m_ratio * lengthB;
	float32 linearError = b2Abs(C);

	float32 impulse = -mass * C;

	b2Vec2 PA = -impulse * uA;
	b2Vec2 PB = -m_ratio * impulse * uB;

	cA += m_invMassA * PA;
	aA += m_invIA * b2Cross(rA, PA);
	cB += m_invMassB * PB;
	aB += m_invIB * b2Cross(rB, PB);

	data.positions[m_indexA].c = cA;
	data.positions[m_indexA].a = aA;
	data.positions[m_indexB].c = cB;
	data.positions[m_indexB].a = aB;

	return linearError < b2_linearSlop;
}

b2Vec2 b2PulleyJoint::GetAnchorA() const
{
	return m_bodyA->GetWorldPoint(m_localAnchorA);
}

b2Vec2 b2PulleyJoint::GetAnchorB() const
{
	return m_bodyB->GetWorldPoint(m_localAnchorB);
}

b2Vec2 b2PulleyJoint::GetReactionForce(float32 inv_dt) const
{
	b2Vec2 P = m_impulse * m_uB;
	return inv_dt * P;
}

float32 b2PulleyJoint::GetReactionTorque(float32 inv_dt) const
{
	B2_NOT_USED(inv_dt);
	return 0.0f;
}

b2Vec2 b2PulleyJoint::GetGroundAnchorA() const
{
	return m_groundAnchorA;
}

b2Vec2 b2PulleyJoint::GetGroundAnchorB() const
{
	return m_groundAnchorB;
}

float32 b2PulleyJoint::GetLengthA() const
{
	return m_lengthA;
}

float32 b2PulleyJoint::GetLengthB() const
{
	return m_lengthB;
}

float32 b2PulleyJoint::GetRatio() const
{
	return m_ratio;
}

float32 b2PulleyJoint::GetCurrentLengthA() const
{
	b2Vec2 p = m_bodyA->GetWorldPoint(m_localAnchorA);
	b2Vec2 s = m_groundAnchorA;
	b2Vec2 d = p - s;
	return d.Length();
}

float32 b2PulleyJoint::GetCurrentLengthB() const
{
	b2Vec2 p = m_bodyB->GetWorldPoint(m_localAnchorB);
	b2Vec2 s = m_groundAnchorB;
	b2Vec2 d = p - s;
	return d.Length();
}

void b2PulleyJoint::Dump()
{
	int32 indexA = m_bodyA->m_islandIndex;
	int32 indexB = m_bodyB->m_islandIndex;

	b2Log("  b2PulleyJointDef jd;\n");
	b2Log("  jd.bodyA = bodies[%d];\n", indexA);
	b2Log("  jd.bodyB = bodies[%d];\n", indexB);
	b2Log("  jd.collideConnected = bool(%d);\n", m_collideConnected);
	b2Log("  jd.groundAnchorA.Set(%.15lef, %.15lef);\n", m_groundAnchorA.x, m_groundAnchorA.y);
	b2Log("  jd.groundAnchorB.Set(%.15lef, %.15lef);\n", m_groundAnchorB.x, m_groundAnchorB.y);
	b2Log("  jd.localAnchorA.Set(%.15lef, %.15lef);\n", m_localAnchorA.x, m_localAnchorA.y);
	b2Log("  jd.localAnchorB.Set(%.15lef, %.15lef);\n", m_localAnchorB.x, m_localAnchorB.y);
	b2Log("  jd.lengthA = %.15lef;\n", m_lengthA);
	b2Log("  jd.lengthB = %.15lef;\n", m_lengthB);
	b2Log("  jd.ratio = %.15lef;\n", m_ratio);
	b2Log("  joints[%d] = m_world->CreateJoint(&jd);\n", m_index);
}

void b2PulleyJoint::ShiftOrigin(const b2Vec2& newOrigin)
{
	m_groundAnchorA -= newOrigin;
	m_groundAnchorB -= newOrigin;
}

// end of PulleyJoint.cpp

b2World::b2World(const b2Vec2& gravity)
{
	Init(gravity);
}

b2World::~b2World()
{
	// Some shapes allocate using b2Alloc.
	b2Body* b = m_bodyList;
	while (b)
	{
		b2Body* bNext = b->m_next;

		b2Fixture* f = b->m_fixtureList;
		while (f)
		{
			b2Fixture* fNext = f->m_next;
			f->m_proxyCount = 0;
			f->Destroy(&m_blockAllocator);
			f = fNext;
		}

		b = bNext;
	}

	while (m_particleSystemList)
	{
		DestroyParticleSystem(m_particleSystemList);
	}

	// Even though the block allocator frees them for us, for safety,
	// we should ensure that all buffers have been freed.
	b2Assert(m_blockAllocator.GetNumGiantAllocations() == 0);
}

void b2World::SetDestructionListener(b2DestructionListener* listener)
{
	m_destructionListener = listener;
}

void b2World::SetContactFilter(b2ContactFilter* filter)
{
	m_contactManager.m_contactFilter = filter;
}

void b2World::SetContactListener(b2ContactListener* listener)
{
	m_contactManager.m_contactListener = listener;
}

void b2World::SetDebugDraw(b2Draw* debugDraw)
{
	m_debugDraw = debugDraw;
}

b2Body* b2World::CreateBody(const b2BodyDef* def)
{
	b2Assert(IsLocked() == false);
	if (IsLocked())
	{
		return NULL;
	}

	void* mem = m_blockAllocator.Allocate(sizeof(b2Body));
	b2Body* b = new (mem) b2Body(def, this);

	// Add to world doubly linked list.
	b->m_prev = NULL;
	b->m_next = m_bodyList;
	if (m_bodyList)
	{
		m_bodyList->m_prev = b;
	}
	m_bodyList = b;
	++m_bodyCount;

	return b;
}

void b2World::DestroyBody(b2Body* b)
{
	b2Assert(m_bodyCount > 0);
	b2Assert(IsLocked() == false);
	if (IsLocked())
	{
		return;
	}

	// Delete the attached joints.
	b2JointEdge* je = b->m_jointList;
	while (je)
	{
		b2JointEdge* je0 = je;
		je = je->next;

		if (m_destructionListener)
		{
			m_destructionListener->SayGoodbye(je0->joint);
		}

		DestroyJoint(je0->joint);

		b->m_jointList = je;
	}
	b->m_jointList = NULL;

	// Delete the attached contacts.
	b2ContactEdge* ce = b->m_contactList;
	while (ce)
	{
		b2ContactEdge* ce0 = ce;
		ce = ce->next;
		m_contactManager.Destroy(ce0->contact);
	}
	b->m_contactList = NULL;

	// Delete the attached fixtures. This destroys broad-phase proxies.
	b2Fixture* f = b->m_fixtureList;
	while (f)
	{
		b2Fixture* f0 = f;
		f = f->m_next;

		if (m_destructionListener)
		{
			m_destructionListener->SayGoodbye(f0);
		}

		f0->DestroyProxies(&m_contactManager.m_broadPhase);
		f0->Destroy(&m_blockAllocator);
		f0->~b2Fixture();
		m_blockAllocator.Free(f0, sizeof(b2Fixture));

		b->m_fixtureList = f;
		b->m_fixtureCount -= 1;
	}
	b->m_fixtureList = NULL;
	b->m_fixtureCount = 0;

	// Remove world body list.
	if (b->m_prev)
	{
		b->m_prev->m_next = b->m_next;
	}

	if (b->m_next)
	{
		b->m_next->m_prev = b->m_prev;
	}

	if (b == m_bodyList)
	{
		m_bodyList = b->m_next;
	}

	--m_bodyCount;
	b->~b2Body();
	m_blockAllocator.Free(b, sizeof(b2Body));
}

b2Joint* b2World::CreateJoint(const b2JointDef* def)
{
	b2Assert(IsLocked() == false);
	if (IsLocked())
	{
		return NULL;
	}

	b2Joint* j = b2Joint::Create(def, &m_blockAllocator);

	// Connect to the world list.
	j->m_prev = NULL;
	j->m_next = m_jointList;
	if (m_jointList)
	{
		m_jointList->m_prev = j;
	}
	m_jointList = j;
	++m_jointCount;

	// Connect to the bodies' doubly linked lists.
	j->m_edgeA.joint = j;
	j->m_edgeA.other = j->m_bodyB;
	j->m_edgeA.prev = NULL;
	j->m_edgeA.next = j->m_bodyA->m_jointList;
	if (j->m_bodyA->m_jointList) j->m_bodyA->m_jointList->prev = &j->m_edgeA;
	j->m_bodyA->m_jointList = &j->m_edgeA;

	j->m_edgeB.joint = j;
	j->m_edgeB.other = j->m_bodyA;
	j->m_edgeB.prev = NULL;
	j->m_edgeB.next = j->m_bodyB->m_jointList;
	if (j->m_bodyB->m_jointList) j->m_bodyB->m_jointList->prev = &j->m_edgeB;
	j->m_bodyB->m_jointList = &j->m_edgeB;

	b2Body* bodyA = def->bodyA;
	b2Body* bodyB = def->bodyB;

	// If the joint prevents collisions, then flag any contacts for filtering.
	if (def->collideConnected == false)
	{
		b2ContactEdge* edge = bodyB->GetContactList();
		while (edge)
		{
			if (edge->other == bodyA)
			{
				// Flag the contact for filtering at the next time step (where either
				// body is awake).
				edge->contact->FlagForFiltering();
			}

			edge = edge->next;
		}
	}

	// Note: creating a joint doesn't wake the bodies.

	return j;
}

void b2World::DestroyJoint(b2Joint* j)
{
	b2Assert(IsLocked() == false);
	if (IsLocked())
	{
		return;
	}

	bool collideConnected = j->m_collideConnected;

	// Remove from the doubly linked list.
	if (j->m_prev)
	{
		j->m_prev->m_next = j->m_next;
	}

	if (j->m_next)
	{
		j->m_next->m_prev = j->m_prev;
	}

	if (j == m_jointList)
	{
		m_jointList = j->m_next;
	}

	// Disconnect from island graph.
	b2Body* bodyA = j->m_bodyA;
	b2Body* bodyB = j->m_bodyB;

	// Wake up connected bodies.
	bodyA->SetAwake(true);
	bodyB->SetAwake(true);

	// Remove from body 1.
	if (j->m_edgeA.prev)
	{
		j->m_edgeA.prev->next = j->m_edgeA.next;
	}

	if (j->m_edgeA.next)
	{
		j->m_edgeA.next->prev = j->m_edgeA.prev;
	}

	if (&j->m_edgeA == bodyA->m_jointList)
	{
		bodyA->m_jointList = j->m_edgeA.next;
	}

	j->m_edgeA.prev = NULL;
	j->m_edgeA.next = NULL;

	// Remove from body 2
	if (j->m_edgeB.prev)
	{
		j->m_edgeB.prev->next = j->m_edgeB.next;
	}

	if (j->m_edgeB.next)
	{
		j->m_edgeB.next->prev = j->m_edgeB.prev;
	}

	if (&j->m_edgeB == bodyB->m_jointList)
	{
		bodyB->m_jointList = j->m_edgeB.next;
	}

	j->m_edgeB.prev = NULL;
	j->m_edgeB.next = NULL;

	b2Joint::Destroy(j, &m_blockAllocator);

	b2Assert(m_jointCount > 0);
	--m_jointCount;

	// If the joint prevents collisions, then flag any contacts for filtering.
	if (collideConnected == false)
	{
		b2ContactEdge* edge = bodyB->GetContactList();
		while (edge)
		{
			if (edge->other == bodyA)
			{
				// Flag the contact for filtering at the next time step (where either
				// body is awake).
				edge->contact->FlagForFiltering();
			}

			edge = edge->next;
		}
	}
}

b2ParticleSystem* b2World::CreateParticleSystem(const b2ParticleSystemDef* def)
{
	b2Assert(IsLocked() == false);
	if (IsLocked())
	{
		return NULL;
	}

	void* mem = m_blockAllocator.Allocate(sizeof(b2ParticleSystem));
	b2ParticleSystem* p = new (mem) b2ParticleSystem(def, this);

	// Add to world doubly linked list.
	p->m_prev = NULL;
	p->m_next = m_particleSystemList;
	if (m_particleSystemList)
	{
		m_particleSystemList->m_prev = p;
	}
	m_particleSystemList = p;

	return p;
}

void b2World::DestroyParticleSystem(b2ParticleSystem* p)
{
	b2Assert(m_particleSystemList != NULL);
	b2Assert(IsLocked() == false);
	if (IsLocked())
	{
		return;
	}

	// Remove world particleSystem list.
	if (p->m_prev)
	{
		p->m_prev->m_next = p->m_next;
	}

	if (p->m_next)
	{
		p->m_next->m_prev = p->m_prev;
	}

	if (p == m_particleSystemList)
	{
		m_particleSystemList = p->m_next;
	}

	p->~b2ParticleSystem();
	m_blockAllocator.Free(p, sizeof(b2ParticleSystem));
}

//
void b2World::SetAllowSleeping(bool flag)
{
	if (flag == m_allowSleep)
	{
		return;
	}

	m_allowSleep = flag;
	if (m_allowSleep == false)
	{
		for (b2Body* b = m_bodyList; b; b = b->m_next)
		{
			b->SetAwake(true);
		}
	}
}

// Initialize the world with a specified gravity.
void b2World::Init(const b2Vec2& gravity)
{
	m_destructionListener = NULL;
	m_debugDraw = NULL;

	m_bodyList = NULL;
	m_jointList = NULL;
	m_particleSystemList = NULL;

	m_bodyCount = 0;
	m_jointCount = 0;

	m_warmStarting = true;
	m_continuousPhysics = true;
	m_subStepping = false;

	m_stepComplete = true;

	m_allowSleep = true;
	m_gravity = gravity;

	m_flags = e_clearForces;

	m_inv_dt0 = 0.0f;

	m_contactManager.m_allocator = &m_blockAllocator;

	m_liquidFunVersion = &b2_liquidFunVersion;
	m_liquidFunVersionString = b2_liquidFunVersionString;

	memset(&m_profile, 0, sizeof(b2Profile));
}

// Find islands, integrate and solve constraints, solve position constraints
void b2World::Solve(const b2TimeStep& step)
{
	// update previous transforms
	for (b2Body* b = m_bodyList; b; b = b->m_next)
	{
		b->m_xf0 = b->m_xf;
	}

	m_profile.solveInit = 0.0f;
	m_profile.solveVelocity = 0.0f;
	m_profile.solvePosition = 0.0f;

	// Size the island for the worst case.
	b2Island island(m_bodyCount,
					m_contactManager.m_contactCount,
					m_jointCount,
					&m_stackAllocator,
					m_contactManager.m_contactListener);

	// Clear all the island flags.
	for (b2Body* b = m_bodyList; b; b = b->m_next)
	{
		b->m_flags &= ~b2Body::e_islandFlag;
	}
	for (b2Contact* c = m_contactManager.m_contactList; c; c = c->m_next)
	{
		c->m_flags &= ~b2Contact::e_islandFlag;
	}
	for (b2Joint* j = m_jointList; j; j = j->m_next)
	{
		j->m_islandFlag = false;
	}

	// Build and simulate all awake islands.
	int32 stackSize = m_bodyCount;
	b2Body** stack = (b2Body**)m_stackAllocator.Allocate(stackSize * sizeof(b2Body*));
	for (b2Body* seed = m_bodyList; seed; seed = seed->m_next)
	{
		if (seed->m_flags & b2Body::e_islandFlag)
		{
			continue;
		}

		if (seed->IsAwake() == false || seed->IsActive() == false)
		{
			continue;
		}

		// The seed can be dynamic or kinematic.
		if (seed->GetType() == b2_staticBody)
		{
			continue;
		}

		// Reset island and stack.
		island.Clear();
		int32 stackCount = 0;
		stack[stackCount++] = seed;
		seed->m_flags |= b2Body::e_islandFlag;

		// Perform a depth first search (DFS) on the constraint graph.
		while (stackCount > 0)
		{
			// Grab the next body off the stack and add it to the island.
			b2Body* b = stack[--stackCount];
			b2Assert(b->IsActive() == true);
			island.Add(b);

			// Make sure the body is awake.
			b->SetAwake(true);

			// To keep islands as small as possible, we don't
			// propagate islands across static bodies.
			if (b->GetType() == b2_staticBody)
			{
				continue;
			}

			// Search all contacts connected to this body.
			for (b2ContactEdge* ce = b->m_contactList; ce; ce = ce->next)
			{
				b2Contact* contact = ce->contact;

				// Has this contact already been added to an island?
				if (contact->m_flags & b2Contact::e_islandFlag)
				{
					continue;
				}

				// Is this contact solid and touching?
				if (contact->IsEnabled() == false ||
					contact->IsTouching() == false)
				{
					continue;
				}

				// Skip sensors.
				bool sensorA = contact->m_fixtureA->m_isSensor;
				bool sensorB = contact->m_fixtureB->m_isSensor;
				if (sensorA || sensorB)
				{
					continue;
				}

				island.Add(contact);
				contact->m_flags |= b2Contact::e_islandFlag;

				b2Body* other = ce->other;

				// Was the other body already added to this island?
				if (other->m_flags & b2Body::e_islandFlag)
				{
					continue;
				}

				b2Assert(stackCount < stackSize);
				stack[stackCount++] = other;
				other->m_flags |= b2Body::e_islandFlag;
			}

			// Search all joints connect to this body.
			for (b2JointEdge* je = b->m_jointList; je; je = je->next)
			{
				if (je->joint->m_islandFlag == true)
				{
					continue;
				}

				b2Body* other = je->other;

				// Don't simulate joints connected to inactive bodies.
				if (other->IsActive() == false)
				{
					continue;
				}

				island.Add(je->joint);
				je->joint->m_islandFlag = true;

				if (other->m_flags & b2Body::e_islandFlag)
				{
					continue;
				}

				b2Assert(stackCount < stackSize);
				stack[stackCount++] = other;
				other->m_flags |= b2Body::e_islandFlag;
			}
		}

		b2Profile profile;
		island.Solve(&profile, step, m_gravity, m_allowSleep);
		m_profile.solveInit += profile.solveInit;
		m_profile.solveVelocity += profile.solveVelocity;
		m_profile.solvePosition += profile.solvePosition;

		// Post solve cleanup.
		for (int32 i = 0; i < island.m_bodyCount; ++i)
		{
			// Allow static bodies to participate in other islands.
			b2Body* b = island.m_bodies[i];
			if (b->GetType() == b2_staticBody)
			{
				b->m_flags &= ~b2Body::e_islandFlag;
			}
		}
	}

	m_stackAllocator.Free(stack);

	{
		b2Timer timer;
		// Synchronize fixtures, check for out of range bodies.
		for (b2Body* b = m_bodyList; b; b = b->GetNext())
		{
			// If a body was not in an island then it did not move.
			if ((b->m_flags & b2Body::e_islandFlag) == 0)
			{
				continue;
			}

			if (b->GetType() == b2_staticBody)
			{
				continue;
			}

			// Update fixtures (for broad-phase).
			b->SynchronizeFixtures();
		}

		// Look for new contacts.
		m_contactManager.FindNewContacts();
		m_profile.broadphase = timer.GetMilliseconds();
	}
}

// Find TOI contacts and solve them.
void b2World::SolveTOI(const b2TimeStep& step)
{
	b2Island island(2 * b2_maxTOIContacts, b2_maxTOIContacts, 0, &m_stackAllocator, m_contactManager.m_contactListener);

	if (m_stepComplete)
	{
		for (b2Body* b = m_bodyList; b; b = b->m_next)
		{
			b->m_flags &= ~b2Body::e_islandFlag;
			b->m_sweep.alpha0 = 0.0f;
		}

		for (b2Contact* c = m_contactManager.m_contactList; c; c = c->m_next)
		{
			// Invalidate TOI
			c->m_flags &= ~(b2Contact::e_toiFlag | b2Contact::e_islandFlag);
			c->m_toiCount = 0;
			c->m_toi = 1.0f;
		}
	}

	// Find TOI events and solve them.
	for (;;)
	{
		// Find the first TOI.
		b2Contact* minContact = NULL;
		float32 minAlpha = 1.0f;

		for (b2Contact* c = m_contactManager.m_contactList; c; c = c->m_next)
		{
			// Is this contact disabled?
			if (c->IsEnabled() == false)
			{
				continue;
			}

			// Prevent excessive sub-stepping.
			if (c->m_toiCount > b2_maxSubSteps)
			{
				continue;
			}

			float32 alpha = 1.0f;
			if (c->m_flags & b2Contact::e_toiFlag)
			{
				// This contact has a valid cached TOI.
				alpha = c->m_toi;
			}
			else
			{
				b2Fixture* fA = c->GetFixtureA();
				b2Fixture* fB = c->GetFixtureB();

				// Is there a sensor?
				if (fA->IsSensor() || fB->IsSensor())
				{
					continue;
				}

				b2Body* bA = fA->GetBody();
				b2Body* bB = fB->GetBody();

				b2BodyType typeA = bA->m_type;
				b2BodyType typeB = bB->m_type;
				b2Assert(typeA == b2_dynamicBody || typeB == b2_dynamicBody);

				bool activeA = bA->IsAwake() && typeA != b2_staticBody;
				bool activeB = bB->IsAwake() && typeB != b2_staticBody;

				// Is at least one body active (awake and dynamic or kinematic)?
				if (activeA == false && activeB == false)
				{
					continue;
				}

				bool collideA = bA->IsBullet() || typeA != b2_dynamicBody;
				bool collideB = bB->IsBullet() || typeB != b2_dynamicBody;

				// Are these two non-bullet dynamic bodies?
				if (collideA == false && collideB == false)
				{
					continue;
				}

				// Compute the TOI for this contact.
				// Put the sweeps onto the same time interval.
				float32 alpha0 = bA->m_sweep.alpha0;

				if (bA->m_sweep.alpha0 < bB->m_sweep.alpha0)
				{
					alpha0 = bB->m_sweep.alpha0;
					bA->m_sweep.Advance(alpha0);
				}
				else if (bB->m_sweep.alpha0 < bA->m_sweep.alpha0)
				{
					alpha0 = bA->m_sweep.alpha0;
					bB->m_sweep.Advance(alpha0);
				}

				b2Assert(alpha0 < 1.0f);

				int32 indexA = c->GetChildIndexA();
				int32 indexB = c->GetChildIndexB();

				// Compute the time of impact in interval [0, minTOI]
				b2TOIInput input;
				input.proxyA.Set(fA->GetShape(), indexA);
				input.proxyB.Set(fB->GetShape(), indexB);
				input.sweepA = bA->m_sweep;
				input.sweepB = bB->m_sweep;
				input.tMax = 1.0f;

				b2TOIOutput output;
				b2TimeOfImpact(&output, &input);

				// Beta is the fraction of the remaining portion of the .
				float32 beta = output.t;
				if (output.state == b2TOIOutput::e_touching)
				{
					alpha = b2Min(alpha0 + (1.0f - alpha0) * beta, 1.0f);
				}
				else
				{
					alpha = 1.0f;
				}

				c->m_toi = alpha;
				c->m_flags |= b2Contact::e_toiFlag;
			}

			if (alpha < minAlpha)
			{
				// This is the minimum TOI found so far.
				minContact = c;
				minAlpha = alpha;
			}
		}

		if (minContact == NULL || 1.0f - 10.0f * b2_epsilon < minAlpha)
		{
			// No more TOI events. Done!
			m_stepComplete = true;
			break;
		}

		// Advance the bodies to the TOI.
		b2Fixture* fA = minContact->GetFixtureA();
		b2Fixture* fB = minContact->GetFixtureB();
		b2Body* bA = fA->GetBody();
		b2Body* bB = fB->GetBody();

		b2Sweep backup1 = bA->m_sweep;
		b2Sweep backup2 = bB->m_sweep;

		bA->Advance(minAlpha);
		bB->Advance(minAlpha);

		// The TOI contact likely has some new contact points.
		minContact->Update(m_contactManager.m_contactListener);
		minContact->m_flags &= ~b2Contact::e_toiFlag;
		++minContact->m_toiCount;

		// Is the contact solid?
		if (minContact->IsEnabled() == false || minContact->IsTouching() == false)
		{
			// Restore the sweeps.
			minContact->SetEnabled(false);
			bA->m_sweep = backup1;
			bB->m_sweep = backup2;
			bA->SynchronizeTransform();
			bB->SynchronizeTransform();
			continue;
		}

		bA->SetAwake(true);
		bB->SetAwake(true);

		// Build the island
		island.Clear();
		island.Add(bA);
		island.Add(bB);
		island.Add(minContact);

		bA->m_flags |= b2Body::e_islandFlag;
		bB->m_flags |= b2Body::e_islandFlag;
		minContact->m_flags |= b2Contact::e_islandFlag;

		// Get contacts on bodyA and bodyB.
		b2Body* bodies[2] = {bA, bB};
		for (int32 i = 0; i < 2; ++i)
		{
			b2Body* body = bodies[i];
			if (body->m_type == b2_dynamicBody)
			{
				for (b2ContactEdge* ce = body->m_contactList; ce; ce = ce->next)
				{
					if (island.m_bodyCount == island.m_bodyCapacity)
					{
						break;
					}

					if (island.m_contactCount == island.m_contactCapacity)
					{
						break;
					}

					b2Contact* contact = ce->contact;

					// Has this contact already been added to the island?
					if (contact->m_flags & b2Contact::e_islandFlag)
					{
						continue;
					}

					// Only add static, kinematic, or bullet bodies.
					b2Body* other = ce->other;
					if (other->m_type == b2_dynamicBody &&
						body->IsBullet() == false && other->IsBullet() == false)
					{
						continue;
					}

					// Skip sensors.
					bool sensorA = contact->m_fixtureA->m_isSensor;
					bool sensorB = contact->m_fixtureB->m_isSensor;
					if (sensorA || sensorB)
					{
						continue;
					}

					// Tentatively advance the body to the TOI.
					b2Sweep backup = other->m_sweep;
					if ((other->m_flags & b2Body::e_islandFlag) == 0)
					{
						other->Advance(minAlpha);
					}

					// Update the contact points
					contact->Update(m_contactManager.m_contactListener);

					// Was the contact disabled by the user?
					if (contact->IsEnabled() == false)
					{
						other->m_sweep = backup;
						other->SynchronizeTransform();
						continue;
					}

					// Are there contact points?
					if (contact->IsTouching() == false)
					{
						other->m_sweep = backup;
						other->SynchronizeTransform();
						continue;
					}

					// Add the contact to the island
					contact->m_flags |= b2Contact::e_islandFlag;
					island.Add(contact);

					// Has the other body already been added to the island?
					if (other->m_flags & b2Body::e_islandFlag)
					{
						continue;
					}

					// Add the other body to the island.
					other->m_flags |= b2Body::e_islandFlag;

					if (other->m_type != b2_staticBody)
					{
						other->SetAwake(true);
					}

					island.Add(other);
				}
			}
		}

		b2TimeStep subStep;
		subStep.dt = (1.0f - minAlpha) * step.dt;
		subStep.inv_dt = 1.0f / subStep.dt;
		subStep.dtRatio = 1.0f;
		subStep.positionIterations = 20;
		subStep.velocityIterations = step.velocityIterations;
		subStep.particleIterations = step.particleIterations;
		subStep.warmStarting = false;
		island.SolveTOI(subStep, bA->m_islandIndex, bB->m_islandIndex);

		// Reset island flags and synchronize broad-phase proxies.
		for (int32 i = 0; i < island.m_bodyCount; ++i)
		{
			b2Body* body = island.m_bodies[i];
			body->m_flags &= ~b2Body::e_islandFlag;

			if (body->m_type != b2_dynamicBody)
			{
				continue;
			}

			body->SynchronizeFixtures();

			// Invalidate all contact TOIs on this displaced body.
			for (b2ContactEdge* ce = body->m_contactList; ce; ce = ce->next)
			{
				ce->contact->m_flags &= ~(b2Contact::e_toiFlag | b2Contact::e_islandFlag);
			}
		}

		// Commit fixture proxy movements to the broad-phase so that new contacts are created.
		// Also, some contacts can be destroyed.
		m_contactManager.FindNewContacts();

		if (m_subStepping)
		{
			m_stepComplete = false;
			break;
		}
	}
}

void b2World::Step(
	float32 dt,
	int32 velocityIterations,
	int32 positionIterations,
	int32 particleIterations)
{
	b2Timer stepTimer;

	// If new fixtures were added, we need to find the new contacts.
	if (m_flags & e_newFixture)
	{
		m_contactManager.FindNewContacts();
		m_flags &= ~e_newFixture;
	}

	m_flags |= e_locked;

	b2TimeStep step;
	step.dt = dt;
	step.velocityIterations	= velocityIterations;
	step.positionIterations = positionIterations;
	step.particleIterations = particleIterations;
	if (dt > 0.0f)
	{
		step.inv_dt = 1.0f / dt;
	}
	else
	{
		step.inv_dt = 0.0f;
	}

	step.dtRatio = m_inv_dt0 * dt;

	step.warmStarting = m_warmStarting;

	// Update contacts. This is where some contacts are destroyed.
	{
		b2Timer timer;
		m_contactManager.Collide();
		m_profile.collide = timer.GetMilliseconds();
	}

	// Integrate velocities, solve velocity constraints, and integrate positions.
	if (m_stepComplete && step.dt > 0.0f)
	{
		b2Timer timer;
		for (b2ParticleSystem* p = m_particleSystemList; p; p = p->GetNext())
		{
			p->Solve(step); // Particle Simulation
		}
		Solve(step);
		m_profile.solve = timer.GetMilliseconds();
	}

	// Handle TOI events.
	if (m_continuousPhysics && step.dt > 0.0f)
	{
		b2Timer timer;
		SolveTOI(step);
		m_profile.solveTOI = timer.GetMilliseconds();
	}

	if (step.dt > 0.0f)
	{
		m_inv_dt0 = step.inv_dt;
	}

	if (m_flags & e_clearForces)
	{
		ClearForces();
	}

	m_flags &= ~e_locked;

	m_profile.step = stepTimer.GetMilliseconds();
}

void b2World::ClearForces()
{
	for (b2Body* body = m_bodyList; body; body = body->GetNext())
	{
		body->m_force.SetZero();
		body->m_torque = 0.0f;
	}
}

struct b2WorldQueryWrapper
{
	bool QueryCallback(int32 proxyId)
	{
		b2FixtureProxy* proxy = (b2FixtureProxy*)broadPhase->GetUserData(proxyId);
		return callback->ReportFixture(proxy->fixture);
	}

	const b2BroadPhase* broadPhase;
	b2QueryCallback* callback;
};

void b2World::QueryAABB(b2QueryCallback* callback, const b2AABB& aabb) const
{
	b2WorldQueryWrapper wrapper;
	wrapper.broadPhase = &m_contactManager.m_broadPhase;
	wrapper.callback = callback;
	m_contactManager.m_broadPhase.Query(&wrapper, aabb);
	for (b2ParticleSystem* p = m_particleSystemList; p; p = p->GetNext())
	{
		if (callback->ShouldQueryParticleSystem(p))
		{
			p->QueryAABB(callback, aabb);
		}
	}
}

void b2World::QueryShapeAABB(b2QueryCallback* callback, const b2Shape& shape,
                             const b2Transform& xf) const
{
	b2AABB aabb;
	shape.ComputeAABB(&aabb, xf, 0);
	QueryAABB(callback, aabb);
}

struct b2WorldRayCastWrapper
{
	float32 RayCastCallback(const b2RayCastInput& input, int32 proxyId)
	{
		void* userData = broadPhase->GetUserData(proxyId);
		b2FixtureProxy* proxy = (b2FixtureProxy*)userData;
		b2Fixture* fixture = proxy->fixture;
		int32 index = proxy->childIndex;
		b2RayCastOutput output;
		bool hit = fixture->RayCast(&output, input, index);

		if (hit)
		{
			float32 fraction = output.fraction;
			b2Vec2 point = (1.0f - fraction) * input.p1 + fraction * input.p2;
			return callback->ReportFixture(fixture, point, output.normal, fraction);
		}

		return input.maxFraction;
	}

	const b2BroadPhase* broadPhase;
	b2RayCastCallback* callback;
};

void b2World::RayCast(b2RayCastCallback* callback, const b2Vec2& point1, const b2Vec2& point2) const
{
	b2WorldRayCastWrapper wrapper;
	wrapper.broadPhase = &m_contactManager.m_broadPhase;
	wrapper.callback = callback;
	b2RayCastInput input;
	input.maxFraction = 1.0f;
	input.p1 = point1;
	input.p2 = point2;
	m_contactManager.m_broadPhase.RayCast(&wrapper, input);
	for (b2ParticleSystem* p = m_particleSystemList; p; p = p->GetNext())
	{
		if (callback->ShouldQueryParticleSystem(p))
		{
			p->RayCast(callback, point1, point2);
		}
	}
}

void b2World::DrawShape(b2Fixture* fixture, const b2Transform& xf, const b2Color& color)
{
	switch (fixture->GetType())
	{
	case b2Shape::e_circle:
		{
			b2CircleShape* circle = (b2CircleShape*)fixture->GetShape();

			b2Vec2 center = b2Mul(xf, circle->m_p);
			float32 radius = circle->m_radius;
			b2Vec2 axis = b2Mul(xf.q, b2Vec2(1.0f, 0.0f));

			m_debugDraw->DrawSolidCircle(center, radius, axis, color);
		}
		break;

	case b2Shape::e_edge:
		{
			b2EdgeShape* edge = (b2EdgeShape*)fixture->GetShape();
			b2Vec2 v1 = b2Mul(xf, edge->m_vertex1);
			b2Vec2 v2 = b2Mul(xf, edge->m_vertex2);
			m_debugDraw->DrawSegment(v1, v2, color);
		}
		break;

	case b2Shape::e_chain:
		{
			b2ChainShape* chain = (b2ChainShape*)fixture->GetShape();
			int32 count = chain->m_count;
			const b2Vec2* vertices = chain->m_vertices;

			b2Vec2 v1 = b2Mul(xf, vertices[0]);
			for (int32 i = 1; i < count; ++i)
			{
				b2Vec2 v2 = b2Mul(xf, vertices[i]);
				m_debugDraw->DrawSegment(v1, v2, color);
				m_debugDraw->DrawCircle(v1, 0.05f, color);
				v1 = v2;
			}
		}
		break;

	case b2Shape::e_polygon:
		{
			b2PolygonShape* poly = (b2PolygonShape*)fixture->GetShape();
			int32 vertexCount = poly->m_count;
			b2Assert(vertexCount <= b2_maxPolygonVertices);
			b2Vec2 vertices[b2_maxPolygonVertices];

			for (int32 i = 0; i < vertexCount; ++i)
			{
				vertices[i] = b2Mul(xf, poly->m_vertices[i]);
			}

			m_debugDraw->DrawSolidPolygon(vertices, vertexCount, color);
		}
		break;

	default:
		break;
	}
}

void b2World::DrawJoint(b2Joint* joint)
{
	b2Body* bodyA = joint->GetBodyA();
	b2Body* bodyB = joint->GetBodyB();
	const b2Transform& xf1 = bodyA->GetTransform();
	const b2Transform& xf2 = bodyB->GetTransform();
	b2Vec2 x1 = xf1.p;
	b2Vec2 x2 = xf2.p;
	b2Vec2 p1 = joint->GetAnchorA();
	b2Vec2 p2 = joint->GetAnchorB();

	b2Color color(0.5f, 0.8f, 0.8f);

	switch (joint->GetType())
	{
	case e_distanceJoint:
		m_debugDraw->DrawSegment(p1, p2, color);
		break;

	case e_pulleyJoint:
		{
			b2PulleyJoint* pulley = (b2PulleyJoint*)joint;
			b2Vec2 s1 = pulley->GetGroundAnchorA();
			b2Vec2 s2 = pulley->GetGroundAnchorB();
			m_debugDraw->DrawSegment(s1, p1, color);
			m_debugDraw->DrawSegment(s2, p2, color);
			m_debugDraw->DrawSegment(s1, s2, color);
		}
		break;

	case e_mouseJoint:
		// don't draw this
		break;

	default:
		m_debugDraw->DrawSegment(x1, p1, color);
		m_debugDraw->DrawSegment(p1, p2, color);
		m_debugDraw->DrawSegment(x2, p2, color);
	}
}

void b2World::DrawParticleSystem(const b2ParticleSystem& system)
{
	int32 particleCount = system.GetParticleCount();
	if (particleCount)
	{
		float32 radius = system.GetRadius();
		const b2Vec2* positionBuffer = system.GetPositionBuffer();
		if (system.m_colorBuffer.data)
		{
			const b2ParticleColor* colorBuffer = system.GetColorBuffer();
			m_debugDraw->DrawParticles(positionBuffer, radius, colorBuffer, particleCount);
		}
		else
		{
			m_debugDraw->DrawParticles(positionBuffer, radius, NULL, particleCount);
		}
	}
}

void b2World::DrawDebugData()
{
	if (m_debugDraw == NULL)
	{
		return;
	}

	uint32 flags = m_debugDraw->GetFlags();

	if (flags & b2Draw::e_shapeBit)
	{
		for (b2Body* b = m_bodyList; b; b = b->GetNext())
		{
			const b2Transform& xf = b->GetTransform();
			for (b2Fixture* f = b->GetFixtureList(); f; f = f->GetNext())
			{
				if (b->IsActive() == false)
				{
					DrawShape(f, xf, b2Color(0.5f, 0.5f, 0.3f));
				}
				else if (b->GetType() == b2_staticBody)
				{
					DrawShape(f, xf, b2Color(0.5f, 0.9f, 0.5f));
				}
				else if (b->GetType() == b2_kinematicBody)
				{
					DrawShape(f, xf, b2Color(0.5f, 0.5f, 0.9f));
				}
				else if (b->IsAwake() == false)
				{
					DrawShape(f, xf, b2Color(0.6f, 0.6f, 0.6f));
				}
				else
				{
					DrawShape(f, xf, b2Color(0.9f, 0.7f, 0.7f));
				}
			}
		}
	}

	if (flags & b2Draw::e_particleBit)
	{
		for (b2ParticleSystem* p = m_particleSystemList; p; p = p->GetNext())
		{
			DrawParticleSystem(*p);
		}
	}

	if (flags & b2Draw::e_jointBit)
	{
		for (b2Joint* j = m_jointList; j; j = j->GetNext())
		{

			DrawJoint(j);
		}
	}

	if (flags & b2Draw::e_pairBit)
	{
		b2Color color(0.3f, 0.9f, 0.9f);
		for (b2Contact* c = m_contactManager.m_contactList; c; c = c->GetNext())
		{
			//b2Fixture* fixtureA = c->GetFixtureA();
			//b2Fixture* fixtureB = c->GetFixtureB();

			//b2Vec2 cA = fixtureA->GetAABB().GetCenter();
			//b2Vec2 cB = fixtureB->GetAABB().GetCenter();

			//m_debugDraw->DrawSegment(cA, cB, color);
		}
	}

	if (flags & b2Draw::e_aabbBit)
	{
		b2Color color(0.9f, 0.3f, 0.9f);
		b2BroadPhase* bp = &m_contactManager.m_broadPhase;

		for (b2Body* b = m_bodyList; b; b = b->GetNext())
		{
			if (b->IsActive() == false)
			{
				continue;
			}

			for (b2Fixture* f = b->GetFixtureList(); f; f = f->GetNext())
			{
				for (int32 i = 0; i < f->m_proxyCount; ++i)
				{
					b2FixtureProxy* proxy = f->m_proxies + i;
					b2AABB aabb = bp->GetFatAABB(proxy->proxyId);
					b2Vec2 vs[4];
					vs[0].Set(aabb.lowerBound.x, aabb.lowerBound.y);
					vs[1].Set(aabb.upperBound.x, aabb.lowerBound.y);
					vs[2].Set(aabb.upperBound.x, aabb.upperBound.y);
					vs[3].Set(aabb.lowerBound.x, aabb.upperBound.y);

					m_debugDraw->DrawPolygon(vs, 4, color);
				}
			}
		}
	}

	if (flags & b2Draw::e_centerOfMassBit)
	{
		for (b2Body* b = m_bodyList; b; b = b->GetNext())
		{
			b2Transform xf = b->GetTransform();
			xf.p = b->GetWorldCenter();
			m_debugDraw->DrawTransform(xf);
		}
	}
}

static float32 GetSmallestRadius(const b2World* world)
{
	float32 smallestRadius = b2_maxFloat;
	for (const b2ParticleSystem* system = world->GetParticleSystemList();
		 system != NULL;
		 system = system->GetNext())
	{
		smallestRadius = b2Min(smallestRadius, system->GetRadius());
	}
	return smallestRadius;
}

int b2World::CalculateReasonableParticleIterations(float32 timeStep) const
{
	if (m_particleSystemList == NULL)
		return 1;

	// Use the smallest radius, since that represents the worst-case.
	return b2CalculateParticleIterations(m_gravity.Length(),
										 GetSmallestRadius(this),
										 timeStep);
}

int32 b2World::GetProxyCount() const
{
	return m_contactManager.m_broadPhase.GetProxyCount();
}

int32 b2World::GetTreeHeight() const
{
	return m_contactManager.m_broadPhase.GetTreeHeight();
}

int32 b2World::GetTreeBalance() const
{
	return m_contactManager.m_broadPhase.GetTreeBalance();
}

float32 b2World::GetTreeQuality() const
{
	return m_contactManager.m_broadPhase.GetTreeQuality();
}

void b2World::ShiftOrigin(const b2Vec2& newOrigin)
{
	b2Assert((m_flags & e_locked) == 0);
	if ((m_flags & e_locked) == e_locked)
	{
		return;
	}

	for (b2Body* b = m_bodyList; b; b = b->m_next)
	{
		b->m_xf.p -= newOrigin;
		b->m_sweep.c0 -= newOrigin;
		b->m_sweep.c -= newOrigin;
	}

	for (b2Joint* j = m_jointList; j; j = j->m_next)
	{
		j->ShiftOrigin(newOrigin);
	}

	m_contactManager.m_broadPhase.ShiftOrigin(newOrigin);
}

void b2World::Dump()
{
	if ((m_flags & e_locked) == e_locked)
	{
		return;
	}

	b2Log("b2Vec2 g(%.15lef, %.15lef);\n", m_gravity.x, m_gravity.y);
	b2Log("m_world->SetGravity(g);\n");

	b2Log("b2Body** bodies = (b2Body**)b2Alloc(%d * sizeof(b2Body*));\n", m_bodyCount);
	b2Log("b2Joint** joints = (b2Joint**)b2Alloc(%d * sizeof(b2Joint*));\n", m_jointCount);
	int32 i = 0;
	for (b2Body* b = m_bodyList; b; b = b->m_next)
	{
		b->m_islandIndex = i;
		b->Dump();
		++i;
	}

	i = 0;
	for (b2Joint* j = m_jointList; j; j = j->m_next)
	{
		j->m_index = i;
		++i;
	}

	// First pass on joints, skip gear joints.
	for (b2Joint* j = m_jointList; j; j = j->m_next)
	{
		if (j->m_type == e_gearJoint)
		{
			continue;
		}

		b2Log("{\n");
		j->Dump();
		b2Log("}\n");
	}

	// Second pass on joints, only gear joints.
	for (b2Joint* j = m_jointList; j; j = j->m_next)
	{
		if (j->m_type != e_gearJoint)
		{
			continue;
		}

		b2Log("{\n");
		j->Dump();
		b2Log("}\n");
	}

	b2Log("b2Free(joints);\n");
	b2Log("b2Free(bodies);\n");
	b2Log("joints = NULL;\n");
	b2Log("bodies = NULL;\n");
}

// end of World.cpp

/// Distance joint definition. This requires defining an
/// anchor point on both bodies and the non-zero length of the
/// distance joint. The definition uses local anchor points
/// so that the initial configuration can violate the constraint
/// slightly. This helps when saving and loading a game.
/// @warning Do not use a zero or short length.
struct b2DistanceJointDef : public b2JointDef
{
	b2DistanceJointDef()
	{
		type = e_distanceJoint;
		localAnchorA.Set(0.0f, 0.0f);
		localAnchorB.Set(0.0f, 0.0f);
		length = 1.0f;
		frequencyHz = 0.0f;
		dampingRatio = 0.0f;
	}

	/// Initialize the bodies, anchors, and length using the world
	/// anchors.
	void Initialize(b2Body* bodyA, b2Body* bodyB,
					const b2Vec2& anchorA, const b2Vec2& anchorB);

	/// The local anchor point relative to bodyA's origin.
	b2Vec2 localAnchorA;

	/// The local anchor point relative to bodyB's origin.
	b2Vec2 localAnchorB;

	/// The natural length between the anchor points.
	float32 length;

	/// The mass-spring-damper frequency in Hertz. A value of 0
	/// disables softness.
	float32 frequencyHz;

	/// The damping ratio. 0 = no damping, 1 = critical damping.
	float32 dampingRatio;
};

/// A distance joint constrains two points on two bodies
/// to remain at a fixed distance from each other. You can view
/// this as a massless, rigid rod.
class b2DistanceJoint : public b2Joint
{
public:

	b2Vec2 GetAnchorA() const;
	b2Vec2 GetAnchorB() const;

	/// Get the reaction force given the inverse time step.
	/// Unit is N.
	b2Vec2 GetReactionForce(float32 inv_dt) const;

	/// Get the reaction torque given the inverse time step.
	/// Unit is N*m. This is always zero for a distance joint.
	float32 GetReactionTorque(float32 inv_dt) const;

	/// The local anchor point relative to bodyA's origin.
	const b2Vec2& GetLocalAnchorA() const { return m_localAnchorA; }

	/// The local anchor point relative to bodyB's origin.
	const b2Vec2& GetLocalAnchorB() const  { return m_localAnchorB; }

	/// Set/get the natural length.
	/// Manipulating the length can lead to non-physical behavior when the frequency is zero.
	void SetLength(float32 length);
	float32 GetLength() const;

	/// Set/get frequency in Hz.
	void SetFrequency(float32 hz);
	float32 GetFrequency() const;

	/// Set/get damping ratio.
	void SetDampingRatio(float32 ratio);
	float32 GetDampingRatio() const;

	/// Dump joint to dmLog
	void Dump();

protected:

	friend class b2Joint;
	b2DistanceJoint(const b2DistanceJointDef* data);

	void InitVelocityConstraints(const b2SolverData& data);
	void SolveVelocityConstraints(const b2SolverData& data);
	bool SolvePositionConstraints(const b2SolverData& data);

	float32 m_frequencyHz;
	float32 m_dampingRatio;
	float32 m_bias;

	// Solver shared
	b2Vec2 m_localAnchorA;
	b2Vec2 m_localAnchorB;
	float32 m_gamma;
	float32 m_impulse;
	float32 m_length;

	// Solver temp
	int32 m_indexA;
	int32 m_indexB;
	b2Vec2 m_u;
	b2Vec2 m_rA;
	b2Vec2 m_rB;
	b2Vec2 m_localCenterA;
	b2Vec2 m_localCenterB;
	float32 m_invMassA;
	float32 m_invMassB;
	float32 m_invIA;
	float32 m_invIB;
	float32 m_mass;
};

inline void b2DistanceJoint::SetLength(float32 length)
{
	m_length = length;
}

inline float32 b2DistanceJoint::GetLength() const
{
	return m_length;
}

inline void b2DistanceJoint::SetFrequency(float32 hz)
{
	m_frequencyHz = hz;
}

inline float32 b2DistanceJoint::GetFrequency() const
{
	return m_frequencyHz;
}

inline void b2DistanceJoint::SetDampingRatio(float32 ratio)
{
	m_dampingRatio = ratio;
}

inline float32 b2DistanceJoint::GetDampingRatio() const
{
	return m_dampingRatio;
}

// end of DistanceJoint.h

// 1-D constrained system
// m (v2 - v1) = lambda
// v2 + (beta/h) * x1 + gamma * lambda = 0, gamma has units of inverse mass.
// x2 = x1 + h * v2

// 1-D mass-damper-spring system
// m (v2 - v1) + h * d * v2 + h * k * 

// C = norm(p2 - p1) - L
// u = (p2 - p1) / norm(p2 - p1)
// Cdot = dot(u, v2 + cross(w2, r2) - v1 - cross(w1, r1))
// J = [-u -cross(r1, u) u cross(r2, u)]
// K = J * invM * JT
//   = invMass1 + invI1 * cross(r1, u)^2 + invMass2 + invI2 * cross(r2, u)^2

void b2DistanceJointDef::Initialize(b2Body* b1, b2Body* b2,
									const b2Vec2& anchor1, const b2Vec2& anchor2)
{
	bodyA = b1;
	bodyB = b2;
	localAnchorA = bodyA->GetLocalPoint(anchor1);
	localAnchorB = bodyB->GetLocalPoint(anchor2);
	b2Vec2 d = anchor2 - anchor1;
	length = d.Length();
}

b2DistanceJoint::b2DistanceJoint(const b2DistanceJointDef* def)
: b2Joint(def)
{
	m_localAnchorA = def->localAnchorA;
	m_localAnchorB = def->localAnchorB;
	m_length = def->length;
	m_frequencyHz = def->frequencyHz;
	m_dampingRatio = def->dampingRatio;
	m_impulse = 0.0f;
	m_gamma = 0.0f;
	m_bias = 0.0f;
}

void b2DistanceJoint::InitVelocityConstraints(const b2SolverData& data)
{
	m_indexA = m_bodyA->m_islandIndex;
	m_indexB = m_bodyB->m_islandIndex;
	m_localCenterA = m_bodyA->m_sweep.localCenter;
	m_localCenterB = m_bodyB->m_sweep.localCenter;
	m_invMassA = m_bodyA->m_invMass;
	m_invMassB = m_bodyB->m_invMass;
	m_invIA = m_bodyA->m_invI;
	m_invIB = m_bodyB->m_invI;

	b2Vec2 cA = data.positions[m_indexA].c;
	float32 aA = data.positions[m_indexA].a;
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;

	b2Vec2 cB = data.positions[m_indexB].c;
	float32 aB = data.positions[m_indexB].a;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	b2Rot qA(aA), qB(aB);

	m_rA = b2Mul(qA, m_localAnchorA - m_localCenterA);
	m_rB = b2Mul(qB, m_localAnchorB - m_localCenterB);
	m_u = cB + m_rB - cA - m_rA;

	// Handle singularity.
	float32 length = m_u.Length();
	if (length > b2_linearSlop)
	{
		m_u *= 1.0f / length;
	}
	else
	{
		m_u.Set(0.0f, 0.0f);
	}

	float32 crAu = b2Cross(m_rA, m_u);
	float32 crBu = b2Cross(m_rB, m_u);
	float32 invMass = m_invMassA + m_invIA * crAu * crAu + m_invMassB + m_invIB * crBu * crBu;

	// Compute the effective mass matrix.
	m_mass = invMass != 0.0f ? 1.0f / invMass : 0.0f;

	if (m_frequencyHz > 0.0f)
	{
		float32 C = length - m_length;

		// Frequency
		float32 omega = 2.0f * b2_pi * m_frequencyHz;

		// Damping coefficient
		float32 d = 2.0f * m_mass * m_dampingRatio * omega;

		// Spring stiffness
		float32 k = m_mass * omega * omega;

		// magic formulas
		float32 h = data.step.dt;
		m_gamma = h * (d + h * k);
		m_gamma = m_gamma != 0.0f ? 1.0f / m_gamma : 0.0f;
		m_bias = C * h * k * m_gamma;

		invMass += m_gamma;
		m_mass = invMass != 0.0f ? 1.0f / invMass : 0.0f;
	}
	else
	{
		m_gamma = 0.0f;
		m_bias = 0.0f;
	}

	if (data.step.warmStarting)
	{
		// Scale the impulse to support a variable time step.
		m_impulse *= data.step.dtRatio;

		b2Vec2 P = m_impulse * m_u;
		vA -= m_invMassA * P;
		wA -= m_invIA * b2Cross(m_rA, P);
		vB += m_invMassB * P;
		wB += m_invIB * b2Cross(m_rB, P);
	}
	else
	{
		m_impulse = 0.0f;
	}

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

void b2DistanceJoint::SolveVelocityConstraints(const b2SolverData& data)
{
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	// Cdot = dot(u, v + cross(w, r))
	b2Vec2 vpA = vA + b2Cross(wA, m_rA);
	b2Vec2 vpB = vB + b2Cross(wB, m_rB);
	float32 Cdot = b2Dot(m_u, vpB - vpA);

	float32 impulse = -m_mass * (Cdot + m_bias + m_gamma * m_impulse);
	m_impulse += impulse;

	b2Vec2 P = impulse * m_u;
	vA -= m_invMassA * P;
	wA -= m_invIA * b2Cross(m_rA, P);
	vB += m_invMassB * P;
	wB += m_invIB * b2Cross(m_rB, P);

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

bool b2DistanceJoint::SolvePositionConstraints(const b2SolverData& data)
{
	if (m_frequencyHz > 0.0f)
	{
		// There is no position correction for soft distance constraints.
		return true;
	}

	b2Vec2 cA = data.positions[m_indexA].c;
	float32 aA = data.positions[m_indexA].a;
	b2Vec2 cB = data.positions[m_indexB].c;
	float32 aB = data.positions[m_indexB].a;

	b2Rot qA(aA), qB(aB);

	b2Vec2 rA = b2Mul(qA, m_localAnchorA - m_localCenterA);
	b2Vec2 rB = b2Mul(qB, m_localAnchorB - m_localCenterB);
	b2Vec2 u = cB + rB - cA - rA;

	float32 length = u.Normalize();
	float32 C = length - m_length;
	C = b2Clamp(C, -b2_maxLinearCorrection, b2_maxLinearCorrection);

	float32 impulse = -m_mass * C;
	b2Vec2 P = impulse * u;

	cA -= m_invMassA * P;
	aA -= m_invIA * b2Cross(rA, P);
	cB += m_invMassB * P;
	aB += m_invIB * b2Cross(rB, P);

	data.positions[m_indexA].c = cA;
	data.positions[m_indexA].a = aA;
	data.positions[m_indexB].c = cB;
	data.positions[m_indexB].a = aB;

	return b2Abs(C) < b2_linearSlop;
}

b2Vec2 b2DistanceJoint::GetAnchorA() const
{
	return m_bodyA->GetWorldPoint(m_localAnchorA);
}

b2Vec2 b2DistanceJoint::GetAnchorB() const
{
	return m_bodyB->GetWorldPoint(m_localAnchorB);
}

b2Vec2 b2DistanceJoint::GetReactionForce(float32 inv_dt) const
{
	b2Vec2 F = (inv_dt * m_impulse) * m_u;
	return F;
}

float32 b2DistanceJoint::GetReactionTorque(float32 inv_dt) const
{
	B2_NOT_USED(inv_dt);
	return 0.0f;
}

void b2DistanceJoint::Dump()
{
	int32 indexA = m_bodyA->m_islandIndex;
	int32 indexB = m_bodyB->m_islandIndex;

	b2Log("  b2DistanceJointDef jd;\n");
	b2Log("  jd.bodyA = bodies[%d];\n", indexA);
	b2Log("  jd.bodyB = bodies[%d];\n", indexB);
	b2Log("  jd.collideConnected = bool(%d);\n", m_collideConnected);
	b2Log("  jd.localAnchorA.Set(%.15lef, %.15lef);\n", m_localAnchorA.x, m_localAnchorA.y);
	b2Log("  jd.localAnchorB.Set(%.15lef, %.15lef);\n", m_localAnchorB.x, m_localAnchorB.y);
	b2Log("  jd.length = %.15lef;\n", m_length);
	b2Log("  jd.frequencyHz = %.15lef;\n", m_frequencyHz);
	b2Log("  jd.dampingRatio = %.15lef;\n", m_dampingRatio);
	b2Log("  joints[%d] = m_world->CreateJoint(&jd);\n", m_index);
}

// end of DistanceJoint.cpp

/// Friction joint definition.
struct b2FrictionJointDef : public b2JointDef
{
	b2FrictionJointDef()
	{
		type = e_frictionJoint;
		localAnchorA.SetZero();
		localAnchorB.SetZero();
		maxForce = 0.0f;
		maxTorque = 0.0f;
	}

	/// Initialize the bodies, anchors, axis, and reference angle using the world
	/// anchor and world axis.
	void Initialize(b2Body* bodyA, b2Body* bodyB, const b2Vec2& anchor);

	/// The local anchor point relative to bodyA's origin.
	b2Vec2 localAnchorA;

	/// The local anchor point relative to bodyB's origin.
	b2Vec2 localAnchorB;

	/// The maximum friction force in N.
	float32 maxForce;

	/// The maximum friction torque in N-m.
	float32 maxTorque;
};

/// Friction joint. This is used for top-down friction.
/// It provides 2D translational friction and angular friction.
class b2FrictionJoint : public b2Joint
{
public:
	b2Vec2 GetAnchorA() const;
	b2Vec2 GetAnchorB() const;

	b2Vec2 GetReactionForce(float32 inv_dt) const;
	float32 GetReactionTorque(float32 inv_dt) const;

	/// The local anchor point relative to bodyA's origin.
	const b2Vec2& GetLocalAnchorA() const { return m_localAnchorA; }

	/// The local anchor point relative to bodyB's origin.
	const b2Vec2& GetLocalAnchorB() const  { return m_localAnchorB; }

	/// Set the maximum friction force in N.
	void SetMaxForce(float32 force);

	/// Get the maximum friction force in N.
	float32 GetMaxForce() const;

	/// Set the maximum friction torque in N*m.
	void SetMaxTorque(float32 torque);

	/// Get the maximum friction torque in N*m.
	float32 GetMaxTorque() const;

	/// Dump joint to dmLog
	void Dump();

protected:

	friend class b2Joint;

	b2FrictionJoint(const b2FrictionJointDef* def);

	void InitVelocityConstraints(const b2SolverData& data);
	void SolveVelocityConstraints(const b2SolverData& data);
	bool SolvePositionConstraints(const b2SolverData& data);

	b2Vec2 m_localAnchorA;
	b2Vec2 m_localAnchorB;

	// Solver shared
	b2Vec2 m_linearImpulse;
	float32 m_angularImpulse;
	float32 m_maxForce;
	float32 m_maxTorque;

	// Solver temp
	int32 m_indexA;
	int32 m_indexB;
	b2Vec2 m_rA;
	b2Vec2 m_rB;
	b2Vec2 m_localCenterA;
	b2Vec2 m_localCenterB;
	float32 m_invMassA;
	float32 m_invMassB;
	float32 m_invIA;
	float32 m_invIB;
	b2Mat22 m_linearMass;
	float32 m_angularMass;
};

// end of FrictionJoint.h

// Point-to-point constraint
// Cdot = v2 - v1
//      = v2 + cross(w2, r2) - v1 - cross(w1, r1)
// J = [-I -r1_skew I r2_skew ]
// Identity used:
// w k % (rx i + ry j) = w * (-ry i + rx j)

// Angle constraint
// Cdot = w2 - w1
// J = [0 0 -1 0 0 1]
// K = invI1 + invI2

void b2FrictionJointDef::Initialize(b2Body* bA, b2Body* bB, const b2Vec2& anchor)
{
	bodyA = bA;
	bodyB = bB;
	localAnchorA = bodyA->GetLocalPoint(anchor);
	localAnchorB = bodyB->GetLocalPoint(anchor);
}

b2FrictionJoint::b2FrictionJoint(const b2FrictionJointDef* def)
: b2Joint(def)
{
	m_localAnchorA = def->localAnchorA;
	m_localAnchorB = def->localAnchorB;

	m_linearImpulse.SetZero();
	m_angularImpulse = 0.0f;

	m_maxForce = def->maxForce;
	m_maxTorque = def->maxTorque;
}

void b2FrictionJoint::InitVelocityConstraints(const b2SolverData& data)
{
	m_indexA = m_bodyA->m_islandIndex;
	m_indexB = m_bodyB->m_islandIndex;
	m_localCenterA = m_bodyA->m_sweep.localCenter;
	m_localCenterB = m_bodyB->m_sweep.localCenter;
	m_invMassA = m_bodyA->m_invMass;
	m_invMassB = m_bodyB->m_invMass;
	m_invIA = m_bodyA->m_invI;
	m_invIB = m_bodyB->m_invI;

	float32 aA = data.positions[m_indexA].a;
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;

	float32 aB = data.positions[m_indexB].a;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	b2Rot qA(aA), qB(aB);

	// Compute the effective mass matrix.
	m_rA = b2Mul(qA, m_localAnchorA - m_localCenterA);
	m_rB = b2Mul(qB, m_localAnchorB - m_localCenterB);

	// J = [-I -r1_skew I r2_skew]
	//     [ 0       -1 0       1]
	// r_skew = [-ry; rx]

	// Matlab
	// K = [ mA+r1y^2*iA+mB+r2y^2*iB,  -r1y*iA*r1x-r2y*iB*r2x,          -r1y*iA-r2y*iB]
	//     [  -r1y*iA*r1x-r2y*iB*r2x, mA+r1x^2*iA+mB+r2x^2*iB,           r1x*iA+r2x*iB]
	//     [          -r1y*iA-r2y*iB,           r1x*iA+r2x*iB,                   iA+iB]

	float32 mA = m_invMassA, mB = m_invMassB;
	float32 iA = m_invIA, iB = m_invIB;

	b2Mat22 K;
	K.ex.x = mA + mB + iA * m_rA.y * m_rA.y + iB * m_rB.y * m_rB.y;
	K.ex.y = -iA * m_rA.x * m_rA.y - iB * m_rB.x * m_rB.y;
	K.ey.x = K.ex.y;
	K.ey.y = mA + mB + iA * m_rA.x * m_rA.x + iB * m_rB.x * m_rB.x;

	m_linearMass = K.GetInverse();

	m_angularMass = iA + iB;
	if (m_angularMass > 0.0f)
	{
		m_angularMass = 1.0f / m_angularMass;
	}

	if (data.step.warmStarting)
	{
		// Scale impulses to support a variable time step.
		m_linearImpulse *= data.step.dtRatio;
		m_angularImpulse *= data.step.dtRatio;

		b2Vec2 P(m_linearImpulse.x, m_linearImpulse.y);
		vA -= mA * P;
		wA -= iA * (b2Cross(m_rA, P) + m_angularImpulse);
		vB += mB * P;
		wB += iB * (b2Cross(m_rB, P) + m_angularImpulse);
	}
	else
	{
		m_linearImpulse.SetZero();
		m_angularImpulse = 0.0f;
	}

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

void b2FrictionJoint::SolveVelocityConstraints(const b2SolverData& data)
{
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	float32 mA = m_invMassA, mB = m_invMassB;
	float32 iA = m_invIA, iB = m_invIB;

	float32 h = data.step.dt;

	// Solve angular friction
	{
		float32 Cdot = wB - wA;
		float32 impulse = -m_angularMass * Cdot;

		float32 oldImpulse = m_angularImpulse;
		float32 maxImpulse = h * m_maxTorque;
		m_angularImpulse = b2Clamp(m_angularImpulse + impulse, -maxImpulse, maxImpulse);
		impulse = m_angularImpulse - oldImpulse;

		wA -= iA * impulse;
		wB += iB * impulse;
	}

	// Solve linear friction
	{
		b2Vec2 Cdot = vB + b2Cross(wB, m_rB) - vA - b2Cross(wA, m_rA);

		b2Vec2 impulse = -b2Mul(m_linearMass, Cdot);
		b2Vec2 oldImpulse = m_linearImpulse;
		m_linearImpulse += impulse;

		float32 maxImpulse = h * m_maxForce;

		if (m_linearImpulse.LengthSquared() > maxImpulse * maxImpulse)
		{
			m_linearImpulse.Normalize();
			m_linearImpulse *= maxImpulse;
		}

		impulse = m_linearImpulse - oldImpulse;

		vA -= mA * impulse;
		wA -= iA * b2Cross(m_rA, impulse);

		vB += mB * impulse;
		wB += iB * b2Cross(m_rB, impulse);
	}

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

bool b2FrictionJoint::SolvePositionConstraints(const b2SolverData& data)
{
	B2_NOT_USED(data);

	return true;
}

b2Vec2 b2FrictionJoint::GetAnchorA() const
{
	return m_bodyA->GetWorldPoint(m_localAnchorA);
}

b2Vec2 b2FrictionJoint::GetAnchorB() const
{
	return m_bodyB->GetWorldPoint(m_localAnchorB);
}

b2Vec2 b2FrictionJoint::GetReactionForce(float32 inv_dt) const
{
	return inv_dt * m_linearImpulse;
}

float32 b2FrictionJoint::GetReactionTorque(float32 inv_dt) const
{
	return inv_dt * m_angularImpulse;
}

void b2FrictionJoint::SetMaxForce(float32 force)
{
	b2Assert(b2IsValid(force) && force >= 0.0f);
	m_maxForce = force;
}

float32 b2FrictionJoint::GetMaxForce() const
{
	return m_maxForce;
}

void b2FrictionJoint::SetMaxTorque(float32 torque)
{
	b2Assert(b2IsValid(torque) && torque >= 0.0f);
	m_maxTorque = torque;
}

float32 b2FrictionJoint::GetMaxTorque() const
{
	return m_maxTorque;
}

void b2FrictionJoint::Dump()
{
	int32 indexA = m_bodyA->m_islandIndex;
	int32 indexB = m_bodyB->m_islandIndex;

	b2Log("  b2FrictionJointDef jd;\n");
	b2Log("  jd.bodyA = bodies[%d];\n", indexA);
	b2Log("  jd.bodyB = bodies[%d];\n", indexB);
	b2Log("  jd.collideConnected = bool(%d);\n", m_collideConnected);
	b2Log("  jd.localAnchorA.Set(%.15lef, %.15lef);\n", m_localAnchorA.x, m_localAnchorA.y);
	b2Log("  jd.localAnchorB.Set(%.15lef, %.15lef);\n", m_localAnchorB.x, m_localAnchorB.y);
	b2Log("  jd.maxForce = %.15lef;\n", m_maxForce);
	b2Log("  jd.maxTorque = %.15lef;\n", m_maxTorque);
	b2Log("  joints[%d] = m_world->CreateJoint(&jd);\n", m_index);
}

// end of FrictionJoint.cpp

/// Gear joint definition. This definition requires two existing
/// revolute or prismatic joints (any combination will work).
struct b2GearJointDef : public b2JointDef
{
	b2GearJointDef()
	{
		type = e_gearJoint;
		joint1 = NULL;
		joint2 = NULL;
		ratio = 1.0f;
	}

	/// The first revolute/prismatic joint attached to the gear joint.
	b2Joint* joint1;

	/// The second revolute/prismatic joint attached to the gear joint.
	b2Joint* joint2;

	/// The gear ratio.
	/// @see b2GearJoint for explanation.
	float32 ratio;
};

/// A gear joint is used to connect two joints together. Either joint
/// can be a revolute or prismatic joint. You specify a gear ratio
/// to bind the motions together:
/// coordinate1 + ratio * coordinate2 = constant
/// The ratio can be negative or positive. If one joint is a revolute joint
/// and the other joint is a prismatic joint, then the ratio will have units
/// of length or units of 1/length.
/// @warning You have to manually destroy the gear joint if joint1 or joint2
/// is destroyed.
class b2GearJoint : public b2Joint
{
public:
	b2Vec2 GetAnchorA() const;
	b2Vec2 GetAnchorB() const;

	b2Vec2 GetReactionForce(float32 inv_dt) const;
	float32 GetReactionTorque(float32 inv_dt) const;

	/// Get the first joint.
	b2Joint* GetJoint1() { return m_joint1; }

	/// Get the second joint.
	b2Joint* GetJoint2() { return m_joint2; }

	/// Set/Get the gear ratio.
	void SetRatio(float32 ratio);
	float32 GetRatio() const;

	/// Dump joint to dmLog
	void Dump();

protected:

	friend class b2Joint;
	b2GearJoint(const b2GearJointDef* data);

	void InitVelocityConstraints(const b2SolverData& data);
	void SolveVelocityConstraints(const b2SolverData& data);
	bool SolvePositionConstraints(const b2SolverData& data);

	b2Joint* m_joint1;
	b2Joint* m_joint2;

	b2JointType m_typeA;
	b2JointType m_typeB;

	// Body A is connected to body C
	// Body B is connected to body D
	b2Body* m_bodyC;
	b2Body* m_bodyD;

	// Solver shared
	b2Vec2 m_localAnchorA;
	b2Vec2 m_localAnchorB;
	b2Vec2 m_localAnchorC;
	b2Vec2 m_localAnchorD;

	b2Vec2 m_localAxisC;
	b2Vec2 m_localAxisD;

	float32 m_referenceAngleA;
	float32 m_referenceAngleB;

	float32 m_constant;
	float32 m_ratio;

	float32 m_impulse;

	// Solver temp
	int32 m_indexA, m_indexB, m_indexC, m_indexD;
	b2Vec2 m_lcA, m_lcB, m_lcC, m_lcD;
	float32 m_mA, m_mB, m_mC, m_mD;
	float32 m_iA, m_iB, m_iC, m_iD;
	b2Vec2 m_JvAC, m_JvBD;
	float32 m_JwA, m_JwB, m_JwC, m_JwD;
	float32 m_mass;
};

// end of GearJoint.h

/// Prismatic joint definition. This requires defining a line of
/// motion using an axis and an anchor point. The definition uses local
/// anchor points and a local axis so that the initial configuration
/// can violate the constraint slightly. The joint translation is zero
/// when the local anchor points coincide in world space. Using local
/// anchors and a local axis helps when saving and loading a game.
struct b2PrismaticJointDef : public b2JointDef
{
	b2PrismaticJointDef()
	{
		type = e_prismaticJoint;
		localAnchorA.SetZero();
		localAnchorB.SetZero();
		localAxisA.Set(1.0f, 0.0f);
		referenceAngle = 0.0f;
		enableLimit = false;
		lowerTranslation = 0.0f;
		upperTranslation = 0.0f;
		enableMotor = false;
		maxMotorForce = 0.0f;
		motorSpeed = 0.0f;
	}

	/// Initialize the bodies, anchors, axis, and reference angle using the world
	/// anchor and unit world axis.
	void Initialize(b2Body* bodyA, b2Body* bodyB, const b2Vec2& anchor, const b2Vec2& axis);

	/// The local anchor point relative to bodyA's origin.
	b2Vec2 localAnchorA;

	/// The local anchor point relative to bodyB's origin.
	b2Vec2 localAnchorB;

	/// The local translation unit axis in bodyA.
	b2Vec2 localAxisA;

	/// The constrained angle between the bodies: bodyB_angle - bodyA_angle.
	float32 referenceAngle;

	/// Enable/disable the joint limit.
	bool enableLimit;

	/// The lower translation limit, usually in meters.
	float32 lowerTranslation;

	/// The upper translation limit, usually in meters.
	float32 upperTranslation;

	/// Enable/disable the joint motor.
	bool enableMotor;

	/// The maximum motor torque, usually in N-m.
	float32 maxMotorForce;

	/// The desired motor speed in radians per second.
	float32 motorSpeed;
};

/// A prismatic joint. This joint provides one degree of freedom: translation
/// along an axis fixed in bodyA. Relative rotation is prevented. You can
/// use a joint limit to restrict the range of motion and a joint motor to
/// drive the motion or to model joint friction.
class b2PrismaticJoint : public b2Joint
{
public:
	b2Vec2 GetAnchorA() const;
	b2Vec2 GetAnchorB() const;

	b2Vec2 GetReactionForce(float32 inv_dt) const;
	float32 GetReactionTorque(float32 inv_dt) const;

	/// The local anchor point relative to bodyA's origin.
	const b2Vec2& GetLocalAnchorA() const { return m_localAnchorA; }

	/// The local anchor point relative to bodyB's origin.
	const b2Vec2& GetLocalAnchorB() const  { return m_localAnchorB; }

	/// The local joint axis relative to bodyA.
	const b2Vec2& GetLocalAxisA() const { return m_localXAxisA; }

	/// Get the reference angle.
	float32 GetReferenceAngle() const { return m_referenceAngle; }

	/// Get the current joint translation, usually in meters.
	float32 GetJointTranslation() const;

	/// Get the current joint translation speed, usually in meters per second.
	float32 GetJointSpeed() const;

	/// Is the joint limit enabled?
	bool IsLimitEnabled() const;

	/// Enable/disable the joint limit.
	void EnableLimit(bool flag);

	/// Get the lower joint limit, usually in meters.
	float32 GetLowerLimit() const;

	/// Get the upper joint limit, usually in meters.
	float32 GetUpperLimit() const;

	/// Set the joint limits, usually in meters.
	void SetLimits(float32 lower, float32 upper);

	/// Is the joint motor enabled?
	bool IsMotorEnabled() const;

	/// Enable/disable the joint motor.
	void EnableMotor(bool flag);

	/// Set the motor speed, usually in meters per second.
	void SetMotorSpeed(float32 speed);

	/// Get the motor speed, usually in meters per second.
	float32 GetMotorSpeed() const;

	/// Set the maximum motor force, usually in N.
	void SetMaxMotorForce(float32 force);
	float32 GetMaxMotorForce() const { return m_maxMotorForce; }

	/// Get the current motor force given the inverse time step, usually in N.
	float32 GetMotorForce(float32 inv_dt) const;

	/// Dump to b2Log
	void Dump();

protected:
	friend class b2Joint;
	friend class b2GearJoint;
	b2PrismaticJoint(const b2PrismaticJointDef* def);

	void InitVelocityConstraints(const b2SolverData& data);
	void SolveVelocityConstraints(const b2SolverData& data);
	bool SolvePositionConstraints(const b2SolverData& data);

	// Solver shared
	b2Vec2 m_localAnchorA;
	b2Vec2 m_localAnchorB;
	b2Vec2 m_localXAxisA;
	b2Vec2 m_localYAxisA;
	float32 m_referenceAngle;
	b2Vec3 m_impulse;
	float32 m_motorImpulse;
	float32 m_lowerTranslation;
	float32 m_upperTranslation;
	float32 m_maxMotorForce;
	float32 m_motorSpeed;
	bool m_enableLimit;
	bool m_enableMotor;
	b2LimitState m_limitState;

	// Solver temp
	int32 m_indexA;
	int32 m_indexB;
	b2Vec2 m_localCenterA;
	b2Vec2 m_localCenterB;
	float32 m_invMassA;
	float32 m_invMassB;
	float32 m_invIA;
	float32 m_invIB;
	b2Vec2 m_axis, m_perp;
	float32 m_s1, m_s2;
	float32 m_a1, m_a2;
	b2Mat33 m_K;
	float32 m_motorMass;
};

inline float32 b2PrismaticJoint::GetMotorSpeed() const
{
	return m_motorSpeed;
}

// end of PrismaticJoint.h

// Linear constraint (point-to-line)
// d = p2 - p1 = x2 + r2 - x1 - r1
// C = dot(perp, d)
// Cdot = dot(d, cross(w1, perp)) + dot(perp, v2 + cross(w2, r2) - v1 - cross(w1, r1))
//      = -dot(perp, v1) - dot(cross(d + r1, perp), w1) + dot(perp, v2) + dot(cross(r2, perp), v2)
// J = [-perp, -cross(d + r1, perp), perp, cross(r2,perp)]
//
// Angular constraint
// C = a2 - a1 + a_initial
// Cdot = w2 - w1
// J = [0 0 -1 0 0 1]
//
// K = J * invM * JT
//
// J = [-a -s1 a s2]
//     [0  -1  0  1]
// a = perp
// s1 = cross(d + r1, a) = cross(p2 - x1, a)
// s2 = cross(r2, a) = cross(p2 - x2, a)


// Motor/Limit linear constraint
// C = dot(ax1, d)
// Cdot = = -dot(ax1, v1) - dot(cross(d + r1, ax1), w1) + dot(ax1, v2) + dot(cross(r2, ax1), v2)
// J = [-ax1 -cross(d+r1,ax1) ax1 cross(r2,ax1)]

// Block Solver
// We develop a block solver that includes the joint limit. This makes the limit stiff (inelastic) even
// when the mass has poor distribution (leading to large torques about the joint anchor points).
//
// The Jacobian has 3 rows:
// J = [-uT -s1 uT s2] // linear
//     [0   -1   0  1] // angular
//     [-vT -a1 vT a2] // limit
//
// u = perp
// v = axis
// s1 = cross(d + r1, u), s2 = cross(r2, u)
// a1 = cross(d + r1, v), a2 = cross(r2, v)

// M * (v2 - v1) = JT * df
// J * v2 = bias
//
// v2 = v1 + invM * JT * df
// J * (v1 + invM * JT * df) = bias
// K * df = bias - J * v1 = -Cdot
// K = J * invM * JT
// Cdot = J * v1 - bias
//
// Now solve for f2.
// df = f2 - f1
// K * (f2 - f1) = -Cdot
// f2 = invK * (-Cdot) + f1
//
// Clamp accumulated limit impulse.
// lower: f2(3) = max(f2(3), 0)
// upper: f2(3) = min(f2(3), 0)
//
// Solve for correct f2(1:2)
// K(1:2, 1:2) * f2(1:2) = -Cdot(1:2) - K(1:2,3) * f2(3) + K(1:2,1:3) * f1
//                       = -Cdot(1:2) - K(1:2,3) * f2(3) + K(1:2,1:2) * f1(1:2) + K(1:2,3) * f1(3)
// K(1:2, 1:2) * f2(1:2) = -Cdot(1:2) - K(1:2,3) * (f2(3) - f1(3)) + K(1:2,1:2) * f1(1:2)
// f2(1:2) = invK(1:2,1:2) * (-Cdot(1:2) - K(1:2,3) * (f2(3) - f1(3))) + f1(1:2)
//
// Now compute impulse to be applied:
// df = f2 - f1

void b2PrismaticJointDef::Initialize(b2Body* bA, b2Body* bB, const b2Vec2& anchor, const b2Vec2& axis)
{
	bodyA = bA;
	bodyB = bB;
	localAnchorA = bodyA->GetLocalPoint(anchor);
	localAnchorB = bodyB->GetLocalPoint(anchor);
	localAxisA = bodyA->GetLocalVector(axis);
	referenceAngle = bodyB->GetAngle() - bodyA->GetAngle();
}

b2PrismaticJoint::b2PrismaticJoint(const b2PrismaticJointDef* def)
: b2Joint(def)
{
	m_localAnchorA = def->localAnchorA;
	m_localAnchorB = def->localAnchorB;
	m_localXAxisA = def->localAxisA;
	m_localXAxisA.Normalize();
	m_localYAxisA = b2Cross(1.0f, m_localXAxisA);
	m_referenceAngle = def->referenceAngle;

	m_impulse.SetZero();
	m_motorMass = 0.0f;
	m_motorImpulse = 0.0f;

	m_lowerTranslation = def->lowerTranslation;
	m_upperTranslation = def->upperTranslation;
	m_maxMotorForce = def->maxMotorForce;
	m_motorSpeed = def->motorSpeed;
	m_enableLimit = def->enableLimit;
	m_enableMotor = def->enableMotor;
	m_limitState = e_inactiveLimit;

	m_axis.SetZero();
	m_perp.SetZero();
}

void b2PrismaticJoint::InitVelocityConstraints(const b2SolverData& data)
{
	m_indexA = m_bodyA->m_islandIndex;
	m_indexB = m_bodyB->m_islandIndex;
	m_localCenterA = m_bodyA->m_sweep.localCenter;
	m_localCenterB = m_bodyB->m_sweep.localCenter;
	m_invMassA = m_bodyA->m_invMass;
	m_invMassB = m_bodyB->m_invMass;
	m_invIA = m_bodyA->m_invI;
	m_invIB = m_bodyB->m_invI;

	b2Vec2 cA = data.positions[m_indexA].c;
	float32 aA = data.positions[m_indexA].a;
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;

	b2Vec2 cB = data.positions[m_indexB].c;
	float32 aB = data.positions[m_indexB].a;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	b2Rot qA(aA), qB(aB);

	// Compute the effective masses.
	b2Vec2 rA = b2Mul(qA, m_localAnchorA - m_localCenterA);
	b2Vec2 rB = b2Mul(qB, m_localAnchorB - m_localCenterB);
	b2Vec2 d = (cB - cA) + rB - rA;

	float32 mA = m_invMassA, mB = m_invMassB;
	float32 iA = m_invIA, iB = m_invIB;

	// Compute motor Jacobian and effective mass.
	{
		m_axis = b2Mul(qA, m_localXAxisA);
		m_a1 = b2Cross(d + rA, m_axis);
		m_a2 = b2Cross(rB, m_axis);

		m_motorMass = mA + mB + iA * m_a1 * m_a1 + iB * m_a2 * m_a2;
		if (m_motorMass > 0.0f)
		{
			m_motorMass = 1.0f / m_motorMass;
		}
	}

	// Prismatic constraint.
	{
		m_perp = b2Mul(qA, m_localYAxisA);

		m_s1 = b2Cross(d + rA, m_perp);
		m_s2 = b2Cross(rB, m_perp);

		float32 k11 = mA + mB + iA * m_s1 * m_s1 + iB * m_s2 * m_s2;
		float32 k12 = iA * m_s1 + iB * m_s2;
		float32 k13 = iA * m_s1 * m_a1 + iB * m_s2 * m_a2;
		float32 k22 = iA + iB;
		if (k22 == 0.0f)
		{
			// For bodies with fixed rotation.
			k22 = 1.0f;
		}
		float32 k23 = iA * m_a1 + iB * m_a2;
		float32 k33 = mA + mB + iA * m_a1 * m_a1 + iB * m_a2 * m_a2;

		m_K.ex.Set(k11, k12, k13);
		m_K.ey.Set(k12, k22, k23);
		m_K.ez.Set(k13, k23, k33);
	}

	// Compute motor and limit terms.
	if (m_enableLimit)
	{
		float32 jointTranslation = b2Dot(m_axis, d);
		if (b2Abs(m_upperTranslation - m_lowerTranslation) < 2.0f * b2_linearSlop)
		{
			m_limitState = e_equalLimits;
		}
		else if (jointTranslation <= m_lowerTranslation)
		{
			if (m_limitState != e_atLowerLimit)
			{
				m_limitState = e_atLowerLimit;
				m_impulse.z = 0.0f;
			}
		}
		else if (jointTranslation >= m_upperTranslation)
		{
			if (m_limitState != e_atUpperLimit)
			{
				m_limitState = e_atUpperLimit;
				m_impulse.z = 0.0f;
			}
		}
		else
		{
			m_limitState = e_inactiveLimit;
			m_impulse.z = 0.0f;
		}
	}
	else
	{
		m_limitState = e_inactiveLimit;
		m_impulse.z = 0.0f;
	}

	if (m_enableMotor == false)
	{
		m_motorImpulse = 0.0f;
	}

	if (data.step.warmStarting)
	{
		// Account for variable time step.
		m_impulse *= data.step.dtRatio;
		m_motorImpulse *= data.step.dtRatio;

		b2Vec2 P = m_impulse.x * m_perp + (m_motorImpulse + m_impulse.z) * m_axis;
		float32 LA = m_impulse.x * m_s1 + m_impulse.y + (m_motorImpulse + m_impulse.z) * m_a1;
		float32 LB = m_impulse.x * m_s2 + m_impulse.y + (m_motorImpulse + m_impulse.z) * m_a2;

		vA -= mA * P;
		wA -= iA * LA;

		vB += mB * P;
		wB += iB * LB;
	}
	else
	{
		m_impulse.SetZero();
		m_motorImpulse = 0.0f;
	}

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

void b2PrismaticJoint::SolveVelocityConstraints(const b2SolverData& data)
{
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	float32 mA = m_invMassA, mB = m_invMassB;
	float32 iA = m_invIA, iB = m_invIB;

	// Solve linear motor constraint.
	if (m_enableMotor && m_limitState != e_equalLimits)
	{
		float32 Cdot = b2Dot(m_axis, vB - vA) + m_a2 * wB - m_a1 * wA;
		float32 impulse = m_motorMass * (m_motorSpeed - Cdot);
		float32 oldImpulse = m_motorImpulse;
		float32 maxImpulse = data.step.dt * m_maxMotorForce;
		m_motorImpulse = b2Clamp(m_motorImpulse + impulse, -maxImpulse, maxImpulse);
		impulse = m_motorImpulse - oldImpulse;

		b2Vec2 P = impulse * m_axis;
		float32 LA = impulse * m_a1;
		float32 LB = impulse * m_a2;

		vA -= mA * P;
		wA -= iA * LA;

		vB += mB * P;
		wB += iB * LB;
	}

	b2Vec2 Cdot1;
	Cdot1.x = b2Dot(m_perp, vB - vA) + m_s2 * wB - m_s1 * wA;
	Cdot1.y = wB - wA;

	if (m_enableLimit && m_limitState != e_inactiveLimit)
	{
		// Solve prismatic and limit constraint in block form.
		float32 Cdot2;
		Cdot2 = b2Dot(m_axis, vB - vA) + m_a2 * wB - m_a1 * wA;
		b2Vec3 Cdot(Cdot1.x, Cdot1.y, Cdot2);

		b2Vec3 f1 = m_impulse;
		b2Vec3 df =  m_K.Solve33(-Cdot);
		m_impulse += df;

		if (m_limitState == e_atLowerLimit)
		{
			m_impulse.z = b2Max(m_impulse.z, 0.0f);
		}
		else if (m_limitState == e_atUpperLimit)
		{
			m_impulse.z = b2Min(m_impulse.z, 0.0f);
		}

		// f2(1:2) = invK(1:2,1:2) * (-Cdot(1:2) - K(1:2,3) * (f2(3) - f1(3))) + f1(1:2)
		b2Vec2 b = -Cdot1 - (m_impulse.z - f1.z) * b2Vec2(m_K.ez.x, m_K.ez.y);
		b2Vec2 f2r = m_K.Solve22(b) + b2Vec2(f1.x, f1.y);
		m_impulse.x = f2r.x;
		m_impulse.y = f2r.y;

		df = m_impulse - f1;

		b2Vec2 P = df.x * m_perp + df.z * m_axis;
		float32 LA = df.x * m_s1 + df.y + df.z * m_a1;
		float32 LB = df.x * m_s2 + df.y + df.z * m_a2;

		vA -= mA * P;
		wA -= iA * LA;

		vB += mB * P;
		wB += iB * LB;
	}
	else
	{
		// Limit is inactive, just solve the prismatic constraint in block form.
		b2Vec2 df = m_K.Solve22(-Cdot1);
		m_impulse.x += df.x;
		m_impulse.y += df.y;

		b2Vec2 P = df.x * m_perp;
		float32 LA = df.x * m_s1 + df.y;
		float32 LB = df.x * m_s2 + df.y;

		vA -= mA * P;
		wA -= iA * LA;

		vB += mB * P;
		wB += iB * LB;
	}

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

bool b2PrismaticJoint::SolvePositionConstraints(const b2SolverData& data)
{
	b2Vec2 cA = data.positions[m_indexA].c;
	float32 aA = data.positions[m_indexA].a;
	b2Vec2 cB = data.positions[m_indexB].c;
	float32 aB = data.positions[m_indexB].a;

	b2Rot qA(aA), qB(aB);

	float32 mA = m_invMassA, mB = m_invMassB;
	float32 iA = m_invIA, iB = m_invIB;

	// Compute fresh Jacobians
	b2Vec2 rA = b2Mul(qA, m_localAnchorA - m_localCenterA);
	b2Vec2 rB = b2Mul(qB, m_localAnchorB - m_localCenterB);
	b2Vec2 d = cB + rB - cA - rA;

	b2Vec2 axis = b2Mul(qA, m_localXAxisA);
	float32 a1 = b2Cross(d + rA, axis);
	float32 a2 = b2Cross(rB, axis);
	b2Vec2 perp = b2Mul(qA, m_localYAxisA);

	float32 s1 = b2Cross(d + rA, perp);
	float32 s2 = b2Cross(rB, perp);

	b2Vec3 impulse;
	b2Vec2 C1;
	C1.x = b2Dot(perp, d);
	C1.y = aB - aA - m_referenceAngle;

	float32 linearError = b2Abs(C1.x);
	float32 angularError = b2Abs(C1.y);

	bool active = false;
	float32 C2 = 0.0f;
	if (m_enableLimit)
	{
		float32 translation = b2Dot(axis, d);
		if (b2Abs(m_upperTranslation - m_lowerTranslation) < 2.0f * b2_linearSlop)
		{
			// Prevent large angular corrections
			C2 = b2Clamp(translation, -b2_maxLinearCorrection, b2_maxLinearCorrection);
			linearError = b2Max(linearError, b2Abs(translation));
			active = true;
		}
		else if (translation <= m_lowerTranslation)
		{
			// Prevent large linear corrections and allow some slop.
			C2 = b2Clamp(translation - m_lowerTranslation + b2_linearSlop, -b2_maxLinearCorrection, 0.0f);
			linearError = b2Max(linearError, m_lowerTranslation - translation);
			active = true;
		}
		else if (translation >= m_upperTranslation)
		{
			// Prevent large linear corrections and allow some slop.
			C2 = b2Clamp(translation - m_upperTranslation - b2_linearSlop, 0.0f, b2_maxLinearCorrection);
			linearError = b2Max(linearError, translation - m_upperTranslation);
			active = true;
		}
	}

	if (active)
	{
		float32 k11 = mA + mB + iA * s1 * s1 + iB * s2 * s2;
		float32 k12 = iA * s1 + iB * s2;
		float32 k13 = iA * s1 * a1 + iB * s2 * a2;
		float32 k22 = iA + iB;
		if (k22 == 0.0f)
		{
			// For fixed rotation
			k22 = 1.0f;
		}
		float32 k23 = iA * a1 + iB * a2;
		float32 k33 = mA + mB + iA * a1 * a1 + iB * a2 * a2;

		b2Mat33 K;
		K.ex.Set(k11, k12, k13);
		K.ey.Set(k12, k22, k23);
		K.ez.Set(k13, k23, k33);

		b2Vec3 C;
		C.x = C1.x;
		C.y = C1.y;
		C.z = C2;

		impulse = K.Solve33(-C);
	}
	else
	{
		float32 k11 = mA + mB + iA * s1 * s1 + iB * s2 * s2;
		float32 k12 = iA * s1 + iB * s2;
		float32 k22 = iA + iB;
		if (k22 == 0.0f)
		{
			k22 = 1.0f;
		}

		b2Mat22 K;
		K.ex.Set(k11, k12);
		K.ey.Set(k12, k22);

		b2Vec2 impulse1 = K.Solve(-C1);
		impulse.x = impulse1.x;
		impulse.y = impulse1.y;
		impulse.z = 0.0f;
	}

	b2Vec2 P = impulse.x * perp + impulse.z * axis;
	float32 LA = impulse.x * s1 + impulse.y + impulse.z * a1;
	float32 LB = impulse.x * s2 + impulse.y + impulse.z * a2;

	cA -= mA * P;
	aA -= iA * LA;
	cB += mB * P;
	aB += iB * LB;

	data.positions[m_indexA].c = cA;
	data.positions[m_indexA].a = aA;
	data.positions[m_indexB].c = cB;
	data.positions[m_indexB].a = aB;

	return linearError <= b2_linearSlop && angularError <= b2_angularSlop;
}

b2Vec2 b2PrismaticJoint::GetAnchorA() const
{
	return m_bodyA->GetWorldPoint(m_localAnchorA);
}

b2Vec2 b2PrismaticJoint::GetAnchorB() const
{
	return m_bodyB->GetWorldPoint(m_localAnchorB);
}

b2Vec2 b2PrismaticJoint::GetReactionForce(float32 inv_dt) const
{
	return inv_dt * (m_impulse.x * m_perp + (m_motorImpulse + m_impulse.z) * m_axis);
}

float32 b2PrismaticJoint::GetReactionTorque(float32 inv_dt) const
{
	return inv_dt * m_impulse.y;
}

float32 b2PrismaticJoint::GetJointTranslation() const
{
	b2Vec2 pA = m_bodyA->GetWorldPoint(m_localAnchorA);
	b2Vec2 pB = m_bodyB->GetWorldPoint(m_localAnchorB);
	b2Vec2 d = pB - pA;
	b2Vec2 axis = m_bodyA->GetWorldVector(m_localXAxisA);

	float32 translation = b2Dot(d, axis);
	return translation;
}

float32 b2PrismaticJoint::GetJointSpeed() const
{
	b2Body* bA = m_bodyA;
	b2Body* bB = m_bodyB;

	b2Vec2 rA = b2Mul(bA->m_xf.q, m_localAnchorA - bA->m_sweep.localCenter);
	b2Vec2 rB = b2Mul(bB->m_xf.q, m_localAnchorB - bB->m_sweep.localCenter);
	b2Vec2 p1 = bA->m_sweep.c + rA;
	b2Vec2 p2 = bB->m_sweep.c + rB;
	b2Vec2 d = p2 - p1;
	b2Vec2 axis = b2Mul(bA->m_xf.q, m_localXAxisA);

	b2Vec2 vA = bA->m_linearVelocity;
	b2Vec2 vB = bB->m_linearVelocity;
	float32 wA = bA->m_angularVelocity;
	float32 wB = bB->m_angularVelocity;

	float32 speed = b2Dot(d, b2Cross(wA, axis)) + b2Dot(axis, vB + b2Cross(wB, rB) - vA - b2Cross(wA, rA));
	return speed;
}

bool b2PrismaticJoint::IsLimitEnabled() const
{
	return m_enableLimit;
}

void b2PrismaticJoint::EnableLimit(bool flag)
{
	if (flag != m_enableLimit)
	{
		m_bodyA->SetAwake(true);
		m_bodyB->SetAwake(true);
		m_enableLimit = flag;
		m_impulse.z = 0.0f;
	}
}

float32 b2PrismaticJoint::GetLowerLimit() const
{
	return m_lowerTranslation;
}

float32 b2PrismaticJoint::GetUpperLimit() const
{
	return m_upperTranslation;
}

void b2PrismaticJoint::SetLimits(float32 lower, float32 upper)
{
	b2Assert(lower <= upper);
	if (lower != m_lowerTranslation || upper != m_upperTranslation)
	{
		m_bodyA->SetAwake(true);
		m_bodyB->SetAwake(true);
		m_lowerTranslation = lower;
		m_upperTranslation = upper;
		m_impulse.z = 0.0f;
	}
}

bool b2PrismaticJoint::IsMotorEnabled() const
{
	return m_enableMotor;
}

void b2PrismaticJoint::EnableMotor(bool flag)
{
	m_bodyA->SetAwake(true);
	m_bodyB->SetAwake(true);
	m_enableMotor = flag;
}

void b2PrismaticJoint::SetMotorSpeed(float32 speed)
{
	m_bodyA->SetAwake(true);
	m_bodyB->SetAwake(true);
	m_motorSpeed = speed;
}

void b2PrismaticJoint::SetMaxMotorForce(float32 force)
{
	m_bodyA->SetAwake(true);
	m_bodyB->SetAwake(true);
	m_maxMotorForce = force;
}

float32 b2PrismaticJoint::GetMotorForce(float32 inv_dt) const
{
	return inv_dt * m_motorImpulse;
}

void b2PrismaticJoint::Dump()
{
	int32 indexA = m_bodyA->m_islandIndex;
	int32 indexB = m_bodyB->m_islandIndex;

	b2Log("  b2PrismaticJointDef jd;\n");
	b2Log("  jd.bodyA = bodies[%d];\n", indexA);
	b2Log("  jd.bodyB = bodies[%d];\n", indexB);
	b2Log("  jd.collideConnected = bool(%d);\n", m_collideConnected);
	b2Log("  jd.localAnchorA.Set(%.15lef, %.15lef);\n", m_localAnchorA.x, m_localAnchorA.y);
	b2Log("  jd.localAnchorB.Set(%.15lef, %.15lef);\n", m_localAnchorB.x, m_localAnchorB.y);
	b2Log("  jd.localAxisA.Set(%.15lef, %.15lef);\n", m_localXAxisA.x, m_localXAxisA.y);
	b2Log("  jd.referenceAngle = %.15lef;\n", m_referenceAngle);
	b2Log("  jd.enableLimit = bool(%d);\n", m_enableLimit);
	b2Log("  jd.lowerTranslation = %.15lef;\n", m_lowerTranslation);
	b2Log("  jd.upperTranslation = %.15lef;\n", m_upperTranslation);
	b2Log("  jd.enableMotor = bool(%d);\n", m_enableMotor);
	b2Log("  jd.motorSpeed = %.15lef;\n", m_motorSpeed);
	b2Log("  jd.maxMotorForce = %.15lef;\n", m_maxMotorForce);
	b2Log("  joints[%d] = m_world->CreateJoint(&jd);\n", m_index);
}

// end of PrismaticJoint.cpp

/// Revolute joint definition. This requires defining an
/// anchor point where the bodies are joined. The definition
/// uses local anchor points so that the initial configuration
/// can violate the constraint slightly. You also need to
/// specify the initial relative angle for joint limits. This
/// helps when saving and loading a game.
/// The local anchor points are measured from the body's origin
/// rather than the center of mass because:
/// 1. you might not know where the center of mass will be.
/// 2. if you add/remove shapes from a body and recompute the mass,
///    the joints will be broken.
struct b2RevoluteJointDef : public b2JointDef
{
	b2RevoluteJointDef()
	{
		type = e_revoluteJoint;
		localAnchorA.Set(0.0f, 0.0f);
		localAnchorB.Set(0.0f, 0.0f);
		referenceAngle = 0.0f;
		lowerAngle = 0.0f;
		upperAngle = 0.0f;
		maxMotorTorque = 0.0f;
		motorSpeed = 0.0f;
		enableLimit = false;
		enableMotor = false;
	}

	/// Initialize the bodies, anchors, and reference angle using a world
	/// anchor point.
	void Initialize(b2Body* bodyA, b2Body* bodyB, const b2Vec2& anchor);

	/// The local anchor point relative to bodyA's origin.
	b2Vec2 localAnchorA;

	/// The local anchor point relative to bodyB's origin.
	b2Vec2 localAnchorB;

	/// The bodyB angle minus bodyA angle in the reference state (radians).
	float32 referenceAngle;

	/// A flag to enable joint limits.
	bool enableLimit;

	/// The lower angle for the joint limit (radians).
	float32 lowerAngle;

	/// The upper angle for the joint limit (radians).
	float32 upperAngle;

	/// A flag to enable the joint motor.
	bool enableMotor;

	/// The desired motor speed. Usually in radians per second.
	float32 motorSpeed;

	/// The maximum motor torque used to achieve the desired motor speed.
	/// Usually in N-m.
	float32 maxMotorTorque;
};

/// A revolute joint constrains two bodies to share a common point while they
/// are free to rotate about the point. The relative rotation about the shared
/// point is the joint angle. You can limit the relative rotation with
/// a joint limit that specifies a lower and upper angle. You can use a motor
/// to drive the relative rotation about the shared point. A maximum motor torque
/// is provided so that infinite forces are not generated.
class b2RevoluteJoint : public b2Joint
{
public:
	b2Vec2 GetAnchorA() const;
	b2Vec2 GetAnchorB() const;

	/// The local anchor point relative to bodyA's origin.
	const b2Vec2& GetLocalAnchorA() const { return m_localAnchorA; }

	/// The local anchor point relative to bodyB's origin.
	const b2Vec2& GetLocalAnchorB() const  { return m_localAnchorB; }

	/// Get the reference angle.
	float32 GetReferenceAngle() const { return m_referenceAngle; }

	/// Get the current joint angle in radians.
	float32 GetJointAngle() const;

	/// Get the current joint angle speed in radians per second.
	float32 GetJointSpeed() const;

	/// Is the joint limit enabled?
	bool IsLimitEnabled() const;

	/// Enable/disable the joint limit.
	void EnableLimit(bool flag);

	/// Get the lower joint limit in radians.
	float32 GetLowerLimit() const;

	/// Get the upper joint limit in radians.
	float32 GetUpperLimit() const;

	/// Set the joint limits in radians.
	void SetLimits(float32 lower, float32 upper);

	/// Is the joint motor enabled?
	bool IsMotorEnabled() const;

	/// Enable/disable the joint motor.
	void EnableMotor(bool flag);

	/// Set the motor speed in radians per second.
	void SetMotorSpeed(float32 speed);

	/// Get the motor speed in radians per second.
	float32 GetMotorSpeed() const;

	/// Set the maximum motor torque, usually in N-m.
	void SetMaxMotorTorque(float32 torque);
	float32 GetMaxMotorTorque() const { return m_maxMotorTorque; }

	/// Get the reaction force given the inverse time step.
	/// Unit is N.
	b2Vec2 GetReactionForce(float32 inv_dt) const;

	/// Get the reaction torque due to the joint limit given the inverse time step.
	/// Unit is N*m.
	float32 GetReactionTorque(float32 inv_dt) const;

	/// Get the current motor torque given the inverse time step.
	/// Unit is N*m.
	float32 GetMotorTorque(float32 inv_dt) const;

	/// Dump to b2Log.
	void Dump();

protected:
	
	friend class b2Joint;
	friend class b2GearJoint;

	b2RevoluteJoint(const b2RevoluteJointDef* def);

	void InitVelocityConstraints(const b2SolverData& data);
	void SolveVelocityConstraints(const b2SolverData& data);
	bool SolvePositionConstraints(const b2SolverData& data);

	// Solver shared
	b2Vec2 m_localAnchorA;
	b2Vec2 m_localAnchorB;
	b2Vec3 m_impulse;
	float32 m_motorImpulse;

	bool m_enableMotor;
	float32 m_maxMotorTorque;
	float32 m_motorSpeed;

	bool m_enableLimit;
	float32 m_referenceAngle;
	float32 m_lowerAngle;
	float32 m_upperAngle;

	// Solver temp
	int32 m_indexA;
	int32 m_indexB;
	b2Vec2 m_rA;
	b2Vec2 m_rB;
	b2Vec2 m_localCenterA;
	b2Vec2 m_localCenterB;
	float32 m_invMassA;
	float32 m_invMassB;
	float32 m_invIA;
	float32 m_invIB;
	b2Mat33 m_mass;			// effective mass for point-to-point constraint.
	float32 m_motorMass;	// effective mass for motor/limit angular constraint.
	b2LimitState m_limitState;
};

inline float32 b2RevoluteJoint::GetMotorSpeed() const
{
	return m_motorSpeed;
}

// end of RevoluteJoint.h

// Point-to-point constraint
// C = p2 - p1
// Cdot = v2 - v1
//      = v2 + cross(w2, r2) - v1 - cross(w1, r1)
// J = [-I -r1_skew I r2_skew ]
// Identity used:
// w k % (rx i + ry j) = w * (-ry i + rx j)

// Motor constraint
// Cdot = w2 - w1
// J = [0 0 -1 0 0 1]
// K = invI1 + invI2

void b2RevoluteJointDef::Initialize(b2Body* bA, b2Body* bB, const b2Vec2& anchor)
{
	bodyA = bA;
	bodyB = bB;
	localAnchorA = bodyA->GetLocalPoint(anchor);
	localAnchorB = bodyB->GetLocalPoint(anchor);
	referenceAngle = bodyB->GetAngle() - bodyA->GetAngle();
}

b2RevoluteJoint::b2RevoluteJoint(const b2RevoluteJointDef* def)
: b2Joint(def)
{
	m_localAnchorA = def->localAnchorA;
	m_localAnchorB = def->localAnchorB;
	m_referenceAngle = def->referenceAngle;

	m_impulse.SetZero();
	m_motorImpulse = 0.0f;

	m_lowerAngle = def->lowerAngle;
	m_upperAngle = def->upperAngle;
	m_maxMotorTorque = def->maxMotorTorque;
	m_motorSpeed = def->motorSpeed;
	m_enableLimit = def->enableLimit;
	m_enableMotor = def->enableMotor;
	m_limitState = e_inactiveLimit;
}

void b2RevoluteJoint::InitVelocityConstraints(const b2SolverData& data)
{
	m_indexA = m_bodyA->m_islandIndex;
	m_indexB = m_bodyB->m_islandIndex;
	m_localCenterA = m_bodyA->m_sweep.localCenter;
	m_localCenterB = m_bodyB->m_sweep.localCenter;
	m_invMassA = m_bodyA->m_invMass;
	m_invMassB = m_bodyB->m_invMass;
	m_invIA = m_bodyA->m_invI;
	m_invIB = m_bodyB->m_invI;

	float32 aA = data.positions[m_indexA].a;
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;

	float32 aB = data.positions[m_indexB].a;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	b2Rot qA(aA), qB(aB);

	m_rA = b2Mul(qA, m_localAnchorA - m_localCenterA);
	m_rB = b2Mul(qB, m_localAnchorB - m_localCenterB);

	// J = [-I -r1_skew I r2_skew]
	//     [ 0       -1 0       1]
	// r_skew = [-ry; rx]

	// Matlab
	// K = [ mA+r1y^2*iA+mB+r2y^2*iB,  -r1y*iA*r1x-r2y*iB*r2x,          -r1y*iA-r2y*iB]
	//     [  -r1y*iA*r1x-r2y*iB*r2x, mA+r1x^2*iA+mB+r2x^2*iB,           r1x*iA+r2x*iB]
	//     [          -r1y*iA-r2y*iB,           r1x*iA+r2x*iB,                   iA+iB]

	float32 mA = m_invMassA, mB = m_invMassB;
	float32 iA = m_invIA, iB = m_invIB;

	bool fixedRotation = (iA + iB == 0.0f);

	m_mass.ex.x = mA + mB + m_rA.y * m_rA.y * iA + m_rB.y * m_rB.y * iB;
	m_mass.ey.x = -m_rA.y * m_rA.x * iA - m_rB.y * m_rB.x * iB;
	m_mass.ez.x = -m_rA.y * iA - m_rB.y * iB;
	m_mass.ex.y = m_mass.ey.x;
	m_mass.ey.y = mA + mB + m_rA.x * m_rA.x * iA + m_rB.x * m_rB.x * iB;
	m_mass.ez.y = m_rA.x * iA + m_rB.x * iB;
	m_mass.ex.z = m_mass.ez.x;
	m_mass.ey.z = m_mass.ez.y;
	m_mass.ez.z = iA + iB;

	m_motorMass = iA + iB;
	if (m_motorMass > 0.0f)
	{
		m_motorMass = 1.0f / m_motorMass;
	}

	if (m_enableMotor == false || fixedRotation)
	{
		m_motorImpulse = 0.0f;
	}

	if (m_enableLimit && fixedRotation == false)
	{
		float32 jointAngle = aB - aA - m_referenceAngle;
		if (b2Abs(m_upperAngle - m_lowerAngle) < 2.0f * b2_angularSlop)
		{
			m_limitState = e_equalLimits;
		}
		else if (jointAngle <= m_lowerAngle)
		{
			if (m_limitState != e_atLowerLimit)
			{
				m_impulse.z = 0.0f;
			}
			m_limitState = e_atLowerLimit;
		}
		else if (jointAngle >= m_upperAngle)
		{
			if (m_limitState != e_atUpperLimit)
			{
				m_impulse.z = 0.0f;
			}
			m_limitState = e_atUpperLimit;
		}
		else
		{
			m_limitState = e_inactiveLimit;
			m_impulse.z = 0.0f;
		}
	}
	else
	{
		m_limitState = e_inactiveLimit;
	}

	if (data.step.warmStarting)
	{
		// Scale impulses to support a variable time step.
		m_impulse *= data.step.dtRatio;
		m_motorImpulse *= data.step.dtRatio;

		b2Vec2 P(m_impulse.x, m_impulse.y);

		vA -= mA * P;
		wA -= iA * (b2Cross(m_rA, P) + m_motorImpulse + m_impulse.z);

		vB += mB * P;
		wB += iB * (b2Cross(m_rB, P) + m_motorImpulse + m_impulse.z);
	}
	else
	{
		m_impulse.SetZero();
		m_motorImpulse = 0.0f;
	}

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

void b2RevoluteJoint::SolveVelocityConstraints(const b2SolverData& data)
{
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	float32 mA = m_invMassA, mB = m_invMassB;
	float32 iA = m_invIA, iB = m_invIB;

	bool fixedRotation = (iA + iB == 0.0f);

	// Solve motor constraint.
	if (m_enableMotor && m_limitState != e_equalLimits && fixedRotation == false)
	{
		float32 Cdot = wB - wA - m_motorSpeed;
		float32 impulse = -m_motorMass * Cdot;
		float32 oldImpulse = m_motorImpulse;
		float32 maxImpulse = data.step.dt * m_maxMotorTorque;
		m_motorImpulse = b2Clamp(m_motorImpulse + impulse, -maxImpulse, maxImpulse);
		impulse = m_motorImpulse - oldImpulse;

		wA -= iA * impulse;
		wB += iB * impulse;
	}

	// Solve limit constraint.
	if (m_enableLimit && m_limitState != e_inactiveLimit && fixedRotation == false)
	{
		b2Vec2 Cdot1 = vB + b2Cross(wB, m_rB) - vA - b2Cross(wA, m_rA);
		float32 Cdot2 = wB - wA;
		b2Vec3 Cdot(Cdot1.x, Cdot1.y, Cdot2);

		b2Vec3 impulse = -m_mass.Solve33(Cdot);

		if (m_limitState == e_equalLimits)
		{
			m_impulse += impulse;
		}
		else if (m_limitState == e_atLowerLimit)
		{
			float32 newImpulse = m_impulse.z + impulse.z;
			if (newImpulse < 0.0f)
			{
				b2Vec2 rhs = -Cdot1 + m_impulse.z * b2Vec2(m_mass.ez.x, m_mass.ez.y);
				b2Vec2 reduced = m_mass.Solve22(rhs);
				impulse.x = reduced.x;
				impulse.y = reduced.y;
				impulse.z = -m_impulse.z;
				m_impulse.x += reduced.x;
				m_impulse.y += reduced.y;
				m_impulse.z = 0.0f;
			}
			else
			{
				m_impulse += impulse;
			}
		}
		else if (m_limitState == e_atUpperLimit)
		{
			float32 newImpulse = m_impulse.z + impulse.z;
			if (newImpulse > 0.0f)
			{
				b2Vec2 rhs = -Cdot1 + m_impulse.z * b2Vec2(m_mass.ez.x, m_mass.ez.y);
				b2Vec2 reduced = m_mass.Solve22(rhs);
				impulse.x = reduced.x;
				impulse.y = reduced.y;
				impulse.z = -m_impulse.z;
				m_impulse.x += reduced.x;
				m_impulse.y += reduced.y;
				m_impulse.z = 0.0f;
			}
			else
			{
				m_impulse += impulse;
			}
		}

		b2Vec2 P(impulse.x, impulse.y);

		vA -= mA * P;
		wA -= iA * (b2Cross(m_rA, P) + impulse.z);

		vB += mB * P;
		wB += iB * (b2Cross(m_rB, P) + impulse.z);
	}
	else
	{
		// Solve point-to-point constraint
		b2Vec2 Cdot = vB + b2Cross(wB, m_rB) - vA - b2Cross(wA, m_rA);
		b2Vec2 impulse = m_mass.Solve22(-Cdot);

		m_impulse.x += impulse.x;
		m_impulse.y += impulse.y;

		vA -= mA * impulse;
		wA -= iA * b2Cross(m_rA, impulse);

		vB += mB * impulse;
		wB += iB * b2Cross(m_rB, impulse);
	}

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

bool b2RevoluteJoint::SolvePositionConstraints(const b2SolverData& data)
{
	b2Vec2 cA = data.positions[m_indexA].c;
	float32 aA = data.positions[m_indexA].a;
	b2Vec2 cB = data.positions[m_indexB].c;
	float32 aB = data.positions[m_indexB].a;

	b2Rot qA(aA), qB(aB);

	float32 angularError = 0.0f;
	float32 positionError = 0.0f;

	bool fixedRotation = (m_invIA + m_invIB == 0.0f);

	// Solve angular limit constraint.
	if (m_enableLimit && m_limitState != e_inactiveLimit && fixedRotation == false)
	{
		float32 angle = aB - aA - m_referenceAngle;
		float32 limitImpulse = 0.0f;

		if (m_limitState == e_equalLimits)
		{
			// Prevent large angular corrections
			float32 C = b2Clamp(angle - m_lowerAngle, -b2_maxAngularCorrection, b2_maxAngularCorrection);
			limitImpulse = -m_motorMass * C;
			angularError = b2Abs(C);
		}
		else if (m_limitState == e_atLowerLimit)
		{
			float32 C = angle - m_lowerAngle;
			angularError = -C;

			// Prevent large angular corrections and allow some slop.
			C = b2Clamp(C + b2_angularSlop, -b2_maxAngularCorrection, 0.0f);
			limitImpulse = -m_motorMass * C;
		}
		else if (m_limitState == e_atUpperLimit)
		{
			float32 C = angle - m_upperAngle;
			angularError = C;

			// Prevent large angular corrections and allow some slop.
			C = b2Clamp(C - b2_angularSlop, 0.0f, b2_maxAngularCorrection);
			limitImpulse = -m_motorMass * C;
		}

		aA -= m_invIA * limitImpulse;
		aB += m_invIB * limitImpulse;
	}

	// Solve point-to-point constraint.
	{
		qA.Set(aA);
		qB.Set(aB);
		b2Vec2 rA = b2Mul(qA, m_localAnchorA - m_localCenterA);
		b2Vec2 rB = b2Mul(qB, m_localAnchorB - m_localCenterB);

		b2Vec2 C = cB + rB - cA - rA;
		positionError = C.Length();

		float32 mA = m_invMassA, mB = m_invMassB;
		float32 iA = m_invIA, iB = m_invIB;

		b2Mat22 K;
		K.ex.x = mA + mB + iA * rA.y * rA.y + iB * rB.y * rB.y;
		K.ex.y = -iA * rA.x * rA.y - iB * rB.x * rB.y;
		K.ey.x = K.ex.y;
		K.ey.y = mA + mB + iA * rA.x * rA.x + iB * rB.x * rB.x;

		b2Vec2 impulse = -K.Solve(C);

		cA -= mA * impulse;
		aA -= iA * b2Cross(rA, impulse);

		cB += mB * impulse;
		aB += iB * b2Cross(rB, impulse);
	}

	data.positions[m_indexA].c = cA;
	data.positions[m_indexA].a = aA;
	data.positions[m_indexB].c = cB;
	data.positions[m_indexB].a = aB;
	
	return positionError <= b2_linearSlop && angularError <= b2_angularSlop;
}

b2Vec2 b2RevoluteJoint::GetAnchorA() const
{
	return m_bodyA->GetWorldPoint(m_localAnchorA);
}

b2Vec2 b2RevoluteJoint::GetAnchorB() const
{
	return m_bodyB->GetWorldPoint(m_localAnchorB);
}

b2Vec2 b2RevoluteJoint::GetReactionForce(float32 inv_dt) const
{
	b2Vec2 P(m_impulse.x, m_impulse.y);
	return inv_dt * P;
}

float32 b2RevoluteJoint::GetReactionTorque(float32 inv_dt) const
{
	return inv_dt * m_impulse.z;
}

float32 b2RevoluteJoint::GetJointAngle() const
{
	b2Body* bA = m_bodyA;
	b2Body* bB = m_bodyB;
	return bB->m_sweep.a - bA->m_sweep.a - m_referenceAngle;
}

float32 b2RevoluteJoint::GetJointSpeed() const
{
	b2Body* bA = m_bodyA;
	b2Body* bB = m_bodyB;
	return bB->m_angularVelocity - bA->m_angularVelocity;
}

bool b2RevoluteJoint::IsMotorEnabled() const
{
	return m_enableMotor;
}

void b2RevoluteJoint::EnableMotor(bool flag)
{
	m_bodyA->SetAwake(true);
	m_bodyB->SetAwake(true);
	m_enableMotor = flag;
}

float32 b2RevoluteJoint::GetMotorTorque(float32 inv_dt) const
{
	return inv_dt * m_motorImpulse;
}

void b2RevoluteJoint::SetMotorSpeed(float32 speed)
{
	m_bodyA->SetAwake(true);
	m_bodyB->SetAwake(true);
	m_motorSpeed = speed;
}

void b2RevoluteJoint::SetMaxMotorTorque(float32 torque)
{
	m_bodyA->SetAwake(true);
	m_bodyB->SetAwake(true);
	m_maxMotorTorque = torque;
}

bool b2RevoluteJoint::IsLimitEnabled() const
{
	return m_enableLimit;
}

void b2RevoluteJoint::EnableLimit(bool flag)
{
	if (flag != m_enableLimit)
	{
		m_bodyA->SetAwake(true);
		m_bodyB->SetAwake(true);
		m_enableLimit = flag;
		m_impulse.z = 0.0f;
	}
}

float32 b2RevoluteJoint::GetLowerLimit() const
{
	return m_lowerAngle;
}

float32 b2RevoluteJoint::GetUpperLimit() const
{
	return m_upperAngle;
}

void b2RevoluteJoint::SetLimits(float32 lower, float32 upper)
{
	b2Assert(lower <= upper);
	
	if (lower != m_lowerAngle || upper != m_upperAngle)
	{
		m_bodyA->SetAwake(true);
		m_bodyB->SetAwake(true);
		m_impulse.z = 0.0f;
		m_lowerAngle = lower;
		m_upperAngle = upper;
	}
}

void b2RevoluteJoint::Dump()
{
	int32 indexA = m_bodyA->m_islandIndex;
	int32 indexB = m_bodyB->m_islandIndex;

	b2Log("  b2RevoluteJointDef jd;\n");
	b2Log("  jd.bodyA = bodies[%d];\n", indexA);
	b2Log("  jd.bodyB = bodies[%d];\n", indexB);
	b2Log("  jd.collideConnected = bool(%d);\n", m_collideConnected);
	b2Log("  jd.localAnchorA.Set(%.15lef, %.15lef);\n", m_localAnchorA.x, m_localAnchorA.y);
	b2Log("  jd.localAnchorB.Set(%.15lef, %.15lef);\n", m_localAnchorB.x, m_localAnchorB.y);
	b2Log("  jd.referenceAngle = %.15lef;\n", m_referenceAngle);
	b2Log("  jd.enableLimit = bool(%d);\n", m_enableLimit);
	b2Log("  jd.lowerAngle = %.15lef;\n", m_lowerAngle);
	b2Log("  jd.upperAngle = %.15lef;\n", m_upperAngle);
	b2Log("  jd.enableMotor = bool(%d);\n", m_enableMotor);
	b2Log("  jd.motorSpeed = %.15lef;\n", m_motorSpeed);
	b2Log("  jd.maxMotorTorque = %.15lef;\n", m_maxMotorTorque);
	b2Log("  joints[%d] = m_world->CreateJoint(&jd);\n", m_index);
}

// end of RevoluteJoint.cpp

// Gear Joint:
// C0 = (coordinate1 + ratio * coordinate2)_initial
// C = (coordinate1 + ratio * coordinate2) - C0 = 0
// J = [J1 ratio * J2]
// K = J * invM * JT
//   = J1 * invM1 * J1T + ratio * ratio * J2 * invM2 * J2T
//
// Revolute:
// coordinate = rotation
// Cdot = angularVelocity
// J = [0 0 1]
// K = J * invM * JT = invI
//
// Prismatic:
// coordinate = dot(p - pg, ug)
// Cdot = dot(v + cross(w, r), ug)
// J = [ug cross(r, ug)]
// K = J * invM * JT = invMass + invI * cross(r, ug)^2

b2GearJoint::b2GearJoint(const b2GearJointDef* def)
: b2Joint(def)
{
	m_joint1 = def->joint1;
	m_joint2 = def->joint2;

	m_typeA = m_joint1->GetType();
	m_typeB = m_joint2->GetType();

	b2Assert(m_typeA == e_revoluteJoint || m_typeA == e_prismaticJoint);
	b2Assert(m_typeB == e_revoluteJoint || m_typeB == e_prismaticJoint);

	float32 coordinateA, coordinateB;

	// TODO_ERIN there might be some problem with the joint edges in b2Joint.

	m_bodyC = m_joint1->GetBodyA();
	m_bodyA = m_joint1->GetBodyB();

	// Get geometry of joint1
	b2Transform xfA = m_bodyA->m_xf;
	float32 aA = m_bodyA->m_sweep.a;
	b2Transform xfC = m_bodyC->m_xf;
	float32 aC = m_bodyC->m_sweep.a;

	if (m_typeA == e_revoluteJoint)
	{
		b2RevoluteJoint* revolute = (b2RevoluteJoint*)def->joint1;
		m_localAnchorC = revolute->m_localAnchorA;
		m_localAnchorA = revolute->m_localAnchorB;
		m_referenceAngleA = revolute->m_referenceAngle;
		m_localAxisC.SetZero();

		coordinateA = aA - aC - m_referenceAngleA;
	}
	else
	{
		b2PrismaticJoint* prismatic = (b2PrismaticJoint*)def->joint1;
		m_localAnchorC = prismatic->m_localAnchorA;
		m_localAnchorA = prismatic->m_localAnchorB;
		m_referenceAngleA = prismatic->m_referenceAngle;
		m_localAxisC = prismatic->m_localXAxisA;

		b2Vec2 pC = m_localAnchorC;
		b2Vec2 pA = b2MulT(xfC.q, b2Mul(xfA.q, m_localAnchorA) + (xfA.p - xfC.p));
		coordinateA = b2Dot(pA - pC, m_localAxisC);
	}

	m_bodyD = m_joint2->GetBodyA();
	m_bodyB = m_joint2->GetBodyB();

	// Get geometry of joint2
	b2Transform xfB = m_bodyB->m_xf;
	float32 aB = m_bodyB->m_sweep.a;
	b2Transform xfD = m_bodyD->m_xf;
	float32 aD = m_bodyD->m_sweep.a;

	if (m_typeB == e_revoluteJoint)
	{
		b2RevoluteJoint* revolute = (b2RevoluteJoint*)def->joint2;
		m_localAnchorD = revolute->m_localAnchorA;
		m_localAnchorB = revolute->m_localAnchorB;
		m_referenceAngleB = revolute->m_referenceAngle;
		m_localAxisD.SetZero();

		coordinateB = aB - aD - m_referenceAngleB;
	}
	else
	{
		b2PrismaticJoint* prismatic = (b2PrismaticJoint*)def->joint2;
		m_localAnchorD = prismatic->m_localAnchorA;
		m_localAnchorB = prismatic->m_localAnchorB;
		m_referenceAngleB = prismatic->m_referenceAngle;
		m_localAxisD = prismatic->m_localXAxisA;

		b2Vec2 pD = m_localAnchorD;
		b2Vec2 pB = b2MulT(xfD.q, b2Mul(xfB.q, m_localAnchorB) + (xfB.p - xfD.p));
		coordinateB = b2Dot(pB - pD, m_localAxisD);
	}

	m_ratio = def->ratio;

	m_constant = coordinateA + m_ratio * coordinateB;

	m_impulse = 0.0f;
}

void b2GearJoint::InitVelocityConstraints(const b2SolverData& data)
{
	m_indexA = m_bodyA->m_islandIndex;
	m_indexB = m_bodyB->m_islandIndex;
	m_indexC = m_bodyC->m_islandIndex;
	m_indexD = m_bodyD->m_islandIndex;
	m_lcA = m_bodyA->m_sweep.localCenter;
	m_lcB = m_bodyB->m_sweep.localCenter;
	m_lcC = m_bodyC->m_sweep.localCenter;
	m_lcD = m_bodyD->m_sweep.localCenter;
	m_mA = m_bodyA->m_invMass;
	m_mB = m_bodyB->m_invMass;
	m_mC = m_bodyC->m_invMass;
	m_mD = m_bodyD->m_invMass;
	m_iA = m_bodyA->m_invI;
	m_iB = m_bodyB->m_invI;
	m_iC = m_bodyC->m_invI;
	m_iD = m_bodyD->m_invI;

	float32 aA = data.positions[m_indexA].a;
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;

	float32 aB = data.positions[m_indexB].a;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	float32 aC = data.positions[m_indexC].a;
	b2Vec2 vC = data.velocities[m_indexC].v;
	float32 wC = data.velocities[m_indexC].w;

	float32 aD = data.positions[m_indexD].a;
	b2Vec2 vD = data.velocities[m_indexD].v;
	float32 wD = data.velocities[m_indexD].w;

	b2Rot qA(aA), qB(aB), qC(aC), qD(aD);

	m_mass = 0.0f;

	if (m_typeA == e_revoluteJoint)
	{
		m_JvAC.SetZero();
		m_JwA = 1.0f;
		m_JwC = 1.0f;
		m_mass += m_iA + m_iC;
	}
	else
	{
		b2Vec2 u = b2Mul(qC, m_localAxisC);
		b2Vec2 rC = b2Mul(qC, m_localAnchorC - m_lcC);
		b2Vec2 rA = b2Mul(qA, m_localAnchorA - m_lcA);
		m_JvAC = u;
		m_JwC = b2Cross(rC, u);
		m_JwA = b2Cross(rA, u);
		m_mass += m_mC + m_mA + m_iC * m_JwC * m_JwC + m_iA * m_JwA * m_JwA;
	}

	if (m_typeB == e_revoluteJoint)
	{
		m_JvBD.SetZero();
		m_JwB = m_ratio;
		m_JwD = m_ratio;
		m_mass += m_ratio * m_ratio * (m_iB + m_iD);
	}
	else
	{
		b2Vec2 u = b2Mul(qD, m_localAxisD);
		b2Vec2 rD = b2Mul(qD, m_localAnchorD - m_lcD);
		b2Vec2 rB = b2Mul(qB, m_localAnchorB - m_lcB);
		m_JvBD = m_ratio * u;
		m_JwD = m_ratio * b2Cross(rD, u);
		m_JwB = m_ratio * b2Cross(rB, u);
		m_mass += m_ratio * m_ratio * (m_mD + m_mB) + m_iD * m_JwD * m_JwD + m_iB * m_JwB * m_JwB;
	}

	// Compute effective mass.
	m_mass = m_mass > 0.0f ? 1.0f / m_mass : 0.0f;

	if (data.step.warmStarting)
	{
		vA += (m_mA * m_impulse) * m_JvAC;
		wA += m_iA * m_impulse * m_JwA;
		vB += (m_mB * m_impulse) * m_JvBD;
		wB += m_iB * m_impulse * m_JwB;
		vC -= (m_mC * m_impulse) * m_JvAC;
		wC -= m_iC * m_impulse * m_JwC;
		vD -= (m_mD * m_impulse) * m_JvBD;
		wD -= m_iD * m_impulse * m_JwD;
	}
	else
	{
		m_impulse = 0.0f;
	}

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
	data.velocities[m_indexC].v = vC;
	data.velocities[m_indexC].w = wC;
	data.velocities[m_indexD].v = vD;
	data.velocities[m_indexD].w = wD;
}

void b2GearJoint::SolveVelocityConstraints(const b2SolverData& data)
{
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;
	b2Vec2 vC = data.velocities[m_indexC].v;
	float32 wC = data.velocities[m_indexC].w;
	b2Vec2 vD = data.velocities[m_indexD].v;
	float32 wD = data.velocities[m_indexD].w;

	float32 Cdot = b2Dot(m_JvAC, vA - vC) + b2Dot(m_JvBD, vB - vD);
	Cdot += (m_JwA * wA - m_JwC * wC) + (m_JwB * wB - m_JwD * wD);

	float32 impulse = -m_mass * Cdot;
	m_impulse += impulse;

	vA += (m_mA * impulse) * m_JvAC;
	wA += m_iA * impulse * m_JwA;
	vB += (m_mB * impulse) * m_JvBD;
	wB += m_iB * impulse * m_JwB;
	vC -= (m_mC * impulse) * m_JvAC;
	wC -= m_iC * impulse * m_JwC;
	vD -= (m_mD * impulse) * m_JvBD;
	wD -= m_iD * impulse * m_JwD;

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
	data.velocities[m_indexC].v = vC;
	data.velocities[m_indexC].w = wC;
	data.velocities[m_indexD].v = vD;
	data.velocities[m_indexD].w = wD;
}

bool b2GearJoint::SolvePositionConstraints(const b2SolverData& data)
{
	b2Vec2 cA = data.positions[m_indexA].c;
	float32 aA = data.positions[m_indexA].a;
	b2Vec2 cB = data.positions[m_indexB].c;
	float32 aB = data.positions[m_indexB].a;
	b2Vec2 cC = data.positions[m_indexC].c;
	float32 aC = data.positions[m_indexC].a;
	b2Vec2 cD = data.positions[m_indexD].c;
	float32 aD = data.positions[m_indexD].a;

	b2Rot qA(aA), qB(aB), qC(aC), qD(aD);

	float32 linearError = 0.0f;

	float32 coordinateA, coordinateB;

	b2Vec2 JvAC, JvBD;
	float32 JwA, JwB, JwC, JwD;
	float32 mass = 0.0f;

	if (m_typeA == e_revoluteJoint)
	{
		JvAC.SetZero();
		JwA = 1.0f;
		JwC = 1.0f;
		mass += m_iA + m_iC;

		coordinateA = aA - aC - m_referenceAngleA;
	}
	else
	{
		b2Vec2 u = b2Mul(qC, m_localAxisC);
		b2Vec2 rC = b2Mul(qC, m_localAnchorC - m_lcC);
		b2Vec2 rA = b2Mul(qA, m_localAnchorA - m_lcA);
		JvAC = u;
		JwC = b2Cross(rC, u);
		JwA = b2Cross(rA, u);
		mass += m_mC + m_mA + m_iC * JwC * JwC + m_iA * JwA * JwA;

		b2Vec2 pC = m_localAnchorC - m_lcC;
		b2Vec2 pA = b2MulT(qC, rA + (cA - cC));
		coordinateA = b2Dot(pA - pC, m_localAxisC);
	}

	if (m_typeB == e_revoluteJoint)
	{
		JvBD.SetZero();
		JwB = m_ratio;
		JwD = m_ratio;
		mass += m_ratio * m_ratio * (m_iB + m_iD);

		coordinateB = aB - aD - m_referenceAngleB;
	}
	else
	{
		b2Vec2 u = b2Mul(qD, m_localAxisD);
		b2Vec2 rD = b2Mul(qD, m_localAnchorD - m_lcD);
		b2Vec2 rB = b2Mul(qB, m_localAnchorB - m_lcB);
		JvBD = m_ratio * u;
		JwD = m_ratio * b2Cross(rD, u);
		JwB = m_ratio * b2Cross(rB, u);
		mass += m_ratio * m_ratio * (m_mD + m_mB) + m_iD * JwD * JwD + m_iB * JwB * JwB;

		b2Vec2 pD = m_localAnchorD - m_lcD;
		b2Vec2 pB = b2MulT(qD, rB + (cB - cD));
		coordinateB = b2Dot(pB - pD, m_localAxisD);
	}

	float32 C = (coordinateA + m_ratio * coordinateB) - m_constant;

	float32 impulse = 0.0f;
	if (mass > 0.0f)
	{
		impulse = -C / mass;
	}

	cA += m_mA * impulse * JvAC;
	aA += m_iA * impulse * JwA;
	cB += m_mB * impulse * JvBD;
	aB += m_iB * impulse * JwB;
	cC -= m_mC * impulse * JvAC;
	aC -= m_iC * impulse * JwC;
	cD -= m_mD * impulse * JvBD;
	aD -= m_iD * impulse * JwD;

	data.positions[m_indexA].c = cA;
	data.positions[m_indexA].a = aA;
	data.positions[m_indexB].c = cB;
	data.positions[m_indexB].a = aB;
	data.positions[m_indexC].c = cC;
	data.positions[m_indexC].a = aC;
	data.positions[m_indexD].c = cD;
	data.positions[m_indexD].a = aD;

	// TODO_ERIN not implemented
	return linearError < b2_linearSlop;
}

b2Vec2 b2GearJoint::GetAnchorA() const
{
	return m_bodyA->GetWorldPoint(m_localAnchorA);
}

b2Vec2 b2GearJoint::GetAnchorB() const
{
	return m_bodyB->GetWorldPoint(m_localAnchorB);
}

b2Vec2 b2GearJoint::GetReactionForce(float32 inv_dt) const
{
	b2Vec2 P = m_impulse * m_JvAC;
	return inv_dt * P;
}

float32 b2GearJoint::GetReactionTorque(float32 inv_dt) const
{
	float32 L = m_impulse * m_JwA;
	return inv_dt * L;
}

void b2GearJoint::SetRatio(float32 ratio)
{
	b2Assert(b2IsValid(ratio));
	m_ratio = ratio;
}

float32 b2GearJoint::GetRatio() const
{
	return m_ratio;
}

void b2GearJoint::Dump()
{
	int32 indexA = m_bodyA->m_islandIndex;
	int32 indexB = m_bodyB->m_islandIndex;

	int32 index1 = m_joint1->m_index;
	int32 index2 = m_joint2->m_index;

	b2Log("  b2GearJointDef jd;\n");
	b2Log("  jd.bodyA = bodies[%d];\n", indexA);
	b2Log("  jd.bodyB = bodies[%d];\n", indexB);
	b2Log("  jd.collideConnected = bool(%d);\n", m_collideConnected);
	b2Log("  jd.joint1 = joints[%d];\n", index1);
	b2Log("  jd.joint2 = joints[%d];\n", index2);
	b2Log("  jd.ratio = %.15lef;\n", m_ratio);
	b2Log("  joints[%d] = m_world->CreateJoint(&jd);\n", m_index);
}

// end of GearJoint.cpp

/// Wheel joint definition. This requires defining a line of
/// motion using an axis and an anchor point. The definition uses local
/// anchor points and a local axis so that the initial configuration
/// can violate the constraint slightly. The joint translation is zero
/// when the local anchor points coincide in world space. Using local
/// anchors and a local axis helps when saving and loading a game.
struct b2WheelJointDef : public b2JointDef
{
	b2WheelJointDef()
	{
		type = e_wheelJoint;
		localAnchorA.SetZero();
		localAnchorB.SetZero();
		localAxisA.Set(1.0f, 0.0f);
		enableMotor = false;
		maxMotorTorque = 0.0f;
		motorSpeed = 0.0f;
		frequencyHz = 2.0f;
		dampingRatio = 0.7f;
	}

	/// Initialize the bodies, anchors, axis, and reference angle using the world
	/// anchor and world axis.
	void Initialize(b2Body* bodyA, b2Body* bodyB, const b2Vec2& anchor, const b2Vec2& axis);

	/// The local anchor point relative to bodyA's origin.
	b2Vec2 localAnchorA;

	/// The local anchor point relative to bodyB's origin.
	b2Vec2 localAnchorB;

	/// The local translation axis in bodyA.
	b2Vec2 localAxisA;

	/// Enable/disable the joint motor.
	bool enableMotor;

	/// The maximum motor torque, usually in N-m.
	float32 maxMotorTorque;

	/// The desired motor speed in radians per second.
	float32 motorSpeed;

	/// Suspension frequency, zero indicates no suspension
	float32 frequencyHz;

	/// Suspension damping ratio, one indicates critical damping
	float32 dampingRatio;
};

/// A wheel joint. This joint provides two degrees of freedom: translation
/// along an axis fixed in bodyA and rotation in the plane. You can use a
/// joint limit to restrict the range of motion and a joint motor to drive
/// the rotation or to model rotational friction.
/// This joint is designed for vehicle suspensions.
class b2WheelJoint : public b2Joint
{
public:
	b2Vec2 GetAnchorA() const;
	b2Vec2 GetAnchorB() const;

	b2Vec2 GetReactionForce(float32 inv_dt) const;
	float32 GetReactionTorque(float32 inv_dt) const;

	/// The local anchor point relative to bodyA's origin.
	const b2Vec2& GetLocalAnchorA() const { return m_localAnchorA; }

	/// The local anchor point relative to bodyB's origin.
	const b2Vec2& GetLocalAnchorB() const  { return m_localAnchorB; }

	/// The local joint axis relative to bodyA.
	const b2Vec2& GetLocalAxisA() const { return m_localXAxisA; }

	/// Get the current joint translation, usually in meters.
	float32 GetJointTranslation() const;

	/// Get the current joint translation speed, usually in meters per second.
	float32 GetJointSpeed() const;

	/// Is the joint motor enabled?
	bool IsMotorEnabled() const;

	/// Enable/disable the joint motor.
	void EnableMotor(bool flag);

	/// Set the motor speed, usually in radians per second.
	void SetMotorSpeed(float32 speed);

	/// Get the motor speed, usually in radians per second.
	float32 GetMotorSpeed() const;

	/// Set/Get the maximum motor force, usually in N-m.
	void SetMaxMotorTorque(float32 torque);
	float32 GetMaxMotorTorque() const;

	/// Get the current motor torque given the inverse time step, usually in N-m.
	float32 GetMotorTorque(float32 inv_dt) const;

	/// Set/Get the spring frequency in hertz. Setting the frequency to zero disables the spring.
	void SetSpringFrequencyHz(float32 hz);
	float32 GetSpringFrequencyHz() const;

	/// Set/Get the spring damping ratio
	void SetSpringDampingRatio(float32 ratio);
	float32 GetSpringDampingRatio() const;

	/// Dump to b2Log
	void Dump();

protected:

	friend class b2Joint;
	b2WheelJoint(const b2WheelJointDef* def);

	void InitVelocityConstraints(const b2SolverData& data);
	void SolveVelocityConstraints(const b2SolverData& data);
	bool SolvePositionConstraints(const b2SolverData& data);

	float32 m_frequencyHz;
	float32 m_dampingRatio;

	// Solver shared
	b2Vec2 m_localAnchorA;
	b2Vec2 m_localAnchorB;
	b2Vec2 m_localXAxisA;
	b2Vec2 m_localYAxisA;

	float32 m_impulse;
	float32 m_motorImpulse;
	float32 m_springImpulse;

	float32 m_maxMotorTorque;
	float32 m_motorSpeed;
	bool m_enableMotor;

	// Solver temp
	int32 m_indexA;
	int32 m_indexB;
	b2Vec2 m_localCenterA;
	b2Vec2 m_localCenterB;
	float32 m_invMassA;
	float32 m_invMassB;
	float32 m_invIA;
	float32 m_invIB;

	b2Vec2 m_ax, m_ay;
	float32 m_sAx, m_sBx;
	float32 m_sAy, m_sBy;

	float32 m_mass;
	float32 m_motorMass;
	float32 m_springMass;

	float32 m_bias;
	float32 m_gamma;
};

inline float32 b2WheelJoint::GetMotorSpeed() const
{
	return m_motorSpeed;
}

inline float32 b2WheelJoint::GetMaxMotorTorque() const
{
	return m_maxMotorTorque;
}

inline void b2WheelJoint::SetSpringFrequencyHz(float32 hz)
{
	m_frequencyHz = hz;
}

inline float32 b2WheelJoint::GetSpringFrequencyHz() const
{
	return m_frequencyHz;
}

inline void b2WheelJoint::SetSpringDampingRatio(float32 ratio)
{
	m_dampingRatio = ratio;
}

inline float32 b2WheelJoint::GetSpringDampingRatio() const
{
	return m_dampingRatio;
}

// end of WheelJoint.h

// Linear constraint (point-to-line)
// d = pB - pA = xB + rB - xA - rA
// C = dot(ay, d)
// Cdot = dot(d, cross(wA, ay)) + dot(ay, vB + cross(wB, rB) - vA - cross(wA, rA))
//      = -dot(ay, vA) - dot(cross(d + rA, ay), wA) + dot(ay, vB) + dot(cross(rB, ay), vB)
// J = [-ay, -cross(d + rA, ay), ay, cross(rB, ay)]

// Spring linear constraint
// C = dot(ax, d)
// Cdot = = -dot(ax, vA) - dot(cross(d + rA, ax), wA) + dot(ax, vB) + dot(cross(rB, ax), vB)
// J = [-ax -cross(d+rA, ax) ax cross(rB, ax)]

// Motor rotational constraint
// Cdot = wB - wA
// J = [0 0 -1 0 0 1]

void b2WheelJointDef::Initialize(b2Body* bA, b2Body* bB, const b2Vec2& anchor, const b2Vec2& axis)
{
	bodyA = bA;
	bodyB = bB;
	localAnchorA = bodyA->GetLocalPoint(anchor);
	localAnchorB = bodyB->GetLocalPoint(anchor);
	localAxisA = bodyA->GetLocalVector(axis);
}

b2WheelJoint::b2WheelJoint(const b2WheelJointDef* def)
: b2Joint(def)
{
	m_localAnchorA = def->localAnchorA;
	m_localAnchorB = def->localAnchorB;
	m_localXAxisA = def->localAxisA;
	m_localYAxisA = b2Cross(1.0f, m_localXAxisA);

	m_mass = 0.0f;
	m_impulse = 0.0f;
	m_motorMass = 0.0f;
	m_motorImpulse = 0.0f;
	m_springMass = 0.0f;
	m_springImpulse = 0.0f;

	m_maxMotorTorque = def->maxMotorTorque;
	m_motorSpeed = def->motorSpeed;
	m_enableMotor = def->enableMotor;

	m_frequencyHz = def->frequencyHz;
	m_dampingRatio = def->dampingRatio;

	m_bias = 0.0f;
	m_gamma = 0.0f;

	m_ax.SetZero();
	m_ay.SetZero();
}

void b2WheelJoint::InitVelocityConstraints(const b2SolverData& data)
{
	m_indexA = m_bodyA->m_islandIndex;
	m_indexB = m_bodyB->m_islandIndex;
	m_localCenterA = m_bodyA->m_sweep.localCenter;
	m_localCenterB = m_bodyB->m_sweep.localCenter;
	m_invMassA = m_bodyA->m_invMass;
	m_invMassB = m_bodyB->m_invMass;
	m_invIA = m_bodyA->m_invI;
	m_invIB = m_bodyB->m_invI;

	float32 mA = m_invMassA, mB = m_invMassB;
	float32 iA = m_invIA, iB = m_invIB;

	b2Vec2 cA = data.positions[m_indexA].c;
	float32 aA = data.positions[m_indexA].a;
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;

	b2Vec2 cB = data.positions[m_indexB].c;
	float32 aB = data.positions[m_indexB].a;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	b2Rot qA(aA), qB(aB);

	// Compute the effective masses.
	b2Vec2 rA = b2Mul(qA, m_localAnchorA - m_localCenterA);
	b2Vec2 rB = b2Mul(qB, m_localAnchorB - m_localCenterB);
	b2Vec2 d = cB + rB - cA - rA;

	// Point to line constraint
	{
		m_ay = b2Mul(qA, m_localYAxisA);
		m_sAy = b2Cross(d + rA, m_ay);
		m_sBy = b2Cross(rB, m_ay);

		m_mass = mA + mB + iA * m_sAy * m_sAy + iB * m_sBy * m_sBy;

		if (m_mass > 0.0f)
		{
			m_mass = 1.0f / m_mass;
		}
	}

	// Spring constraint
	m_springMass = 0.0f;
	m_bias = 0.0f;
	m_gamma = 0.0f;
	if (m_frequencyHz > 0.0f)
	{
		m_ax = b2Mul(qA, m_localXAxisA);
		m_sAx = b2Cross(d + rA, m_ax);
		m_sBx = b2Cross(rB, m_ax);

		float32 invMass = mA + mB + iA * m_sAx * m_sAx + iB * m_sBx * m_sBx;

		if (invMass > 0.0f)
		{
			m_springMass = 1.0f / invMass;

			float32 C = b2Dot(d, m_ax);

			// Frequency
			float32 omega = 2.0f * b2_pi * m_frequencyHz;

			// Damping coefficient
			float32 d = 2.0f * m_springMass * m_dampingRatio * omega;

			// Spring stiffness
			float32 k = m_springMass * omega * omega;

			// magic formulas
			float32 h = data.step.dt;
			m_gamma = h * (d + h * k);
			if (m_gamma > 0.0f)
			{
				m_gamma = 1.0f / m_gamma;
			}

			m_bias = C * h * k * m_gamma;

			m_springMass = invMass + m_gamma;
			if (m_springMass > 0.0f)
			{
				m_springMass = 1.0f / m_springMass;
			}
		}
	}
	else
	{
		m_springImpulse = 0.0f;
	}

	// Rotational motor
	if (m_enableMotor)
	{
		m_motorMass = iA + iB;
		if (m_motorMass > 0.0f)
		{
			m_motorMass = 1.0f / m_motorMass;
		}
	}
	else
	{
		m_motorMass = 0.0f;
		m_motorImpulse = 0.0f;
	}

	if (data.step.warmStarting)
	{
		// Account for variable time step.
		m_impulse *= data.step.dtRatio;
		m_springImpulse *= data.step.dtRatio;
		m_motorImpulse *= data.step.dtRatio;

		b2Vec2 P = m_impulse * m_ay + m_springImpulse * m_ax;
		float32 LA = m_impulse * m_sAy + m_springImpulse * m_sAx + m_motorImpulse;
		float32 LB = m_impulse * m_sBy + m_springImpulse * m_sBx + m_motorImpulse;

		vA -= m_invMassA * P;
		wA -= m_invIA * LA;

		vB += m_invMassB * P;
		wB += m_invIB * LB;
	}
	else
	{
		m_impulse = 0.0f;
		m_springImpulse = 0.0f;
		m_motorImpulse = 0.0f;
	}

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

void b2WheelJoint::SolveVelocityConstraints(const b2SolverData& data)
{
	float32 mA = m_invMassA, mB = m_invMassB;
	float32 iA = m_invIA, iB = m_invIB;

	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	// Solve spring constraint
	{
		float32 Cdot = b2Dot(m_ax, vB - vA) + m_sBx * wB - m_sAx * wA;
		float32 impulse = -m_springMass * (Cdot + m_bias + m_gamma * m_springImpulse);
		m_springImpulse += impulse;

		b2Vec2 P = impulse * m_ax;
		float32 LA = impulse * m_sAx;
		float32 LB = impulse * m_sBx;

		vA -= mA * P;
		wA -= iA * LA;

		vB += mB * P;
		wB += iB * LB;
	}

	// Solve rotational motor constraint
	{
		float32 Cdot = wB - wA - m_motorSpeed;
		float32 impulse = -m_motorMass * Cdot;

		float32 oldImpulse = m_motorImpulse;
		float32 maxImpulse = data.step.dt * m_maxMotorTorque;
		m_motorImpulse = b2Clamp(m_motorImpulse + impulse, -maxImpulse, maxImpulse);
		impulse = m_motorImpulse - oldImpulse;

		wA -= iA * impulse;
		wB += iB * impulse;
	}

	// Solve point to line constraint
	{
		float32 Cdot = b2Dot(m_ay, vB - vA) + m_sBy * wB - m_sAy * wA;
		float32 impulse = -m_mass * Cdot;
		m_impulse += impulse;

		b2Vec2 P = impulse * m_ay;
		float32 LA = impulse * m_sAy;
		float32 LB = impulse * m_sBy;

		vA -= mA * P;
		wA -= iA * LA;

		vB += mB * P;
		wB += iB * LB;
	}

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

bool b2WheelJoint::SolvePositionConstraints(const b2SolverData& data)
{
	b2Vec2 cA = data.positions[m_indexA].c;
	float32 aA = data.positions[m_indexA].a;
	b2Vec2 cB = data.positions[m_indexB].c;
	float32 aB = data.positions[m_indexB].a;

	b2Rot qA(aA), qB(aB);

	b2Vec2 rA = b2Mul(qA, m_localAnchorA - m_localCenterA);
	b2Vec2 rB = b2Mul(qB, m_localAnchorB - m_localCenterB);
	b2Vec2 d = (cB - cA) + rB - rA;

	b2Vec2 ay = b2Mul(qA, m_localYAxisA);

	float32 sAy = b2Cross(d + rA, ay);
	float32 sBy = b2Cross(rB, ay);

	float32 C = b2Dot(d, ay);

	float32 k = m_invMassA + m_invMassB + m_invIA * m_sAy * m_sAy + m_invIB * m_sBy * m_sBy;

	float32 impulse;
	if (k != 0.0f)
	{
		impulse = - C / k;
	}
	else
	{
		impulse = 0.0f;
	}

	b2Vec2 P = impulse * ay;
	float32 LA = impulse * sAy;
	float32 LB = impulse * sBy;

	cA -= m_invMassA * P;
	aA -= m_invIA * LA;
	cB += m_invMassB * P;
	aB += m_invIB * LB;

	data.positions[m_indexA].c = cA;
	data.positions[m_indexA].a = aA;
	data.positions[m_indexB].c = cB;
	data.positions[m_indexB].a = aB;

	return b2Abs(C) <= b2_linearSlop;
}

b2Vec2 b2WheelJoint::GetAnchorA() const
{
	return m_bodyA->GetWorldPoint(m_localAnchorA);
}

b2Vec2 b2WheelJoint::GetAnchorB() const
{
	return m_bodyB->GetWorldPoint(m_localAnchorB);
}

b2Vec2 b2WheelJoint::GetReactionForce(float32 inv_dt) const
{
	return inv_dt * (m_impulse * m_ay + m_springImpulse * m_ax);
}

float32 b2WheelJoint::GetReactionTorque(float32 inv_dt) const
{
	return inv_dt * m_motorImpulse;
}

float32 b2WheelJoint::GetJointTranslation() const
{
	b2Body* bA = m_bodyA;
	b2Body* bB = m_bodyB;

	b2Vec2 pA = bA->GetWorldPoint(m_localAnchorA);
	b2Vec2 pB = bB->GetWorldPoint(m_localAnchorB);
	b2Vec2 d = pB - pA;
	b2Vec2 axis = bA->GetWorldVector(m_localXAxisA);

	float32 translation = b2Dot(d, axis);
	return translation;
}

float32 b2WheelJoint::GetJointSpeed() const
{
	float32 wA = m_bodyA->m_angularVelocity;
	float32 wB = m_bodyB->m_angularVelocity;
	return wB - wA;
}

bool b2WheelJoint::IsMotorEnabled() const
{
	return m_enableMotor;
}

void b2WheelJoint::EnableMotor(bool flag)
{
	m_bodyA->SetAwake(true);
	m_bodyB->SetAwake(true);
	m_enableMotor = flag;
}

void b2WheelJoint::SetMotorSpeed(float32 speed)
{
	m_bodyA->SetAwake(true);
	m_bodyB->SetAwake(true);
	m_motorSpeed = speed;
}

void b2WheelJoint::SetMaxMotorTorque(float32 torque)
{
	m_bodyA->SetAwake(true);
	m_bodyB->SetAwake(true);
	m_maxMotorTorque = torque;
}

float32 b2WheelJoint::GetMotorTorque(float32 inv_dt) const
{
	return inv_dt * m_motorImpulse;
}

void b2WheelJoint::Dump()
{
	int32 indexA = m_bodyA->m_islandIndex;
	int32 indexB = m_bodyB->m_islandIndex;

	b2Log("  b2WheelJointDef jd;\n");
	b2Log("  jd.bodyA = bodies[%d];\n", indexA);
	b2Log("  jd.bodyB = bodies[%d];\n", indexB);
	b2Log("  jd.collideConnected = bool(%d);\n", m_collideConnected);
	b2Log("  jd.localAnchorA.Set(%.15lef, %.15lef);\n", m_localAnchorA.x, m_localAnchorA.y);
	b2Log("  jd.localAnchorB.Set(%.15lef, %.15lef);\n", m_localAnchorB.x, m_localAnchorB.y);
	b2Log("  jd.localAxisA.Set(%.15lef, %.15lef);\n", m_localXAxisA.x, m_localXAxisA.y);
	b2Log("  jd.enableMotor = bool(%d);\n", m_enableMotor);
	b2Log("  jd.motorSpeed = %.15lef;\n", m_motorSpeed);
	b2Log("  jd.maxMotorTorque = %.15lef;\n", m_maxMotorTorque);
	b2Log("  jd.frequencyHz = %.15lef;\n", m_frequencyHz);
	b2Log("  jd.dampingRatio = %.15lef;\n", m_dampingRatio);
	b2Log("  joints[%d] = m_world->CreateJoint(&jd);\n", m_index);
}

// end of WheelJoint.cpp

/// Mouse joint definition. This requires a world target point,
/// tuning parameters, and the time step.
struct b2MouseJointDef : public b2JointDef
{
	b2MouseJointDef()
	{
		type = e_mouseJoint;
		target.Set(0.0f, 0.0f);
		maxForce = 0.0f;
		frequencyHz = 5.0f;
		dampingRatio = 0.7f;
	}

	/// The initial world target point. This is assumed
	/// to coincide with the body anchor initially.
	b2Vec2 target;

	/// The maximum constraint force that can be exerted
	/// to move the candidate body. Usually you will express
	/// as some multiple of the weight (multiplier * mass * gravity).
	float32 maxForce;

	/// The response speed.
	float32 frequencyHz;

	/// The damping ratio. 0 = no damping, 1 = critical damping.
	float32 dampingRatio;
};

/// A mouse joint is used to make a point on a body track a
/// specified world point. This a soft constraint with a maximum
/// force. This allows the constraint to stretch and without
/// applying huge forces.
/// NOTE: this joint is not documented in the manual because it was
/// developed to be used in the testbed. If you want to learn how to
/// use the mouse joint, look at the testbed.
class b2MouseJoint : public b2Joint
{
public:

	/// Implements b2Joint.
	b2Vec2 GetAnchorA() const;

	/// Implements b2Joint.
	b2Vec2 GetAnchorB() const;

	/// Implements b2Joint.
	b2Vec2 GetReactionForce(float32 inv_dt) const;

	/// Implements b2Joint.
	float32 GetReactionTorque(float32 inv_dt) const;

	/// Use this to update the target point.
	void SetTarget(const b2Vec2& target);
	const b2Vec2& GetTarget() const;

	/// Set/get the maximum force in Newtons.
	void SetMaxForce(float32 force);
	float32 GetMaxForce() const;

	/// Set/get the frequency in Hertz.
	void SetFrequency(float32 hz);
	float32 GetFrequency() const;

	/// Set/get the damping ratio (dimensionless).
	void SetDampingRatio(float32 ratio);
	float32 GetDampingRatio() const;

	/// The mouse joint does not support dumping.
	void Dump() { b2Log("Mouse joint dumping is not supported.\n"); }

	/// Implement b2Joint::ShiftOrigin
	void ShiftOrigin(const b2Vec2& newOrigin);

protected:
	friend class b2Joint;

	b2MouseJoint(const b2MouseJointDef* def);

	void InitVelocityConstraints(const b2SolverData& data);
	void SolveVelocityConstraints(const b2SolverData& data);
	bool SolvePositionConstraints(const b2SolverData& data);

	b2Vec2 m_localAnchorB;
	b2Vec2 m_targetA;
	float32 m_frequencyHz;
	float32 m_dampingRatio;
	float32 m_beta;
	
	// Solver shared
	b2Vec2 m_impulse;
	float32 m_maxForce;
	float32 m_gamma;

	// Solver temp
	int32 m_indexA;
	int32 m_indexB;
	b2Vec2 m_rB;
	b2Vec2 m_localCenterB;
	float32 m_invMassB;
	float32 m_invIB;
	b2Mat22 m_mass;
	b2Vec2 m_C;
};

// end of MouseJoint.h

// p = attached point, m = mouse point
// C = p - m
// Cdot = v
//      = v + cross(w, r)
// J = [I r_skew]
// Identity used:
// w k % (rx i + ry j) = w * (-ry i + rx j)

b2MouseJoint::b2MouseJoint(const b2MouseJointDef* def)
: b2Joint(def)
{
	b2Assert(def->target.IsValid());
	b2Assert(b2IsValid(def->maxForce) && def->maxForce >= 0.0f);
	b2Assert(b2IsValid(def->frequencyHz) && def->frequencyHz >= 0.0f);
	b2Assert(b2IsValid(def->dampingRatio) && def->dampingRatio >= 0.0f);

	m_targetA = def->target;
	m_localAnchorB = b2MulT(m_bodyB->GetTransform(), m_targetA);

	m_maxForce = def->maxForce;
	m_impulse.SetZero();

	m_frequencyHz = def->frequencyHz;
	m_dampingRatio = def->dampingRatio;

	m_beta = 0.0f;
	m_gamma = 0.0f;
}

void b2MouseJoint::SetTarget(const b2Vec2& target)
{
	if (m_bodyB->IsAwake() == false)
	{
		m_bodyB->SetAwake(true);
	}
	m_targetA = target;
}

const b2Vec2& b2MouseJoint::GetTarget() const
{
	return m_targetA;
}

void b2MouseJoint::SetMaxForce(float32 force)
{
	m_maxForce = force;
}

float32 b2MouseJoint::GetMaxForce() const
{
	return m_maxForce;
}

void b2MouseJoint::SetFrequency(float32 hz)
{
	m_frequencyHz = hz;
}

float32 b2MouseJoint::GetFrequency() const
{
	return m_frequencyHz;
}

void b2MouseJoint::SetDampingRatio(float32 ratio)
{
	m_dampingRatio = ratio;
}

float32 b2MouseJoint::GetDampingRatio() const
{
	return m_dampingRatio;
}

void b2MouseJoint::InitVelocityConstraints(const b2SolverData& data)
{
	m_indexB = m_bodyB->m_islandIndex;
	m_localCenterB = m_bodyB->m_sweep.localCenter;
	m_invMassB = m_bodyB->m_invMass;
	m_invIB = m_bodyB->m_invI;

	b2Vec2 cB = data.positions[m_indexB].c;
	float32 aB = data.positions[m_indexB].a;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	b2Rot qB(aB);

	float32 mass = m_bodyB->GetMass();

	// Frequency
	float32 omega = 2.0f * b2_pi * m_frequencyHz;

	// Damping coefficient
	float32 d = 2.0f * mass * m_dampingRatio * omega;

	// Spring stiffness
	float32 k = mass * (omega * omega);

	// magic formulas
	// gamma has units of inverse mass.
	// beta has units of inverse time.
	float32 h = data.step.dt;
	b2Assert(d + h * k > b2_epsilon);
	m_gamma = h * (d + h * k);
	if (m_gamma != 0.0f)
	{
		m_gamma = 1.0f / m_gamma;
	}
	m_beta = h * k * m_gamma;

	// Compute the effective mass matrix.
	m_rB = b2Mul(qB, m_localAnchorB - m_localCenterB);

	// K    = [(1/m1 + 1/m2) * eye(2) - skew(r1) * invI1 * skew(r1) - skew(r2) * invI2 * skew(r2)]
	//      = [1/m1+1/m2     0    ] + invI1 * [r1.y*r1.y -r1.x*r1.y] + invI2 * [r1.y*r1.y -r1.x*r1.y]
	//        [    0     1/m1+1/m2]           [-r1.x*r1.y r1.x*r1.x]           [-r1.x*r1.y r1.x*r1.x]
	b2Mat22 K;
	K.ex.x = m_invMassB + m_invIB * m_rB.y * m_rB.y + m_gamma;
	K.ex.y = -m_invIB * m_rB.x * m_rB.y;
	K.ey.x = K.ex.y;
	K.ey.y = m_invMassB + m_invIB * m_rB.x * m_rB.x + m_gamma;

	m_mass = K.GetInverse();

	m_C = cB + m_rB - m_targetA;
	m_C *= m_beta;

	// Cheat with some damping
	wB *= 0.98f;

	if (data.step.warmStarting)
	{
		m_impulse *= data.step.dtRatio;
		vB += m_invMassB * m_impulse;
		wB += m_invIB * b2Cross(m_rB, m_impulse);
	}
	else
	{
		m_impulse.SetZero();
	}

	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

void b2MouseJoint::SolveVelocityConstraints(const b2SolverData& data)
{
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	// Cdot = v + cross(w, r)
	b2Vec2 Cdot = vB + b2Cross(wB, m_rB);
	b2Vec2 impulse = b2Mul(m_mass, -(Cdot + m_C + m_gamma * m_impulse));

	b2Vec2 oldImpulse = m_impulse;
	m_impulse += impulse;
	float32 maxImpulse = data.step.dt * m_maxForce;
	if (m_impulse.LengthSquared() > maxImpulse * maxImpulse)
	{
		m_impulse *= maxImpulse / m_impulse.Length();
	}
	impulse = m_impulse - oldImpulse;

	vB += m_invMassB * impulse;
	wB += m_invIB * b2Cross(m_rB, impulse);

	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

bool b2MouseJoint::SolvePositionConstraints(const b2SolverData& data)
{
	B2_NOT_USED(data);
	return true;
}

b2Vec2 b2MouseJoint::GetAnchorA() const
{
	return m_targetA;
}

b2Vec2 b2MouseJoint::GetAnchorB() const
{
	return m_bodyB->GetWorldPoint(m_localAnchorB);
}

b2Vec2 b2MouseJoint::GetReactionForce(float32 inv_dt) const
{
	return inv_dt * m_impulse;
}

float32 b2MouseJoint::GetReactionTorque(float32 inv_dt) const
{
	return inv_dt * 0.0f;
}

void b2MouseJoint::ShiftOrigin(const b2Vec2& newOrigin)
{
	m_targetA -= newOrigin;
}

// end of MouseJoint.cpp


/// Weld joint definition. You need to specify local anchor points
/// where they are attached and the relative body angle. The position
/// of the anchor points is important for computing the reaction torque.
struct b2WeldJointDef : public b2JointDef
{
	b2WeldJointDef()
	{
		type = e_weldJoint;
		localAnchorA.Set(0.0f, 0.0f);
		localAnchorB.Set(0.0f, 0.0f);
		referenceAngle = 0.0f;
		frequencyHz = 0.0f;
		dampingRatio = 0.0f;
	}

	/// Initialize the bodies, anchors, and reference angle using a world
	/// anchor point.
	void Initialize(b2Body* bodyA, b2Body* bodyB, const b2Vec2& anchor);

	/// The local anchor point relative to bodyA's origin.
	b2Vec2 localAnchorA;

	/// The local anchor point relative to bodyB's origin.
	b2Vec2 localAnchorB;

	/// The bodyB angle minus bodyA angle in the reference state (radians).
	float32 referenceAngle;
	
	/// The mass-spring-damper frequency in Hertz. Rotation only.
	/// Disable softness with a value of 0.
	float32 frequencyHz;

	/// The damping ratio. 0 = no damping, 1 = critical damping.
	float32 dampingRatio;
};

/// A weld joint essentially glues two bodies together. A weld joint may
/// distort somewhat because the island constraint solver is approximate.
class b2WeldJoint : public b2Joint
{
public:
	b2Vec2 GetAnchorA() const;
	b2Vec2 GetAnchorB() const;

	b2Vec2 GetReactionForce(float32 inv_dt) const;
	float32 GetReactionTorque(float32 inv_dt) const;

	/// The local anchor point relative to bodyA's origin.
	const b2Vec2& GetLocalAnchorA() const { return m_localAnchorA; }

	/// The local anchor point relative to bodyB's origin.
	const b2Vec2& GetLocalAnchorB() const  { return m_localAnchorB; }

	/// Get the reference angle.
	float32 GetReferenceAngle() const { return m_referenceAngle; }

	/// Set/get frequency in Hz.
	void SetFrequency(float32 hz) { m_frequencyHz = hz; }
	float32 GetFrequency() const { return m_frequencyHz; }

	/// Set/get damping ratio.
	void SetDampingRatio(float32 ratio) { m_dampingRatio = ratio; }
	float32 GetDampingRatio() const { return m_dampingRatio; }

	/// Dump to b2Log
	void Dump();

protected:

	friend class b2Joint;

	b2WeldJoint(const b2WeldJointDef* def);

	void InitVelocityConstraints(const b2SolverData& data);
	void SolveVelocityConstraints(const b2SolverData& data);
	bool SolvePositionConstraints(const b2SolverData& data);

	float32 m_frequencyHz;
	float32 m_dampingRatio;
	float32 m_bias;

	// Solver shared
	b2Vec2 m_localAnchorA;
	b2Vec2 m_localAnchorB;
	float32 m_referenceAngle;
	float32 m_gamma;
	b2Vec3 m_impulse;

	// Solver temp
	int32 m_indexA;
	int32 m_indexB;
	b2Vec2 m_rA;
	b2Vec2 m_rB;
	b2Vec2 m_localCenterA;
	b2Vec2 m_localCenterB;
	float32 m_invMassA;
	float32 m_invMassB;
	float32 m_invIA;
	float32 m_invIB;
	b2Mat33 m_mass;
};

// end of WeldJoint.h

// Point-to-point constraint
// C = p2 - p1
// Cdot = v2 - v1
//      = v2 + cross(w2, r2) - v1 - cross(w1, r1)
// J = [-I -r1_skew I r2_skew ]
// Identity used:
// w k % (rx i + ry j) = w * (-ry i + rx j)

// Angle constraint
// C = angle2 - angle1 - referenceAngle
// Cdot = w2 - w1
// J = [0 0 -1 0 0 1]
// K = invI1 + invI2

void b2WeldJointDef::Initialize(b2Body* bA, b2Body* bB, const b2Vec2& anchor)
{
	bodyA = bA;
	bodyB = bB;
	localAnchorA = bodyA->GetLocalPoint(anchor);
	localAnchorB = bodyB->GetLocalPoint(anchor);
	referenceAngle = bodyB->GetAngle() - bodyA->GetAngle();
}

b2WeldJoint::b2WeldJoint(const b2WeldJointDef* def)
: b2Joint(def)
{
	m_localAnchorA = def->localAnchorA;
	m_localAnchorB = def->localAnchorB;
	m_referenceAngle = def->referenceAngle;
	m_frequencyHz = def->frequencyHz;
	m_dampingRatio = def->dampingRatio;

	m_impulse.SetZero();
}

void b2WeldJoint::InitVelocityConstraints(const b2SolverData& data)
{
	m_indexA = m_bodyA->m_islandIndex;
	m_indexB = m_bodyB->m_islandIndex;
	m_localCenterA = m_bodyA->m_sweep.localCenter;
	m_localCenterB = m_bodyB->m_sweep.localCenter;
	m_invMassA = m_bodyA->m_invMass;
	m_invMassB = m_bodyB->m_invMass;
	m_invIA = m_bodyA->m_invI;
	m_invIB = m_bodyB->m_invI;

	float32 aA = data.positions[m_indexA].a;
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;

	float32 aB = data.positions[m_indexB].a;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	b2Rot qA(aA), qB(aB);

	m_rA = b2Mul(qA, m_localAnchorA - m_localCenterA);
	m_rB = b2Mul(qB, m_localAnchorB - m_localCenterB);

	// J = [-I -r1_skew I r2_skew]
	//     [ 0       -1 0       1]
	// r_skew = [-ry; rx]

	// Matlab
	// K = [ mA+r1y^2*iA+mB+r2y^2*iB,  -r1y*iA*r1x-r2y*iB*r2x,          -r1y*iA-r2y*iB]
	//     [  -r1y*iA*r1x-r2y*iB*r2x, mA+r1x^2*iA+mB+r2x^2*iB,           r1x*iA+r2x*iB]
	//     [          -r1y*iA-r2y*iB,           r1x*iA+r2x*iB,                   iA+iB]

	float32 mA = m_invMassA, mB = m_invMassB;
	float32 iA = m_invIA, iB = m_invIB;

	b2Mat33 K;
	K.ex.x = mA + mB + m_rA.y * m_rA.y * iA + m_rB.y * m_rB.y * iB;
	K.ey.x = -m_rA.y * m_rA.x * iA - m_rB.y * m_rB.x * iB;
	K.ez.x = -m_rA.y * iA - m_rB.y * iB;
	K.ex.y = K.ey.x;
	K.ey.y = mA + mB + m_rA.x * m_rA.x * iA + m_rB.x * m_rB.x * iB;
	K.ez.y = m_rA.x * iA + m_rB.x * iB;
	K.ex.z = K.ez.x;
	K.ey.z = K.ez.y;
	K.ez.z = iA + iB;

	if (m_frequencyHz > 0.0f)
	{
		K.GetInverse22(&m_mass);

		float32 invM = iA + iB;
		float32 m = invM > 0.0f ? 1.0f / invM : 0.0f;

		float32 C = aB - aA - m_referenceAngle;

		// Frequency
		float32 omega = 2.0f * b2_pi * m_frequencyHz;

		// Damping coefficient
		float32 d = 2.0f * m * m_dampingRatio * omega;

		// Spring stiffness
		float32 k = m * omega * omega;

		// magic formulas
		float32 h = data.step.dt;
		m_gamma = h * (d + h * k);
		m_gamma = m_gamma != 0.0f ? 1.0f / m_gamma : 0.0f;
		m_bias = C * h * k * m_gamma;

		invM += m_gamma;
		m_mass.ez.z = invM != 0.0f ? 1.0f / invM : 0.0f;
	}
	else
	{
		K.GetSymInverse33(&m_mass);
		m_gamma = 0.0f;
		m_bias = 0.0f;
	}

	if (data.step.warmStarting)
	{
		// Scale impulses to support a variable time step.
		m_impulse *= data.step.dtRatio;

		b2Vec2 P(m_impulse.x, m_impulse.y);

		vA -= mA * P;
		wA -= iA * (b2Cross(m_rA, P) + m_impulse.z);

		vB += mB * P;
		wB += iB * (b2Cross(m_rB, P) + m_impulse.z);
	}
	else
	{
		m_impulse.SetZero();
	}

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

void b2WeldJoint::SolveVelocityConstraints(const b2SolverData& data)
{
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	float32 mA = m_invMassA, mB = m_invMassB;
	float32 iA = m_invIA, iB = m_invIB;

	if (m_frequencyHz > 0.0f)
	{
		float32 Cdot2 = wB - wA;

		float32 impulse2 = -m_mass.ez.z * (Cdot2 + m_bias + m_gamma * m_impulse.z);
		m_impulse.z += impulse2;

		wA -= iA * impulse2;
		wB += iB * impulse2;

		b2Vec2 Cdot1 = vB + b2Cross(wB, m_rB) - vA - b2Cross(wA, m_rA);

		b2Vec2 impulse1 = -b2Mul22(m_mass, Cdot1);
		m_impulse.x += impulse1.x;
		m_impulse.y += impulse1.y;

		b2Vec2 P = impulse1;

		vA -= mA * P;
		wA -= iA * b2Cross(m_rA, P);

		vB += mB * P;
		wB += iB * b2Cross(m_rB, P);
	}
	else
	{
		b2Vec2 Cdot1 = vB + b2Cross(wB, m_rB) - vA - b2Cross(wA, m_rA);
		float32 Cdot2 = wB - wA;
		b2Vec3 Cdot(Cdot1.x, Cdot1.y, Cdot2);

		b2Vec3 impulse = -b2Mul(m_mass, Cdot);
		m_impulse += impulse;

		b2Vec2 P(impulse.x, impulse.y);

		vA -= mA * P;
		wA -= iA * (b2Cross(m_rA, P) + impulse.z);

		vB += mB * P;
		wB += iB * (b2Cross(m_rB, P) + impulse.z);
	}

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

bool b2WeldJoint::SolvePositionConstraints(const b2SolverData& data)
{
	b2Vec2 cA = data.positions[m_indexA].c;
	float32 aA = data.positions[m_indexA].a;
	b2Vec2 cB = data.positions[m_indexB].c;
	float32 aB = data.positions[m_indexB].a;

	b2Rot qA(aA), qB(aB);

	float32 mA = m_invMassA, mB = m_invMassB;
	float32 iA = m_invIA, iB = m_invIB;

	b2Vec2 rA = b2Mul(qA, m_localAnchorA - m_localCenterA);
	b2Vec2 rB = b2Mul(qB, m_localAnchorB - m_localCenterB);

	float32 positionError, angularError;

	b2Mat33 K;
	K.ex.x = mA + mB + rA.y * rA.y * iA + rB.y * rB.y * iB;
	K.ey.x = -rA.y * rA.x * iA - rB.y * rB.x * iB;
	K.ez.x = -rA.y * iA - rB.y * iB;
	K.ex.y = K.ey.x;
	K.ey.y = mA + mB + rA.x * rA.x * iA + rB.x * rB.x * iB;
	K.ez.y = rA.x * iA + rB.x * iB;
	K.ex.z = K.ez.x;
	K.ey.z = K.ez.y;
	K.ez.z = iA + iB;

	if (m_frequencyHz > 0.0f)
	{
		b2Vec2 C1 =  cB + rB - cA - rA;

		positionError = C1.Length();
		angularError = 0.0f;

		b2Vec2 P = -K.Solve22(C1);

		cA -= mA * P;
		aA -= iA * b2Cross(rA, P);

		cB += mB * P;
		aB += iB * b2Cross(rB, P);
	}
	else
	{
		b2Vec2 C1 =  cB + rB - cA - rA;
		float32 C2 = aB - aA - m_referenceAngle;

		positionError = C1.Length();
		angularError = b2Abs(C2);

		b2Vec3 C(C1.x, C1.y, C2);
	
		b2Vec3 impulse = -K.Solve33(C);
		b2Vec2 P(impulse.x, impulse.y);

		cA -= mA * P;
		aA -= iA * (b2Cross(rA, P) + impulse.z);

		cB += mB * P;
		aB += iB * (b2Cross(rB, P) + impulse.z);
	}

	data.positions[m_indexA].c = cA;
	data.positions[m_indexA].a = aA;
	data.positions[m_indexB].c = cB;
	data.positions[m_indexB].a = aB;

	return positionError <= b2_linearSlop && angularError <= b2_angularSlop;
}

b2Vec2 b2WeldJoint::GetAnchorA() const
{
	return m_bodyA->GetWorldPoint(m_localAnchorA);
}

b2Vec2 b2WeldJoint::GetAnchorB() const
{
	return m_bodyB->GetWorldPoint(m_localAnchorB);
}

b2Vec2 b2WeldJoint::GetReactionForce(float32 inv_dt) const
{
	b2Vec2 P(m_impulse.x, m_impulse.y);
	return inv_dt * P;
}

float32 b2WeldJoint::GetReactionTorque(float32 inv_dt) const
{
	return inv_dt * m_impulse.z;
}

void b2WeldJoint::Dump()
{
	int32 indexA = m_bodyA->m_islandIndex;
	int32 indexB = m_bodyB->m_islandIndex;

	b2Log("  b2WeldJointDef jd;\n");
	b2Log("  jd.bodyA = bodies[%d];\n", indexA);
	b2Log("  jd.bodyB = bodies[%d];\n", indexB);
	b2Log("  jd.collideConnected = bool(%d);\n", m_collideConnected);
	b2Log("  jd.localAnchorA.Set(%.15lef, %.15lef);\n", m_localAnchorA.x, m_localAnchorA.y);
	b2Log("  jd.localAnchorB.Set(%.15lef, %.15lef);\n", m_localAnchorB.x, m_localAnchorB.y);
	b2Log("  jd.referenceAngle = %.15lef;\n", m_referenceAngle);
	b2Log("  jd.frequencyHz = %.15lef;\n", m_frequencyHz);
	b2Log("  jd.dampingRatio = %.15lef;\n", m_dampingRatio);
	b2Log("  joints[%d] = m_world->CreateJoint(&jd);\n", m_index);
}

// end of WeldJoint.cpp

/// Rope joint definition. This requires two body anchor points and
/// a maximum lengths.
/// Note: by default the connected objects will not collide.
/// see collideConnected in b2JointDef.
struct b2RopeJointDef : public b2JointDef
{
	b2RopeJointDef()
	{
		type = e_ropeJoint;
		localAnchorA.Set(-1.0f, 0.0f);
		localAnchorB.Set(1.0f, 0.0f);
		maxLength = 0.0f;
	}

	/// The local anchor point relative to bodyA's origin.
	b2Vec2 localAnchorA;

	/// The local anchor point relative to bodyB's origin.
	b2Vec2 localAnchorB;

	/// The maximum length of the rope.
	/// Warning: this must be larger than b2_linearSlop or
	/// the joint will have no effect.
	float32 maxLength;
};

/// A rope joint enforces a maximum distance between two points
/// on two bodies. It has no other effect.
/// Warning: if you attempt to change the maximum length during
/// the simulation you will get some non-physical behavior.
/// A model that would allow you to dynamically modify the length
/// would have some sponginess, so I chose not to implement it
/// that way. See b2DistanceJoint if you want to dynamically
/// control length.
class b2RopeJoint : public b2Joint
{
public:
	b2Vec2 GetAnchorA() const;
	b2Vec2 GetAnchorB() const;

	b2Vec2 GetReactionForce(float32 inv_dt) const;
	float32 GetReactionTorque(float32 inv_dt) const;

	/// The local anchor point relative to bodyA's origin.
	const b2Vec2& GetLocalAnchorA() const { return m_localAnchorA; }

	/// The local anchor point relative to bodyB's origin.
	const b2Vec2& GetLocalAnchorB() const  { return m_localAnchorB; }

	/// Set/Get the maximum length of the rope.
	void SetMaxLength(float32 length) { m_maxLength = length; }
	float32 GetMaxLength() const;

	b2LimitState GetLimitState() const;

	/// Dump joint to dmLog
	void Dump();

protected:

	friend class b2Joint;
	b2RopeJoint(const b2RopeJointDef* data);

	void InitVelocityConstraints(const b2SolverData& data);
	void SolveVelocityConstraints(const b2SolverData& data);
	bool SolvePositionConstraints(const b2SolverData& data);

	// Solver shared
	b2Vec2 m_localAnchorA;
	b2Vec2 m_localAnchorB;
	float32 m_maxLength;
	float32 m_length;
	float32 m_impulse;

	// Solver temp
	int32 m_indexA;
	int32 m_indexB;
	b2Vec2 m_u;
	b2Vec2 m_rA;
	b2Vec2 m_rB;
	b2Vec2 m_localCenterA;
	b2Vec2 m_localCenterB;
	float32 m_invMassA;
	float32 m_invMassB;
	float32 m_invIA;
	float32 m_invIB;
	float32 m_mass;
	b2LimitState m_state;
};

// end of RopeJoint.h

// Limit:
// C = norm(pB - pA) - L
// u = (pB - pA) / norm(pB - pA)
// Cdot = dot(u, vB + cross(wB, rB) - vA - cross(wA, rA))
// J = [-u -cross(rA, u) u cross(rB, u)]
// K = J * invM * JT
//   = invMassA + invIA * cross(rA, u)^2 + invMassB + invIB * cross(rB, u)^2

b2RopeJoint::b2RopeJoint(const b2RopeJointDef* def)
: b2Joint(def)
{
	m_localAnchorA = def->localAnchorA;
	m_localAnchorB = def->localAnchorB;

	m_maxLength = def->maxLength;

	m_mass = 0.0f;
	m_impulse = 0.0f;
	m_state = e_inactiveLimit;
	m_length = 0.0f;
}

void b2RopeJoint::InitVelocityConstraints(const b2SolverData& data)
{
	m_indexA = m_bodyA->m_islandIndex;
	m_indexB = m_bodyB->m_islandIndex;
	m_localCenterA = m_bodyA->m_sweep.localCenter;
	m_localCenterB = m_bodyB->m_sweep.localCenter;
	m_invMassA = m_bodyA->m_invMass;
	m_invMassB = m_bodyB->m_invMass;
	m_invIA = m_bodyA->m_invI;
	m_invIB = m_bodyB->m_invI;

	b2Vec2 cA = data.positions[m_indexA].c;
	float32 aA = data.positions[m_indexA].a;
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;

	b2Vec2 cB = data.positions[m_indexB].c;
	float32 aB = data.positions[m_indexB].a;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	b2Rot qA(aA), qB(aB);

	m_rA = b2Mul(qA, m_localAnchorA - m_localCenterA);
	m_rB = b2Mul(qB, m_localAnchorB - m_localCenterB);
	m_u = cB + m_rB - cA - m_rA;

	m_length = m_u.Length();

	float32 C = m_length - m_maxLength;
	if (C > 0.0f)
	{
		m_state = e_atUpperLimit;
	}
	else
	{
		m_state = e_inactiveLimit;
	}

	if (m_length > b2_linearSlop)
	{
		m_u *= 1.0f / m_length;
	}
	else
	{
		m_u.SetZero();
		m_mass = 0.0f;
		m_impulse = 0.0f;
		return;
	}

	// Compute effective mass.
	float32 crA = b2Cross(m_rA, m_u);
	float32 crB = b2Cross(m_rB, m_u);
	float32 invMass = m_invMassA + m_invIA * crA * crA + m_invMassB + m_invIB * crB * crB;

	m_mass = invMass != 0.0f ? 1.0f / invMass : 0.0f;

	if (data.step.warmStarting)
	{
		// Scale the impulse to support a variable time step.
		m_impulse *= data.step.dtRatio;

		b2Vec2 P = m_impulse * m_u;
		vA -= m_invMassA * P;
		wA -= m_invIA * b2Cross(m_rA, P);
		vB += m_invMassB * P;
		wB += m_invIB * b2Cross(m_rB, P);
	}
	else
	{
		m_impulse = 0.0f;
	}

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

void b2RopeJoint::SolveVelocityConstraints(const b2SolverData& data)
{
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	// Cdot = dot(u, v + cross(w, r))
	b2Vec2 vpA = vA + b2Cross(wA, m_rA);
	b2Vec2 vpB = vB + b2Cross(wB, m_rB);
	float32 C = m_length - m_maxLength;
	float32 Cdot = b2Dot(m_u, vpB - vpA);

	// Predictive constraint.
	if (C < 0.0f)
	{
		Cdot += data.step.inv_dt * C;
	}

	float32 impulse = -m_mass * Cdot;
	float32 oldImpulse = m_impulse;
	m_impulse = b2Min(0.0f, m_impulse + impulse);
	impulse = m_impulse - oldImpulse;

	b2Vec2 P = impulse * m_u;
	vA -= m_invMassA * P;
	wA -= m_invIA * b2Cross(m_rA, P);
	vB += m_invMassB * P;
	wB += m_invIB * b2Cross(m_rB, P);

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

bool b2RopeJoint::SolvePositionConstraints(const b2SolverData& data)
{
	b2Vec2 cA = data.positions[m_indexA].c;
	float32 aA = data.positions[m_indexA].a;
	b2Vec2 cB = data.positions[m_indexB].c;
	float32 aB = data.positions[m_indexB].a;

	b2Rot qA(aA), qB(aB);

	b2Vec2 rA = b2Mul(qA, m_localAnchorA - m_localCenterA);
	b2Vec2 rB = b2Mul(qB, m_localAnchorB - m_localCenterB);
	b2Vec2 u = cB + rB - cA - rA;

	float32 length = u.Normalize();
	float32 C = length - m_maxLength;

	C = b2Clamp(C, 0.0f, b2_maxLinearCorrection);

	float32 impulse = -m_mass * C;
	b2Vec2 P = impulse * u;

	cA -= m_invMassA * P;
	aA -= m_invIA * b2Cross(rA, P);
	cB += m_invMassB * P;
	aB += m_invIB * b2Cross(rB, P);

	data.positions[m_indexA].c = cA;
	data.positions[m_indexA].a = aA;
	data.positions[m_indexB].c = cB;
	data.positions[m_indexB].a = aB;

	return length - m_maxLength < b2_linearSlop;
}

b2Vec2 b2RopeJoint::GetAnchorA() const
{
	return m_bodyA->GetWorldPoint(m_localAnchorA);
}

b2Vec2 b2RopeJoint::GetAnchorB() const
{
	return m_bodyB->GetWorldPoint(m_localAnchorB);
}

b2Vec2 b2RopeJoint::GetReactionForce(float32 inv_dt) const
{
	b2Vec2 F = (inv_dt * m_impulse) * m_u;
	return F;
}

float32 b2RopeJoint::GetReactionTorque(float32 inv_dt) const
{
	B2_NOT_USED(inv_dt);
	return 0.0f;
}

float32 b2RopeJoint::GetMaxLength() const
{
	return m_maxLength;
}

b2LimitState b2RopeJoint::GetLimitState() const
{
	return m_state;
}

void b2RopeJoint::Dump()
{
	int32 indexA = m_bodyA->m_islandIndex;
	int32 indexB = m_bodyB->m_islandIndex;

	b2Log("  b2RopeJointDef jd;\n");
	b2Log("  jd.bodyA = bodies[%d];\n", indexA);
	b2Log("  jd.bodyB = bodies[%d];\n", indexB);
	b2Log("  jd.collideConnected = bool(%d);\n", m_collideConnected);
	b2Log("  jd.localAnchorA.Set(%.15lef, %.15lef);\n", m_localAnchorA.x, m_localAnchorA.y);
	b2Log("  jd.localAnchorB.Set(%.15lef, %.15lef);\n", m_localAnchorB.x, m_localAnchorB.y);
	b2Log("  jd.maxLength = %.15lef;\n", m_maxLength);
	b2Log("  joints[%d] = m_world->CreateJoint(&jd);\n", m_index);
}

// end of RopeJoint.cpp

/// Motor joint definition.
struct b2MotorJointDef : public b2JointDef
{
	b2MotorJointDef()
	{
		type = e_motorJoint;
		linearOffset.SetZero();
		angularOffset = 0.0f;
		maxForce = 1.0f;
		maxTorque = 1.0f;
		correctionFactor = 0.3f;
	}

	/// Initialize the bodies and offsets using the current transforms.
	void Initialize(b2Body* bodyA, b2Body* bodyB);

	/// Position of bodyB minus the position of bodyA, in bodyA's frame, in meters.
	b2Vec2 linearOffset;

	/// The bodyB angle minus bodyA angle in radians.
	float32 angularOffset;
	
	/// The maximum motor force in N.
	float32 maxForce;

	/// The maximum motor torque in N-m.
	float32 maxTorque;

	/// Position correction factor in the range [0,1].
	float32 correctionFactor;
};

/// A motor joint is used to control the relative motion
/// between two bodies. A typical usage is to control the movement
/// of a dynamic body with respect to the ground.
class b2MotorJoint : public b2Joint
{
public:
	b2Vec2 GetAnchorA() const;
	b2Vec2 GetAnchorB() const;

	b2Vec2 GetReactionForce(float32 inv_dt) const;
	float32 GetReactionTorque(float32 inv_dt) const;

	/// Set/get the target linear offset, in frame A, in meters.
	void SetLinearOffset(const b2Vec2& linearOffset);
	const b2Vec2& GetLinearOffset() const;

	/// Set/get the target angular offset, in radians.
	void SetAngularOffset(float32 angularOffset);
	float32 GetAngularOffset() const;

	/// Set the maximum friction force in N.
	void SetMaxForce(float32 force);

	/// Get the maximum friction force in N.
	float32 GetMaxForce() const;

	/// Set the maximum friction torque in N*m.
	void SetMaxTorque(float32 torque);

	/// Get the maximum friction torque in N*m.
	float32 GetMaxTorque() const;

	/// Set the position correction factor in the range [0,1].
	void SetCorrectionFactor(float32 factor);

	/// Get the position correction factor in the range [0,1].
	float32 GetCorrectionFactor() const;

	/// Dump to b2Log
	void Dump();

protected:

	friend class b2Joint;

	b2MotorJoint(const b2MotorJointDef* def);

	void InitVelocityConstraints(const b2SolverData& data);
	void SolveVelocityConstraints(const b2SolverData& data);
	bool SolvePositionConstraints(const b2SolverData& data);

	// Solver shared
	b2Vec2 m_linearOffset;
	float32 m_angularOffset;
	b2Vec2 m_linearImpulse;
	float32 m_angularImpulse;
	float32 m_maxForce;
	float32 m_maxTorque;
	float32 m_correctionFactor;

	// Solver temp
	int32 m_indexA;
	int32 m_indexB;
	b2Vec2 m_rA;
	b2Vec2 m_rB;
	b2Vec2 m_localCenterA;
	b2Vec2 m_localCenterB;
	b2Vec2 m_linearError;
	float32 m_angularError;
	float32 m_invMassA;
	float32 m_invMassB;
	float32 m_invIA;
	float32 m_invIB;
	b2Mat22 m_linearMass;
	float32 m_angularMass;
};

// end of MotorJoint.h

// Point-to-point constraint
// Cdot = v2 - v1
//      = v2 + cross(w2, r2) - v1 - cross(w1, r1)
// J = [-I -r1_skew I r2_skew ]
// Identity used:
// w k % (rx i + ry j) = w * (-ry i + rx j)

// Angle constraint
// Cdot = w2 - w1
// J = [0 0 -1 0 0 1]
// K = invI1 + invI2

void b2MotorJointDef::Initialize(b2Body* bA, b2Body* bB)
{
	bodyA = bA;
	bodyB = bB;
	b2Vec2 xB = bodyB->GetPosition();
	linearOffset = bodyA->GetLocalPoint(xB);

	float32 angleA = bodyA->GetAngle();
	float32 angleB = bodyB->GetAngle();
	angularOffset = angleB - angleA;
}

b2MotorJoint::b2MotorJoint(const b2MotorJointDef* def)
: b2Joint(def)
{
	m_linearOffset = def->linearOffset;
	m_angularOffset = def->angularOffset;

	m_linearImpulse.SetZero();
	m_angularImpulse = 0.0f;

	m_maxForce = def->maxForce;
	m_maxTorque = def->maxTorque;
	m_correctionFactor = def->correctionFactor;
}

void b2MotorJoint::InitVelocityConstraints(const b2SolverData& data)
{
	m_indexA = m_bodyA->m_islandIndex;
	m_indexB = m_bodyB->m_islandIndex;
	m_localCenterA = m_bodyA->m_sweep.localCenter;
	m_localCenterB = m_bodyB->m_sweep.localCenter;
	m_invMassA = m_bodyA->m_invMass;
	m_invMassB = m_bodyB->m_invMass;
	m_invIA = m_bodyA->m_invI;
	m_invIB = m_bodyB->m_invI;

	b2Vec2 cA = data.positions[m_indexA].c;
	float32 aA = data.positions[m_indexA].a;
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;

	b2Vec2 cB = data.positions[m_indexB].c;
	float32 aB = data.positions[m_indexB].a;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	b2Rot qA(aA), qB(aB);

	// Compute the effective mass matrix.
	m_rA = b2Mul(qA, -m_localCenterA);
	m_rB = b2Mul(qB, -m_localCenterB);

	// J = [-I -r1_skew I r2_skew]
	//     [ 0       -1 0       1]
	// r_skew = [-ry; rx]

	// Matlab
	// K = [ mA+r1y^2*iA+mB+r2y^2*iB,  -r1y*iA*r1x-r2y*iB*r2x,          -r1y*iA-r2y*iB]
	//     [  -r1y*iA*r1x-r2y*iB*r2x, mA+r1x^2*iA+mB+r2x^2*iB,           r1x*iA+r2x*iB]
	//     [          -r1y*iA-r2y*iB,           r1x*iA+r2x*iB,                   iA+iB]

	float32 mA = m_invMassA, mB = m_invMassB;
	float32 iA = m_invIA, iB = m_invIB;

	b2Mat22 K;
	K.ex.x = mA + mB + iA * m_rA.y * m_rA.y + iB * m_rB.y * m_rB.y;
	K.ex.y = -iA * m_rA.x * m_rA.y - iB * m_rB.x * m_rB.y;
	K.ey.x = K.ex.y;
	K.ey.y = mA + mB + iA * m_rA.x * m_rA.x + iB * m_rB.x * m_rB.x;

	m_linearMass = K.GetInverse();

	m_angularMass = iA + iB;
	if (m_angularMass > 0.0f)
	{
		m_angularMass = 1.0f / m_angularMass;
	}

	m_linearError = cB + m_rB - cA - m_rA - b2Mul(qA, m_linearOffset);
	m_angularError = aB - aA - m_angularOffset;

	if (data.step.warmStarting)
	{
		// Scale impulses to support a variable time step.
		m_linearImpulse *= data.step.dtRatio;
		m_angularImpulse *= data.step.dtRatio;

		b2Vec2 P(m_linearImpulse.x, m_linearImpulse.y);
		vA -= mA * P;
		wA -= iA * (b2Cross(m_rA, P) + m_angularImpulse);
		vB += mB * P;
		wB += iB * (b2Cross(m_rB, P) + m_angularImpulse);
	}
	else
	{
		m_linearImpulse.SetZero();
		m_angularImpulse = 0.0f;
	}

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

void b2MotorJoint::SolveVelocityConstraints(const b2SolverData& data)
{
	b2Vec2 vA = data.velocities[m_indexA].v;
	float32 wA = data.velocities[m_indexA].w;
	b2Vec2 vB = data.velocities[m_indexB].v;
	float32 wB = data.velocities[m_indexB].w;

	float32 mA = m_invMassA, mB = m_invMassB;
	float32 iA = m_invIA, iB = m_invIB;

	float32 h = data.step.dt;
	float32 inv_h = data.step.inv_dt;

	// Solve angular friction
	{
		float32 Cdot = wB - wA + inv_h * m_correctionFactor * m_angularError;
		float32 impulse = -m_angularMass * Cdot;

		float32 oldImpulse = m_angularImpulse;
		float32 maxImpulse = h * m_maxTorque;
		m_angularImpulse = b2Clamp(m_angularImpulse + impulse, -maxImpulse, maxImpulse);
		impulse = m_angularImpulse - oldImpulse;

		wA -= iA * impulse;
		wB += iB * impulse;
	}

	// Solve linear friction
	{
		b2Vec2 Cdot = vB + b2Cross(wB, m_rB) - vA - b2Cross(wA, m_rA) + inv_h * m_correctionFactor * m_linearError;

		b2Vec2 impulse = -b2Mul(m_linearMass, Cdot);
		b2Vec2 oldImpulse = m_linearImpulse;
		m_linearImpulse += impulse;

		float32 maxImpulse = h * m_maxForce;

		if (m_linearImpulse.LengthSquared() > maxImpulse * maxImpulse)
		{
			m_linearImpulse.Normalize();
			m_linearImpulse *= maxImpulse;
		}

		impulse = m_linearImpulse - oldImpulse;

		vA -= mA * impulse;
		wA -= iA * b2Cross(m_rA, impulse);

		vB += mB * impulse;
		wB += iB * b2Cross(m_rB, impulse);
	}

	data.velocities[m_indexA].v = vA;
	data.velocities[m_indexA].w = wA;
	data.velocities[m_indexB].v = vB;
	data.velocities[m_indexB].w = wB;
}

bool b2MotorJoint::SolvePositionConstraints(const b2SolverData& data)
{
	B2_NOT_USED(data);

	return true;
}

b2Vec2 b2MotorJoint::GetAnchorA() const
{
	return m_bodyA->GetPosition();
}

b2Vec2 b2MotorJoint::GetAnchorB() const
{
	return m_bodyB->GetPosition();
}

b2Vec2 b2MotorJoint::GetReactionForce(float32 inv_dt) const
{
	return inv_dt * m_linearImpulse;
}

float32 b2MotorJoint::GetReactionTorque(float32 inv_dt) const
{
	return inv_dt * m_angularImpulse;
}

void b2MotorJoint::SetMaxForce(float32 force)
{
	b2Assert(b2IsValid(force) && force >= 0.0f);
	m_maxForce = force;
}

float32 b2MotorJoint::GetMaxForce() const
{
	return m_maxForce;
}

void b2MotorJoint::SetMaxTorque(float32 torque)
{
	b2Assert(b2IsValid(torque) && torque >= 0.0f);
	m_maxTorque = torque;
}

float32 b2MotorJoint::GetMaxTorque() const
{
	return m_maxTorque;
}

void b2MotorJoint::SetCorrectionFactor(float32 factor)
{
	b2Assert(b2IsValid(factor) && 0.0f <= factor && factor <= 1.0f);
	m_correctionFactor = factor;
}

float32 b2MotorJoint::GetCorrectionFactor() const
{
	return m_correctionFactor;
}

void b2MotorJoint::SetLinearOffset(const b2Vec2& linearOffset)
{
	if (linearOffset.x != m_linearOffset.x || linearOffset.y != m_linearOffset.y)
	{
		m_bodyA->SetAwake(true);
		m_bodyB->SetAwake(true);
		m_linearOffset = linearOffset;
	}
}

const b2Vec2& b2MotorJoint::GetLinearOffset() const
{
	return m_linearOffset;
}

void b2MotorJoint::SetAngularOffset(float32 angularOffset)
{
	if (angularOffset != m_angularOffset)
	{
		m_bodyA->SetAwake(true);
		m_bodyB->SetAwake(true);
		m_angularOffset = angularOffset;
	}
}

float32 b2MotorJoint::GetAngularOffset() const
{
	return m_angularOffset;
}

void b2MotorJoint::Dump()
{
	int32 indexA = m_bodyA->m_islandIndex;
	int32 indexB = m_bodyB->m_islandIndex;

	b2Log("  b2MotorJointDef jd;\n");
	b2Log("  jd.bodyA = bodies[%d];\n", indexA);
	b2Log("  jd.bodyB = bodies[%d];\n", indexB);
	b2Log("  jd.collideConnected = bool(%d);\n", m_collideConnected);
	b2Log("  jd.linearOffset.Set(%.15lef, %.15lef);\n", m_linearOffset.x, m_linearOffset.y);
	b2Log("  jd.angularOffset = %.15lef;\n", m_angularOffset);
	b2Log("  jd.maxForce = %.15lef;\n", m_maxForce);
	b2Log("  jd.maxTorque = %.15lef;\n", m_maxTorque);
	b2Log("  jd.correctionFactor = %.15lef;\n", m_correctionFactor);
	b2Log("  joints[%d] = m_world->CreateJoint(&jd);\n", m_index);
}

// end of MotorJoint.cpp

b2Joint* b2Joint::Create(const b2JointDef* def, b2BlockAllocator* allocator)
{
	b2Joint* joint = NULL;

	switch (def->type)
	{
	case e_distanceJoint:
		{
			void* mem = allocator->Allocate(sizeof(b2DistanceJoint));
			joint = new (mem) b2DistanceJoint(static_cast<const b2DistanceJointDef*>(def));
		}
		break;

	case e_mouseJoint:
		{
			void* mem = allocator->Allocate(sizeof(b2MouseJoint));
			joint = new (mem) b2MouseJoint(static_cast<const b2MouseJointDef*>(def));
		}
		break;

	case e_prismaticJoint:
		{
			void* mem = allocator->Allocate(sizeof(b2PrismaticJoint));
			joint = new (mem) b2PrismaticJoint(static_cast<const b2PrismaticJointDef*>(def));
		}
		break;

	case e_revoluteJoint:
		{
			void* mem = allocator->Allocate(sizeof(b2RevoluteJoint));
			joint = new (mem) b2RevoluteJoint(static_cast<const b2RevoluteJointDef*>(def));
		}
		break;

	case e_pulleyJoint:
		{
			void* mem = allocator->Allocate(sizeof(b2PulleyJoint));
			joint = new (mem) b2PulleyJoint(static_cast<const b2PulleyJointDef*>(def));
		}
		break;

	case e_gearJoint:
		{
			void* mem = allocator->Allocate(sizeof(b2GearJoint));
			joint = new (mem) b2GearJoint(static_cast<const b2GearJointDef*>(def));
		}
		break;

	case e_wheelJoint:
		{
			void* mem = allocator->Allocate(sizeof(b2WheelJoint));
			joint = new (mem) b2WheelJoint(static_cast<const b2WheelJointDef*>(def));
		}
		break;

	case e_weldJoint:
		{
			void* mem = allocator->Allocate(sizeof(b2WeldJoint));
			joint = new (mem) b2WeldJoint(static_cast<const b2WeldJointDef*>(def));
		}
		break;
        
	case e_frictionJoint:
		{
			void* mem = allocator->Allocate(sizeof(b2FrictionJoint));
			joint = new (mem) b2FrictionJoint(static_cast<const b2FrictionJointDef*>(def));
		}
		break;

	case e_ropeJoint:
		{
			void* mem = allocator->Allocate(sizeof(b2RopeJoint));
			joint = new (mem) b2RopeJoint(static_cast<const b2RopeJointDef*>(def));
		}
		break;

	case e_motorJoint:
		{
			void* mem = allocator->Allocate(sizeof(b2MotorJoint));
			joint = new (mem) b2MotorJoint(static_cast<const b2MotorJointDef*>(def));
		}
		break;

	default:
		b2Assert(false);
		break;
	}

	return joint;
}

void b2Joint::Destroy(b2Joint* joint, b2BlockAllocator* allocator)
{
	joint->~b2Joint();
	switch (joint->m_type)
	{
	case e_distanceJoint:
		allocator->Free(joint, sizeof(b2DistanceJoint));
		break;

	case e_mouseJoint:
		allocator->Free(joint, sizeof(b2MouseJoint));
		break;

	case e_prismaticJoint:
		allocator->Free(joint, sizeof(b2PrismaticJoint));
		break;

	case e_revoluteJoint:
		allocator->Free(joint, sizeof(b2RevoluteJoint));
		break;

	case e_pulleyJoint:
		allocator->Free(joint, sizeof(b2PulleyJoint));
		break;

	case e_gearJoint:
		allocator->Free(joint, sizeof(b2GearJoint));
		break;

	case e_wheelJoint:
		allocator->Free(joint, sizeof(b2WheelJoint));
		break;
    
	case e_weldJoint:
		allocator->Free(joint, sizeof(b2WeldJoint));
		break;

	case e_frictionJoint:
		allocator->Free(joint, sizeof(b2FrictionJoint));
		break;

	case e_ropeJoint:
		allocator->Free(joint, sizeof(b2RopeJoint));
		break;

	case e_motorJoint:
		allocator->Free(joint, sizeof(b2MotorJoint));
		break;

	default:
		b2Assert(false);
		break;
	}
}

b2Joint::b2Joint(const b2JointDef* def)
{
	b2Assert(def->bodyA != def->bodyB);

	m_type = def->type;
	m_prev = NULL;
	m_next = NULL;
	m_bodyA = def->bodyA;
	m_bodyB = def->bodyB;
	m_index = 0;
	m_collideConnected = def->collideConnected;
	m_islandFlag = false;
	m_userData = def->userData;

	m_edgeA.joint = NULL;
	m_edgeA.other = NULL;
	m_edgeA.prev = NULL;
	m_edgeA.next = NULL;

	m_edgeB.joint = NULL;
	m_edgeB.other = NULL;
	m_edgeB.prev = NULL;
	m_edgeB.next = NULL;
}

bool b2Joint::IsActive() const
{
	return m_bodyA->IsActive() && m_bodyB->IsActive();
}

// end of Joint.cpp

class b2BlockAllocator;

class b2ChainAndCircleContact : public b2Contact
{
public:
	static b2Contact* Create(	b2Fixture* fixtureA, int32 indexA,
								b2Fixture* fixtureB, int32 indexB, b2BlockAllocator* allocator);
	static void Destroy(b2Contact* contact, b2BlockAllocator* allocator);

	b2ChainAndCircleContact(b2Fixture* fixtureA, int32 indexA, b2Fixture* fixtureB, int32 indexB);
	~b2ChainAndCircleContact() {}

	void Evaluate(b2Manifold* manifold, const b2Transform& xfA, const b2Transform& xfB);
};

// end of ChainAndCircleContact.h

b2Contact* b2ChainAndCircleContact::Create(b2Fixture* fixtureA, int32 indexA, b2Fixture* fixtureB, int32 indexB, b2BlockAllocator* allocator)
{
	void* mem = allocator->Allocate(sizeof(b2ChainAndCircleContact));
	return new (mem) b2ChainAndCircleContact(fixtureA, indexA, fixtureB, indexB);
}

void b2ChainAndCircleContact::Destroy(b2Contact* contact, b2BlockAllocator* allocator)
{
	((b2ChainAndCircleContact*)contact)->~b2ChainAndCircleContact();
	allocator->Free(contact, sizeof(b2ChainAndCircleContact));
}

b2ChainAndCircleContact::b2ChainAndCircleContact(b2Fixture* fixtureA, int32 indexA, b2Fixture* fixtureB, int32 indexB)
: b2Contact(fixtureA, indexA, fixtureB, indexB)
{
	b2Assert(m_fixtureA->GetType() == b2Shape::e_chain);
	b2Assert(m_fixtureB->GetType() == b2Shape::e_circle);
}

void b2ChainAndCircleContact::Evaluate(b2Manifold* manifold, const b2Transform& xfA, const b2Transform& xfB)
{
	b2ChainShape* chain = (b2ChainShape*)m_fixtureA->GetShape();
	b2EdgeShape edge;
	chain->GetChildEdge(&edge, m_indexA);
	b2CollideEdgeAndCircle(	manifold, &edge, xfA,
							(b2CircleShape*)m_fixtureB->GetShape(), xfB);
}

// end of ChainAndCircleContact.cpp

class b2BlockAllocator;

class b2ChainAndPolygonContact : public b2Contact
{
public:
	static b2Contact* Create(	b2Fixture* fixtureA, int32 indexA,
								b2Fixture* fixtureB, int32 indexB, b2BlockAllocator* allocator);
	static void Destroy(b2Contact* contact, b2BlockAllocator* allocator);

	b2ChainAndPolygonContact(b2Fixture* fixtureA, int32 indexA, b2Fixture* fixtureB, int32 indexB);
	~b2ChainAndPolygonContact() {}

	void Evaluate(b2Manifold* manifold, const b2Transform& xfA, const b2Transform& xfB);
};

// end of ChainAndPolygonContact.h

b2Contact* b2ChainAndPolygonContact::Create(b2Fixture* fixtureA, int32 indexA, b2Fixture* fixtureB, int32 indexB, b2BlockAllocator* allocator)
{
	void* mem = allocator->Allocate(sizeof(b2ChainAndPolygonContact));
	return new (mem) b2ChainAndPolygonContact(fixtureA, indexA, fixtureB, indexB);
}

void b2ChainAndPolygonContact::Destroy(b2Contact* contact, b2BlockAllocator* allocator)
{
	((b2ChainAndPolygonContact*)contact)->~b2ChainAndPolygonContact();
	allocator->Free(contact, sizeof(b2ChainAndPolygonContact));
}

b2ChainAndPolygonContact::b2ChainAndPolygonContact(b2Fixture* fixtureA, int32 indexA, b2Fixture* fixtureB, int32 indexB)
: b2Contact(fixtureA, indexA, fixtureB, indexB)
{
	b2Assert(m_fixtureA->GetType() == b2Shape::e_chain);
	b2Assert(m_fixtureB->GetType() == b2Shape::e_polygon);
}

void b2ChainAndPolygonContact::Evaluate(b2Manifold* manifold, const b2Transform& xfA, const b2Transform& xfB)
{
	b2ChainShape* chain = (b2ChainShape*)m_fixtureA->GetShape();
	b2EdgeShape edge;
	chain->GetChildEdge(&edge, m_indexA);
	b2CollideEdgeAndPolygon(	manifold, &edge, xfA,
								(b2PolygonShape*)m_fixtureB->GetShape(), xfB);
}

// end of ChainAndPolygonContact.cpp

class b2BlockAllocator;

class b2CircleContact : public b2Contact
{
public:
	static b2Contact* Create(	b2Fixture* fixtureA, int32 indexA,
								b2Fixture* fixtureB, int32 indexB, b2BlockAllocator* allocator);
	static void Destroy(b2Contact* contact, b2BlockAllocator* allocator);

	b2CircleContact(b2Fixture* fixtureA, b2Fixture* fixtureB);
	~b2CircleContact() {}

	void Evaluate(b2Manifold* manifold, const b2Transform& xfA, const b2Transform& xfB);
};

// end of CircleContact.h

b2Contact* b2CircleContact::Create(b2Fixture* fixtureA, int32, b2Fixture* fixtureB, int32, b2BlockAllocator* allocator)
{
	void* mem = allocator->Allocate(sizeof(b2CircleContact));
	return new (mem) b2CircleContact(fixtureA, fixtureB);
}

void b2CircleContact::Destroy(b2Contact* contact, b2BlockAllocator* allocator)
{
	((b2CircleContact*)contact)->~b2CircleContact();
	allocator->Free(contact, sizeof(b2CircleContact));
}

b2CircleContact::b2CircleContact(b2Fixture* fixtureA, b2Fixture* fixtureB)
	: b2Contact(fixtureA, 0, fixtureB, 0)
{
	b2Assert(m_fixtureA->GetType() == b2Shape::e_circle);
	b2Assert(m_fixtureB->GetType() == b2Shape::e_circle);
}

void b2CircleContact::Evaluate(b2Manifold* manifold, const b2Transform& xfA, const b2Transform& xfB)
{
	b2CollideCircles(manifold,
					(b2CircleShape*)m_fixtureA->GetShape(), xfA,
					(b2CircleShape*)m_fixtureB->GetShape(), xfB);
}

// end of CircleContact.cpp

#define B2_DEBUG_SOLVER 0

struct b2ContactPositionConstraint
{
	b2Vec2 localPoints[b2_maxManifoldPoints];
	b2Vec2 localNormal;
	b2Vec2 localPoint;
	int32 indexA;
	int32 indexB;
	float32 invMassA, invMassB;
	b2Vec2 localCenterA, localCenterB;
	float32 invIA, invIB;
	b2Manifold::Type type;
	float32 radiusA, radiusB;
	int32 pointCount;
};

b2ContactSolver::b2ContactSolver(b2ContactSolverDef* def)
{
	m_step = def->step;
	m_allocator = def->allocator;
	m_count = def->count;
	m_positionConstraints = (b2ContactPositionConstraint*)m_allocator->Allocate(m_count * sizeof(b2ContactPositionConstraint));
	m_velocityConstraints = (b2ContactVelocityConstraint*)m_allocator->Allocate(m_count * sizeof(b2ContactVelocityConstraint));
	m_positions = def->positions;
	m_velocities = def->velocities;
	m_contacts = def->contacts;

	// Initialize position independent portions of the constraints.
	for (int32 i = 0; i < m_count; ++i)
	{
		b2Contact* contact = m_contacts[i];

		b2Fixture* fixtureA = contact->m_fixtureA;
		b2Fixture* fixtureB = contact->m_fixtureB;
		b2Shape* shapeA = fixtureA->GetShape();
		b2Shape* shapeB = fixtureB->GetShape();
		float32 radiusA = shapeA->m_radius;
		float32 radiusB = shapeB->m_radius;
		b2Body* bodyA = fixtureA->GetBody();
		b2Body* bodyB = fixtureB->GetBody();
		b2Manifold* manifold = contact->GetManifold();

		int32 pointCount = manifold->pointCount;
		b2Assert(pointCount > 0);

		b2ContactVelocityConstraint* vc = m_velocityConstraints + i;
		vc->friction = contact->m_friction;
		vc->restitution = contact->m_restitution;
		vc->tangentSpeed = contact->m_tangentSpeed;
		vc->indexA = bodyA->m_islandIndex;
		vc->indexB = bodyB->m_islandIndex;
		vc->invMassA = bodyA->m_invMass;
		vc->invMassB = bodyB->m_invMass;
		vc->invIA = bodyA->m_invI;
		vc->invIB = bodyB->m_invI;
		vc->contactIndex = i;
		vc->pointCount = pointCount;
		vc->K.SetZero();
		vc->normalMass.SetZero();

		b2ContactPositionConstraint* pc = m_positionConstraints + i;
		pc->indexA = bodyA->m_islandIndex;
		pc->indexB = bodyB->m_islandIndex;
		pc->invMassA = bodyA->m_invMass;
		pc->invMassB = bodyB->m_invMass;
		pc->localCenterA = bodyA->m_sweep.localCenter;
		pc->localCenterB = bodyB->m_sweep.localCenter;
		pc->invIA = bodyA->m_invI;
		pc->invIB = bodyB->m_invI;
		pc->localNormal = manifold->localNormal;
		pc->localPoint = manifold->localPoint;
		pc->pointCount = pointCount;
		pc->radiusA = radiusA;
		pc->radiusB = radiusB;
		pc->type = manifold->type;

		for (int32 j = 0; j < pointCount; ++j)
		{
			b2ManifoldPoint* cp = manifold->points + j;
			b2VelocityConstraintPoint* vcp = vc->points + j;
	
			if (m_step.warmStarting)
			{
				vcp->normalImpulse = m_step.dtRatio * cp->normalImpulse;
				vcp->tangentImpulse = m_step.dtRatio * cp->tangentImpulse;
			}
			else
			{
				vcp->normalImpulse = 0.0f;
				vcp->tangentImpulse = 0.0f;
			}

			vcp->rA.SetZero();
			vcp->rB.SetZero();
			vcp->normalMass = 0.0f;
			vcp->tangentMass = 0.0f;
			vcp->velocityBias = 0.0f;

			pc->localPoints[j] = cp->localPoint;
		}
	}
}

b2ContactSolver::~b2ContactSolver()
{
	m_allocator->Free(m_velocityConstraints);
	m_allocator->Free(m_positionConstraints);
}

// Initialize position dependent portions of the velocity constraints.
void b2ContactSolver::InitializeVelocityConstraints()
{
	for (int32 i = 0; i < m_count; ++i)
	{
		b2ContactVelocityConstraint* vc = m_velocityConstraints + i;
		b2ContactPositionConstraint* pc = m_positionConstraints + i;

		float32 radiusA = pc->radiusA;
		float32 radiusB = pc->radiusB;
		b2Manifold* manifold = m_contacts[vc->contactIndex]->GetManifold();

		int32 indexA = vc->indexA;
		int32 indexB = vc->indexB;

		float32 mA = vc->invMassA;
		float32 mB = vc->invMassB;
		float32 iA = vc->invIA;
		float32 iB = vc->invIB;
		b2Vec2 localCenterA = pc->localCenterA;
		b2Vec2 localCenterB = pc->localCenterB;

		b2Vec2 cA = m_positions[indexA].c;
		float32 aA = m_positions[indexA].a;
		b2Vec2 vA = m_velocities[indexA].v;
		float32 wA = m_velocities[indexA].w;

		b2Vec2 cB = m_positions[indexB].c;
		float32 aB = m_positions[indexB].a;
		b2Vec2 vB = m_velocities[indexB].v;
		float32 wB = m_velocities[indexB].w;

		b2Assert(manifold->pointCount > 0);

		b2Transform xfA, xfB;
		xfA.q.Set(aA);
		xfB.q.Set(aB);
		xfA.p = cA - b2Mul(xfA.q, localCenterA);
		xfB.p = cB - b2Mul(xfB.q, localCenterB);

		b2WorldManifold worldManifold;
		worldManifold.Initialize(manifold, xfA, radiusA, xfB, radiusB);

		vc->normal = worldManifold.normal;

		int32 pointCount = vc->pointCount;
		for (int32 j = 0; j < pointCount; ++j)
		{
			b2VelocityConstraintPoint* vcp = vc->points + j;

			vcp->rA = worldManifold.points[j] - cA;
			vcp->rB = worldManifold.points[j] - cB;

			float32 rnA = b2Cross(vcp->rA, vc->normal);
			float32 rnB = b2Cross(vcp->rB, vc->normal);

			float32 kNormal = mA + mB + iA * rnA * rnA + iB * rnB * rnB;

			vcp->normalMass = kNormal > 0.0f ? 1.0f / kNormal : 0.0f;

			b2Vec2 tangent = b2Cross(vc->normal, 1.0f);

			float32 rtA = b2Cross(vcp->rA, tangent);
			float32 rtB = b2Cross(vcp->rB, tangent);

			float32 kTangent = mA + mB + iA * rtA * rtA + iB * rtB * rtB;

			vcp->tangentMass = kTangent > 0.0f ? 1.0f /  kTangent : 0.0f;

			// Setup a velocity bias for restitution.
			vcp->velocityBias = 0.0f;
			float32 vRel = b2Dot(vc->normal, vB + b2Cross(wB, vcp->rB) - vA - b2Cross(wA, vcp->rA));
			if (vRel < -b2_velocityThreshold)
			{
				vcp->velocityBias = -vc->restitution * vRel;
			}
		}

		// If we have two points, then prepare the block solver.
		if (vc->pointCount == 2)
		{
			b2VelocityConstraintPoint* vcp1 = vc->points + 0;
			b2VelocityConstraintPoint* vcp2 = vc->points + 1;

			float32 rn1A = b2Cross(vcp1->rA, vc->normal);
			float32 rn1B = b2Cross(vcp1->rB, vc->normal);
			float32 rn2A = b2Cross(vcp2->rA, vc->normal);
			float32 rn2B = b2Cross(vcp2->rB, vc->normal);

			float32 k11 = mA + mB + iA * rn1A * rn1A + iB * rn1B * rn1B;
			float32 k22 = mA + mB + iA * rn2A * rn2A + iB * rn2B * rn2B;
			float32 k12 = mA + mB + iA * rn1A * rn2A + iB * rn1B * rn2B;

			// Ensure a reasonable condition number.
			const float32 k_maxConditionNumber = 1000.0f;
			if (k11 * k11 < k_maxConditionNumber * (k11 * k22 - k12 * k12))
			{
				// K is safe to invert.
				vc->K.ex.Set(k11, k12);
				vc->K.ey.Set(k12, k22);
				vc->normalMass = vc->K.GetInverse();
			}
			else
			{
				// The constraints are redundant, just use one.
				// TODO_ERIN use deepest?
				vc->pointCount = 1;
			}
		}
	}
}

void b2ContactSolver::WarmStart()
{
	// Warm start.
	for (int32 i = 0; i < m_count; ++i)
	{
		b2ContactVelocityConstraint* vc = m_velocityConstraints + i;

		int32 indexA = vc->indexA;
		int32 indexB = vc->indexB;
		float32 mA = vc->invMassA;
		float32 iA = vc->invIA;
		float32 mB = vc->invMassB;
		float32 iB = vc->invIB;
		int32 pointCount = vc->pointCount;

		b2Vec2 vA = m_velocities[indexA].v;
		float32 wA = m_velocities[indexA].w;
		b2Vec2 vB = m_velocities[indexB].v;
		float32 wB = m_velocities[indexB].w;

		b2Vec2 normal = vc->normal;
		b2Vec2 tangent = b2Cross(normal, 1.0f);

		for (int32 j = 0; j < pointCount; ++j)
		{
			b2VelocityConstraintPoint* vcp = vc->points + j;
			b2Vec2 P = vcp->normalImpulse * normal + vcp->tangentImpulse * tangent;
			wA -= iA * b2Cross(vcp->rA, P);
			vA -= mA * P;
			wB += iB * b2Cross(vcp->rB, P);
			vB += mB * P;
		}

		m_velocities[indexA].v = vA;
		m_velocities[indexA].w = wA;
		m_velocities[indexB].v = vB;
		m_velocities[indexB].w = wB;
	}
}

void b2ContactSolver::SolveVelocityConstraints()
{
	for (int32 i = 0; i < m_count; ++i)
	{
		b2ContactVelocityConstraint* vc = m_velocityConstraints + i;

		int32 indexA = vc->indexA;
		int32 indexB = vc->indexB;
		float32 mA = vc->invMassA;
		float32 iA = vc->invIA;
		float32 mB = vc->invMassB;
		float32 iB = vc->invIB;
		int32 pointCount = vc->pointCount;

		b2Vec2 vA = m_velocities[indexA].v;
		float32 wA = m_velocities[indexA].w;
		b2Vec2 vB = m_velocities[indexB].v;
		float32 wB = m_velocities[indexB].w;

		b2Vec2 normal = vc->normal;
		b2Vec2 tangent = b2Cross(normal, 1.0f);
		float32 friction = vc->friction;

		b2Assert(pointCount == 1 || pointCount == 2);

		// Solve tangent constraints first because non-penetration is more important
		// than friction.
		for (int32 j = 0; j < pointCount; ++j)
		{
			b2VelocityConstraintPoint* vcp = vc->points + j;

			// Relative velocity at contact
			b2Vec2 dv = vB + b2Cross(wB, vcp->rB) - vA - b2Cross(wA, vcp->rA);

			// Compute tangent force
			float32 vt = b2Dot(dv, tangent) - vc->tangentSpeed;
			float32 lambda = vcp->tangentMass * (-vt);

			// b2Clamp the accumulated force
			float32 maxFriction = friction * vcp->normalImpulse;
			float32 newImpulse = b2Clamp(vcp->tangentImpulse + lambda, -maxFriction, maxFriction);
			lambda = newImpulse - vcp->tangentImpulse;
			vcp->tangentImpulse = newImpulse;

			// Apply contact impulse
			b2Vec2 P = lambda * tangent;

			vA -= mA * P;
			wA -= iA * b2Cross(vcp->rA, P);

			vB += mB * P;
			wB += iB * b2Cross(vcp->rB, P);
		}

		// Solve normal constraints
		if (vc->pointCount == 1)
		{
			b2VelocityConstraintPoint* vcp = vc->points + 0;

			// Relative velocity at contact
			b2Vec2 dv = vB + b2Cross(wB, vcp->rB) - vA - b2Cross(wA, vcp->rA);

			// Compute normal impulse
			float32 vn = b2Dot(dv, normal);
			float32 lambda = -vcp->normalMass * (vn - vcp->velocityBias);

			// b2Clamp the accumulated impulse
			float32 newImpulse = b2Max(vcp->normalImpulse + lambda, 0.0f);
			lambda = newImpulse - vcp->normalImpulse;
			vcp->normalImpulse = newImpulse;

			// Apply contact impulse
			b2Vec2 P = lambda * normal;
			vA -= mA * P;
			wA -= iA * b2Cross(vcp->rA, P);

			vB += mB * P;
			wB += iB * b2Cross(vcp->rB, P);
		}
		else
		{
			// Block solver developed in collaboration with Dirk Gregorius (back in 01/07 on Box2D_Lite).
			// Build the mini LCP for this contact patch
			//
			// vn = A * x + b, vn >= 0, , vn >= 0, x >= 0 and vn_i * x_i = 0 with i = 1..2
			//
			// A = J * W * JT and J = ( -n, -r1 x n, n, r2 x n )
			// b = vn0 - velocityBias
			//
			// The system is solved using the "Total enumeration method" (s. Murty). The complementary constraint vn_i * x_i
			// implies that we must have in any solution either vn_i = 0 or x_i = 0. So for the 2D contact problem the cases
			// vn1 = 0 and vn2 = 0, x1 = 0 and x2 = 0, x1 = 0 and vn2 = 0, x2 = 0 and vn1 = 0 need to be tested. The first valid
			// solution that satisfies the problem is chosen.
			// 
			// In order to account of the accumulated impulse 'a' (because of the iterative nature of the solver which only requires
			// that the accumulated impulse is clamped and not the incremental impulse) we change the impulse variable (x_i).
			//
			// Substitute:
			// 
			// x = a + d
			// 
			// a := old total impulse
			// x := new total impulse
			// d := incremental impulse 
			//
			// For the current iteration we extend the formula for the incremental impulse
			// to compute the new total impulse:
			//
			// vn = A * d + b
			//    = A * (x - a) + b
			//    = A * x + b - A * a
			//    = A * x + b'
			// b' = b - A * a;

			b2VelocityConstraintPoint* cp1 = vc->points + 0;
			b2VelocityConstraintPoint* cp2 = vc->points + 1;

			b2Vec2 a(cp1->normalImpulse, cp2->normalImpulse);
			b2Assert(a.x >= 0.0f && a.y >= 0.0f);

			// Relative velocity at contact
			b2Vec2 dv1 = vB + b2Cross(wB, cp1->rB) - vA - b2Cross(wA, cp1->rA);
			b2Vec2 dv2 = vB + b2Cross(wB, cp2->rB) - vA - b2Cross(wA, cp2->rA);

			// Compute normal velocity
			float32 vn1 = b2Dot(dv1, normal);
			float32 vn2 = b2Dot(dv2, normal);

			b2Vec2 b;
			b.x = vn1 - cp1->velocityBias;
			b.y = vn2 - cp2->velocityBias;

			// Compute b'
			b -= b2Mul(vc->K, a);

			const float32 k_errorTol = 1e-3f;
			B2_NOT_USED(k_errorTol);

			for (;;)
			{
				//
				// Case 1: vn = 0
				//
				// 0 = A * x + b'
				//
				// Solve for x:
				//
				// x = - inv(A) * b'
				//
				b2Vec2 x = - b2Mul(vc->normalMass, b);

				if (x.x >= 0.0f && x.y >= 0.0f)
				{
					// Get the incremental impulse
					b2Vec2 d = x - a;

					// Apply incremental impulse
					b2Vec2 P1 = d.x * normal;
					b2Vec2 P2 = d.y * normal;
					vA -= mA * (P1 + P2);
					wA -= iA * (b2Cross(cp1->rA, P1) + b2Cross(cp2->rA, P2));

					vB += mB * (P1 + P2);
					wB += iB * (b2Cross(cp1->rB, P1) + b2Cross(cp2->rB, P2));

					// Accumulate
					cp1->normalImpulse = x.x;
					cp2->normalImpulse = x.y;

#if B2_DEBUG_SOLVER == 1
					// Postconditions
					dv1 = vB + b2Cross(wB, cp1->rB) - vA - b2Cross(wA, cp1->rA);
					dv2 = vB + b2Cross(wB, cp2->rB) - vA - b2Cross(wA, cp2->rA);

					// Compute normal velocity
					vn1 = b2Dot(dv1, normal);
					vn2 = b2Dot(dv2, normal);

					b2Assert(b2Abs(vn1 - cp1->velocityBias) < k_errorTol);
					b2Assert(b2Abs(vn2 - cp2->velocityBias) < k_errorTol);
#endif
					break;
				}

				//
				// Case 2: vn1 = 0 and x2 = 0
				//
				//   0 = a11 * x1 + a12 * 0 + b1' 
				// vn2 = a21 * x1 + a22 * 0 + b2'
				//
				x.x = - cp1->normalMass * b.x;
				x.y = 0.0f;
				vn2 = vc->K.ex.y * x.x + b.y;

				if (x.x >= 0.0f && vn2 >= 0.0f)
				{
					// Get the incremental impulse
					b2Vec2 d = x - a;

					// Apply incremental impulse
					b2Vec2 P1 = d.x * normal;
					b2Vec2 P2 = d.y * normal;
					vA -= mA * (P1 + P2);
					wA -= iA * (b2Cross(cp1->rA, P1) + b2Cross(cp2->rA, P2));

					vB += mB * (P1 + P2);
					wB += iB * (b2Cross(cp1->rB, P1) + b2Cross(cp2->rB, P2));

					// Accumulate
					cp1->normalImpulse = x.x;
					cp2->normalImpulse = x.y;

#if B2_DEBUG_SOLVER == 1
					// Postconditions
					dv1 = vB + b2Cross(wB, cp1->rB) - vA - b2Cross(wA, cp1->rA);

					// Compute normal velocity
					vn1 = b2Dot(dv1, normal);

					b2Assert(b2Abs(vn1 - cp1->velocityBias) < k_errorTol);
#endif
					break;
				}


				//
				// Case 3: vn2 = 0 and x1 = 0
				//
				// vn1 = a11 * 0 + a12 * x2 + b1' 
				//   0 = a21 * 0 + a22 * x2 + b2'
				//
				x.x = 0.0f;
				x.y = - cp2->normalMass * b.y;
				vn1 = vc->K.ey.x * x.y + b.x;

				if (x.y >= 0.0f && vn1 >= 0.0f)
				{
					// Resubstitute for the incremental impulse
					b2Vec2 d = x - a;

					// Apply incremental impulse
					b2Vec2 P1 = d.x * normal;
					b2Vec2 P2 = d.y * normal;
					vA -= mA * (P1 + P2);
					wA -= iA * (b2Cross(cp1->rA, P1) + b2Cross(cp2->rA, P2));

					vB += mB * (P1 + P2);
					wB += iB * (b2Cross(cp1->rB, P1) + b2Cross(cp2->rB, P2));

					// Accumulate
					cp1->normalImpulse = x.x;
					cp2->normalImpulse = x.y;

#if B2_DEBUG_SOLVER == 1
					// Postconditions
					dv2 = vB + b2Cross(wB, cp2->rB) - vA - b2Cross(wA, cp2->rA);

					// Compute normal velocity
					vn2 = b2Dot(dv2, normal);

					b2Assert(b2Abs(vn2 - cp2->velocityBias) < k_errorTol);
#endif
					break;
				}

				//
				// Case 4: x1 = 0 and x2 = 0
				// 
				// vn1 = b1
				// vn2 = b2;
				x.x = 0.0f;
				x.y = 0.0f;
				vn1 = b.x;
				vn2 = b.y;

				if (vn1 >= 0.0f && vn2 >= 0.0f )
				{
					// Resubstitute for the incremental impulse
					b2Vec2 d = x - a;

					// Apply incremental impulse
					b2Vec2 P1 = d.x * normal;
					b2Vec2 P2 = d.y * normal;
					vA -= mA * (P1 + P2);
					wA -= iA * (b2Cross(cp1->rA, P1) + b2Cross(cp2->rA, P2));

					vB += mB * (P1 + P2);
					wB += iB * (b2Cross(cp1->rB, P1) + b2Cross(cp2->rB, P2));

					// Accumulate
					cp1->normalImpulse = x.x;
					cp2->normalImpulse = x.y;

					break;
				}

				// No solution, give up. This is hit sometimes, but it doesn't seem to matter.
				break;
			}
		}

		m_velocities[indexA].v = vA;
		m_velocities[indexA].w = wA;
		m_velocities[indexB].v = vB;
		m_velocities[indexB].w = wB;
	}
}

void b2ContactSolver::StoreImpulses()
{
	for (int32 i = 0; i < m_count; ++i)
	{
		b2ContactVelocityConstraint* vc = m_velocityConstraints + i;
		b2Manifold* manifold = m_contacts[vc->contactIndex]->GetManifold();

		for (int32 j = 0; j < vc->pointCount; ++j)
		{
			manifold->points[j].normalImpulse = vc->points[j].normalImpulse;
			manifold->points[j].tangentImpulse = vc->points[j].tangentImpulse;
		}
	}
}

struct b2PositionSolverManifold
{
	void Initialize(b2ContactPositionConstraint* pc, const b2Transform& xfA, const b2Transform& xfB, int32 index)
	{
		b2Assert(pc->pointCount > 0);

		switch (pc->type)
		{
		case b2Manifold::e_circles:
			{
				b2Vec2 pointA = b2Mul(xfA, pc->localPoint);
				b2Vec2 pointB = b2Mul(xfB, pc->localPoints[0]);
				normal = pointB - pointA;
				normal.Normalize();
				point = 0.5f * (pointA + pointB);
				separation = b2Dot(pointB - pointA, normal) - pc->radiusA - pc->radiusB;
			}
			break;

		case b2Manifold::e_faceA:
			{
				normal = b2Mul(xfA.q, pc->localNormal);
				b2Vec2 planePoint = b2Mul(xfA, pc->localPoint);

				b2Vec2 clipPoint = b2Mul(xfB, pc->localPoints[index]);
				separation = b2Dot(clipPoint - planePoint, normal) - pc->radiusA - pc->radiusB;
				point = clipPoint;
			}
			break;

		case b2Manifold::e_faceB:
			{
				normal = b2Mul(xfB.q, pc->localNormal);
				b2Vec2 planePoint = b2Mul(xfB, pc->localPoint);

				b2Vec2 clipPoint = b2Mul(xfA, pc->localPoints[index]);
				separation = b2Dot(clipPoint - planePoint, normal) - pc->radiusA - pc->radiusB;
				point = clipPoint;

				// Ensure normal points from A to B
				normal = -normal;
			}
			break;
		default:
			{
				// This shouldn't be executed if pc->type is valid.
				separation = 0.0f;
				normal = b2Vec2_zero;
				point = b2Vec2_zero;
				b2Assert(false);
			}
			break;
		}
	}

	b2Vec2 normal;
	b2Vec2 point;
	float32 separation;
};

// Sequential solver.
bool b2ContactSolver::SolvePositionConstraints()
{
	float32 minSeparation = 0.0f;

	for (int32 i = 0; i < m_count; ++i)
	{
		b2ContactPositionConstraint* pc = m_positionConstraints + i;

		int32 indexA = pc->indexA;
		int32 indexB = pc->indexB;
		b2Vec2 localCenterA = pc->localCenterA;
		float32 mA = pc->invMassA;
		float32 iA = pc->invIA;
		b2Vec2 localCenterB = pc->localCenterB;
		float32 mB = pc->invMassB;
		float32 iB = pc->invIB;
		int32 pointCount = pc->pointCount;

		b2Vec2 cA = m_positions[indexA].c;
		float32 aA = m_positions[indexA].a;

		b2Vec2 cB = m_positions[indexB].c;
		float32 aB = m_positions[indexB].a;

		// Solve normal constraints
		for (int32 j = 0; j < pointCount; ++j)
		{
			b2Transform xfA, xfB;
			xfA.q.Set(aA);
			xfB.q.Set(aB);
			xfA.p = cA - b2Mul(xfA.q, localCenterA);
			xfB.p = cB - b2Mul(xfB.q, localCenterB);

			b2PositionSolverManifold psm;
			psm.Initialize(pc, xfA, xfB, j);
			b2Vec2 normal = psm.normal;

			b2Vec2 point = psm.point;
			float32 separation = psm.separation;

			b2Vec2 rA = point - cA;
			b2Vec2 rB = point - cB;

			// Track max constraint error.
			minSeparation = b2Min(minSeparation, separation);

			// Prevent large corrections and allow slop.
			float32 C = b2Clamp(b2_baumgarte * (separation + b2_linearSlop), -b2_maxLinearCorrection, 0.0f);

			// Compute the effective mass.
			float32 rnA = b2Cross(rA, normal);
			float32 rnB = b2Cross(rB, normal);
			float32 K = mA + mB + iA * rnA * rnA + iB * rnB * rnB;

			// Compute normal impulse
			float32 impulse = K > 0.0f ? - C / K : 0.0f;

			b2Vec2 P = impulse * normal;

			cA -= mA * P;
			aA -= iA * b2Cross(rA, P);

			cB += mB * P;
			aB += iB * b2Cross(rB, P);
		}

		m_positions[indexA].c = cA;
		m_positions[indexA].a = aA;

		m_positions[indexB].c = cB;
		m_positions[indexB].a = aB;
	}

	// We can't expect minSpeparation >= -b2_linearSlop because we don't
	// push the separation above -b2_linearSlop.
	return minSeparation >= -3.0f * b2_linearSlop;
}

// Sequential position solver for position constraints.
bool b2ContactSolver::SolveTOIPositionConstraints(int32 toiIndexA, int32 toiIndexB)
{
	float32 minSeparation = 0.0f;

	for (int32 i = 0; i < m_count; ++i)
	{
		b2ContactPositionConstraint* pc = m_positionConstraints + i;

		int32 indexA = pc->indexA;
		int32 indexB = pc->indexB;
		b2Vec2 localCenterA = pc->localCenterA;
		b2Vec2 localCenterB = pc->localCenterB;
		int32 pointCount = pc->pointCount;

		float32 mA = 0.0f;
		float32 iA = 0.0f;
		if (indexA == toiIndexA || indexA == toiIndexB)
		{
			mA = pc->invMassA;
			iA = pc->invIA;
		}

		float32 mB = 0.0f;
		float32 iB = 0.;
		if (indexB == toiIndexA || indexB == toiIndexB)
		{
			mB = pc->invMassB;
			iB = pc->invIB;
		}

		b2Vec2 cA = m_positions[indexA].c;
		float32 aA = m_positions[indexA].a;

		b2Vec2 cB = m_positions[indexB].c;
		float32 aB = m_positions[indexB].a;

		// Solve normal constraints
		for (int32 j = 0; j < pointCount; ++j)
		{
			b2Transform xfA, xfB;
			xfA.q.Set(aA);
			xfB.q.Set(aB);
			xfA.p = cA - b2Mul(xfA.q, localCenterA);
			xfB.p = cB - b2Mul(xfB.q, localCenterB);

			b2PositionSolverManifold psm;
			psm.Initialize(pc, xfA, xfB, j);
			b2Vec2 normal = psm.normal;

			b2Vec2 point = psm.point;
			float32 separation = psm.separation;

			b2Vec2 rA = point - cA;
			b2Vec2 rB = point - cB;

			// Track max constraint error.
			minSeparation = b2Min(minSeparation, separation);

			// Prevent large corrections and allow slop.
			float32 C = b2Clamp(b2_toiBaugarte * (separation + b2_linearSlop), -b2_maxLinearCorrection, 0.0f);

			// Compute the effective mass.
			float32 rnA = b2Cross(rA, normal);
			float32 rnB = b2Cross(rB, normal);
			float32 K = mA + mB + iA * rnA * rnA + iB * rnB * rnB;

			// Compute normal impulse
			float32 impulse = K > 0.0f ? - C / K : 0.0f;

			b2Vec2 P = impulse * normal;

			cA -= mA * P;
			aA -= iA * b2Cross(rA, P);

			cB += mB * P;
			aB += iB * b2Cross(rB, P);
		}

		m_positions[indexA].c = cA;
		m_positions[indexA].a = aA;

		m_positions[indexB].c = cB;
		m_positions[indexB].a = aB;
	}

	// We can't expect minSpeparation >= -b2_linearSlop because we don't
	// push the separation above -b2_linearSlop.
	return minSeparation >= -1.5f * b2_linearSlop;
}

// end of ContactSolver.cpp

class b2BlockAllocator;

class b2EdgeAndCircleContact : public b2Contact
{
public:
	static b2Contact* Create(	b2Fixture* fixtureA, int32 indexA,
								b2Fixture* fixtureB, int32 indexB, b2BlockAllocator* allocator);
	static void Destroy(b2Contact* contact, b2BlockAllocator* allocator);

	b2EdgeAndCircleContact(b2Fixture* fixtureA, b2Fixture* fixtureB);
	~b2EdgeAndCircleContact() {}

	void Evaluate(b2Manifold* manifold, const b2Transform& xfA, const b2Transform& xfB);
};

// end of EdgeAndCircleContact.h

b2Contact* b2EdgeAndCircleContact::Create(b2Fixture* fixtureA, int32, b2Fixture* fixtureB, int32, b2BlockAllocator* allocator)
{
	void* mem = allocator->Allocate(sizeof(b2EdgeAndCircleContact));
	return new (mem) b2EdgeAndCircleContact(fixtureA, fixtureB);
}

void b2EdgeAndCircleContact::Destroy(b2Contact* contact, b2BlockAllocator* allocator)
{
	((b2EdgeAndCircleContact*)contact)->~b2EdgeAndCircleContact();
	allocator->Free(contact, sizeof(b2EdgeAndCircleContact));
}

b2EdgeAndCircleContact::b2EdgeAndCircleContact(b2Fixture* fixtureA, b2Fixture* fixtureB)
: b2Contact(fixtureA, 0, fixtureB, 0)
{
	b2Assert(m_fixtureA->GetType() == b2Shape::e_edge);
	b2Assert(m_fixtureB->GetType() == b2Shape::e_circle);
}

void b2EdgeAndCircleContact::Evaluate(b2Manifold* manifold, const b2Transform& xfA, const b2Transform& xfB)
{
	b2CollideEdgeAndCircle(	manifold,
								(b2EdgeShape*)m_fixtureA->GetShape(), xfA,
								(b2CircleShape*)m_fixtureB->GetShape(), xfB);
}

// end of EdgeAndCircleContact.cpp

class b2BlockAllocator;

class b2EdgeAndPolygonContact : public b2Contact
{
public:
	static b2Contact* Create(	b2Fixture* fixtureA, int32 indexA,
								b2Fixture* fixtureB, int32 indexB, b2BlockAllocator* allocator);
	static void Destroy(b2Contact* contact, b2BlockAllocator* allocator);

	b2EdgeAndPolygonContact(b2Fixture* fixtureA, b2Fixture* fixtureB);
	~b2EdgeAndPolygonContact() {}

	void Evaluate(b2Manifold* manifold, const b2Transform& xfA, const b2Transform& xfB);
};

// end of EdgeAndPolygonContact.h 

b2Contact* b2EdgeAndPolygonContact::Create(b2Fixture* fixtureA, int32, b2Fixture* fixtureB, int32, b2BlockAllocator* allocator)
{
	void* mem = allocator->Allocate(sizeof(b2EdgeAndPolygonContact));
	return new (mem) b2EdgeAndPolygonContact(fixtureA, fixtureB);
}

void b2EdgeAndPolygonContact::Destroy(b2Contact* contact, b2BlockAllocator* allocator)
{
	((b2EdgeAndPolygonContact*)contact)->~b2EdgeAndPolygonContact();
	allocator->Free(contact, sizeof(b2EdgeAndPolygonContact));
}

b2EdgeAndPolygonContact::b2EdgeAndPolygonContact(b2Fixture* fixtureA, b2Fixture* fixtureB)
: b2Contact(fixtureA, 0, fixtureB, 0)
{
	b2Assert(m_fixtureA->GetType() == b2Shape::e_edge);
	b2Assert(m_fixtureB->GetType() == b2Shape::e_polygon);
}

void b2EdgeAndPolygonContact::Evaluate(b2Manifold* manifold, const b2Transform& xfA, const b2Transform& xfB)
{
	b2CollideEdgeAndPolygon(	manifold,
								(b2EdgeShape*)m_fixtureA->GetShape(), xfA,
								(b2PolygonShape*)m_fixtureB->GetShape(), xfB);
}

// end of EdgeAndPolygonContact.cpp

class b2BlockAllocator;

class b2PolygonAndCircleContact : public b2Contact
{
public:
	static b2Contact* Create(b2Fixture* fixtureA, int32 indexA, b2Fixture* fixtureB, int32 indexB, b2BlockAllocator* allocator);
	static void Destroy(b2Contact* contact, b2BlockAllocator* allocator);

	b2PolygonAndCircleContact(b2Fixture* fixtureA, b2Fixture* fixtureB);
	~b2PolygonAndCircleContact() {}

	void Evaluate(b2Manifold* manifold, const b2Transform& xfA, const b2Transform& xfB);
};

// end of PolygonAndCircleContact.h

b2Contact* b2PolygonAndCircleContact::Create(b2Fixture* fixtureA, int32, b2Fixture* fixtureB, int32, b2BlockAllocator* allocator)
{
	void* mem = allocator->Allocate(sizeof(b2PolygonAndCircleContact));
	return new (mem) b2PolygonAndCircleContact(fixtureA, fixtureB);
}

void b2PolygonAndCircleContact::Destroy(b2Contact* contact, b2BlockAllocator* allocator)
{
	((b2PolygonAndCircleContact*)contact)->~b2PolygonAndCircleContact();
	allocator->Free(contact, sizeof(b2PolygonAndCircleContact));
}

b2PolygonAndCircleContact::b2PolygonAndCircleContact(b2Fixture* fixtureA, b2Fixture* fixtureB)
: b2Contact(fixtureA, 0, fixtureB, 0)
{
	b2Assert(m_fixtureA->GetType() == b2Shape::e_polygon);
	b2Assert(m_fixtureB->GetType() == b2Shape::e_circle);
}

void b2PolygonAndCircleContact::Evaluate(b2Manifold* manifold, const b2Transform& xfA, const b2Transform& xfB)
{
	b2CollidePolygonAndCircle(	manifold,
								(b2PolygonShape*)m_fixtureA->GetShape(), xfA,
								(b2CircleShape*)m_fixtureB->GetShape(), xfB);
}

// end of PolygonAndCircleContact.cpp

class b2BlockAllocator;

class b2PolygonContact : public b2Contact
{
public:
	static b2Contact* Create(	b2Fixture* fixtureA, int32 indexA,
								b2Fixture* fixtureB, int32 indexB, b2BlockAllocator* allocator);
	static void Destroy(b2Contact* contact, b2BlockAllocator* allocator);

	b2PolygonContact(b2Fixture* fixtureA, b2Fixture* fixtureB);
	~b2PolygonContact() {}

	void Evaluate(b2Manifold* manifold, const b2Transform& xfA, const b2Transform& xfB);
};

// end of PolygonContact.h

b2Contact* b2PolygonContact::Create(b2Fixture* fixtureA, int32, b2Fixture* fixtureB, int32, b2BlockAllocator* allocator)
{
	void* mem = allocator->Allocate(sizeof(b2PolygonContact));
	return new (mem) b2PolygonContact(fixtureA, fixtureB);
}

void b2PolygonContact::Destroy(b2Contact* contact, b2BlockAllocator* allocator)
{
	((b2PolygonContact*)contact)->~b2PolygonContact();
	allocator->Free(contact, sizeof(b2PolygonContact));
}

b2PolygonContact::b2PolygonContact(b2Fixture* fixtureA, b2Fixture* fixtureB)
	: b2Contact(fixtureA, 0, fixtureB, 0)
{
	b2Assert(m_fixtureA->GetType() == b2Shape::e_polygon);
	b2Assert(m_fixtureB->GetType() == b2Shape::e_polygon);
}

void b2PolygonContact::Evaluate(b2Manifold* manifold, const b2Transform& xfA, const b2Transform& xfB)
{
	b2CollidePolygons(	manifold,
						(b2PolygonShape*)m_fixtureA->GetShape(), xfA,
						(b2PolygonShape*)m_fixtureB->GetShape(), xfB);
}

// end of PolygonContact.cpp

b2ContactRegister b2Contact::s_registers[b2Shape::e_typeCount][b2Shape::e_typeCount];
bool b2Contact::s_initialized = false;

void b2Contact::InitializeRegisters()
{
	AddType(b2CircleContact::Create, b2CircleContact::Destroy, b2Shape::e_circle, b2Shape::e_circle);
	AddType(b2PolygonAndCircleContact::Create, b2PolygonAndCircleContact::Destroy, b2Shape::e_polygon, b2Shape::e_circle);
	AddType(b2PolygonContact::Create, b2PolygonContact::Destroy, b2Shape::e_polygon, b2Shape::e_polygon);
	AddType(b2EdgeAndCircleContact::Create, b2EdgeAndCircleContact::Destroy, b2Shape::e_edge, b2Shape::e_circle);
	AddType(b2EdgeAndPolygonContact::Create, b2EdgeAndPolygonContact::Destroy, b2Shape::e_edge, b2Shape::e_polygon);
	AddType(b2ChainAndCircleContact::Create, b2ChainAndCircleContact::Destroy, b2Shape::e_chain, b2Shape::e_circle);
	AddType(b2ChainAndPolygonContact::Create, b2ChainAndPolygonContact::Destroy, b2Shape::e_chain, b2Shape::e_polygon);
}

void b2Contact::AddType(b2ContactCreateFcn* createFcn, b2ContactDestroyFcn* destoryFcn,
						b2Shape::Type type1, b2Shape::Type type2)
{
	b2Assert(0 <= type1 && type1 < b2Shape::e_typeCount);
	b2Assert(0 <= type2 && type2 < b2Shape::e_typeCount);
	
	s_registers[type1][type2].createFcn = createFcn;
	s_registers[type1][type2].destroyFcn = destoryFcn;
	s_registers[type1][type2].primary = true;

	if (type1 != type2)
	{
		s_registers[type2][type1].createFcn = createFcn;
		s_registers[type2][type1].destroyFcn = destoryFcn;
		s_registers[type2][type1].primary = false;
	}
}

b2Contact* b2Contact::Create(b2Fixture* fixtureA, int32 indexA, b2Fixture* fixtureB, int32 indexB, b2BlockAllocator* allocator)
{
	if (s_initialized == false)
	{
		InitializeRegisters();
		s_initialized = true;
	}

	b2Shape::Type type1 = fixtureA->GetType();
	b2Shape::Type type2 = fixtureB->GetType();

	b2Assert(0 <= type1 && type1 < b2Shape::e_typeCount);
	b2Assert(0 <= type2 && type2 < b2Shape::e_typeCount);
	
	b2ContactCreateFcn* createFcn = s_registers[type1][type2].createFcn;
	if (createFcn)
	{
		if (s_registers[type1][type2].primary)
		{
			return createFcn(fixtureA, indexA, fixtureB, indexB, allocator);
		}
		else
		{
			return createFcn(fixtureB, indexB, fixtureA, indexA, allocator);
		}
	}
	else
	{
		return NULL;
	}
}

void b2Contact::Destroy(b2Contact* contact, b2BlockAllocator* allocator)
{
	b2Assert(s_initialized == true);

	b2Fixture* fixtureA = contact->m_fixtureA;
	b2Fixture* fixtureB = contact->m_fixtureB;

	if (contact->m_manifold.pointCount > 0 &&
		fixtureA->IsSensor() == false &&
		fixtureB->IsSensor() == false)
	{
		fixtureA->GetBody()->SetAwake(true);
		fixtureB->GetBody()->SetAwake(true);
	}

	b2Shape::Type typeA = fixtureA->GetType();
	b2Shape::Type typeB = fixtureB->GetType();

	b2Assert(0 <= typeA && typeB < b2Shape::e_typeCount);
	b2Assert(0 <= typeA && typeB < b2Shape::e_typeCount);

	b2ContactDestroyFcn* destroyFcn = s_registers[typeA][typeB].destroyFcn;
	destroyFcn(contact, allocator);
}

b2Contact::b2Contact(b2Fixture* fA, int32 indexA, b2Fixture* fB, int32 indexB)
{
	m_flags = e_enabledFlag;

	m_fixtureA = fA;
	m_fixtureB = fB;

	m_indexA = indexA;
	m_indexB = indexB;

	m_manifold.pointCount = 0;

	m_prev = NULL;
	m_next = NULL;

	m_nodeA.contact = NULL;
	m_nodeA.prev = NULL;
	m_nodeA.next = NULL;
	m_nodeA.other = NULL;

	m_nodeB.contact = NULL;
	m_nodeB.prev = NULL;
	m_nodeB.next = NULL;
	m_nodeB.other = NULL;

	m_toiCount = 0;

	m_friction = b2MixFriction(m_fixtureA->m_friction, m_fixtureB->m_friction);
	m_restitution = b2MixRestitution(m_fixtureA->m_restitution, m_fixtureB->m_restitution);

	m_tangentSpeed = 0.0f;
}

// Update the contact manifold and touching status.
// Note: do not assume the fixture AABBs are overlapping or are valid.
void b2Contact::Update(b2ContactListener* listener)
{
	b2Manifold oldManifold = m_manifold;

	// Re-enable this contact.
	m_flags |= e_enabledFlag;

	bool touching = false;
	bool wasTouching = (m_flags & e_touchingFlag) == e_touchingFlag;

	bool sensorA = m_fixtureA->IsSensor();
	bool sensorB = m_fixtureB->IsSensor();
	bool sensor = sensorA || sensorB;

	b2Body* bodyA = m_fixtureA->GetBody();
	b2Body* bodyB = m_fixtureB->GetBody();
	const b2Transform& xfA = bodyA->GetTransform();
	const b2Transform& xfB = bodyB->GetTransform();

	// Is this contact a sensor?
	if (sensor)
	{
		const b2Shape* shapeA = m_fixtureA->GetShape();
		const b2Shape* shapeB = m_fixtureB->GetShape();
		touching = b2TestOverlap(shapeA, m_indexA, shapeB, m_indexB, xfA, xfB);

		// Sensors don't generate manifolds.
		m_manifold.pointCount = 0;
	}
	else
	{
		Evaluate(&m_manifold, xfA, xfB);
		touching = m_manifold.pointCount > 0;

		// Match old contact ids to new contact ids and copy the
		// stored impulses to warm start the solver.
		for (int32 i = 0; i < m_manifold.pointCount; ++i)
		{
			b2ManifoldPoint* mp2 = m_manifold.points + i;
			mp2->normalImpulse = 0.0f;
			mp2->tangentImpulse = 0.0f;
			b2ContactID id2 = mp2->id;

			for (int32 j = 0; j < oldManifold.pointCount; ++j)
			{
				b2ManifoldPoint* mp1 = oldManifold.points + j;

				if (mp1->id.key == id2.key)
				{
					mp2->normalImpulse = mp1->normalImpulse;
					mp2->tangentImpulse = mp1->tangentImpulse;
					break;
				}
			}
		}

		if (touching != wasTouching)
		{
			bodyA->SetAwake(true);
			bodyB->SetAwake(true);
		}
	}

	if (touching)
	{
		m_flags |= e_touchingFlag;
	}
	else
	{
		m_flags &= ~e_touchingFlag;
	}

	if (wasTouching == false && touching == true && listener)
	{
		listener->BeginContact(this);
	}

	if (wasTouching == true && touching == false && listener)
	{
		listener->EndContact(this);
	}

	if (sensor == false && touching && listener)
	{
		listener->PreSolve(this, &oldManifold);
	}
}

// end of Contact.cpp

struct b2ParticleContact;

struct FindContactCheck
{
    uint16 particleIndex;
    uint16 comparatorIndex;
};

struct FindContactInput
{
    uint32 proxyIndex;
    b2Vec2 position;
};

enum { NUM_V32_SLOTS = 4 };

#ifdef __cplusplus
extern "C" {
#endif

extern int CalculateTags_Simd(const b2Vec2* positions,
                              int count,
                              const float& inverseDiameter,
                              uint32* outTags);

extern void FindContactsFromChecks_Simd(
	const FindContactInput* reordered,
	const FindContactCheck* checks,
	int numChecks,
  const float& particleDiameterSq,
  const float& particleDiameterInv,
  const uint32* flags,
	b2GrowableBuffer<b2ParticleContact>& contacts);

#ifdef __cplusplus
} // extern "C"
#endif

// end of ParticleAssembly.h

extern "C" {

// Helper function, called from assembly routine.
void GrowParticleContactBuffer(
	b2GrowableBuffer<b2ParticleContact>& contacts)
{
	// Set contacts.count = capacity instead of count because there are
	// items past the end of the array waiting to be post-processed.
	// We must maintain the entire contacts array.
	// TODO: It would be better to have the items awaiting post-processing
	// in their own array on the stack.
	contacts.SetCount(contacts.GetCapacity());
	contacts.Grow();
}

} // extern "C"

// end of ParticleAssembly.cpp

#define B2PARTICLECOLOR_BITS_PER_COMPONENT (sizeof(uint8) << 3)
// Maximum value of a b2ParticleColor component.
#define B2PARTICLECOLOR_MAX_VALUE \
	((1U << B2PARTICLECOLOR_BITS_PER_COMPONENT) - 1)

/// Number of bits used to store each b2ParticleColor component.
const uint8 b2ParticleColor::k_bitsPerComponent =
	B2PARTICLECOLOR_BITS_PER_COMPONENT;
const float32 b2ParticleColor::k_maxValue = (float)B2PARTICLECOLOR_MAX_VALUE;
const float32 b2ParticleColor::k_inverseMaxValue =
	1.0f / (float)B2PARTICLECOLOR_MAX_VALUE;

b2ParticleColor b2ParticleColor_zero(0, 0, 0, 0);

b2ParticleColor::b2ParticleColor(const b2Color& color)
{
	Set(color);
}

b2Color b2ParticleColor::GetColor() const
{
	return b2Color(k_inverseMaxValue * r,
				   k_inverseMaxValue * g,
				   k_inverseMaxValue * b);
}

void b2ParticleColor::Set(const b2Color& color)
{
	Set((uint8)(k_maxValue * color.r),
		(uint8)(k_maxValue * color.g),
		(uint8)(k_maxValue * color.b),
		B2PARTICLECOLOR_MAX_VALUE);
}

int32 b2CalculateParticleIterations(
	float32 gravity, float32 radius, float32 timeStep)
{
	// In some situations you may want more particle iterations than this,
	// but to avoid excessive cycle cost, don't recommend more than this.
	const int32 B2_MAX_RECOMMENDED_PARTICLE_ITERATIONS = 8;
	const float32 B2_RADIUS_THRESHOLD = 0.01f;
	int32 iterations =
		(int32) ceilf(b2Sqrt(gravity / (B2_RADIUS_THRESHOLD * radius)) * timeStep);
	return b2Clamp(iterations, 1, B2_MAX_RECOMMENDED_PARTICLE_ITERATIONS);
}
 
// end of Particle.cpp

class b2Shape;
class b2World;
class b2ParticleSystem;
class b2ParticleGroup;
class b2ParticleColor;
#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
class b2CircleShape;
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

/// @file

/// The particle group type.  Can be combined with the | operator.
enum b2ParticleGroupFlag
{
	/// Prevents overlapping or leaking.
	b2_solidParticleGroup = 1 << 0,
	/// Keeps its shape.
	b2_rigidParticleGroup = 1 << 1,
	/// Won't be destroyed if it gets empty.
	b2_particleGroupCanBeEmpty = 1 << 2,
	/// Will be destroyed on next simulation step.
	b2_particleGroupWillBeDestroyed = 1 << 3,
	/// Updates depth data on next simulation step.
	b2_particleGroupNeedsUpdateDepth = 1 << 4,
	b2_particleGroupInternalMask =
		b2_particleGroupWillBeDestroyed |
		b2_particleGroupNeedsUpdateDepth,
};

/// A particle group definition holds all the data needed to construct a
/// particle group.  You can safely re-use these definitions.
struct b2ParticleGroupDef
{

	b2ParticleGroupDef()
	{
		flags = 0;
		groupFlags = 0;
		position = b2Vec2_zero;
		angle = 0;
		linearVelocity = b2Vec2_zero;
		angularVelocity = 0;
		color = b2ParticleColor_zero;
		strength = 1;
		shape = NULL;
		shapes = NULL;
		shapeCount = 0;
		stride = 0;
		particleCount = 0;
		positionData = NULL;
		lifetime = 0.0f;
		userData = NULL;
		group = NULL;

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
		circleShapes = NULL;
		ownShapesArray = false;
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API
	}

	~b2ParticleGroupDef()
	{
#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
		FreeShapesMemory();
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API
	}

	/// The particle-behavior flags (See #b2ParticleFlag).
	uint32 flags;

	/// The group-construction flags (See #b2ParticleGroupFlag).
	uint32 groupFlags;

	/// The world position of the group.
	/// Moves the group's shape a distance equal to the value of position.
	b2Vec2 position;

	/// The world angle of the group in radians.
	/// Rotates the shape by an angle equal to the value of angle.
	float32 angle;

	/// The linear velocity of the group's origin in world co-ordinates.
	b2Vec2 linearVelocity;

	/// The angular velocity of the group.
	float32 angularVelocity;

	/// The color of all particles in the group.
	b2ParticleColor color;

	/// The strength of cohesion among the particles in a group with flag
	/// b2_elasticParticle or b2_springParticle.
	float32 strength;

	/// The shape where particles will be added.
	const b2Shape* shape;

	/// A array of shapes where particles will be added.
	const b2Shape* const* shapes;

	/// The number of shapes.
	int32 shapeCount;

	/// The interval of particles in the shape.
	/// If it is 0, b2_particleStride * particleDiameter is used instead.
	float32 stride;

	/// The number of particles in addition to ones added in the shape.
	int32 particleCount;

	/// The initial positions of the particleCount particles.
	const b2Vec2* positionData;

	/// Lifetime of the particle group in seconds.  A value <= 0.0f indicates a
	/// particle group with infinite lifetime.
	float32 lifetime;

	/// Use this to store application-specific group data.
	void* userData;

	/// An existing particle group to which the particles will be added.
	b2ParticleGroup* group;

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
	/// Storage for constructed CircleShapes from an incoming vertex list
	const b2CircleShape* circleShapes;

	/// True if we create the shapes array internally.
	bool ownShapesArray;

	/// Clean up all memory associated with SetCircleShapesFromVertexList
	void FreeShapesMemory();

	/// From a vertex list created by an external language API, construct
	/// a list of circle shapes that can be used to create a b2ParticleGroup
	/// This eliminates cumbersome array-interfaces between languages.
	void SetCircleShapesFromVertexList(void* inBuf,
									   int numShapes,
									   float radius);

	/// Set position with direct floats.
	void SetPosition(float32 x, float32 y);

	/// Set color with direct ints.
	void SetColor(int32 r, int32 g, int32 b, int32 a);
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API
};

/// A group of particles. b2ParticleGroup::CreateParticleGroup creates these.
class b2ParticleGroup
{

public:

	/// Get the next particle group from the list in b2_World.
	b2ParticleGroup* GetNext();
	const b2ParticleGroup* GetNext() const;

	/// Get the particle system that holds this particle group.
	b2ParticleSystem* GetParticleSystem();
	const b2ParticleSystem* GetParticleSystem() const;

	/// Get the number of particles.
	int32 GetParticleCount() const;

	/// Get the offset of this group in the global particle buffer
	int32 GetBufferIndex() const;

	/// Does this group contain the particle.
	bool ContainsParticle(int32 index) const;

	/// Get the logical sum of particle flags.
	uint32 GetAllParticleFlags() const;

	/// Get the construction flags for the group.
	uint32 GetGroupFlags() const;

	/// Set the construction flags for the group.
	void SetGroupFlags(uint32 flags);

	/// Get the total mass of the group: the sum of all particles in it.
	float32 GetMass() const;

	/// Get the moment of inertia for the group.
	float32 GetInertia() const;

	/// Get the center of gravity for the group.
	b2Vec2 GetCenter() const;

	/// Get the linear velocity of the group.
	b2Vec2 GetLinearVelocity() const;

	/// Get the angular velocity of the group.
	float32 GetAngularVelocity() const;

	/// Get the position of the group's origin and rotation.
	/// Used only with groups of rigid particles.
	const b2Transform& GetTransform() const;

	/// Get position of the particle group as a whole.
	/// Used only with groups of rigid particles.
	const b2Vec2& GetPosition() const;

	/// Get the rotational angle of the particle group as a whole.
	/// Used only with groups of rigid particles.
	float32 GetAngle() const;

	/// Get the world linear velocity of a world point, from the average linear
	/// and angular velocities of the particle group.
	/// @param a point in world coordinates.
	/// @return the world velocity of a point.
	b2Vec2 GetLinearVelocityFromWorldPoint(const b2Vec2& worldPoint) const;

	/// Get the user data pointer that was provided in the group definition.
	void* GetUserData() const;

	/// Set the user data. Use this to store your application specific data.
	void SetUserData(void* data);

	/// Call b2ParticleSystem::ApplyForce for every particle in the group.
	void ApplyForce(const b2Vec2& force);

	/// Call b2ParticleSystem::ApplyLinearImpulse for every particle in the
	/// group.
	void ApplyLinearImpulse(const b2Vec2& impulse);

	/// Destroy all the particles in this group.
	/// This function is locked during callbacks.
	/// @param Whether to call the world b2DestructionListener for each
	/// particle is destroyed.
	/// @warning This function is locked during callbacks.
	void DestroyParticles(bool callDestructionListener);

	/// Destroy all particles in this group without enabling the destruction
	/// callback for destroyed particles.
	/// This function is locked during callbacks.
	/// @warning This function is locked during callbacks.
	void DestroyParticles();

private:

	friend class b2ParticleSystem;

	b2ParticleSystem* m_system;
	int32 m_firstIndex, m_lastIndex;
	uint32 m_groupFlags;
	float32 m_strength;
	b2ParticleGroup* m_prev;
	b2ParticleGroup* m_next;

	mutable int32 m_timestamp;
	mutable float32 m_mass;
	mutable float32 m_inertia;
	mutable b2Vec2 m_center;
	mutable b2Vec2 m_linearVelocity;
	mutable float32 m_angularVelocity;
	mutable b2Transform m_transform;

	void* m_userData;

	b2ParticleGroup();
	~b2ParticleGroup();
	void UpdateStatistics() const;

};

inline b2ParticleGroup* b2ParticleGroup::GetNext()
{
	return m_next;
}

inline const b2ParticleGroup* b2ParticleGroup::GetNext() const
{
	return m_next;
}

inline b2ParticleSystem* b2ParticleGroup::GetParticleSystem()
{
	return m_system;
}

inline const b2ParticleSystem* b2ParticleGroup::GetParticleSystem() const
{
	return m_system;
}

inline int32 b2ParticleGroup::GetParticleCount() const
{
	return m_lastIndex - m_firstIndex;
}

inline bool b2ParticleGroup::ContainsParticle(int32 index) const
{
	return m_firstIndex <= index && index < m_lastIndex;
}

inline b2ParticleGroup::~b2ParticleGroup()
{
}

inline int32 b2ParticleGroup::GetBufferIndex() const
{
  return m_firstIndex;
}

inline uint32 b2ParticleGroup::GetGroupFlags() const
{
	return m_groupFlags & ~b2_particleGroupInternalMask;
}

inline float32 b2ParticleGroup::GetMass() const
{
	UpdateStatistics();
	return m_mass;
}

inline float32 b2ParticleGroup::GetInertia() const
{
	UpdateStatistics();
	return m_inertia;
}

inline b2Vec2 b2ParticleGroup::GetCenter() const
{
	UpdateStatistics();
	return m_center;
}

inline b2Vec2 b2ParticleGroup::GetLinearVelocity() const
{
	UpdateStatistics();
	return m_linearVelocity;
}

inline float32 b2ParticleGroup::GetAngularVelocity() const
{
	UpdateStatistics();
	return m_angularVelocity;
}

inline const b2Transform& b2ParticleGroup::GetTransform() const
{
	return m_transform;
}

inline const b2Vec2& b2ParticleGroup::GetPosition() const
{
	return m_transform.p;
}

inline float32 b2ParticleGroup::GetAngle() const
{
	return m_transform.q.GetAngle();
}

inline b2Vec2 b2ParticleGroup::GetLinearVelocityFromWorldPoint(
												const b2Vec2& worldPoint) const
{
	UpdateStatistics();
	return m_linearVelocity + b2Cross(m_angularVelocity, worldPoint - m_center);
}

inline void* b2ParticleGroup::GetUserData() const
{
	return m_userData;
}

inline void b2ParticleGroup::SetUserData(void* data)
{
	m_userData = data;
}

inline void b2ParticleGroup::DestroyParticles()
{
	DestroyParticles(false);
}

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
inline void b2ParticleGroupDef::SetPosition(float32 x, float32 y)
{
	position.Set(x, y);
}

inline void b2ParticleGroupDef::SetColor(int32 r, int32 g, int32 b, int32 a)
{
	color.Set((uint8)r, (uint8)g, (uint8)b, (uint8)a);
}
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

// end of ParticleGroup.h

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
#include <Box2D/Collision/Shapes/b2CircleShape.h>
#endif //LIQUIDFUN_EXTERNAL_LANGUAGE_API

b2ParticleGroup::b2ParticleGroup()
{

	m_system = NULL;
	m_firstIndex = 0;
	m_lastIndex = 0;
	m_groupFlags = 0;
	m_strength = 1.0f;
	m_prev = NULL;
	m_next = NULL;

	m_timestamp = -1;
	m_mass = 0;
	m_inertia = 0;
	m_center = b2Vec2_zero;
	m_linearVelocity = b2Vec2_zero;
	m_angularVelocity = 0;
	m_transform.SetIdentity();

	m_userData = NULL;

}

uint32 b2ParticleGroup::GetAllParticleFlags() const
{
	uint32 flags = 0;
	for (int32 i = m_firstIndex; i < m_lastIndex; i++)
	{
		flags |= m_system->m_flagsBuffer.data[i];
	}
	return flags;
}

void b2ParticleGroup::SetGroupFlags(uint32 flags)
{
	b2Assert((flags & b2_particleGroupInternalMask) == 0);
	flags |= m_groupFlags & b2_particleGroupInternalMask;
	m_system->SetGroupFlags(this, flags);
}

void b2ParticleGroup::UpdateStatistics() const
{
	if (m_timestamp != m_system->m_timestamp)
	{
		float32 m = m_system->GetParticleMass();
		m_mass = 0;
		m_center.SetZero();
		m_linearVelocity.SetZero();
		for (int32 i = m_firstIndex; i < m_lastIndex; i++)
		{
			m_mass += m;
			m_center += m * m_system->m_positionBuffer.data[i];
			m_linearVelocity += m * m_system->m_velocityBuffer.data[i];
		}
		if (m_mass > 0)
		{
			m_center *= 1 / m_mass;
			m_linearVelocity *= 1 / m_mass;
		}
		m_inertia = 0;
		m_angularVelocity = 0;
		for (int32 i = m_firstIndex; i < m_lastIndex; i++)
		{
			b2Vec2 p = m_system->m_positionBuffer.data[i] - m_center;
			b2Vec2 v = m_system->m_velocityBuffer.data[i] - m_linearVelocity;
			m_inertia += m * b2Dot(p, p);
			m_angularVelocity += m * b2Cross(p, v);
		}
		if (m_inertia > 0)
		{
			m_angularVelocity *= 1 / m_inertia;
		}
		m_timestamp = m_system->m_timestamp;
	}
}

void b2ParticleGroup::ApplyForce(const b2Vec2& force)
{
	m_system->ApplyForce(m_firstIndex, m_lastIndex, force);
}

void b2ParticleGroup::ApplyLinearImpulse(const b2Vec2& impulse)
{
	m_system->ApplyLinearImpulse(m_firstIndex, m_lastIndex, impulse);
}

void b2ParticleGroup::DestroyParticles(bool callDestructionListener)
{
	b2Assert(m_system->m_world->IsLocked() == false);
	if (m_system->m_world->IsLocked())
	{
		return;
	}

	for (int32 i = m_firstIndex; i < m_lastIndex; i++) {
		m_system->DestroyParticle(i, callDestructionListener);
	}
}

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API
void b2ParticleGroupDef::FreeShapesMemory() {
	if (circleShapes)
	{
		delete[] circleShapes;
		circleShapes = NULL;
	}
	if (ownShapesArray && shapes)
	{
		delete[] shapes;
		shapes = NULL;
		ownShapesArray = false;
	}
}

void b2ParticleGroupDef::SetCircleShapesFromVertexList(void* inBuf,
													   int numShapes,
													   float radius)
{
	float* points = (float*) inBuf;
	// Create circle shapes from vertex list and radius
	b2CircleShape* pCircleShapes = new b2CircleShape[numShapes];
	b2Shape** pShapes = new b2Shape*[numShapes];
	for (int i = 0; i < numShapes; ++i) {
		pCircleShapes[i].m_radius = radius;
		pCircleShapes[i].m_p = b2Vec2(points[i*2], points[i*2+1]);
		pShapes[i] = &pCircleShapes[i];
	}

	// Clean up existing buffers
	FreeShapesMemory();

	// Assign to newly created buffers
	ownShapesArray = true;
	circleShapes = pCircleShapes;
	shapes = pShapes;
	shapeCount = numShapes;
}
#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

// end of ParticleGroup.cpp

template <typename T>
class b2StackQueue
{

public:

	b2StackQueue(b2StackAllocator *allocator, int32 capacity)
	{
		m_allocator = allocator;
		m_buffer = (T*) m_allocator->Allocate(sizeof(T) * capacity);
		m_front = 0;
		m_back = 0;
		m_capacity = capacity;
	}

	~b2StackQueue()
	{
		m_allocator->Free(m_buffer);
	}

	void Push(const T &item)
	{
		if (m_back >= m_capacity)
		{
			for (int32 i = m_front; i < m_back; i++)
			{
				m_buffer[i - m_front] = m_buffer[i];
			}
			m_back -= m_front;
			m_front = 0;
			if (m_back >= m_capacity)
			{
				if (m_capacity > 0)
				{
					m_capacity *= 2;
				}
				else
				{
					m_capacity = 1;
				}
				m_buffer = (T*) m_allocator->Reallocate(m_buffer,
														sizeof(T) * m_capacity);
			}
		}
		m_buffer[m_back] = item;
		m_back++;
	}

	void Pop()
	{
		b2Assert(m_front < m_back);
		m_front++;
	}

	bool Empty() const
	{
		b2Assert(m_front <= m_back);
		return m_front == m_back;
	}

	const T &Front() const
	{
		return m_buffer[m_front];
	}

private:

	b2StackAllocator *m_allocator;
	T* m_buffer;
	int32 m_front;
	int32 m_back;
	int32 m_capacity;

};

// end of StackQueue.h

class b2StackAllocator;
struct b2AABB;

/// A field representing the nearest generator from each point.
class b2VoronoiDiagram
{

public:

	b2VoronoiDiagram(b2StackAllocator* allocator, int32 generatorCapacity);
	~b2VoronoiDiagram();

	/// Add a generator.
	/// @param the position of the generator.
	/// @param a tag used to identify the generator in callback functions.
	/// @param whether to callback for nodes associated with the generator.
	void AddGenerator(const b2Vec2& center, int32 tag, bool necessary);

	/// Generate the Voronoi diagram. It is rasterized with a given interval
	/// in the same range as the necessary generators exist.
	/// @param the interval of the diagram.
	/// @param margin for which the range of the diagram is extended.
	void Generate(float32 radius, float32 margin);

	/// Callback used by GetNodes().
	class NodeCallback
	{
	public:
		virtual ~NodeCallback() {}
		/// Receive tags for generators associated with a node.
		virtual void operator()(int32 a, int32 b, int32 c) = 0;
	};

	/// Enumerate all nodes that contain at least one necessary generator.
	/// @param a callback function object called for each node.
	void GetNodes(NodeCallback& callback) const;

private:

	struct Generator
	{
		b2Vec2 center;
		int32 tag;
		bool necessary;
	};

	struct b2VoronoiDiagramTask
	{
		int32 m_x, m_y, m_i;
		Generator* m_generator;

		b2VoronoiDiagramTask() {}
		b2VoronoiDiagramTask(int32 x, int32 y, int32 i, Generator* g)
		{
			m_x = x;
			m_y = y;
			m_i = i;
			m_generator = g;
		}
	};

	b2StackAllocator *m_allocator;
	Generator* m_generatorBuffer;
	int32 m_generatorCapacity;
	int32 m_generatorCount;
	int32 m_countX, m_countY;
	Generator** m_diagram;

};

// end of VoronoiDiagram.h

b2VoronoiDiagram::b2VoronoiDiagram(
	b2StackAllocator* allocator, int32 generatorCapacity)
{
	m_allocator = allocator;
	m_generatorBuffer =
		(Generator*) allocator->Allocate(
			sizeof(Generator) * generatorCapacity);
	m_generatorCapacity = generatorCapacity;
	m_generatorCount = 0;
	m_countX = 0;
	m_countY = 0;
	m_diagram = NULL;
}

b2VoronoiDiagram::~b2VoronoiDiagram()
{
	if (m_diagram)
	{
		m_allocator->Free(m_diagram);
	}
	m_allocator->Free(m_generatorBuffer);
}

void b2VoronoiDiagram::AddGenerator(
	const b2Vec2& center, int32 tag, bool necessary)
{
	b2Assert(m_generatorCount < m_generatorCapacity);
	Generator& g = m_generatorBuffer[m_generatorCount++];
	g.center = center;
	g.tag = tag;
	g.necessary = necessary;
}

void b2VoronoiDiagram::Generate(float32 radius, float32 margin)
{
	b2Assert(m_diagram == NULL);
	float32 inverseRadius = 1 / radius;
	b2Vec2 lower(+b2_maxFloat, +b2_maxFloat);
	b2Vec2 upper(-b2_maxFloat, -b2_maxFloat);
	for (int32 k = 0; k < m_generatorCount; k++)
	{
		Generator& g = m_generatorBuffer[k];
		if (g.necessary)
		{
			lower = b2Min(lower, g.center);
			upper = b2Max(upper, g.center);
		}
	}
	lower.x -= margin;
	lower.y -= margin;
	upper.x += margin;
	upper.y += margin;
	m_countX = 1 + (int32) (inverseRadius * (upper.x - lower.x));
	m_countY = 1 + (int32) (inverseRadius * (upper.y - lower.y));
	m_diagram = (Generator**)
		m_allocator->Allocate(sizeof(Generator*) * m_countX * m_countY);
	for (int32 i = 0; i < m_countX * m_countY; i++)
	{
		m_diagram[i] = NULL;
	}
	// (4 * m_countX * m_countY) is the queue capacity that is experimentally
	// known to be necessary and sufficient for general particle distributions.
	b2StackQueue<b2VoronoiDiagramTask> queue(
		m_allocator, 4 * m_countX * m_countY);
	for (int32 k = 0; k < m_generatorCount; k++)
	{
		Generator& g = m_generatorBuffer[k];
		g.center = inverseRadius * (g.center - lower);
		int32 x = (int32) g.center.x;
		int32 y = (int32) g.center.y;
		if (x >=0 && y >= 0 && x < m_countX && y < m_countY)
		{
			queue.Push(b2VoronoiDiagramTask(x, y, x + y * m_countX, &g));
		}
	}
	while (!queue.Empty())
	{
		int32 x = queue.Front().m_x;
		int32 y = queue.Front().m_y;
		int32 i = queue.Front().m_i;
		Generator* g = queue.Front().m_generator;
		queue.Pop();
		if (!m_diagram[i])
		{
			m_diagram[i] = g;
			if (x > 0)
			{
				queue.Push(b2VoronoiDiagramTask(x - 1, y, i - 1, g));
			}
			if (y > 0)
			{
				queue.Push(b2VoronoiDiagramTask(x, y - 1, i - m_countX, g));
			}
			if (x < m_countX - 1)
			{
				queue.Push(b2VoronoiDiagramTask(x + 1, y, i + 1, g));
			}
			if (y < m_countY - 1)
			{
				queue.Push(b2VoronoiDiagramTask(x, y + 1, i + m_countX, g));
			}
		}
	}
	for (int32 y = 0; y < m_countY; y++)
	{
		for (int32 x = 0; x < m_countX - 1; x++)
		{
			int32 i = x + y * m_countX;
			Generator* a = m_diagram[i];
			Generator* b = m_diagram[i + 1];
			if (a != b)
			{
				queue.Push(b2VoronoiDiagramTask(x, y, i, b));
				queue.Push(b2VoronoiDiagramTask(x + 1, y, i + 1, a));
			}
		}
	}
	for (int32 y = 0; y < m_countY - 1; y++)
	{
		for (int32 x = 0; x < m_countX; x++)
		{
			int32 i = x + y * m_countX;
			Generator* a = m_diagram[i];
			Generator* b = m_diagram[i + m_countX];
			if (a != b)
			{
				queue.Push(b2VoronoiDiagramTask(x, y, i, b));
				queue.Push(b2VoronoiDiagramTask(x, y + 1, i + m_countX, a));
			}
		}
	}
	while (!queue.Empty())
	{
		const b2VoronoiDiagramTask& task = queue.Front();
		int32 x = task.m_x;
		int32 y = task.m_y;
		int32 i = task.m_i;
		Generator* k = task.m_generator;
		queue.Pop();
		Generator* a = m_diagram[i];
		Generator* b = k;
		if (a != b)
		{
			float32 ax = a->center.x - x;
			float32 ay = a->center.y - y;
			float32 bx = b->center.x - x;
			float32 by = b->center.y - y;
			float32 a2 = ax * ax + ay * ay;
			float32 b2 = bx * bx + by * by;
			if (a2 > b2)
			{
				m_diagram[i] = b;
				if (x > 0)
				{
					queue.Push(b2VoronoiDiagramTask(x - 1, y, i - 1, b));
				}
				if (y > 0)
				{
					queue.Push(b2VoronoiDiagramTask(x, y - 1, i - m_countX, b));
				}
				if (x < m_countX - 1)
				{
					queue.Push(b2VoronoiDiagramTask(x + 1, y, i + 1, b));
				}
				if (y < m_countY - 1)
				{
					queue.Push(b2VoronoiDiagramTask(x, y + 1, i + m_countX, b));
				}
			}
		}
	}
}

void b2VoronoiDiagram::GetNodes(NodeCallback& callback) const
{
	for (int32 y = 0; y < m_countY - 1; y++)
	{
		for (int32 x = 0; x < m_countX - 1; x++)
		{
			int32 i = x + y * m_countX;
			const Generator* a = m_diagram[i];
			const Generator* b = m_diagram[i + 1];
			const Generator* c = m_diagram[i + m_countX];
			const Generator* d = m_diagram[i + 1 + m_countX];
			if (b != c)
			{
				if (a != b && a != c &&
					(a->necessary || b->necessary || c->necessary))
				{
					callback(a->tag, b->tag, c->tag);
				}
				if (d != b && d != c &&
					(b->necessary || d->necessary || c->necessary))
				{
					callback(b->tag, d->tag, c->tag);
				}
			}
		}
	}
}

// end of VoronoiDiagram.cpp

// Define LIQUIDFUN_SIMD_TEST_VS_REFERENCE to run both SIMD and reference
// versions, and assert that the results are identical. This is useful when
// modifying one of the functions, to help verify correctness.
// #define LIQUIDFUN_SIMD_TEST_VS_REFERENCE

// For ease of debugging, remove 'inline'. Then, when an assert hits in the
// test-vs-reference functions, you can easily jump the instruction pointer
// to the top of the function to re-run the test.
#define LIQUIDFUN_SIMD_INLINE inline


static const uint32 xTruncBits = 12;
static const uint32 yTruncBits = 12;
static const uint32 tagBits = 8u * sizeof(uint32);
static const uint32 yOffset = 1u << (yTruncBits - 1u);
static const uint32 yShift = tagBits - yTruncBits;
static const uint32 xShift = tagBits - yTruncBits - xTruncBits;
static const uint32 xScale = 1u << xShift;
static const uint32 xOffset = xScale * (1u << (xTruncBits - 1u));
static const uint32 yMask = ((1u << yTruncBits) - 1u) << yShift;
static const uint32 xMask = ~yMask;
static const uint32 relativeTagRight = 1u << xShift;
static const uint32 relativeTagBottomLeft = (uint32)((1 << yShift) +
                                                    ((~uint32(0)) << xShift));

static const uint32 relativeTagBottomRight = (1u << yShift) + (1u << xShift);

// This functor is passed to std::remove_if in RemoveSpuriousBodyContacts
// to implement the algorithm described there.  It was hoisted out and friended
// as it would not compile with g++ 4.6.3 as a local class.  It is only used in
// that function.
class b2ParticleBodyContactRemovePredicate
{
public:
	b2ParticleBodyContactRemovePredicate(b2ParticleSystem* system,
										 int32* discarded)
		: m_system(system), m_lastIndex(-1), m_currentContacts(0),
		  m_discarded(discarded) {}

	bool operator()(const b2ParticleBodyContact& contact)
	{
		// This implements the selection criteria described in
		// RemoveSpuriousBodyContacts().
		// This functor is iterating through a list of Body contacts per
		// Particle, ordered from near to far.  For up to the maximum number of
		// contacts we allow per point per step, we verify that the contact
		// normal of the Body that genenerated the contact makes physical sense
		// by projecting a point back along that normal and seeing if it
		// intersects the fixture generating the contact.

		if (contact.index != m_lastIndex)
		{
			m_currentContacts = 0;
			m_lastIndex = contact.index;
		}

		if (m_currentContacts++ > k_maxContactsPerPoint)
		{
			++(*m_discarded);
			return true;
		}

		// Project along inverse normal (as returned in the contact) to get the
		// point to check.
		b2Vec2 n = contact.normal;
		// weight is 1-(inv(diameter) * distance)
		n *= m_system->m_particleDiameter * (1 - contact.weight);
		b2Vec2 pos = m_system->m_positionBuffer.data[contact.index] + n;

		// pos is now a point projected back along the contact normal to the
		// contact distance. If the surface makes sense for a contact, pos will
		// now lie on or in the fixture generating
		if (!contact.fixture->TestPoint(pos))
		{
			int32 childCount = contact.fixture->GetShape()->GetChildCount();
			for (int32 childIndex = 0; childIndex < childCount; childIndex++)
			{
				float32 distance;
				b2Vec2 normal;
				contact.fixture->ComputeDistance(pos, &distance, &normal,
																	childIndex);
				if (distance < b2_linearSlop)
				{
					return false;
				}
			}
			++(*m_discarded);
			return true;
		}

		return false;
	}
private:
	// Max number of contacts processed per particle, from nearest to farthest.
	// This must be at least 2 for correctness with concave shapes; 3 was
	// experimentally arrived at as looking reasonable.
	static const int32 k_maxContactsPerPoint = 3;
	const b2ParticleSystem* m_system;
	// Index of last particle processed.
	int32 m_lastIndex;
	// Number of contacts processed for the current particle.
	int32 m_currentContacts;
	// Output the number of discarded contacts.
	int32* m_discarded;
};

namespace {

// Compares the expiration time of two particle indices.
class ExpirationTimeComparator
{
public:
	// Initialize the class with a pointer to an array of particle
	// lifetimes.
	ExpirationTimeComparator(const int32* const expirationTimes) :
		m_expirationTimes(expirationTimes)
	{
	}
	// Empty destructor.
	~ExpirationTimeComparator() { }

	// Compare the lifetime of particleIndexA and particleIndexB
	// returning true if the lifetime of A is greater than B for particles
	// that will expire.  If either particle's lifetime is infinite (<= 0.0f)
	// this function return true if the lifetime of A is lesser than B.
	// When used with std::sort() this results in an array of particle
	// indicies sorted in reverse order by particle lifetime.
	// For example, the set of lifetimes
	// (1.0, 0.7, 0.3, 0.0, -1.0, -2.0)
	// would be sorted as
	// (0.0, -1.0, -2.0, 1.0, 0.7, 0.3)
	bool operator() (const int32 particleIndexA,
					 const int32 particleIndexB) const
	{
		const int32 expirationTimeA = m_expirationTimes[particleIndexA];
		const int32 expirationTimeB = m_expirationTimes[particleIndexB];
		const bool infiniteExpirationTimeA = expirationTimeA <= 0.0f;
		const bool infiniteExpirationTimeB = expirationTimeB <= 0.0f;
		return infiniteExpirationTimeA == infiniteExpirationTimeB ?
			expirationTimeA > expirationTimeB : infiniteExpirationTimeA;
	}

private:
	const int32* m_expirationTimes;
};

// *Very* lightweight pair implementation.
template<typename A, typename B>
struct LightweightPair
{
	A first;
	B second;

	// Compares the value of two FixtureParticle objects returning
	// true if left is a smaller value than right.
	static bool Compare(const LightweightPair& left,
						const LightweightPair& right)
	{
		return left.first < right.first &&
			left.second < right.second;
	}

};

// Allocator for a fixed set of items.
class FixedSetAllocator
{
public:
	// Associate a memory allocator with this object.
	FixedSetAllocator(b2StackAllocator* allocator);
	// Deallocate storage for this class.
	~FixedSetAllocator()
	{
		Clear();
	}

	// Allocate internal storage for this object returning the size.
	int32 Allocate(const int32 itemSize, const int32 count);

	// Deallocate the internal buffer if it's allocated.
	void Clear();

	// Get the number of items in the set.
	int32 GetCount() const { return m_count; }

	// Invalidate an item from the set by index.
	void Invalidate(const int32 itemIndex)
	{
		b2Assert(m_valid);
		m_valid[itemIndex] = 0;
	}

	// Get the buffer which indicates whether items are valid in the set.
	const int8* GetValidBuffer() const { return m_valid; }

protected:
	// Get the internal buffer.
	void* GetBuffer() const { return m_buffer; }
	void* GetBuffer() { return m_buffer; }

	// Reduce the number of items in the set.
	void SetCount(int32 count)
	{
		b2Assert(count <= m_count);
		m_count = count;
	}

private:
	// Set buffer.
	void* m_buffer;
	// Array of size m_count which indicates whether an item is in the
	// corresponding index of m_set (1) or the item is invalid (0).
	int8* m_valid;
	// Number of items in m_set.
	int32 m_count;
	// Allocator used to allocate / free the set.
	b2StackAllocator* m_allocator;
};

// Allocator for a fixed set of objects.
template<typename T>
class TypedFixedSetAllocator : public FixedSetAllocator
{
public:
	// Initialize members of this class.
	TypedFixedSetAllocator(b2StackAllocator* allocator) :
		FixedSetAllocator(allocator) { }

	// Allocate a set of objects, returning the new size of the set.
	int32 Allocate(const int32 numberOfObjects)
	{
		Clear();
		return FixedSetAllocator::Allocate(sizeof(T), numberOfObjects);
	}

	// Get the index of an item in the set if it's valid return an index
	// >= 0, -1 otherwise.
	int32 GetIndex(const T* item) const
	{
		if (item)
		{
			b2Assert(item >= GetBuffer() &&
					 item < GetBuffer() + GetCount());
			const int32 index =
				(int32)(((uint8*)item - (uint8*)GetBuffer()) /
						sizeof(*item));
			if (GetValidBuffer()[index])
			{
				return index;
			}
		}
		return -1;
	}

	// Get the internal buffer.
	const T* GetBuffer() const
	{
		return (const T*)FixedSetAllocator::GetBuffer();
	}
	T* GetBuffer() { return (T*)FixedSetAllocator::GetBuffer(); }
};

// Associates a fixture with a particle index.
typedef LightweightPair<b2Fixture*,int32> FixtureParticle;

// Associates a fixture with a particle index.
typedef LightweightPair<int32,int32> ParticlePair;

}  // namespace

// Set of fixture / particle indices.
class FixtureParticleSet :
	public TypedFixedSetAllocator<FixtureParticle>
{
public:
	// Initialize members of this class.
	FixtureParticleSet(b2StackAllocator* allocator) :
		TypedFixedSetAllocator<FixtureParticle>(allocator) { }


	// Initialize from a set of particle / body contacts for particles
	// that have the b2_fixtureContactListenerParticle flag set.
	void Initialize(const b2ParticleBodyContact * const bodyContacts,
					const int32 numBodyContacts,
					const uint32 * const particleFlagsBuffer);

	// Find the index of a particle / fixture pair in the set or -1
	// if it's not present.
	// NOTE: This was not written as a template function to avoid
	// exposing any dependencies via this header.
	int32 Find(const FixtureParticle& fixtureParticle) const;
};

// Set of particle / particle pairs.
class b2ParticlePairSet : public TypedFixedSetAllocator<ParticlePair>
{
public:
	// Initialize members of this class.
	b2ParticlePairSet(b2StackAllocator* allocator) :
		TypedFixedSetAllocator<ParticlePair>(allocator) { }

	// Initialize from a set of particle contacts.
	void Initialize(const b2ParticleContact * const contacts,
					const int32 numContacts,
					const uint32 * const particleFlagsBuffer);

	// Find the index of a particle pair in the set or -1
	// if it's not present.
	// NOTE: This was not written as a template function to avoid
	// exposing any dependencies via this header.
	int32 Find(const ParticlePair& pair) const;
};

static inline uint32 computeTag(float32 x, float32 y)
{
	return ((uint32)(y + yOffset) << yShift) + (uint32)(xScale * x + xOffset);
}

static inline uint32 computeRelativeTag(uint32 tag, int32 x, int32 y)
{
	return tag + (y << yShift) + (x << xShift);
}

b2ParticleSystem::InsideBoundsEnumerator::InsideBoundsEnumerator(
	uint32 lower, uint32 upper, const Proxy* first, const Proxy* last)
{
	m_xLower = lower & xMask;
	m_xUpper = upper & xMask;
	m_yLower = lower & yMask;
	m_yUpper = upper & yMask;
	m_first = first;
	m_last = last;
	b2Assert(m_first <= m_last);
}

int32 b2ParticleSystem::InsideBoundsEnumerator::GetNext()
{
	while (m_first < m_last)
	{
		uint32 xTag = m_first->tag & xMask;
#if B2_ASSERT_ENABLED
		uint32 yTag = m_first->tag & yMask;
		b2Assert(yTag >= m_yLower);
		b2Assert(yTag <= m_yUpper);
#endif
		if (xTag >= m_xLower && xTag <= m_xUpper)
		{
			return (m_first++)->index;
		}
		m_first++;
	}
	return b2_invalidParticleIndex;
}

b2ParticleSystem::b2ParticleSystem(const b2ParticleSystemDef* def,
								   b2World* world) :
	m_handleAllocator(b2_minParticleSystemBufferCapacity),
	m_stuckParticleBuffer(world->m_blockAllocator),
	m_proxyBuffer(world->m_blockAllocator),
	m_contactBuffer(world->m_blockAllocator),
	m_bodyContactBuffer(world->m_blockAllocator),
	m_pairBuffer(world->m_blockAllocator),
	m_triadBuffer(world->m_blockAllocator)
{
	b2Assert(def);
	m_paused = false;
	m_timestamp = 0;
	m_allParticleFlags = 0;
	m_needsUpdateAllParticleFlags = false;
	m_allGroupFlags = 0;
	m_needsUpdateAllGroupFlags = false;
	m_hasForce = false;
	m_iterationIndex = 0;

	SetStrictContactCheck(def->strictContactCheck);
	SetDensity(def->density);
	SetGravityScale(def->gravityScale);
	SetRadius(def->radius);
	SetMaxParticleCount(def->maxCount);

	m_count = 0;
	m_internalAllocatedCapacity = 0;
	m_forceBuffer = NULL;
	m_weightBuffer = NULL;
	m_staticPressureBuffer = NULL;
	m_accumulationBuffer = NULL;
	m_accumulation2Buffer = NULL;
	m_depthBuffer = NULL;
	m_groupBuffer = NULL;

	m_groupCount = 0;
	m_groupList = NULL;

	b2Assert(def->lifetimeGranularity > 0.0f);
	m_def = *def;

	m_world = world;

	m_stuckThreshold = 0;

	m_timeElapsed = 0;
	m_expirationTimeBufferRequiresSorting = false;

	SetDestructionByAge(m_def.destroyByAge);
}

b2ParticleSystem::~b2ParticleSystem()
{
	while (m_groupList)
	{
		DestroyParticleGroup(m_groupList);
	}

	FreeUserOverridableBuffer(&m_handleIndexBuffer);
	FreeUserOverridableBuffer(&m_flagsBuffer);
	FreeUserOverridableBuffer(&m_lastBodyContactStepBuffer);
	FreeUserOverridableBuffer(&m_bodyContactCountBuffer);
	FreeUserOverridableBuffer(&m_consecutiveContactStepsBuffer);
	FreeUserOverridableBuffer(&m_positionBuffer);
	FreeUserOverridableBuffer(&m_velocityBuffer);
	FreeUserOverridableBuffer(&m_colorBuffer);
	FreeUserOverridableBuffer(&m_userDataBuffer);
	FreeUserOverridableBuffer(&m_expirationTimeBuffer);
	FreeUserOverridableBuffer(&m_indexByExpirationTimeBuffer);
	FreeBuffer(&m_forceBuffer, m_internalAllocatedCapacity);
	FreeBuffer(&m_weightBuffer, m_internalAllocatedCapacity);
	FreeBuffer(&m_staticPressureBuffer, m_internalAllocatedCapacity);
	FreeBuffer(&m_accumulationBuffer, m_internalAllocatedCapacity);
	FreeBuffer(&m_accumulation2Buffer, m_internalAllocatedCapacity);
	FreeBuffer(&m_depthBuffer, m_internalAllocatedCapacity);
	FreeBuffer(&m_groupBuffer, m_internalAllocatedCapacity);
}

template <typename T> void b2ParticleSystem::FreeBuffer(T** b, int capacity)
{
	if (*b == NULL)
		return;

	m_world->m_blockAllocator.Free(*b, sizeof(**b) * capacity);
	*b = NULL;
}

// Free buffer, if it was allocated with b2World's block allocator
template <typename T> void b2ParticleSystem::FreeUserOverridableBuffer(
	UserOverridableBuffer<T>* b)
{
	if (b->userSuppliedCapacity == 0)
	{
		FreeBuffer(&b->data, m_internalAllocatedCapacity);
	}
}

// Reallocate a buffer
template <typename T> T* b2ParticleSystem::ReallocateBuffer(
	T* oldBuffer, int32 oldCapacity, int32 newCapacity)
{
	b2Assert(newCapacity > oldCapacity);
	T* newBuffer = (T*) m_world->m_blockAllocator.Allocate(
		sizeof(T) * newCapacity);
	if (oldBuffer)
	{
		memcpy(newBuffer, oldBuffer, sizeof(T) * oldCapacity);
		m_world->m_blockAllocator.Free(oldBuffer, sizeof(T) * oldCapacity);
	}
	return newBuffer;
}

// Reallocate a buffer
template <typename T> T* b2ParticleSystem::ReallocateBuffer(
	T* buffer, int32 userSuppliedCapacity, int32 oldCapacity,
	int32 newCapacity, bool deferred)
{
	b2Assert(newCapacity > oldCapacity);
	// A 'deferred' buffer is reallocated only if it is not NULL.
	// If 'userSuppliedCapacity' is not zero, buffer is user supplied and must
	// be kept.
	b2Assert(!userSuppliedCapacity || newCapacity <= userSuppliedCapacity);
	if ((!deferred || buffer) && !userSuppliedCapacity)
	{
		buffer = ReallocateBuffer(buffer, oldCapacity, newCapacity);
	}
	return buffer;
}

// Reallocate a buffer
template <typename T> T* b2ParticleSystem::ReallocateBuffer(
	UserOverridableBuffer<T>* buffer, int32 oldCapacity, int32 newCapacity,
	bool deferred)
{
	b2Assert(newCapacity > oldCapacity);
	return ReallocateBuffer(buffer->data, buffer->userSuppliedCapacity,
							oldCapacity, newCapacity, deferred);
}

/// Reallocate the handle / index map and schedule the allocation of a new
/// pool for handle allocation.
void b2ParticleSystem::ReallocateHandleBuffers(int32 newCapacity)
{
	b2Assert(newCapacity > m_internalAllocatedCapacity);
	// Reallocate a new handle / index map buffer, copying old handle pointers
	// is fine since they're kept around.
	m_handleIndexBuffer.data = ReallocateBuffer(
		&m_handleIndexBuffer, m_internalAllocatedCapacity, newCapacity,
		true);
	// Set the size of the next handle allocation.
	m_handleAllocator.SetItemsPerSlab(newCapacity -
									  m_internalAllocatedCapacity);
}

template <typename T> T* b2ParticleSystem::RequestBuffer(T* buffer)
{
	if (!buffer)
	{
		if (m_internalAllocatedCapacity == 0)
		{
			ReallocateInternalAllocatedBuffers(
				b2_minParticleSystemBufferCapacity);
		}
		buffer = (T*) (m_world->m_blockAllocator.Allocate(
						   sizeof(T) * m_internalAllocatedCapacity));
		b2Assert(buffer);
		memset(buffer, 0, sizeof(T) * m_internalAllocatedCapacity);
	}
	return buffer;
}

b2ParticleColor* b2ParticleSystem::GetColorBuffer()
{
	m_colorBuffer.data = RequestBuffer(m_colorBuffer.data);
	return m_colorBuffer.data;
}

void** b2ParticleSystem::GetUserDataBuffer()
{
	m_userDataBuffer.data = RequestBuffer(m_userDataBuffer.data);
	return m_userDataBuffer.data;
}

static int32 LimitCapacity(int32 capacity, int32 maxCount)
{
	return maxCount && capacity > maxCount ? maxCount : capacity;
}

void b2ParticleSystem::ReallocateInternalAllocatedBuffers(int32 capacity)
{
	// Don't increase capacity beyond the smallest user-supplied buffer size.
	capacity = LimitCapacity(capacity, m_def.maxCount);
	capacity = LimitCapacity(capacity, m_flagsBuffer.userSuppliedCapacity);
	capacity = LimitCapacity(capacity, m_positionBuffer.userSuppliedCapacity);
	capacity = LimitCapacity(capacity, m_velocityBuffer.userSuppliedCapacity);
	capacity = LimitCapacity(capacity, m_colorBuffer.userSuppliedCapacity);
	capacity = LimitCapacity(capacity, m_userDataBuffer.userSuppliedCapacity);
	if (m_internalAllocatedCapacity < capacity)
	{
		ReallocateHandleBuffers(capacity);
		m_flagsBuffer.data = ReallocateBuffer(
			&m_flagsBuffer, m_internalAllocatedCapacity, capacity, false);

		// Conditionally defer these as they are optional if the feature is
		// not enabled.
		const bool stuck = m_stuckThreshold > 0;
		m_lastBodyContactStepBuffer.data = ReallocateBuffer(
			&m_lastBodyContactStepBuffer, m_internalAllocatedCapacity,
			capacity, stuck);
		m_bodyContactCountBuffer.data = ReallocateBuffer(
			&m_bodyContactCountBuffer, m_internalAllocatedCapacity, capacity,
			stuck);
		m_consecutiveContactStepsBuffer.data = ReallocateBuffer(
			&m_consecutiveContactStepsBuffer, m_internalAllocatedCapacity,
			capacity, stuck);
		m_positionBuffer.data = ReallocateBuffer(
			&m_positionBuffer, m_internalAllocatedCapacity, capacity, false);
		m_velocityBuffer.data = ReallocateBuffer(
			&m_velocityBuffer, m_internalAllocatedCapacity, capacity, false);
		m_forceBuffer = ReallocateBuffer(
			m_forceBuffer, 0, m_internalAllocatedCapacity, capacity, false);
		m_weightBuffer = ReallocateBuffer(
			m_weightBuffer, 0, m_internalAllocatedCapacity, capacity, false);
		m_staticPressureBuffer = ReallocateBuffer(
			m_staticPressureBuffer, 0, m_internalAllocatedCapacity, capacity,
			true);
		m_accumulationBuffer = ReallocateBuffer(
			m_accumulationBuffer, 0, m_internalAllocatedCapacity, capacity,
			false);
		m_accumulation2Buffer = ReallocateBuffer(
			m_accumulation2Buffer, 0, m_internalAllocatedCapacity, capacity,
			true);
		m_depthBuffer = ReallocateBuffer(
			m_depthBuffer, 0, m_internalAllocatedCapacity, capacity, true);
		m_colorBuffer.data = ReallocateBuffer(
			&m_colorBuffer, m_internalAllocatedCapacity, capacity, true);
		m_groupBuffer = ReallocateBuffer(
			m_groupBuffer, 0, m_internalAllocatedCapacity, capacity, false);
		m_userDataBuffer.data = ReallocateBuffer(
			&m_userDataBuffer, m_internalAllocatedCapacity, capacity, true);
		m_expirationTimeBuffer.data = ReallocateBuffer(
			&m_expirationTimeBuffer, m_internalAllocatedCapacity, capacity,
			true);
		m_indexByExpirationTimeBuffer.data = ReallocateBuffer(
			&m_indexByExpirationTimeBuffer, m_internalAllocatedCapacity,
			capacity, true);
		m_internalAllocatedCapacity = capacity;
	}
}

int32 b2ParticleSystem::CreateParticle(const b2ParticleDef& def)
{
	b2Assert(m_world->IsLocked() == false);
	if (m_world->IsLocked())
	{
		return 0;
	}

	if (m_count >= m_internalAllocatedCapacity)
	{
		// Double the particle capacity.
		int32 capacity =
			m_count ? 2 * m_count : b2_minParticleSystemBufferCapacity;
		ReallocateInternalAllocatedBuffers(capacity);
	}
	if (m_count >= m_internalAllocatedCapacity)
	{
		// If the oldest particle should be destroyed...
		if (m_def.destroyByAge)
		{
			DestroyOldestParticle(0, false);
			// Need to destroy this particle *now* so that it's possible to
			// create a new particle.
			SolveZombie();
		}
		else
		{
			return b2_invalidParticleIndex;
		}
	}
	int32 index = m_count++;
	m_flagsBuffer.data[index] = 0;
	if (m_lastBodyContactStepBuffer.data)
	{
		m_lastBodyContactStepBuffer.data[index] = 0;
	}
	if (m_bodyContactCountBuffer.data)
	{
		m_bodyContactCountBuffer.data[index] = 0;
	}
	if (m_consecutiveContactStepsBuffer.data)
	{
		m_consecutiveContactStepsBuffer.data[index] = 0;
	}
	m_positionBuffer.data[index] = def.position;
	m_velocityBuffer.data[index] = def.velocity;
	m_weightBuffer[index] = 0;
	m_forceBuffer[index] = b2Vec2_zero;
	if (m_staticPressureBuffer)
	{
		m_staticPressureBuffer[index] = 0;
	}
	if (m_depthBuffer)
	{
		m_depthBuffer[index] = 0;
	}
	if (m_colorBuffer.data || !def.color.IsZero())
	{
		m_colorBuffer.data = RequestBuffer(m_colorBuffer.data);
		m_colorBuffer.data[index] = def.color;
	}
	if (m_userDataBuffer.data || def.userData)
	{
		m_userDataBuffer.data= RequestBuffer(m_userDataBuffer.data);
		m_userDataBuffer.data[index] = def.userData;
	}
	if (m_handleIndexBuffer.data)
	{
		m_handleIndexBuffer.data[index] = NULL;
	}
	Proxy& proxy = m_proxyBuffer.Append();

	// If particle lifetimes are enabled or the lifetime is set in the particle
	// definition, initialize the lifetime.
	const bool finiteLifetime = def.lifetime > 0;
	if (m_expirationTimeBuffer.data || finiteLifetime)
	{
		SetParticleLifetime(index, finiteLifetime ? def.lifetime :
								ExpirationTimeToLifetime(
									-GetQuantizedTimeElapsed()));
		// Add a reference to the newly added particle to the end of the
		// queue.
		m_indexByExpirationTimeBuffer.data[index] = index;
	}

	proxy.index = index;
	b2ParticleGroup* group = def.group;
	m_groupBuffer[index] = group;
	if (group)
	{
		if (group->m_firstIndex < group->m_lastIndex)
		{
			// Move particles in the group just before the new particle.
			RotateBuffer(group->m_firstIndex, group->m_lastIndex, index);
			b2Assert(group->m_lastIndex == index);
			// Update the index range of the group to contain the new particle.
			group->m_lastIndex = index + 1;
		}
		else
		{
			// If the group is empty, reset the index range to contain only the
			// new particle.
			group->m_firstIndex = index;
			group->m_lastIndex = index + 1;
		}
	}
	SetParticleFlags(index, def.flags);
	return index;
}

/// Retrieve a handle to the particle at the specified index.
const b2ParticleHandle* b2ParticleSystem::GetParticleHandleFromIndex(
	const int32 index)
{
	b2Assert(index >= 0 && index < GetParticleCount() &&
			 index != b2_invalidParticleIndex);
	m_handleIndexBuffer.data = RequestBuffer(m_handleIndexBuffer.data);
	b2ParticleHandle* handle = m_handleIndexBuffer.data[index];
	if (handle)
	{
		return handle;
	}
	// Create a handle.
	handle = m_handleAllocator.Allocate();
	b2Assert(handle);
	handle->SetIndex(index);
	m_handleIndexBuffer.data[index] = handle;
	return handle;
}


void b2ParticleSystem::DestroyParticle(
	int32 index, bool callDestructionListener)
{
	uint32 flags = b2_zombieParticle;
	if (callDestructionListener)
	{
		flags |= b2_destructionListenerParticle;
	}
	SetParticleFlags(index, m_flagsBuffer.data[index] | flags);
}

void b2ParticleSystem::DestroyOldestParticle(
	const int32 index, const bool callDestructionListener)
{
	const int32 particleCount = GetParticleCount();
	b2Assert(index >= 0 && index < particleCount);
	// Make sure particle lifetime tracking is enabled.
	b2Assert(m_indexByExpirationTimeBuffer.data);
	// Destroy the oldest particle (preferring to destroy finite
	// lifetime particles first) to free a slot in the buffer.
	const int32 oldestFiniteLifetimeParticle =
		m_indexByExpirationTimeBuffer.data[particleCount - (index + 1)];
	const int32 oldestInfiniteLifetimeParticle =
		m_indexByExpirationTimeBuffer.data[index];
	DestroyParticle(
		m_expirationTimeBuffer.data[oldestFiniteLifetimeParticle] > 0.0f ?
			oldestFiniteLifetimeParticle : oldestInfiniteLifetimeParticle,
		callDestructionListener);
}

int32 b2ParticleSystem::DestroyParticlesInShape(
	const b2Shape& shape, const b2Transform& xf,
	bool callDestructionListener)
{
	b2Assert(m_world->IsLocked() == false);
	if (m_world->IsLocked())
	{
		return 0;
	}

	class DestroyParticlesInShapeCallback : public b2QueryCallback
	{
	public:
		DestroyParticlesInShapeCallback(
			b2ParticleSystem* system, const b2Shape& shape,
			const b2Transform& xf, bool callDestructionListener)
		{
			m_system = system;
			m_shape = &shape;
			m_xf = xf;
			m_callDestructionListener = callDestructionListener;
			m_destroyed = 0;
		}

		bool ReportFixture(b2Fixture* fixture)
		{
			B2_NOT_USED(fixture);
			return false;
		}

		bool ReportParticle(const b2ParticleSystem* particleSystem, int32 index)
		{
			if (particleSystem != m_system)
				return false;

			b2Assert(index >=0 && index < m_system->m_count);
			if (m_shape->TestPoint(m_xf,
								   m_system->m_positionBuffer.data[index]))
			{
				m_system->DestroyParticle(index, m_callDestructionListener);
				m_destroyed++;
			}
			return true;
		}

		int32 Destroyed() { return m_destroyed; }

	private:
		b2ParticleSystem* m_system;
		const b2Shape* m_shape;
		b2Transform m_xf;
		bool m_callDestructionListener;
		int32 m_destroyed;
	} callback(this, shape, xf, callDestructionListener);
	b2AABB aabb;
	shape.ComputeAABB(&aabb, xf, 0);
	m_world->QueryAABB(&callback, aabb);
	return callback.Destroyed();
}

int32 b2ParticleSystem::CreateParticleForGroup(
	const b2ParticleGroupDef& groupDef, const b2Transform& xf, const b2Vec2& p)
{
	b2ParticleDef particleDef;
	particleDef.flags = groupDef.flags;
	particleDef.position = b2Mul(xf, p);
	particleDef.velocity =
		groupDef.linearVelocity +
		b2Cross(groupDef.angularVelocity,
				particleDef.position - groupDef.position);
	particleDef.color = groupDef.color;
	particleDef.lifetime = groupDef.lifetime;
	particleDef.userData = groupDef.userData;
	return CreateParticle(particleDef);
}

void b2ParticleSystem::CreateParticlesStrokeShapeForGroup(
	const b2Shape *shape,
	const b2ParticleGroupDef& groupDef, const b2Transform& xf)
{
	float32 stride = groupDef.stride;
	if (stride == 0)
	{
		stride = GetParticleStride();
	}
	float32 positionOnEdge = 0;
	int32 childCount = shape->GetChildCount();
	for (int32 childIndex = 0; childIndex < childCount; childIndex++)
	{
		b2EdgeShape edge;
		if (shape->GetType() == b2Shape::e_edge)
		{
			edge = *(b2EdgeShape*) shape;
		}
		else
		{
			b2Assert(shape->GetType() == b2Shape::e_chain);
			((b2ChainShape*) shape)->GetChildEdge(&edge, childIndex);
		}
		b2Vec2 d = edge.m_vertex2 - edge.m_vertex1;
		float32 edgeLength = d.Length();
		while (positionOnEdge < edgeLength)
		{
			b2Vec2 p = edge.m_vertex1 + positionOnEdge / edgeLength * d;
			CreateParticleForGroup(groupDef, xf, p);
			positionOnEdge += stride;
		}
		positionOnEdge -= edgeLength;
	}
}

void b2ParticleSystem::CreateParticlesFillShapeForGroup(
	const b2Shape *shape,
	const b2ParticleGroupDef& groupDef, const b2Transform& xf)
{
	float32 stride = groupDef.stride;
	if (stride == 0)
	{
		stride = GetParticleStride();
	}
	b2Transform identity;
	identity.SetIdentity();
	b2AABB aabb;
	b2Assert(shape->GetChildCount() == 1);
	shape->ComputeAABB(&aabb, identity, 0);
	for (float32 y = floorf(aabb.lowerBound.y / stride) * stride;
		y < aabb.upperBound.y; y += stride)
	{
		for (float32 x = floorf(aabb.lowerBound.x / stride) * stride;
			x < aabb.upperBound.x; x += stride)
		{
			b2Vec2 p(x, y);
			if (shape->TestPoint(identity, p))
			{
				CreateParticleForGroup(groupDef, xf, p);
			}
		}
	}
}

void b2ParticleSystem::CreateParticlesWithShapeForGroup(
	const b2Shape* shape,
	const b2ParticleGroupDef& groupDef, const b2Transform& xf)
{
	switch (shape->GetType()) {
	case b2Shape::e_edge:
	case b2Shape::e_chain:
		CreateParticlesStrokeShapeForGroup(shape, groupDef, xf);
		break;
	case b2Shape::e_polygon:
	case b2Shape::e_circle:
		CreateParticlesFillShapeForGroup(shape, groupDef, xf);
		break;
	default:
		b2Assert(false);
		break;
	}
}

void b2ParticleSystem::CreateParticlesWithShapesForGroup(
	const b2Shape* const* shapes, int32 shapeCount,
	const b2ParticleGroupDef& groupDef, const b2Transform& xf)
{
	class CompositeShape : public b2Shape
	{
	public:
		CompositeShape(const b2Shape* const* shapes, int32 shapeCount)
		{
			m_shapes = shapes;
			m_shapeCount = shapeCount;
		}
		b2Shape* Clone(b2BlockAllocator* allocator) const
		{
			b2Assert(false);
			B2_NOT_USED(allocator);
			return NULL;
		}
		int32 GetChildCount() const
		{
			return 1;
		}
		bool TestPoint(const b2Transform& xf, const b2Vec2& p) const
		{
			for (int32 i = 0; i < m_shapeCount; i++)
			{
				if (m_shapes[i]->TestPoint(xf, p))
				{
					return true;
				}
			}
			return false;
		}
		void ComputeDistance(const b2Transform& xf, const b2Vec2& p,
					float32* distance, b2Vec2* normal, int32 childIndex) const
		{
			b2Assert(false);
			B2_NOT_USED(xf);
			B2_NOT_USED(p);
			B2_NOT_USED(distance);
			B2_NOT_USED(normal);
			B2_NOT_USED(childIndex);
		}
		bool RayCast(b2RayCastOutput* output, const b2RayCastInput& input,
						const b2Transform& transform, int32 childIndex) const
		{
			b2Assert(false);
			B2_NOT_USED(output);
			B2_NOT_USED(input);
			B2_NOT_USED(transform);
			B2_NOT_USED(childIndex);
			return false;
		}
		void ComputeAABB(
				b2AABB* aabb, const b2Transform& xf, int32 childIndex) const
		{
			B2_NOT_USED(childIndex);
			aabb->lowerBound.x = +FLT_MAX;
			aabb->lowerBound.y = +FLT_MAX;
			aabb->upperBound.x = -FLT_MAX;
			aabb->upperBound.y = -FLT_MAX;
			b2Assert(childIndex == 0);
			for (int32 i = 0; i < m_shapeCount; i++)
			{
				int32 childCount = m_shapes[i]->GetChildCount();
				for (int32 j = 0; j < childCount; j++)
				{
					b2AABB subaabb;
					m_shapes[i]->ComputeAABB(&subaabb, xf, j);
					aabb->Combine(subaabb);
				}
			}
		}
		void ComputeMass(b2MassData* massData, float32 density) const
		{
			b2Assert(false);
			B2_NOT_USED(massData);
			B2_NOT_USED(density);
		}
	private:
		const b2Shape* const* m_shapes;
		int32 m_shapeCount;
	} compositeShape(shapes, shapeCount);
	CreateParticlesFillShapeForGroup(&compositeShape, groupDef, xf);
}

b2ParticleGroup* b2ParticleSystem::CreateParticleGroup(
	const b2ParticleGroupDef& groupDef)
{
	b2Assert(m_world->IsLocked() == false);
	if (m_world->IsLocked())
	{
		return 0;
	}

	b2Transform transform;
	transform.Set(groupDef.position, groupDef.angle);
	int32 firstIndex = m_count;
	if (groupDef.shape)
	{
		CreateParticlesWithShapeForGroup(groupDef.shape, groupDef, transform);
	}
	if (groupDef.shapes)
	{
		CreateParticlesWithShapesForGroup(
					groupDef.shapes, groupDef.shapeCount, groupDef, transform);
	}
	if (groupDef.particleCount)
	{
		b2Assert(groupDef.positionData);
		for (int32 i = 0; i < groupDef.particleCount; i++)
		{
			b2Vec2 p = groupDef.positionData[i];
			CreateParticleForGroup(groupDef, transform, p);
		}
	}
	int32 lastIndex = m_count;

	void* mem = m_world->m_blockAllocator.Allocate(sizeof(b2ParticleGroup));
	b2ParticleGroup* group = new (mem) b2ParticleGroup();
	group->m_system = this;
	group->m_firstIndex = firstIndex;
	group->m_lastIndex = lastIndex;
	group->m_strength = groupDef.strength;
	group->m_userData = groupDef.userData;
	group->m_transform = transform;
	group->m_prev = NULL;
	group->m_next = m_groupList;
	if (m_groupList)
	{
		m_groupList->m_prev = group;
	}
	m_groupList = group;
	++m_groupCount;
	for (int32 i = firstIndex; i < lastIndex; i++)
	{
		m_groupBuffer[i] = group;
	}
	SetGroupFlags(group, groupDef.groupFlags);

	// Create pairs and triads between particles in the group.
	ConnectionFilter filter;
	UpdateContacts(true);
	UpdatePairsAndTriads(firstIndex, lastIndex, filter);

	if (groupDef.group)
	{
		JoinParticleGroups(groupDef.group, group);
		group = groupDef.group;
	}

	return group;
}

void b2ParticleSystem::JoinParticleGroups(b2ParticleGroup* groupA,
										  b2ParticleGroup* groupB)
{
	b2Assert(m_world->IsLocked() == false);
	if (m_world->IsLocked())
	{
		return;
	}

	b2Assert(groupA != groupB);
	RotateBuffer(groupB->m_firstIndex, groupB->m_lastIndex, m_count);
	b2Assert(groupB->m_lastIndex == m_count);
	RotateBuffer(groupA->m_firstIndex, groupA->m_lastIndex,
				 groupB->m_firstIndex);
	b2Assert(groupA->m_lastIndex == groupB->m_firstIndex);

	// Create pairs and triads connecting groupA and groupB.
	class JoinParticleGroupsFilter : public ConnectionFilter
	{
		bool ShouldCreatePair(int32 a, int32 b) const
		{
			return
				(a < m_threshold && m_threshold <= b) ||
				(b < m_threshold && m_threshold <= a);
		}
		bool ShouldCreateTriad(int32 a, int32 b, int32 c) const
		{
			return
				(a < m_threshold || b < m_threshold || c < m_threshold) &&
				(m_threshold <= a || m_threshold <= b || m_threshold <= c);
		}
		int32 m_threshold;
	public:
		JoinParticleGroupsFilter(int32 threshold)
		{
			m_threshold = threshold;
		}
	} filter(groupB->m_firstIndex);
	UpdateContacts(true);
	UpdatePairsAndTriads(groupA->m_firstIndex, groupB->m_lastIndex, filter);

	for (int32 i = groupB->m_firstIndex; i < groupB->m_lastIndex; i++)
	{
		m_groupBuffer[i] = groupA;
	}
	uint32 groupFlags = groupA->m_groupFlags | groupB->m_groupFlags;
	SetGroupFlags(groupA, groupFlags);
	groupA->m_lastIndex = groupB->m_lastIndex;
	groupB->m_firstIndex = groupB->m_lastIndex;
	DestroyParticleGroup(groupB);
}

void b2ParticleSystem::SplitParticleGroup(b2ParticleGroup* group)
{
	UpdateContacts(true);
	int32 particleCount = group->GetParticleCount();
	// We create several linked lists. Each list represents a set of connected
	// particles.
	ParticleListNode* nodeBuffer =
		(ParticleListNode*) m_world->m_stackAllocator.Allocate(
									sizeof(ParticleListNode) * particleCount);
	InitializeParticleLists(group, nodeBuffer);
	MergeParticleListsInContact(group, nodeBuffer);
	ParticleListNode* survivingList =
									FindLongestParticleList(group, nodeBuffer);
	MergeZombieParticleListNodes(group, nodeBuffer, survivingList);
	CreateParticleGroupsFromParticleList(group, nodeBuffer, survivingList);
	UpdatePairsAndTriadsWithParticleList(group, nodeBuffer);
	m_world->m_stackAllocator.Free(nodeBuffer);
}

void b2ParticleSystem::InitializeParticleLists(
	const b2ParticleGroup* group, ParticleListNode* nodeBuffer)
{
	int32 bufferIndex = group->GetBufferIndex();
	int32 particleCount = group->GetParticleCount();
	for (int32 i = 0; i < particleCount; i++)
	{
		ParticleListNode* node = &nodeBuffer[i];
		node->list = node;
		node->next = NULL;
		node->count = 1;
		node->index = i + bufferIndex;
	}
}

void b2ParticleSystem::MergeParticleListsInContact(
	const b2ParticleGroup* group, ParticleListNode* nodeBuffer) const
{
	int32 bufferIndex = group->GetBufferIndex();
	for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
	{
		const b2ParticleContact& contact = m_contactBuffer[k];
		int32 a = contact.GetIndexA();
		int32 b = contact.GetIndexB();
		if (!group->ContainsParticle(a) || !group->ContainsParticle(b)) {
			continue;
		}
		ParticleListNode* listA = nodeBuffer[a - bufferIndex].list;
		ParticleListNode* listB = nodeBuffer[b - bufferIndex].list;
		if (listA == listB) {
			continue;
		}
		// To minimize the cost of insertion, make sure listA is longer than
		// listB.
		if (listA->count < listB->count)
		{
			b2Swap(listA, listB);
		}
		b2Assert(listA->count >= listB->count);
		MergeParticleLists(listA, listB);
	}
}

void b2ParticleSystem::MergeParticleLists(
	ParticleListNode* listA, ParticleListNode* listB)
{
	// Insert listB between index 0 and 1 of listA
	// Example:
	//     listA => a1 => a2 => a3 => NULL
	//     listB => b1 => b2 => NULL
	// to
	//     listA => listB => b1 => b2 => a1 => a2 => a3 => NULL
	b2Assert(listA != listB);
	for (ParticleListNode* b = listB;;)
	{
		b->list = listA;
		ParticleListNode* nextB = b->next;
		if (nextB)
		{
			b = nextB;
		}
		else
		{
			b->next = listA->next;
			break;
		}
	}
	listA->next = listB;
	listA->count += listB->count;
	listB->count = 0;
}

b2ParticleSystem::ParticleListNode* b2ParticleSystem::FindLongestParticleList(
	const b2ParticleGroup* group, ParticleListNode* nodeBuffer)
{
	int32 particleCount = group->GetParticleCount();
	ParticleListNode* result = nodeBuffer;
	for (int32 i = 0; i < particleCount; i++)
	{
		ParticleListNode* node = &nodeBuffer[i];
		if (result->count < node->count)
		{
			result = node;
		}
	}
	return result;
}

void b2ParticleSystem::MergeZombieParticleListNodes(
	const b2ParticleGroup* group, ParticleListNode* nodeBuffer,
	ParticleListNode* survivingList) const
{
	int32 particleCount = group->GetParticleCount();
	for (int32 i = 0; i < particleCount; i++)
	{
		ParticleListNode* node = &nodeBuffer[i];
		if (node != survivingList &&
			(m_flagsBuffer.data[node->index] & b2_zombieParticle))
		{
			MergeParticleListAndNode(survivingList, node);
		}
	}
}

void b2ParticleSystem::MergeParticleListAndNode(
	ParticleListNode* list, ParticleListNode* node)
{
	// Insert node between index 0 and 1 of list
	// Example:
	//     list => a1 => a2 => a3 => NULL
	//     node => NULL
	// to
	//     list => node => a1 => a2 => a3 => NULL
	b2Assert(node != list);
	b2Assert(node->list == node);
	b2Assert(node->count == 1);
	node->list = list;
	node->next = list->next;
	list->next = node;
	list->count++;
	node->count = 0;
}

void b2ParticleSystem::CreateParticleGroupsFromParticleList(
	const b2ParticleGroup* group, ParticleListNode* nodeBuffer,
	const ParticleListNode* survivingList)
{
	int32 particleCount = group->GetParticleCount();
	b2ParticleGroupDef def;
	def.groupFlags = group->GetGroupFlags();
	def.userData = group->GetUserData();
	for (int32 i = 0; i < particleCount; i++)
	{
		ParticleListNode* list = &nodeBuffer[i];
		if (!list->count || list == survivingList)
		{
			continue;
		}
		b2Assert(list->list == list);
		b2ParticleGroup* newGroup = CreateParticleGroup(def);
		for (ParticleListNode* node = list; node; node = node->next)
		{
			int32 oldIndex = node->index;
			b2Assert(!(m_flagsBuffer.data[oldIndex] & b2_zombieParticle));
			int32 newIndex = CloneParticle(oldIndex, newGroup);
			m_flagsBuffer.data[oldIndex] |= b2_zombieParticle;
			node->index = newIndex;
		}
	}
}

void b2ParticleSystem::UpdatePairsAndTriadsWithParticleList(
	const b2ParticleGroup* group, const ParticleListNode* nodeBuffer)
{
	int32 bufferIndex = group->GetBufferIndex();
	// Update indices in pairs and triads. If an index belongs to the group,
	// replace it with the corresponding value in nodeBuffer.
	// Note that nodeBuffer is allocated only for the group and the index should
	// be shifted by bufferIndex.
	for (int32 k = 0; k < m_pairBuffer.GetCount(); k++)
	{
		b2ParticlePair& pair = m_pairBuffer[k];
		int32 a = pair.indexA;
		int32 b = pair.indexB;
		if (group->ContainsParticle(a))
		{
			pair.indexA = nodeBuffer[a - bufferIndex].index;
		}
		if (group->ContainsParticle(b))
		{
			pair.indexB = nodeBuffer[b - bufferIndex].index;
		}
	}
	for (int32 k = 0; k < m_triadBuffer.GetCount(); k++)
	{
		b2ParticleTriad& triad = m_triadBuffer[k];
		int32 a = triad.indexA;
		int32 b = triad.indexB;
		int32 c = triad.indexC;
		if (group->ContainsParticle(a))
		{
			triad.indexA = nodeBuffer[a - bufferIndex].index;
		}
		if (group->ContainsParticle(b))
		{
			triad.indexB = nodeBuffer[b - bufferIndex].index;
		}
		if (group->ContainsParticle(c))
		{
			triad.indexC = nodeBuffer[c - bufferIndex].index;
		}
	}
}

int32 b2ParticleSystem::CloneParticle(int32 oldIndex, b2ParticleGroup* group)
{
	b2ParticleDef def;
	def.flags = m_flagsBuffer.data[oldIndex];
	def.position = m_positionBuffer.data[oldIndex];
	def.velocity = m_velocityBuffer.data[oldIndex];
	if (m_colorBuffer.data)
	{
		def.color = m_colorBuffer.data[oldIndex];
	}
	if (m_userDataBuffer.data)
	{
		def.userData = m_userDataBuffer.data[oldIndex];
	}
	def.group = group;
	int32 newIndex = CreateParticle(def);
	if (m_handleIndexBuffer.data)
	{
		b2ParticleHandle* handle = m_handleIndexBuffer.data[oldIndex];
		if (handle) handle->SetIndex(newIndex);
		m_handleIndexBuffer.data[newIndex] = handle;
		m_handleIndexBuffer.data[oldIndex] = NULL;
	}
	if (m_lastBodyContactStepBuffer.data)
	{
		m_lastBodyContactStepBuffer.data[newIndex] =
			m_lastBodyContactStepBuffer.data[oldIndex];
	}
	if (m_bodyContactCountBuffer.data)
	{
		m_bodyContactCountBuffer.data[newIndex] =
			m_bodyContactCountBuffer.data[oldIndex];
	}
	if (m_consecutiveContactStepsBuffer.data)
	{
		m_consecutiveContactStepsBuffer.data[newIndex] =
			m_consecutiveContactStepsBuffer.data[oldIndex];
	}
	if (m_hasForce)
	{
		m_forceBuffer[newIndex] = m_forceBuffer[oldIndex];
	}
	if (m_staticPressureBuffer)
	{
		m_staticPressureBuffer[newIndex] = m_staticPressureBuffer[oldIndex];
	}
	if (m_depthBuffer)
	{
		m_depthBuffer[newIndex] = m_depthBuffer[oldIndex];
	}
	if (m_expirationTimeBuffer.data)
	{
		m_expirationTimeBuffer.data[newIndex] =
			m_expirationTimeBuffer.data[oldIndex];
	}
	return newIndex;
}

void b2ParticleSystem::UpdatePairsAndTriadsWithReactiveParticles()
{
	class ReactiveFilter : public ConnectionFilter
	{
		bool IsNecessary(int32 index) const
		{
			return (m_flagsBuffer[index] & b2_reactiveParticle) != 0;
		}
		const uint32* m_flagsBuffer;
	public:
		ReactiveFilter(uint32* flagsBuffer)
		{
			m_flagsBuffer = flagsBuffer;
		}
	} filter(m_flagsBuffer.data);
	UpdatePairsAndTriads(0, m_count, filter);

	for (int32 i = 0; i < m_count; i++)
	{
		m_flagsBuffer.data[i] &= ~b2_reactiveParticle;
	}
	m_allParticleFlags &= ~b2_reactiveParticle;
}

static bool ParticleCanBeConnected(
	uint32 flags, b2ParticleGroup* group)
{
	return
		(flags & (b2_wallParticle | b2_springParticle | b2_elasticParticle)) ||
		(group && group->GetGroupFlags() & b2_rigidParticleGroup);
}

void b2ParticleSystem::UpdatePairsAndTriads(
	int32 firstIndex, int32 lastIndex, const ConnectionFilter& filter)
{
	// Create pairs or triads.
	// All particles in each pair/triad should satisfy the following:
	// * firstIndex <= index < lastIndex
	// * don't have b2_zombieParticle
	// * ParticleCanBeConnected returns true
	// * ShouldCreatePair/ShouldCreateTriad returns true
	// Any particles in each pair/triad should satisfy the following:
	// * filter.IsNeeded returns true
	// * have one of k_pairFlags/k_triadsFlags
	b2Assert(firstIndex <= lastIndex);
	uint32 particleFlags = 0;
	for (int32 i = firstIndex; i < lastIndex; i++)
	{
		particleFlags |= m_flagsBuffer.data[i];
	}
	if (particleFlags & k_pairFlags)
	{
		for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
		{
			const b2ParticleContact& contact = m_contactBuffer[k];
			int32 a = contact.GetIndexA();
			int32 b = contact.GetIndexB();
			uint32 af = m_flagsBuffer.data[a];
			uint32 bf = m_flagsBuffer.data[b];
			b2ParticleGroup* groupA = m_groupBuffer[a];
			b2ParticleGroup* groupB = m_groupBuffer[b];
			if (a >= firstIndex && a < lastIndex &&
				b >= firstIndex && b < lastIndex &&
				!((af | bf) & b2_zombieParticle) &&
				((af | bf) & k_pairFlags) &&
				(filter.IsNecessary(a) || filter.IsNecessary(b)) &&
				ParticleCanBeConnected(af, groupA) &&
				ParticleCanBeConnected(bf, groupB) &&
				filter.ShouldCreatePair(a, b))
			{
				b2ParticlePair& pair = m_pairBuffer.Append();
				pair.indexA = a;
				pair.indexB = b;
				pair.flags = contact.GetFlags();
				pair.strength = b2Min(
					groupA ? groupA->m_strength : 1,
					groupB ? groupB->m_strength : 1);
				pair.distance = b2Distance(m_positionBuffer.data[a],
										   m_positionBuffer.data[b]);
			}
		}
		std::stable_sort(
			m_pairBuffer.Begin(), m_pairBuffer.End(), ComparePairIndices);
		m_pairBuffer.Unique(MatchPairIndices);
	}
	if (particleFlags & k_triadFlags)
	{
		b2VoronoiDiagram diagram(
			&m_world->m_stackAllocator, lastIndex - firstIndex);
		for (int32 i = firstIndex; i < lastIndex; i++)
		{
			uint32 flags = m_flagsBuffer.data[i];
			b2ParticleGroup* group = m_groupBuffer[i];
			if (!(flags & b2_zombieParticle) &&
				ParticleCanBeConnected(flags, group))
			{
				diagram.AddGenerator(
					m_positionBuffer.data[i], i, filter.IsNecessary(i));
			}
		}
		float32 stride = GetParticleStride();
		diagram.Generate(stride / 2, stride * 2);
		class UpdateTriadsCallback : public b2VoronoiDiagram::NodeCallback
		{
			void operator()(int32 a, int32 b, int32 c)
			{
				uint32 af = m_system->m_flagsBuffer.data[a];
				uint32 bf = m_system->m_flagsBuffer.data[b];
				uint32 cf = m_system->m_flagsBuffer.data[c];
				if (((af | bf | cf) & k_triadFlags) &&
					m_filter->ShouldCreateTriad(a, b, c))
				{
					const b2Vec2& pa = m_system->m_positionBuffer.data[a];
					const b2Vec2& pb = m_system->m_positionBuffer.data[b];
					const b2Vec2& pc = m_system->m_positionBuffer.data[c];
					b2Vec2 dab = pa - pb;
					b2Vec2 dbc = pb - pc;
					b2Vec2 dca = pc - pa;
					float32 maxDistanceSquared = b2_maxTriadDistanceSquared *
												 m_system->m_squaredDiameter;
					if (b2Dot(dab, dab) > maxDistanceSquared ||
						b2Dot(dbc, dbc) > maxDistanceSquared ||
						b2Dot(dca, dca) > maxDistanceSquared)
					{
						return;
					}
					b2ParticleGroup* groupA = m_system->m_groupBuffer[a];
					b2ParticleGroup* groupB = m_system->m_groupBuffer[b];
					b2ParticleGroup* groupC = m_system->m_groupBuffer[c];
					b2ParticleTriad& triad = m_system->m_triadBuffer.Append();
					triad.indexA = a;
					triad.indexB = b;
					triad.indexC = c;
					triad.flags = af | bf | cf;
					triad.strength = b2Min(b2Min(
						groupA ? groupA->m_strength : 1,
						groupB ? groupB->m_strength : 1),
						groupC ? groupC->m_strength : 1);
					b2Vec2 midPoint = (float32) 1 / 3 * (pa + pb + pc);
					triad.pa = pa - midPoint;
					triad.pb = pb - midPoint;
					triad.pc = pc - midPoint;
					triad.ka = -b2Dot(dca, dab);
					triad.kb = -b2Dot(dab, dbc);
					triad.kc = -b2Dot(dbc, dca);
					triad.s = b2Cross(pa, pb) + b2Cross(pb, pc) + b2Cross(pc, pa);
				}
			}
			b2ParticleSystem* m_system;
			const ConnectionFilter* m_filter;
		public:
			UpdateTriadsCallback(
				b2ParticleSystem* system, const ConnectionFilter* filter)
			{
				m_system = system;
				m_filter = filter;
			}
		} callback(this, &filter);
		diagram.GetNodes(callback);
		std::stable_sort(
			m_triadBuffer.Begin(), m_triadBuffer.End(), CompareTriadIndices);
		m_triadBuffer.Unique(MatchTriadIndices);
	}
}

bool b2ParticleSystem::ComparePairIndices(
							const b2ParticlePair& a, const b2ParticlePair& b)
{
	int32 diffA = a.indexA - b.indexA;
	if (diffA != 0) return diffA < 0;
	return a.indexB < b.indexB;
}

bool b2ParticleSystem::MatchPairIndices(
							const b2ParticlePair& a, const b2ParticlePair& b)
{
	return a.indexA == b.indexA && a.indexB == b.indexB;
}

bool b2ParticleSystem::CompareTriadIndices(
							const b2ParticleTriad& a, const b2ParticleTriad& b)
{
	int32 diffA = a.indexA - b.indexA;
	if (diffA != 0) return diffA < 0;
	int32 diffB = a.indexB - b.indexB;
	if (diffB != 0) return diffB < 0;
	return a.indexC < b.indexC;
}

bool b2ParticleSystem::MatchTriadIndices(
							const b2ParticleTriad& a, const b2ParticleTriad& b)
{
	return a.indexA == b.indexA && a.indexB == b.indexB && a.indexC == b.indexC;
}

// Only called from SolveZombie() or JoinParticleGroups().
void b2ParticleSystem::DestroyParticleGroup(b2ParticleGroup* group)
{
	b2Assert(m_groupCount > 0);
	b2Assert(group);

	if (m_world->m_destructionListener)
	{
		m_world->m_destructionListener->SayGoodbye(group);
	}

	SetGroupFlags(group, 0);
	for (int32 i = group->m_firstIndex; i < group->m_lastIndex; i++)
	{
		m_groupBuffer[i] = NULL;
	}

	if (group->m_prev)
	{
		group->m_prev->m_next = group->m_next;
	}
	if (group->m_next)
	{
		group->m_next->m_prev = group->m_prev;
	}
	if (group == m_groupList)
	{
		m_groupList = group->m_next;
	}

	--m_groupCount;
	group->~b2ParticleGroup();
	m_world->m_blockAllocator.Free(group, sizeof(b2ParticleGroup));
}

void b2ParticleSystem::ComputeWeight()
{
	// calculates the sum of contact-weights for each particle
	// that means dimensionless density
	memset(m_weightBuffer, 0, sizeof(*m_weightBuffer) * m_count);
	for (int32 k = 0; k < m_bodyContactBuffer.GetCount(); k++)
	{
		const b2ParticleBodyContact& contact = m_bodyContactBuffer[k];
		int32 a = contact.index;
		float32 w = contact.weight;
		m_weightBuffer[a] += w;
	}
	for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
	{
		const b2ParticleContact& contact = m_contactBuffer[k];
		int32 a = contact.GetIndexA();
		int32 b = contact.GetIndexB();
		float32 w = contact.GetWeight();
		m_weightBuffer[a] += w;
		m_weightBuffer[b] += w;
	}
}

void b2ParticleSystem::ComputeDepth()
{
	b2ParticleContact* contactGroups = (b2ParticleContact*) m_world->
		m_stackAllocator.Allocate(sizeof(b2ParticleContact) * m_contactBuffer.GetCount());
	int32 contactGroupsCount = 0;
	for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
	{
		const b2ParticleContact& contact = m_contactBuffer[k];
		int32 a = contact.GetIndexA();
		int32 b = contact.GetIndexB();
		const b2ParticleGroup* groupA = m_groupBuffer[a];
		const b2ParticleGroup* groupB = m_groupBuffer[b];
		if (groupA && groupA == groupB &&
			(groupA->m_groupFlags & b2_particleGroupNeedsUpdateDepth))
		{
			contactGroups[contactGroupsCount++] = contact;
		}
	}
	b2ParticleGroup** groupsToUpdate = (b2ParticleGroup**) m_world->
		m_stackAllocator.Allocate(sizeof(b2ParticleGroup*) * m_groupCount);
	int32 groupsToUpdateCount = 0;
	for (b2ParticleGroup* group = m_groupList; group; group = group->GetNext())
	{
		if (group->m_groupFlags & b2_particleGroupNeedsUpdateDepth)
		{
			groupsToUpdate[groupsToUpdateCount++] = group;
			SetGroupFlags(group,
						  group->m_groupFlags &
						  ~b2_particleGroupNeedsUpdateDepth);
			for (int32 i = group->m_firstIndex; i < group->m_lastIndex; i++)
			{
				m_accumulationBuffer[i] = 0;
			}
		}
	}
	// Compute sum of weight of contacts except between different groups.
	for (int32 k = 0; k < contactGroupsCount; k++)
	{
		const b2ParticleContact& contact = contactGroups[k];
		int32 a = contact.GetIndexA();
		int32 b = contact.GetIndexB();
		float32 w = contact.GetWeight();
		m_accumulationBuffer[a] += w;
		m_accumulationBuffer[b] += w;
	}
	b2Assert(m_depthBuffer);
	for (int32 i = 0; i < groupsToUpdateCount; i++)
	{
		const b2ParticleGroup* group = groupsToUpdate[i];
		for (int32 i = group->m_firstIndex; i < group->m_lastIndex; i++)
		{
			float32 w = m_accumulationBuffer[i];
			m_depthBuffer[i] = w < 0.8f ? 0 : b2_maxFloat;
		}
	}
	// The number of iterations is equal to particle number from the deepest
	// particle to the nearest surface particle, and in general it is smaller
	// than sqrt of total particle number.
	int32 iterationCount = (int32)b2Sqrt((float)m_count);
	for (int32 t = 0; t < iterationCount; t++)
	{
		bool updated = false;
		for (int32 k = 0; k < contactGroupsCount; k++)
		{
			const b2ParticleContact& contact = contactGroups[k];
			int32 a = contact.GetIndexA();
			int32 b = contact.GetIndexB();
			float32 r = 1 - contact.GetWeight();
			float32& ap0 = m_depthBuffer[a];
			float32& bp0 = m_depthBuffer[b];
			float32 ap1 = bp0 + r;
			float32 bp1 = ap0 + r;
			if (ap0 > ap1)
			{
				ap0 = ap1;
				updated = true;
			}
			if (bp0 > bp1)
			{
				bp0 = bp1;
				updated = true;
			}
		}
		if (!updated)
		{
			break;
		}
	}
	for (int32 i = 0; i < groupsToUpdateCount; i++)
	{
		const b2ParticleGroup* group = groupsToUpdate[i];
		for (int32 i = group->m_firstIndex; i < group->m_lastIndex; i++)
		{
			float32& p = m_depthBuffer[i];
			if (p < b2_maxFloat)
			{
				p *= m_particleDiameter;
			}
			else
			{
				p = 0;
			}
		}
	}
	m_world->m_stackAllocator.Free(groupsToUpdate);
	m_world->m_stackAllocator.Free(contactGroups);
}

b2ParticleSystem::InsideBoundsEnumerator
b2ParticleSystem::GetInsideBoundsEnumerator(const b2AABB& aabb) const
{
	uint32 lowerTag = computeTag(m_inverseDiameter * aabb.lowerBound.x - 1,
								 m_inverseDiameter * aabb.lowerBound.y - 1);
	uint32 upperTag = computeTag(m_inverseDiameter * aabb.upperBound.x + 1,
								 m_inverseDiameter * aabb.upperBound.y + 1);
	const Proxy* beginProxy = m_proxyBuffer.Begin();
	const Proxy* endProxy = m_proxyBuffer.End();
	const Proxy* firstProxy = std::lower_bound(beginProxy, endProxy, lowerTag);
	const Proxy* lastProxy = std::upper_bound(firstProxy, endProxy, upperTag);
	return InsideBoundsEnumerator(lowerTag, upperTag, firstProxy, lastProxy);
}

inline void b2ParticleSystem::AddContact(int32 a, int32 b,
	b2GrowableBuffer<b2ParticleContact>& contacts) const
{
	b2Vec2 d = m_positionBuffer.data[b] - m_positionBuffer.data[a];
	float32 distBtParticlesSq = b2Dot(d, d);
	if (distBtParticlesSq < m_squaredDiameter)
	{
		float32 invD = b2InvSqrt(distBtParticlesSq);
		b2ParticleContact& contact = contacts.Append();
		contact.SetIndices(a, b);
		contact.SetFlags(m_flagsBuffer.data[a] | m_flagsBuffer.data[b]);
		// 1 - distBtParticles / diameter
		contact.SetWeight(1 - distBtParticlesSq * invD * m_inverseDiameter);
		contact.SetNormal(invD * d);
	}
}

void b2ParticleSystem::FindContacts_Reference(
	b2GrowableBuffer<b2ParticleContact>& contacts) const
{
	const Proxy* beginProxy = m_proxyBuffer.Begin();
	const Proxy* endProxy = m_proxyBuffer.End();

	contacts.SetCount(0);
	for (const Proxy *a = beginProxy, *c = beginProxy; a < endProxy; a++)
	{
		uint32 rightTag = computeRelativeTag(a->tag, 1, 0);
		for (const Proxy* b = a + 1; b < endProxy; b++)
		{
			if (rightTag < b->tag) break;
			AddContact(a->index, b->index, contacts);
		}
		uint32 bottomLeftTag = computeRelativeTag(a->tag, -1, 1);
		for (; c < endProxy; c++)
		{
			if (bottomLeftTag <= c->tag) break;
		}
		uint32 bottomRightTag = computeRelativeTag(a->tag, 1, 1);
		for (const Proxy* b = c; b < endProxy; b++)
		{
			if (bottomRightTag < b->tag) break;
			AddContact(a->index, b->index, contacts);
		}
	}
}

// Put the positions and indices in proxy-order. This allows us to process
// particles with SIMD, since adjacent particles are adjacent in memory.
void b2ParticleSystem::ReorderForFindContact(FindContactInput* reordered,
	                                         int alignedCount) const
{
	int i = 0;
	for (; i < m_count; ++i)
	{
		const int proxyIndex = m_proxyBuffer[i].index;
		FindContactInput& r = reordered[i];
		r.proxyIndex = proxyIndex;
		r.position = m_positionBuffer.data[proxyIndex];
	}

	// We process multiple elements at a time, so we may read off the end of
	// the array. Pad the array with a few elements, so we don't end up
	// outputing spurious contacts.
	for (; i < alignedCount; ++i)
	{
		FindContactInput& r = reordered[i];
		r.proxyIndex = 0;
		r.position = b2Vec2(b2_maxFloat, b2_maxFloat);
	}
}

// Check particles to the right of 'startIndex', outputing FindContactChecks
// until we find an index that is greater than 'bound'. We skip over the
// indices NUM_V32_SLOTS at a time, because they are processed in groups
// in the SIMD function.
inline void b2ParticleSystem::GatherChecksOneParticle(
	const uint32 bound,
	const int startIndex,
	const int particleIndex,
	int* nextUncheckedIndex,
	b2GrowableBuffer<FindContactCheck>& checks) const
{
	// The particles have to be heavily packed together in order for this
	// loop to iterate more than once. In almost all situations, it will
	// iterate less than twice.
	for (int comparatorIndex = startIndex;
		 comparatorIndex < m_count;
	     comparatorIndex += NUM_V32_SLOTS)
	{
		if (m_proxyBuffer[comparatorIndex].tag > bound)
			break;

		FindContactCheck& out = checks.Append();
		out.particleIndex = (uint16)particleIndex;
		out.comparatorIndex = (uint16)comparatorIndex;

		// This is faster inside the 'for' since there are so few iterations.
		if (nextUncheckedIndex != NULL)
		{
			*nextUncheckedIndex = comparatorIndex + NUM_V32_SLOTS;
		}
	}
}

void b2ParticleSystem::GatherChecks(
	b2GrowableBuffer<FindContactCheck>& checks) const
{
	int bottomLeftIndex = 0;
	for (int particleIndex = 0; particleIndex < m_count; ++particleIndex)
	{
		const uint32 particleTag = m_proxyBuffer[particleIndex].tag;

		// Add checks for particles to the right.
		const uint32 rightBound = particleTag + relativeTagRight;
		int nextUncheckedIndex = particleIndex + 1;
		GatherChecksOneParticle(rightBound,
								particleIndex + 1,
								particleIndex,
								&nextUncheckedIndex,
								checks);

		// Find comparator index below and to left of particle.
		const uint32 bottomLeftTag = particleTag + relativeTagBottomLeft;
		for (; bottomLeftIndex < m_count; ++bottomLeftIndex)
		{
			if (bottomLeftTag <= m_proxyBuffer[bottomLeftIndex].tag)
				break;
		}

		// Add checks for particles below.
		const uint32 bottomRightBound = particleTag + relativeTagBottomRight;
		const int bottomStartIndex = b2Max(bottomLeftIndex, nextUncheckedIndex);
		GatherChecksOneParticle(bottomRightBound,
								bottomStartIndex,
								particleIndex,
								NULL,
								checks);
	}
}

#if defined(LIQUIDFUN_SIMD_NEON)
void b2ParticleSystem::FindContacts_Simd(
	b2GrowableBuffer<b2ParticleContact>& contacts) const
{
	contacts.SetCount(0);

	const int alignedCount = m_count + NUM_V32_SLOTS;
	FindContactInput* reordered = (FindContactInput*)
		m_world->m_stackAllocator.Allocate(
			sizeof(FindContactInput) * alignedCount);

	// Put positions and indices into proxy-order.
	// This allows us to efficiently check for contacts using SIMD.
	ReorderForFindContact(reordered, alignedCount);

	// Perform broad-band contact check using tags to approximate
	// positions. This reduces the number of narrow-band contact checks
	// that use actual positions.
	static const int MAX_EXPECTED_CHECKS_PER_PARTICLE = 3;
	b2GrowableBuffer<FindContactCheck> checks(m_world->m_blockAllocator);
	checks.Reserve(MAX_EXPECTED_CHECKS_PER_PARTICLE * m_count);
	GatherChecks(checks);

	// Perform narrow-band contact checks using actual positions.
	// Any particles whose centers are within one diameter of each other are
	// considered contacting.
	FindContactsFromChecks_Simd(reordered, checks.Data(), checks.GetCount(),
								m_squaredDiameter, m_inverseDiameter,
								m_flagsBuffer.data, contacts);

	m_world->m_stackAllocator.Free(reordered);
}
#endif // defined(LIQUIDFUN_SIMD_NEON)

LIQUIDFUN_SIMD_INLINE
void b2ParticleSystem::FindContacts(
	b2GrowableBuffer<b2ParticleContact>& contacts) const
{
	#if defined(LIQUIDFUN_SIMD_NEON)
		FindContacts_Simd(contacts);
	#else
		FindContacts_Reference(contacts);
	#endif

	#if defined(LIQUIDFUN_SIMD_TEST_VS_REFERENCE)
		b2GrowableBuffer<b2ParticleContact>
			reference(m_world->m_blockAllocator);
		FindContacts_Reference(reference);

		b2Assert(contacts.GetCount() == reference.GetCount());
		for (int32 i = 0; i < contacts.GetCount(); ++i)
		{
			b2Assert(contacts[i].ApproximatelyEqual(reference[i]));
		}
	#endif // defined(LIQUIDFUN_SIMD_TEST_VS_REFERENCE)
}

static inline bool b2ParticleContactIsZombie(const b2ParticleContact& contact)
{
	return (contact.GetFlags() & b2_zombieParticle) == b2_zombieParticle;
}

// Get the world's contact filter if any particles with the
// b2_particleContactFilterParticle flag are present in the system.
inline b2ContactFilter* b2ParticleSystem::GetParticleContactFilter() const
{
	return (m_allParticleFlags & b2_particleContactFilterParticle) ?
		m_world->m_contactManager.m_contactFilter : NULL;
}

// Get the world's contact listener if any particles with the
// b2_particleContactListenerParticle flag are present in the system.
inline b2ContactListener* b2ParticleSystem::GetParticleContactListener() const
{
	return (m_allParticleFlags & b2_particleContactListenerParticle) ?
		m_world->m_contactManager.m_contactListener : NULL;
}

// Recalculate 'tag' in proxies using m_positionBuffer.
// The 'tag' is an approximation of position, in left-right, top-bottom order.
void b2ParticleSystem::UpdateProxies_Reference(
	b2GrowableBuffer<Proxy>& proxies) const
{
	const Proxy* const endProxy = proxies.End();
	for (Proxy* proxy = proxies.Begin(); proxy < endProxy; ++proxy)
	{
		int32 i = proxy->index;
		b2Vec2 p = m_positionBuffer.data[i];
		proxy->tag = computeTag(m_inverseDiameter * p.x,
								m_inverseDiameter * p.y);
	}
}

#if defined(LIQUIDFUN_SIMD_NEON)
// static
void b2ParticleSystem::UpdateProxyTags(
	const uint32* const tags,
	b2GrowableBuffer<Proxy>& proxies)
{
	const Proxy* const endProxy = proxies.End();
	for (Proxy* proxy = proxies.Begin(); proxy < endProxy; ++proxy)
	{
		proxy->tag = tags[proxy->index];
	}
}

void b2ParticleSystem::UpdateProxies_Simd(
	b2GrowableBuffer<Proxy>& proxies) const
{
	uint32* tags = (uint32*)
		m_world->m_stackAllocator.Allocate(m_count * sizeof(uint32));

	// Calculate tag for every position.
	// 'tags' array is in position-order.
	CalculateTags_Simd(m_positionBuffer.data, m_count,
					   m_inverseDiameter, tags);

	// Update 'tag' element in the 'proxies' array to the new values.
	UpdateProxyTags(tags, proxies);

	m_world->m_stackAllocator.Free(tags);
}
#endif // defined(LIQUIDFUN_SIMD_NEON)

// static
bool b2ParticleSystem::ProxyBufferHasIndex(
	int32 index, const Proxy* const a, int count)
{
	for (int j = 0; j < count; ++j)
	{
		if (a[j].index == index)
			return true;
	}
	return false;
}

// static
int b2ParticleSystem::NumProxiesWithSameTag(
	const Proxy* const a, const Proxy* const b, int count)
{
	const uint32 tag = a[0].tag;
	for (int num = 0; num < count; ++num)
	{
		if (a[num].tag != tag || b[num].tag != tag)
			return num;
	}
	return count;
}

// Precondition: both 'a' and 'b' should be sorted by tag, but don't need to be
// sorted by index.
// static
bool b2ParticleSystem::AreProxyBuffersTheSame(const b2GrowableBuffer<Proxy>& a,
								   			  const b2GrowableBuffer<Proxy>& b)
{
	if (a.GetCount() != b.GetCount())
		return false;

	// A given tag may have several indices. The order of these indices is
	// not important, but the set must be equivalent.
	for (int i = 0; i < a.GetCount();)
	{
		const int numWithSameTag = NumProxiesWithSameTag(
			&a[i], &b[i], a.GetCount() - i);
		if (numWithSameTag == 0)
			return false;

		for (int j = 0; j < numWithSameTag; ++j)
		{
			const bool hasIndex = ProxyBufferHasIndex(
				a[i + j].index, &b[i], numWithSameTag);
			if (!hasIndex)
				return false;
		}

		i += numWithSameTag;
	}
	return true;
}

LIQUIDFUN_SIMD_INLINE
void b2ParticleSystem::UpdateProxies(
	b2GrowableBuffer<Proxy>& proxies) const
{
	#if defined(LIQUIDFUN_SIMD_TEST_VS_REFERENCE)
		b2GrowableBuffer<Proxy> reference(proxies);
	#endif

	#if defined(LIQUIDFUN_SIMD_NEON)
		UpdateProxies_Simd(proxies);
	#else
		UpdateProxies_Reference(proxies);
	#endif

	#if defined(LIQUIDFUN_SIMD_TEST_VS_REFERENCE)
		UpdateProxies_Reference(reference);
		b2Assert(AreProxyBuffersTheSame(proxies, reference));
	#endif
}


// Sort the proxy array by 'tag'. This orders the particles into rows that
// run left-to-right, top-to-bottom. The rows are spaced m_particleDiameter
// apart, such that a particle in one row can only collide with the rows
// immediately above and below it. This ordering makes collision computation
// tractable.
//
// TODO OPT: The sort is a hot spot on the profiles. We could use SIMD to
// speed this up. See http://www.vldb.org/pvldb/1/1454171.pdf for an excellent
// explanation of a SIMD mergesort algorithm.
void b2ParticleSystem::SortProxies(b2GrowableBuffer<Proxy>& proxies) const
{
	std::sort(proxies.Begin(), proxies.End());
}

class b2ParticleContactRemovePredicate
{
public:
	b2ParticleContactRemovePredicate(
		b2ParticleSystem* system,
		b2ContactFilter* contactFilter) :
		m_system(system),
		m_contactFilter(contactFilter)
	{}

	bool operator()(const b2ParticleContact& contact)
	{
	    return (contact.GetFlags() & b2_particleContactFilterParticle)
	        && !m_contactFilter->ShouldCollide(m_system, contact.GetIndexA(),
	        								   contact.GetIndexB());
	}

private:
	b2ParticleSystem* m_system;
	b2ContactFilter* m_contactFilter;
};

// Only changes 'contacts', but the contact filter has a non-const 'this'
// pointer, so this member function cannot be const.
void b2ParticleSystem::FilterContacts(
	b2GrowableBuffer<b2ParticleContact>& contacts)
{
	// Optionally filter the contact.
	b2ContactFilter* const contactFilter = GetParticleContactFilter();
	if (contactFilter == NULL)
		return;

	contacts.RemoveIf(b2ParticleContactRemovePredicate(this, contactFilter));
}

void b2ParticleSystem::NotifyContactListenerPreContact(
	b2ParticlePairSet* particlePairs) const
{
	b2ContactListener* const contactListener = GetParticleContactListener();
	if (contactListener == NULL)
		return;

	particlePairs->Initialize(m_contactBuffer.Begin(),
							  m_contactBuffer.GetCount(),
						      GetFlagsBuffer());
}

// Note: This function is not const because 'this' in BeginContact and
// EndContact callbacks must be non-const. However, this function itself
// does not change any internal data (though the callbacks might).
void b2ParticleSystem::NotifyContactListenerPostContact(
	b2ParticlePairSet& particlePairs)
{
	b2ContactListener* const contactListener = GetParticleContactListener();
	if (contactListener == NULL)
		return;

	// Loop through all new contacts, reporting any new ones, and
	// "invalidating" the ones that still exist.
	const b2ParticleContact* const endContact = m_contactBuffer.End();
	for (b2ParticleContact* contact = m_contactBuffer.Begin();
		 contact < endContact; ++contact)
	{
		ParticlePair pair;
		pair.first = contact->GetIndexA();
		pair.second = contact->GetIndexB();
		const int32 itemIndex = particlePairs.Find(pair);
		if (itemIndex >= 0)
		{
			// Already touching, ignore this contact.
			particlePairs.Invalidate(itemIndex);
		}
		else
		{
			// Just started touching, inform the listener.
			contactListener->BeginContact(this, contact);
		}
	}

	// Report particles that are no longer touching.
	// That is, any pairs that were not invalidated above.
	const int32 pairCount = particlePairs.GetCount();
	const ParticlePair* const pairs = particlePairs.GetBuffer();
	const int8* const valid = particlePairs.GetValidBuffer();
	for (int32 i = 0; i < pairCount; ++i)
	{
		if (valid[i])
		{
			contactListener->EndContact(this, pairs[i].first,
										pairs[i].second);
		}
	}
}

void b2ParticleSystem::UpdateContacts(bool exceptZombie)
{
	UpdateProxies(m_proxyBuffer);
	SortProxies(m_proxyBuffer);

	b2ParticlePairSet particlePairs(&m_world->m_stackAllocator);
	NotifyContactListenerPreContact(&particlePairs);

	FindContacts(m_contactBuffer);
	FilterContacts(m_contactBuffer);

	NotifyContactListenerPostContact(particlePairs);

	if (exceptZombie)
	{
		m_contactBuffer.RemoveIf(b2ParticleContactIsZombie);
	}
}

void b2ParticleSystem::DetectStuckParticle(int32 particle)
{
	// Detect stuck particles
	//
	// The basic algorithm is to allow the user to specify an optional
	// threshold where we detect whenever a particle is contacting
	// more than one fixture for more than threshold consecutive
	// steps. This is considered to be "stuck", and these are put
	// in a list the user can query per step, if enabled, to deal with
	// such particles.

	if (m_stuckThreshold <= 0)
	{
		return;
	}

	// Get the state variables for this particle.
	int32 * const consecutiveCount =
			&m_consecutiveContactStepsBuffer.data[particle];
	int32 * const lastStep = &m_lastBodyContactStepBuffer.data[particle];
	int32 * const bodyCount = &m_bodyContactCountBuffer.data[particle];

	// This is only called when there is a body contact for this particle.
	++(*bodyCount);

	// We want to only trigger detection once per step, the first time we
	// contact more than one fixture in a step for a given particle.
	if (*bodyCount == 2)
	{
		++(*consecutiveCount);
		if (*consecutiveCount > m_stuckThreshold)
		{
			int32& newStuckParticle = m_stuckParticleBuffer.Append();
			newStuckParticle = particle;
		}
	}
	*lastStep = m_timestamp;
}

// Get the world's contact listener if any particles with the
// b2_fixtureContactListenerParticle flag are present in the system.
inline b2ContactListener* b2ParticleSystem::GetFixtureContactListener() const
{
	return (m_allParticleFlags & b2_fixtureContactListenerParticle) ?
		m_world->m_contactManager.m_contactListener : NULL;
}

// Get the world's contact filter if any particles with the
// b2_fixtureContactFilterParticle flag are present in the system.
inline b2ContactFilter* b2ParticleSystem::GetFixtureContactFilter() const
{
	return (m_allParticleFlags & b2_fixtureContactFilterParticle) ?
		m_world->m_contactManager.m_contactFilter : NULL;
}

/// Compute the axis-aligned bounding box for all particles contained
/// within this particle system.
/// @param aabb Returns the axis-aligned bounding box of the system.
void b2ParticleSystem::ComputeAABB(b2AABB* const aabb) const
{
	const int32 particleCount = GetParticleCount();
	b2Assert(aabb);
	aabb->lowerBound.x = +b2_maxFloat;
	aabb->lowerBound.y = +b2_maxFloat;
	aabb->upperBound.x = -b2_maxFloat;
	aabb->upperBound.y = -b2_maxFloat;

	for (int32 i = 0; i < particleCount; i++)
	{
		b2Vec2 p = m_positionBuffer.data[i];
		aabb->lowerBound = b2Min(aabb->lowerBound, p);
		aabb->upperBound = b2Max(aabb->upperBound, p);
	}
	aabb->lowerBound.x -= m_particleDiameter;
	aabb->lowerBound.y -= m_particleDiameter;
	aabb->upperBound.x += m_particleDiameter;
	aabb->upperBound.y += m_particleDiameter;
}

// Associate a memory allocator with this object.
FixedSetAllocator::FixedSetAllocator(
		b2StackAllocator* allocator) :
	m_buffer(NULL), m_valid(NULL), m_count(0), m_allocator(allocator)
{
	b2Assert(allocator);
}

// Allocate internal storage for this object.
int32 FixedSetAllocator::Allocate(
	const int32 itemSize, const int32 count)
{
	Clear();
	if (count)
	{
		m_buffer = m_allocator->Allocate(
			(itemSize + sizeof(*m_valid)) * count);
		b2Assert(m_buffer);
		m_valid = (int8*)m_buffer + (itemSize * count);
		memset(m_valid, 1, sizeof(*m_valid) * count);
		m_count = count;
	}
	return m_count;
}

// Deallocate the internal buffer if it's allocated.
void FixedSetAllocator::Clear()
{
	if (m_buffer)
	{
		m_allocator->Free(m_buffer);
		m_buffer = NULL;
        m_count = 0;
	}
}

// Search set for item returning the index of the item if it's found, -1
// otherwise.
template<typename T>
static int32 FindItemIndexInFixedSet(const TypedFixedSetAllocator<T>& set,
									 const T& item)
{
	if (set.GetCount())
	{
		const T* buffer = set.GetBuffer();
		const T* last = buffer + set.GetCount();
		const T* found = std::lower_bound( buffer, buffer + set.GetCount(),
											item, T::Compare);
		if( found != last )
		{
			return set.GetIndex( found );
		}
	}
	return -1;
}

// Initialize from a set of particle / body contacts for particles
// that have the b2_fixtureContactListenerParticle flag set.
void FixtureParticleSet::Initialize(
	const b2ParticleBodyContact * const bodyContacts,
	const int32 numBodyContacts,
	const uint32 * const particleFlagsBuffer)
{
	Clear();
	if (Allocate(numBodyContacts))
	{
		FixtureParticle* set = GetBuffer();
		int32 insertedContacts = 0;
		for (int32 i = 0; i < numBodyContacts; ++i)
		{
			FixtureParticle* const fixtureParticle = &set[i];
			const b2ParticleBodyContact& bodyContact = bodyContacts[i];
			if (bodyContact.index == b2_invalidParticleIndex ||
				!(particleFlagsBuffer[bodyContact.index] &
				  b2_fixtureContactListenerParticle))
			{
				continue;
			}
			fixtureParticle->first = bodyContact.fixture;
			fixtureParticle->second = bodyContact.index;
			insertedContacts++;
		}
		SetCount(insertedContacts);
		std::sort(set, set + insertedContacts, FixtureParticle::Compare);
	}
}

// Find the index of a particle / fixture pair in the set or -1 if it's not
// present.
int32 FixtureParticleSet::Find(
	const FixtureParticle& fixtureParticle) const
{
	return FindItemIndexInFixedSet(*this, fixtureParticle);
}

// Initialize from a set of particle contacts.
void b2ParticlePairSet::Initialize(
	const b2ParticleContact * const contacts, const int32 numContacts,
	const uint32 * const particleFlagsBuffer)
{
	Clear();
	if (Allocate(numContacts))
	{
		ParticlePair* set = GetBuffer();
		int32 insertedContacts = 0;
		for (int32 i = 0; i < numContacts; ++i)
		{
			ParticlePair* const pair = &set[i];
			const b2ParticleContact& contact = contacts[i];
			if (contact.GetIndexA() == b2_invalidParticleIndex ||
				contact.GetIndexB() == b2_invalidParticleIndex ||
				!((particleFlagsBuffer[contact.GetIndexA()] |
				   particleFlagsBuffer[contact.GetIndexB()]) &
				  b2_particleContactListenerParticle))
			{
				continue;
			}
			pair->first = contact.GetIndexA();
			pair->second = contact.GetIndexB();
			insertedContacts++;
		}
		SetCount(insertedContacts);
		std::sort(set, set + insertedContacts, ParticlePair::Compare);
	}
}

// Find the index of a particle pair in the set or -1 if it's not present.
int32 b2ParticlePairSet::Find(const ParticlePair& pair) const
{
	int32 index = FindItemIndexInFixedSet(*this, pair);
	if (index < 0)
	{
		ParticlePair swapped;
		swapped.first = pair.second;
		swapped.second = pair.first;
		index = FindItemIndexInFixedSet(*this, swapped);
	}
	return index;
}

/// Callback class to receive pairs of fixtures and particles which may be
/// overlapping. Used as an argument of b2World::QueryAABB.
class b2FixtureParticleQueryCallback : public b2QueryCallback
{
public:
	explicit b2FixtureParticleQueryCallback(b2ParticleSystem* system)
	{
		m_system = system;
	}

private:
	// Skip reporting particles.
	bool ShouldQueryParticleSystem(const b2ParticleSystem* system)
	{
		B2_NOT_USED(system);
		return false;
	}

	// Receive a fixture and call ReportFixtureAndParticle() for each particle
	// inside aabb of the fixture.
	bool ReportFixture(b2Fixture* fixture)
	{
		if (fixture->IsSensor())
		{
			return true;
		}
		const b2Shape* shape = fixture->GetShape();
		int32 childCount = shape->GetChildCount();
		for (int32 childIndex = 0; childIndex < childCount; childIndex++)
		{
			b2AABB aabb = fixture->GetAABB(childIndex);
			b2ParticleSystem::InsideBoundsEnumerator enumerator =
								m_system->GetInsideBoundsEnumerator(aabb);
			int32 index;
			while ((index = enumerator.GetNext()) >= 0)
			{
				ReportFixtureAndParticle(fixture, childIndex, index);
			}
		}
		return true;
	}

	// Receive a fixture and a particle which may be overlapping.
	virtual void ReportFixtureAndParticle(
						b2Fixture* fixture, int32 childIndex, int32 index) = 0;

protected:
	b2ParticleSystem* m_system;
};

void b2ParticleSystem::NotifyBodyContactListenerPreContact(
	FixtureParticleSet* fixtureSet) const
{
	b2ContactListener* const contactListener = GetFixtureContactListener();
	if (contactListener == NULL)
		return;

	fixtureSet->Initialize(m_bodyContactBuffer.Begin(),
						   m_bodyContactBuffer.GetCount(),
						   GetFlagsBuffer());
}

// If a contact listener is present and the contact is just starting
// report the contact.  If the contact is already in progress invalid
// the contact from m_fixtureSet.
void b2ParticleSystem::NotifyBodyContactListenerPostContact(
	FixtureParticleSet& fixtureSet)
{
	b2ContactListener* const contactListener = GetFixtureContactListener();
	if (contactListener == NULL)
		return;

	// Loop through all new contacts, reporting any new ones, and
	// "invalidating" the ones that still exist.
	for (b2ParticleBodyContact* contact = m_bodyContactBuffer.Begin();
		 contact != m_bodyContactBuffer.End(); ++contact)
	{
		b2Assert(contact);
		FixtureParticle fixtureParticleToFind;
		fixtureParticleToFind.first = contact->fixture;
		fixtureParticleToFind.second = contact->index;
		const int32 index = fixtureSet.Find(fixtureParticleToFind);
		if (index >= 0)
		{
			// Already touching remove this from the set.
			fixtureSet.Invalidate(index);
		}
		else
		{
			// Just started touching, report it!
			contactListener->BeginContact(this, contact);
		}
	}

	// If the contact listener is enabled, report all fixtures that are no
	// longer in contact with particles.
	const FixtureParticle* const fixtureParticles = fixtureSet.GetBuffer();
	const int8* const fixtureParticlesValid = fixtureSet.GetValidBuffer();
	const int32 fixtureParticleCount = fixtureSet.GetCount();
	for (int32 i = 0; i < fixtureParticleCount; ++i)
	{
		if (fixtureParticlesValid[i])
		{
			const FixtureParticle* const fixtureParticle =
				&fixtureParticles[i];
			contactListener->EndContact(fixtureParticle->first, this,
										fixtureParticle->second);
		}
	}
}


void b2ParticleSystem::UpdateBodyContacts()
{
	// If the particle contact listener is enabled, generate a set of
	// fixture / particle contacts.
	FixtureParticleSet fixtureSet(&m_world->m_stackAllocator);
	NotifyBodyContactListenerPreContact(&fixtureSet);

	if (m_stuckThreshold > 0)
	{
		const int32 particleCount = GetParticleCount();
		for (int32 i = 0; i < particleCount; i++)
		{
			// Detect stuck particles, see comment in
			// b2ParticleSystem::DetectStuckParticle()
			m_bodyContactCountBuffer.data[i] = 0;
			if (m_timestamp > (m_lastBodyContactStepBuffer.data[i] + 1))
			{
				m_consecutiveContactStepsBuffer.data[i] = 0;
			}
		}
	}
	m_bodyContactBuffer.SetCount(0);
	m_stuckParticleBuffer.SetCount(0);

	class UpdateBodyContactsCallback : public b2FixtureParticleQueryCallback
	{
		// Call the contact filter if it's set, to determine whether to
		// filter this contact.  Returns true if contact calculations should
		// be performed, false otherwise.
		inline bool ShouldCollide(b2Fixture * const fixture,
								  int32 particleIndex)
		{
			if (m_contactFilter)
			{
				const uint32* const flags = m_system->GetFlagsBuffer();
				if (flags[particleIndex] & b2_fixtureContactFilterParticle)
				{
					return m_contactFilter->ShouldCollide(fixture, m_system,
														  particleIndex);
				}
			}
			return true;
		}

		void ReportFixtureAndParticle(
								b2Fixture* fixture, int32 childIndex, int32 a)
		{
			b2Vec2 ap = m_system->m_positionBuffer.data[a];
			float32 d;
			b2Vec2 n;
			fixture->ComputeDistance(ap, &d, &n, childIndex);
			if (d < m_system->m_particleDiameter && ShouldCollide(fixture, a))
			{
				b2Body* b = fixture->GetBody();
				b2Vec2 bp = b->GetWorldCenter();
				float32 bm = b->GetMass();
				float32 bI =
					b->GetInertia() - bm * b->GetLocalCenter().LengthSquared();
				float32 invBm = bm > 0 ? 1 / bm : 0;
				float32 invBI = bI > 0 ? 1 / bI : 0;
				float32 invAm =
					m_system->m_flagsBuffer.data[a] &
					b2_wallParticle ? 0 : m_system->GetParticleInvMass();
				b2Vec2 rp = ap - bp;
				float32 rpn = b2Cross(rp, n);
				float32 invM = invAm + invBm + invBI * rpn * rpn;

				b2ParticleBodyContact& contact =
					m_system->m_bodyContactBuffer.Append();
				contact.index = a;
				contact.body = b;
				contact.fixture = fixture;
				contact.weight = 1 - d * m_system->m_inverseDiameter;
				contact.normal = -n;
				contact.mass = invM > 0 ? 1 / invM : 0;
				m_system->DetectStuckParticle(a);
			}
		}

		b2ContactFilter* m_contactFilter;

	public:
		UpdateBodyContactsCallback(
			b2ParticleSystem* system, b2ContactFilter* contactFilter):
			b2FixtureParticleQueryCallback(system)
		{
			m_contactFilter = contactFilter;
		}
	} callback(this, GetFixtureContactFilter());

	b2AABB aabb;
	ComputeAABB(&aabb);
	m_world->QueryAABB(&callback, aabb);

	if (m_def.strictContactCheck)
	{
		RemoveSpuriousBodyContacts();
	}

	NotifyBodyContactListenerPostContact(fixtureSet);
}

void b2ParticleSystem::RemoveSpuriousBodyContacts()
{
	// At this point we have a list of contact candidates based on AABB
	// overlap.The AABB query that  generated this returns all collidable
	// fixtures overlapping particle bounding boxes.  This breaks down around
	// vertices where two shapes intersect, such as a "ground" surface made
	// of multiple b2PolygonShapes; it potentially applies a lot of spurious
	// impulses from normals that should not actually contribute.  See the
	// Ramp example in Testbed.
	//
	// To correct for this, we apply this algorithm:
	//   * sort contacts by particle and subsort by weight (nearest to farthest)
	//   * for each contact per particle:
	//      - project a point at the contact distance along the inverse of the
	//        contact normal
	//      - if this intersects the fixture that generated the contact, apply
	//         it, otherwise discard as impossible
	//      - repeat for up to n nearest contacts, currently we get good results
	//        from n=3.
	std::sort(m_bodyContactBuffer.Begin(), m_bodyContactBuffer.End(),
				b2ParticleSystem::BodyContactCompare);

	int32 discarded = 0;
	std::remove_if(m_bodyContactBuffer.Begin(),
					m_bodyContactBuffer.End(),
					b2ParticleBodyContactRemovePredicate(this, &discarded));

	m_bodyContactBuffer.SetCount(m_bodyContactBuffer.GetCount() - discarded);
}

bool b2ParticleSystem::BodyContactCompare(const b2ParticleBodyContact &lhs,
										  const b2ParticleBodyContact &rhs)
{
	if (lhs.index == rhs.index)
	{
		// Subsort by weight, decreasing.
		return lhs.weight > rhs.weight;
	}
	return lhs.index < rhs.index;
}


void b2ParticleSystem::SolveCollision(const b2TimeStep& step)
{
	// This function detects particles which are crossing boundary of bodies
	// and modifies velocities of them so that they will move just in front of
	// boundary. This function function also applies the reaction force to
	// bodies as precisely as the numerical stability is kept.
	b2AABB aabb;
	aabb.lowerBound.x = +b2_maxFloat;
	aabb.lowerBound.y = +b2_maxFloat;
	aabb.upperBound.x = -b2_maxFloat;
	aabb.upperBound.y = -b2_maxFloat;
	for (int32 i = 0; i < m_count; i++)
	{
		b2Vec2 v = m_velocityBuffer.data[i];
		b2Vec2 p1 = m_positionBuffer.data[i];
		b2Vec2 p2 = p1 + step.dt * v;
		aabb.lowerBound = b2Min(aabb.lowerBound, b2Min(p1, p2));
		aabb.upperBound = b2Max(aabb.upperBound, b2Max(p1, p2));
	}
	class SolveCollisionCallback : public b2FixtureParticleQueryCallback
	{
		// Call the contact filter if it's set, to determine whether to
		// filter this contact.  Returns true if contact calculations should
		// be performed, false otherwise.
		inline bool ShouldCollide(b2Fixture * const fixture,
								  int32 particleIndex)
		{
			if (m_contactFilter) {
				const uint32* const flags = m_system->GetFlagsBuffer();
				if (flags[particleIndex] & b2_fixtureContactFilterParticle) {
					return m_contactFilter->ShouldCollide(fixture, m_system,
														  particleIndex);
				}
			}
			return true;
		}

		void ReportFixtureAndParticle(
								b2Fixture* fixture, int32 childIndex, int32 a)
		{
			if (ShouldCollide(fixture, a)) {
				b2Body* body = fixture->GetBody();
				b2Vec2 ap = m_system->m_positionBuffer.data[a];
				b2Vec2 av = m_system->m_velocityBuffer.data[a];
				b2RayCastOutput output;
				b2RayCastInput input;
				if (m_system->m_iterationIndex == 0)
				{
					// Put 'ap' in the local space of the previous frame
					b2Vec2 p1 = b2MulT(body->m_xf0, ap);
					if (fixture->GetShape()->GetType() == b2Shape::e_circle)
					{
						// Make relative to the center of the circle
						p1 -= body->GetLocalCenter();
						// Re-apply rotation about the center of the
						// circle
						p1 = b2Mul(body->m_xf0.q, p1);
						// Subtract rotation of the current frame
						p1 = b2MulT(body->m_xf.q, p1);
						// Return to local space
						p1 += body->GetLocalCenter();
					}
					// Return to global space and apply rotation of current frame
					input.p1 = b2Mul(body->m_xf, p1);
				}
				else
				{
					input.p1 = ap;
				}
				input.p2 = ap + m_step.dt * av;
				input.maxFraction = 1;
				if (fixture->RayCast(&output, input, childIndex))
				{
					b2Vec2 n = output.normal;
					b2Vec2 p =
						(1 - output.fraction) * input.p1 +
						output.fraction * input.p2 +
						b2_linearSlop * n;
					b2Vec2 v = m_step.inv_dt * (p - ap);
					m_system->m_velocityBuffer.data[a] = v;
					b2Vec2 f = m_step.inv_dt *
						m_system->GetParticleMass() * (av - v);
					m_system->ParticleApplyForce(a, f);
				}
			}
		}

		b2TimeStep m_step;
		b2ContactFilter* m_contactFilter;

	public:
		SolveCollisionCallback(
			b2ParticleSystem* system, const b2TimeStep& step, b2ContactFilter* contactFilter) :
			b2FixtureParticleQueryCallback(system)
		{
			m_step = step;
			m_contactFilter = contactFilter;
		}
	} callback(this, step, GetFixtureContactFilter());
	m_world->QueryAABB(&callback, aabb);
}

void b2ParticleSystem::SolveBarrier(const b2TimeStep& step)
{
	// If a particle is passing between paired barrier particles,
	// its velocity will be decelerated to avoid passing.
	for (int32 i = 0; i < m_count; i++)
	{
		uint32 flags = m_flagsBuffer.data[i];
		static const uint32 k_barrierWallFlags =
										b2_barrierParticle | b2_wallParticle;
		if ((flags & k_barrierWallFlags) == k_barrierWallFlags)
		{
			m_velocityBuffer.data[i].SetZero();
		}
	}
	float32 tmax = b2_barrierCollisionTime * step.dt;
	for (int32 k = 0; k < m_pairBuffer.GetCount(); k++)
	{
		const b2ParticlePair& pair = m_pairBuffer[k];
		if (pair.flags & b2_barrierParticle)
		{
			int32 a = pair.indexA;
			int32 b = pair.indexB;
			b2Vec2 pa = m_positionBuffer.data[a];
			b2Vec2 pb = m_positionBuffer.data[b];
			b2AABB aabb;
			aabb.lowerBound = b2Min(pa, pb);
			aabb.upperBound = b2Max(pa, pb);
			b2ParticleGroup *aGroup = m_groupBuffer[a];
			b2ParticleGroup *bGroup = m_groupBuffer[b];
			b2Vec2 va = GetLinearVelocity(aGroup, a, pa);
			b2Vec2 vb = GetLinearVelocity(bGroup, b, pb);
			b2Vec2 pba = pb - pa;
			b2Vec2 vba = vb - va;
			InsideBoundsEnumerator enumerator = GetInsideBoundsEnumerator(aabb);
			int32 c;
			while ((c = enumerator.GetNext()) >= 0)
			{
				b2Vec2 pc = m_positionBuffer.data[c];
				b2ParticleGroup *cGroup = m_groupBuffer[c];
				if (aGroup != cGroup && bGroup != cGroup)
				{
					b2Vec2 vc = GetLinearVelocity(cGroup, c, pc);
					// Solve the equation below:
					//   (1-s)*(pa+t*va)+s*(pb+t*vb) = pc+t*vc
					// which expresses that the particle c will pass a line
					// connecting the particles a and b at the time of t.
					// if s is between 0 and 1, c will pass between a and b.
					b2Vec2 pca = pc - pa;
					b2Vec2 vca = vc - va;
					float32 e2 = b2Cross(vba, vca);
					float32 e1 = b2Cross(pba, vca) - b2Cross(pca, vba);
					float32 e0 = b2Cross(pba, pca);
					float32 s, t;
					b2Vec2 qba, qca;
					if (e2 == 0)
					{
						if (e1 == 0) continue;
						t = - e0 / e1;
						if (!(t >= 0 && t < tmax)) continue;
						qba = pba + t * vba;
						qca = pca + t * vca;
						s = b2Dot(qba, qca) / b2Dot(qba, qba);
						if (!(s >= 0 && s <= 1)) continue;
					}
					else
					{
						float32 det = e1 * e1 - 4 * e0 * e2;
						if (det < 0) continue;
						float32 sqrtDet = b2Sqrt(det);
						float32 t1 = (- e1 - sqrtDet) / (2 * e2);
						float32 t2 = (- e1 + sqrtDet) / (2 * e2);
						if (t1 > t2) b2Swap(t1, t2);
						t = t1;
						qba = pba + t * vba;
						qca = pca + t * vca;
						s = b2Dot(qba, qca) / b2Dot(qba, qba);
						if (!(t >= 0 && t < tmax && s >= 0 && s <= 1))
						{
							t = t2;
							if (!(t >= 0 && t < tmax)) continue;
							qba = pba + t * vba;
							qca = pca + t * vca;
							s = b2Dot(qba, qca) / b2Dot(qba, qba);
							if (!(s >= 0 && s <= 1)) continue;
						}
					}
					// Apply a force to particle c so that it will have the
					// interpolated velocity at the collision point on line ab.
					b2Vec2 dv = va + s * vba - vc;
					b2Vec2 f = GetParticleMass() * dv;
					if (IsRigidGroup(cGroup))
					{
						// If c belongs to a rigid group, the force will be
						// distributed in the group.
						float32 mass = cGroup->GetMass();
						float32 inertia = cGroup->GetInertia();
						if (mass > 0)
						{
							cGroup->m_linearVelocity += 1 / mass * f;
						}
						if (inertia > 0)
						{
							cGroup->m_angularVelocity +=
								b2Cross(pc - cGroup->GetCenter(), f) / inertia;
						}
					}
					else
					{
						m_velocityBuffer.data[c] += dv;
					}
					// Apply a reversed force to particle c after particle
					// movement so that momentum will be preserved.
					ParticleApplyForce(c, -step.inv_dt * f);
				}
			}
		}
	}
}

void b2ParticleSystem::Solve(const b2TimeStep& step)
{
	if (m_count == 0)
	{
		return;
	}
	// If particle lifetimes are enabled, destroy particles that are too old.
	if (m_expirationTimeBuffer.data)
	{
		SolveLifetimes(step);
	}
	if (m_allParticleFlags & b2_zombieParticle)
	{
		SolveZombie();
	}
	if (m_needsUpdateAllParticleFlags)
	{
		UpdateAllParticleFlags();
	}
	if (m_needsUpdateAllGroupFlags)
	{
		UpdateAllGroupFlags();
	}
	if (m_paused)
	{
		return;
	}
	for (m_iterationIndex = 0;
		m_iterationIndex < step.particleIterations;
		m_iterationIndex++)
	{
		++m_timestamp;
		b2TimeStep subStep = step;
		subStep.dt /= step.particleIterations;
		subStep.inv_dt *= step.particleIterations;
		UpdateContacts(false);
		UpdateBodyContacts();
		ComputeWeight();
		if (m_allGroupFlags & b2_particleGroupNeedsUpdateDepth)
		{
			ComputeDepth();
		}
		if (m_allParticleFlags & b2_reactiveParticle)
		{
			UpdatePairsAndTriadsWithReactiveParticles();
		}
		if (m_hasForce)
		{
			SolveForce(subStep);
		}
		if (m_allParticleFlags & b2_viscousParticle)
		{
			SolveViscous();
		}
		if (m_allParticleFlags & b2_repulsiveParticle)
		{
			SolveRepulsive(subStep);
		}
		if (m_allParticleFlags & b2_powderParticle)
		{
			SolvePowder(subStep);
		}
		if (m_allParticleFlags & b2_tensileParticle)
		{
			SolveTensile(subStep);
		}
		if (m_allGroupFlags & b2_solidParticleGroup)
		{
			SolveSolid(subStep);
		}
		if (m_allParticleFlags & b2_colorMixingParticle)
		{
			SolveColorMixing();
		}
		SolveGravity(subStep);
		if (m_allParticleFlags & b2_staticPressureParticle)
		{
			SolveStaticPressure(subStep);
		}
		SolvePressure(subStep);
		SolveDamping(subStep);
		if (m_allParticleFlags & k_extraDampingFlags)
		{
			SolveExtraDamping();
		}
		// SolveElastic and SolveSpring refer the current velocities for
		// numerical stability, they should be called as late as possible.
		if (m_allParticleFlags & b2_elasticParticle)
		{
			SolveElastic(subStep);
		}
		if (m_allParticleFlags & b2_springParticle)
		{
			SolveSpring(subStep);
		}
		LimitVelocity(subStep);
		if (m_allGroupFlags & b2_rigidParticleGroup)
		{
			SolveRigidDamping();
		}
		if (m_allParticleFlags & b2_barrierParticle)
		{
			SolveBarrier(subStep);
		}
		// SolveCollision, SolveRigid and SolveWall should be called after
		// other force functions because they may require particles to have
		// specific velocities.
		SolveCollision(subStep);
		if (m_allGroupFlags & b2_rigidParticleGroup)
		{
			SolveRigid(subStep);
		}
		if (m_allParticleFlags & b2_wallParticle)
		{
			SolveWall();
		}
		// The particle positions can be updated only at the end of substep.
		for (int32 i = 0; i < m_count; i++)
		{
			m_positionBuffer.data[i] += subStep.dt * m_velocityBuffer.data[i];
		}
	}
}

void b2ParticleSystem::UpdateAllParticleFlags()
{
	m_allParticleFlags = 0;
	for (int32 i = 0; i < m_count; i++)
	{
		m_allParticleFlags |= m_flagsBuffer.data[i];
	}
	m_needsUpdateAllParticleFlags = false;
}

void b2ParticleSystem::UpdateAllGroupFlags()
{
	m_allGroupFlags = 0;
	for (const b2ParticleGroup* group = m_groupList; group;
		 group = group->GetNext())
	{
		m_allGroupFlags |= group->m_groupFlags;
	}
	m_needsUpdateAllGroupFlags = false;
}

void b2ParticleSystem::LimitVelocity(const b2TimeStep& step)
{
	float32 criticalVelocitySquared = GetCriticalVelocitySquared(step);
	for (int32 i = 0; i < m_count; i++)
	{
		b2Vec2& v = m_velocityBuffer.data[i];
		float32 v2 = b2Dot(v, v);
		if (v2 > criticalVelocitySquared)
		{
			v *= b2Sqrt(criticalVelocitySquared / v2);
		}
	}
}

void b2ParticleSystem::SolveGravity(const b2TimeStep& step)
{
	b2Vec2 gravity = step.dt * m_def.gravityScale * m_world->GetGravity();
	for (int32 i = 0; i < m_count; i++)
	{
		m_velocityBuffer.data[i] += gravity;
	}
}

void b2ParticleSystem::SolveStaticPressure(const b2TimeStep& step)
{
	m_staticPressureBuffer = RequestBuffer(m_staticPressureBuffer);
	float32 criticalPressure = GetCriticalPressure(step);
	float32 pressurePerWeight = m_def.staticPressureStrength * criticalPressure;
	float32 maxPressure = b2_maxParticlePressure * criticalPressure;
	float32 relaxation = m_def.staticPressureRelaxation;
	/// Compute pressure satisfying the modified Poisson equation:
	///     Sum_for_j((p_i - p_j) * w_ij) + relaxation * p_i =
	///     pressurePerWeight * (w_i - b2_minParticleWeight)
	/// by iterating the calculation:
	///     p_i = (Sum_for_j(p_j * w_ij) + pressurePerWeight *
	///           (w_i - b2_minParticleWeight)) / (w_i + relaxation)
	/// where
	///     p_i and p_j are static pressure of particle i and j
	///     w_ij is contact weight between particle i and j
	///     w_i is sum of contact weight of particle i
	for (int32 t = 0; t < m_def.staticPressureIterations; t++)
	{
		memset(m_accumulationBuffer, 0,
			   sizeof(*m_accumulationBuffer) * m_count);
		for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
		{
			const b2ParticleContact& contact = m_contactBuffer[k];
			if (contact.GetFlags() & b2_staticPressureParticle)
			{
				int32 a = contact.GetIndexA();
				int32 b = contact.GetIndexB();
				float32 w = contact.GetWeight();
				m_accumulationBuffer[a] +=
					w * m_staticPressureBuffer[b]; // a <- b
				m_accumulationBuffer[b] +=
					w * m_staticPressureBuffer[a]; // b <- a
			}
		}
		for (int32 i = 0; i < m_count; i++)
		{
			float32 w = m_weightBuffer[i];
			if (m_flagsBuffer.data[i] & b2_staticPressureParticle)
			{
				float32 wh = m_accumulationBuffer[i];
				float32 h =
					(wh + pressurePerWeight * (w - b2_minParticleWeight)) /
					(w + relaxation);
				m_staticPressureBuffer[i] = b2Clamp(h, 0.0f, maxPressure);
			}
			else
			{
				m_staticPressureBuffer[i] = 0;
			}
		}
	}
}

void b2ParticleSystem::SolvePressure(const b2TimeStep& step)
{
	// calculates pressure as a linear function of density
	float32 criticalPressure = GetCriticalPressure(step);
	float32 pressurePerWeight = m_def.pressureStrength * criticalPressure;
	float32 maxPressure = b2_maxParticlePressure * criticalPressure;
	for (int32 i = 0; i < m_count; i++)
	{
		float32 w = m_weightBuffer[i];
		float32 h = pressurePerWeight * b2Max(0.0f, w - b2_minParticleWeight);
		m_accumulationBuffer[i] = b2Min(h, maxPressure);
	}
	// ignores particles which have their own repulsive force
	if (m_allParticleFlags & k_noPressureFlags)
	{
		for (int32 i = 0; i < m_count; i++)
		{
			if (m_flagsBuffer.data[i] & k_noPressureFlags)
			{
				m_accumulationBuffer[i] = 0;
			}
		}
	}
	// static pressure
	if (m_allParticleFlags & b2_staticPressureParticle)
	{
		b2Assert(m_staticPressureBuffer);
		for (int32 i = 0; i < m_count; i++)
		{
			if (m_flagsBuffer.data[i] & b2_staticPressureParticle)
			{
				m_accumulationBuffer[i] += m_staticPressureBuffer[i];
			}
		}
	}
	// applies pressure between each particles in contact
	float32 velocityPerPressure = step.dt / (m_def.density * m_particleDiameter);
	for (int32 k = 0; k < m_bodyContactBuffer.GetCount(); k++)
	{
		const b2ParticleBodyContact& contact = m_bodyContactBuffer[k];
		int32 a = contact.index;
		b2Body* b = contact.body;
		float32 w = contact.weight;
		float32 m = contact.mass;
		b2Vec2 n = contact.normal;
		b2Vec2 p = m_positionBuffer.data[a];
		float32 h = m_accumulationBuffer[a] + pressurePerWeight * w;
		b2Vec2 f = velocityPerPressure * w * m * h * n;
		m_velocityBuffer.data[a] -= GetParticleInvMass() * f;
		b->ApplyLinearImpulse(f, p, true);
	}
	for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
	{
		const b2ParticleContact& contact = m_contactBuffer[k];
		int32 a = contact.GetIndexA();
		int32 b = contact.GetIndexB();
		float32 w = contact.GetWeight();
		b2Vec2 n = contact.GetNormal();
		float32 h = m_accumulationBuffer[a] + m_accumulationBuffer[b];
		b2Vec2 f = velocityPerPressure * w * h * n;
		m_velocityBuffer.data[a] -= f;
		m_velocityBuffer.data[b] += f;
	}
}

void b2ParticleSystem::SolveDamping(const b2TimeStep& step)
{
	// reduces normal velocity of each contact
	float32 linearDamping = m_def.dampingStrength;
	float32 quadraticDamping = 1 / GetCriticalVelocity(step);
	for (int32 k = 0; k < m_bodyContactBuffer.GetCount(); k++)
	{
		const b2ParticleBodyContact& contact = m_bodyContactBuffer[k];
		int32 a = contact.index;
		b2Body* b = contact.body;
		float32 w = contact.weight;
		float32 m = contact.mass;
		b2Vec2 n = contact.normal;
		b2Vec2 p = m_positionBuffer.data[a];
		b2Vec2 v = b->GetLinearVelocityFromWorldPoint(p) -
				   m_velocityBuffer.data[a];
		float32 vn = b2Dot(v, n);
		if (vn < 0)
		{
			float32 damping =
				b2Max(linearDamping * w, b2Min(- quadraticDamping * vn, 0.5f));
			b2Vec2 f = damping * m * vn * n;
			m_velocityBuffer.data[a] += GetParticleInvMass() * f;
			b->ApplyLinearImpulse(-f, p, true);
		}
	}
	for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
	{
		const b2ParticleContact& contact = m_contactBuffer[k];
		int32 a = contact.GetIndexA();
		int32 b = contact.GetIndexB();
		float32 w = contact.GetWeight();
		b2Vec2 n = contact.GetNormal();
		b2Vec2 v = m_velocityBuffer.data[b] - m_velocityBuffer.data[a];
		float32 vn = b2Dot(v, n);
		if (vn < 0)
		{
			float32 damping =
				b2Max(linearDamping * w, b2Min(- quadraticDamping * vn, 0.5f));
			b2Vec2 f = damping * vn * n;
			m_velocityBuffer.data[a] += f;
			m_velocityBuffer.data[b] -= f;
		}
	}
}

inline bool b2ParticleSystem::IsRigidGroup(b2ParticleGroup *group) const
{
	return group && (group->m_groupFlags & b2_rigidParticleGroup);
}

inline b2Vec2 b2ParticleSystem::GetLinearVelocity(
	b2ParticleGroup *group, int32 particleIndex,
	const b2Vec2 &point) const
{
	if (IsRigidGroup(group))
	{
		return group->GetLinearVelocityFromWorldPoint(point);
	}
	else
	{
		return m_velocityBuffer.data[particleIndex];
	}
}

inline void b2ParticleSystem::InitDampingParameter(
	float32* invMass, float32* invInertia, float32* tangentDistance,
	float32 mass, float32 inertia, const b2Vec2& center,
	const b2Vec2& point, const b2Vec2& normal) const
{
	*invMass = mass > 0 ? 1 / mass : 0;
	*invInertia = inertia > 0 ? 1 / inertia : 0;
	*tangentDistance = b2Cross(point - center, normal);
}

inline void b2ParticleSystem::InitDampingParameterWithRigidGroupOrParticle(
	float32* invMass, float32* invInertia, float32* tangentDistance,
	bool isRigidGroup, b2ParticleGroup* group, int32 particleIndex,
	const b2Vec2& point, const b2Vec2& normal) const
{
	if (isRigidGroup)
	{
		InitDampingParameter(
			invMass, invInertia, tangentDistance,
			group->GetMass(), group->GetInertia(), group->GetCenter(),
			point, normal);
	}
	else
	{
		uint32 flags = m_flagsBuffer.data[particleIndex];
		InitDampingParameter(
			invMass, invInertia, tangentDistance,
			flags & b2_wallParticle ? 0 : GetParticleMass(), 0, point,
			point, normal);
	}
}

inline float32 b2ParticleSystem::ComputeDampingImpulse(
	float32 invMassA, float32 invInertiaA, float32 tangentDistanceA,
	float32 invMassB, float32 invInertiaB, float32 tangentDistanceB,
	float32 normalVelocity) const
{
	float32 invMass =
		invMassA + invInertiaA * tangentDistanceA * tangentDistanceA +
		invMassB + invInertiaB * tangentDistanceB * tangentDistanceB;
	return invMass > 0 ? normalVelocity / invMass : 0;
}

inline void b2ParticleSystem::ApplyDamping(
	float32 invMass, float32 invInertia, float32 tangentDistance,
	bool isRigidGroup, b2ParticleGroup* group, int32 particleIndex,
	float32 impulse, const b2Vec2& normal)
{
	if (isRigidGroup)
	{
		group->m_linearVelocity += impulse * invMass * normal;
		group->m_angularVelocity += impulse * tangentDistance * invInertia;
	}
	else
	{
		m_velocityBuffer.data[particleIndex] += impulse * invMass * normal;
	}
}

void b2ParticleSystem::SolveRigidDamping()
{
	// Apply impulse to rigid particle groups colliding with other objects
	// to reduce relative velocity at the colliding point.
	float32 damping = m_def.dampingStrength;
	for (int32 k = 0; k < m_bodyContactBuffer.GetCount(); k++)
	{
		const b2ParticleBodyContact& contact = m_bodyContactBuffer[k];
		int32 a = contact.index;
		b2ParticleGroup* aGroup = m_groupBuffer[a];
		if (IsRigidGroup(aGroup))
		{
			b2Body* b = contact.body;
			b2Vec2 n = contact.normal;
			float32 w = contact.weight;
			b2Vec2 p = m_positionBuffer.data[a];
			b2Vec2 v = b->GetLinearVelocityFromWorldPoint(p) -
					   aGroup->GetLinearVelocityFromWorldPoint(p);
			float32 vn = b2Dot(v, n);
			if (vn < 0)
			// The group's average velocity at particle position 'p' is pushing
			// the particle into the body.
			{
				float32 invMassA, invInertiaA, tangentDistanceA;
				float32 invMassB, invInertiaB, tangentDistanceB;
				InitDampingParameterWithRigidGroupOrParticle(
					&invMassA, &invInertiaA, &tangentDistanceA,
					true, aGroup, a, p, n);
				InitDampingParameter(
					&invMassB, &invInertiaB, &tangentDistanceB,
					b->GetMass(),
					// Calculate b->m_I from public functions of b2Body.
					b->GetInertia() -
							b->GetMass() * b->GetLocalCenter().LengthSquared(),
					b->GetWorldCenter(),
					p, n);
				float32 f = damping * b2Min(w, 1.0f) * ComputeDampingImpulse(
					invMassA, invInertiaA, tangentDistanceA,
					invMassB, invInertiaB, tangentDistanceB,
					vn);
				ApplyDamping(
					invMassA, invInertiaA, tangentDistanceA,
					true, aGroup, a, f, n);
				b->ApplyLinearImpulse(-f * n, p, true);
			}
		}
	}
	for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
	{
		const b2ParticleContact& contact = m_contactBuffer[k];
		int32 a = contact.GetIndexA();
		int32 b = contact.GetIndexB();
		b2Vec2 n = contact.GetNormal();
		float32 w = contact.GetWeight();
		b2ParticleGroup* aGroup = m_groupBuffer[a];
		b2ParticleGroup* bGroup = m_groupBuffer[b];
		bool aRigid = IsRigidGroup(aGroup);
		bool bRigid = IsRigidGroup(bGroup);
		if (aGroup != bGroup && (aRigid || bRigid))
		{
			b2Vec2 p =
				0.5f * (m_positionBuffer.data[a] + m_positionBuffer.data[b]);
			b2Vec2 v =
				GetLinearVelocity(bGroup, b, p) -
				GetLinearVelocity(aGroup, a, p);
			float32 vn = b2Dot(v, n);
			if (vn < 0)
			{
				float32 invMassA, invInertiaA, tangentDistanceA;
				float32 invMassB, invInertiaB, tangentDistanceB;
				InitDampingParameterWithRigidGroupOrParticle(
					&invMassA, &invInertiaA, &tangentDistanceA,
					aRigid, aGroup, a,
					p, n);
				InitDampingParameterWithRigidGroupOrParticle(
					&invMassB, &invInertiaB, &tangentDistanceB,
					bRigid, bGroup, b,
					p, n);
				float32 f = damping * w * ComputeDampingImpulse(
					invMassA, invInertiaA, tangentDistanceA,
					invMassB, invInertiaB, tangentDistanceB,
					vn);
				ApplyDamping(
					invMassA, invInertiaA, tangentDistanceA,
					aRigid, aGroup, a, f, n);
				ApplyDamping(
					invMassB, invInertiaB, tangentDistanceB,
					bRigid, bGroup, b, -f, n);
			}
		}
	}
}

void b2ParticleSystem::SolveExtraDamping()
{
	// Applies additional damping force between bodies and particles which can
	// produce strong repulsive force. Applying damping force multiple times
	// is effective in suppressing vibration.
	for (int32 k = 0; k < m_bodyContactBuffer.GetCount(); k++)
	{
		const b2ParticleBodyContact& contact = m_bodyContactBuffer[k];
		int32 a = contact.index;
		if (m_flagsBuffer.data[a] & k_extraDampingFlags)
		{
			b2Body* b = contact.body;
			float32 m = contact.mass;
			b2Vec2 n = contact.normal;
			b2Vec2 p = m_positionBuffer.data[a];
			b2Vec2 v =
				b->GetLinearVelocityFromWorldPoint(p) -
				m_velocityBuffer.data[a];
			float32 vn = b2Dot(v, n);
			if (vn < 0)
			{
				b2Vec2 f = 0.5f * m * vn * n;
				m_velocityBuffer.data[a] += GetParticleInvMass() * f;
				b->ApplyLinearImpulse(-f, p, true);
			}
		}
	}
}

void b2ParticleSystem::SolveWall()
{
	for (int32 i = 0; i < m_count; i++)
	{
		if (m_flagsBuffer.data[i] & b2_wallParticle)
		{
			m_velocityBuffer.data[i].SetZero();
		}
	}
}

void b2ParticleSystem::SolveRigid(const b2TimeStep& step)
{
	for (b2ParticleGroup* group = m_groupList; group; group = group->GetNext())
	{
		if (group->m_groupFlags & b2_rigidParticleGroup)
		{
			group->UpdateStatistics();
			b2Rot rotation(step.dt * group->m_angularVelocity);
			b2Transform transform(
				group->m_center + step.dt * group->m_linearVelocity -
				b2Mul(rotation, group->m_center), rotation);
			group->m_transform = b2Mul(transform, group->m_transform);
			b2Transform velocityTransform;
			velocityTransform.p.x = step.inv_dt * transform.p.x;
			velocityTransform.p.y = step.inv_dt * transform.p.y;
			velocityTransform.q.s = step.inv_dt * transform.q.s;
			velocityTransform.q.c = step.inv_dt * (transform.q.c - 1);
			for (int32 i = group->m_firstIndex; i < group->m_lastIndex; i++)
			{
				m_velocityBuffer.data[i] = b2Mul(velocityTransform,
												 m_positionBuffer.data[i]);
			}
		}
	}
}

void b2ParticleSystem::SolveElastic(const b2TimeStep& step)
{
	float32 elasticStrength = step.inv_dt * m_def.elasticStrength;
	for (int32 k = 0; k < m_triadBuffer.GetCount(); k++)
	{
		const b2ParticleTriad& triad = m_triadBuffer[k];
		if (triad.flags & b2_elasticParticle)
		{
			int32 a = triad.indexA;
			int32 b = triad.indexB;
			int32 c = triad.indexC;
			const b2Vec2& oa = triad.pa;
			const b2Vec2& ob = triad.pb;
			const b2Vec2& oc = triad.pc;
			b2Vec2 pa = m_positionBuffer.data[a];
			b2Vec2 pb = m_positionBuffer.data[b];
			b2Vec2 pc = m_positionBuffer.data[c];
			b2Vec2& va = m_velocityBuffer.data[a];
			b2Vec2& vb = m_velocityBuffer.data[b];
			b2Vec2& vc = m_velocityBuffer.data[c];
			pa += step.dt * va;
			pb += step.dt * vb;
			pc += step.dt * vc;
			b2Vec2 midPoint = (float32) 1 / 3 * (pa + pb + pc);
			pa -= midPoint;
			pb -= midPoint;
			pc -= midPoint;
			b2Rot r;
			r.s = b2Cross(oa, pa) + b2Cross(ob, pb) + b2Cross(oc, pc);
			r.c = b2Dot(oa, pa) + b2Dot(ob, pb) + b2Dot(oc, pc);
			float32 r2 = r.s * r.s + r.c * r.c;
			float32 invR = b2InvSqrt(r2);
			r.s *= invR;
			r.c *= invR;
			float32 strength = elasticStrength * triad.strength;
			va += strength * (b2Mul(r, oa) - pa);
			vb += strength * (b2Mul(r, ob) - pb);
			vc += strength * (b2Mul(r, oc) - pc);
		}
	}
}

void b2ParticleSystem::SolveSpring(const b2TimeStep& step)
{
	float32 springStrength = step.inv_dt * m_def.springStrength;
	for (int32 k = 0; k < m_pairBuffer.GetCount(); k++)
	{
		const b2ParticlePair& pair = m_pairBuffer[k];
		if (pair.flags & b2_springParticle)
		{
			int32 a = pair.indexA;
			int32 b = pair.indexB;
			b2Vec2 pa = m_positionBuffer.data[a];
			b2Vec2 pb = m_positionBuffer.data[b];
			b2Vec2& va = m_velocityBuffer.data[a];
			b2Vec2& vb = m_velocityBuffer.data[b];
			pa += step.dt * va;
			pb += step.dt * vb;
			b2Vec2 d = pb - pa;
			float32 r0 = pair.distance;
			float32 r1 = d.Length();
			float32 strength = springStrength * pair.strength;
			b2Vec2 f = strength * (r0 - r1) / r1 * d;
			va -= f;
			vb += f;
		}
	}
}

void b2ParticleSystem::SolveTensile(const b2TimeStep& step)
{
	b2Assert(m_accumulation2Buffer);
	for (int32 i = 0; i < m_count; i++)
	{
		m_accumulation2Buffer[i] = b2Vec2_zero;
	}
	for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
	{
		const b2ParticleContact& contact = m_contactBuffer[k];
		if (contact.GetFlags() & b2_tensileParticle)
		{
			int32 a = contact.GetIndexA();
			int32 b = contact.GetIndexB();
			float32 w = contact.GetWeight();
			b2Vec2 n = contact.GetNormal();
			b2Vec2 weightedNormal = (1 - w) * w * n;
			m_accumulation2Buffer[a] -= weightedNormal;
			m_accumulation2Buffer[b] += weightedNormal;
		}
	}
	float32 criticalVelocity = GetCriticalVelocity(step);
	float32 pressureStrength = m_def.surfaceTensionPressureStrength
							 * criticalVelocity;
	float32 normalStrength = m_def.surfaceTensionNormalStrength
						   * criticalVelocity;
	float32 maxVelocityVariation = b2_maxParticleForce * criticalVelocity;
	for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
	{
		const b2ParticleContact& contact = m_contactBuffer[k];
		if (contact.GetFlags() & b2_tensileParticle)
		{
			int32 a = contact.GetIndexA();
			int32 b = contact.GetIndexB();
			float32 w = contact.GetWeight();
			b2Vec2 n = contact.GetNormal();
			float32 h = m_weightBuffer[a] + m_weightBuffer[b];
			b2Vec2 s = m_accumulation2Buffer[b] - m_accumulation2Buffer[a];
			float32 fn = b2Min(
					pressureStrength * (h - 2) + normalStrength * b2Dot(s, n),
					maxVelocityVariation) * w;
			b2Vec2 f = fn * n;
			m_velocityBuffer.data[a] -= f;
			m_velocityBuffer.data[b] += f;
		}
	}
}

void b2ParticleSystem::SolveViscous()
{
	float32 viscousStrength = m_def.viscousStrength;
	for (int32 k = 0; k < m_bodyContactBuffer.GetCount(); k++)
	{
		const b2ParticleBodyContact& contact = m_bodyContactBuffer[k];
		int32 a = contact.index;
		if (m_flagsBuffer.data[a] & b2_viscousParticle)
		{
			b2Body* b = contact.body;
			float32 w = contact.weight;
			float32 m = contact.mass;
			b2Vec2 p = m_positionBuffer.data[a];
			b2Vec2 v = b->GetLinearVelocityFromWorldPoint(p) -
					   m_velocityBuffer.data[a];
			b2Vec2 f = viscousStrength * m * w * v;
			m_velocityBuffer.data[a] += GetParticleInvMass() * f;
			b->ApplyLinearImpulse(-f, p, true);
		}
	}
	for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
	{
		const b2ParticleContact& contact = m_contactBuffer[k];
		if (contact.GetFlags() & b2_viscousParticle)
		{
			int32 a = contact.GetIndexA();
			int32 b = contact.GetIndexB();
			float32 w = contact.GetWeight();
			b2Vec2 v = m_velocityBuffer.data[b] - m_velocityBuffer.data[a];
			b2Vec2 f = viscousStrength * w * v;
			m_velocityBuffer.data[a] += f;
			m_velocityBuffer.data[b] -= f;
		}
	}
}

void b2ParticleSystem::SolveRepulsive(const b2TimeStep& step)
{
	float32 repulsiveStrength =
		m_def.repulsiveStrength * GetCriticalVelocity(step);
	for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
	{
		const b2ParticleContact& contact = m_contactBuffer[k];
		if (contact.GetFlags() & b2_repulsiveParticle)
		{
			int32 a = contact.GetIndexA();
			int32 b = contact.GetIndexB();
			if (m_groupBuffer[a] != m_groupBuffer[b])
			{
				float32 w = contact.GetWeight();
				b2Vec2 n = contact.GetNormal();
				b2Vec2 f = repulsiveStrength * w * n;
				m_velocityBuffer.data[a] -= f;
				m_velocityBuffer.data[b] += f;
			}
		}
	}
}

void b2ParticleSystem::SolvePowder(const b2TimeStep& step)
{
	float32 powderStrength = m_def.powderStrength * GetCriticalVelocity(step);
	float32 minWeight = 1.0f - b2_particleStride;
	for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
	{
		const b2ParticleContact& contact = m_contactBuffer[k];
		if (contact.GetFlags() & b2_powderParticle)
		{
			float32 w = contact.GetWeight();
			if (w > minWeight)
			{
				int32 a = contact.GetIndexA();
				int32 b = contact.GetIndexB();
				b2Vec2 n = contact.GetNormal();
				b2Vec2 f = powderStrength * (w - minWeight) * n;
				m_velocityBuffer.data[a] -= f;
				m_velocityBuffer.data[b] += f;
			}
		}
	}
}

void b2ParticleSystem::SolveSolid(const b2TimeStep& step)
{
	// applies extra repulsive force from solid particle groups
	b2Assert(m_depthBuffer);
	float32 ejectionStrength = step.inv_dt * m_def.ejectionStrength;
	for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
	{
		const b2ParticleContact& contact = m_contactBuffer[k];
		int32 a = contact.GetIndexA();
		int32 b = contact.GetIndexB();
		if (m_groupBuffer[a] != m_groupBuffer[b])
		{
			float32 w = contact.GetWeight();
			b2Vec2 n = contact.GetNormal();
			float32 h = m_depthBuffer[a] + m_depthBuffer[b];
			b2Vec2 f = ejectionStrength * h * w * n;
			m_velocityBuffer.data[a] -= f;
			m_velocityBuffer.data[b] += f;
		}
	}
}

void b2ParticleSystem::SolveForce(const b2TimeStep& step)
{
	float32 velocityPerForce = step.dt * GetParticleInvMass();
	for (int32 i = 0; i < m_count; i++)
	{
		m_velocityBuffer.data[i] += velocityPerForce * m_forceBuffer[i];
	}
	m_hasForce = false;
}

void b2ParticleSystem::SolveColorMixing()
{
	// mixes color between contacting particles
	b2Assert(m_colorBuffer.data);
	const int32 colorMixing128 = (int32) (128 * m_def.colorMixingStrength);
	if (colorMixing128) {
		for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
		{
			const b2ParticleContact& contact = m_contactBuffer[k];
			int32 a = contact.GetIndexA();
			int32 b = contact.GetIndexB();
			if (m_flagsBuffer.data[a] & m_flagsBuffer.data[b] &
				b2_colorMixingParticle)
			{
				b2ParticleColor& colorA = m_colorBuffer.data[a];
				b2ParticleColor& colorB = m_colorBuffer.data[b];
				// Use the static method to ensure certain compilers inline
				// this correctly.
				b2ParticleColor::MixColors(&colorA, &colorB, colorMixing128);
			}
		}
	}
}

void b2ParticleSystem::SolveZombie()
{
	// removes particles with zombie flag
	int32 newCount = 0;
	int32* newIndices = (int32*) m_world->m_stackAllocator.Allocate(
		sizeof(int32) * m_count);
	uint32 allParticleFlags = 0;
	for (int32 i = 0; i < m_count; i++)
	{
		int32 flags = m_flagsBuffer.data[i];
		if (flags & b2_zombieParticle)
		{
			b2DestructionListener * const destructionListener =
				m_world->m_destructionListener;
			if ((flags & b2_destructionListenerParticle) &&
				destructionListener)
			{
				destructionListener->SayGoodbye(this, i);
			}
			// Destroy particle handle.
			if (m_handleIndexBuffer.data)
			{
				b2ParticleHandle * const handle = m_handleIndexBuffer.data[i];
				if (handle)
				{
					handle->SetIndex(b2_invalidParticleIndex);
					m_handleIndexBuffer.data[i] = NULL;
					m_handleAllocator.Free(handle);
				}
			}
			newIndices[i] = b2_invalidParticleIndex;
		}
		else
		{
			newIndices[i] = newCount;
			if (i != newCount)
			{
				// Update handle to reference new particle index.
				if (m_handleIndexBuffer.data)
				{
					b2ParticleHandle * const handle =
						m_handleIndexBuffer.data[i];
					if (handle) handle->SetIndex(newCount);
					m_handleIndexBuffer.data[newCount] = handle;
				}
				m_flagsBuffer.data[newCount] = m_flagsBuffer.data[i];
				if (m_lastBodyContactStepBuffer.data)
				{
					m_lastBodyContactStepBuffer.data[newCount] =
						m_lastBodyContactStepBuffer.data[i];
				}
				if (m_bodyContactCountBuffer.data)
				{
					m_bodyContactCountBuffer.data[newCount] =
						m_bodyContactCountBuffer.data[i];
				}
				if (m_consecutiveContactStepsBuffer.data)
				{
					m_consecutiveContactStepsBuffer.data[newCount] =
						m_consecutiveContactStepsBuffer.data[i];
				}
				m_positionBuffer.data[newCount] = m_positionBuffer.data[i];
				m_velocityBuffer.data[newCount] = m_velocityBuffer.data[i];
				m_groupBuffer[newCount] = m_groupBuffer[i];
				if (m_hasForce)
				{
					m_forceBuffer[newCount] = m_forceBuffer[i];
				}
				if (m_staticPressureBuffer)
				{
					m_staticPressureBuffer[newCount] =
						m_staticPressureBuffer[i];
				}
				if (m_depthBuffer)
				{
					m_depthBuffer[newCount] = m_depthBuffer[i];
				}
				if (m_colorBuffer.data)
				{
					m_colorBuffer.data[newCount] = m_colorBuffer.data[i];
				}
				if (m_userDataBuffer.data)
				{
					m_userDataBuffer.data[newCount] = m_userDataBuffer.data[i];
				}
				if (m_expirationTimeBuffer.data)
				{
					m_expirationTimeBuffer.data[newCount] =
						m_expirationTimeBuffer.data[i];
				}
			}
			newCount++;
			allParticleFlags |= flags;
		}
	}

	// predicate functions
	struct Test
	{
		static bool IsProxyInvalid(const Proxy& proxy)
		{
			return proxy.index < 0;
		}
		static bool IsContactInvalid(const b2ParticleContact& contact)
		{
			return contact.GetIndexA() < 0 || contact.GetIndexB() < 0;
		}
		static bool IsBodyContactInvalid(const b2ParticleBodyContact& contact)
		{
			return contact.index < 0;
		}
		static bool IsPairInvalid(const b2ParticlePair& pair)
		{
			return pair.indexA < 0 || pair.indexB < 0;
		}
		static bool IsTriadInvalid(const b2ParticleTriad& triad)
		{
			return triad.indexA < 0 || triad.indexB < 0 || triad.indexC < 0;
		}
	};

	// update proxies
	for (int32 k = 0; k < m_proxyBuffer.GetCount(); k++)
	{
		Proxy& proxy = m_proxyBuffer.Begin()[k];
		proxy.index = newIndices[proxy.index];
	}
	m_proxyBuffer.RemoveIf(Test::IsProxyInvalid);

	// update contacts
	for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
	{
		b2ParticleContact& contact = m_contactBuffer[k];
		contact.SetIndices(newIndices[contact.GetIndexA()],
						   newIndices[contact.GetIndexB()]);
	}
	m_contactBuffer.RemoveIf(Test::IsContactInvalid);

	// update particle-body contacts
	for (int32 k = 0; k < m_bodyContactBuffer.GetCount(); k++)
	{
		b2ParticleBodyContact& contact = m_bodyContactBuffer[k];
		contact.index = newIndices[contact.index];
	}
	m_bodyContactBuffer.RemoveIf(Test::IsBodyContactInvalid);

	// update pairs
	for (int32 k = 0; k < m_pairBuffer.GetCount(); k++)
	{
		b2ParticlePair& pair = m_pairBuffer[k];
		pair.indexA = newIndices[pair.indexA];
		pair.indexB = newIndices[pair.indexB];
	}
	m_pairBuffer.RemoveIf(Test::IsPairInvalid);

	// update triads
	for (int32 k = 0; k < m_triadBuffer.GetCount(); k++)
	{
		b2ParticleTriad& triad = m_triadBuffer[k];
		triad.indexA = newIndices[triad.indexA];
		triad.indexB = newIndices[triad.indexB];
		triad.indexC = newIndices[triad.indexC];
	}
	m_triadBuffer.RemoveIf(Test::IsTriadInvalid);

	// Update lifetime indices.
	if (m_indexByExpirationTimeBuffer.data)
	{
		int32 writeOffset = 0;
		for (int32 readOffset = 0; readOffset < m_count; readOffset++)
		{
			const int32 newIndex = newIndices[
				m_indexByExpirationTimeBuffer.data[readOffset]];
			if (newIndex != b2_invalidParticleIndex)
			{
				m_indexByExpirationTimeBuffer.data[writeOffset++] = newIndex;
			}
		}
	}

	// update groups
	for (b2ParticleGroup* group = m_groupList; group; group = group->GetNext())
	{
		int32 firstIndex = newCount;
		int32 lastIndex = 0;
		bool modified = false;
		for (int32 i = group->m_firstIndex; i < group->m_lastIndex; i++)
		{
			int32 j = newIndices[i];
			if (j >= 0) {
				firstIndex = b2Min(firstIndex, j);
				lastIndex = b2Max(lastIndex, j + 1);
			} else {
				modified = true;
			}
		}
		if (firstIndex < lastIndex)
		{
			group->m_firstIndex = firstIndex;
			group->m_lastIndex = lastIndex;
			if (modified)
			{
				if (group->m_groupFlags & b2_solidParticleGroup)
				{
					SetGroupFlags(group,
								  group->m_groupFlags |
								  b2_particleGroupNeedsUpdateDepth);
				}
			}
		}
		else
		{
			group->m_firstIndex = 0;
			group->m_lastIndex = 0;
			if (!(group->m_groupFlags & b2_particleGroupCanBeEmpty))
			{
				SetGroupFlags(group,
					group->m_groupFlags | b2_particleGroupWillBeDestroyed);
			}
		}
	}

	// update particle count
	m_count = newCount;
	m_world->m_stackAllocator.Free(newIndices);
	m_allParticleFlags = allParticleFlags;
	m_needsUpdateAllParticleFlags = false;

	// destroy bodies with no particles
	for (b2ParticleGroup* group = m_groupList; group;)
	{
		b2ParticleGroup* next = group->GetNext();
		if (group->m_groupFlags & b2_particleGroupWillBeDestroyed)
		{
			DestroyParticleGroup(group);
		}
		group = next;
	}
}

/// Destroy all particles which have outlived their lifetimes set by
/// SetParticleLifetime().
void b2ParticleSystem::SolveLifetimes(const b2TimeStep& step)
{
	b2Assert(m_expirationTimeBuffer.data);
	b2Assert(m_indexByExpirationTimeBuffer.data);
	// Update the time elapsed.
	m_timeElapsed = LifetimeToExpirationTime(step.dt);
	// Get the floor (non-fractional component) of the elapsed time.
	const int32 quantizedTimeElapsed = GetQuantizedTimeElapsed();

	const int32* const expirationTimes = m_expirationTimeBuffer.data;
	int32* const expirationTimeIndices = m_indexByExpirationTimeBuffer.data;
	const int32 particleCount = GetParticleCount();
	// Sort the lifetime buffer if it's required.
	if (m_expirationTimeBufferRequiresSorting)
	{
		const ExpirationTimeComparator expirationTimeComparator(
			expirationTimes);
		std::sort(expirationTimeIndices,
				  expirationTimeIndices + particleCount,
				  expirationTimeComparator);
		m_expirationTimeBufferRequiresSorting = false;
	}

	// Destroy particles which have expired.
	for (int32 i = particleCount - 1; i >= 0; --i)
	{
		const int32 particleIndex = expirationTimeIndices[i];
		const int32 expirationTime = expirationTimes[particleIndex];
		// If no particles need to be destroyed, skip this.
		if (quantizedTimeElapsed < expirationTime || expirationTime <= 0)
		{
			break;
		}
		// Destroy this particle.
		DestroyParticle(particleIndex);
	}
}

void b2ParticleSystem::RotateBuffer(int32 start, int32 mid, int32 end)
{
	// move the particles assigned to the given group toward the end of array
	if (start == mid || mid == end)
	{
		return;
	}
	b2Assert(mid >= start && mid <= end);
	struct NewIndices
	{
		int32 operator[](int32 i) const
		{
			if (i < start)
			{
				return i;
			}
			else if (i < mid)
			{
				return i + end - mid;
			}
			else if (i < end)
			{
				return i + start - mid;
			}
			else
			{
				return i;
			}
		}
		int32 start, mid, end;
	} newIndices;
	newIndices.start = start;
	newIndices.mid = mid;
	newIndices.end = end;

	std::rotate(m_flagsBuffer.data + start, m_flagsBuffer.data + mid,
				m_flagsBuffer.data + end);
	if (m_lastBodyContactStepBuffer.data)
	{
		std::rotate(m_lastBodyContactStepBuffer.data + start,
					m_lastBodyContactStepBuffer.data + mid,
					m_lastBodyContactStepBuffer.data + end);
	}
	if (m_bodyContactCountBuffer.data)
	{
		std::rotate(m_bodyContactCountBuffer.data + start,
					m_bodyContactCountBuffer.data + mid,
					m_bodyContactCountBuffer.data + end);
	}
	if (m_consecutiveContactStepsBuffer.data)
	{
		std::rotate(m_consecutiveContactStepsBuffer.data + start,
					m_consecutiveContactStepsBuffer.data + mid,
					m_consecutiveContactStepsBuffer.data + end);
	}
	std::rotate(m_positionBuffer.data + start, m_positionBuffer.data + mid,
				m_positionBuffer.data + end);
	std::rotate(m_velocityBuffer.data + start, m_velocityBuffer.data + mid,
				m_velocityBuffer.data + end);
	std::rotate(m_groupBuffer + start, m_groupBuffer + mid,
				m_groupBuffer + end);
	if (m_hasForce)
	{
		std::rotate(m_forceBuffer + start, m_forceBuffer + mid,
					m_forceBuffer + end);
	}
	if (m_staticPressureBuffer)
	{
		std::rotate(m_staticPressureBuffer + start,
					m_staticPressureBuffer + mid,
					m_staticPressureBuffer + end);
	}
	if (m_depthBuffer)
	{
		std::rotate(m_depthBuffer + start, m_depthBuffer + mid,
					m_depthBuffer + end);
	}
	if (m_colorBuffer.data)
	{
		std::rotate(m_colorBuffer.data + start,
					m_colorBuffer.data + mid, m_colorBuffer.data + end);
	}
	if (m_userDataBuffer.data)
	{
		std::rotate(m_userDataBuffer.data + start,
					m_userDataBuffer.data + mid, m_userDataBuffer.data + end);
	}

	// Update handle indices.
	if (m_handleIndexBuffer.data)
	{
		std::rotate(m_handleIndexBuffer.data + start,
					m_handleIndexBuffer.data + mid,
					m_handleIndexBuffer.data + end);
		for (int32 i = start; i < end; ++i)
		{
			b2ParticleHandle * const handle = m_handleIndexBuffer.data[i];
			if (handle) handle->SetIndex(newIndices[handle->GetIndex()]);
		}
	}

	if (m_expirationTimeBuffer.data)
	{
		std::rotate(m_expirationTimeBuffer.data + start,
					m_expirationTimeBuffer.data + mid,
					m_expirationTimeBuffer.data + end);
		// Update expiration time buffer indices.
		const int32 particleCount = GetParticleCount();
		int32* const indexByExpirationTime =
			m_indexByExpirationTimeBuffer.data;
		for (int32 i = 0; i < particleCount; ++i)
		{
			indexByExpirationTime[i] = newIndices[indexByExpirationTime[i]];
		}
	}

	// update proxies
	for (int32 k = 0; k < m_proxyBuffer.GetCount(); k++)
	{
		Proxy& proxy = m_proxyBuffer.Begin()[k];
		proxy.index = newIndices[proxy.index];
	}

	// update contacts
	for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
	{
		b2ParticleContact& contact = m_contactBuffer[k];
		contact.SetIndices(newIndices[contact.GetIndexA()],
						   newIndices[contact.GetIndexB()]);
	}

	// update particle-body contacts
	for (int32 k = 0; k < m_bodyContactBuffer.GetCount(); k++)
	{
		b2ParticleBodyContact& contact = m_bodyContactBuffer[k];
		contact.index = newIndices[contact.index];
	}

	// update pairs
	for (int32 k = 0; k < m_pairBuffer.GetCount(); k++)
	{
		b2ParticlePair& pair = m_pairBuffer[k];
		pair.indexA = newIndices[pair.indexA];
		pair.indexB = newIndices[pair.indexB];
	}

	// update triads
	for (int32 k = 0; k < m_triadBuffer.GetCount(); k++)
	{
		b2ParticleTriad& triad = m_triadBuffer[k];
		triad.indexA = newIndices[triad.indexA];
		triad.indexB = newIndices[triad.indexB];
		triad.indexC = newIndices[triad.indexC];
	}

	// update groups
	for (b2ParticleGroup* group = m_groupList; group; group = group->GetNext())
	{
		group->m_firstIndex = newIndices[group->m_firstIndex];
		group->m_lastIndex = newIndices[group->m_lastIndex - 1] + 1;
	}
}

/// Set the lifetime (in seconds) of a particle relative to the current
/// time.
void b2ParticleSystem::SetParticleLifetime(const int32 index,
										   const float32 lifetime)
{
	b2Assert(ValidateParticleIndex(index));
	const bool initializeExpirationTimes =
		m_indexByExpirationTimeBuffer.data == NULL;
	m_expirationTimeBuffer.data = RequestBuffer(
		m_expirationTimeBuffer.data);
	m_indexByExpirationTimeBuffer.data = RequestBuffer(
		m_indexByExpirationTimeBuffer.data);

	// Initialize the inverse mapping buffer.
	if (initializeExpirationTimes)
	{
		const int32 particleCount = GetParticleCount();
		for (int32 i = 0; i < particleCount; ++i)
		{
			m_indexByExpirationTimeBuffer.data[i] = i;
		}
	}
	const int32 quantizedLifetime = (int32)(lifetime /
											m_def.lifetimeGranularity);
	// Use a negative lifetime so that it's possible to track which
	// of the infinite lifetime particles are older.
	const int32 newExpirationTime = quantizedLifetime > 0 ?
		GetQuantizedTimeElapsed() + quantizedLifetime : quantizedLifetime;
	if (newExpirationTime != m_expirationTimeBuffer.data[index])
	{
		m_expirationTimeBuffer.data[index] = newExpirationTime;
		m_expirationTimeBufferRequiresSorting = true;
	}
}


/// Convert a lifetime value in returned by GetExpirationTimeBuffer()
/// to a value in seconds relative to the current simulation time.
float32 b2ParticleSystem::ExpirationTimeToLifetime(
	const int32 expirationTime) const
{
	return (float32)(expirationTime > 0 ?
					 	expirationTime - GetQuantizedTimeElapsed() :
					 	expirationTime) * m_def.lifetimeGranularity;
}

/// Get the lifetime (in seconds) of a particle relative to the current
/// time.
float32 b2ParticleSystem::GetParticleLifetime(const int32 index)
{
	b2Assert(ValidateParticleIndex(index));
	return ExpirationTimeToLifetime(GetExpirationTimeBuffer()[index]);
}

/// Get the array of particle lifetimes indexed by particle index.
/// GetParticleCount() items are in the returned array.
const int32* b2ParticleSystem::GetExpirationTimeBuffer()
{
	m_expirationTimeBuffer.data = RequestBuffer(
		m_expirationTimeBuffer.data);
	return m_expirationTimeBuffer.data;
}

/// Get the array of particle indices ordered by lifetime.
/// GetExpirationTimeBuffer(
///    GetIndexByExpirationTimeBuffer()[index])
/// is equivalent to GetParticleLifetime(index).
/// GetParticleCount() items are in the returned array.
const int32* b2ParticleSystem::GetIndexByExpirationTimeBuffer()
{
	// If particles are present, initialize / reinitialize the lifetime buffer.
	if (GetParticleCount())
	{
		SetParticleLifetime(0, GetParticleLifetime(0));
	}
	else
	{
		m_indexByExpirationTimeBuffer.data = RequestBuffer(
			m_indexByExpirationTimeBuffer.data);
	}
	return m_indexByExpirationTimeBuffer.data;
}

void b2ParticleSystem::SetDestructionByAge(const bool enable)
{
	if (enable)
	{
		GetExpirationTimeBuffer();
	}
	m_def.destroyByAge = enable;
}

/// Get the time elapsed in b2ParticleSystemDef::lifetimeGranularity.
int32 b2ParticleSystem::GetQuantizedTimeElapsed() const
{
	return (int32)(m_timeElapsed >> 32);
}

/// Convert a lifetime in seconds to an expiration time.
int64 b2ParticleSystem::LifetimeToExpirationTime(const float32 lifetime) const
{
	return m_timeElapsed + (int64)((lifetime / m_def.lifetimeGranularity) *
								   (float32)(1LL << 32));
}

template <typename T> void b2ParticleSystem::SetUserOverridableBuffer(
	UserOverridableBuffer<T>* buffer, T* newData, int32 newCapacity)
{
	b2Assert((newData && newCapacity) || (!newData && !newCapacity));
	if (!buffer->userSuppliedCapacity && buffer->data)
	{
		m_world->m_blockAllocator.Free(
			buffer->data, sizeof(T) * m_internalAllocatedCapacity);
	}
	buffer->data = newData;
	buffer->userSuppliedCapacity = newCapacity;
}

void b2ParticleSystem::SetFlagsBuffer(uint32* buffer, int32 capacity)
{
	SetUserOverridableBuffer(&m_flagsBuffer, buffer, capacity);
}

void b2ParticleSystem::SetPositionBuffer(b2Vec2* buffer,
												 int32 capacity)
{
	SetUserOverridableBuffer(&m_positionBuffer, buffer, capacity);
}

void b2ParticleSystem::SetVelocityBuffer(b2Vec2* buffer,
												 int32 capacity)
{
	SetUserOverridableBuffer(&m_velocityBuffer, buffer, capacity);
}

void b2ParticleSystem::SetColorBuffer(b2ParticleColor* buffer,
											  int32 capacity)
{
	SetUserOverridableBuffer(&m_colorBuffer, buffer, capacity);
}

void b2ParticleSystem::SetUserDataBuffer(void** buffer, int32 capacity)
{
	SetUserOverridableBuffer(&m_userDataBuffer, buffer, capacity);
}

void b2ParticleSystem::SetParticleFlags(int32 index, uint32 newFlags)
{
	uint32* oldFlags = &m_flagsBuffer.data[index];
	if (*oldFlags & ~newFlags)
	{
		// If any flags might be removed
		m_needsUpdateAllParticleFlags = true;
	}
	if (~m_allParticleFlags & newFlags)
	{
		// If any flags were added
		if (newFlags & b2_tensileParticle)
		{
			m_accumulation2Buffer = RequestBuffer(
				m_accumulation2Buffer);
		}
		if (newFlags & b2_colorMixingParticle)
		{
			m_colorBuffer.data = RequestBuffer(m_colorBuffer.data);
		}
		m_allParticleFlags |= newFlags;
	}
	*oldFlags = newFlags;
}

void b2ParticleSystem::SetGroupFlags(
	b2ParticleGroup* group, uint32 newFlags)
{
	uint32* oldFlags = &group->m_groupFlags;
	if ((*oldFlags ^ newFlags) & b2_solidParticleGroup)
	{
		// If the b2_solidParticleGroup flag changed schedule depth update.
		newFlags |= b2_particleGroupNeedsUpdateDepth;
	}
	if (*oldFlags & ~newFlags)
	{
		// If any flags might be removed
		m_needsUpdateAllGroupFlags = true;
	}
	if (~m_allGroupFlags & newFlags)
	{
		// If any flags were added
		if (newFlags & b2_solidParticleGroup)
		{
			m_depthBuffer = RequestBuffer(m_depthBuffer);
		}
		m_allGroupFlags |= newFlags;
	}
	*oldFlags = newFlags;
}

static inline bool IsSignificantForce(const b2Vec2& force)
{
	return force.x != 0 || force.y != 0;
}

inline bool b2ParticleSystem::ForceCanBeApplied(uint32 flags) const
{
	return !(flags & b2_wallParticle);
}

inline void b2ParticleSystem::PrepareForceBuffer()
{
	if (!m_hasForce)
	{
		memset(m_forceBuffer, 0, sizeof(*m_forceBuffer) * m_count);
		m_hasForce = true;
	}
}

void b2ParticleSystem::ApplyForce(int32 firstIndex, int32 lastIndex,
								  const b2Vec2& force)
{
	// Ensure we're not trying to apply force to particles that can't move,
	// such as wall particles.
#if B2_ASSERT_ENABLED
	uint32 flags = 0;
	for (int32 i = firstIndex; i < lastIndex; i++)
	{
		flags |= m_flagsBuffer.data[i];
	}
	b2Assert(ForceCanBeApplied(flags));
#endif

	// Early out if force does nothing (optimization).
	const b2Vec2 distributedForce = force / (float32)(lastIndex - firstIndex);
	if (IsSignificantForce(distributedForce))
	{
		PrepareForceBuffer();

		// Distribute the force over all the particles.
		for (int32 i = firstIndex; i < lastIndex; i++)
		{
			m_forceBuffer[i] += distributedForce;
		}
	}
}

void b2ParticleSystem::ParticleApplyForce(int32 index, const b2Vec2& force)
{
	if (IsSignificantForce(force) &&
		ForceCanBeApplied(m_flagsBuffer.data[index]))
	{
		PrepareForceBuffer();
		m_forceBuffer[index] += force;
	}
}

void b2ParticleSystem::ApplyLinearImpulse(int32 firstIndex, int32 lastIndex,
										  const b2Vec2& impulse)
{
	const float32 numParticles = (float32)(lastIndex - firstIndex);
	const float32 totalMass = numParticles * GetParticleMass();
	const b2Vec2 velocityDelta = impulse / totalMass;
	for (int32 i = firstIndex; i < lastIndex; i++)
	{
		m_velocityBuffer.data[i] += velocityDelta;
	}
}

void b2ParticleSystem::QueryAABB(b2QueryCallback* callback,
								 const b2AABB& aabb) const
{
	if (m_proxyBuffer.GetCount() == 0)
	{
		return;
	}
	const Proxy* beginProxy = m_proxyBuffer.Begin();
	const Proxy* endProxy = m_proxyBuffer.End();
	const Proxy* firstProxy = std::lower_bound(
		beginProxy, endProxy,
		computeTag(
			m_inverseDiameter * aabb.lowerBound.x,
			m_inverseDiameter * aabb.lowerBound.y));
	const Proxy* lastProxy = std::upper_bound(
		firstProxy, endProxy,
		computeTag(
			m_inverseDiameter * aabb.upperBound.x,
			m_inverseDiameter * aabb.upperBound.y));
	for (const Proxy* proxy = firstProxy; proxy < lastProxy; ++proxy)
	{
		int32 i = proxy->index;
		const b2Vec2& p = m_positionBuffer.data[i];
		if (aabb.lowerBound.x < p.x && p.x < aabb.upperBound.x &&
			aabb.lowerBound.y < p.y && p.y < aabb.upperBound.y)
		{
			if (!callback->ReportParticle(this, i))
			{
				break;
			}
		}
	}
}

void b2ParticleSystem::QueryShapeAABB(b2QueryCallback* callback,
									  const b2Shape& shape,
									  const b2Transform& xf) const
{
	b2AABB aabb;
	shape.ComputeAABB(&aabb, xf, 0);
	QueryAABB(callback, aabb);
}

void b2ParticleSystem::RayCast(b2RayCastCallback* callback,
							   const b2Vec2& point1,
							   const b2Vec2& point2) const
{
	if (m_proxyBuffer.GetCount() == 0)
	{
		return;
	}
	b2AABB aabb;
	aabb.lowerBound = b2Min(point1, point2);
	aabb.upperBound = b2Max(point1, point2);
	float32 fraction = 1;
	// solving the following equation:
	// ((1-t)*point1+t*point2-position)^2=diameter^2
	// where t is a potential fraction
	b2Vec2 v = point2 - point1;
	float32 v2 = b2Dot(v, v);
	InsideBoundsEnumerator enumerator = GetInsideBoundsEnumerator(aabb);
	int32 i;
	while ((i = enumerator.GetNext()) >= 0)
	{
		b2Vec2 p = point1 - m_positionBuffer.data[i];
		float32 pv = b2Dot(p, v);
		float32 p2 = b2Dot(p, p);
		float32 determinant = pv * pv - v2 * (p2 - m_squaredDiameter);
		if (determinant >= 0)
		{
			float32 sqrtDeterminant = b2Sqrt(determinant);
			// find a solution between 0 and fraction
			float32 t = (-pv - sqrtDeterminant) / v2;
			if (t > fraction)
			{
				continue;
			}
			if (t < 0)
			{
				t = (-pv + sqrtDeterminant) / v2;
				if (t < 0 || t > fraction)
				{
					continue;
				}
			}
			b2Vec2 n = p + t * v;
			n.Normalize();
			float32 f = callback->ReportParticle(this, i, point1 + t * v, n, t);
			fraction = b2Min(fraction, f);
			if (fraction <= 0)
			{
				break;
			}
		}
	}
}

float32 b2ParticleSystem::ComputeCollisionEnergy() const
{
	float32 sum_v2 = 0;
	for (int32 k = 0; k < m_contactBuffer.GetCount(); k++)
	{
		const b2ParticleContact& contact = m_contactBuffer[k];
		int32 a = contact.GetIndexA();
		int32 b = contact.GetIndexB();
		b2Vec2 n = contact.GetNormal();
		b2Vec2 v = m_velocityBuffer.data[b] - m_velocityBuffer.data[a];
		float32 vn = b2Dot(v, n);
		if (vn < 0)
		{
			sum_v2 += vn * vn;
		}
	}
	return 0.5f * GetParticleMass() * sum_v2;
}

void b2ParticleSystem::SetStuckThreshold(int32 steps)
{
	m_stuckThreshold = steps;

	if (steps > 0)
	{
		m_lastBodyContactStepBuffer.data = RequestBuffer(
			m_lastBodyContactStepBuffer.data);
		m_bodyContactCountBuffer.data = RequestBuffer(
			m_bodyContactCountBuffer.data);
		m_consecutiveContactStepsBuffer.data = RequestBuffer(
			m_consecutiveContactStepsBuffer.data);
	}
}

#if LIQUIDFUN_EXTERNAL_LANGUAGE_API

b2ParticleSystem::b2ExceptionType b2ParticleSystem::IsBufCopyValid(
	int startIndex, int numParticles, int copySize, int bufSize) const
{
	const int maxNumParticles = GetParticleCount();

	// are we actually copying?
	if (copySize == 0)
	{
		return b2_noExceptions;
	}

	// is the index out of bounds?
	if (startIndex < 0 ||
		startIndex >= maxNumParticles ||
		numParticles < 0 ||
		numParticles + startIndex > maxNumParticles)
	{
		return b2_particleIndexOutOfBounds;
	}

	// are we copying within the boundaries?
	if (copySize > bufSize)
	{
		return b2_bufferTooSmall;
	}

	return b2_noExceptions;
}

#endif // LIQUIDFUN_EXTERNAL_LANGUAGE_API

// end of ParticleSystem.cpp

#endif
