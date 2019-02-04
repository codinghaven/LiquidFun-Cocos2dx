#ifndef LIQUID_FUN_H_DEFINE
#define LIQUID_FUN_H_DEFINE

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
#endif
