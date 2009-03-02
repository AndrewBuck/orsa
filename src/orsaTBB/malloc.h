#ifndef _ORSA_TBB_MALLOC_
#define _ORSA_TBB_MALLOC_

#ifdef ORSA_USE_TBB

// use this to check that every compiled file includes this header
/* 
   #warning ***************************************************
   #warning ***************************************************
   #warning ******** FILE OK INCLUDES TBB MALLOC.H ************
   #warning ***************************************************
   #warning ***************************************************
*/

#include <tbb/scalable_allocator.h>
#include <tbb/cache_aligned_allocator.h>

void * operator new(std::size_t size) throw(std::bad_alloc);
void * operator new(std::size_t size, const std::nothrow_t &) throw();
void * operator new[](std::size_t size) throw(std::bad_alloc);
void * operator new[](std::size_t size, const std::nothrow_t &) throw();

void operator delete(void * ptr) throw();
void operator delete(void * ptr, const std::nothrow_t &) throw();
void operator delete[](void * ptr) throw();
void operator delete[](void * ptr, const std::nothrow_t &) throw();

#endif // ORSA_USE_TBB

#endif // _ORSA_TBB_MALLOC_
