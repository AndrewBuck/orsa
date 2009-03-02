#include <orsaTBB/malloc.h>

#ifdef ORSA_USE_TBB

// #include <iostream>

void * operator new(std::size_t size) throw(std::bad_alloc) {
  // fprintf(stderr,"--MARK-- (void * operator new(std::size_t size))\n");
  if (size == 0) size = 1;
  if (void * ptr = scalable_malloc(size)) {
    return ptr;
  }
  throw std::bad_alloc();
}

void * operator new(std::size_t size, const std::nothrow_t &) throw() {
  // fprintf(stderr,"--MARK-- (void * operator new(std::size_t size, const std::nothrow_t &))\n");
  if (size == 0) size = 1;
  if (void * ptr = scalable_malloc(size)) {
    return ptr;
  }
  return NULL;
}

void * operator new[](std::size_t size) throw(std::bad_alloc) {
  // fprintf(stderr,"--MARK-- (void * operator new[](std::size_t size))\n");
  return operator new(size);
}

void * operator new[](std::size_t size, const std::nothrow_t &) throw() {
  // fprintf(stderr,"--MARK-- (void * operator new[](std::size_t size, const std::nothrow_t &))\n");
  return operator new(size, std::nothrow);
}

void operator delete(void * ptr) throw() {
  // fprintf(stderr,"--MARK-- (void operator delete(void * ptr))\n");
  if (ptr != 0) scalable_free(ptr);
}

void operator delete(void * ptr, const std::nothrow_t &) throw() {
  // fprintf(stderr,"--MARK-- (void operator delete(void * ptr, const std::nothrow_t &))\n");
  if (ptr != 0) scalable_free(ptr);
}

void operator delete[](void * ptr) throw() {
  // fprintf(stderr,"--MARK-- (void operator delete[](void * ptr))\n");
  operator delete(ptr);
}

void operator delete[](void * ptr, const std::nothrow_t &) throw() {
  // fprintf(stderr,"--MARK-- (void operator delete[](void * ptr, const std::nothrow_t &))\n");
  operator delete(ptr, std::nothrow);
}

#endif // ORSA_USE_TBB
