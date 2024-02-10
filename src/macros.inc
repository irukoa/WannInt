#ifdef _OPENMP
#define _OMPTGT_(x) $OMP x
#else
#define _OMPTGT_(x) disabled
#endif
