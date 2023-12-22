
#ifndef CLMULK_H
#define CLMULK_H

//************************************************
//

#if !defined(__x86_64__) && !defined(_MSC_VER)
static inline uint32_t clmulk_ctz(uint64_t x) { return (x!=0) ? (uint32_t)__builtin_ctzl(x) : 64; }
static inline uint32_t clmulk_pop(uint64_t x) { return (uint32_t)__builtin_popcountl(x); }

#else

#if defined(_MSC_VER)
#include <wmmintrin.h>
#else
#include <x86intrin.h>
#endif

static inline uint32_t clmulk_ctz(uint64_t x) { return (uint32_t)_tzcnt_u64(x); }
static inline uint32_t clmulk_clz(uint64_t x) { return (uint32_t)__lzcnt64(x);  }
static inline uint32_t clmulk_pop(uint64_t x) { return (uint32_t)_mm_popcnt_u64(x); }
#endif

static inline uint32_t clmulk_log2(uint64_t x)      { return (63-clmulk_clz(x  )); }
static inline uint32_t clmulk_log2_ceil(uint64_t x) { return (64-clmulk_clz(x-1)); }
static inline uint32_t clmulk_parity(uint64_t x)    { return clmulk_pop(x) & 1; }

typedef struct { uint64_t q, r;} clmulk_pair_t; 


//************************************************
// tiny register machine 
//
// * only two operations:
//   1) r_{n+1} = r_a ^ (r_b << s)
//   2) return r_n << s
// * destintation register is implict. each is
//   is only assigned to once.
// * this is more general than I'm currently using.
//   it would allow for additive factoring. As of now
//   there's only:
//   0) returns & explicit short cuts
//      * 
//      * 
//      * 
//      * pop(k) >= 4 
//   1) multiply by 2-bit number: (1+2^a)
//      r[rn+1] = r[rn] ^ (r[rn  ] << a); rn++;
//   2) multiply by 3-bit number: (1+2^a+2^b)
//      r[rn+1] = r[rn] ^ (r[rn  ] << a); rn++;
//      r[rn+1] = r[rn] ^ (r[rn-1] << b); rn++;
//   3) multiply by 2-bit number and add 1
//      2^s(1+2^a)+1
//      r[rn+1] = r[rn] ^ (r[rn  ] << a); rn++;
//      r[rn+1] = r[ 0] ^ (r[rn  ] << s); rn++;
//   4) multiply by 2-bit number and add 2^s
//      

// warning this could be passed around
typedef struct {
  uint8_t op;    // only low bit used ATM
  uint8_t a;     // register of add
  uint8_t b;     // register of shifted
  uint8_t s;     // [0,63]
} clmulk_op_t;

// enough space for schoolbook. can do
// better than this if (as an example)
// encode via (~k + ~0) support is added for
// large popcount vs sig digits...so low
// transition counts.
#define CLMULK_MAXOPS 65

// quick wrapper
typedef struct {
  uint64_t    k;
  clmulk_op_t op[CLMULK_MAXOPS];
} clmulk_vm_t;


#define CLMULK_DEF(OP,A,B,S) (clmulk_op_t){.op=(uint8_t)(OP),.a=(uint8_t)(A),.b=(uint8_t)(B),.s=(uint8_t)(S)}

#define CLMULK_INVALID(A,B,S)  CLMULK_DEF(0xff,A,B,S)
#define CLMULK_MUL(A,B,S)      CLMULK_DEF(1,A,B,S)
#define CLMULK_ADD(A,B)        CLMULK_DEF(1,A,B,0)
#define CLMULK_RET             CLMULK_DEF(0,0,0,0)
#define CLMULK_SHIFT_RET(S)    CLMULK_DEF(0,0,0,S)

// 
typedef struct {
  // change to uint32_t? probably
  uint8_t rn;
  uint8_t limit;
  clmulk_op_t* op;
} clmulk_frag_t;


// should hold memorization info
typedef struct {
  uint32_t foo;
} clmulk_env_t;

// strategy level. should this be different than builder?
// has two working fragments to be able to swap when
// a lower cost sequence has been found.
typedef struct {
  uint8_t limit;
  clmulk_frag_t* f[2];  // f[0]: best, f[1]: current working
} clmulk_build_t;

// should be user interface to creating a constant
typedef struct {
  clmulk_env_t* env;
  uint64_t      k;

  // not really..humm...
  clmulk_op_t   op1[CLMULK_MAXOPS];
  clmulk_op_t   op2[CLMULK_MAXOPS];
} clmulk_builder_t;

extern void clmulk_print_c(uint64_t k, clmulk_op_t* op);

extern uint64_t clmulk_eval(clmulk_op_t* op, uint64_t x);

extern uint32_t clmulk_g2br1(clmulk_frag_t* f, uint64_t k);

extern uint64_t clmulk_sanity(clmulk_op_t* op, const uint64_t k);

#endif
