// Marc B. Reynolds, 2023
// Public Domain under http://unlicense.org, see link for details.


#include <stdint.h>
//#include <smmintrin.h>
// 
#if defined(_MSC_VER)
#include <wmmintrin.h>
#else
#include <x86intrin.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

/************************************************

  s-significant bit number:
  * even parity has (at least one) factor of 3
    more generally there's some 2-bit number that is
    a divisor.
  * 2-bit  irreducible numbers have form 1 {0}_n 1
    where 'n' is an even number (including zero).
    first few are 0x3,0x9,0x21,0x81..

 ************************************************/

// * initial cost limit should be pop(k)-1 & not pop(k)
// * check on: ceil(log2(transition_count)) as low bound
//   estimate. this is ad-hoc based on my current thinking
//   of the "space" that's explored.

#define CARRYLESS_IMPLEMENTATION
#define INTOPS_IMPLEMENTATION

#include "carryless.h"
#include "intops.h"
#include "bitops.h"

#include "clmulk.h"


// just for print_c: usually we know it...but don't care atm
uint32_t clmulk_len(clmulk_op_t* op)
{
  uint32_t i = 0;

  while(op[i].op != 0) i++;

  return i;
}

void clmulk_print_c(uint64_t k, clmulk_op_t* op)
{
  uint32_t len = clmulk_len(op);
  
  printf("uint64_t cl_mulk_%lx(uint64_t x)\n{\n  uint64_t r[%u];\n\n  r[ 0] = x;\n", k,len+1);

  for(uint32_t i=0; i<CLMULK_VM_MAXOPS; i++) {
    switch(op[i].op) {
    case 1:
      printf("  r[%2u] = r[%2u] ^ ", i+1, op[i].a);
      if (op[i].s != 0)
	printf("(r[%2u] << %2u);\n", op[i].b, op[i].s);
      else
	printf("r[%2u];\n", op[i].b);
      break;

    case 0:
      if (op[i].s == 0)
	printf("\n  return r[%u];\n}\n", i);
      else
	printf("\n  return r[%u] << %u;\n}\n", i, op[i].s);
      return;
 
    default:
      printf("error: add stuff here\n"); return;
      // error to reach here
    }
  }

  // error to reach here
}

void clmulk_vm_print_c(clmulk_vm_t* vm)
{
  clmulk_print_c(vm->k, vm->op);
}

//*************************************************************

// dumb but don't care
void printb(uint64_t v, uint32_t b)
{
  uint64_t m = UINT64_C(1)<<(b-1);
  do { printf("%c", (v & m)!=0 ? '1':'.'); m >>= 1; } while (m != 0);
}

void lf(void) { printf("\n"); }

void dump_line_64(char* prefix, uint64_t x)
{
  printf("%s 0x%016lx : ", prefix, x);
  printb(x,64);
  printf(" %2u ", pop_64(x));
}

//*************************************************************

uint64_t clmulk_eval(clmulk_op_t* op, uint64_t x)
{
  uint64_t reg[CLMULK_VM_MAXOPS];

  reg[0] = x;

  for(uint32_t i=0; i<CLMULK_VM_MAXOPS; i++) {
    switch(op[i].op) {
      case 1:
	reg[i+1] = reg[op[i].a] ^ (reg[op[i].b] << op[i].s);
	continue;
      case 0:
	return reg[i] << op[i].s;
      default:
	goto error; 
    }
  }

 error:
  printf("error: add stuff here\n");
  return 0;
}

// kind of pointless wrappers
uint64_t clmulk_vm_eval(clmulk_vm_t* vm, uint64_t x)
{
  return clmulk_eval(vm->op, x);
}

uint64_t clmulk_frag_eval(clmulk_frag_t* f, uint64_t x)
{
  return clmulk_eval(f->op, x);
}


// validate the produced sequence is correct.
// (it's an internal bug to fail)

uint64_t clmulk_sanity(clmulk_op_t* op, const uint64_t k)
{ 
  uint64_t r = 0;

  // If defined perform a full matrix
  // compare. Not needed unless the code
  // starts producing anything other
  // than left xorshifts.
#if defined(CLMULK_VM_INTERNAL_SANITY)
  uint64_t x = 1;
  uint64_t t = k;

  do {
    uint64_t e = clmulk_eval(op, x);
    r  |= t ^ e;
    x <<= 1;
    t <<= 1;
  } while(x);

  r ^= k;

#else
  r = clmulk_eval(op, 1);
#endif  
  
  if (r != k) {
    printf("error: clmulk_sanity failed.\n");
    printf("// produced this incorrect addition chain\n");
    dump_line_64("// desired  k ", k); lf();
    dump_line_64("// produced k ", r); lf();
    clmulk_print_c(k,op);
  }
  
  return r^k;
}


uint64_t clmulk_vm_sanity(clmulk_vm_t* vm)
{
  return clmulk_sanity(vm->op, vm->k);
}


//*************************************************************

static inline void clmulk_mul_assert(uint64_t k, uint32_t pop)
{
  if (k & 1) {
    uint32_t p = pop_64(k);
    if (p == pop)
      return;
    else {
      fprintf(stderr, "error: k=0x%lx : pop=%u when required=%u\n", k, p,pop);
    }
  }
  else
    fprintf(stderr, "error: k=0x%lx : must be odd\n", k);
}


static inline void clmulk_emit_return(clmulk_frag_t* f, uint32_t s)
{
  clmulk_op_t* op = f->op;
  uint32_t     rn = f->rn;

  op[rn] = CLMULK_SHIFT_RET(s);

  //f->cost += (s != 0);
}


// multiply rx by (1 + 2^a)
static inline void clmulk_emit_mul2_i(clmulk_frag_t* f, uint32_t rx, uint32_t a)
{
  clmulk_op_t* op = f->op;
  uint32_t     rn = f->rn;

  op[rn]   = CLMULK_MUL(rx, rx, a); rn++;
  //f->cost += 1;
  
  f->rn = (uint8_t)rn;
}

// multiply rn by (1 + 2^a)
static inline void clmulk_rn_mul2_i(clmulk_frag_t* f, uint32_t a)
{
  clmulk_emit_mul2_i(f, f->rn, a);
}


// multiply rx by odd k with pop(k)=2
static inline void clmulk_emit_mul2(clmulk_frag_t* f, uint32_t rx, uint64_t k)
{
  clmulk_mul_assert(k,2);

  uint64_t n = k;         n = n ^ 1;
  uint32_t a = ctz_64(n); n = n & (n-1);

  clmulk_emit_mul2_i(f,rx,a);
}


// multiply rx by (1 + 2^a + 2^b)
void clmulk_emit_mul3_i(clmulk_frag_t* f, uint32_t rx, uint32_t a, uint32_t b)
{
  clmulk_op_t* op = f->op;
  uint32_t     rn = f->rn;

  op[rn]   = CLMULK_MUL(rx, rx, a); rn++;
  op[rn]   = CLMULK_MUL(rn, rx, b); rn++;
  //f->cost += 2;
  
  f->rn = (uint8_t)rn;
}

// multiply rx by odd k with pop(k)=3
static inline void clmulk_emit_mul3(clmulk_frag_t* f, uint32_t rx, uint64_t k)
{
  clmulk_mul_assert(k,3);

  uint64_t n = k;         n = n ^ 1;
  uint32_t a = ctz_64(n); n = n & (n-1);
  uint32_t b = ctz_64(n); n = n & (n-1);

  clmulk_emit_mul3_i(f,rx,a,b);
}



// NOPE
void clmulk_emit_mul4_i(clmulk_frag_t* f, uint32_t rx, uint32_t a, uint32_t b, uint32_t c)
{
  clmulk_op_t* op = f->op;
  uint32_t     rn = f->rn;

  // can we factor into (1+2^a)(1+2^b)?
  if (c != a+b) {
    op[rn]   = CLMULK_MUL(rn,rx,a); rn++;
    op[rn]   = CLMULK_MUL(rn,rx,b); rn++;
    op[rn]   = CLMULK_MUL(rn,rx,c); rn++;
    //f->cost += 3;
  }
  else {
    op[rn]   = CLMULK_MUL(rn,rx,a); rn++;
    op[rn]   = CLMULK_MUL(rn,rn,b); rn++;
    //f->cost += 2;
  }
  
  f->rn = (uint8_t)rn;
}

// helper function: directly generates 'k' which
// is odd and has exactly 4 bits set into either
// a 2 or 3 op sequence.
void clmulk_emit_mul4(clmulk_frag_t* f, uint32_t rx, uint64_t k)
{
  uint64_t n;
  uint32_t a,b,c;

  clmulk_mul_assert(k,4);

  // k = 1 + 2^a + 2^b + 2^c (a <= b <= c)
  n = k;         n = n ^ 1;
  a = ctz_64(n); n = n & (n-1);
  b = ctz_64(n); n = n & (n-1);
  c = ctz_64(n); n = n & (n-1);

  clmulk_emit_mul4_i(f,rx,a,b,c);  
}

static inline void clmulk_rn_mul2(clmulk_frag_t* f, uint64_t k)
{
  clmulk_emit_mul2(f,f->rn,k);
}

static inline void clmulk_rn_mul3(clmulk_frag_t* f, uint64_t k)
{
  clmulk_emit_mul3(f,f->rn,k);
}

static inline void clmulk_rn_mul4(clmulk_frag_t* f, uint64_t k)
{
  clmulk_emit_mul4(f,f->rn,k);
}

//************************************************

// build schoolbook: pop(k)-1 steps
// TODO: add building ~k + ~0 for large pop
void clmulk_vm_schoolbook(clmulk_vm_t* vm, uint64_t k)
{
  uint32_t i = 0;
  uint32_t s = 0;

  vm->k = k;

  if (k != 0) {
    uint32_t bit;

    
    // make odd and kill the low one.
    // the first op handles the two lowest bits
    s   = ctz_64(k);
    k >>= s;
    k   = k & (k-1);
    
    // for the rest add one at a time
    while(k) {
      bit       = ctz_64(k);
      k         = k & (k-1);
      vm->op[i] = CLMULK_MUL(i,0,bit);
      i++;
    };
  }
  else {
    // zero is just: return x^x
    vm->op[0] = CLMULK_MUL(0,0,0);
    i++;
  }
  
  // reintroduce any initial shift and return
  vm->op[i] = CLMULK_SHIFT_RET(s);
}




pair_u64_t clmul_find_2(uint64_t k)
{
  pair_u64_t divrem;
  
  uint32_t bit = (uint32_t)log2_ceil_u64(k);

  while(bit != 0) {
    uint64_t d = (UINT64_C(1)<<bit) ^ 1;

    divrem = cl_divrem_64(k,d);

    if (divrem.r == 0) {
      return divrem;
    }
    
    bit--;
  }

  return cl_divrem_64(k,3);
}

#if 0
uint64_t clmul_find_3(uint64_t k)
{
}
#endif


//************************************************
// greedy 2-bit/remainder-1 method


static const uint32_t clmulk_tail_max = 4;

#if 0
static inline void
clmulk_small_tail(clmulk_frag_t* f, uint64_t k, uint32_t pop, uint32_t rn)
{
  switch(pop) {
    case 4: rn = clmulk_rn_mul4(f, k); break;
    case 3: rn = clmulk_rn_mul3(f, k); break;
    case 2: rn = clmulk_rn_mul2(f, k); break;
  }
  return rn;
}
#endif


//************************************************
// greedy 2-bit/remainder-1 method

uint32_t clmulk_g2br1_i(clmulk_frag_t* f, uint64_t k, uint32_t rn)
{
  // temp hack
  if (k > 1) {
    uint32_t pop = pop_64(k);
    uint32_t s   = 0;
    
    if (pop & 1) {
      k  ^= 1;
      s   = ctz_64(k);
      k >>= s;
    }
    
    // here: k is even so divisible by some 2-bit number
    // this is "naive" version. change to determine small
    // divisor and find any multiple
    pair_u64_t divrem;
    uint32_t   pos = (uint32_t)log2_ceil_u64(k);
    
    // make a find 2-bit function?
    do {
      divrem = cl_divrem_64(k, (1ul<<pos)^1);
      if (divrem.r == 0) break;
      pos--;
    } while(pos != 0);
    
    if (s == 0) {
      rn        = clmulk_g2br1_i(f, divrem.q, rn);
      f->op[rn] = CLMULK_MUL(rn, rn, pos); rn++;
    }
    else {
      rn        = clmulk_g2br1_i(f, divrem.q, rn);
      f->op[rn] = CLMULK_MUL(rn, rn, pos); rn++;
      f->op[rn] = CLMULK_MUL(0,  rn, s);   rn++;
    }
  }

  return rn;
}


void clmulk_g2br1_new(clmulk_frag_t* f, uint64_t k)
{
  uint32_t s  = ctz_64(k);

  //f->cost += (s != 0);
  
  uint32_t rn = clmulk_g2br1_i(f,k>>s,0);

  f->op[rn] = CLMULK_SHIFT_RET(s);
}



//************************************************
//************************************************

