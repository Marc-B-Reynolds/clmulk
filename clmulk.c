// Marc B. Reynolds, 2023
// Public Domain under http://unlicense.org, see link for details.


#include <stdint.h>

#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <stdbool.h>
#include <assert.h>

#if 1
#if defined(__GNUC__) && !defined(__clang__)
static_assert(0, "problem under GCC atm. infinite loop. opps!");
#endif
#endif

/************************************************

  All of this is a major mess & eyesore. I'm hacking
  in small bursts. XXX

  -----------------------------------------------

  s-significant bit number:
  * even parity has (at least one) factor of 3
    more generally there's some 2-bit number that is
    a divisor.
  * 2-bit  irreducible numbers have form 1 {0}_n 1
    where 'n' is an even number (including zero).
    first few are 0x3,0x9,0x21,0x81..

 ************************************************/

// * recursive structure is intended as a dev hack
//   can be converted to non-recursive.
// * initial cost limit should be pop(k)-1. can be
//   tighter than that in high pop vs. num of sig bits
//   cases (low transition count).
// * check on: ceil(log2(transition_count)) as low bound
//   estimate. this is ad-hoc based on my current thinking
//   of the "space" that's explored.
// * put asserts back in and add tracing
// * memorization:
//   * internal: look-up sub-expression during building.
//     Don't think it's useful outside of empirical
//     testing.
//   * root: keeping a small LRU could be useful for
//     a compiler.
//   * humm...but a searching version will hit the same
//     terms frequently?? 
// * cost modeling. currently there is none as I'm only
//   producing one dep-chain. So it's almost the number
//   of vm-ops produced (the exception being if the
//   return op is shifted or not). This would need to
//   change if I think of a useful additive factoring.

//*************************************************************
// method of attack thinking
//
// a possible method (not doing this yet)
// 1. reduce to odd (low clear)
// 2. while even parity
//    * find largest 2-bit divisor. can be done by looking at
//      irreduciables and if found determine the largest power
//      of two of that which is a divisor.
//    * quotient is alway odd
// 3. while odd parity
//    * see if any small 3-bit primitives are divisors. if not
//      break to 4
//    * find the largest power-of-two of found and factor
//    * quotient is alway odd
//    * SIGH! sometime much better and sometimes much worse
//      with this step ATM. need to figure it out for greedy.
// 4. punt and factor k into (k+1)+1. sometimes sucks.
//    sometimes works really well. sigh.
//    need to consider other "punt" methods

//*************************************************************
// big problems
// * greed methods aren't currently limiting an can produce
//   longer chain than schoolbook! Doh!



/*
 first 8 irreduciables of odd n-bits.
  
 3-bit :
....111   67109828  25.0004%
...1.11   33558864  12.5017
...11.1   33552933  12.4994
..1..11   16777772   6.2502
..11..1   16777312   6.2500
.1..1.1    8389934   3.1255
.1.1..1    8394884   3.1273
1....11    4195571   1.5630
         125147420  46.6210 (none)

5-bit : 
..11111 16785389  6.2530%
.1.1111  8394103  3.1270
.11.111  8389880  3.1255
.111.11  8390993  3.1259
.1111.1  8391199  3.1260
1.1.111  4192705  1.5619
1.11.11  4190816  1.5612
11..111  4193318  1.5621
                 78.7543 (none)

7-bit (waste of time)
..1.111111  2098432  0.7817%
..111.1111  2099451  0.7821
..1111.111  2097352  0.7813
..111111.1  2098649  0.7818
.1..111111  1049297  0.3909
.1.1.11111  1048021  0.3904
.1.111.111  1048324  0.3905
.1.1111.11  1048535  0.3906
                    95.4044 (none)		
*/


//*************************************************************


#include "clmulk.h"


//*************************************************************

static inline uint64_t clmulk_pop_next(uint64_t x)
{
  uint64_t t = x + (x & -x);
  x = x & ~t;
  x = (uint64_t)((int64_t)x >> clmulk_ctz(x));
  x = (uint64_t)((int64_t)x >> 1);

  return t^x;
}

static inline uint64_t clmulk_pop_prev(uint64_t x)
{
  uint64_t a = ~x & (x+1); 
  uint64_t b =  x - a;
  uint64_t c = ~x & b;
  uint64_t d = (uint64_t)((int64_t)c >> clmulk_ctz(c));
  d = (uint64_t)((int64_t)c >> 1);

  return b^d;
}

static inline uint64_t clmulk_pop_next_odd(uint64_t x)
{
  x = clmulk_pop_next((uint64_t)((int64_t)x >> 1));
  x = (x<<1) ^ 1;
  return x;
}

static inline uint64_t clmulk_pop_prev_odd(uint64_t x)
{
  x = clmulk_pop_prev((uint64_t)((int64_t)x >> 1));
  x = (x<<1) ^ 1;
  return x;
}


clmulk_pair_t clmulk_divrem(uint64_t a, uint64_t b)
{
  uint64_t q  = 0;

  // nuke compare. illegal. add assert
  
  if (b != 0) {
    uint64_t lb = clmulk_clz(b);
    uint64_t t  = lb - clmulk_clz(a);
    
    while ((int64_t)t >= 0) {
      a ^= b << t;
      q ^= UINT64_C(1) << t;
      t  = lb - clmulk_clz(a);
    }
  }
    
  return (clmulk_pair_t){.q=q, .r=a};
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
  printf(" %2u ", clmulk_pop(x));
}


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

  // dumping line doesn't need a comment (needed?)
  uint32_t skip = 0;
  
  dump_line_64("//", k);
  printf("\nuint64_t cl_mulk_%lx(uint64_t x)\n{\n  uint64_t r[%u];\n\n  r[ 0] = x;\n", k,len+1);

  for(uint32_t i=0; i<CLMULK_MAXOPS; i++) {
    switch(op[i].op) {
    case 1:
      printf("  r[%2u] = r[%2u] ^ ", i+1, op[i].a);
      if (op[i].s != 0)
	printf("(r[%2u] << %2u);", op[i].b, op[i].s);
      else
	printf("r[%2u];", op[i].b);

      // attempt to add comment (silly to deduce...but hacking yo!)
      if (skip == 0) {
	printf("  // ");

	// is r[rn] ^ (rn[rn] << s)
	if (op[i].a == op[i].b) {

	  // is next op not a return
	  if (op[i+1].op == 1) {
	    
	    if (op[i+1].a == 0) {
	      uint64_t n = 1ul ^ (1ul << op[i].s);
	      printf("0x%lx*rn+1", n << op[i+1].s);
	      skip = 1;
	    }
	    else if (op[i+1].b == i) {
	      uint64_t n = 1ul ^ (1ul << op[i].s) ^ (1ul << op[i+1].s);
	      printf("0x%lx*rn (3-bit)",n);
	      skip = 1;
	    }
	    else {
	      uint64_t n = 1ul ^ (1ul << op[i].s);
	      printf("0x%lx*rn (2-bit)",n);
	    }
	  }
	  else {
	    uint64_t n = 1ul ^ (1ul << op[i].s);
	    printf("0x%lx*rn (2-bit)",n);
	  }
	}
      }
      else skip = 0;
      
      printf("\n");
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

// doesn't general work ATM. for internal hacking
#if 1
void clmulk_print_rn_c(uint64_t k, clmulk_op_t* op)
{
  uint32_t len = clmulk_len(op);

  printf("uint64_t cl_mulk_%lx(uint64_t x)\n{\n  uint64_t r[%u];\n  uint32_t rn = 0;\n\n  r[ 0] = x;\n", k,len+1);

  for(uint32_t i=0; i<CLMULK_MAXOPS; i++) {
    switch(op[i].op) {

    case 1:

      if ((op[i].a != 0) && (i != 0)) 
	printf("  r[rn+1] = r[rn] ^ ");
      else
	printf("  r[rn+1] = r[ 0] ^ ");

      if (op[i].s != 0)
	printf("(r[rn] << %2u); rn++;\n", op[i].s);
      else
	printf("r[rn]; rn++;\n");

      break;

    case 0:
      if (op[i].s == 0)
	printf("\n  return r[rn];\n}\n");
      else
	printf("\n  return r[rn] << %u;\n}\n", op[i].s);
      return;
 
    default:
      printf("error: add stuff here\n"); return;
      // error to reach here
    }
  }

  // error to reach here
}
#endif

void clmulk_vm_print_c(clmulk_vm_t* vm)
{
  clmulk_print_c(vm->k, vm->op);
}


//*************************************************************

uint64_t clmulk_eval(clmulk_op_t* op, uint64_t x)
{
  uint64_t reg[CLMULK_MAXOPS];

  reg[0] = x;

  for(uint32_t i=0; i<CLMULK_MAXOPS; i++) {
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
// dev hack
uint32_t clmulk_trace = 0;

// temp name
static inline void trace_str(char* str)
{
  if (clmulk_trace) {
    printf("%s\n", str);
  }
}

static inline void trace_line(char* str, uint64_t n)
{
  if (clmulk_trace) {
    dump_line_64(str,n);
    printf("\n");
  }
}


//*************************************************************

static inline void clmulk_mul_assert(uint64_t k, uint32_t pop)
{
  if (k & 1) {
    uint32_t p = clmulk_pop(k);
    if (p == pop)
      return;
    else {
      fprintf(stderr, "error: k=0x%lx : pop=%u when required=%u\n", k, p,pop);
    }
  }
  else
    fprintf(stderr, "error: k=0x%lx : must be odd\n", k);
}



//************************************************

static inline uint32_t clmulk_limit(uint32_t limit, uint64_t k)
{
  uint32_t t = clmulk_pop(k)-1;
  return limit <= t ? limit : t;
}


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
    s   = clmulk_ctz(k);
    k >>= s;
    k   = k & (k-1);
    
    // for the rest add one at a time
    while(k) {
      bit       = clmulk_ctz(k);
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




clmulk_pair_t clmulk_find_2(uint64_t k)
{
  clmulk_pair_t divrem;
  
  uint32_t bit = clmulk_log2_ceil(k);

  while(bit != 0) {
    uint64_t d = (UINT64_C(1)<<bit) ^ 1;

    divrem = clmulk_divrem(k,d);

    if (divrem.r == 0) {
      return divrem;
    }
    
    bit--;
  }

  return clmulk_divrem(k,3);
}

// expensive. for experimenting
#if 0
clmulk_pair_t clmulk_find_3(uint64_t k)
{
  clmulk_pair_t divrem;

  uint64_t d = 0x7;
  uint64_t e = k >> 2;

  while(d <= e) {
    divrem = clmulk_divrem(k,d);

    if (divrem.r == 0) {

#if 1
      dump_line_64("k", k);        lf();
      dump_line_64("d", d);        lf();
      dump_line_64("q", divrem.q); lf();
#endif      

      divrem.r = d;

      return divrem;
    }
    
    d = clmulk_pop_next_odd(d);
    
    if (clmulk_pop(d) != 3) {printf("what?\n");}
  }

  divrem.r = 0;

  return divrem;
}
#elif 0
clmulk_pair_t clmulk_find_3(uint64_t k)
{
  clmulk_pair_t divrem;

  uint64_t d = (UINT64_C(3)<<62) >> clmulk_clz(k);
  
  d = (d>>2)|1;

#if 0
  dump_line_64("k  ", k);  lf();
  dump_line_64("div", d);  lf();
#endif  

  while(d != 0) {
    if (divrem.r == 0) {
      
#if 1
      dump_line_64("k", k);        lf();
      dump_line_64("d", d);        lf();
      dump_line_64("q", divrem.q); lf();
#endif      
      
      divrem.r = d;
      
      return divrem;
    }
    
    uint64_t d2 = clmulk_pop_prev_odd(d);

    if (d == d2) { printf("no\n"); }

    d = d2;
    
    if (clmulk_pop(d) != 3) {printf("what?\n");}
  }
  
  divrem.r = 0;
  
  return divrem;
}

#else

// dumb version
clmulk_pair_t clmulk_find_3(uint64_t k)
{
  //uint32_t pop = clmulk_pop(k);

  // greedy popcount filtering hurts sometimes...sigh
  
  // powers of 7
  uint64_t      div    = 0x7;
  clmulk_pair_t divrem = clmulk_divrem(k, div);

  if (divrem.r == 0) {
    clmulk_pair_t dr2 = clmulk_divrem(k,0x15);
    
    if (dr2.r == 0) {
      clmulk_pair_t dr3 = clmulk_divrem(k,0x111);
      
      if (dr3.r == 0) {
	dr3.r = 0x111;
	
	//if (clmulk_pop(dr3.q) < pop)
	  return dr3;
      }
      
      dr2.r = 0x15;
      //if (clmulk_pop(dr2.q) < pop)
	return dr2;
    }
    
    divrem.r = 0x7;
    
    //if (clmulk_pop(dr2.q) < pop)
      return divrem;
  }
  
  divrem.r = 0;
  
  return divrem;
}

#endif
  

//************************************************
// goofying around. these are tots inconsistent.

static inline void clmulk_emit_return(clmulk_frag_t* f, uint32_t s)
{
  clmulk_op_t* op = f->op;
  uint32_t     rn = f->rn;

  op[rn] = CLMULK_SHIFT_RET(s);
}


static inline uint32_t clmulk_mul2(clmulk_frag_t* f, uint32_t rn, uint32_t a)
{
  f->op[rn] = CLMULK_MUL(rn,rn,a); rn++;
  return rn;
}

static inline uint32_t clmulk_add1(clmulk_frag_t* f, uint32_t rn, uint32_t a)
{
  f->op[rn] = CLMULK_MUL(0,rn,a); rn++;
  return rn;
}

static inline uint32_t
clmulk_rx_mul3(clmulk_frag_t* f, uint32_t rx, uint32_t rn, uint32_t a, uint32_t b)
{
  clmulk_op_t* op = f->op;

  op[rn] = CLMULK_MUL(rx, rx, a); rn++;
  op[rn] = CLMULK_MUL(rn, rx, b); rn++;
  
  return rn;
}

static inline uint32_t
clmulk_rx_mul4(clmulk_frag_t* f, uint32_t rn, uint32_t rx, uint32_t a, uint32_t b, uint32_t c)
{
  clmulk_op_t* op = f->op;

  // can we factor into (1+2^a)(1+2^b)?
  if (c != a+b) {
    op[rn] = CLMULK_MUL(rn,rx,a); rn++;
    op[rn] = CLMULK_MUL(rn,rx,b); rn++;
    op[rn] = CLMULK_MUL(rn,rx,c); rn++;
  }
  else {
    op[rn] = CLMULK_MUL(rn,rx,a); rn++;
    op[rn] = CLMULK_MUL(rn,rn,b); rn++;
  }
  
  return rn;
}

static inline uint32_t
clmulk_mul3(clmulk_frag_t* f, uint32_t rn, uint32_t a, uint32_t b)
{
  return clmulk_rx_mul3(f,rn,rn,a,b);
}

static inline uint32_t
clmulk_mul4(clmulk_frag_t* f, uint32_t rn, uint32_t a, uint32_t b, uint32_t c)
{
  return clmulk_rx_mul4(f,rn,rn,a,b,c);
}

static inline uint32_t
clmulk_emit_mul2(clmulk_frag_t* f, uint32_t rn, uint64_t k)
{
  uint64_t n = k;             n = n ^ 1;
  uint32_t a = clmulk_ctz(n); n = n & (n-1);

  return clmulk_mul2(f,rn,a);
}

static inline uint32_t
clmulk_emit_mul3(clmulk_frag_t* f, uint32_t rn, uint64_t k)
{
  uint64_t n = k;             n = n ^ 1;
  uint32_t a = clmulk_ctz(n); n = n & (n-1);
  uint32_t b = clmulk_ctz(n); n = n & (n-1);

  return clmulk_mul3(f,rn,a,b);
}

static inline uint32_t
clmulk_emit_mul4(clmulk_frag_t* f, uint32_t rn, uint64_t k)
{
  uint64_t n;
  uint32_t a,b,c;

  // k = 1 + 2^a + 2^b + 2^c (a <= b <= c)
  n = k;             n = n ^ 1;
  a = clmulk_ctz(n); n = n & (n-1);
  b = clmulk_ctz(n); n = n & (n-1);
  c = clmulk_ctz(n); n = n & (n-1);

  return clmulk_mul4(f,rn,a,b,c);  
}


// short-cut constants with popcount smaller
// than `clmulk_tail_max`.

const uint32_t clmulk_tail_max = 4;

static inline uint32_t
clmulk_tail(clmulk_frag_t* f, uint64_t k, uint32_t rn)
{
  // caller has this. compiler should remove.
  uint32_t pop = clmulk_pop(k);

  switch(pop) {
    case 4: rn = clmulk_emit_mul4(f, rn, k); break;
    case 3: rn = clmulk_emit_mul3(f, rn, k); break;
    case 2: rn = clmulk_emit_mul2(f, rn, k); break;
  }
  return rn;
}


static void foo(void)
{
}


//************************************************
// greedy 2-bit/remainder-1 method


uint32_t clmulk_g2br1_i(clmulk_frag_t* f, uint64_t k, uint32_t rn)
{
  clmulk_pair_t divrem;

  uint32_t pop = clmulk_pop(k);
  
  if (pop > clmulk_tail_max) {
    uint32_t s   = 0;
    
    if (pop & 1) {
      k  ^= 1;
      s   = clmulk_ctz(k);
      k >>= s;
    }
    // here: k is even so divisible by some 2-bit number
    // this is "naive" version. change to determine small
    // divisor and find any multiple
    uint32_t   pos = clmulk_log2_ceil(k);
    
    // make a find 2-bit function?
    // - return q & pos
    // - but a try would go in the center...but that's
    //   a different function.
    do {
      divrem = clmulk_divrem(k, (1ul<<pos)^1);
      if (divrem.r == 0) break;
      pos--;
    } while(pos != 0);

    if (s != 0) {
      trace_str("+1 factor");
      trace_line("  k ", k);
    //trace_line("  k'", divrem.q);
      trace_line("  f ", (1ul<<pos)^1);
    }
    else {
      trace_str("2-bit factor");
      trace_line("  k ", k);
    //trace_line("  k'", divrem.q);
      trace_line("  f ", (1ul<<pos)^1);
    }

    rn = clmulk_g2br1_i(f, divrem.q, rn);
    rn = clmulk_mul2(f, rn, pos);
    
    if (s != 0)
      rn = clmulk_add1(f,rn,s);

    return rn;
  }

  return clmulk_tail(f,k,rn);
}

// temp wrapper
uint32_t clmulk_g2br1(clmulk_frag_t* f, uint64_t k)
{
  uint32_t s  = clmulk_ctz(k);
  uint32_t rn = clmulk_g2br1_i(f,k>>s,0);
  
  f->op[rn] = CLMULK_SHIFT_RET(s);

  return rn;
}


//************************************************

#define recurse clmulk_greedy_i

uint32_t clmulk_greedy_i(clmulk_frag_t* f, uint64_t k, uint32_t rn)
{
  clmulk_pair_t divrem;

  uint32_t pop = clmulk_pop(k);
  
  if (pop > clmulk_tail_max) {
    uint32_t s   = 0;
    
    if (pop & 1) {
#if 1
      divrem = clmulk_find_3(k);

      if (divrem.r) {
	trace_str("3-bit factor");
	trace_line("  k ", k);
      //trace_line("  k'", divrem.q);
	trace_line("  f ", divrem.r);
	rn = recurse(f, divrem.q, rn);
	rn = clmulk_emit_mul3(f,rn,divrem.r);
	return rn;
      }
#endif      
      
      k  ^= 1;
      s   = clmulk_ctz(k);
      k >>= s;
    }
    
    // here: k is even so divisible by some 2-bit number
    // this is "naive" version. change to determine small
    // divisor and find any multiple
    uint32_t   pos = clmulk_log2_ceil(k)-1;

    // make a find 2-bit function?
    // - return q & pos
    // - but a try would go in the center...but that's
    //   a different function.
    do {
      divrem = clmulk_divrem(k, (1ul<<pos)^1);
      if (divrem.r == 0) break;
      pos--;
    } while(pos != 0);

    if (s != 0) {
      trace_str("+1 factor");
      trace_line("  k ", k);
    //trace_line("  k'", divrem.q);
      trace_line("  f ", (1ul<<pos)^1);
    }
    else {
      trace_str("2-bit factor");
      trace_line("  k ", k);
    //trace_line("  k'", divrem.q);
      trace_line("  f ", (1ul<<pos)^1);
    }
    
    rn = recurse(f, divrem.q, rn);
    rn = clmulk_mul2(f, rn, pos);
    
    if (s != 0) {
      rn = clmulk_add1(f,rn,s);
    }
    return rn;
  }

  return clmulk_tail(f,k,rn);
}

#undef recurse



// temp wrapper
uint32_t clmulk_greedy(clmulk_frag_t* f, uint64_t k)
{
  uint32_t s  = clmulk_ctz(k);
  uint32_t rn = clmulk_greedy_i(f,k>>s,0);
  
  f->op[rn] = CLMULK_SHIFT_RET(s);

  return rn;
}


//************************************************

uint32_t clmulk_g23b_i(clmulk_frag_t* f, uint64_t k, uint32_t rn);


uint32_t clmulk_g23b_o(clmulk_frag_t* f, uint64_t k, uint32_t rn)
{
  clmulk_pair_t divrem;

  uint32_t s = 0;
  
  k  ^= 1;
  s   = clmulk_ctz(k);
  k >>= s;
  
  uint32_t pos = clmulk_log2_ceil(k);
  
  do {
    divrem = clmulk_divrem(k, (1ul<<pos)^1);
    if (divrem.r == 0) break;
    pos--;
  } while(pos != 0);
  
  rn = clmulk_g23b_i(f, divrem.q, rn);
  rn = clmulk_mul2(f, rn, pos);
  rn = clmulk_add1(f,rn,s);
  
  return rn;
}

uint32_t clmulk_g23b_e(clmulk_frag_t* f, uint64_t k, uint32_t rn)
{
  uint32_t pos = clmulk_log2_ceil(k);
  
  clmulk_pair_t divrem;
  
  do {
    divrem = clmulk_divrem(k, (1ul<<pos)^1);
    if (divrem.r == 0) break;
    pos--;
  } while(pos != 0);
  
  rn = clmulk_g23b_i(f, divrem.q, rn);
  rn = clmulk_mul2(f, rn, pos);

  return rn;
}

uint32_t clmulk_g23b_i(clmulk_frag_t* f, uint64_t k, uint32_t rn)
{
  if (k > 1) {
    uint32_t pop = clmulk_pop(k);

    if (pop & 1)
      return clmulk_g23b_o(f,k,rn);

    return clmulk_g23b_e(f,k,rn);
  }
  
  return rn;
}

//************************************************

#if 0
static inline uint32_t clmulk_build_try(uint32_t rn)
{
  uint32_t r0 = f(build, rn, k);
  
  
  
  return rn;
}
#endif
