
#include <stdint.h>
#include <stdio.h>

#include "clmulk.h"

#include "primitives_u16.h"

#define LENGTHOF(X) (sizeof(X)/sizeof(X[0]))

// scaled 1/golden-ratio rounded to odd
const uint64_t weyl_add_u64 = UINT64_C(0x9e3779b97f4a7c15);

// keep the low bit the same. to insure all odd or all even
// walks.
static inline uint64_t weyl_next(uint64_t k)
{
  return k + (weyl_add_u64 ^ 1);
}

typedef uint32_t (frag_func_t)(clmulk_frag_t*, uint64_t);

// dev helper. compare sequence lens of two greedy
// funcs.

#define TRACK_MAX

// temp hack
extern uint32_t clmulk_trace;
			      
// for checking two frag methods:
// 1) produce identical results
// 2) trying to locate "failure" terms for non-indentical
//    (I should make a test suite out of hard terms)			      
void weyl_compare(frag_func_t* f1, frag_func_t* f2)
{
  static const uint32_t len = 1<<20;

  uint32_t cnts[3] = {0};
  
  uint64_t      k = 1;

  uint64_t f1_k    = 0;
  uint32_t f1_diff = 0;
  uint32_t f1_f1_l = 0;
  uint32_t f1_f2_l = 0;
  uint32_t f1_sum  = 0;
  
  uint64_t f2_k    = 0;
  uint32_t f2_diff = 0;
  uint32_t f2_f1_l = 0;
  uint32_t f2_f2_l = 0;
  uint32_t f2_sum  = 0;


  clmulk_vm_t   vm1;
  clmulk_vm_t   vm2;
  clmulk_frag_t frag1 = {.op = (clmulk_op_t*)&vm1.op};
  clmulk_frag_t frag2 = {.op = (clmulk_op_t*)&vm2.op};
  
  for(uint32_t i=0; i<len; i++) {

#if !defined(TRACK_MAX)    
    if ((i & 0xff)==0) { printf("."); fflush(stdout); }
#endif    
    
    frag1.rn = 0;
    frag2.rn = 0;
    
    k = weyl_next(k);

    uint32_t r1 = f1(&frag1, k);
    uint32_t r2 = f2(&frag2, k);

    if (r1 != r2) {
      uint32_t pop = clmulk_pop(k);
      
      if (r1 > r2) {
	// f1 is longer
	uint32_t diff = r1-r2;
	cnts[1]++;
	f1_sum += diff;
	
	if (diff > f1_diff) {
	  f1_diff = diff; f1_k = k; f1_f1_l = r1; f1_f2_l = r2;
#if defined(TRACK_MAX)
	  printf("f1 max: 0x%016lx : r1 = %2u : r2 = %2u : diff %2u : ave diff %f : pop %u\n",
		 f1_k, f1_f1_l, f1_f2_l, f1_diff, (double)f1_sum / (double)cnts[1], pop);
#endif	  
	  
	}
      }
      else {
	// f2 is longer
	uint32_t diff = r2-r1;
	cnts[0]++;
	f2_sum += diff;
	
	if (diff > f2_diff) {
	  f2_diff = diff; f2_k = k; f2_f1_l = r1; f2_f2_l = r2;
#if defined(TRACK_MAX)
	  printf("f2 max: 0x%016lx : r1 = %2u : r2 = %2u : diff %2u : ave diff %f : pop %u\n",
		 f2_k, f2_f1_l, f2_f2_l, f2_diff, (double)f2_sum / (double)cnts[0], pop);
#endif	  
	}
      }
    }
    else cnts[2]++;
  }

  printf("\n");
  printf("f1 max: 0x%016lx : r1 = %2u : r2 = %2u : diff %u : ave diff %f\n",
	 f1_k, f1_f1_l, f1_f2_l, f1_diff, (double)f1_sum / (double)cnts[1]);
  printf("f2 max: 0x%016lx : r1 = %2u : r2 = %2u : diff %u : ave diff %f\n",
	 f2_k, f2_f1_l, f2_f2_l, f2_diff, (double)f2_sum / (double)cnts[0]);


  printf("r1<r2 %u : r2<r1 %u : r1==r2 %u\n", cnts[0],cnts[1],cnts[2]);
}

			      
uint32_t weyl_frag(frag_func_t* f)
{
  static const uint32_t len = 0xffff;
  
  double   t = 0;
  uint64_t k = 1;
  uint32_t e = 0;

  clmulk_vm_t   vm;
  clmulk_frag_t frag = {.op = (clmulk_op_t*)&vm.op};
  
  for(uint32_t i=0; i<len; i++) {
    frag.rn = 0;
    
    k = weyl_next(k);

    if ((i & 0xff)==0) { printf("."); fflush(stdout); }
    
    uint32_t l = f(&frag, k)-1;
    uint32_t p = clmulk_pop(k);
    double   v = ((double)p-(double)l)/(double)p;

#if 0
    if (l < (p >> 1)) {
      clmulk_print_c(k, frag.op);
    }
#endif    
    
    t += v;

    if (clmulk_sanity(frag.op, k) == 0)
      continue;

    if (e++ > 300) {
      printf("**********bailing : i=%u e=%u\n",i,e); return 0;
    }
  }

  printf("done: %f\n", t*(1.0/(double)len));
  return 1;
}

uint32_t simple_frag(frag_func_t* f)
{
  uint32_t e = 0;
  double   t = 0;

  printf("foo\n");
  
  clmulk_vm_t   vm;
  clmulk_frag_t frag = {.op = (clmulk_op_t*)&vm.op};
  
  for(uint32_t i=0; i<LENGTHOF(gf2_primitives_u16); i++) {
    uint64_t k = gf2_primitives_u16[i];
    
    frag.rn = 0;

    uint32_t l = f(&frag, k)-1;
    double   p = (double)clmulk_pop(k);
    double   v = (p-(double)l)/p;

    t += v;

    if (clmulk_sanity(frag.op, k) == 0)
      continue;

    if (e++ > 300) {
      printf("**********bailing : i=%u e=%u\n",i,e); return 0;
    }
  }
  printf("done: %f\n", t*(1.0/(double)LENGTHOF(gf2_primitives_u16)));
  return 1;
}


extern uint32_t clmulk_greedy(clmulk_frag_t* f, uint64_t k);

int main(void)
{
#if 0  
  if (!simple_frag(&clmulk_g2br1))  return -1;
  if (!simple_frag(&clmulk_greedy)) return -1;
  if (!weyl_frag(&clmulk_g2br1))    return -1;
#endif  
  if (!weyl_frag(&clmulk_greedy))   return -1;

  weyl_compare(&clmulk_g2br1, &clmulk_greedy);
  
  return 0;
}
