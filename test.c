
#include <stdint.h>
#include <stdio.h>

#include "clmulk.h"

#include "primitives_u16.h"

#define LENGTHOF(X) (sizeof(X)/sizeof(X[0]))

const uint64_t weyl_add_u64 = UINT64_C(0x9e3779b97f4a7c15);

typedef uint32_t (frag_func_t)(clmulk_frag_t*, uint64_t);

void weyl_frag(frag_func_t* f)
{
  static const uint32_t len = 0xffff;
  
  double   t = 0;
  uint64_t k = 1;
  uint32_t e = 0;

  clmulk_vm_t   vm;
  clmulk_frag_t frag = {.op = (clmulk_op_t*)&vm.op};
  
  for(uint32_t i=0; i<len; i++) {
    frag.rn = 0;
    
    k += weyl_add_u64;

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
      printf("**********bailing : i=%u e=%u\n",i,e); return;
    }
  }

  printf("done: %f\n", t*(1.0/(double)len));
}

void simple_frag(frag_func_t* f)
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
      printf("**********bailing : i=%u e=%u\n",i,e); return;
    }
  }
  printf("done: %f\n", t*(1.0/(double)LENGTHOF(gf2_primitives_u16)));
}


int main(void)
{
  simple_frag(&clmulk_g2br1);
  weyl_frag(&clmulk_g2br1);
  return 0;
}
