
#include <stdint.h>
#include <stdio.h>

#include "clmulk.h"

#include "primitives_u16.h"

// temp hack
extern void clmulk_g2br1_new(clmulk_frag_t* f, uint64_t k);


#define LENGTHOF(X) (sizeof(X)/sizeof(X[0]))

const uint64_t weyl_add_u64 = UINT64_C(0x9e3779b97f4a7c15);

typedef void (frag_func_t)(clmulk_frag_t*, uint64_t);

void weyl_frag(frag_func_t* f)
{
  uint64_t k = 1;
  uint32_t e = 0;

  clmulk_vm_t   vm;
  clmulk_frag_t frag = {.op = (clmulk_op_t*)&vm.op};
  
  for(uint32_t i=0; i<0xfffff; i++) {
    frag.rn = 0;
    
    k += weyl_add_u64;

    uint64_t k0 = k ;//& 0xffff;

    f(&frag, k0);

    if (clmulk_sanity(frag.op, k0) == 0)
      continue;

    if (e++ > 300) {
      printf("**********bailing : i=%u e=%u\n",i,e); return;
    }
  }
}

void simple_frag(frag_func_t* f)
{
  uint32_t e = 0;

  clmulk_vm_t   vm;
  clmulk_frag_t frag = {.op = (clmulk_op_t*)&vm.op};
  
  for(uint32_t i=0; i<LENGTHOF(gf2_primitives_u16); i++) {
    uint64_t k = gf2_primitives_u16[i];
    
    frag.rn = 0;

    f(&frag, k);

    if (clmulk_sanity(frag.op, k) == 0)
      continue;

    if (e++ > 300) {
      printf("**********bailing : i=%u e=%u\n",i,e); return;
    }
  }
}


int main(void)
{
  weyl_frag(&clmulk_g2br1_new);
  simple_frag(&clmulk_g2br1_new);
  return 0;
}
