
#include <stdint.h>
#include <stdio.h>

#include "clmulk.h"

// hack
extern void clmulk_print_c(uint64_t k, clmulk_op_t* op);

int main(int argc, char** argv)
{
  clmulk_vm_t vm;
  clmulk_frag_t frag = {.op = (clmulk_op_t*)&vm.op, .rn=0};
  
  for(int i=1; i<argc; i++) {
    uint64_t k = strtoul(argv[i], (char**)NULL, 0);

    if (k != 0) {
      vm.k = k;
      clmulk_g2br1(&frag, k);

      // check the sequence is correct
      if (clmulk_sanity(frag.op, k) == 0) {
	clmulk_print_c(k, frag.op);
	printf("\n");
      }
    }
    else
      printf(" skipping %s\n", argv[i]);
  }
  
  return 0;
}
