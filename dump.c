
#include <stdint.h>
#include <stdio.h>

#include "clmulk.h"

// hacks
extern uint32_t clmulk_greedy(clmulk_frag_t* f, uint64_t k);
extern void clmulk_print_c(uint64_t k, clmulk_op_t* op);

  
extern uint32_t clmulk_trace;

int main(int argc, char** argv)
{
  clmulk_vm_t vm;
  clmulk_frag_t frag = {.op = (clmulk_op_t*)&vm.op, .rn=0};

  clmulk_trace = 1;
  
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

      clmulk_greedy(&frag, k);
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
