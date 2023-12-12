
#ifndef CLMULK_H
#define CLMULK_H

//************************************************
// tiny register machine 
//
// * only two operations:
//   1) r_{n+1} = r_a ^ (r_b << s)
//   2) return r_n << s
// * destintation register is implict. each is
//   is only assigned to once.

typedef struct {
  uint8_t op;
  uint8_t a;
  uint8_t b;
  uint8_t s;
} clmulk_op_t;

// enough space for schoolbook (even that can do
// better than this if encode via (~k + ~0)
// support is added for large popcount)
#define CLMULK_VM_MAXOPS 65

// quick wrapper
typedef struct {
  uint64_t    k;
  clmulk_op_t op[CLMULK_VM_MAXOPS];
} clmulk_vm_t;


#define CLMULK_DEF(OP,A,B,S) (clmulk_op_t){.op=(uint8_t)(OP),.a=(uint8_t)(A),.b=(uint8_t)(B),.s=(uint8_t)(S)}

#define CLMULK_INVALID(A,B,S)  CLMULK_DEF(0xff,A,B,S)
#define CLMULK_MUL(A,B,S)      CLMULK_DEF(1,A,B,S)
#define CLMULK_ADD(A,B)        CLMULK_DEF(1,A,B,0)
#define CLMULK_RET             CLMULK_DEF(0,0,0,0)
#define CLMULK_SHIFT_RET(S)    CLMULK_DEF(0,0,0,S)

// 
typedef struct {
  uint8_t rx;  // kill? only in hacky test ATM. belongs on logical stack anyway
//uint8_t r0;  // kill? ditto

  uint8_t rn;
//uint8_t shift;
//uint8_t cost;
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
  clmulk_op_t   op1[CLMULK_VM_MAXOPS];
  clmulk_op_t   op2[CLMULK_VM_MAXOPS];
} clmulk_builder_t;


extern void     clmulk_g2br1(clmulk_frag_t* f, uint64_t k);

extern uint64_t clmulk_sanity(clmulk_op_t* op, const uint64_t k);

#endif
