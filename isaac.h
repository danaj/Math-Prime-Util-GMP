#ifndef MPU_ISAAC_H
#define MPU_ISAAC_H

#include "ptypes.h"

extern void isaac_init(uint32_t bytes, const unsigned char* data);
extern uint32_t isaac_rand32(void);
extern uint32_t isaac_rand(uint32_t n);
extern void isaac_rand_bytes(uint32_t bytes, unsigned char* data);
extern int isaac_seeded(void);

#endif
