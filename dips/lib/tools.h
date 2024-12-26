#pragma once
#include "SFMT/dSFMT/dSFMT.h"
void init_random_seed() {
    // Randomize the seed for generating random numbers
    dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
}