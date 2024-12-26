#pragma once
// Randomize the seed for generating random numbers
void init_random_seed() {
    dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
}