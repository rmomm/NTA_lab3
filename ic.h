#pragma once
#include <cstdint>
#include <vector>

using namespace std;

uint64_t pow_mod(uint64_t a, uint64_t b, uint64_t mod);
uint64_t gcd(uint64_t a, uint64_t b);
int64_t ext_gcd(int64_t a, int64_t b, int64_t& x, int64_t& y);
uint64_t inv_mod(uint64_t a, uint64_t mod);
vector<uint64_t> buildFactorBase(uint64_t bound);
bool factorOverBase(uint64_t v, const vector<uint64_t>& base, vector<int>& exp);


