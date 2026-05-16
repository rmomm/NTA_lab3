#pragma once
#include <iostream>
#include <vector>
#include <cstdint>
#include <chrono>
#include <random>
#include <algorithm>

using namespace std;

uint64_t pow_mod(uint64_t a, uint64_t b, uint64_t mod);
uint64_t gcd(uint64_t a, uint64_t b);
int64_t ext_gcd(int64_t a, int64_t b, int64_t& x, int64_t& y);
uint64_t inv_mod(uint64_t a, uint64_t mod);
uint64_t computeBound(uint64_t n);
vector<uint64_t> buildFactorBase(uint64_t B);
bool factorOverBase(uint64_t v, const vector<uint64_t>& base, vector<int>& exp);

struct Relation {
    uint64_t k;
    vector<int> exponents;
}; 

bool tryRelation(uint64_t alpha, uint64_t k, uint64_t p, const vector<uint64_t>& base, Relation& r);
vector<Relation> collectRelations(uint64_t alpha, uint64_t p, const vector<uint64_t>& base, size_t n);

vector<uint64_t> factor(uint64_t n);
bool is_generator(uint64_t a, uint64_t p);
uint64_t findGenerator(uint64_t p);
uint64_t solveCongruence(uint64_t u, uint64_t v, uint64_t mod);


vector<uint64_t> SLSMod(vector<vector<int>> A, vector<uint64_t> b, uint64_t mod);
uint64_t computeLogBeta(uint64_t alpha, uint64_t beta, uint64_t p, const vector<uint64_t>& base, const vector<uint64_t>& log_p);

uint64_t index_calculus(uint64_t alpha, uint64_t beta, uint64_t p);
uint64_t index_calculus_general(uint64_t alpha, uint64_t beta, uint64_t p);