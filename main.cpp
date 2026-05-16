#include <iostream>
#include <chrono>
#include <iomanip>

#include "ic.h"

using namespace std;
using namespace std::chrono;

int main()
{
    uint64_t alpha, beta, p;

    cout << "Enter alpha: ";
    cin >> alpha;

    cout << "Enter beta: ";
    cin >> beta;

    cout << "Enter prime p: ";
    cin >> p;

    cout << fixed << setprecision(4);

    auto start1 = high_resolution_clock::now();

    uint64_t x1;

    while (true) {
        x1 = index_calculus_general(alpha, beta, p);

        if (pow_mod(alpha, x1, p) == beta)
            break;
    }

    auto end1 = high_resolution_clock::now();

    double duration1 = duration_cast<microseconds>(end1 - start1).count();

    cout << "\n index-calculus \n";
    cout << " x = " << x1 << endl;
    uint64_t check1 = pow_mod(alpha, x1, p);

    cout << "Check: " << alpha << "^" << x1   << " mod " << p   << " = " << check1 << endl;

    if (check1 == beta)
        cout << "Correct solution!\n";
    else
        cout << "Wrong result\n";

    cout << "Execution time: "  << duration1 << " microseconds\n";


    auto start2 = high_resolution_clock::now();

    uint64_t x2;

    while (true) {
        x2 = index_calculus_general_parallel(alpha, beta, p);

        if (pow_mod(alpha, x2, p) == beta)
            break;
    }

    auto end2 = high_resolution_clock::now();

    double duration2 = duration_cast<microseconds>(end2 - start2).count();

    cout << "\n index-calculus (parallel)\n";

    cout << " x = " << x2 << endl;

    uint64_t check2 = pow_mod(alpha, x2, p);

    cout << "Check: " << alpha << "^" << x2 << " mod " << p << " = " << check2 << endl;

    if (check2 == beta)
        cout << "Correct solution!\n";
    else
        cout << "Wrong result\n";

    cout << "Execution time: " << duration2 << " microseconds\n";

    return 0;
}