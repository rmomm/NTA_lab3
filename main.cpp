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


    auto start = high_resolution_clock::now();
    uint64_t x;
    while (true) {
        x = index_calculus_general(alpha, beta, p);

        if (pow_mod(alpha, x, p) == beta)
            break;
    }

    auto end = high_resolution_clock::now();

    double duration = duration_cast<microseconds>(end - start).count();

    cout << "\nCorrect x = " << x << endl;

    uint64_t check = pow_mod(alpha, x, p);

    cout << "Check: " << alpha << "^" << x<< " mod "<< p << " = " << check << endl;

    if (check == beta)
        cout << "Correct solution!\n";
    else
        cout << "Wrong result\n";

    cout << fixed << setprecision(4);

    cout << "\nExecution time: " << duration << " microseconds\n";

    return 0;
}