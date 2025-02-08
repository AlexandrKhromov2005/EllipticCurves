#include <iostream>
#include "curve.hpp"
#include <chrono>

using namespace std::chrono;

int main() {
    mpz_class a("-2"), b("1"), p("41");

    Curve curve(a, b, p);

    std::cout << "Finding points on the curve: y^2 = x^3 + " 
              << a << "x + " << b << " (mod " << p << ")" << std::endl;

    auto start = steady_clock::now();

    curve.find_points();
    curve.calculate_orders();

    auto end = steady_clock::now();

    std::cout << "\nFound points and their orders:" << std::endl;
    curve.print_points();

    auto elapsed = duration_cast<milliseconds>(end - start).count();
    std::cout << "Elapsed time: " << elapsed << " milliseconds." << std::endl;

    return 0;
}
