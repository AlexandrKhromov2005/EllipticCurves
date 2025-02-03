#include <iostream>
#include "curve.hpp"
#include "point.hpp"
#include "algorithms_for_primes.hpp"

int main() {
    // Example 1: Generate and display all points
    {
        std::cout << "=== Example 1: Generating all points ===" << std::endl;
        mpz_class a = 1, b = 1, p = 23; // Curve y² = x³ + x + 1 over GF(23)
        Curve curve(a, b, p);

        curve.find_points(); 

        std::cout << "Curve: y^2 = x^3 + " 
                  << curve.get_a() << "x + " 
                  << curve.get_b() << " over GF(" 
                  << curve.get_p() << ")\n";

        std::cout << "\nTotal points: " << curve.points.size() << "\n";
        curve.print_points();
        std::cout << "\n";
    }

    // Example 2: Validate generated points
    {
        std::cout << "\n=== Example 2: Point validation ===" << std::endl;
        mpz_class a = 0, b = 7, p = 17; // Curve y² = x³ + 7 over GF(17)
        Curve curve(a, b, p);
        curve.find_points();

        int invalid_count = 0;
        for (const Point& p : curve.points) {
            if (!curve.is_on_curve(p)) {
                std::cerr << "ERROR: Point (" 
                          << p.get_x() << ", " 
                          << p.get_y() << ") is invalid!\n";
                invalid_count++;
            }
        }

        if (invalid_count == 0) {
            std::cout << "All " << curve.points.size() 
                      << " points are valid!\n";
        } else {
            std::cout << "Found " << invalid_count 
                      << " invalid points!\n";
        }
    }

    // Пример 3: Операции с точками (исправленный вывод)
{
    std::cout << "\n=== Example 3: Point operations ===" << std::endl;
    mpz_class a = 1, b = 1, p = 23;
    Curve curve(a, b, p);
    curve.find_points();

    Point P = curve.points[1]; // (0, 1)
    Point Q = curve.points[2]; // (0, 22)

    Point sum = P + Q;
    Point double_P = P * 2;

    // Вывод с проверкой на бесконечность
    auto print_point = [](const std::string& name, const Point& p) {
        std::cout << name << ": ";
        if (p.isInfinity()) {
            std::cout << "INF";
        } else {
            std::cout << "(" << p.get_x() << ", " << p.get_y() << ")";
        }
        std::cout << "\n";
    };

    print_point("P", P);
    print_point("Q", Q);
    print_point("P + Q", sum);
    print_point("2P", double_P);
}

    return 0;
}