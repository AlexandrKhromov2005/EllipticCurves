#include "curve.hpp"
#include <iostream>

Curve::Curve(mpz_class a, mpz_class b, mpz_class p) : a(a), b(b), p(p) {}
Curve::Curve() : a(0), b(0), p(0) {}

mpz_class Curve::get_a() { return a; }
mpz_class Curve::get_b() { return b; }
mpz_class Curve::get_p() { return p; }

void Curve::find_points() {
    points.push_back(Point(this));  
    for (mpz_class x = 0; x < p; ++x) {
        mpz_class y_pow_2 = (x * x * x + a * x + b) % p;
        if (y_pow_2 == 0) {
            points.push_back(Point(x, 0, this));  // добавляем точку (x, 0)
            continue;
        }

        if (legendre_symbol(y_pow_2, p) == 1) {
            points.push_back(Point(x, mod_sqrt(y_pow_2, p), this));  // добавляем точку (x, y)
            points.push_back(Point(x, p - mod_sqrt(y_pow_2, p), this));  // добавляем точку (x, -y)
        }
    }
}

void Curve::print_points() {
    mpz_class i = 0;
    for (Point point : points) {
        if (point.isInfinity()) {
            std::cout << "Infinity point" << std::endl;
            ++i;
            continue;
        }
        std::cout << "Point " << i << ": (" << point.get_x() << " , " << point.get_y() << ")" << std::endl;
        ++i;
    }
}

bool Curve::is_on_curve(const Point& P) const {
    if (P.isInfinity()) return true; // Бесконечно удаленная точка всегда на кривой
    mpz_class x = P.get_x();
    mpz_class y = P.get_y();

    // Проверяем уравнение кривой: y² ≡ x³ + a*x + b mod p
    mpz_class lhs = (y * y) % p;         // Левая часть: y²
    mpz_class rhs = (x*x*x + a*x + b) % p; // Правая часть: x³ + ax + b
    return lhs == rhs;
}
