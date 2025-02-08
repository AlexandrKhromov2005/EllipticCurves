#include "curve.hpp"
#include <iostream>

Curve::Curve(mpz_class a, mpz_class b, mpz_class p) : a(a), b(b), p(p) {}
Curve::Curve() : a(0), b(0), p(0) {}

mpz_class Curve::get_a() { return a; }
mpz_class Curve::get_b() { return b; }
mpz_class Curve::get_p() { return p; }

void Curve::find_points() {
    // Добавляем точку на бесконечности
    points.push_back(Point(this)); 
    orders.push_back(mpz_class(-1));  // Здесь можно поставить особое значение, например, -1

    for (mpz_class x = 0; x < p; ++x) {
        mpz_class y_pow_2 = (x * x * x + a * x + b) % p;
        if (y_pow_2 == 0) {
            Point pt(x, 0, this);
            points.push_back(pt);
            orders.push_back(pt.calculate_order()); 
            continue;
        }

        if (legendre_symbol(y_pow_2, p) == 1) {
            Point pt1(x, mod_sqrt(y_pow_2, p), this);
            points.push_back(pt1);
            orders.push_back(pt1.calculate_order());

            Point pt2(x, p - mod_sqrt(y_pow_2, p), this);
            points.push_back(pt2);
            orders.push_back(pt2.calculate_order());
        }
    }
}

void Curve::print_points() {
    mpz_class i = 0;
    for (const Point& point : points) {
        if (point.isInfinity()) {
            std::cout << "Infinity point" << std::endl;
            ++i;
            continue;
        }
        std::cout << "Point " << i << ": (" << point.get_x() << " , " << point.get_y() << ")"
                  << " Order: " << orders[i.get_ui()] << std::endl;
        ++i;
    }
}

bool Curve::is_on_curve(const Point& P) const {
    if (P.isInfinity()) return true; 
    mpz_class x = P.get_x();
    mpz_class y = P.get_y();

    mpz_class lhs = (y * y) % p;         
    mpz_class rhs = (x * x * x + a * x + b) % p; 
    return lhs == rhs;
}

void Curve::calculate_orders() {
    orders.clear(); 
    for (const Point& pt : points) {
        orders.push_back(pt.calculate_order());
    }
}
