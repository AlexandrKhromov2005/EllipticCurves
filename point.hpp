#ifndef POINT_HPP
#define POINT_HPP

#include <gmpxx.h>
#include <stdexcept>
#include "algorithms_for_primes.hpp"

class Curve;

class Point {
private:
    mpz_class x;
    mpz_class y;
    bool is_infinity;
    Curve* crv;

public:
    Point(Curve* crv);
    Point(const mpz_class& x, const mpz_class& y, Curve* crv);
    
    mpz_class get_x() const;
    mpz_class get_y() const;
    bool isInfinity() const;
    Curve* get_curve() const;

    Point operator+(const Point& other) const;
    Point operator*(const mpz_class& k) const;
    bool operator==(const Point& other) const;
};

#endif // POINT_HPP