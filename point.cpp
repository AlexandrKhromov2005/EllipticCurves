#include "point.hpp"
#include "curve.hpp"
#include "algorithms_for_primes.hpp"

// Конструкторы
Point::Point(Curve* crv) : x(0), y(0), is_infinity(true), crv(crv) {}
Point::Point(const mpz_class& x, const mpz_class& y, Curve* crv) 
    : x(mod(x, crv->get_p())), 
      y(mod(y, crv->get_p())), 
      is_infinity(false), 
      crv(crv) {}

// Геттеры
mpz_class Point::get_x() const { return x; }
mpz_class Point::get_y() const { return y; }
bool Point::isInfinity() const { return is_infinity; }
Curve* Point::get_curve() const { return crv; }

bool Point::operator==(const Point& other) const {
    if (isInfinity() && other.isInfinity()) return true;
    if (isInfinity() || other.isInfinity()) return false;
    return (x == other.x) && (y == other.y) && (crv == other.crv);
}

Point Point::operator+(const Point& other) const {
    if (crv != other.crv) throw std::invalid_argument("Points are on different curves");

    if (isInfinity()) return other;
    if (other.isInfinity()) return *this;

    mpz_class p = crv->get_p();
    mpz_class x1 = x, y1 = y;
    mpz_class x2 = other.x, y2 = other.y;

    if (x1 == x2 && y1 != y2) return Point(crv);

    mpz_class lambda;
    if (*this == other) { 
        if (y1 == 0) return Point(crv); 
        mpz_class numerator = mod(3 * mod(x1 * x1, p) + crv->get_a(), p);
        mpz_class denominator = mod(2 * y1, p);
        mpz_class inv_denominator = mod_inverse(denominator, p);
        lambda = mod(numerator * inv_denominator, p);
    } else { 
        mpz_class numerator = mod(y2 - y1, p);
        mpz_class denominator = mod(x2 - x1, p);
        mpz_class inv_denominator = mod_inverse(denominator, p);
        lambda = mod(numerator * inv_denominator, p);
    }

    mpz_class x3 = mod(lambda * lambda - x1 - x2, p);
    mpz_class y3 = mod(lambda * (x1 - x3) - y1, p);

    return Point(x3, y3, crv);
}

Point Point::operator*(const mpz_class& k) const {
    Point result(crv);
    Point base = *this;
    mpz_class exponent = k;

    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result = result + base;
        }
        base = base + base;
        exponent /= 2;
    }

    return result;
}