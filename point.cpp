#include "point.hpp"
#include "curve.hpp"
#include <stdexcept>

int extendedEuclid(int a, int b, int& x, int& y) {
    if (b == 0) { x = 1; y = 0; return a; }
    int gcd = extendedEuclid(b, a % b, y, x);
    y -= (a / b) * x;
    return gcd;
}

int modInverse(int a, int p) {
    int x, y, gcd = extendedEuclid(a, p, x, y);
    return (gcd != 1) ? -1 : (x % p + p) % p;
}

Point::Point(int x, int y, bool is_inf, const Curve* crv) : x(x), y(y), is_infinity(is_inf), curve(crv) {}

bool Point::isInfinity() const { return is_infinity; }
int Point::getX() const { return x; }
int Point::getY() const { return y; }

Point Point::operator+(const Point& Q) const {
    if (is_infinity) return Q;
    if (Q.is_infinity) return *this;
    int p = curve->getP();
    if (x == Q.x && (y + Q.y) % p == 0) return Point(0, 0, true, curve);
    int s, numerator, denominator, inv_denominator;
    if (*this == Q) {
        numerator = (3 * x * x + curve->getA()) % p;
        denominator = (2 * y) % p;
    } else {
        numerator = (Q.y - y) % p;
        denominator = (Q.x - x) % p;
    }
    inv_denominator = modInverse(denominator, p);
    if (inv_denominator == -1) throw std::logic_error("No inverse");
    s = (numerator * inv_denominator) % p;
    int x3 = (s * s - x - Q.x) % p;
    int y3 = (s * (x - x3) - y) % p;
    return Point((x3 + p) % p, (y3 + p) % p, false, curve);
}

Point Point::multiply(int k) const {
    if (is_infinity) return *this;
    int order = curve->getOrder();
    k = k % order;
    Point result(0, 0, true, curve);
    Point current = *this;
    while (k > 0) {
        if (k % 2 == 1) result = result + current;
        current = current + current;
        k /= 2;
    }
    return result;
}

bool operator==(const Point& P, const Point& Q) {
    return (P.is_infinity && Q.is_infinity) || (P.x == Q.x && P.y == Q.y && P.curve == Q.curve);
}