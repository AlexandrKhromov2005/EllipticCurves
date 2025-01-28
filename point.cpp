#include "point.hpp"
#include "curve.hpp"
#include <stdexcept>
#include <cmath>

int modInverse(int a, int p) {
    a = a % p;
    if (a < 0) a += p;
    int gcd, x;
    int old_r = p, r = a;
    int old_s = 0, s = 1;

    while (r != 0) {
        int quotient = old_r / r;
        int temp = r;
        r = old_r - quotient * r;
        old_r = temp;
        temp = s;
        s = old_s - quotient * s;
        old_s = temp;
    }

    if (old_r != 1) throw std::runtime_error("Inverse does not exist");
    return (old_s % p + p) % p;
}

Point::Point(int x, int y, bool is_inf, const Curve* crv) : x(x), y(y), is_infinity(is_inf), curve(crv) {}

bool Point::isInfinity() const { return is_infinity; }
int Point::getX() const { return x; }
int Point::getY() const { return y; }

Point Point::operator+(const Point& other) const {
    if (this->isInfinity()) return other;
    if (other.isInfinity()) return *this;
    if (this->x == other.x && this->y == other.y) {
        if (this->y == 0) return Point(0, 0, true, curve);
        int numerator = (3 * x * x + curve->getA()) % curve->getP();
        int denominator = (2 * y) % curve->getP();
        int s = (numerator * modInverse(denominator, curve->getP())) % curve->getP();
        int x3 = (s * s - 2 * x) % curve->getP();
        x3 = (x3 + curve->getP()) % curve->getP();
        int y3 = (s * (x - x3) - y) % curve->getP();
        y3 = (y3 + curve->getP()) % curve->getP();
        return Point(x3, y3, false, curve);
    } else {
        if (this->x == other.x) return Point(0, 0, true, curve);
        int numerator = (other.y - y) % curve->getP();
        int denominator = (other.x - x) % curve->getP();
        int s = (numerator * modInverse(denominator, curve->getP())) % curve->getP();
        int x3 = (s * s - x - other.x) % curve->getP();
        x3 = (x3 + curve->getP()) % curve->getP();
        int y3 = (s * (x - x3) - y) % curve->getP();
        y3 = (y3 + curve->getP()) % curve->getP();
        return Point(x3, y3, false, curve);
    }
}

Point Point::multiply(int k) const {
    Point result(0, 0, true, curve);
    Point addend = *this;
    k = k % (curve->getOrderOfGroup());
    if (k < 0) k += curve->getOrderOfGroup();
    while (k > 0) {
        if (k % 2 == 1) {
            result = result + addend;
        }
        addend = addend + addend;
        k /= 2;
    }
    return result;
}

bool operator==(const Point& P, const Point& Q) {
    return P.x == Q.x && P.y == Q.y && P.is_infinity == Q.is_infinity;
}