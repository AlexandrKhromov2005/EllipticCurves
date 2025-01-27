#ifndef POINT_HPP
#define POINT_HPP

#include "curve.hpp"
class Curve;

class Point {
private:
    int x;
    int y;
    bool is_infinity;
    const Curve* curve;

public:
    Point(int x, int y, bool is_inf, const Curve* crv);
    bool isInfinity() const;
    int getX() const;
    int getY() const;
    Point operator+(const Point& other) const;
    Point multiply(int k) const;
    friend bool operator==(const Point& P, const Point& Q);
};

#endif