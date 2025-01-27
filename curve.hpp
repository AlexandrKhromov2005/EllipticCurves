#ifndef CURVE_HPP
#define CURVE_HPP

#include <vector>
class Point;

class Curve {
private:
    int a;
    int b;
    int p;
    std::vector<Point> points;
    int curveEquation(int x) const;
    std::vector<int> findYValues(int value) const;

public:
    Curve(int a, int b, int p);
    int getA() const;
    int getB() const;
    int getP() const;
    const std::vector<Point>& getPoints() const;
    int getOrder() const;
};

#endif