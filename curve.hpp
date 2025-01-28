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
    mutable std::vector<int> orders;
    mutable std::vector<int> primeSubgroups; // Добавлено
    int curveEquation(int x) const;
    std::vector<int> findYValues(int value) const;

public:
    Curve(int a, int b, int p);
    int getA() const;
    int getB() const;
    int getP() const;
    const std::vector<Point>& getPoints() const;
    int getOrderOfGroup() const;
    void getOrderOfPoints() const;
    const std::vector<int>& getOrders() const;
    void findPrimeSubgroups() const; // Добавлено
    const std::vector<int>& getPrimeSubgroups() const; // Добавлено
    bool isCyclic() const; // Добавлено
};

#endif