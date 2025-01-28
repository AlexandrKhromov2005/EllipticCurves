#include "curve.hpp"
#include "point.hpp"
#include <stdexcept>
#include <set>

bool isPrime(int n, int k = 5) {
    if (n <= 1) return false;
    if (n <= 3) return true;
    int d = n - 1, s = 0;
    while (d % 2 == 0) { d /= 2; s++; }
    for (int i = 0; i < k; ++i) {
        int a = 2 + rand() % (n - 3);
        int x = 1, pow = d;
        while (pow--) x = (x * a) % n; 
        if (x == 1 || x == n - 1) continue;
        for (int r = 1; r < s; ++r) {
            x = (x * x) % n;
            if (x == n - 1) break;
        }
        if (x != n - 1) return false;
    }
    return true;
}

int Curve::curveEquation(int x) const {
    x = (x % p + p) % p; 
    return (x * x * x + a * x + b) % p;
}

std::vector<int> Curve::findYValues(int value) const {
    std::vector<int> yValues;
    value = (value % p + p) % p;
    for (int y = 0; y < p; ++y) {
        if ((y * y) % p == value) yValues.push_back(y);
    }
    return yValues;
}

Curve::Curve(int a, int b, int p) : a(a % p), b(b % p), p(p) {
    if (p <= 3) throw std::invalid_argument("p > 3 required");
    if (!isPrime(p)) throw std::invalid_argument("p must be prime");
    int discriminant = (-16 * (4 * a * a * a + 27 * b * b)) % p;
    if (discriminant == 0) throw std::invalid_argument("Singular curve");

    std::set<std::pair<int, int>> unique_points;
    for (int x = 0; x < p; ++x) {
        int rhs = curveEquation(x);
        std::vector<int> yCandidates = findYValues(rhs);
        for (int y : yCandidates) {
            if (y != 0 && y != p - y) {
                unique_points.insert({x, std::min(y, p - y)});
                unique_points.insert({x, std::max(y, p - y)});
            } else {
                unique_points.insert({x, y});
            }
        }
    }

    for (const auto& pt : unique_points) {
        points.emplace_back(pt.first, pt.second, false, this);
    }
    points.emplace_back(0, 0, true, this);
}

int Curve::getA() const { return a; }
int Curve::getB() const { return b; }
int Curve::getP() const { return p; }
const std::vector<Point>& Curve::getPoints() const { return points; }
int Curve::getOrderOfGroup() const { return points.size(); }

void Curve::getOrderOfPoints() const {
    orders.clear();
    for (const auto& point : points) {
        if (point.isInfinity()) {
            orders.push_back(1);
            continue;
        }
        int order = 1;
        Point current = point;
        while (!current.isInfinity()) {
            current = current + point;
            order++;
        }
        orders.push_back(order);
    }
}

const std::vector<int>& Curve::getOrders() const {
    return orders;
}

bool Curve::isCyclic() const {
    for (const auto& order : orders) {
        if (order == points.size()) return true;
    }
    return false;
}

void Curve::findPrimeSubgroups() const {
    primeSubgroups.clear();
    int groupOrder = getOrderOfGroup();

    // Факторизация порядка группы
    std::vector<int> primes;
    int temp = groupOrder;
    for (int i = 2; i * i <= temp; ++i) {
        if (temp % i == 0) {
            primes.push_back(i);
            while (temp % i == 0) temp /= i;
        }
    }
    if (temp > 1) primes.push_back(temp);

    // Поиск точек с простым порядком
    for (int prime : primes) {
        for (size_t i = 0; i < points.size(); ++i) {
            if (!points[i].isInfinity() && orders[i] == prime) {
                primeSubgroups.push_back(prime);
                break;
            }
        }
    }
}

const std::vector<int>& Curve::getPrimeSubgroups() const {
    return primeSubgroups;
}