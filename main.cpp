#include <iostream>
#include "curve.hpp"
#include "point.hpp"

int main() {
    try {
        int p, a, b;
        std::cout << "Enter curve parameters (p a b): ";
        std::cin >> p >> a >> b;
        Curve curve(a, b, p);
        std::cout << "Curve created. Order: " << curve.getOrder() << "\n";
        
        // Вывод всех точек
        std::cout << "All points:\n";
        for (const Point& pt : curve.getPoints()) {
            if (pt.isInfinity()) {
                std::cout << "O (infinity)\n";
            } else {
                std::cout << "(" << pt.getX() << ", " << pt.getY() << ")\n";
            }
        }

        if (!curve.getPoints().empty()) {
            Point P = curve.getPoints()[1];
            int k;
            std::cout << "Enter scalar k: ";
            std::cin >> k;
            Point Q = P.multiply(k);
            std::cout << "Result " << k << "P: (" << Q.getX() << ", " << Q.getY() << ")\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
    return 0;
}