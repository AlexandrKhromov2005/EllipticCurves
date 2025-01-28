#include <iostream>
#include <limits>
#include "curve.hpp"
#include "point.hpp"

int main() {
    int a, b, p;
    
    std::cout << "Enter elliptic curve parameters (a, b, p):\n";
    std::cout << "a: ";
    std::cin >> a;
    std::cout << "b: ";
    std::cin >> b;
    std::cout << "p (prime number > 3): ";
    std::cin >> p;

    try {
        Curve curve(a, b, p);
        curve.getOrderOfPoints();
        curve.findPrimeSubgroups();

        const std::vector<Point>& points = curve.getPoints();
        const std::vector<int>& orders = curve.getOrders();
        const std::vector<int>& primes = curve.getPrimeSubgroups();

        // Group info
        std::cout << "\nGroup structure:\n";
        std::cout << "Order: " << curve.getOrderOfGroup() << "\n";
        std::cout << "Is cyclic: " << (curve.isCyclic() ? "Yes" : "No") << "\n";
        if (!primes.empty()) {
            std::cout << "Prime-order subgroups: ";
            for (int p : primes) std::cout << p << " ";
            std::cout << "\n";
        }

        // Points list
        std::cout << "\nCurve points and their orders:\n";
        for (size_t i = 0; i < points.size(); ++i) {
            const Point& pt = points[i];
            if (pt.isInfinity()) {
                std::cout << "[" << i << "] Infinity point : order " << orders[i] << "\n";
            } else {
                std::cout << "[" << i << "] (" << pt.getX() << ", " << pt.getY() << ") : order " << orders[i] << "\n";
            }
        }

        // Menu
        while (true) {
            std::cout << "\nChoose an action:\n"
                      << "1. Multiply a point\n"
                      << "2. Show group info\n"
                      << "3. Exit\n"
                      << "> ";
            int choice;
            std::cin >> choice;
            if (choice == 3) break;

            if (choice == 1) {
                // Point multiplication
                std::cout << "Enter point index (0-" << points.size() - 1 << "): ";
                int index;
                std::cin >> index;
                
                if (std::cin.fail() || index < 0 || index >= static_cast<int>(points.size())) {
                    std::cerr << "Invalid point index!\n";
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    continue;
                }

                std::cout << "Enter multiplier: ";
                int k;
                std::cin >> k;
                if (std::cin.fail()) {
                    std::cerr << "Invalid multiplier!\n";
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    continue;
                }

                try {
                    Point result = points[index].multiply(k);
                    std::cout << "Result: ";
                    if (result.isInfinity()) {
                        std::cout << "Infinity point\n";
                    } else {
                        std::cout << "(" << result.getX() << ", " << result.getY() << ")\n";
                    }
                } catch (const std::exception& e) {
                    std::cerr << "Error: " << e.what() << "\n";
                }

            } else if (choice == 2) {
                // Group info
                std::cout << "Group order: " << curve.getOrderOfGroup() << "\n";
                std::cout << "Cyclic: " << (curve.isCyclic() ? "Yes" : "No") << "\n";
                std::cout << "Prime subgroups: ";
                for (int p : primes) std::cout << p << " ";
                std::cout << "\n";
            } else {
                std::cerr << "Invalid choice!\n";
            }
        }
    } 
    catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
    catch (...) {
        std::cerr << "Unknown error!\n";
    }

    return 0;
}