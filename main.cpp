#include <iostream>
#include <string>
#include <stdexcept>
#include "curve.hpp"

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <a> <b> <p>" << std::endl;
        std::cerr << "Example: " << argv[0] << " 2 3 97" << std::endl;
        return 1;
    }

    try {
        mpz_class a(argv[1]);
        mpz_class b(argv[2]);
        mpz_class p(argv[3]);

        Curve curve(a, b, p);
        curve.find_points();
        for (Point &pt : curve.points){
            pt.calculate_order();
        }
        std::cout << "Elliptic curve: y^2 = x^3 + " << a << "x + " << b << " (mod " << p << ")" << std::endl;
        std::cout << "Found points:" << std::endl;
        curve.print_points();

        std::string command;
        while (true) {
            std::cout << "\nEnter command (add - addition, mul - multiplication, exit - exit): ";
            std::cin >> command;

            if (command == "exit") {
                break;
            } 
            else if (command == "add") {
                std::size_t idx1, idx2;
                std::cout << "Enter index of the first point: ";
                std::cin >> idx1;
                std::cout << "Enter index of the second point: ";
                std::cin >> idx2;

                if (idx1 >= curve.points.size() || idx2 >= curve.points.size()) {
                    std::cout << "Invalid point index." << std::endl;
                    continue;
                }

                Point res = curve.points[idx1] + curve.points[idx2];
                std::cout << "Addition result: ";
                if (res.isInfinity()) {
                    std::cout << "Point at infinity" << std::endl;
                } else {
                    std::cout << "(" << res.get_x() << ", " << res.get_y() << ")" << std::endl;
                }
            } 
            else if (command == "mul") {
                std::size_t idx;
                std::string scalar_str;
                std::cout << "Enter point index: ";
                std::cin >> idx;
                if (idx >= curve.points.size()) {
                    std::cout << "Invalid point index." << std::endl;
                    continue;
                }
                std::cout << "Enter scalar: ";
                std::cin >> scalar_str;
                mpz_class k(scalar_str);

                Point res = curve.points[idx] * k;
                std::cout << "Multiplication result: ";
                if (res.isInfinity()) {
                    std::cout << "Point at infinity" << std::endl;
                } else {
                    std::cout << "(" << res.get_x() << ", " << res.get_y() << ")" << std::endl;
                }
            } 
            else {
                std::cout << "Unknown command. Try again." << std::endl;
            }
        }
    }
    catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
