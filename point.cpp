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

Point Point::operator*(const mpz_class &k) const {
    if (k == mpz_class(-1)){return Point(x, -y, crv);}
    if (k == mpz_class(0)){return  Point(this->get_curve());}
    if (k == mpz_class(1)){return *this;}
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



/*mpz_class Point::calculate_order() const {
    if (this->isInfinity()) {
        return mpz_class(1);  
    }

    if (this->get_y() == 0) {
        return mpz_class(2);  
    }
    mpz_class n(static_cast<long>(this->crv->points.size()));
    mpz_class m = sqrt(n) + 1;
    Point temp = *this * m;
    Point qpoint(temp.get_x(), -temp.get_y(), temp.get_curve());
    std::vector<Point> p_points, q_points;
    p_points.push_back(*this);
    q_points.push_back(qpoint);
    for (mpz_class i = 1; i < m; ++i) {
    size_t idx = i.get_ui(); 
    p_points.push_back(p_points[idx - 1] + *this);
    q_points.push_back(q_points[idx - 1] + qpoint);
    }

    for (mpz_class i = 0; i < m; ++i) {
        for (mpz_class j = 0; j < m; ++j) {
            if (p_points[i.get_ui()] == q_points[j.get_ui()]) {
                return i * m + j + 1;
            }
        }
    }
    return -1;
}*/

/*mpz_class Point::calculate_order() const {
    // Если точка на бесконечности, её порядок равен 1
    if (this->isInfinity())
        return mpz_class(1);
    
    // Если y == 0, то при удвоении получаем нейтральный элемент
    if (this->get_y() == 0)
        return mpz_class(2);
    
    mpz_class N = static_cast<unsigned long>(this->crv->points.size());
    
    mpz_class sqrt_N;
    mpz_sqrt(sqrt_N.get_mpz_t(), N.get_mpz_t());
    mpz_class m;
    if (sqrt_N * sqrt_N < N)
        m = sqrt_N + 1;
    else
        m = sqrt_N;
    
    auto point_to_string = [this](const Point& pt) -> std::string {
        if (pt.isInfinity())
            return "inf";
        return pt.get_x().get_str() + "," + pt.get_y().get_str();
    };
    
    // Шаг 1. Вычисляем "маленькие шаги": для j = 0, 1, ..., m-1 вычисляем jP
    // и сохраняем в map: ключ – строковое представление точки, значение – j.
    std::map<std::string, mpz_class> babySteps;
    for (mpz_class j = 0; j < m; j++) {
        Point baby = (*this) * j;  // 0*P, 1*P, 2*P, ..., (m-1)*P
        std::string key = point_to_string(baby);
        // Сохраним только первый (наименьший) индекс j для данного представления точки
        if (babySteps.find(key) == babySteps.end()) {
            babySteps[key] = j;
        }
    }
    
    Point factor = (*this) * m;
    for (mpz_class i = 1; i < m; i++) {
        Point giant = factor * i;
        Point neg_giant = giant.isInfinity()
                            ? giant
                            : Point(giant.get_x(), mod(-giant.get_y(), this->crv->get_p()), this->crv);
        
        std::string key = point_to_string(neg_giant);
        if (babySteps.find(key) != babySteps.end()) {
            mpz_class j = babySteps[key];
            mpz_class order = i * m + j;
            if (order > 0)
                return order;
        }
    }
    
    // Если ничего не найдено – выполняем запасной перебор (хотя такой случай встречается редко)
    mpz_class order = 1;
    Point current = *this;
    while (!current.isInfinity()) {
        order++;
        current = current + *this;
        if (order > N)
            break;
    }
    return order;
}*/

mpz_class Point::calculate_order() const {
    // Если точка на бесконечности, её порядок равен 1
    if (this->isInfinity())
        return mpz_class(1);
    
    // Если y == 0, то удвоение даёт нейтральный элемент
    if (this->get_y() == 0)
        return mpz_class(2);
    
    // Порядок всей группы (количество точек на кривой)
    mpz_class N = static_cast<unsigned long>(this->crv->points.size());
    
    // Определяем параметр m для baby-step giant-step (m ≈ √N)
    mpz_class sqrt_N;
    mpz_sqrt(sqrt_N.get_mpz_t(), N.get_mpz_t());
    mpz_class m = (sqrt_N * sqrt_N < N) ? sqrt_N + 1 : sqrt_N;
    
    // Лямбда для получения строкового представления точки
    auto point_to_string = [this](const Point& pt) -> std::string {
        if (pt.isInfinity())
            return "inf";
        return pt.get_x().get_str() + "," + pt.get_y().get_str();
    };
    
    // Лямбда для уточнения (рефайна) кандидата на порядок.
    // Идея: если для некоторого простого делителя f от candidate
    // оказывается, что (candidate / f)*P = O, то можно уменьшить кандидат.
    auto refine_order = [this](mpz_class candidate) -> mpz_class {
        mpz_class order = candidate;
        // Перебираем возможные делители начиная с 2.
        // (Простой перебор; предполагается, что order не слишком велик.)
        for (mpz_class f = 2; f <= order; f++) {
            // Пока f делит order, пробуем уменьшить порядок
            while (order % f == 0) {
                mpz_class new_order = order / f;
                // Если new_order уже даёт нейтральный элемент, то уточняем порядок
                if (((*this) * new_order).isInfinity()) {
                    order = new_order;
                } else {
                    break;
                }
            }
        }
        return order;
    };
    
    // Этап 1. Вычисляем маленькие шаги:
    // Для j = 0, 1, …, m-1 вычисляем jP и сохраняем в map
    std::map<std::string, mpz_class> babySteps;
    for (mpz_class j = 0; j < m; j++) {
        Point baby = (*this) * j;  // 0*P, 1*P, 2*P, …, (m-1)*P
        std::string key = point_to_string(baby);
        if (babySteps.find(key) == babySteps.end()) {
            babySteps[key] = j;
        }
    }
    
    // Этап 2. Вычисляем giant-steps:
    Point factor = (*this) * m;
    for (mpz_class i = 1; i < m; i++) {
        Point giant = factor * i;
        // Берём -giant (для сравнения с baby steps)
        Point neg_giant = giant.isInfinity()
                          ? giant
                          : Point(giant.get_x(), mod(-giant.get_y(), this->crv->get_p()), this->crv);
        
        std::string key = point_to_string(neg_giant);
        if (babySteps.find(key) != babySteps.end()) {
            mpz_class j = babySteps[key];
            mpz_class candidate_order = i * m + j;
            // Уточняем найденный порядок, пытаясь «сдернуть» лишние множители
            candidate_order = refine_order(candidate_order);
            // Если найденный порядок является делителем общего порядка группы, возвращаем его
            if (candidate_order > 0 && (N % candidate_order == 0)) {
                return candidate_order;
            }
            // Если кандидат не удовлетворяет условию делимости, продолжаем поиск
        }
    }
    
    // Запасной перебор: если алгоритм baby-step/giant-step не дал корректного кандидата,
    // перебираем последовательное умножение, пока не найдём порядок.
    mpz_class order = 1;
    Point current = *this;
    while (!current.isInfinity()) {
        order++;
        current = current + *this;
        // Если порядок превышает общий порядок группы, прерываем цикл (это маловероятно)
        if (order > N)
            break;
    }
    order = refine_order(order);
    if (order > 0 && (N % order == 0)) {
        return order;
    }
    
    // Если ни один из методов не дал корректного результата, возвращаем полученное значение.
    return order;
}




