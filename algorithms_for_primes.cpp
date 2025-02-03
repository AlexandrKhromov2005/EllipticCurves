#include "algorithms_for_primes.hpp"

mpz_class mod(const mpz_class &a, const mpz_class &p) {
    return (((a % p) + p) % p);
}

mpz_class mod_inverse(const mpz_class& a, const mpz_class& m) {
    mpz_class gcd, x, y;
    gcd = extended_gcd(a, m, x, y); // Добавлено присваивание gcd
    if (gcd != 1) throw std::runtime_error("Inverse does not exist");
    return mod(x, m);
}

mpz_class pow_mod(const mpz_class& base, const mpz_class& exponent, const mpz_class& modulus) {
    if (modulus == 1) return 0; // Ноль для modulus = 1
    mpz_class result = 1;
    mpz_class b = base % modulus;
    mpz_class e = exponent;

    while (e > 0) {
        if (e % 2 == 1) {
            result = (result * b) % modulus;
        }
        b = (b * b) % modulus;
        e /= 2;
    }
    return result;
}

mpz_class mod_sqrt(const mpz_class& a, const mpz_class& p) {
    if (a == 0) return 0; // Корень из нуля

    // Используем нашу функцию legendre_symbol вместо mpz_legendre
    if (legendre_symbol(a, p) != 1) return 0; // Не квадратичный вычет

    // Случай p = 2
    if (p == 2) return a;

    mpz_class q = p - 1;
    mpz_class s = 0;
    while (q % 2 == 0) {
        q /= 2;
        s++;
    }

    // Поиск невычета z с помощью нашей legendre_symbol
    mpz_class z = 2;
    while (legendre_symbol(z, p) != -1) {
        z++;
    }

    mpz_class c = pow_mod(z, q, p);
    mpz_class x = pow_mod(a, (q + 1) / 2, p);
    mpz_class t = pow_mod(a, q, p);
    mpz_class m = s;

    while (t != 1) {
        mpz_class tt = t;
        mpz_class i = 0;
        // Находим наименьшее i: t^(2^i) ≡ 1 mod p
        while (tt != 1 && i < m) {
            tt = pow_mod(tt, 2, p);
            i++;
        }

        // Вычисляем 2^(m-i-1) через GMP-функции
        mpz_class exponent = 1;
        mpz_class shift = m - i - 1;
        mpz_mul_2exp(exponent.get_mpz_t(), exponent.get_mpz_t(), mpz_get_ui(shift.get_mpz_t()));
        
        mpz_class b = pow_mod(c, exponent, p);
        x = (x * b) % p;
        t = (t * b % p * b) % p;
        c = pow_mod(b, 2, p);
        m = i;
    }

    return x;
}

mpz_class extended_gcd(const mpz_class& a, const mpz_class& b, mpz_class& x, mpz_class& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    mpz_class q = a / b;
    mpz_class r = a % b;
    mpz_class x1, y1;
    mpz_class gcd = extended_gcd(b, r, x1, y1);
    x = y1;
    y = x1 - q * y1;
    return gcd;
}

bool invert(mpz_class& result, const mpz_class& a, const mpz_class& m) {
    if (m <= 1) {
        return false; 
    }
    mpz_class gcd, x, y;
    gcd = extended_gcd(a, m, x, y);
    if (gcd != 1) {
        return false; 
    }

    x %= m;
    if (x < 0) {
        x += m;
    }
    result = x;
    return true;
}


int legendre_symbol(const mpz_class& a, const mpz_class& p) {
    // Обработка p = 2
    if (p == 2) {
        mpz_class tmp = a % 2;
        return (tmp == 0) ? 0 : 1;
    }

    // Проверка на четность p (символ Лежандра требует нечетного простого p)
    if (p % 2 == 0) {
        return 0; // Некорректный ввод
    }

    mpz_class tmp = a % p;
    if (tmp == 0) return 0;

    int result = 1;
    mpz_class n = tmp;
    mpz_class d = p;

    while (n != 0) {
        // Извлечение множителей 2 из n
        unsigned long t = mpz_scan1(n.get_mpz_t(), 0); // Найти степень 2
        if (t > 0) {
            mpz_class divisor;
            mpz_ui_pow_ui(divisor.get_mpz_t(), 2, t); // divisor = 2^t
            n /= divisor; // Явное деление вместо n >>= t

            // Проверка условия для множителя 2
            mpz_class mod8 = d % 8;
            if ((t % 2 == 1) && (mod8 == 3 || mod8 == 5)) {
                result *= -1;
            }
        }

        // Квадратичный закон взаимности
        mpz_class mod4_n = n % 4;
        mpz_class mod4_d = d % 4;
        if (mod4_d == 3 && mod4_n == 3) {
            result *= -1;
        }

        // Обмен n и d с взятием модуля
        std::swap(n, d);
        n %= d;
    }

    return (d == 1) ? result : 0;
}