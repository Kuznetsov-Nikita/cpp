#include <cmath>
#include <iostream>
#include <string>
#include <vector>

const int POW = 6;
const int BASE = 1000000;

class BigInteger {
private:
    int length = 0;
    int sign = 1;
    std::vector<int> number;

    friend std::istream& operator>>(std::istream&, BigInteger&);
    friend bool operator==(const BigInteger&, const BigInteger&);
    friend bool operator!=(const BigInteger&, const BigInteger&);
    friend bool operator<(const BigInteger&, const BigInteger&);
    friend bool operator>(const BigInteger&, const BigInteger&);
    friend bool operator>=(const BigInteger&, const BigInteger&);

    void clear() {
        number = {};
        length = 0;
        sign = 1;
    }

public:
    BigInteger() {
        number.push_back(0);
    }

    BigInteger(int src_number) {
        sign = (src_number < 0) ? -1 : 1;
        src_number = abs(src_number);

        do {
            number.push_back(src_number % BASE);
            src_number /= BASE;
            ++length;
        } while (src_number != 0);
    }

    BigInteger(size_t src_number) {
        sign = 1;

        do {
            number.push_back(src_number % BASE);
            src_number /= BASE;
            ++length;
        } while (src_number != 0);
    }

    BigInteger(unsigned long long int src_number) {
        sign = 1;

        do {
            number.push_back(src_number % BASE);
            src_number /= BASE;
            ++length;
        } while (src_number != 0);
    }

    BigInteger(const BigInteger& src_big_int): length(src_big_int.length), sign(src_big_int.sign), number(src_big_int.number) {}

    BigInteger& operator=(const BigInteger& src_big_int) {
        BigInteger copy = src_big_int;
        swap(copy);
        return *this;
    }

    void swap(BigInteger& big_int) {
        std::swap(length, big_int.length);
        std::swap(sign, big_int.sign);
        std::swap(number, big_int.number);
    }

    explicit operator bool() {
        return *this != 0;
    }

    std::string toString() const {
        std::string str_big_int = "";

        for (int i = length - 1; i >= 0; --i) {
            int j = 0;
            int tmp = number[i];
            std::string tmp_str_big_int = "";
            while (j < POW) {
                char dig = tmp % 10 + '0';
                tmp_str_big_int = dig + tmp_str_big_int;
                tmp /= 10;
                ++j;
            }

            str_big_int += tmp_str_big_int;
        }

        std::string ans = "";

        if (sign == -1) {
            ans += '-';
        }

        int i = 0;
        while (i < str_big_int.size() - 1 && str_big_int[i] == '0') {
            ++i;
        }

        while (i < str_big_int.size()) {
            ans += str_big_int[i];
            ++i;
        }

        return ans;
    }

    BigInteger& operator+=(const BigInteger& big_int) {
        if (big_int != 0 && *this != 0) {
            if (sign == big_int.sign) {
                int carry = 0;

                for (int i = 0; i < std::max(length, big_int.length) || carry != 0; ++i) {
                    int new_digit = carry;
                    if (i < length) {
                        new_digit += number[i];
                    }
                    if (i < big_int.length) {
                        new_digit += big_int.number[i];
                    }

                    if (i >= length) {
                        number.push_back(0);
                        ++length;
                    }

                    number[i] = new_digit % BASE;
                    carry = new_digit / BASE;
                }
            } else {
                sign *= -1;
                *this -= big_int;
                sign *= -1;
            }
        } else {
            if (*this == 0) {
                *this = big_int;
            }
        }

        return *this;
    }
    
    BigInteger& operator-=(const BigInteger& big_int) {
        if (big_int != 0 && *this != 0) {
            if (big_int == *this) {
                *this = 0;
                return *this;
            }

            if (sign == big_int.sign) {
                if((sign >= 0 && *this < big_int) || (sign < 0 && *this > big_int)) {
                    BigInteger tmp = big_int;
                    tmp -= *this;
                    *this = tmp;
                    sign = -sign;

                    return *this;
                }

                int carry = 0;
                int i;

                for (i = 0; i < std::max(length, big_int.length); ++i) {
                    int new_digit = carry;

                    if (i < length) {
                        new_digit += number[i];
                    }
                    if (i < big_int.length) {
                        new_digit -= big_int.number[i];
                    }
                    if (new_digit < 0) {
                        new_digit += BASE;
                        carry = -1;
                    } else {
                        carry = 0;
                    }

                    if (length <= i) {
                        number.push_back(0);
                        ++length;
                    }
                    number[i] = new_digit;
                }

                if (carry != 0) {
                    if (i != 0) {
                        number[0] = BASE - number[0];
                    }

                    length = (i != 0) ? 1 : 0;

                    for (int j = 1; j < i; ++j) {
                        number[j] = BASE - 1 - number[j];

                        if (number[i] != 0) {
                            length = j + 1;
                        }
                    }

                    sign *= -1;
                }

                while (length > 1 && number[length - 1] == 0) {
                    --length;
                    number.pop_back();
                }
            } else {
                sign *= -1;
                *this += big_int;
                sign *= -1;
            }
        } else {
            if (*this == 0) {
                *this = big_int;
                sign = -big_int.sign;
            }
        }

        return *this;
    }

    BigInteger& operator*=(const BigInteger& big_int) {
        if (*this == 0 || big_int == 0) {
            *this = 0;
        } else {
            std::vector<int> tmp(length + big_int.length);

            sign *= big_int.sign;

            int i, j = 0;

            for (i = 0; i < big_int.length; ++i) {
                if (big_int.number[i] != 0) {
                    int carry = 0;

                    for (j = 0; j < length || carry != 0; ++j) {
                        int new_digit = tmp[i + j] + ((j < length) ? big_int.number[i] * number[j] : 0) + carry;

                        tmp[i + j] = new_digit % BASE;
                        carry = new_digit / BASE;
                    }
                }
            }

            length = i + j - 1;
            number = tmp;
        }

        return *this;
    }

    BigInteger& operator/=(const BigInteger& big_int) {
        if (*this != 0) {
            sign *= big_int.sign;

            BigInteger tmp = 0;
            BigInteger n = ((big_int.sign > 0) ? big_int : -big_int);

            for (int i = length - 1; i >= 0; --i) {
                tmp *= BASE;
                tmp += number[i];
                number[i] = 0;

                while (tmp >= n) {
                    tmp -= n;
                    ++number[i];
                }
            }

            while (length > 1 && number[length - 1] == 0) {
                --length;
                number.pop_back();
            }
        }

        if (length == 1 && number[0] == 0) {
            sign = 1;
        }

        return *this;
    }

    BigInteger& operator%=(const BigInteger& big_int) {
        if (*this != 0) {
            int s1 = sign, s2 = big_int.sign;

            sign *= big_int.sign;

            BigInteger tmp = 0;
            BigInteger n = ((big_int.sign > 0) ? big_int : -big_int);

            for (int i = length - 1; i >= 0; --i) {
                tmp *= BASE;
                tmp += number[i];
                number[i] = 0;

                while (tmp >= n) {
                    tmp -= n;
                    ++number[i];
                }
            }

            while (length > 1 && number[length - 1] == 0) {
                --length;
                number.pop_back();
            }

            *this = tmp;

            if ((s1 == -1 && s2 == -1) || (s1 == -1 && s2 == 1)) {
                sign = -1;
            }

            if (length == 1 && number[0] == 0) {
                sign = 1;
            }
        }

        return *this;
    }

    BigInteger operator-() const {
        if (*this != 0) {
            BigInteger copy = *this;
            copy.sign = (sign == 1) ? -1 : 1;
            return copy;
        }
        return *this;
    }

    BigInteger& operator++() {
        *this += 1;
        return *this;
    }

    BigInteger& operator--() {
        *this -= 1;
        return *this;
    }

    BigInteger operator++(int) {
        BigInteger copy = *this;
        ++*this;
        return copy;
    }

    BigInteger operator--(int) {
        BigInteger copy = *this;
        --*this;
        return copy;
    }
};

BigInteger operator+(const BigInteger& big_int_1, const BigInteger& big_int_2) {
    BigInteger copy = big_int_1;
    copy += big_int_2;
    return copy;
}

BigInteger operator-(const BigInteger& big_int_1, const BigInteger& big_int_2) {
    BigInteger copy = big_int_1;
    copy -= big_int_2;
    return copy;
}

BigInteger operator*(const BigInteger& big_int_1, const BigInteger& big_int_2) {
    BigInteger copy = big_int_1;
    copy *= big_int_2;
    return copy;
}

BigInteger operator/(const BigInteger& big_int_1, const BigInteger& big_int_2) {
    BigInteger copy = big_int_1;
    copy /= big_int_2;
    return copy;
}

BigInteger operator%(const BigInteger& big_int_1, const BigInteger& big_int_2) {
    BigInteger copy = big_int_1;
    copy %= big_int_2;
    return copy;
}

bool operator==(const BigInteger& big_int_1, const BigInteger& big_int_2) {
    if (big_int_1.sign != big_int_2.sign || big_int_1.length != big_int_2.length) {
        return false;
    } else {
        bool flag = true;

        for (int i = 0; i < big_int_1.length; ++i) {
            if (big_int_1.number[i] != big_int_2.number[i]) {
                flag = false;
                break;
            }
        }

        return flag;
    }
}

bool operator!=(const BigInteger& big_int_1, const BigInteger& big_int_2) {
    return !(big_int_1 == big_int_2);
}

bool operator<(const BigInteger& big_int_1, const BigInteger& big_int_2) {
    if (big_int_1.sign != big_int_2.sign) {
        return big_int_1.sign < big_int_2.sign;
    } else {
        if (big_int_1.length != big_int_2.length) {
            return (big_int_1.length < big_int_2.length) == (big_int_1.sign == 1);
        } else {
            const BigInteger* big_int_pointer_1;
            const BigInteger* big_int_pointer_2;
            if (big_int_1.sign > 0) {
                big_int_pointer_1 = &big_int_1;
                big_int_pointer_2 = &big_int_2;
            } else {
                big_int_pointer_1 = &big_int_2;
                big_int_pointer_2 = &big_int_1;
            }

            bool flag = false;

            for (int i = big_int_1.length - 1; i >= 0; --i) {
                if (big_int_pointer_1->number[i] < big_int_pointer_2->number[i]) {
                    flag = true;
                    break;
                } else if (big_int_pointer_1->number[i] > big_int_pointer_2->number[i]) {
                    break;
                }
            }

            return flag;
        }
    }
}

bool operator>(const BigInteger& big_int_1, const BigInteger& big_int_2) {
    return big_int_2 < big_int_1;
}

bool operator<=(const BigInteger& big_int_1, const BigInteger& big_int_2) {
    return (big_int_1 == big_int_2) || (big_int_1 < big_int_2);
}

bool operator>=(const BigInteger& big_int_1, const BigInteger& big_int_2) {
    return big_int_2 <= big_int_1;
}

std::istream& operator>>(std::istream& in, BigInteger& big_int) {
    big_int.clear();

    std::string in_big_int;
    in >> in_big_int;

    if (!in_big_int.empty()) {
        int first_no_null, size = in_big_int.size();// start position

        if (!isdigit(in_big_int[0])) {
            big_int.sign = (in_big_int[0] == '+') ? 1 : -1;
            first_no_null = 1;
        } else {
            first_no_null = 0;
        }

        while (first_no_null < size && in_big_int[first_no_null] == '0') {
            ++first_no_null;
        }

        if (first_no_null == size) {
            big_int.number.push_back(0);
            ++big_int.length;
        } else {
            for (int i = size - 1; i >= first_no_null; i = i - POW) {
                int next_digit = 0;
                int k = 0, j = i, d = 1;
                while (j - k >= first_no_null && k < POW) {
                    next_digit = next_digit + (in_big_int[j - k] - '0') * d;
                    ++k;
                    d *= 10;
                }

                big_int.number.push_back(next_digit);
                ++big_int.length;
            }
        }
    }

    if (big_int.length == 1 && big_int.number[0] == 0) {
        big_int.sign = 1;
    }

    return in;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& big_int) {
    out << big_int.toString();
    return out;
}

BigInteger operator""_bi(unsigned long long x) {
    return BigInteger(x);
}

BigInteger gcd(BigInteger a, BigInteger b) {
    while (a != 0 && b != 0) {
        if (a > b) {
            a = a % b;
        } else {
            b = b % a;
        }
    }

    return a + b;
}

class Rational {
private:
    int sign = 1;
    BigInteger numerator, denominator;

    friend bool operator==(const Rational&, const Rational&);
    friend bool operator<(const Rational&, const Rational&);

    void normalize() {
        BigInteger gcd_num_den = gcd(numerator, denominator);
        numerator /= gcd_num_den;
        denominator /= gcd_num_den;
    }

public:
    Rational() = default;

    Rational(const BigInteger& numerator, const BigInteger& denominator = 1): sign(((numerator < 0 && denominator < 0) || (numerator >= 0 && denominator >= 0)) ? 1 : -1),
                                                                              numerator(((numerator >= 0) ? numerator : -numerator)),
                                                                              denominator(((denominator >= 0) ? denominator : -denominator)) {
        normalize();
    }

    Rational(int numerator, int denominator = 1): sign(((numerator < 0 && denominator < 0) || (numerator >= 0 && denominator >= 0)) ? 1 : -1),
                                                  numerator(((numerator >= 0) ? numerator : -numerator)),
                                                  denominator(((denominator >= 0) ? denominator : -denominator)) {
        normalize();
    }

    Rational& operator+=(const Rational& rational) {
        if (sign == rational.sign) {
            BigInteger gcd_denominator = gcd(denominator, rational.denominator);
            numerator = numerator * (rational.denominator / gcd_denominator) + rational.numerator * (denominator / gcd_denominator);
            denominator = denominator * (rational.denominator / gcd_denominator);

            normalize();
        } else {
            *this -= -rational;
        }

        return *this;
    }

    Rational& operator-=(const Rational& rational) {
        if (sign == rational.sign) {
            BigInteger gcd_denominator = gcd(denominator, rational.denominator);
            numerator = numerator * (rational.denominator / gcd_denominator) - rational.numerator * (denominator / gcd_denominator);
            denominator = denominator * (rational.denominator / gcd_denominator);

            if (numerator < 0) {
                numerator = -numerator;
                sign *= -1;
            }
            if (numerator == 0) {
                sign = 1;
            }

            normalize();

        } else {
            *this += -rational;
        }

        return *this;
    }

    Rational& operator*=(const Rational& rational) {
        sign = (sign == rational.sign) ? 1 : -1;

        numerator *= rational.numerator;
        denominator *= rational.denominator;

        normalize();

        if (numerator == 0) {
            sign = 1;
        }

        return *this;
    }

    Rational& operator/=(const Rational& rational) {
        sign = (sign == rational.sign) ? 1 : -1;

        numerator *= rational.denominator;
        denominator *= rational.numerator;

        normalize();

        if (numerator == 0) {
            sign = 1;
        }

        return *this;
    }

    Rational operator-() const {
        Rational copy = *this;
        copy.sign = (sign == 1) ? -1 : 1;
        return copy;
    }

    std::string toString() const {
        std::string str_rational = "";

        if (sign == -1) {
            str_rational += '-';
        }

        if (denominator != 1) {
            str_rational += numerator.toString() + '/' + denominator.toString();
        } else {
            str_rational += numerator.toString();
        }

        return str_rational;
    }

    std::string asDecimal(size_t precision = 0) const {
        std::string decimal = ((sign == -1) ? "-" : "") + (numerator / denominator).toString();

        if (precision != 0) {
            decimal += '.';
            unsigned long long int p = pow(10, precision);

            decimal += (numerator % denominator * p / denominator).toString();
        }

        return decimal;
    }

    explicit operator double() const {
        std::string decimal = this->asDecimal(18);

        double answer = decimal[0];
        for (size_t i = 1; i < decimal.length(); ++i) {
            if (decimal[i] != '.') {
                answer = answer * 10 + decimal[i];
            }
        }

        if (decimal.find('.') != std::string::npos) {
            for (size_t i = 0; i < decimal.length() - decimal.find('.'); ++i) {
                answer *= 0.1;
            }
        }

        return answer;
    }
};

Rational operator+(const Rational& rational_1, const Rational& rational_2) {
    Rational copy = rational_1;
    copy += rational_2;
    return copy;
}

Rational operator-(const Rational& rational_1, const Rational& rational_2) {
    Rational copy = rational_1;
    copy -= rational_2;
    return copy;
}

Rational operator*(const Rational& rational_1, const Rational& rational_2) {
    Rational copy = rational_1;
    copy *= rational_2;
    return copy;
}

Rational operator/(const Rational& rational_1, const Rational& rational_2) {
    Rational copy = rational_1;
    copy /= rational_2;
    return copy;
}

bool operator==(const Rational& rational_1, const Rational& rational_2) {
    return (rational_1.sign == rational_2.sign) && (rational_1.numerator == rational_2.numerator) && (rational_1.denominator == rational_2.denominator);
}

bool operator!=(const Rational& rational_1, const Rational& rational_2) {
    return !(rational_1 == rational_2);
}

bool operator<(const Rational& rational_1, const Rational& rational_2) {
    if (rational_1.sign < rational_2.sign) {
        return true;
    } else if (rational_1.sign > rational_2.sign) {
        return false;
    } else {
        BigInteger gcd_denominator = gcd(rational_1.denominator, rational_2.denominator);
        return rational_1.numerator * (rational_2.denominator / gcd_denominator) < rational_2.numerator * (rational_1.denominator / gcd_denominator);
    }
}

bool operator>(const Rational& rational_1, const Rational& rational_2) {
    return rational_2 < rational_1;
}

bool operator>=(const Rational& rational_1, const Rational& rational_2) {
    return rational_1 > rational_2 || rational_1 == rational_2;
}

bool operator<=(const Rational& rational_1, const Rational& rational_2) {
    return rational_2 >= rational_1;
}
