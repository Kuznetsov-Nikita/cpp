#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "biginteger.h"

template <int N>
class Finite {
private:
    long long int number;

    template<int M>
    friend bool operator==(const Finite<M>&, const Finite<M>&);
    template <int M>
    friend bool operator==(const Finite<M>&, int);
    template<int M>
    friend bool operator!=(const Finite<M>&, const Finite<M>&);
    template<int M>
    friend bool operator<(const Finite<M>&, const Finite<M>&);
    template<int M>
    friend bool operator>(const Finite<M>&, const Finite<M>&);
    template<int M>
    friend bool operator<=(const Finite<M>&, const Finite<M>&);
    template<int M>
    friend bool operator>=(const Finite<M>&, const Finite<M>&);

    bool is_prime(int n) {
        int i = 2;

        while (i <= sqrt(n)) {
            if (n % i == 0) {
                return false;
            }
            ++i;
        }

        return true;
    }

public:
    Finite<N>(int n = 0) {
        if (n < 0) {
            number = N + n % N;
        } else {
            number = n % N;
        }
    }

    Finite<N>& operator+=(const Finite<N>& another) {
        number = (number + another.number) % N;
        return *this;
    }

    Finite<N>& operator-=(const Finite<N>& another) {
        number -= another.number;

        if (number < 0) {
            number = N + number;
        }

        return *this;
    }

    Finite<N>& operator*=(const Finite<N>& another) {
        number = number * another.number % N;
        return *this;
    }

    Finite<N>& operator++() {
        *this += Finite<N> (1);
        return *this;
    }

    Finite<N> operator++(int) const {
        Finite<N> copy = this;
        copy += Finite<N> (1);
        return copy;
    }

    Finite<N>& operator--() {
        *this -= Finite<N> (1);
        return *this;
    }

    Finite<N> operator--(int) const {
        Finite<N> copy = this;
        copy -= Finite<N> (1);
        return copy;
    }

    Finite<N>& operator/=(const Finite<N>& another) {
        long long int a = N, b = another.number;
        std::vector<long long int> a_div_b;
        while (a % b != 0) {
            a_div_b.push_back(a / b);
            long long int t = b;
            b = a % b;
            a = t;
        }

        long long int x = 0, y = 1;

        for (int i = a_div_b.size() - 1; i >= 0; --i) {
            int t = x;
            x = y;
            y = t - y * a_div_b[i];
        }

        if (y < 0) {
            y += N;
        }

        *this *= Finite<N>(y);

        return *this;
    }
};

template <int N>
Finite<N> operator+(const Finite<N>& finite_1, const Finite<N>& finite_2) {
    Finite<N> copy = finite_1;
    copy += finite_2;
    return copy;
}

template <int N>
Finite<N> operator-(const Finite<N>& finite_1, const Finite<N>& finite_2) {
    Finite<N> copy = finite_1;
    copy -= finite_2;
    return copy;
}

template <int N>
Finite<N> operator*(const Finite<N>& finite_1, const Finite<N>& finite_2) {
    Finite<N> copy = finite_1;
    copy *= finite_2;
    return copy;
}

template <int N>
Finite<N> operator/(const Finite<N>& finite_1, const Finite<N>& finite_2) {
    Finite<N> copy = finite_1;
    copy /= finite_2;
    return copy;
}

template <int N>
bool operator==(const Finite<N>& finite_1, const Finite<N>& finite_2) {
    return finite_1.number == finite_2.number;
}

template <int N>
bool operator==(const Finite<N>& finite, int another) {
    return finite.number == another;
}

template <int N>
bool operator!=(const Finite<N>& finite_1, const Finite<N>& finite_2) {
    return !(finite_1 == finite_2);
}

template <int N>
bool operator<(const Finite<N>& finite_1, const Finite<N>& finite_2) {
    return finite_1.number < finite_2.number;
}

template <int N>
bool operator>(const Finite<N>& finite_1, const Finite<N>& finite_2) {
    return finite_2 < finite_1;
}

template <int N>
bool operator<=(const Finite<N>& finite_1, const Finite<N>& finite_2) {
    return finite_1.number <= finite_2.number;
}

template <int N>
bool operator>=(const Finite<N>& finite_1, const Finite<N>& finite_2) {
    return finite_2 <= finite_1;
}

template <unsigned int N, unsigned int M = N, typename Field = Rational>
class Matrix {
protected:
    std::vector<std::vector<Field>> matrix;

public:
    Matrix<N, M, Field>() {
        if (N == M) {
            for (unsigned int i = 0; i < N; ++i) {
                matrix.push_back(std::vector<Field> ());

                for (unsigned int j = 0; j < N; ++j) {
                    if (i == j) {
                        matrix[i].push_back(1);
                    } else {
                        matrix[i].push_back(0);
                    }
                }
            }
        }
    }

    Matrix<N, M, Field>(std::vector<std::vector<Field>> src_matrix) {
        matrix = src_matrix;

        for (unsigned int i = 0; i < src_matrix.size(); ++i) {
            for (unsigned int j = src_matrix[i].size(); j < M; ++j) {
                matrix[i].push_back(0);
            }
        }

        for (unsigned int i = src_matrix.size(); i < N; ++i) {
            matrix.push_back(std::vector<Field> (0));
        }
    }

    Matrix<N, M, Field>(std::vector<std::vector<int>> src_matrix) {
        for (size_t i = 0; i < src_matrix.size(); ++i) {
            std::vector<Field> row;

            for (size_t j = 0; j < src_matrix[i].size(); ++j) {
                row.push_back(src_matrix[i][j]);
            }
            matrix.push_back(row);
        }

        for (unsigned int i = 0; i < src_matrix.size(); ++i) {
            for (unsigned int j = src_matrix[i].size(); j < M; ++j) {
                matrix[i].push_back(static_cast<Field>(0));
            }
        }

        for (unsigned int i = src_matrix.size(); i < N; ++i) {
            std::vector<Field> row;

            for (size_t j = 0; j < M; ++j) {
                row.push_back(static_cast<Field>(0));
            }
            matrix.push_back(row);
        }
    }

    std::vector<Field>& operator[](unsigned int index) {
        return matrix[index];
    }

    const std::vector<Field>& operator[](unsigned int index) const {
        return matrix[index];
    }

    std::vector<Field> getRow(unsigned int i) const {
        return matrix[i];
    }

    std::vector<Field> getColumn(unsigned int i) const {
        std::vector<Field> ans;

        for (unsigned int j = 0; j < N; ++j) {
            ans.push_back(matrix[j][i]);
        }

        return ans;
    }

    std::vector<std::vector<Field>> getMatrix() const {
        return matrix;
    }

    Matrix<M, N, Field> transposed() const {
        std::vector<std::vector<Field>> ans_matrix;

        for (unsigned int i = 0; i < M; ++i) {
            ans_matrix.push_back(this->getColumn(i));
        }

        Matrix<M, N, Field> tr = ans_matrix;
        return tr;
    }

    Matrix<N, M, Field>& operator+=(const Matrix<N, M, Field>& another) {
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                matrix[i][j] += another.matrix[i][j];
            }
        }

        return *this;
    }

    Matrix<N, M, Field>& operator-=(const Matrix<N, M, Field>& another) {
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                matrix[i][j] -= another.matrix[i][j];
            }
        }

        return *this;
    }

    Matrix<N, M, Field>& operator*=(int x) {
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                matrix[i][j] *= x;
            }
        }

        return *this;
    }

    Matrix<N, M, Field>& operator*=(Field x) {
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                matrix[i][j] *= x;
            }
        }

        return *this;
    }

    Matrix<N, N, Field>& operator*=(const Matrix<N, N, Field>& another) {
        Matrix<N, N, Field> tmp = *this;

        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < N; ++j) {
                Field value = 0;

                for (unsigned int k = 0; k < N; ++k) {
                    value += matrix[i][k] * another[k][j];
                }

                tmp.matrix[i][j] = value;
            }
        }

        *this = tmp;
        return *this;
    }

    unsigned int rank() const {
        std::vector<std::vector<Field>> copy;

        for (size_t i = 0; i < N; ++i) {
            std::vector<Field> line;
            for (size_t j = 0; j < M; ++j) {
                line.push_back(matrix[i][j]);
            }
            copy.push_back(line);
        }

        unsigned int rank = M;

        for (size_t row = 0; row < rank; row++) {
            if (copy[row][row] != static_cast<Field>(0)) {
                for (size_t col = 0; col < N; col++) {
                    if (col != row) {
                        Field mult = copy[col][row] / copy[row][row];

                        for (size_t i = 0; i < rank; i++) {
                            copy[col][i] -= mult * copy[row][i];
                        }
                    }
                }
            } else {
                bool reduce = true;

                for (size_t i = row + 1; i < N;  i++) {
                    if (copy[i][row] != static_cast<Field>(0)) {

                        for (size_t j = 0; j < rank; j++) {
                            Field temp = copy[row][j];
                            copy[row][j] = copy[i][j];
                            copy[i][j] = temp;
                        }

                        reduce = false;
                        break ;
                    }
                }

                if (reduce) {
                    rank--;

                    for (size_t i = 0; i < N; i ++) {
                        copy[i][row] = copy[i][rank];
                    }
                }

                row--;
            }
        }

        return rank;
    }

    Field trace() const {
        static_assert (N == M);

        Field trace = 0;

        for (size_t i = 0; i < N; ++i) {
            trace += matrix[i][i];
        }

        return trace;
    }

    Field det() const {
        static_assert (N == M);

        Field determinant = 1;

        if (N < 1) {
            determinant = 0;
        } else if (N == 1) {
            determinant = matrix[0][0];
        } else if (N == 2) {
            determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        } else {
            std::vector<std::vector<Field>> copy;

            for (size_t i = 0; i < N; ++i) {
                std::vector<Field> line;
                for (size_t j = 0; j < M; ++j) {
                    line.push_back(matrix[i][j]);
                }
                copy.push_back(line);
            }

            for (size_t i = 0; i < N; ++i) {
                if (copy[i][i] != static_cast<Field>(0)) {
                    for (size_t j = i + 1; j < N; ++j) {
                        Field mult = copy[j][i] / copy[i][i];

                        for (size_t k = i; k < N; ++k) {
                            copy[j][k] -= mult * copy[i][k];
                        }
                    }
                    determinant *= copy[i][i];
                } else {
                    determinant = 0;
                    break;
                }
            }
        }

        return determinant;
    }

    Matrix<N, N, Field> inverted() const {
        static_assert (N == M);

        std::vector<std::vector<Field>> copy;

        for (size_t i = 0; i < N; ++i) {
            std::vector<Field> line;
            for (size_t j = 0; j < M; ++j) {
                line.push_back(matrix[i][j]);
            }
            copy.push_back(line);
        }

        Field temp;

        for (size_t i = 0; i < N; i++) {
            for (size_t j = N; j < 2 * N; j++) {
                if (j == i + N) {
                    copy[i].push_back(1);
                } else {
                    copy[i].push_back(0);
                }
            }
        }

        for (size_t i = N - 1; i > 0; i--) {
            if (copy[i - 1][0] < copy[i][0]) {
                std::vector<Field> temp = copy[i];
                copy[i] = copy[i - 1];
                copy[i - 1] = temp;
            }
        }

        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {
                if (j != i) {
                    temp = copy[j][i] / copy[i][i];

                    for (size_t k = 0; k < 2 * N; k++) {
                        copy[j][k] -= copy[i][k] * temp;
                    }
                }
            }
        }

        for (size_t i = 0; i < N; i++) {
            temp = copy[i][i];

            for (size_t j = 0; j < 2 * N; j++) {
                copy[i][j] = copy[i][j] / temp;
            }
        }

        std::vector<std::vector<Field>> inverted;
        for (size_t i = 0; i < N; ++i) {
            std::vector<Field> row;

            for (size_t j = N; j < 2 * N; ++j) {
                row.push_back(copy[i][j]);
            }
            inverted.push_back(row);
        }

        Matrix<N, N, Field> ans = inverted;
        return ans;
    }

    void invert() {
        static_assert (N == M);
        *this = this->inverted();
    }

    template <unsigned int U, unsigned int T>
    bool operator==(const Matrix<U, T, Field>& another) const {
        if (N != U || M != T) {
            return false;
        } else {
            std::vector<std::vector<Field>> another_matrix = another.getMatrix();

            for (size_t i = 0; i < N; ++i) {
                for (size_t j = 0; j < M; ++j) {
                    if (matrix[i][j] != another_matrix[i][j]) {
                        return false;
                    }
                }
            }

            return true;
        }
    }

    template <unsigned int U, unsigned int T>
    bool operator!=(const Matrix<U, T, Field> another) const {
        return !(*this == another);
    }
};

template <unsigned int N, unsigned int M, typename Field>
Matrix<N, M, Field> operator+(const Matrix<N, M, Field>& matrix_1, const Matrix<N, M, Field>& matrix_2) {
    Matrix<N, M, Field> copy = matrix_1;
    copy += matrix_2;
    return copy;
}

template <unsigned int N, unsigned int M, typename Field>
Matrix<N, M, Field> operator-(const Matrix<N, M, Field>& matrix_1, const Matrix<N, M, Field>& matrix_2) {
    Matrix<N, M, Field> copy = matrix_1;
    copy -= matrix_2;
    return copy;
}

template <unsigned int N, unsigned int M, typename Field>
Matrix<N, M, Field> operator*(const Matrix<N, M, Field>& matrix, int x) {
    Matrix<N, M, Field> copy = matrix;
    copy *= x;
    return copy;
}

template <unsigned int N, unsigned int M, typename Field>
Matrix<N, M, Field> operator*(int x, const Matrix<N, M, Field>& matrix) {
    Matrix<N, M, Field> copy = matrix;
    copy *= x;
    return copy;
}

template <unsigned int N, unsigned int M, typename Field>
Matrix<N, M, Field> operator*(const Matrix<N, M, Field>& matrix, Field x) {
    Matrix<N, M, Field> copy = matrix;
    copy *= x;
    return copy;
}

template <unsigned int N, unsigned int M, typename Field>
Matrix<N, M, Field> operator*(Field x, const Matrix<N, M, Field>& matrix) {
    Matrix<N, M, Field> copy = matrix;
    copy *= x;
    return copy;
}

template <unsigned int N, unsigned int M, unsigned int L, typename Field>
Matrix<N, L, Field> operator*(const Matrix<N, M, Field>& matrix_1, const Matrix<M, L, Field>& matrix_2) {
    std::vector<std::vector<Field>> ans_matrix;

    for (unsigned int i = 0; i < N; ++i) {
        std::vector<Field> row;
        for (unsigned int j = 0; j < L; ++j) {
            Field value = 0;

            for (unsigned int k = 0; k < M; ++k) {
                value += matrix_1[i][k] * matrix_2[k][j];
            }
            row.push_back(value);
        }
        ans_matrix.push_back(row);
    }

    return Matrix<N, L, Field> (ans_matrix);
}

template <unsigned int N, typename Field = Rational>
class SquareMatrix: public Matrix<N, N, Field> {
public:
    SquareMatrix<N, Field>(): Matrix<N, N, Field> () {}

    SquareMatrix<N, Field>(std::vector<std::vector<Field>> src_matrix): Matrix<N, N, Field>(src_matrix) {}

    SquareMatrix<N, Field>(std::vector<std::vector<int>> src_matrix): Matrix<N, N, Field>(src_matrix) {}
};
