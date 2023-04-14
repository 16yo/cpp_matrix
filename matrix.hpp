#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <cinttypes>
#include <iostream>

template <typename T>
class Matrix;

template <typename T>
class Matrix {
public:

    Matrix() {
        n = m = 1;
        M = std::vector<std::vector<T>>(1, std::vector<T>(1, T{}));
    }

    Matrix(const uint64_t n, const uint64_t m, const T& fill = T{}) {
        this->n = n;
        this->m = m;
        this->M = std::vector<std::vector<T>>(n, std::vector<T>(m, fill));
    }

    //Matrix(const uint64_t n, const T& fill = T{}) : Matrix(n, n, fill) {}

    Matrix(const std::vector<std::vector<T>>& mx) {
        n = mx.size();
        m = mx.end() - mx.begin();
        if (n == 0 || m == 0) {
            Matrix();
            return;
        }
        M = mx;
    }

    Matrix(const Matrix& M_) : Matrix(M_.M) {}

    static Matrix O(const uint64_t n) {
        return Matrix(n);
    }

    static Matrix E(const uint64_t n) {
        Matrix e(n);
        for (uint64_t i = 0; i < n; ++i) {
            e[i][i] = 1;
        }
        return e;
    }

    static T algebraic_addition(const Matrix& mx, const uint64_t r, uint64_t c) {
        return T(((r + c) % 2 == 0 ? 1 : -1)) * det(mx.erase_row_n_colomn(r, c));
    }

    static Matrix addition_matrix(const Matrix& mx) {
        Matrix a(mx.n, mx.m);
        for (uint64_t i = 0; i < mx.n; ++i) {
            for (uint64_t j = 0; j < mx.m; ++j) {
                a[i][j] = algebraic_addition(mx, i, j);
            }
        }
        return a.transposition();
    }

    static T det(const Matrix& mx) {
        if (!is_square(mx)) {
            throw std::invalid_argument("Matrix is not square");
        }
        if (mx.n == 1) {
            return mx[0][0];
        }
        if (mx.n == 2) {
            return mx[0][0] * mx[1][1] - mx[0][1] * mx[1][0];
        }
        T d{};
        for (uint64_t i = 0; i < mx.n; ++i) {
            d += (i % 2 == 0 ? T(1) : T(-1)) * mx[0][i] * det(mx.erase_row_n_colomn(0, i));
        }
        return d;
    }

    Matrix& operator=(const Matrix& M_) {
        this->n = M_.n;
        this->m = M_.m;
        this->M = M_.M;
    }    

    const std::vector<T>& operator[](const uint64_t n) const {
        return this->M[n];
    }

    std::vector<T>& operator[](const uint64_t n) {
        return this->M[n];
    }

    Matrix& operator+=(const Matrix& M_) {
       if (! (this->n == M_.n && this->m == M_.n) ) {
           throw std::invalid_argument("Matrixs dimentions do not match");
       }
       for (uint64_t i = 0; i < this->n; ++i) {
           for (uint64_t j = 0; j < this->m; ++j) {
               this->M[i][j] += M_[i][j];
           }
       }
       return *this;
    }

    Matrix& operator-=(const Matrix& M_) {
        if (! (this->n == M_.n && this->m == M_.n) ) {
           throw std::invalid_argument("Matrixs dimentions do not match");
       }
       for (uint64_t i = 0; i < this->n; ++i) {
           for (uint64_t j = 0; j < this->m; ++j) {
               this->M[i][j] -= M_[i][j];
           }
       }
       return *this;
    }

    Matrix& operator*=(const T k) {
        for (auto& i : M)
            for (auto& j : i)
                j *= k;
        return *this;
    }

    Matrix erase_row(const uint64_t r) const {
        Matrix a = *this;
        a.M.erase(std::next(a.M.begin(), r));
        --a.n;
        return a;
    }

    Matrix erase_colomn(const uint64_t c) const {
        Matrix a(n, m - 1);
        for (uint64_t i = 0; i < n; ++i) {
            for (uint64_t j = 0; j < m; ++j) {
                if (j == c) 
                    continue;
                a[i][j - (j < c ? 0 : 1)] = M[i][j];
            }
        }
        return a;
    }

    Matrix erase_row_n_colomn(const uint64_t r, const uint64_t c) const {
        Matrix a(n - 1, m - 1);
        uint8_t rs = 0, cs = 0;
        for (uint64_t i = 0; i < n; ++i) {
            cs = 0;
            if (i == r) {
                rs = 1;
                continue;
            }
            for (uint64_t j = 0; j < m; ++j) {
                if (j == c) {
                    cs = 1;
                    continue;
                }
                a[i - rs][j - cs] = M[i][j];
            }
        }
        return a;
    }

    Matrix operator*(const Matrix& M_) const {
        if (! (this->m == M_.n)) {
            throw std::invalid_argument("Matrixs dimentions do not match");
        }
        Matrix t(this->n, M_.m);
        for (uint64_t i = 0; i < this->n; ++i) {
            for (uint64_t j = 0; j < M_.m; ++j) {
                for (uint64_t k = 0; k < this->m; ++k) {
                    t[i][j] += this->M[i][k] * M_[k][j];
                }
            }
        }
        return t;
    }

    Matrix& operator*=(const Matrix& M_) {
        return *this = *this * M_;
    }

    Matrix operator+(const Matrix& M_) const {
       Matrix M__ = *this;
       return M__ += M_;
    }

    Matrix operator-(const Matrix& M_) const {
       Matrix M__ = *this;
       return M__ -= M_;
    }

    Matrix transposition() {
        Matrix t(n, m);
        for (uint64_t i = 0; i < n; ++i) {
            for (uint64_t j = 0; j < m; ++j) {
                t[i][j] = (*this)[j][i];
            }
        }
        return t;
    } 

    Matrix reverse() const {
        T d = det(*this);
        if (d == T{}) {
            throw std::invalid_argument("Matrix is singular");
        }
        return addition_matrix(*this) * (T(1) / d);
    }

    Matrix pow(const int64_t k) const {
        if (!is_square(*this)) {
            throw std::invalid_argument("Matrix is not square");
        }
        if (k == -1) {
            return this->reverse();
        }
        if (k == 0) {
            return E(this->n);
        }
        if (k == 1) {
            return *this;
        }
        Matrix k_2 = this->pow(k/2);
        if (k % 2 == 0) {
            return k_2 * k_2;
        }
        return k_2*k_2*(*this);
    }

    static bool is_square(const Matrix& mx) {
        return mx.n == mx.m;
    }

    uint64_t rows() const {
        return this->n;
    }

    uint64_t colomns() const {
        return this->m;
    }

    std::vector<std::vector<T>> to_vector() const {
        return M;
    }

protected:
    uint64_t n;
    uint64_t m;
    std::vector<std::vector<T>> M;

};

template<typename T>
inline Matrix<T> operator * (const T& k, const Matrix<T>& mx) {
    Matrix<T> a = mx;
    return a *= k;
}

template<typename T>
inline Matrix<T> operator * (const Matrix<T>& mx, const T& k) {
    Matrix<T> a = mx;
    return a *= k;
}

template <typename T>
std::ostream& operator << (std::ostream& out, const Matrix<T>& M) {
    for (uint64_t i = 0; i < M.rows(); ++i) {
        for (int j = 0; j < M.colomns(); ++j)
            out << M[i][j] << " ";
        out << std::endl;
    }
    return out;
}

template <typename T>
std::istream& operator >> (std::istream& in, Matrix<T>& M) {
    for (uint64_t i = 0; i < M.rows(); ++i)
        for (int j = 0; j < M.colomns(); ++j)
            in >> M[i][j];
    return in;
}

#endif