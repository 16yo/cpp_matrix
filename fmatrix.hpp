#ifndef FMATRIX_H

#define FMATRIX_H

#include <vector>
#include <cinttypes>
#include <iostream>
#include "func.hpp"
#include "matrix.hpp"


#define mapping func<out, in>

template <class out, class in>
class FMatrix : public Matrix<mapping> {
public: 
    FMatrix() : Matrix<mapping>() {}

    FMatrix(const uint64_t n, const uint64_t m, const mapping& f = mapping{}) : \
     Matrix<mapping>(n, m, f) {}

    FMatrix(const uint64_t n, const mapping& f = mapping{}) : \
     Matrix<mapping>(n, f) {}

    FMatrix(const std::vector<std::vector<mapping>>& mx) : Matrix<mapping>(mx) {}

    FMatrix (const FMatrix& M_) : Matrix<mapping>(M_) {}

    FMatrix (const Matrix<in>& M_) : FMatrix(cast_from_const_to_mapping(M_)) {}

    FMatrix (const Matrix<mapping>& M_) {
        this->n = M_.rows();
        this->m = M_.colomns();
        this->M = M_.to_vector();
    }

    Matrix<out> operator()(in value = in{}) {
        Matrix<out> R(this->rows(), this->colomns());
        for (uint64_t i = 0; i < this->rows(); ++i)
            for (uint64_t j = 0; j < this->colomns(); ++j)
                R[i][j] = (*this)[i][j](value);
        return R;
    }


private:
    std::vector<std::vector<mapping>> cast_from_const_to_mapping(const Matrix<in>& M_) {
        std::vector<std::vector<mapping>> mx(M_.rows(), std::vector<mapping>(M_.colomns()));
        for (uint64_t i = 0; i < M_.rows(); ++i) {
            for (uint64_t j = 0; j < M_.colomns(); ++j) {
                mx[i][j] = mapping(M_[i][j]);
            }
        }
        return mx;
    }

};


#endif