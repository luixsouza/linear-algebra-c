#ifndef LINALG_H
#define LINALG_H

#include <stddef.h>

typedef struct {
    size_t rows, cols;
    double *data;
} Matrix;

typedef struct {
    size_t n;
    double *data;
} Vector;

typedef enum {
    SYS_UNIQUE=0,
    SYS_INFINITE=1,
    SYS_INCONSISTENT=2
} SystemType;

typedef struct {
    SystemType type;
    Vector x;
    Matrix nullspace_basis;
    size_t rankA, rankAug;
} LinearSystemSolution;

Matrix mat_create(size_t r, size_t c);
Matrix mat_identity(size_t n);
void   mat_free(Matrix *A);
double* mat_at(Matrix *A, size_t i, size_t j);
double  mat_get(const Matrix *A, size_t i, size_t j);
void    mat_set(Matrix *A, size_t i, size_t j, double v);
Matrix  mat_copy(const Matrix *A);
void    mat_print(const Matrix *A, const char *name);
Vector  vec_create(size_t n);
void    vec_free(Vector *v);
void    vec_print(const Vector *v, const char *name);

Matrix mat_mul(const Matrix *A, const Matrix *B);
Matrix mat_sub(const Matrix *A, const Matrix *B);
Matrix mat_add(const Matrix *A, const Matrix *B);
Matrix mat_scale(const Matrix *A, double s);
int    mat_equal_size(const Matrix *A, const Matrix *B);
Matrix mat_transpose(const Matrix *A);

size_t mat_rank_rref(Matrix *A, double eps);
size_t mat_rank(const Matrix *A, double eps);
int    mat_rref(Matrix *A, double eps);
double mat_det_2x2(const Matrix *A);
double mat_det_3x3(const Matrix *A);

LinearSystemSolution solve_linear_system(const Matrix *A, const Vector *b, double eps);

int is_injective(const Matrix *A, double eps);
int is_surjective(const Matrix *A, double eps);
int is_bijective(const Matrix *A, double eps);

int forms_basis(const Matrix *V, double eps);

Vector eigenvalues_qr(const Matrix *A, size_t max_iters, double eps);
Matrix eigenvectors_for_lambda(const Matrix *A, double lambda, double eps);

int diagonalize(const Matrix *A, Matrix *P_out, Matrix *D_out, double eps, size_t max_iters);

#endif