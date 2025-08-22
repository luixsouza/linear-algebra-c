#include "linalg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double dabs(double x){ return x<0?-x:x; }
static int nearly_zero(double x, double eps){ return dabs(x) < eps; }

Matrix mat_create(size_t r, size_t c){
    Matrix A; A.rows=r; A.cols=c; A.data=(double*)calloc(r*c, sizeof(double)); return A;
}
Matrix mat_identity(size_t n){
    Matrix I=mat_create(n,n);
    for(size_t i=0;i<n;i++) *mat_at(&I,i,i)=1.0;
    return I;
}
void mat_free(Matrix *A){
    if(A && A->data){ free(A->data); A->data=NULL; A->rows=A->cols=0; }
}
double* mat_at(Matrix *A, size_t i, size_t j){ return &A->data[i*A->cols + j]; }
double  mat_get(const Matrix *A, size_t i, size_t j){ return A->data[i*A->cols + j]; }
void    mat_set(Matrix *A, size_t i, size_t j, double v){ A->data[i*A->cols + j]=v; }

Matrix  mat_copy(const Matrix *A){
    Matrix B=mat_create(A->rows,A->cols);
    memcpy(B.data, A->data, sizeof(double)*A->rows*A->cols);
    return B;
}

void mat_print(const Matrix *A, const char *name){
    if(name) printf("%s =\n", name);
    for(size_t i=0;i<A->rows;i++){
        for(size_t j=0;j<A->cols;j++){
            printf("%10.6f ", mat_get(A,i,j));
        }
        printf("\n");
    }
}

Vector vec_create(size_t n){
    Vector v; v.n=n; v.data=(double*)calloc(n,sizeof(double)); return v;
}
void vec_free(Vector *v){
    if(v && v->data){ free(v->data); v->data=NULL; v->n=0; }
}
void vec_print(const Vector *v, const char *name){
    if(name) printf("%s^T = ", name);
    printf("[");
    for(size_t i=0;i<v->n;i++){
        printf("%10.6f", v->data[i]);
        if(i+1<v->n) printf(", ");
    }
    printf("]\n");
}

Matrix mat_mul(const Matrix *A, const Matrix *B){
    Matrix C=mat_create(A->rows, B->cols);
    for(size_t i=0;i<A->rows;i++){
        for(size_t k=0;k<A->cols;k++){
            double aik = mat_get(A,i,k);
            for(size_t j=0;j<B->cols;j++){
                C.data[i*C.cols+j] += aik * mat_get(B,k,j);
            }
        }
    }
    return C;
}
Matrix mat_sub(const Matrix *A, const Matrix *B){
    Matrix C=mat_create(A->rows, A->cols);
    for(size_t i=0;i<A->rows*A->cols;i++) C.data[i]=A->data[i]-B->data[i];
    return C;
}
Matrix mat_add(const Matrix *A, const Matrix *B){
    Matrix C=mat_create(A->rows, A->cols);
    for(size_t i=0;i<A->rows*A->cols;i++) C.data[i]=A->data[i]+B->data[i];
    return C;
}
Matrix mat_scale(const Matrix *A, double s){
    Matrix C=mat_create(A->rows, A->cols);
    for(size_t i=0;i<A->rows*A->cols;i++) C.data[i]=A->data[i]*s;
    return C;
}
int mat_equal_size(const Matrix *A, const Matrix *B){
    return A->rows==B->rows && A->cols==B->cols;
}
Matrix mat_transpose(const Matrix *A){
    Matrix T=mat_create(A->cols, A->rows);
    for(size_t i=0;i<A->rows;i++)
        for(size_t j=0;j<A->cols;j++)
            mat_set(&T, j, i, mat_get(A,i,j));
    return T;
}

int mat_rref(Matrix *A, double eps){
    size_t lead=0;
    for(size_t r=0; r<A->rows; r++){
        if(lead >= A->cols) return 0;
        size_t i=r;
        while(nearly_zero(mat_get(A,i,lead), eps)){
            i++;
            if(i==A->rows){
                i=r; lead++;
                if(lead==A->cols) return 0;
            }
        }
        if(i!=r){
            for(size_t j=0;j<A->cols;j++){
                double tmp=mat_get(A,r,j);
                mat_set(A,r,j, mat_get(A,i,j));
                mat_set(A,i,j, tmp);
            }
        }
        double lv = mat_get(A,r,lead);
        if(!nearly_zero(lv,eps)){
            for(size_t j=0;j<A->cols;j++) mat_set(A,r,j, mat_get(A,r,j)/lv);
        }
        for(size_t i2=0;i2<A->rows;i2++){
            if(i2!=r){
                double lv2 = mat_get(A,i2,lead);
                if(!nearly_zero(lv2,eps)){
                    for(size_t j=0;j<A->cols;j++)
                        mat_set(A,i2,j, mat_get(A,i2,j) - lv2*mat_get(A,r,j));
                }
            }
        }
        lead++;
    }
    return 0;
}

size_t mat_rank_rref(Matrix *A, double eps){
    mat_rref(A, eps);
    size_t rank=0;
    for(size_t i=0;i<A->rows;i++){
        int nonzero=0;
        for(size_t j=0;j<A->cols;j++){
            if(!nearly_zero(mat_get(A,i,j), eps)){ nonzero=1; break; }
        }
        if(nonzero) rank++;
    }
    return rank;
}

size_t mat_rank(const Matrix *A, double eps){
    Matrix B=mat_copy(A);
    size_t r = mat_rank_rref(&B, eps);
    mat_free(&B);
    return r;
}

double mat_det_2x2(const Matrix *A){
    if(A->rows!=2 || A->cols!=2) return 0.0;
    return mat_get(A,0,0)*mat_get(A,1,1) - mat_get(A,0,1)*mat_get(A,1,0);
}
double mat_det_3x3(const Matrix *A){
    if(A->rows!=3 || A->cols!=3) return 0.0;
    double a=mat_get(A,0,0), b=mat_get(A,0,1), c=mat_get(A,0,2);
    double d=mat_get(A,1,0), e=mat_get(A,1,1), f=mat_get(A,1,2);
    double g=mat_get(A,2,0), h=mat_get(A,2,1), i=mat_get(A,2,2);
    return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
}

LinearSystemSolution solve_linear_system(const Matrix *Ain, const Vector *b, double eps){
    LinearSystemSolution S;
    S.type=SYS_UNIQUE;
    S.rankA=S.rankAug=0;
    S.x=vec_create(Ain->cols);
    S.nullspace_basis=mat_create(Ain->cols, 0);

    Matrix Ab = mat_create(Ain->rows, Ain->cols+1);
    for(size_t i=0;i<Ain->rows;i++){
        for(size_t j=0;j<Ain->cols;j++) mat_set(&Ab,i,j, mat_get(Ain,i,j));
        mat_set(&Ab,i,Ain->cols, b->data[i]);
    }
    Matrix A = mat_copy(Ain);
    S.rankA = mat_rank(&A, eps);
    mat_free(&A);
    Matrix Ab_copy = mat_copy(&Ab);
    S.rankAug = mat_rank(&Ab_copy, eps);
    mat_free(&Ab_copy);

    mat_rref(&Ab, eps);

    size_t m=Ab.rows, n=Ain->cols;
    int inconsistent=0;
    for(size_t i=0;i<m;i++){
        int allzero=1;
        for(size_t j=0;j<n;j++){
            if(!nearly_zero(mat_get(&Ab,i,j), eps)){ allzero=0; break; }
        }
        if(allzero && !nearly_zero(mat_get(&Ab,i,n), eps)){ inconsistent=1; break; }
    }
    if(inconsistent){
        S.type=SYS_INCONSISTENT;
        mat_free(&Ab);
        return S;
    }

    int *pivot_col = (int*)malloc(n*sizeof(int));
    for(size_t j=0;j<n;j++) pivot_col[j]=-1;
    size_t row=0;
    for(size_t col=0; col<n && row<m; col++){
        size_t rlead = (size_t)-1;
        for(size_t i=row;i<m;i++){
            if(nearly_zero(mat_get(&Ab,i,col)-1.0, 1e-7)){
                rlead=i; break;
            }
        }
        if(rlead!=(size_t)-1){
            pivot_col[col]=(int)rlead;
            row=rlead+1;
        }
    }
    for(size_t j=0;j<n;j++) S.x.data[j]=0.0;
    for(size_t j=0;j<n;j++){
        if(pivot_col[j]>=0){
            S.x.data[j]=mat_get(&Ab, pivot_col[j], n);
        }
    }

    size_t free_count=0;
    for(size_t j=0;j<n;j++) if(pivot_col[j]<0) free_count++;
    if(free_count>0){
        Matrix N=mat_create(n, free_count);
        size_t k=0;
        for(size_t j=0;j<n;j++){
            if(pivot_col[j]<0){
                for(size_t t=0;t<n;t++) mat_set(&N,t,k, 0.0);
                mat_set(&N,j,k, 1.0);
                for(size_t pcol=0;pcol<n;pcol++){
                    if(pivot_col[pcol]>=0){
                        double coeff = -mat_get(&Ab, pivot_col[pcol], j);
                        mat_set(&N, pcol, k, coeff);
                    }
                }
                k++;
            }
        }
        S.nullspace_basis=N;
        S.type=SYS_INFINITE;
    }else{
        S.type=SYS_UNIQUE;
    }

    free(pivot_col);
    mat_free(&Ab);
    return S;
}

int is_injective(const Matrix *A, double eps){
    size_t r = mat_rank(A, eps);
    return r == A->cols;
}
int is_surjective(const Matrix *A, double eps){
    size_t r = mat_rank(A, eps);
    return r == A->rows;
}
int is_bijective(const Matrix *A, double eps){
    if(A->rows != A->cols) return 0;
    size_t r = mat_rank(A, eps);
    return r == A->rows;
}

int forms_basis(const Matrix *V, double eps){
    if(V->cols != V->rows) return 0;
    size_t r = mat_rank(V, eps);
    return r == V->rows;
}

static void qr_decompose(const Matrix *A, Matrix *Q, Matrix *R, double eps){
    size_t m=A->rows, n=A->cols;
    *Q = mat_create(m,n);
    *R = mat_create(n,n);
    for(size_t j=0;j<n;j++){
        for(size_t i=0;i<m;i++) mat_set(Q,i,j, mat_get(A,i,j));
        for(size_t k=0;k<j;k++){
            double r=0.0;
            for(size_t i=0;i<m;i++) r += mat_get(Q,i,k)*mat_get(A,i,j);
            mat_set(R,k,j,r);
            for(size_t i=0;i<m;i++)
                mat_set(Q,i,j, mat_get(Q,i,j) - r*mat_get(Q,i,k));
        }
        double norm=0.0;
        for(size_t i=0;i<m;i++){ double v=mat_get(Q,i,j); norm+=v*v; }
        norm = sqrt(norm);
        if(nearly_zero(norm, eps)){
            for(size_t i=0;i<m;i++) mat_set(Q,i,j, 0.0);
            mat_set(R,j,j, 0.0);
        }else{
            for(size_t i=0;i<m;i++) mat_set(Q,i,j, mat_get(Q,i,j)/norm);
            mat_set(R,j,j, norm);
        }
    }
}

Vector eigenvalues_qr(const Matrix *Ain, size_t max_iters, double eps){
    Matrix A = mat_copy(Ain);
    for(size_t it=0; it<max_iters; it++){
        Matrix Q,R;
        qr_decompose(&A, &Q, &R, eps);
        Matrix AQ = mat_mul(&R, &Q);
        mat_free(&A);
        A = AQ;
        mat_free(&Q); mat_free(&R);
    }
    Vector evals = vec_create(A.rows);
    for(size_t i=0;i<A.rows;i++) evals.data[i]=mat_get(&A,i,i);
    mat_free(&A);
    return evals;
}

Matrix eigenvectors_for_lambda(const Matrix *A, double lambda, double eps){
    Matrix M = mat_copy(A);
    for(size_t i=0;i<M.rows;i++) mat_set(&M,i,i, mat_get(&M,i,i)-lambda);
    Matrix R = mat_copy(&M);
    mat_rref(&R, eps);
    size_t m=R.rows, n=R.cols;
    int *pivot = (int*)malloc(n*sizeof(int));
    for(size_t j=0;j<n;j++) pivot[j]=-1;
    size_t row=0;
    for(size_t col=0; col<n && row<m; col++){
        size_t rlead=(size_t)-1;
        for(size_t i=row;i<m;i++){
            if(nearly_zero(mat_get(&R,i,col)-1.0, 1e-7)){ rlead=i; break; }
        }
        if(rlead!=(size_t)-1){ pivot[col]=(int)rlead; row=rlead+1; }
    }
    size_t free_count=0; for(size_t j=0;j<n;j++) if(pivot[j]<0) free_count++;
    Matrix N = mat_create(n, free_count>0?free_count:1);
    if(free_count==0){
        for(size_t i=0;i<n;i++) mat_set(&N,i,0, 0.0);
    }else{
        size_t k=0;
        for(size_t j=0;j<n;j++){
            if(pivot[j]<0){
                for(size_t t=0;t<n;t++) mat_set(&N,t,k, 0.0);
                mat_set(&N,j,k, 1.0);
                for(size_t p=0;p<n;p++){
                    if(pivot[p]>=0){
                        double coeff = -mat_get(&R, pivot[p], j);
                        mat_set(&N, p, k, coeff);
                    }
                }
                k++;
            }
        }
    }
    free(pivot);
    mat_free(&M); mat_free(&R);
    return N;
}

int diagonalize(const Matrix *A, Matrix *P_out, Matrix *D_out, double eps, size_t max_iters){
    if(A->rows != A->cols) return 0;
    size_t n=A->rows;
    Vector evals = eigenvalues_qr(A, max_iters, eps);
    Matrix P = mat_create(n,n);
    for(size_t j=0;j<n;j++){
        Matrix V = eigenvectors_for_lambda(A, evals.data[j], eps);
        for(size_t i=0;i<n;i++) mat_set(&P,i,j, mat_get(&V,i,0));
        mat_free(&V);
    }
    size_t r = mat_rank(&P, eps);
    if(r<n){
        vec_free(&evals);
        mat_free(&P);
        return 0;
    }
    Matrix D = mat_create(n,n);
    for(size_t i=0;i<n;i++) mat_set(&D,i,i, evals.data[i]);
    vec_free(&evals);
    *P_out=P; *D_out=D;
    return 1;
}