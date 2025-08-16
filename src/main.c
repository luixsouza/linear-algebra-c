#include <stdio.h>
#include <stdlib.h>
#include "linalg.h"
#include "parser.h"
#include "io_utils.h"

// === Leitura simplificada de matriz e vetor ===
static Matrix read_matrix() {
    size_t m, n;
    printf("Linhas Colunas: ");
    scanf("%zu %zu", &m, &n);
    Matrix A = mat_create(m, n);
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n; j++)
            scanf("%lf", mat_at(&A, i, j));
    return A;
}

static Vector read_vector_input(size_t n) {
    Vector v = vec_create(n);
    for (size_t i = 0; i < n; i++)
        scanf("%lf", &v.data[i]);
    return v;
}

// === Funções minimalistas ===
static void solve_from_console() {
    Matrix A = read_matrix();
    Vector b = read_vector_input(A.rows);
    Vector x = gauss_jordan(A, b);
    for (size_t i = 0; i < x.size; i++) printf("%.6lf ", x.data[i]);
    printf("\n");
    mat_free(&A); vec_free(&b); vec_free(&x);
}

static void check_transform() {
    Matrix A = read_matrix();
    int inj = is_injective(A);
    int surj = is_surjective(A);
    if (inj && surj) printf("Bijetiva\n");
    else if (inj) printf("Injetiva\n");
    else if (surj) printf("Sobrejetiva\n");
    else printf("Nenhuma\n");
    mat_free(&A);
}

static void check_basis() {
    Matrix A = read_matrix();
    printf("%s\n", is_basis(A) ? "Forma base" : "Não forma base");
    mat_free(&A);
}

static void eigen() {
    Matrix A = read_matrix();
    Vector vals = eigenvalues(A);
    for (size_t i = 0; i < vals.size; i++) printf("%.6lf ", vals.data[i]);
    printf("\n");
    mat_free(&A); vec_free(&vals);
}

static void diag() {
    Matrix A = read_matrix();
    Matrix D = diagonalize(A);
    for (size_t i = 0; i < D.rows; i++) {
        for (size_t j = 0; j < D.cols; j++)
            printf("%.6lf ", mat_at(&D, i, j)[0]);
        printf("\n");
    }
    mat_free(&A); mat_free(&D);
}

static void solve_from_file() {
    Matrix A; Vector b;
    read_system_from_file(&A, &b);
    Vector x = gauss_jordan(A, b);
    for (size_t i = 0; i < x.size; i++) printf("%.6lf ", x.data[i]);
    printf("\n");
    mat_free(&A); vec_free(&b); vec_free(&x);
}

// === Menu visual ===
static void menu() {
    printf("\n==== Álgebra Linear em C ====\n");
    printf("1) Resolver sistema linear\n");
    printf("2) Verificar injetividade/sobrejetividade/bijetividade\n");
    printf("3) Verificar se um conjunto forma base\n");
    printf("4) Autovalores\n");
    printf("5) Diagonalizar matriz\n");
    printf("6) Resolver sistema de arquivo\n");
    printf("0) Sair\n");
    printf("Escolha: ");
}

int main() {
    int op;
    while (1) {
        menu();
        if (scanf("%d", &op) != 1) break;
        int c; while ((c = getchar()) != '\n' && c != EOF) {} // limpa buffer
        switch (op) {
            case 1: solve_from_console(); break;
            case 2: check_transform(); break;
            case 3: check_basis(); break;
            case 4: eigen(); break;
            case 5: diag(); break;
            case 6: solve_from_file(); break;
            case 0: return 0;
            default: printf("Opção inválida\n"); break;
        }
    }
    return 0;
}
