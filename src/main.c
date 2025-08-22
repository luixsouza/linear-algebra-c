#include <stdio.h>
#include <stdlib.h>
#include "linalg.h"
#include "parser.h"
#include "io_utils.h"

#define EPS 1e-9

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

static void solve_from_console() {
    Matrix A = read_matrix();
    Vector b = read_vector_input(A.rows);
    LinearSystemSolution sol = solve_linear_system(&A, &b, EPS);
    
    printf("Tipo da Solucao: ");
    if (sol.type == SYS_UNIQUE) {
        printf("Unica\n");
        vec_print(&sol.x, "x");
    } else if (sol.type == SYS_INFINITE) {
        printf("Infinita\n");
        vec_print(&sol.x, "Solucao Particular");
        mat_print(&sol.nullspace_basis, "Base do Espaco Nulo");
    } else {
        printf("Inconsistente\n");
    }

    mat_free(&A);
    vec_free(&b);
    vec_free(&sol.x);
    mat_free(&sol.nullspace_basis);
}

static void check_transform() {
    Matrix A = read_matrix();
    int inj = is_injective(&A, EPS);
    int surj = is_surjective(&A, EPS);
    if (inj && surj) printf("Bijetiva\n");
    else if (inj) printf("Injetiva\n");
    else if (surj) printf("Sobrejetiva\n");
    else printf("Nenhuma\n");
    mat_free(&A);
}

static void check_basis() {
    Matrix A = read_matrix();
    printf("%s\n", forms_basis(&A, EPS) ? "Forma base" : "Nao forma base");
    mat_free(&A);
}

static void eigen() {
    Matrix A = read_matrix();
    Vector vals = eigenvalues_qr(&A, 1000, EPS);
    
    printf("Autovalores: ");
    for (size_t i = 0; i < vals.n; i++) printf("%.6lf ", vals.data[i]);
    printf("\n");

    mat_free(&A);
    vec_free(&vals);
}

static void diag() {
    Matrix A = read_matrix();
    Matrix P, D;
    if (diagonalize(&A, &P, &D, EPS, 1000)) {
        printf("Matriz P (Autovetores):\n");
        mat_print(&P, NULL);
        printf("\nMatriz D (Autovalores):\n");
        mat_print(&D, NULL);
        mat_free(&P);
        mat_free(&D);
    } else {
        printf("A matriz nao pode ser diagonalizada.\n");
    }
    mat_free(&A);
}

static void solve_from_file() {
    char filepath[256];
    printf("Digite o caminho do arquivo: ");
    scanf("%s", filepath);

    ParsedSystem ps = parse_from_file(filepath, NULL);
    if (ps.vars == NULL) {
        printf("Erro ao ler o arquivo ou arquivo vazio.\n");
        return;
    }

    LinearSystemSolution sol = solve_linear_system(&ps.A, &ps.b, EPS);
    
    printf("Tipo da Solucao: ");
    if (sol.type == SYS_UNIQUE) {
        printf("Unica\n");
        vec_print(&sol.x, "x");
    } else if (sol.type == SYS_INFINITE) {
        printf("Infinita\n");
        vec_print(&sol.x, "Solucao Particular");
        mat_print(&sol.nullspace_basis, "Base do Espaco Nulo");
    } else {
        printf("Inconsistente\n");
    }

    free_parsed_system(&ps);
    vec_free(&sol.x);
    mat_free(&sol.nullspace_basis);
}

static void menu() {
    printf("\n==== Algebra Linear ====\n");
    printf("1) Resolver sistema linear (Console)\n");
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
        int c; while ((c = getchar()) != '\n' && c != EOF) {}
        switch (op) {
            case 1: solve_from_console(); break;
            case 2: check_transform(); break;
            case 3: check_basis(); break;
            case 4: eigen(); break;
            case 5: diag(); break;
            case 6: solve_from_file(); break;
            case 0: return 0;
            default: printf("Opcao invalida\n"); break;
        }
    }
    return 0;
}