#include "io_utils.h"
#include <stdio.h>

void save_solution_to_file(const char *filepath, const ParsedSystem *ps, const LinearSystemSolution *sol){
    FILE *f=fopen(filepath, "a");
    if(!f) return;
    fprintf(f, "Sistema com variaveis: %s\n", ps->vars?ps->vars:"");
    fprintf(f, "A (%zux%zu) e b:\n", ps->A.rows, ps->A.cols);
    for(size_t i=0;i<ps->A.rows;i++){
      for(size_t j=0;j<ps->A.cols;j++) fprintf(f, "%8.4f ", ps->A.data[i*ps->A.cols+j]);
      fprintf(f, " | %8.4f\n", ps->b.data[i]);
    }
    fprintf(f, "rank(A)=%zu, rank([A|b])=%zu\n", sol->rankA, sol->rankAug);
    if(sol->type==SYS_INCONSISTENT){
        fprintf(f, ">> Sistema inconsistente (sem solucao).\n\n");
    }else if(sol->type==SYS_UNIQUE){
        fprintf(f, ">> Solucao unica:\n");
        for(size_t i=0;i<sol->x.n;i++){
            fprintf(f, "  %c = %.6f\n", ps->vars?ps->vars[i]:'x', sol->x.data[i]);
        }
        fprintf(f, "\n");
    }else{
        fprintf(f, ">> Infinas solucoes: x = x_particular + span(N)\n");
        fprintf(f, "x_particular:\n");
        for(size_t i=0;i<sol->x.n;i++){
            fprintf(f, "  %c = %.6f\n", ps->vars?ps->vars[i]:'x', sol->x.data[i]);
        }
        fprintf(f, "Base do nucleo (colunas): %zux%zu\n", sol->nullspace_basis.rows, sol->nullspace_basis.cols);
        for(size_t i=0;i<sol->nullspace_basis.rows;i++){
            for(size_t j=0;j<sol->nullspace_basis.cols;j++){
                fprintf(f, "%8.4f ", sol->nullspace_basis.data[i*sol->nullspace_basis.cols + j]);
            }
            fprintf(f, "\n");
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void save_matrix_to_file(const char *filepath, const Matrix *A, const char *name){
    FILE *f=fopen(filepath, "a");
    if(!f) return;
    fprintf(f, "%s (%zux%zu):\n", name?name:"Matriz", A->rows, A->cols);
    for(size_t i=0;i<A->rows;i++){
        for(size_t j=0;j<A->cols;j++)
            fprintf(f, "%10.6f ", A->data[i*A->cols+j]);
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
    fclose(f);
}

void save_vectors_to_file(const char *filepath, const Matrix *V, const char *name){
    save_matrix_to_file(filepath, V, name);
}
