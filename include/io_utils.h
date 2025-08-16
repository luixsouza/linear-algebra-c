#ifndef IO_UTILS_H
#define IO_UTILS_H

#include "linalg.h"
#include "parser.h"

void save_solution_to_file(const char *filepath, const ParsedSystem *ps, const LinearSystemSolution *sol);
void save_matrix_to_file(const char *filepath, const Matrix *A, const char *name);
void save_vectors_to_file(const char *filepath, const Matrix *V, const char *name);

#endif
