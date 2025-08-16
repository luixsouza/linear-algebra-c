#ifndef PARSER_H
#define PARSER_H

#include "linalg.h"

// Parses equations like: "4x + 2y + z - 2t = 3"
// Supports variable names that are single letters (a-z or A-Z).
// Returns matrix A (m x n) and vector b. 'vars_out' string will contain the variable
// order found (e.g., "xyzt"). Caller must free A, b, and vars_out.
typedef struct {
    Matrix A;
    Vector b;
    char *vars; // null-terminated order of variables
} ParsedSystem;

// Parse from an array of lines (m lines). If vars_hint is non-NULL, it defines the variable order to expect.
ParsedSystem parse_equations(char **lines, size_t m, const char *vars_hint);

// Parse from file: each non-empty line is an equation.
ParsedSystem parse_from_file(const char *filepath, const char *vars_hint);

// Free ParsedSystem allocations
void free_parsed_system(ParsedSystem *ps);

#endif
