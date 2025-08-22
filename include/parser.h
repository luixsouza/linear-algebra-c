#ifndef PARSER_H
#define PARSER_H

#include "linalg.h"

typedef struct {
    Matrix A;
    Vector b;
    char *vars;
} ParsedSystem;

ParsedSystem parse_equations(char **lines, size_t m, const char *vars_hint);

ParsedSystem parse_from_file(const char *filepath, const char *vars_hint);

void free_parsed_system(ParsedSystem *ps);

#endif