#include "parser.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int is_var_char(char c){ return (c>='a'&&c<='z')||(c>='A'&&c<='Z'); }

static void add_var(char *vars, size_t *len, char v){
    for(size_t i=0;i<*len;i++) if(vars[i]==v) return;
    vars[*len]=v; (*len)++;
}

static void sort_vars(char *vars, size_t len){
    for(size_t i=0;i<len;i++)
        for(size_t j=i+1;j<len;j++)
            if(vars[j]<vars[i]){ char t=vars[i]; vars[i]=vars[j]; vars[j]=t; }
}

static void trim(char *s){
    // remove spaces around
    size_t n=strlen(s), i=0;
    while(i<n && isspace((unsigned char)s[i])) i++;
    size_t start=i;
    size_t end=n;
    while(end>start && isspace((unsigned char)s[end-1])) end--;
    memmove(s, s+start, end-start);
    s[end-start]='\0';
}

ParsedSystem parse_equations(char **lines, size_t m, const char *vars_hint){
    ParsedSystem ps;
    ps.A = mat_create(0,0);
    ps.b = vec_create(0);
    ps.vars=NULL;

    // Discover variable set
    char vset[64]; size_t vlen=0;
    if(vars_hint && *vars_hint){
        for(const char *p=vars_hint; *p; ++p) add_var(vset,&vlen,*p);
    }else{
        for(size_t i=0;i<m;i++){
            const char *p=lines[i];
            while(*p){
                if(is_var_char(*p)) add_var(vset,&vlen,*p);
                p++;
            }
        }
        sort_vars(vset, vlen);
    }
    ps.vars = (char*)malloc(vlen+1);
    memcpy(ps.vars, vset, vlen);
    ps.vars[vlen]='\0';

    ps.A = mat_create(m, vlen);
    ps.b = vec_create(m);

    // Parse each line
    for(size_t r=0;r<m;r++){
        char *buf = strdup(lines[r]);
        trim(buf);
        // split at '='
        char *eq = strchr(buf, '=');
        if(!eq){ free(buf); continue; }
        *eq='\0';
        char *left=buf;
        char *right=eq+1;
        trim(left); trim(right);

        // parse right constant
        ps.b.data[r] = atof(right);

        // parse left sum of terms
        // tokens: number*var or +/- var, handle signs
        int sign=+1;
        for(size_t i=0; left[i]; ){
            if(left[i]=='+'){ sign=+1; i++; continue; }
            if(left[i]=='-'){ sign=-1; i++; continue; }
            // read coefficient (optional)
            int hasnum=0;
            double coeff=0.0;
            size_t j=i;
            // allow spaces
            while(left[j] && isspace((unsigned char)left[j])) j++;
            i=j;
            // read number (possibly with decimal)
            char numbuf[64]; size_t nb=0;
            if((left[i]>='0'&&left[i]<='9') || left[i]=='.'){
                hasnum=1;
                while((left[i]>='0'&&left[i]<='9')||left[i]=='.'){
                    if(nb<63) numbuf[nb++]=left[i];
                    i++;
                }
                numbuf[nb]='\0';
                coeff = atof(numbuf);
                // optional '*'
                while(left[i] && isspace((unsigned char)left[i])) i++;
                if(left[i]=='*') i++;
            }else{
                coeff = 1.0;
            }
            while(left[i] && isspace((unsigned char)left[i])) i++;
            // expect variable
            if(!is_var_char(left[i])){
                // skip unknown char
                i++; continue;
            }
            char var = left[i++];
            // find column
            size_t col=0;
            for(; col<vlen; col++) if(ps.vars[col]==var) break;
            if(col<vlen){
                ps.A.data[r*ps.A.cols + col] += sign*coeff;
            }
        }
        free(buf);
    }

    return ps;
}

ParsedSystem parse_from_file(const char *filepath, const char *vars_hint){
    FILE *f = fopen(filepath, "r");
    if(!f){
        ParsedSystem empty; empty.A=mat_create(0,0); empty.b=vec_create(0); empty.vars=NULL; return empty;
    }
    char **lines=NULL; size_t m=0; size_t cap=0;
    char buf[512];
    while(fgets(buf, sizeof(buf), f)){
        // skip empty
        int only_space=1;
        for(char *p=buf; *p; ++p){ if(!isspace((unsigned char)*p)){ only_space=0; break; } }
        if(only_space) continue;
        if(m==cap){ cap=cap?cap*2:8; lines=(char**)realloc(lines, cap*sizeof(char*)); }
        lines[m++]=strdup(buf);
    }
    fclose(f);
    ParsedSystem ps = parse_equations(lines, m, vars_hint);
    for(size_t i=0;i<m;i++) free(lines[i]);
    free(lines);
    return ps;
}

void free_parsed_system(ParsedSystem *ps){
    mat_free(&ps->A);
    vec_free(&ps->b);
    if(ps->vars) free(ps->vars);
    ps->vars=NULL;
}
