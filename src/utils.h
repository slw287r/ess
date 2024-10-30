#pragma once
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <unistd.h>
#include <inttypes.h>
#include <libgen.h>
#include <errno.h>
#include <time.h>
#include <math.h>

#include <ftw.h>
#include <sys/stat.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <sys/ioctl.h>
#ifdef __linux__
#include <sys/auxv.h>
#endif

#include "khashl.h"

extern const char *__progname;

#define BUL "\e[90m\xE2\x97\x8F\e[0m"
#define ARR "\e[90m\xE2\x97\x82\e[0m"
#define BUL "\e[90m\xE2\x97\x8F\e[0m"
#define PP fprintf(stderr, "%s\t%d\t<%s>\n", __FILE__, __LINE__, __func__);

#define error(msg, ...) do {                                 \
    fputs("\e[31m\xE2\x9C\x97\e[0m ", stderr);               \
    fprintf(stderr, msg, ##__VA_ARGS__);                     \
    if (errno) fprintf(stderr, ": %s", strerror(errno));     \
    fflush(stderr);                                          \
    exit(EXIT_FAILURE);                                      \
} while (0)

// map of unsigned long longs
KHASHL_MAP_INIT(KH_LOCAL, kh_t, kh, uint64_t, uint64_t, kh_hash_uint64, kh_eq_generic)

typedef struct
{
	double step, peak;
} sp_t;

// kmer functions
void seqtok(char *seq, int len, int k, kh_t *h);
void u64s(uint64_t x, char* s, size_t k);
uint64_t s64u(const char* s);
void print_kmer(const kh_t *h, const int kl, const int n);
void step_and_peak(const double n, sp_t *sp);
int strlen_wo_esc(const char *str);
bool ends_with(const char *str, const char *sfx);
void horiz(const int n);
