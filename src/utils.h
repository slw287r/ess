#pragma once
#ifndef _GNU_SOURCE
  #define _GNU_SOURCE 1
#endif

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

#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))
#define BUL "\e[90m\xE2\x97\x8F\e[0m"
#define ARR "\e[90m\xE2\x97\x82\e[0m"
#define BUL "\e[90m\xE2\x97\x8F\e[0m"
#define INF "\e[1;34m\xE2\x84\xb9\e[0;0m "
#define PP fprintf(stderr, "%s\t%d\t<%s>\n", __FILE__, __LINE__, __func__);

#define error(msg, ...) do {                                 \
    fputs("\e[31m\xE2\x9C\x97\e[0m ", stderr);               \
    fprintf(stderr, msg, ##__VA_ARGS__);                     \
    if (errno) fprintf(stderr, ": %s", strerror(errno));     \
    fflush(stderr);                                          \
    exit(EXIT_FAILURE);                                      \
} while (0)

#define warning(msg, ...) do {                               \
    fputs("\e[31m\xE2\x9A\xA0\e[0m ", stderr);               \
    fprintf(stderr, msg, ##__VA_ARGS__);                     \
    if (errno) fprintf(stderr, ": %s", strerror(errno));     \
    fflush(stderr);                                          \
} while (0)

typedef const char *cstr_t;
// set of unsigned ints
KHASHL_SET_INIT(KH_LOCAL, u32_t, u32, uint32_t, kh_hash_uint32, kh_eq_generic)
KHASHL_SET_INIT(KH_LOCAL, u64_t, u64, uint64_t, kh_hash_uint64, kh_eq_generic)
// map of unsigned long longs
KHASHL_MAP_INIT(KH_LOCAL, kh_t, kh, uint64_t, uint64_t, kh_hash_uint64, kh_eq_generic)
// unsigned to string hash table
KHASHL_MAP_INIT(KH_LOCAL, us_t, us, uint32_t, cstr_t, kh_hash_uint32, kh_eq_generic)
// unsigned to string array hash table
KHASHL_MAP_INIT(KH_LOCAL, ua_t, ua, uint32_t, char **, kh_hash_uint32, kh_eq_generic)
// unsigned to hash array hash table
KHASHL_MAP_INIT(KH_LOCAL, hh_t, hh, uint64_t, kh_t *, kh_hash_uint64, kh_eq_generic)

// set of strings
KHASHL_CSET_INIT(KH_LOCAL, ss_t, ss, cstr_t, kh_hash_str, kh_eq_str)
// string to uint32_t map: keys need to be freed
KHASHL_CMAP_INIT(KH_LOCAL, su_t, su, cstr_t, uint32_t, kh_hash_str, kh_eq_str)
// string to unsigned array hash table
KHASHL_CMAP_INIT(KH_LOCAL, sh_t, sh, cstr_t, u32_t *, kh_hash_str, kh_eq_str)
// string to string map: keys and values need to be freed
KHASHL_CMAP_INIT(KH_LOCAL, s2_t, s2, cstr_t, cstr_t, kh_hash_str, kh_eq_str)

// unsigned to struct array hash table
KHASHL_MAP_INIT(KH_LOCAL, hs_t, hs, uint32_t, u64_t *, kh_hash_uint32, kh_eq_generic)

// insert key-val pair into hash
void kh_ins(kh_t *h, const uint64_t n, const int64_t v);
void hh_ins(hh_t *h, const uint32_t j, const uint32_t n, const int32_t v);
// fetch value by key with check
uint64_t kh_xval(const kh_t *h, const uint64_t n);
bool kh_has(const kh_t *h, const uint64_t a);
void kh_asn(kh_t *h, const uint64_t a, const uint64_t b);
void ss_ins(ss_t *h, const char *a);
void s2_ins(s2_t *h, const char *a, const char *b);
void su_ins(su_t *h, const char *a, const uint32_t n);
uint64_t su_val(const su_t *h, const char *a);
void su_asn(su_t *h, const char *a, const uint32_t n);
void su_des(su_t *h);
void ss_des(ss_t *h);
void s2_des(s2_t *h);
void us_ins(us_t *h, const uint32_t n, const char *a);
void us_des(us_t *h);
void u32_ins(u32_t *h, uint32_t n);
void u64_ins(u64_t *h, uint64_t n);
bool u64_has(const u64_t *h, const uint64_t n);
uint32_t u32_val(const u32_t *h);
bool u32_has(const u32_t *h, const uint32_t n);

// kmer functions
void seqtok(char *seq, int len, int k, kh_t *h);
void u64s(uint64_t x, char* s, size_t k);
uint64_t s64u(const char* s);
void print_kmer(const kh_t *h, const int kl, const int n);

void chkfile(int num, ...);
bool begins_with(const char *str, const char *prefix);
bool ends_with(const char *str, const char *sfx);
void swap(char *xp, char *yp);
void reverse(char *str);
void complement(char *str);
char *revcom(const char *str);
char *strrstr(const char *haystack, const char *needle);
void horiz(const int n);
