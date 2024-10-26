#define _GNU_SOURCE 1
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <assert.h>
#include <stdbool.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

#include "kseq.h"
#include "ketopt.h"
#include "khashl.h"
#include "utils.h"
#include "version.h"
#include "iplot.h"
#include "izlib.h"

#define DM_CC 5
#define DM_MAX 12
#define MIN_RL 35

#define DEF_IS 500
#define MAX_IS 0xFFFF
#define MAX_TRIES 0xFFFF

KSEQ_INIT(gzFile, gzread)
// argument struct
typedef struct
{
	char *in, *out, *ref, *plot, *dep;
} arg_t;

int8_t seq_comp_table[16] = {
	0,  8, 4, 12,
	2, 10, 6, 14,
	1,  9, 5, 13,
	3, 11, 7, 15
};

typedef struct
{
	int8_t left, right;
} sc_t;

int is_gzip(const char *fn);
// sequence accession ID to sequence length from fa index
int acc2len(const faidx_t *fai, const char *acc);

void get_sname(const char *bam, char *sname);
void bam_get_ref(bam_hdr_t *h, char *ref);

// check if the 5'-end is the CC signature
void get_sc(const bam1_t *b, sc_t *sc);
bool bam_is_cc(const bam1_t *b);
int get_nm(const bam1_t *b);
// expected and observed dimmer frequencies
float exp_dmf(const char *fa);
float obs_dmf(const char *bam, const faidx_t *fai);
// insert size stats
void isize(const char *bam, const faidx_t *fai, int *is);
void lrsd(const int *is, const int n, sd_t *sd);
// parse command line arguments
void prs_arg(int argc, char **argv, arg_t *arg);
// print ESS equation
void dump_ess_fn();
void usage();
