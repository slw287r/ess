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
#include "izlib.h"

#define DM_CC 5
#define DM_MAX 12
#define MAX_RL 35

KSEQ_INIT(gzFile, gzread)
// argument struct
typedef struct
{
	char *in, *out, *ref;
} arg_t;

// bam auxillary struct
typedef struct
{
	samFile *fp;
	bam_hdr_t *hdr;
} aux_t;

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

void bam_get_ref(bam_hdr_t *h, char *ref);

/**
 * @brief  get query sequence in base characters
 *
 * @param b  pointer to bam1_t struct
 *
 * @return   query sequence; needs to be freed by the invoker
 */
char *bam_get_seq_str(const bam1_t *b);
// check if the 5'-end is the CC signature
void get_sc(const bam1_t *b, sc_t *sc);
bool bam_is_cc(const bam1_t *b);
void bam_get_cigar_str(bam1_t *b, kstring_t *s);
int get_nm(const bam1_t *b);

float exp_dmf(const char *fa);
float obs_dmf(const char *bam, const faidx_t *fai);
void prs_arg(int argc, char **argv, arg_t *arg);
void dump_ess_fn();
void usage();
