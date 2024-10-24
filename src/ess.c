#include "ess.h"

int main(int argc, char *argv[])
{
	setenv("FONTCONFIG_PATH", "/etc/fonts", 1);
	arg_t *arg = calloc(1, sizeof(arg_t));
	if (argc == 1)
		usage();
	prs_arg(argc, argv, arg);
	// load reference index
	faidx_t *fai = NULL;
	if (!(fai = fai_load(arg->ref)))
		error("Error loading reference index of [%s]\n", arg->ref);
	float exp = exp_dmf(arg->ref), obs = obs_dmf(arg->in, fai);
	if (fai)
		fai_destroy(fai);
	// output
	FILE *fo = arg->out ? fopen(arg->out, "w") : stdout;
	if (!fo)
		error("Failed to create output file [%s]\n", arg->out ? arg->out : "stdout");
	fprintf(fo, "obs\t%f\nexp\t%f\no/e\t%f\n", exp, obs, exp ? obs / exp : 0);
	fclose(fo);
	free(arg);
	return 0;
}

static ko_longopt_t long_options[] = {
	{ "in",                        ko_required_argument, 'i' },
	{ "out",                       ko_required_argument, 'o' },
	{ "ref",                       ko_required_argument, 'r' },
	{ "help",                      ko_no_argument, 'h' },
	{ "version",                   ko_no_argument, 'v' },
	{ NULL, 0, 0 }
};

void prs_arg(int argc, char **argv, arg_t *arg)
{
	int c = 0;
	ketopt_t opt = KETOPT_INIT;
	const char *opt_str = "i:o:r:hv";
	while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >= 0)
	{
		switch (c)
		{
			case 'i': arg->in = opt.arg; break;
			case 'o': arg->out = opt.arg; break;
			case 'r': arg->ref = opt.arg; break;
			case 'h': usage(); break;
			case 'v':
				if (strlen(BRANCH_COMMIT))
					printf("%s (%s)\n", VERSION, BRANCH_COMMIT);
				else
					puts(VERSION);
				exit(EXIT_SUCCESS);
			case '?':
				printf("Invalid option: [%c]\n", opt.opt); exit(EXIT_SUCCESS);
			default:
				printf("Invalid option: [%s]", opt.arg); exit(EXIT_SUCCESS);
		}
	}
	if (!arg->in || access(arg->in, R_OK))
		error("Error: input bam is unspecified or inaccessible!\n");
	if (!(ends_with(arg->in, ".bam") && is_gzip(arg->in)))
		error("Oops! only bam input is supported. Invalid bam file [%s]\n", arg->in);
	if (!arg->ref)
	{
		static char ref[PATH_MAX] = {'\0'};
		samFile *fp = sam_open(arg->in, "r");
		bam_hdr_t *hdr = sam_hdr_read(fp);
		if (!hdr)
			error("Error reading input bam header!\n");
		bam_get_ref(hdr, ref);
		bam_hdr_destroy(hdr);
		hts_close(fp);
		if (!strlen(ref))
			error("Failed getting unspecified reference from bam file!\n");
		if (access(ref, R_OK))
			error("Reference [%s] in bam is inaccessible!\n", ref);
		arg->ref = ref;
	}
	else if (access(arg->ref, R_OK))
		error("Error: specified reference [%s] is inaccessible!\n", arg->ref);
	char *bai, *fai;
	asprintf(&bai, "%s.bai", arg->in);
	if (access(bai, R_OK))
	{
		free(bai);
		error("Error: bam's index file (.bai) is required, please use samtools sort and index to create it.\n");
	}
	asprintf(&fai, "%s.fai", arg->ref);
	if (access(fai, R_OK))
	{
		free(fai);
		error("Error: fasta's index file (.fai) is required, please use samtools faidx to create it.\n");
	}
}

int is_gzip(const char *fn)
{
	char buf[2];
	FILE *fp;
	int gzip = 0;
	if ((fp = fopen(fn, "rb")) == NULL)
		error("[ERROR] Unable to open file: %s\n", fn);
	if (fread(buf, 1, 2, fp) == 2)
		if (((int)buf[0] == 0x1f) && ((int)(buf[1]&0xFF) == 0x8b))
			gzip = 1;
	fclose(fp);
	return gzip;
}

int get_nm(const bam1_t *b)
{
	int32_t nm_i = 0;
	uint8_t *nm = bam_aux_get(b, "NM");
	if (nm) nm_i = bam_aux2i(nm);
	return nm_i;
}

void bam_get_cigar_str(bam1_t *b, kstring_t *s)
{
	int i;
	s->l = 0;
	uint32_t *c = bam_get_cigar(b);
	int n = b->core.n_cigar;
	for (i = 0; i < n; ++i)
	{
		int op = bam_cigar_op(c[i]);
		int ol = bam_cigar_oplen(c[i]);
		switch (op)
		{
			case BAM_CMATCH: // M
				ksprintf(s, "%dM", ol);
				break;
			case BAM_CHARD_CLIP: // H
				ksprintf(s, "%dH", ol);
				break;
			case BAM_CREF_SKIP: // S
			case BAM_CSOFT_CLIP: // S
				ksprintf(s, "%dS", ol);
				break;
			case BAM_CDEL: // D
				ksprintf(s, "%dD", ol);
				break;
			case BAM_CPAD: // P
				ksprintf(s, "%dP", ol);
				break;
			case BAM_CINS: // I
				ksprintf(s, "%dI", ol);
				break;
			 default:
				break;
		}
	}
}

int get_match_base(const bam1_t *b)
{
	int i = 0, match_base = 0;
	uint32_t *c = bam_get_cigar(b);
	int o = b->core.n_cigar;
	for(i = 0; i < o; ++i)
		if (bam_cigar_op(c[i]) == BAM_CMATCH)
			match_base += bam_cigar_oplen(c[i]);
	return match_base;
}

int get_mm(const bam1_t *b, const unsigned a, const unsigned z)
{
	int pos = 0, countout = 0;
	const bam1_core_t *c = &b->core;
	bool fwd = !(c->flag & BAM_FREVERSE);
	char *md = NULL;
	uint8_t *md_p = bam_aux_get(b, "MD");
	if (md_p) md = bam_aux2Z(md_p);
	while (*md)
	{
		if (isdigit(*md))
		{
			uint8_t *endptr;
			long i = strtol((char *)md, (char **)&endptr, 10);
			md = (char *)endptr;
			pos += i;
			continue;
		}
		if (*md == '^') // deletion.
		{
			while (*++md && !isdigit(*md))
				continue;
			continue;
		}
		// substitution
		if ((fwd && (pos <= a || pos >= c->l_qseq - z)) ||
			(!fwd && (pos <= z || pos >= c->l_qseq - a)))
		{
			char *mm = md;
			while (!isdigit(*mm++))
				++countout;
		}
		md++;
	}
	return get_nm(b) - countout;
}

float get_gc(const bam1_t *b)
{
	int i, l_qseq = b->core.l_qseq, bases[16] = {0};
	char *seq = (char *)bam_get_seq(b);
	for (i = 0; i < l_qseq; ++i)
		++bases[bam_seqi(seq, i)];
	return l_qseq ? 1.0 * (bases[2] + bases[4]) / l_qseq : 0.0f;
}

/*
 * @PG  ID:bwa       PN:bwa       VN:0.7.18  CL:bwa mem -at8 /path/to/ref.fa /path/to/fq.gz -R
 * @PG  ID:bwa-mem2  PN:bwa-mem2  VN:2.2.1a  CL:bwa-mem2 mem -v1 -zt8 -h64 /path/to/ref.fa /path/to/fq.gz
 */
void bam_get_ref(bam_hdr_t *h, char *ref)
{
	int i;
	kstring_t ks = {0, 0, 0};
	if (sam_hdr_find_line_id(h, "PG", "ID", "bwa-mem2", &ks) &&
			sam_hdr_find_line_id(h, "PG", "ID", "bwa", &ks))
		return;
	const char sfx[][16] = {".fasta.gz", ".fasta", ".fa.gz", ".fa", ".fna.gz", ".fna"};
	char *p = NULL, *q = NULL;
	for (i = 0; i < 6; ++i)
	{
		if ((p = strstr(ks.s, sfx[i])))
		{
			q = p;
			while (!isspace(*q--));
			strncpy(ref, q + 2, p - q + strlen(sfx[i]) - 2);
			break;
		}
	}
	if (ks.m) free(ks.s);
}

/* return the read, reverse complemented if necessary */
char *bam_get_seq_str(const bam1_t *b)
{
	int n, len = b->core.l_qseq + 1;
	char *read = calloc(1, len);
	char *seq = (char *)bam_get_seq(b);
	if (!read) return NULL;
	for (n = 0; n < b->core.l_qseq; ++n)
	{
		if (b->core.flag & BAM_FREVERSE)
			read[n] = seq_nt16_str[seq_comp_table[bam_seqi(seq,n)]];
		else
			read[n] = seq_nt16_str[bam_seqi(seq,n)];
	}
	if (b->core.flag & BAM_FREVERSE)
		reverse(read);
	return read;
}

void get_sc(const bam1_t *b, sc_t *sc)
{
	int k;
	for (k = 0; k < b->core.n_cigar; ++k)
	{
		int o = bam_get_cigar(b)[k] & BAM_CIGAR_MASK;
		int l = bam_get_cigar(b)[k] >> BAM_CIGAR_SHIFT;
		switch(o)
		{
			case BAM_CMATCH:
			case BAM_CHARD_CLIP:
			case BAM_CREF_SKIP:
			case BAM_CPAD:
			case BAM_CDEL:
			case BAM_CINS:
				break;
			case BAM_CSOFT_CLIP:
				if (k == 0) // leading SC
					sc->left = l;
				else // trailing SC
					sc->right = l;
				break;
			default:
				assert(0);
				break;
		}
	}
}

/* return the read, reverse complemented if necessary */
bool bam_is_cc(const bam1_t *b)
{
	int i, j = 0;
	char dm[3] = {'\0'};
	char *seq = (char *)bam_get_seq(b);
	sc_t sc = {0, 0};
	get_sc(b, &sc);
	// TODO check first in pair
	// ignore fragments with 5'-softclips
	if (b->core.flag & BAM_FREVERSE)
	{
		if (sc.right)
			return false;
		for (i = b->core.l_qseq - 1; i > b->core.l_qseq - 3; --i)
			dm[j++] = seq_nt16_str[seq_comp_table[bam_seqi(seq, i)]];
	}
	else
	{
		if (sc.left)
			return false;
		for (i = 0; i < 2; ++i)
			dm[j++] = seq_nt16_str[bam_seqi(seq, i)];
	}
	if (!strncmp(dm, "CC", 2) || !strncmp(dm, "GG", 2))
		return true;
	else
		return false;
}

float exp_dmf(const char *fa)
{
	float dmf = 0.0f;
	uint64_t tot = 0;
	kh_t *h = kh_init();
	gzFile fp = gzopen(fa, "r");
	kseq_t *ks = kseq_init(fp);
	while (kseq_read(ks) >= 0)
		seqtok(ks->seq.s, ks->seq.l, 2, h);
	kseq_destroy(ks);
	gzclose(fp);
	khint_t k;
	kh_foreach(h, k)
		tot += kh_val(h, k);
	if ((k = kh_get(h, DM_CC)) != kh_end(h) && tot)
		dmf = 1.0 * kh_val(h, k) / tot;
	kh_destroy(h);
	return dmf;
}

float obs_dmf(const char *bam, const faidx_t *fai)
{
	float dmf = 0.0f;
	uint64_t mtf = 0llu, tot = 0llu;
	// read bam for obs dmf
	samFile *fp = sam_open(bam, "r");
	bam_hdr_t *hdr = sam_hdr_read(fp);
	bam1_t *b = bam_init1();
	while(sam_read1(fp, hdr, b) >= 0)
	{
		bam1_core_t *c = &b->core;
		if (c->flag & (BAM_FQCFAIL | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
			continue;
		if (!faidx_has_seq(fai, sam_hdr_tid2name(hdr, c->tid)))
			continue;
		if (get_nm(b) >= 10) // FIXME
			continue;
		++tot;
		mtf += bam_is_cc(b);
	}
	bam_destroy1(b);
	bam_hdr_destroy(hdr);
	hts_close(fp);
	dmf = tot ? 1.0 * mtf / tot : 0.0f;
	return dmf;
}

void dump_ess_fn()
{
	printf("%50s\n", "obs(CC_end_freq + GG_end_freq)");
	printf("%19s", "ESS (O/E) = ");
	horiz(32);
	printf("%50s\n", "exp(CC_end_freq + GG_end_freq)");
}

void usage()
{
	struct winsize w;
	ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
	int wx = min(60, w.ws_col);
	horiz(wx);
	char title[] = "\e[1mCalculate End-Signature-Score (ESS) from BAM file\e[0m";
	int title_len = strlen("Calculate End-Signature-Score (ESS) from BAM file");
	printf("%*.*s\n", (int)((wx - title_len) / 2 + strlen(title)), (int)strlen(title), title);
	putchar('\n');
	dump_ess_fn();
	putchar('\n');
	horiz(wx);
	printf("%s \e[1mUsage\e[0m: \e[1;31m%s\e[0;0m \e[1;90m[options]\e[0;0m -i <bam> -o <tsv>\n", BUL, __progname);
	putchar('\n');
	puts(BUL " \e[1mOptions\e[0m:");
	puts("  -i, --in  \e[3mFILE\e[0m   Input BAM file with bai index (\e[31mrequired\e[0m)");
	puts("  -o, --out \e[3mSTR\e[0m    Output ESS value to file \e[90m[stdout]\e[0m");
	puts("  -r, --ref \e[3mFILE\e[0m   Reference fasta with fai index \e[90m[auto]\e[0m");
	putchar('\n');
	puts("  -h, --help       Display this message");
	puts("  -v, --version    Display program version");
	putchar('\n');
	puts(BUL " \e[1mContact\e[0m: \e[4mmeta@geneplus.cn\e[0m for support and bug report");
	horiz(wx);
	exit(1);
}
