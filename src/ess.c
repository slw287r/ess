#include "ess.h"

int main(int argc, char *argv[])
{
	setenv("FONTCONFIG_PATH", "/etc/fonts", 1);
	arg_t *arg = calloc(1, sizeof(arg_t));
	arg->mis = MIN_IS;
	arg->xis = DEF_IS;
	if (argc == 1)
		usage();
	prs_arg(argc, argv, arg);
	// load reference index
	faidx_t *fai = NULL;
	if (!(fai = fai_load(arg->ref)))
		error("Error loading reference index of [%s]\n", arg->ref);
	float exp = exp_dmf(arg->ref), obs = obs_dmf(arg->in, arg->mis, MAX_IS, fai);
	if (arg->plot)
	{
		int i, tot = 0;
		int *is = calloc(MAX_IS + 1, sizeof(int));
		double *cis = calloc(MAX_IS + 1, sizeof(double));
		sd_t sd = {0};
		isize(arg->in, fai, arg->mis, MAX_IS + 1, is);
		/* dbg is
		for (i = 0; i <= MAX_IS; ++i)
			printf("%d\t%d\n", i, is[i]);
		*/
		for (i = 0; i <= MAX_IS; ++i)
			cis[i] = (tot += is[i]);
		for (i = 0; i <= MAX_IS; ++i)
			cis[i] = 1 - (cis[i] /= cis[MAX_IS - 1]);
		lrsd(is, MAX_IS, &sd);
		cairo_surface_t *sf = NULL;
		if (ends_with(arg->plot, ".svg"))
			sf = cairo_svg_surface_create(arg->plot, WIDTH * 1.02, HEIGHT);
		else if (ends_with(arg->plot, ".png"))
			sf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, WIDTH, HEIGHT);
		cairo_t *cr = cairo_create(sf);
		cairo_set_antialias(cr, CAIRO_ANTIALIAS_BEST);
		char sname[NAME_MAX];
		get_sname(arg->in, sname);
		i = MAX_IS;
		while (!is[i--]);
		++i;
		do_drawing(cr, is, cis, fmin(arg->xis, i), &sd, sname, arg->sub);
		// clean canvas
		if (ends_with(arg->plot, ".png"))
			cairo_surface_write_to_png(cairo_get_target(cr), arg->plot);
		cairo_surface_destroy(sf);
		cairo_destroy(cr);
		// clean sizes
		free(is);
		free(cis);
	}
	if (fai)
		fai_destroy(fai);
	// output
	FILE *fo = arg->out ? fopen(arg->out, "w") : stdout;
	if (!fo)
		error("Failed to create output file [%s]\n", arg->out ? arg->out : "stdout");
	fprintf(fo, "exp\t%f\nobs\t%f\no/e\t%f\n", exp, obs, exp ? obs / exp: 0);
	fclose(fo);
	free(arg);
	return 0;
}

static ko_longopt_t long_options[] = {
	{ "in",                        ko_required_argument, 'i' },
	{ "out",                       ko_required_argument, 'o' },
	{ "ref",                       ko_required_argument, 'r' },
	{ "mis",                       ko_required_argument, 'm' },
	{ "plot",                      ko_required_argument, 'p' },
	{ "sub",                       ko_required_argument, 's' },
	{ "help",                      ko_no_argument, 'h' },
	{ "version",                   ko_no_argument, 'v' },
	{ NULL, 0, 0 }
};

void prs_arg(int argc, char **argv, arg_t *arg)
{
	int c = 0;
	char *p = NULL;
	ketopt_t opt = KETOPT_INIT;
	const char *opt_str = "i:o:r:m:p:s:hv";
	while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >= 0)
	{
		switch (c)
		{
			case 'i': arg->in = opt.arg; break;
			case 'o': arg->out = opt.arg; break;
			case 'r': arg->ref = opt.arg; break;
			case 'm':
				if ((p = strchr(opt.arg, ',')))
				{
					if (p != opt.arg)
						arg->mis = atoi(opt.arg);
					arg->xis = atoi(p + 1);
				}
				else
					arg->mis = atoi(opt.arg);
				break;
			case 'p': arg->plot = opt.arg; break;
			case 's': arg->sub = opt.arg; break;
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
	if (arg->mis < 0)
		error("Error: invalid mismatch value specified [%d]\n", arg->mis);
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
	free(bai);
	free(fai);
	if (arg->plot && !(ends_with(arg->plot, ".png") || ends_with(arg->plot, ".svg")))
		error("Error: unsupported plot format: [%s]\n", arg->plot);
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

void get_sname(const char *bam, char *sam)
{
	char *p = NULL;
	strncpy(sam, (p = strrchr(bam, '/')) ? p + 1 : bam, NAME_MAX);
	*strrchr(sam, '.') = '\0';
}

int get_nm(const bam1_t *b)
{
	int32_t nm_i = 0;
	uint8_t *nm = bam_aux_get(b, "NM");
	if (nm) nm_i = bam_aux2i(nm);
	return nm_i;
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
	const char sfx[][16] = {".fasta.gz ", ".fasta ", ".fa.gz ", ".fa ", ".fna.gz ", ".fna "};
	char *p = NULL, *q = NULL;
	for (i = 0; i < 6; ++i)
	{
		if ((p = strstr(ks.s, sfx[i])))
		{
			q = p;
			while (!isspace(*q--));
			strncpy(ref, q + 2, p - q + strlen(sfx[i]) - 3);
			break;
		}
	}
	if (ks.m) free(ks.s);
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

bool bam_is_cc(const bam1_t *b, bool *skip)
{
	int i, j = 0;
	*skip = false;
	char dm[3] = {'\0'};
	char *seq = (char *)bam_get_seq(b);
	sc_t sc = {0, 0};
	get_sc(b, &sc);
	// ignore fragments with 5'-softclips
	if (b->core.flag & BAM_FREVERSE)
	{
		if (sc.right)
		{
			*skip = true;
			return false;
		}
		for (i = b->core.l_qseq - 1; i > b->core.l_qseq - 3; --i)
			dm[j++] = seq_nt16_str[seq_comp_table[bam_seqi(seq, i)]];
	}
	else
	{
		if (sc.left)
		{
			*skip = true;
			return false;
		}
		for (i = 0; i < 2; ++i)
			dm[j++] = seq_nt16_str[bam_seqi(seq, i)];
	}
	return !strncmp(dm, "CC", 2);
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

float obs_dmf(const char *bam, const int mis, const int xis, const faidx_t *fai)
{
	bool skip;
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
		if (llabs(c->isize) > xis)
			continue;
		if (c->flag & BAM_FPAIRED)
		{
			if (c->mtid != c->tid || llabs(c->isize) < mis)
				continue;
		}
		else
		{
			sc_t sc = {0, 0};
			get_sc(b, &sc);
			int is = c->l_qseq - sc.left - sc.right;
			if (is < mis)
				continue;
		}
		if (get_nm(b) > c->l_qseq * 0.94)
			continue;
		mtf += bam_is_cc(b, &skip) * (c->flag & BAM_FPAIRED ? 1 : 2);
		tot += !skip * (bool)(c->flag & BAM_FPAIRED ? c->flag & BAM_FREAD1 : true);
	}
	bam_destroy1(b);
	bam_hdr_destroy(hdr);
	hts_close(fp);
	dmf = tot ? 1.0 * mtf / tot : 0.0f;
	return dmf;
}

void isize(const char *bam, const faidx_t *fai, const int mis, const int xis, int *is)
{
	int is1, tries = MAX_TRIES;
	samFile *fp = sam_open(bam, "r");
	bam_hdr_t *hdr = sam_hdr_read(fp);
	bam1_t *b = bam_init1();
	while(sam_read1(fp, hdr, b) >= 0)
	{
		bam1_core_t *c = &b->core;
		if (c->flag & (BAM_FQCFAIL | BAM_FREAD2 | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
			continue;
		if (!faidx_has_seq(fai, sam_hdr_tid2name(hdr, c->tid)))
			continue;
		if (llabs(c->isize) > xis)
			continue;
		if (c->flag & BAM_FPAIRED)
		{
			if (c->mtid != c->tid || llabs(c->isize) < mis)
				continue;
			if ((is1 = c->isize) <= 0 || is1 > xis)
				continue;
			++is[(int)fmin(is1, xis)];
			if (!--tries) break;
		}
		else
		{
			sc_t sc = {0, 0};
			get_sc(b, &sc);
			if ((is1 = c->l_qseq - sc.left - sc.right) <= xis && is1 >= mis)
			{
				++is[is1];
				if (!--tries) break;
			}
		}
	}
	bam_destroy1(b);
	bam_hdr_destroy(hdr);
	hts_close(fp);
}

void lrsd(const int *is, const int n, sd_t *sd)
{
	int i, j = 0, t = 0, pk = 0, cm = 0, rm = 0;
	// record insert sizes
	double *x = calloc(n, sizeof(double)), *y = calloc(n, sizeof(double));
	for (i = 1; i < n; ++i)
	{
		x[i] = i;
		y[i] = is[i];
		t += is[i];
		if (is[i] > cm)
		{
			cm = is[i];
			pk = i;
		}
	}
	t /= 2;
	for (i = 1; i < n; ++i)
	{
		if(is[i])
		{
			if(j <= t && t <= j + is[i])
				break;
			j += is[i];
		}
	}
	int lc = 0, rc = 0;
	double lsd = 0, rsd = 0;
	for (i = 1, j = 0; i < n; ++i)
	{
		rm = i;
		j += is[i];
		if (j >= t * 2 * 0.98) break;
	}
	for (i = 1; i <= rm; ++i)
	{
		int delta = i - pk;
		if (delta < 0)
		{
			lsd += is[i] * pow(delta, 2);
			lc += is[i];
		}
		else
		{
			rsd += is[i] * pow(delta, 2);
			rc += is[i];
		}
	}
	lc = lc > 1 ? lc : 1;
	rc = rc > 1 ? rc : 1;
	lsd = sqrt(lsd / lc);
	rsd = sqrt(rsd / rc);
	sd->pk = pk;
	sd->lsd = lsd;
	sd->rsd = rsd;
}

void dump_ess_fn()
{
	printf("%50s\n", "obs(CC_end_freq + GG_end_freq)");
	printf("%29s", "\e[3mESS\e[0m (O/E) = ");
	horiz(32);
	printf("%50s\n", "exp(CC_end_freq + GG_end_freq)");
}

void usage()
{
	struct winsize w;
	ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
	int wx = fmin(65, w.ws_col);
	horiz(wx);
	char title[] = "\e[1mCalculate End-Signature-Score (ESS) from BAM file\e[0m";
	int title_len = strlen_wo_esc(title);
	printf("%*.*s\n", (int)((wx - title_len) / 2 + strlen(title)), (int)strlen(title), title);
	putchar('\n');
	dump_ess_fn();
	putchar('\n');
	horiz(wx);
	printf("%s \e[1mUsage\e[0m: \e[1;31m%s\e[0;0m \e[1;90m[options]\e[0;0m -i <bam> -o <tsv>\n", BUL, __progname);
	putchar('\n');
	puts(BUL " \e[1mOptions\e[0m:");
	puts("  -i, --in  \e[3mFILE\e[0m     Input BAM file with bai index");
	puts("  -o, --out \e[3mSTR\e[0m      Output ESS value to file \e[90m[stdout]\e[0m");
	puts("  -r, --ref \e[3mFILE\e[0m     Reference fasta with fai index \e[90m[auto]\e[0m");
	puts("  -p, --plot \e[3mFILE\e[0m    Insert size plot png file \e[90m[none]\e[0m");
	printf("  -m, --mis \e[3mINT\e[90m,INT\e[0m\e[0m  Minimum\e[90m,Maximum\e[0m insert size to plot \e[90m[%d,%d]\e[0m\n", MIN_IS, DEF_IS);
	puts("  -s, --sub \e[3mFILE\e[0m     Sub-title of insert size plot \e[90m[none]\e[0m");
	putchar('\n');
	puts("  -h, --help         Display this message");
	puts("  -v, --version      Display program version");
	putchar('\n');
	puts(BUL " \e[1mContact\e[0m: \e[4mmeta@geneplus.cn\e[0m for support and bug report");
	horiz(wx);
	exit(1);
}
