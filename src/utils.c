#include "utils.h"

const uint8_t kmertochar[5] = { 'A', 'C', 'G', 'T', 'N' };

const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void seqtok(char *seq, int len, int k, kh_t *h)
{
	int i, l;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int absent, c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				uint64_t y = x[0] < x[1]? x[0] : x[1];
				khint_t itr = kh_put(h, y, &absent); // only add one strand!
				if (absent) kh_val(h, itr) = 0;
				++kh_val(h, itr);
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

void u64s(uint64_t x, char* s, size_t k)
{
	size_t i = 0;
	for (i = 0; i < k; ++i) {
		s[k - i - 1] = "ACGT"[(uint8_t) x & 0x3];
		x >>= 2;
	}
	s[i] = '\0';
}

uint64_t s64u(const char* s)
{
	int i, k = strlen(s);
	
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = 0, x[0] = x[1] = 0; i < k; ++i)
		{
			int c = seq_nt4_table[(uint8_t)s[i]];
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
		}
	return x[0] < x[1]? x[0] : x[1];
}

void print_kmer(const kh_t *h, const int kl, const int n)
{
	khint_t k;
	char s[kl + 1];
	kh_foreach(h, k)
	{
		u64s(kh_key(h, k), s, kl);
		//printf("%"PRIu64"\t%d\t%s\n", kh_key(h, k), kh_val(h, k), s);
		uint64_t cnt = kh_val(h, k);
		if (cnt >= n)
			printf("%"PRIu64"\t%s\t%"PRIu64"\n", kh_key(h, k), s, cnt);
	}
}

bool ends_with(const char *str, const char *sfx)
{
	int ret = 0;
	int str_len = strlen(str);
	int sfx_len = strlen(sfx);
	if ((str_len >= sfx_len) && (0 == strcasecmp(str + (str_len-sfx_len), sfx)))
		ret = 1;
	return ret;
}

void step_and_peak(const double n, sp_t *sp)
{
	int i;
	int m = log10(n);
	sp->peak = pow(10, m);
	printf("%d\t%f\n", m, sp->peak);
	for (i = 0; i < 10; ++i)
		if ((sp->peak += pow(10, m - 1)) >= n)
			break;
	sp->step = 5 * pow(10, m - 1) / 2;
}

void horiz(const int _n)
{
	struct winsize w;
	ioctl(0, TIOCGWINSZ, &w);
	int i, n = (w.ws_col >= _n) ? _n : w.ws_col;
	for (i = 0; i < n; ++i) fputs("\e[90m\xe2\x94\x80\e[0m", stdout);
	fputc('\n', stdout);
}
