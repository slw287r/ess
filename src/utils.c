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

int strlen_wo_esc(const char *str)
{
	int length = 0;
	while (*str != '\0')
	{
		// Check if this is the start of an ANSI escape sequence (ESC [ ... m)
		if (*str == '\033' && *(str + 1) == '[')
		{
			// Move past \033[
			str += 2;
			// Skip until we reach 'm', which ends the ANSI sequence
			while (*str != '\0' && *str != 'm')
				str++;
			// Move past 'm' if we found it
			if (*str == 'm')
				str++;
			continue;
		}
		// Count visible characters
		length++;
		str++;
	}
	return length;
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

double nice_interval(const double x, const double n)
{
	double rough_intv = x / n;
	double magnitude = pow(10, floor(log10(rough_intv)));
	double nice_intv;
	if (rough_intv / magnitude <= 1)
		nice_intv = 1 * magnitude;
	else if (rough_intv / magnitude <= 2)
		nice_intv = 2 * magnitude;
	else if (rough_intv / magnitude <= 5)
		nice_intv = 5 * magnitude;
	else
		nice_intv = 10 * magnitude;
	return nice_intv;
}

void step_and_peak(const double n, sp_t *sp)
{
	double intv = nice_interval(n, 10);
	sp->step = intv;
	double i = intv, j = 0;
	while (i < n)
	{
		i += intv;
		++j;
	}
	if (fmod(j, 2))
		++j;
	sp->peak = intv * j;
}

void horiz(const int _n)
{
	struct winsize w;
	ioctl(0, TIOCGWINSZ, &w);
	int i, n = (w.ws_col >= _n) ? _n : w.ws_col;
	for (i = 0; i < n; ++i) fputs("\e[90m\xe2\x94\x80\e[0m", stdout);
	fputc('\n', stdout);
}
