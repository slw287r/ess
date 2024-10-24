#include "utils.h"

void kh_ins(kh_t *h, const uint64_t n, const int64_t v)
{
	khint_t k;
	int absent;
	if ((k = kh_get(h, n)) == kh_end(h))
	{
		k = kh_put(h, n, &absent);
		kh_val(h, k) = v;
	}
	else
		kh_val(h, k) += v;
}

void kh_asn(kh_t *h, const uint64_t a, const uint64_t b)
{
	khint_t k;
	int absent;
	if ((k = kh_get(h, a)) == kh_end(h))
	{
		k = kh_put(h, a, &absent);
		kh_val(h, k) = b;
	}
	else
		kh_val(h, k) = b;
}

bool kh_has(const kh_t *h, const uint64_t a)
{
	return kh_get(h, a) != kh_end(h);
}

uint64_t kh_xval(const kh_t *h, const uint64_t n)
{
	khint_t k;
	if ((k = kh_get(h, n)) != kh_end(h))
		return kh_val(h, k);
	else
		return 0;
}

void ss_ins(ss_t *h, const char *a)
{
	khint_t k;
	int absent;
	if ((k = ss_get(h, a)) == kh_end(h))
		k = ss_put(h, strdup(a), &absent);
}

void s2_ins(s2_t *h, const char *a, const char *b)
{
	khint_t k;
	int absent;
	if ((k = s2_get(h, a)) == kh_end(h))
	{
		k = s2_put(h, strdup(a), &absent);
		kh_val(h, k) = strdup(b);
	}
}

void su_ins(su_t *h, const char *a, const uint32_t n)
{
	khint_t k;
	int absent;
	if ((k = su_get(h, a)) == kh_end(h))
	{
		k = su_put(h, strdup(a), &absent);
		kh_val(h, k) = n;
	}
	else
		kh_val(h, k) += n;
}

uint64_t su_val(const su_t *h, const char *a)
{
	khint_t k;
	if ((k = su_get(h, a)) == kh_end(h))
		return 0;
	else
		return kh_val(h, k);
}

void su_asn(su_t *h, const char *a, const uint32_t n)
{
	khint_t k;
	int absent;
	if ((k = su_get(h, a)) == kh_end(h))
	{
		k = su_put(h, strdup(a), &absent);
		kh_val(h, k) = n;
	}
	else
		kh_val(h, k) = n;
}

void su_des(su_t *h)
{
	khint_t k;
	kh_foreach(h, k)
		free((char *)kh_key(h, k));
	su_destroy(h);
}

void ss_des(ss_t *h)
{
	khint_t k;
	kh_foreach(h, k)
		free((char *)kh_key(h, k));
	ss_destroy(h);
}

void s2_des(s2_t *h)
{
	khint_t k;
	kh_foreach(h, k)
	{
		free((char *)kh_key(h, k));
		free((char *)kh_val(h, k));
	}
	s2_destroy(h);
}

void us_ins(us_t *h, const uint32_t n, const char *a)
{
	khint_t k;
	int absent;
	if ((k = us_get(h, n)) == kh_end(h))
	{
		k = us_put(h, n, &absent);
		kh_val(h, k) = strdup(a);
	}
}

void us_des(us_t *h)
{
	khint_t k;
	kh_foreach(h, k)
		free((char *)kh_val(h, k));
	us_destroy(h);
}

void u32_ins(u32_t *h, uint32_t n)
{
	khint_t k;
	int absent;
	if ((k = u32_get(h, n)) == kh_end(h))
		u32_put(h, n, &absent);
}

void u64_ins(u64_t *h, uint64_t n)
{
	khint_t k;
	int absent;
	if ((k = u64_get(h, n)) == kh_end(h))
		u64_put(h, n, &absent);
}

bool u64_has(const u64_t *h, const uint64_t n)
{
	return u64_get(h, n) != kh_end(h);
}

uint32_t u32_val(const u32_t *h)
{
	khint_t k;
	uint32_t r = 0;
	kh_foreach(h, k)
	{
		r = kh_key(h, k);
		break;
	}
	return r;
}

bool u32_has(const u32_t *h, const uint32_t n)
{
	return u32_get(h, n) != kh_end(h);
}

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

void exec(void *_a)
{
	const char *a = (char *)_a;
	char commit_suicide[PATH_MAX];
	int r = system(a);
	if (WEXITSTATUS(r))
	{
		fprintf(stderr, "CMD error [%s]\n", a);
		if (errno)
			fprintf(stderr, "\e[31m%s\e[0m\n", strerror(errno));
		sprintf(commit_suicide, "kill -9 %d &>/dev/null", getpid());
		system(commit_suicide);
		exit(EXIT_FAILURE);
	}
}

void chomp(char *str)
{
	if (strchr(str, '\n'))
		str[strchr(str, '\n') - str] = '\0';
}

void chkfile(int num, ...)
{
	int i;
	va_list valist;
	va_start(valist, num);
	for (i = 0; i < num; ++i)
	{
		char *fn = va_arg(valist, char *);
		if (access(fn, R_OK))
		{
			va_end(valist);
			error("Error accessing required file [%s]\n", fn);
		}
	}
	va_end(valist);
}

bool begins_with(const char *str, const char *prefix)
{
	int str_len = strlen(str);
	int prefix_len = strlen(prefix);
	return (str_len >= prefix_len) && (0 == strncmp(str, prefix, prefix_len));
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

void swap(char *xp, char *yp)
{
	*xp = *xp ^ *yp;
	*yp = *xp ^ *yp;
	*xp = *xp ^ *yp;
}

void reverse(char *str)
{
	int i = 0, j = strlen(str) - 1;
	while (i < j)
		swap(str + i++, str + j--);
}

char comptab[] = {
   0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
  64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
 'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
  64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
 'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

void complement(char *str)
{
	while(*str)
	{
		*str = comptab[(int)*str];
		++str;
	}
}

char *revcom(const char *str)
{
	char *str_cpy = strdup(str);
	reverse(str_cpy);
	complement(str_cpy);
	return str_cpy;
}

char *strrstr(const char *haystack, const char *needle)
{
	char *r = NULL;
	if (!needle[0])
		return (char*)haystack + strlen(haystack);
	while (1)
	{
		char *p = strstr(haystack, needle);
		if (!p)
			return r;
		r = p;
		haystack = p + 1;
	}
}

void horiz(const int _n)
{
	struct winsize w;
	ioctl(0, TIOCGWINSZ, &w);
	int i, n = (w.ws_col >= _n) ? _n : w.ws_col;
	for (i = 0; i < n; ++i) fputs("\e[90m\xe2\x94\x80\e[0m", stdout);
	fputc('\n', stdout);
}
