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

void hh_ins(hh_t *h, const uint32_t j, const uint32_t n, const int32_t v)
{
	khint_t k;
	int absent;
	kh_t *h1 = NULL;
	if ((k = hh_get(h, j)) == kh_end(h))
	{
		h1 = kh_init();
		kh_ins(h1, n, v);
		k = hh_put(h, j, &absent);
		kh_val(h, k) = h1;
	}
	else
	{
		kh_t *h1 = kh_val(h, k);
		kh_ins(h1, n, v);
	}
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

int unlink_cb(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf)
{
	return remove(fpath);
}

int clean(const char *path)
{
	return nftw(path, unlink_cb, FOPEN_MAX, FTW_DEPTH | FTW_PHYS);
}

void chomp(char *str)
{
	if (strchr(str, '\n'))
		str[strchr(str, '\n') - str] = '\0';
}

int is_empty(const char *fn)
{
	int z = 0;
	if (!access(fn, R_OK))
	{
		FILE *fp = fopen(fn, "r");
		if (!fp) ++z;
		else if (fgetc(fp) == EOF) ++z;
		fclose(fp);
	}
	return z;
}

int is_dir(const char *path)
{
	struct stat statbuf;
	if (stat(path, &statbuf) != 0)
		return 0;
	return S_ISDIR(statbuf.st_mode);
}

int mkdir_p(const char *path)
{
	// check for regular file
	if (!access(path, F_OK) && !is_dir(path))
		error("Not a directory [%s]\n", path);
	if (!access(path, F_OK) && is_dir(path))
		return EXIT_SUCCESS;
	const size_t len = strlen(path);
	char _path[PATH_MAX];
	errno = 0;
	if (len > sizeof(_path)-1)
	{
		errno = ENAMETOOLONG;
		return EXIT_FAILURE;
	}
	strcpy(_path, path);
	char *p = 0;
	for (p = _path + 1; *p; p++)
	{
		if (*p == '/')
		{
			*p = '\0';
			if (mkdir(_path, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) != 0)
				if (errno != EEXIST)
					return EXIT_FAILURE;
			*p = '/';
		}
	}
	if (mkdir(_path, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) != 0)
		if (errno != EEXIST)
			return EXIT_FAILURE;
	return EXIT_SUCCESS;
}

bool is_sql3(const char *fn)
{
	bool sql3 = false;
	char line[17] = {0};
	FILE *fp = fopen(fn, "rb");
	if (!fp)
	{
		fprintf(stderr, "Error opening [%s]\n", fn);
		exit(1);
	}
	fread(line, sizeof(char), 16, fp);
	sql3 = (bool)(strncmp(line, "SQLite format 3", 16) == 0 ||
			strncmp(line, "\xd5\x78\x9e\xed\x58\x1d\x32\x12", 8) == 0);
	fclose(fp);
	return sql3;
}

char *decode(char *a)
{
	size_t i, s = 16;
	static const char password[16] = "invalid pointer";
	char *r = strdup(a);
	for (i = 0; i < strlen(a); ++i)
		r[i] ^= ~password[i % s];
	return r;
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

// http://www.concentric.net/~Ttwang/tech/inthash.htm
unsigned long mix(unsigned long a, unsigned long b, unsigned long c)
{
	a=a-b;  a=a-c;  a=a^(c >> 13);
	b=b-c;  b=b-a;  b=b^(a << 8);
	c=c-a;  c=c-b;  c=c^(b >> 13);
	a=a-b;  a=a-c;  a=a^(c >> 12);
	b=b-c;  b=b-a;  b=b^(a << 16);
	c=c-a;  c=c-b;  c=c^(b >> 5);
	a=a-b;  a=a-c;  a=a^(c >> 3);
	b=b-c;  b=b-a;  b=b^(a << 10);
	c=c-a;  c=c-b;  c=c^(b >> 15);
	return c;
}

/**
 * @brief Generate thread-safe random string of specified length
 *
 * @param s   the string to be returned
 * @param len the length of the returned string
 */
void gen_random(char *s, const int len)
{
	int i = 0;
	static const char ascii[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
	for (i = 0; i < len; ++i)
		s[i] = ascii[rand() % (sizeof(ascii) - 1)];
	s[len] = '\0';
}

char *tmpstr(const char *pfx)
{
	char randstr[] = "XXXXXX", *a = NULL;
	srand(mix(clock(), time(0), getpid()));
	gen_random(randstr, 6);
	if (pfx)
	{
		char *od = strdup(pfx);
		asprintf(&a, "%s/%s.%s", dirname(od), __progname, randstr);
		free(od);
	}
	else
		asprintf(&a, "./%s.%s", __progname, randstr);
	return a;
}

void flt_round(const unsigned n, char *flt)
{
	unsigned i = n;
	char *p = strchr(flt, '.');
	if (p)
	{
		if (!i)
			*p = '\0';
		else
		{
			while (*p++ && i--);
			printf("[%d]\t[%s]\n", i, p);
			if (i == -1)
				*p = '\0';
		}
	}
}

char *now(void)
{
	char buf[80];
	time_t _now = time(0);
	strftime(buf, sizeof(buf), "\e[35m%D %X\e[0m", localtime(&_now));
	return strdup(buf);
}

void info(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	char *right_now = now();
	fputs(INF, stderr);
	fprintf(stderr, "%s \e[35m", right_now);
	vfprintf(stderr, format, ap);
	fputs("\e[0m", stderr);
	fflush(stderr);
	va_end(ap);
	free(right_now);
}

void logging(const char *func, const char *msg)
{
	char *right_now = now();
	fputs(INF, stderr);
	fprintf(stderr, "%s ", right_now);
	fprintf(stderr, "[\e[36m%s\e[0m%s\e[36m%s]\e[0m %s\n", __progname, ARR, func, msg);
	fflush(stderr);
	free(right_now);
}

void horiz(const int _n)
{
	struct winsize w;
	ioctl(0, TIOCGWINSZ, &w);
	int i, n = (w.ws_col >= _n) ? _n : w.ws_col;
	for (i = 0; i < n; ++i) fputs("\e[90m\xe2\x94\x80\e[0m", stdout);
	fputc('\n', stdout);
}

void get_ip(char *ip)
{
	int fd, suc = 0;
	int if_num = 0;
	struct ifconf ifc;
	struct ifreq buf[MAX_IF];	
	ifc.ifc_len = sizeof(buf);
	ifc.ifc_buf = (caddr_t) buf;
	fd = socket(AF_INET, SOCK_DGRAM, 0);
	if (fd > 0)
	{
		if (ioctl(fd, SIOCGIFCONF, (char*)&ifc))
		{
			strncpy(ip, "n/a", MAX_IP - 1);
			return;
		}
		if_num = ifc.ifc_len / sizeof(struct ifreq);
		while (if_num-- > 0)
		{
			if (!strncmp(buf[if_num].ifr_name, "en", 2) &&
					!ioctl(fd, SIOCGIFFLAGS, (char*)&buf[if_num]) &&
					(buf[if_num].ifr_flags & IFF_RUNNING))
			{
				ioctl(fd, SIOCGIFADDR, (char*)&buf[if_num]);
				strncpy(ip, inet_ntoa(((struct sockaddr_in *)(&buf[if_num].ifr_addr))->sin_addr), MAX_IP - 1);
				suc = 1;
				break;
			}
		}
		close(fd);
	}
	if (!suc)
		strncpy(ip, "n/a", MAX_IP - 1);
}

void runtime(const time_t start, const char *func)
{
	time_t current;
	time(&current);
	int diff = (int)difftime(current, start);
	timer run_time;
	run_time.d = diff / 86400;
	run_time.h = diff % 86400 / 3600;
	run_time.m = diff % 86400 % 3600 / 60;
	run_time.s = diff % 86400 % 3600 % 60;
	if (func)
	{
		if (run_time.d)
			info("[\e[36m%s\e[0m%s\e[36m%s]\e[0m runtime: %02d:%02d:%02d:%02d\n",
					__progname, ARR, func, run_time.d, run_time.h, run_time.m, run_time.s);
		else
			info("[\e[36m%s\e[0m%s\e[36m%s]\e[0m runtime: %02d:%02d:%02d\n",
					__progname, ARR, func, run_time.h, run_time.m, run_time.s);
	}
	else
	{
		if (run_time.d)
			info("[\e[36m%s\e[0m] runtime: %02d:%02d:%02d:%02d\n",
					__progname, run_time.d, run_time.h, run_time.m, run_time.s);
		else
			info("[\e[36m%s\e[0m] runtime: %02d:%02d:%02d\n",
					 __progname, run_time.h, run_time.m, run_time.s);
	}
}

void dmp_cmd(const int argc, char **argv, const char *version, const char *prefix)
{
	int i = 0;
	FILE *fp = stderr;
	char host[NAME_MAX], fn[PATH_MAX], *p = NULL;
	if (prefix)
	{
		if ((p = strrchr(prefix, '/')))
		{
			if ((p = strrchr(p + 1, '.')))
			{
				memcpy(fn, prefix, p - prefix);
				memcpy(fn + (p - prefix), ".cmd\0", 5);
			}
			else
				sprintf(fn, "%s.cmd", prefix);
		}
		else if ((p = strrchr(prefix, '.')))
		{
			memcpy(fn, prefix, p - prefix);
			memcpy(fn + (p - prefix), ".cmd\0", 5);
		}
		else
			sprintf(fn, "%s.cmd", prefix);
		fp = fopen(fn, "w");
	}
	char _date[64], _time[64], last_cmd_char = 0, ip[MAX_IP] = {'\0'};
	time_t now = time(0);
	strftime(_date, sizeof(_date), "%D", localtime(&now));
	strftime(_time, sizeof(_time), "%X", localtime(&now));
	gethostname(host, NAME_MAX);
	if ((p = strstr(host, ".local")))
		*p = '\0';
	get_ip(ip);
	fprintf(fp, "#<cmd host=%s, ip=%s, version=%s, date=%s, time=%s>\n", host, ip, version, _date, _time);
	for (; i < argc; ++i)
	{
		if (*(argv[i]) == '-')
		{
			if (i != argc - 1 && *argv[i + 1] == '-')
			{
				if (*(argv[i] + 1) == '-')
				{
					fprintf(fp, " %s \\\n", argv[i]);
					last_cmd_char = '\n';
				}
				else
				{
					fprintf(fp, "  %s \\\n", argv[i]);
					last_cmd_char = '\n';
				}
			}
			else
			{
				fprintf(fp, "  %s ", argv[i]);
				last_cmd_char = ' ';
			}
		}
		else
		{
			if (!access(argv[i], F_OK))
			{
				char rp[PATH_MAX + 1];
				if (*(argv[i]) == '~')
				{
					if (strlen(argv[i]) == 1 || (strlen(argv[i]) == 2 && *(argv[i] + 1) == '/'))
						strcpy(rp, getenv("HOME"));
					else
						snprintf(rp, PATH_MAX, "%s/%s", getenv("HOME"), argv[i] + 2);
				}
				else
					realpath(strdup(argv[i]), rp);
				fprintf(fp, "%s %s\n", rp, (i + 1 == argc) ? "" : "\\");
				last_cmd_char = '\n';
			}
			else
			{
				if (*(argv[i]) == '~')
				{
					char rp[PATH_MAX + 1];
					if (strlen(argv[i]) == 1 || (strlen(argv[i]) == 2 && *(argv[i] + 1) == '/'))
						strcpy(rp, getenv("HOME"));
					else
						snprintf(rp, PATH_MAX, "%s/%s", getenv("HOME"), argv[i] + 2);
					fprintf(fp, "%s %s\n", rp, (i + 1 == argc) ? "" : "\\");
					last_cmd_char = '\n';
				}
				else
				{
					fprintf(fp, "%s %s\n", argv[i], (i + 1 == argc) ? "" : "\\");
					last_cmd_char = '\n';
				}
			}
		}
	}
	if (last_cmd_char == ' ')
		fputc('\n', fp);
	fprintf(fp, "#</cmd>\n");
	fclose(fp);
}
