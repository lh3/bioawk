#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "awk.h"

int bio_flag = 0, bio_fmt = BIO_NULL;

static const char *col_defs[][15] = { /* FIXME: this is convenient, but not memory efficient. Shouldn't matter. */
	{"header", NULL},
	{"bed", "chrom", "start", "end", "name", "score", "strand", "thickstart", "thickend", "rgb", "blockcount", "blocksizes", "blockstarts", NULL},
	{"sam", "qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual", NULL},
	{"vcf", "chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", NULL},
	{"gff", "seqname", "source", "feature", "start", "end", "score", "filter", "strand", "group", "attribute", NULL},
	{"fastx", "name", "seq", "qual", "comment", NULL},
	{NULL}
};

static const char *tab_delim = "nyyyyyn", *hdr_chr = "\0#@##\0\0";

/************************
 * Setting column names *
 ************************/

static void set_colnm_aux(const char *p, int col)
{
	const char *q;
	char *r = 0;
	Cell *x;
	for (q = p; *q; ++q) /* test if there are punctuations */
		if (ispunct(*q) && *q != '_') break;
	if (*q || isdigit(*p)) { /* there are punctuations or the first is digit */
		char *qq;
		r = malloc(strlen(p) + 2);
		if (isdigit(*p)) {
			*r = '_';
			strcpy(r + 1, p);
		} else strcpy(r, p);
		for (qq = r; *qq; ++qq)
			if (ispunct(*qq)) *qq = '_';
		q = r;
	} else q = p;
	if ((x = lookup(q, symtab)) != NULL) /* do not add if not appear in the program */
		setfval(x, (Awkfloat)col);
	if (r) free(r);
}

int bio_get_fmt(const char *s)
{
	int i, j;
	if (strcmp(s, "hdr") == 0) return BIO_HDR;
	for (i = 0; col_defs[i][0]; ++i)
		if (strcmp(s, col_defs[i][0]) == 0) return i;
	for (i = 1; col_defs[i][0]; ++i) {
		printf("%s:\n\t", col_defs[i][0]);
		for (j = 1; col_defs[i][j]; ++j)
			printf("%d:%s ", j, col_defs[i][j]);
		putchar('\n');
	}
	return BIO_NULL;
}

int bio_skip_hdr(const char *r)
{
	if (bio_fmt <= BIO_HDR) return 0;
	if (*r && *r == hdr_chr[bio_fmt]) {
		if (bio_flag & BIO_SHOW_HDR) puts(r);
		return 1;
	} else return 0;
}

void bio_set_colnm()
{
	int i;
	if (bio_fmt == BIO_NULL) {
		return;
	} else if (bio_fmt == BIO_HDR) {
		char *p, *q, c;
		for (p = record; *p && isspace(*p); ++p); /* skip leading spaces */
		for (i = 1, q = p; *q; ++q) {
			if (!isspace(*q)) continue;
			c = *q; /* backup the space */
			*q = 0; /* terminate the field */
			set_colnm_aux(p, i);
			*q = c; /* change back */
			++i;
			for (p = q + 1; *p && isspace(*p); ++p); /* skip contiguous spaces */
			q = p;
		}
		set_colnm_aux(p, i); /* the last column */
	} else {
		for (i = 0; col_defs[bio_fmt][i] != NULL; ++i)
			set_colnm_aux(col_defs[bio_fmt][i], i);
		if (tab_delim[bio_fmt] == 'y') *FS = *OFS = "\t";
	}
}

/**********************
 * Built-in functions *
 **********************/

static char comp_tab[] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

static float q_int2real[128];

#define tempfree(x)	  if (istemp(x)) tfree(x); else

Cell *bio_func(int f, Cell *x, Node **a)
{
	Cell *y, *z;
	y = gettemp();
	if (f == BIO_FAND) {
		if (a[1]->nnext == 0) {
			WARNING("and requires two arguments; returning 0.0");
			setfval(y, 0.0);
		} else {
			z = execute(a[1]->nnext);
			setfval(y, (Awkfloat)((long)getfval(x) & (long)getfval(z))); /* FIXME: does (long) always work??? */
			tempfree(z);
		}
	} else if (f == BIO_FOR) {
		if (a[1]->nnext == 0) {
			WARNING("or requires two arguments; returning 0.0");
			setfval(y, 0.0);
		} else {
			z = execute(a[1]->nnext);
			setfval(y, (Awkfloat)((long)getfval(x) | (long)getfval(z)));
			tempfree(z);
		}
	} else if (f == BIO_FXOR) {
		if (a[1]->nnext == 0) {
			WARNING("xor requires two arguments; returning 0.0");
			setfval(y, 0.0);
		} else {
			z = execute(a[1]->nnext);
			setfval(y, (Awkfloat)((long)getfval(x) ^ (long)getfval(z)));
			tempfree(z);
		}
	} else if (f == BIO_FREVERSE) {
		char *buf = getsval(x);
		int i, l, tmp;
		l = strlen(buf);
		for (i = 0; i < l>>1; ++i)
			tmp = buf[i], buf[i] = buf[l-1-i], buf[l-1-i] = tmp;
		setsval(y, buf);
	} else if (f == BIO_FREVCOMP) {
		char *buf;
		int i, l, tmp;
		buf = getsval(x);
		l = strlen(buf);
		for (i = 0; i < l>>1; ++i)
			tmp = comp_tab[(int)buf[i]], buf[i] = comp_tab[(int)buf[l-1-i]], buf[l-1-i] = tmp;
		if (l&1) buf[l>>1] = comp_tab[(int)buf[l>>1]];
		setsval(y, buf);
	} else if (f == BIO_FGC) {
		char *buf;
		int i, l, gc = 0;
		buf = getsval(x);
		l = strlen(buf);
		if (l) { /* don't try for empty strings */
			for (i = 0; i < l; ++i)
				if (buf[i] == 'g' || buf[i] == 'c' ||
					buf[i] == 'G' || buf[i] == 'C')
					gc++;
			setfval(y, (Awkfloat)gc / l);
		}
	} else if (f == BIO_FMEANQUAL) {
		char *buf;
		int i, l, total_qual = 0;
		buf = getsval(x);
		l = strlen(buf);
		if (l) { /* don't try for empty strings */
			for (i = 0; i < l; ++i)
				total_qual += buf[i] - 33;
			setfval(y, (Awkfloat)total_qual / l);
		}
	} else if (f == BIO_FTRIMQ) {
		char *buf;
		double thres = 0.05, s, max;
		int i, l, tmp, beg, end;
		Cell *u = 0, *v = 0;
		if (a[1]->nnext) {
			u = execute(a[1]->nnext); /* begin */
			if (a[1]->nnext->nnext) {
				v = execute(a[1]->nnext->nnext); /* end */
				if (a[1]->nnext->nnext->nnext) {
					z = execute(a[1]->nnext->nnext->nnext);
					thres = getfval(z); /* user defined threshold */
					tempfree(z);
				}
			}
		}
		buf = getsval(x);
		l = strlen(buf);
		if (q_int2real[0] == 0.) /* to initialize */
			for (i = 0; i < 128; ++i)
				q_int2real[i] = pow(10., -(i - 33) / 10.);
		for (i = 0, beg = tmp = 0, end = l, s = max = 0.; i < l; ++i) {
			int q = buf[i];
			if (q < 36) q = 36;
			if (q > 127) q = 127;
			s += thres - q_int2real[q];
			if (s > max) max = s, beg = tmp, end = i + 1;
			if (s < 0) s = 0, tmp = i + 1;
		}
		if (u) { setfval(u, beg); tempfree(u); } /* 1-based position; as substr() is 1-based. */
		if (v) { setfval(v, end); tempfree(v); }
	} else if (f == BIO_FQUALCOUNT) {
		if (a[1]->nnext == 0) {
			WARNING("xor requires two arguments; returning 0.0");
			setfval(y, 0.0);
		} else {
			char *buf;
			int i, l, thres, cnt = 0;
			buf = getsval(x);
			l = strlen(buf);
			z = execute(a[1]->nnext); /* threshold */
			thres = (int)(getfval(z) + .499);
			for (i = 0; i < l; ++i)
				if (buf[i] - 33 >= thres) ++cnt;
			setfval(y, (Awkfloat)cnt);
		}
	} else if (f == BIO_FFASTX) {
		if (a[1]->nnext == 0) {
			FATAL("fastx requires at least two arguments");
		} else {
			char *buf, *name, *seq, *qual;
			int bufsz=3*recsize;
			int has_qual;
			z = execute(a[1]->nnext);
			if ((has_qual = a[1]->nnext->nnext != 0)) {
				y = execute(a[1]->nnext->nnext);
				qual = getsval(y);
			}

			if ((buf = (char *) malloc(bufsz)) == NULL)
				FATAL("out of memory in fastx");

			name = getsval(x);
			seq = getsval(z);
			if (!has_qual) {
				sprintf(buf, ">%s\n%s", name, seq);
			} else {
				if (strlen(seq) != strlen(qual))
					WARNING("fastx arguments seq and qual are not same length");
				sprintf(buf, "@%s\n%s\n+\n%s", name, seq, qual); 
			}
			setsval(y, buf);
			free(buf);
		}
	} /* else: never happens */
	return y;
}

/************************
 * getrec() replacement *
 ************************/

#include <zlib.h> /* FIXME: it would be better to drop this dependency... */
#include "kseq.h"
KSEQ_INIT2(, gzFile, gzread)

static gzFile g_fp;
static kseq_t *g_kseq;
static int g_firsttime = 1, g_is_stdin = 0;
static kstring_t g_str;

int bio_getrec(char **pbuf, int *psize, int isrecord)
{
	extern Awkfloat *ARGC;
	extern int argno, recsize;
	extern char *file;
	extern Cell **fldtab;

	int i, c, saveb0, dret, bufsize = *psize, savesize = *psize;
	char *p, *buf = *pbuf;
	if (g_firsttime) { /* mimicing initgetrec() in lib.c */
		g_firsttime = 0;
		for (i = 1; i < *ARGC; i++) {
			p = getargv(i); /* find 1st real filename */
			if (p == NULL || *p == '\0') {	/* deleted or zapped */
				argno++;
				continue;
			}
			if (!isclvar(p)) {
				setsval(lookup("FILENAME", symtab), p);
				goto getrec_start;
			}
			setclvar(p);	/* a commandline assignment before filename */
			argno++;
		}
		g_fp = gzdopen(fileno(stdin), "r"); /* no filenames, so use stdin */
		g_kseq = kseq_init(g_fp);
		g_is_stdin = 1;
	}

getrec_start:
	if (isrecord) {
		donefld = 0; /* these are defined in lib.c */
		donerec = 1;
	}
	saveb0 = buf[0];
	buf[0] = 0; /* this is effective at the end of file */
	while (argno < *ARGC || g_is_stdin) {
		if (g_kseq == 0) { /* have to open a new file */
			file = getargv(argno);
			if (file == NULL || *file == '\0') { /* deleted or zapped */
				argno++;
				continue;
			}
			if (isclvar(file)) {	/* a var=value arg */
				setclvar(file);
				argno++;
				continue;
			}
			*FILENAME = file;
			if (*file == '-' && *(file+1) == '\0') {
				g_fp = gzdopen(fileno(stdin), "r");
				g_kseq = kseq_init(g_fp);
				g_is_stdin = 1;
			} else {
				if ((g_fp = gzopen(file, "r")) == NULL)
					FATAL("can't open file %s", file);
				g_kseq = kseq_init(g_fp);
				g_is_stdin = 0;
			}
			setfval(fnrloc, 0.0);
		}
		if (bio_fmt != BIO_FASTX) {
			c = ks_getuntil(g_kseq->f, **RS, &g_str, &dret);
		} else {
			c = kseq_read(g_kseq);
			if (c >= 0) {
				g_str.l = 0;
				g_str.m = g_kseq->name.l + g_kseq->comment.l + g_kseq->seq.l + g_kseq->qual.l + 4;
				kroundup32(g_str.m);
				g_str.s = (char*)realloc(g_str.s, g_str.m);
				for (i = 0; i < g_kseq->name.l; ++i)
					g_str.s[g_str.l++] = g_kseq->name.s[i];
				g_str.s[g_str.l++] = '\t';
				for (i = 0; i < g_kseq->seq.l; ++i)
					g_str.s[g_str.l++] = g_kseq->seq.s[i];
				g_str.s[g_str.l++] = '\t';
				for (i = 0; i < g_kseq->qual.l; ++i)
					g_str.s[g_str.l++] = g_kseq->qual.s[i];
				g_str.s[g_str.l++] = '\t';
				for (i = 0; i < g_kseq->comment.l; ++i)
					g_str.s[g_str.l++] = g_kseq->comment.s[i];
				g_str.s[g_str.l++] = '\0';
			} else {
				g_str.l = 0;
				if (g_str.s) g_str.s[0] = '\0';
			}
		}
		adjbuf(&buf, &bufsize, g_str.l + 1, recsize, 0, "bio_getrec");
		memcpy(buf, g_str.s, g_str.l + 1);
		if (c >= 0) {	/* normal record */
			if (isrecord) {
				if (freeable(fldtab[0]))
					xfree(fldtab[0]->sval);
				fldtab[0]->sval = buf;	/* buf == record */
				fldtab[0]->tval = REC | STR | DONTFREE;
				if (is_number(fldtab[0]->sval)) {
					fldtab[0]->fval = atof(fldtab[0]->sval);
					fldtab[0]->tval |= NUM;
				}
			}
			setfval(nrloc, nrloc->fval+1);
			setfval(fnrloc, fnrloc->fval+1);
			*pbuf = buf;
			*psize = bufsize;
			return 1;
		}
		/* EOF arrived on this file; set up next */
		if (!g_is_stdin) {
			kseq_destroy(g_kseq);
			gzclose(g_fp);
		}
		g_fp = 0; g_kseq = 0; g_is_stdin = 0;
		argno++;
	}
	buf[0] = saveb0;
	*pbuf = buf;
	*psize = savesize;
	return 0;	/* true end of file */
}
