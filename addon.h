#ifndef ADDON_H
#define ADDON_H

/* do not change these values! */
#define BIO_NULL -1
#define BIO_HDR  0
#define BIO_BED  1
#define BIO_SAM  2
#define BIO_VCF  3
#define BIO_GFF  4

#define BIO_SHOW_HDR 0x1

extern int bio_fmt, bio_flag;
extern char *bio_hdr_chr;

int bio_get_fmt(const char *s);
int bio_skip_hdr(const char *r);
void bio_set_colnm(void);

#endif
