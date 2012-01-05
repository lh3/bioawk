#include <ctype.h>
#include <stdio.h>
#include "awk.h"
#include "addon.h"

int lh3_has_colnm = 0;

static void set_colnm_aux(const char *p, int col)
{
	const char *q;
	Cell *x;
	for (q = p; *q; ++q)
		if (!isdigit(*q)) break;
	if (*q == 0) return; /* do not set if string p is an integer */
	if ((x = lookup(p, symtab)) != NULL)
		x->tval = NUM, x->fval = col;
}

void lh3_set_colnm()
{
	char *p, *q, c;
	int i;
	if (lh3_has_colnm == 0) return;
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
}
