#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "awk.h"
#include "addon.h"

const char *lh3_col_defn = NULL;
const char *valid_coldefs[] = {"header", "bed", "bedgraph", "sam", "vcf", "gff", NULL};

/*BED*/
const char *bed_coldefs[] = {"chrom", "start", "end", 
                             "name", "score", "strand", 
                             "thickstart", "thickend", "rgb",
                             "blockcount", "blocksizes", "blockstarts",
                              NULL};
/*BEDGRAPH*/
const char *bedgraph_coldefs[] = {"chrom", "start", "end", "score", NULL};
/*SAM*/
const char *sam_coldefs[] = {"qname", "flag", "rname",
                             "pos", "mapq", "cigar",
                             "rnext", "pnext", "tlen",
                             "seq", "qual", NULL};
/*VCF*/
const char *vcf_coldefs[] = {"chrom", "pos", "id",
                             "ref", "alt", "qual",
                             "filter" "info", NULL};
/*GFF/GTF*/
const char *gff_coldefs[] = {"seqname", "source", "feature", 
                             "start", "end", "score",
                             "filter", "strand", "group",
                             "attribute", NULL};

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

int isvalid_coldef(const char *request)
{
    int i;
    for (i = 0; valid_coldefs[i] != NULL; ++i)
        if (strcmp(request, valid_coldefs[i]) == 0)
            return 1;
    return 0;
}

void print_valid_coldefs() 
{
    printf("valid -c options include:\n");
    int i;
    for (i = 0; valid_coldefs[i] != NULL; ++i)
    {
        const char *option = valid_coldefs[i];
        printf("  %d. \"%s\"\n", i+1,option);

        if (strcmp(option, "header") == 0)
            printf("    input should contain a col. defn header as first line\n");
        else if (strcmp(option, "bed") == 0) {
            int j;
            for (j = 0; bed_coldefs[j] != NULL; ++j)
                printf("    %s: column $%d\n", bed_coldefs[j], j+1);
        }
        else if (strcmp(option, "bedgraph") == 0) {
            int j;
            for (j = 0; bedgraph_coldefs[j] != NULL; ++j)
                printf("    %s: column $%d\n", bedgraph_coldefs[j], j+1);
        }
        else if (strcmp(option, "sam") == 0) {
            int j;
            for (j = 0; sam_coldefs[j] != NULL; ++j)
                printf("    %s: column $%d\n", sam_coldefs[j], j+1);
        }
        else if (strcmp(option, "vcf") == 0) {
            int j;
            for (j = 0; vcf_coldefs[j] != NULL; ++j)
                printf("    %s: column $%d\n", vcf_coldefs[j], j+1);
        }
        else if (strcmp(option, "gff") == 0) {
            int j;
            for (j = 0; gff_coldefs[j] != NULL; ++j)
                if (strcmp(gff_coldefs[j], "attribute") != 0)
                    printf("    %s: column $%d\n", gff_coldefs[j], j+1);
                else
                    printf("    %s: column $%d\n", gff_coldefs[j], j);
        }
    }
}

void lh3_set_colnm()
{
    if (lh3_col_defn == NULL) return;
    
    if (strcmp(lh3_col_defn, "header") == 0) {
        char *p, *q, c;
        int i;
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
    else if (strcmp(lh3_col_defn, "bed") == 0) {
        int i;
        for (i = 0; bed_coldefs[i] != NULL; ++i)
            set_colnm_aux(bed_coldefs[i], i+1);
        // force tab delimited input and output
        *FS = "\t";
        *OFS = "\t";
    }
    else if (strcmp(lh3_col_defn, "bedgraph") == 0) {
        int i;
        for (i = 0; bedgraph_coldefs[i] != NULL; ++i)
            set_colnm_aux(bedgraph_coldefs[i], i+1);
        // force tab delimited input and output
        *FS = "\t";
        *OFS = "\t";
    }
    else if (strcmp(lh3_col_defn, "sam") == 0) {
        int i;
        for (i = 0; sam_coldefs[i] != NULL; ++i)
            set_colnm_aux(sam_coldefs[i], i+1);
        // force tab delimited input and output
        *FS = "\t";
        *OFS = "\t";
        // auto-report any header lines
        while (getrec(&record, &recsize, 1) > 0 && record[0] == '@') {
            printf("%s\n", record);
        }
    }
    else if (strcmp(lh3_col_defn, "vcf") == 0) {
        int i;
        for (i = 0; vcf_coldefs[i] != NULL; ++i)
            set_colnm_aux(vcf_coldefs[i], i+1);
        // todo: any intelligent way to handle genotypes?
        // force tab delimited input and output
        *FS = "\t";
        *OFS = "\t";
        // auto-report any header lines
        while (getrec(&record, &recsize, 1) > 0 && record[0] == '#') {
            printf("%s\n", record);
        }
    }
    else if (strcmp(lh3_col_defn, "gff") == 0 || strcmp(lh3_col_defn, "gtf") == 0) {
        int i;
        for (i = 0; gff_coldefs[i] != NULL; ++i)
            // allow "group" and "attribute" to be the ninth column
            if (strcmp(lh3_col_defn, "attribute") != 0)
                set_colnm_aux(gff_coldefs[i], i+1);
            else
                set_colnm_aux(gff_coldefs[i], i);
        // force tab delimited input and output
        *FS = "\t";
        *OFS = "\t";
        // auto-report any header lines
        while (getrec(&record, &recsize, 1) > 0 && record[0] == '#') {
            printf("%s\n", record);
        }
    }
}
