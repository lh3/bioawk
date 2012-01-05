#ifndef ADDON_H
#define ADDON_H

//extern int lh3_has_colnm;

/* column definition, "header" reads frist line  
   other current options include:
     bed, bedgraph, sam, vcf
*/
extern const char *lh3_col_defn;   

void lh3_set_colnm();             /* assoc. col names w/ numbers */ 
int isvalid_coldef(const char *); /* is the req. col defn. valid? */
void print_valid_coldefs();       /* print a list of supported col defns. */

#endif
