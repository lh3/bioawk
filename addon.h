#ifndef ADDON_H
#define ADDON_H

/* column definition, "header" reads frist line  
   other current options include:
     bed, bedgraph, sam, vcf
*/
extern const char *bio_col_defn;   

void bio_set_colnm(void);             /* assoc. col names w/ numbers */ 
int bio_isvalid_coldef(const char *); /* is the req. col defn. valid? */
void bio_print_valid_coldefs(void);   /* print a list of supported col defns. */

#endif
