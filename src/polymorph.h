/* MKT 
 * For the codons encoded here, they are in order of "ACGT"
 * e.g., AAA,AAC,AAG,AAT,ACA,ACC,ACG,ACT,AGA....
 *
 * Based on this order, we encode codon_tbl to a string,
 * and ns_tbl in a 4096*2 array 
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define DIM 3
#define NUMCODON  64
#define codon_tbl "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLF"

#define gap '-'
#define gaps "---"
#define idx_gap 64

#define MINLEN 500
#define MAXLEN 8192
#define CN 5555 // maximum populations
#define STRLEN 180000 // maximum cDNA length (longest protein 34351*3<105000)
#define MIN_ALPHA -99.999
#define MIN_CUTOFF 0.00000001

#define max(x,y) (((x) >= (y)) ? (x) : (y))
#define min(x,y) (((x) <= (y)) ? (x) : (y))

/* polymorphism of each sequence */
typedef struct
{
    /* polymorph nonsyn */
    int *pn;  

    /* fixed nonsyn */
    int *dn;

    /* polymorph syn */
    int *ps;

    /* fixed syn */
    int *ds;

    /* corresponding allele frequencies */
    int *freq;

    /* total polymorphism  */
    int Pn,Dn,Ps,Ds;

    /* length of the sequence */
    int ncodon; 

    /* population size */
    int npop; 

    /* name of the sequence */
    char *name;

} Polym;

typedef struct
{
    float *x;
    float *ax;
    int num;
    char *name;
} SFS;

typedef char CODON[4];
/*https://stackoverflow.com/questions/33587981/how-to-create-array-of-fixed-length-strings-in-c 
 According to statckoverflow, if DIM is 3, we need to declare CODON to be length of 4, so that the compiler will add '\0' to CODON[3];
*/

/* the two table are from MK.pl, in order of ACGT */
extern char NT[];
extern int n_tbl[];
extern int s_tbl[];

/* compute polymorphism of a sequence from fasta alignments */
//int polymorph_site(CODON *pop_codon,CODON *out_codon, int npop, int nout, float cutoff, int *ns);
//Polym polymorph_fasta(char *fasta_nm,char *ingrp_nm, char **outgrp_nm, int npop, int nout, float cutoff);
//void polymorph_fasta(char *fasta_nm,char *ingrp_nm, char **outgrp_nm, int npop, int nout, float cutoff, Polym p);
void cMK_from_fasta(char *fasta_nm,char *name,char **sequence,char **headers,char *ingrp_nm,char **outgrp_nm,int nout,float cutoff,int *MKpolym,float *MKalpha,bool debug);
