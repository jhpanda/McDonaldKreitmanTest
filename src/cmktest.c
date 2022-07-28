/*
 * Read in Fasta
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <getopt.h>
#include <time.h>
#include "fasta.h"
#include "polymorph.h"
#include "fisher_exact.h"

static void usage()
{
    fprintf(stderr,"\n"
"***********************************************************\n"
"             Standard McDonald–Kreitman test\n"
"                   Single thread mode\n"
"***********************************************************\n"
"Reference:\n"
"   1. McDonald, J. H. Kreitman (1991).\n"
"   2. Begun et al 2007.\n"
"   3. Nucleic Acids Res. 2008, W157-62, 10.1093/nar/gkn337.\n"
"\n"
"Usage:\n"
"cmktest_list\n"
"   -f, --fasta      input MSA including outgroup species and ingroup populations.\n"
"   -i, --ingrp      ingroup species (e.g. 'dmel').\n"
"   -o, --outgrp     outgroup species (e.g. 'dsim,dyak').\n"
"   -c, --cutoff     optional, default 0.00000001.\n"
"   -v, --verbose    optional, print polymorphism of each site.\n"
"   -h, --help       optional, print this info.\n"
"\n");
}


int main(int argc, char *argv[])
{

    /*
    char *fasta_nm   = argv[1];
    */

    clock_t t1,t2;
    float time_used;

    int i,ch,verbose_flag=0;
    char *token;
    char fasta_nm[MINLEN];
    char ingrp_nm[MINLEN];
    char outnm[MINLEN];
    char outnm2[MINLEN];
    float cutoff = MIN_CUTOFF;

    static struct option longopts[] =
    {
        {"fasta", required_argument, NULL, 'f'},
        {"ingrp", required_argument, NULL, 'i'},
        {"outgrp", required_argument, NULL, 'o'},
        {"cutoff", optional_argument, NULL, 'c'},
        {"verbose", no_argument, NULL, 'v'},
        {"help", no_argument, NULL, 'h'},
    };

    
    while ((ch = getopt_long(argc, argv, "f:i:o:c:vh", longopts, NULL)) != -1) {
        switch (ch) {
            case 'f':
                strcpy(fasta_nm,optarg);
                break;
            case 'i':
                strcpy(ingrp_nm,optarg);
                break;
            case 'o':
                strcpy(outnm,optarg);
                strcpy(outnm2,optarg);
                break;
            case 'c':
                cutoff = atof(optarg);
                break;
            case 'v':
                verbose_flag = 1;
                break;
            case 'h':
                usage();
                exit(1);
            default:
                usage();
                exit(1);
        }
    }

    if (argc==1)
    {
        fprintf(stderr,"ERROR: Required arguments missing!\n");
        usage();
        exit(1);
    }

    bool debug = false;
    if (verbose_flag==1)
    {
        debug = true;
    }    

    char name[MINLEN];
    basename(fasta_nm,name);
    fprintf(stderr,"%s!\n",fasta_nm);
    fprintf(stderr,"%s!\n",name);

    /* extract outgrp_nm and nout from input */
    int nout = 0;
    
    token = strtok(outnm,",");
    while (token!=NULL)
    {
        nout ++;
        token = strtok(NULL, ",");
    }
    
    char **outgrp_nm;
    snew(outgrp_nm,nout);
    for (i=0;i<nout;i++)
    {
        snew(outgrp_nm[i],MINLEN);
    }

    i = 0;
    token = strtok(outnm2,",");
    while (token!=NULL)
    {
        strcpy(outgrp_nm[i],token);
        i ++;
        token = strtok(NULL, ",");
    }
    /*
    char *ingrp_nm      = "hg19";
    char *outgrp_nm[]   = {"panTro6"};
    */

    int  MKpolym[6];
    float MKalpha[2];

    char **sequence,**headers;
    snew(sequence,CN);
    snew(headers,CN);

    for (i=0;i<CN;i++)
    {
        snew(sequence[i],STRLEN);
        snew(headers[i],MINLEN);
    }
    /*
    t1 = clock();
    */
    cMK_from_fasta(fasta_nm,name,sequence,headers,ingrp_nm,outgrp_nm,nout,cutoff,MKpolym,MKalpha,debug);
    /*
    */
    printf("\nName,Ncodon,Npopulation,Pn,Dn,Ps,Ds,Alpha,Pvalue\n");

    printf("%s,%d,%d,%d,%d,%d,%d,%f,%f\n",name,MKpolym[4],MKpolym[5],MKpolym[0],MKpolym[1],MKpolym[2],MKpolym[3],MKalpha[0],MKalpha[1]);
    /*
    t2 = clock();
    time_used = (float) (t2-t1)/CLOCKS_PER_SEC;
    printf("\ncMK_from_fasta finished in %.3fs\n",time_used);
    */
}
