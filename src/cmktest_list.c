/*
 * Read in Fasta
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>
#include <getopt.h>
#include "fasta.h"
#include "polymorph.h"
#include "fisher_exact.h"


static void usage()
{
    fprintf(stderr,"\n"
"***********************************************************\n"
"             Standard McDonald–Kreitman test\n"
"                  Multi-threading mode\n"
"***********************************************************\n"
"Reference:\n"
"   1. McDonald, J. H. Kreitman (1991).\n"
"   2. Begun et al 2007.\n"
"   3. Nucleic Acids Res. 2008, W157-62, 10.1093/nar/gkn337.\n"
"\n"
"Usage:\n"
"cmktest_list\n"
"   -l, --list       a list of MSA in fasta format.\n"
"   -i, --ingrp      ingroup species (e.g. 'dmel').\n"
"   -o, --outgrp     outgroup species (e.g. 'dsim,dyak').\n"
"   -c, --cutoff     optional, default 0.00000001.\n"
"   -n, --ncpu       optional, default 4.\n"
"   -v, --verbose    optional, print polymorphism of each site.\n"
"   -h, --help       optional, print this info.\n"
"\n");
}

int main(int argc, char *argv[])
{

    /*
    if (argc!=7)
    {
        fprintf(stderr,"Usage: cmktest_list <fasta_list> <ingrp_nm> <outgrp_nm seperated by comma [out1,out2]> <cutoff> <debug,nodebug> <ncpu>\n");
        exit(1);
    }
    char *fasta_nm   = argv[1];
    */

    int  ch,verbose_flag = 0;
    long i;
    char fasta_list_nm[MINLEN];
    char ingrp_nm[MINLEN];
    char outnm[MINLEN];
    char outnm2[MINLEN];
    char *token;
    float cutoff = MIN_CUTOFF;
    int  ncpu = 4;

    static struct option longopts[] =
    {
        {"list", required_argument, NULL, 'l'},
        {"ingrp", required_argument, NULL, 'i'},
        {"outgrp", required_argument, NULL, 'o'},
        {"cutoff", optional_argument, NULL, 'c'},
        {"ncpu", optional_argument, NULL, 'n'},
        {"verbose", no_argument, NULL, 'v'},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    
    while ((ch = getopt_long(argc, argv, "hvn:c:o:i:l:", longopts, NULL)) != -1) {
        switch (ch) {
            case 'l':
                strcpy(fasta_list_nm,optarg);
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
            case 'n':
                ncpu = atoi(optarg);
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

    char  output[MINLEN];
    sprintf(output,"%s_mktest_cutoff_%.2f.out",ingrp_nm,cutoff);

    /* extract outgrp_nm and nout from input */
    long nout = 0;

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
    
    fprintf(stderr,"Read in list: %s\n",fasta_list_nm);

    LIST list = read_fasta_list(fasta_list_nm);
    char name[MINLEN];
    int  MKpolym[6];
    float MKalpha[2];

    int *Pn,*Dn,*Ps,*Ds,*length,*pop;
    char **sequence,**headers;
    float *alpha,*pvalue;

    snew(Pn,list.nfasta);
    snew(Ps,list.nfasta);
    snew(Ds,list.nfasta);
    snew(Dn,list.nfasta);
    snew(length,list.nfasta);
    snew(pop,list.nfasta);
    snew(alpha,list.nfasta);
    snew(pvalue,list.nfasta);


    omp_set_dynamic(0);
    fprintf(stderr,"Initializing OMP...\n");
    //#pragma omp parallel shared(Pn,Dn,Ps,Ds,length,alpha,pvalue,list,ingrp_nm,outgrp_nm,nout,cutoff,debug) private(i,MKpolym,MKalpha,ncodon,npop,name,sequence,headers) num_threads(8)
    #pragma omp parallel private(i,MKpolym,MKalpha,name,sequence,headers) num_threads(ncpu)
    {
        snew(sequence,CN);
        snew(headers,CN);
 
        for (i=0;i<CN;i++)
        {
            snew(sequence[i],STRLEN);
            snew(headers[i],MINLEN);
        }
        #pragma omp for nowait
        for (i=0;i<list.nfasta;i++)
        //for (i=0;i<4;i++)
        {
            basename(list.fasta_list[i],name);
            cMK_from_fasta(list.fasta_list[i],name,sequence,headers,ingrp_nm,outgrp_nm,nout,cutoff,MKpolym,MKalpha,debug);
            Pn[i]     = MKpolym[0];
            Dn[i]     = MKpolym[1];
            Ps[i]     = MKpolym[2];
            Ds[i]     = MKpolym[3];
            length[i] = MKpolym[4];
            pop[i]    = MKpolym[5];
            alpha[i]  = MKalpha[0];
            pvalue[i] = MKalpha[1];
            printf("%d,%s,%d,%d,%d,%d,%d,%d,%.3f,%.3f\n",i,name,MKpolym[4],MKpolym[5],MKpolym[0],MKpolym[1],MKpolym[2],MKpolym[3],MKalpha[0],MKalpha[1]);
        }
    }

    fprintf(stderr,"Finished, now write to file %s\n",output);
    FILE *fp = fopen(output,"w");
    fprintf(fp,"Name,Ncodon,Npopulation,Pn,Dn,Ps,Ds,Alpha,Pvalue\n");
    for (i=0;i<list.nfasta;i++)
    {
        basename(list.fasta_list[i],name);
        fprintf(fp,"%s,%d,%d,%d,%d,%d,%d,%.3f,%.3lf\n",name,length[i],pop[i],Pn[i],Dn[i],Ps[i],Ds[i],alpha[i],pvalue[i]);
    }
    close(fp);

}
