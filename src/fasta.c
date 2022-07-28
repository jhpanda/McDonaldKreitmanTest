/*
 * Read in Fasta
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fasta.h"

int empty_line(char *str)
{
    while (*str) {
        if (!isspace(*str++)) {
            return 0;
        }
    }
    return 1;
}

void *save_calloc(char *name,unsigned nelem,unsigned elsize)
{
  void *p;

  p = NULL;
  if ((nelem == 0) || (elsize == 0))
    p = NULL;
  else
    p = malloc((size_t)nelem * (size_t)elsize);
  return p;
}

char *my_strdup(const char *src)
{
    size_t len = strlen(src) + 1;       // String plus '\0'
    char *dst = malloc(len);            // Allocate space
    if (dst == NULL) return NULL;       // No memory
    memcpy (dst, src, len);             // Copy the block
    return dst;                         // Return the new string
}

void basename(char *nm, char *base_nm)
{
    int i = 0;
    int m = 0;
    int idxf = 0;
    int idxe = 0;

    /* find out index for last '/' and last '.' */
    i = 0;
    while (nm[i]!='\0')
    {
        if (nm[i]=='/')
        {
            idxf = i;
        }
        if (nm[i]=='.')
        {
            idxe = i;
        }
        i ++;
    }

    if (idxe == 0 && idxf !=0 )
    {
        i = idxf;
        m = 0;
        while (nm[i]!='\0')
        {
            base_nm[m] = nm[i];
            m ++;
            i ++;
        }
        base_nm[m] = '\0';
    }

    else if (idxe !=0 && idxf !=0 )
    {
        m = 0;
        for (i=idxf+1;i<idxe;i++)
        {
            base_nm[m] = nm[i];
            m ++;
        }
        base_nm[m] = '\0';
    }

    else
    {
        m = 0;
        for (i=idxf;i<idxe;i++)
        {
            base_nm[m] = nm[i];
            m ++;
        }
        base_nm[m] = '\0';
    }
}

SEQ read_fasta(char *fasta_nm) 
{
    size_t len;
    int  i,read,idx,iseq,seqlen;
    char identifier;
    char header[MINLEN];
    char *str = NULL;
    char *headers[MINLEN];

    FILE *fp;

    SEQ s;

    char base_nm[MAXLEN];

    basename(fasta_nm,base_nm);
    snew(s.name,MAXLEN);
    strcpy(s.name,base_nm);
    /*
    fprintf(stderr,"%s %s\n",fasta_nm,s.name);
    */
    s.nseq = 0;
    fp = fopen(fasta_nm,"r");
    while ((read = getline(&str, &len, fp))!=-1){
        if (!empty_line(str)) 
        {
            identifier = str[0];
            if (identifier == '>') 
            {
                s.nseq ++;
            }
        }
    }
    snew(s.headers,s.nseq);
    snew(s.sequence,s.nseq);
    snew(s.length,s.nseq);
    snew(s.realLen,s.nseq);
    snew(s.index,s.nseq);

     
    /*make sure seq is empty*/
    fp = fopen(fasta_nm,"r");
    header[0] = '\0';
    iseq      = 0;
    while ((read = getline(&str, &len, fp))!=-1){
        if (!empty_line(str)) 
        {
            identifier = str[0];
            if (identifier == '>') 
            {
 
                /* a new header means to summarize previous sequence */
                if (iseq>0) 
                {
                    s.sequence[iseq-1][idx] = '\0';// make sure seq ends
                    s.length[iseq-1] = idx;
                    s.realLen[iseq-1] = idx;
                }
 
                /* a new header means a new sequence */
                seqlen = strlen(str);
                snew(s.headers[iseq],seqlen);
                for (i=0;i<seqlen-2;i++)
                {
                    s.headers[iseq][i] = str[i+1];
                }
                s.headers[iseq][i] = '\0';
 
                idx = 0; // initialize the index of each sequence
                snew(s.sequence[iseq],STRLEN);
                iseq ++;
 
            }
            else if (identifier != '>')
            {
                seqlen = strlen(str);
                for (i=0;i<seqlen-1;i++)
                {
                    s.sequence[iseq-1][idx] = str[i];
                    idx ++;
                } 
            }
        }
    }
    /* don't forget the last sequence */
    s.sequence[iseq-1][idx] = '\0';// make sure seq ends
    s.length[iseq-1] = idx;
    s.realLen[iseq-1] = idx;
    
    //fprintf(stderr, "%6d sequences in %s\n",iseq,fasta_nm);
    /* re_index, so that gaps are not counted in new index */
    /*
    int j,m;
    for (i=0;i<s.nseq;i++)
    {
        snew(s.index[i],s.length[i]);
        m = 0;
        for (j=0;j<s.length[i];j++)
        {
            if (s.sequence[i][j]!='-')
            {
                s.index[i][m] = j;
                m ++;
            }
        }
        s.index[i][m] = '\0';
        s.realLen[i] = m;
    }
    */
    /* close file */
    fclose(fp);
    return s;
}

LIST read_fasta_list(char *fasta_list_nm)
{
    FILE *fp;

    size_t len;
    long i,read,flag;
    char *str,*tok;
    LIST list;

    list.nfasta = 0;
    fp = fopen(fasta_list_nm,"r");
    while ((read = getline(&str, &len, fp))!=-1)
    {
        if (!empty_line(str)) 
        {
            flag   = 0;
            len = strlen(str);
            str[len-1] = '\0';
            tok = strrchr(str, '.');
            if (strcmp(tok,".fa")==0 || strcmp(tok,".fasta")==0)
            {
                flag = 1;
            }
            else
            {
                fprintf(stderr,">WARN: not a valid format %s",str);
                exit(1);
            }
            if (flag)
            {
                list.nfasta ++;
            }
        }
    }

    snew(list.fasta_list,list.nfasta);
    for (i=0;i<list.nfasta;i++)
    {
        snew(list.fasta_list[i],MINLEN);
    }

    fp = fopen(fasta_list_nm,"r");
    i = 0;
    while ((read = getline(&str, &len, fp))!=-1)
    {
        if (!empty_line(str)) 
        {
            len = strlen(str);
            strcpy(list.fasta_list[i],str);
            list.fasta_list[i][len-1] = '\0';
            i ++;
        }
    }

    return list;
}
