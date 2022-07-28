#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "fisher_exact.h"

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        fprintf(stderr,"Usage: fisher_exact_test <n11> <n12> <n21> <n22>\n");
        exit(1);

    }
    int n11 = atoi(argv[1]);
    int n12 = atoi(argv[2]);
    int n21 = atoi(argv[3]);
    int n22 = atoi(argv[4]);
    double _left_p,_right_p,pvalue;
    fisher_exact(n11,n12,n21,n22,&_left_p,&_right_p,&pvalue);
    printf("left  tailed Pvalue: %14.9f\n",_left_p);
    printf("right tailed Pvalue: %14.9f\n",_right_p);
    printf("two   tailed Pvalue: %14.9f\n",pvalue);
}
