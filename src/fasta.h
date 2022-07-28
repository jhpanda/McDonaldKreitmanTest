
#define MINLEN 500
#define MAXLEN 8192
#define CN 5555 // maximum populations
#define STRLEN 180000 // maximum cDNA length (longest protein 34351*3<105000)
#define snew(ptr,nelem) (ptr)=save_calloc(#ptr,(nelem),sizeof(*(ptr)))

typedef struct
{
    int  nseq;
    int  *length;
    int  *realLen;
    int  **index;
    char **headers;
    char **sequence;
    char *name;
} SEQ;

typedef struct
{
    char **fasta_list;
    long nfasta;
} LIST;

/* file io and memory options */
void *save_calloc(char *name,unsigned nelem,unsigned elsize);
int empty_line(char *str);
char *my_strdup(const char *src);

/* read in fasta */
void basename(char *nm, char *base_nm);
SEQ read_fasta(char *fasta_nm);
LIST read_fasta_list(char *fasta_list_nm);
