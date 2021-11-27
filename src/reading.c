#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <../include/hts.h>
#include <../include/vcf.h>
#include <../include/synced_bcf_reader.h>


int reading(st, format, res) 
  char ** st, * format;
  void * res; {
  int i;
  char * rs;
  rs = * st;
  for (i = 0; isspace(rs[i]); i++);
  if (!rs[i]) return 0;
  for (; !isspace(rs[i]); i++);
  if (rs[i]) rs[i++] = 0;
  if (!sscanf( * st, format, res)) return 0;
  * st += i;
  return 1;
}

#define WHERE do { fprintf(stderr,"[%s:%d]",__FILE__,__LINE__);} while(0)
#define WARNING(...) do { fputs("[WARNING]",stderr);WHERE;fprintf(stderr,__VA_ARGS__);fputc('\n',stderr);} while(0)
#define ERROR(...) do { fputs("[ERROR]",stderr);WHERE;fprintf(stderr,__VA_ARGS__);fputc('\n',stderr);abort();} while(0)

int reading_htslib(int argc, char const *argv[]) {

  htsFile *in = bcf_open(argc!=2 ?"-":argv[1],"r");

  if (in == NULL) {
    ERROR("Cannot open input vcf %s.\n", strerror(errno));
    return EXIT_FAILURE;
  }
  
  bcf_hdr_t *header = bcf_hdr_read(in);

  if (header == NULL) {
    ERROR("Cannot open input header %s.\n", strerror(errno));
    return EXIT_FAILURE;
  }

  hts_idx_t *index = bcf_index_load(argc!=2 ?"-":argv[1]);

  if (index == NULL) {
    ERROR("Cannot open input index %s.\n", strerror(errno));
    return EXIT_FAILURE;
  }

  bcf1_t* line = bcf_init();

  int ngt, nsmpl;

  nsmpl = bcf_hdr_nsamples(header);
  int nsnps = 2;

  printf("Processing %d samples\n", nsmpl);

  int32_t *gt_arr = NULL, ngt_arr = 0;

  int (*al_arr_2d)[nsmpl] = malloc(sizeof(int[(nsmpl * 2)-1][nsnps-1]));

  while (bcf_read(in, header, line)==0) {

    ngt = bcf_get_genotypes(header, line, &gt_arr, &ngt_arr);

    int ploidy = ngt/nsmpl;

    if (ploidy != 2) {
      ERROR("Ploidy identified as %d. Exiting.", ploidy);
      return EXIT_FAILURE;  
    }

    for (int i=0; i < nsmpl; i++) {

      int32_t *ptr = gt_arr + i*2;

      for (int j=0; j<=1; j++) {

        int a0 = bcf_gt_allele(ptr[j]);
        al_arr_2d[1][1] = a0;

        printf("Allele 0 is: %d\n", a0);
      }

    }
  }

  printf("Test allele is %d\n", al_arr_2d[1][1]);

  bcf_hdr_destroy(header);
  bcf_destroy(line);
  bcf_close(in);

  return 0;
}
