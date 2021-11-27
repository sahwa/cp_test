#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <../include/hts.h>
#include <../include/vcf.h>
#include <../include/synced_bcf_reader.h>

#include "zlib.h"

#define PI 3.141593

#define WHERE do {fprintf(stderr,"[%s:%d]",__FILE__,__LINE__);} while(0)
#define WARNING(...) do { fputs("[WARNING]",stderr);WHERE;fprintf(stderr,__VA_ARGS__);fputc('\n',stderr);} while(0)
#define ERROR(...) do { fputs("[ERROR]",stderr);WHERE;fprintf(stderr,__VA_ARGS__);fputc('\n',stderr);abort();} while(0)


int reading(char ** st, char *format, void * res);

struct data_t {
  int nhaps;
  int ndonorhaps;
  int nsnps;
  double * positions;
  double * lambda;
  double * copy_prob;
  double * copy_probSTART;
  double * MutProb_vec;
  int * pop_vec;
  int * hap_label_vec;
  int ** cond_chromosomes;
  int ** ind_chromosomes;
  double * copy_prob_new;
  double * copy_prob_newSTART;
  double * MutProb_vec_new;
  double ** back_prob;
  int * ndonorhaps_vec;
};

struct data_t * ReadData(FILE * fd, int ind_val, int ploidy, int * include_ind_vec, char ** pop_label_vec, int num_donor_pops, char ** donor_pop_vec, int * pop_vec_tot, double * copy_prob_tot, double * copy_probSTART_tot, double * MutProb_vec_tot, int * ndonorhaps_tot, int all_versus_all_ind) {
  //int line_max=100000000;
  int line_max;
  struct data_t * dat;
  char * firstline = malloc(1000 * sizeof(char));
  char * step;
  char waste[400];
  int i, j, k;
  int nind, num_cond_haps;
  int cond_hap_count, ind_hap_count, include_hap;

  dat = malloc(sizeof(struct data_t));
  /* if dat==NULL etc ... */

  /* Number of haplotypes */

  if (fgets(firstline, sizeof firstline, fd) == NULL) {
    printf("error with PHASE-style input file\n");
    exit(1);
  }
  
  sscanf(firstline, "%d", & dat -> nhaps);
  nind = (int) dat -> nhaps / ploidy;
  if ((dat -> nhaps <= 0) || (((int) dat -> nhaps) != dat -> nhaps)) {
    printf("Number of total haplotypes must be an integer value and > 0. Exiting...\n");
    exit(1);
  }

  /* Number of SNPs */
  if (fgets(firstline, sizeof firstline, fd) == NULL) {
    printf("error with PHASE-style input file\n");
    exit(1);
  }

  sscanf(firstline, "%d", & dat -> nsnps);
  if ((dat -> nsnps <= 0) || (((int) dat -> nsnps) != dat -> nsnps)) {
    printf("Number of sites must be an integer value and > 0. Exiting...\n");
    exit(1);
  }
  
	line_max = dat -> nsnps * 8 * 2;
  char * line = malloc(line_max * sizeof(char));

  dat -> ndonorhaps_vec = malloc(num_donor_pops * sizeof(int));
  for (k = 0; k < num_donor_pops; k++)
    dat -> ndonorhaps_vec[k] = 0;
  num_cond_haps = 0;
  
  for (i = 0; i < nind; i++) {
    for (k = 0; k < num_donor_pops; k++) {
      if (include_ind_vec[i] != 0 && i != ind_val && all_versus_all_ind == 0 && strcmp(donor_pop_vec[k], pop_label_vec[i]) == 0) {
        num_cond_haps = num_cond_haps + ploidy;
        dat -> ndonorhaps_vec[k] = dat -> ndonorhaps_vec[k] + ploidy;
        break;
      }
      if (include_ind_vec[i] != 0 && i != ind_val && all_versus_all_ind == 1) {
        num_cond_haps = num_cond_haps + ploidy;
        dat -> ndonorhaps_vec[k] = dat -> ndonorhaps_vec[k] + ploidy;
        break;
      }
    }
  }

  dat -> ndonorhaps = num_cond_haps;
  dat -> positions = malloc(dat -> nsnps * sizeof(double));
  dat -> lambda = malloc((dat -> nsnps - 1) * sizeof(double));
  dat -> cond_chromosomes = malloc(dat -> ndonorhaps * sizeof(int * ));
  dat -> ind_chromosomes = malloc(ploidy * sizeof(int * ));
  dat -> copy_prob = malloc(dat -> ndonorhaps * sizeof(double));
  dat -> copy_probSTART = malloc(dat -> ndonorhaps * sizeof(double));
  dat -> MutProb_vec = malloc(dat -> ndonorhaps * sizeof(double));
  dat -> pop_vec = malloc(dat -> ndonorhaps * sizeof(int));
  dat -> hap_label_vec = malloc(dat -> ndonorhaps * sizeof(int));
  
  for (i = 0; i < dat -> ndonorhaps; i++) {
    dat -> cond_chromosomes[i] = malloc(dat -> nsnps * sizeof(int));
  }
  
  for (i = 0; i < ploidy; i++) {
    dat -> ind_chromosomes[i] = malloc(dat -> nsnps * sizeof(int));
  }
  
  dat -> copy_prob_new = malloc(dat -> ndonorhaps * sizeof(double));
  dat -> copy_prob_newSTART = malloc(dat -> ndonorhaps * sizeof(double));
  dat -> MutProb_vec_new = malloc(dat -> ndonorhaps * sizeof(double));
  dat -> back_prob = malloc(8 * sizeof(double * ));
  
  for (i = 0; i < 8; i++) {
    dat -> back_prob[i] = malloc((dat -> ndonorhaps + 1) * sizeof(double));
  }

  /* Positions */

  if (fgets(line, line_max, fd) == NULL) {
    printf("error with PHASE-style input file\n");
    exit(1);
  }
  
  step = line;
  
  reading( & step, "%s", waste);
  
  for (i = 0; i < dat -> nsnps; i++) {
    reading( & step, "%lf", & dat -> positions[i]);
    if (dat -> positions[i] < 0) {
      printf("Basepair positions must be >= 0. Exiting...\n");
      exit(1);
    }
    if (i < (dat -> nsnps - 1)) dat -> lambda[i] = 1.0;
  }

  cond_hap_count = 0;
  ind_hap_count = 0;

  for (i = 0; i < dat -> nhaps; i++) {

    if ((fgets(line, line_max, fd) == NULL) || (strlen(line) != (dat -> nsnps + 1))) {
      printf("error with PHASE-style input file\n");
      exit(1);
    }
    
    include_hap = 0;

    for (k = 0; k < num_donor_pops; k++) {

      if (include_ind_vec[((int) floor(i / ploidy))] != 0 && (((int) floor(i / ploidy)) != ind_val) && all_versus_all_ind == 0 && strcmp(donor_pop_vec[k], pop_label_vec[((int) floor(i / ploidy))]) == 0) {
        include_hap = 1;
        break;
      }

      if (include_ind_vec[((int) floor(i / ploidy))] != 0 && (((int) floor(i / ploidy)) != ind_val) && all_versus_all_ind == 1) {
        include_hap = 1;
        break;
      }
    }
    
    if (include_hap == 1) {

      dat -> hap_label_vec[cond_hap_count] = i + 1;
      dat -> copy_prob[cond_hap_count] = copy_prob_tot[i];
      dat -> copy_probSTART[cond_hap_count] = copy_probSTART_tot[i];
      dat -> MutProb_vec[cond_hap_count] = MutProb_vec_tot[i];
      
      if (i < (ploidy * ind_val)) dat -> pop_vec[cond_hap_count] = pop_vec_tot[i];
      if (i > (ploidy * ind_val)) dat -> pop_vec[cond_hap_count] = pop_vec_tot[i] - all_versus_all_ind;

      for (j = 0; j < dat -> nsnps; j++) {

        if (line[j] == '0')
          dat -> cond_chromosomes[cond_hap_count][j] = 0;
        if (line[j] == '1')
          dat -> cond_chromosomes[cond_hap_count][j] = 1;
        if (line[j] == 'A')
          dat -> cond_chromosomes[cond_hap_count][j] = 2;
        if (line[j] == 'C')
          dat -> cond_chromosomes[cond_hap_count][j] = 3;
        if (line[j] == 'G')
          dat -> cond_chromosomes[cond_hap_count][j] = 4;
        if (line[j] == 'T')
          dat -> cond_chromosomes[cond_hap_count][j] = 5;
        if (line[j] == '?')
          dat -> cond_chromosomes[cond_hap_count][j] = 9;

        if ((line[j] != '0') && (line[j] != '1') && (line[j] != 'A') && (line[j] != 'G') && (line[j] != 'C') && (line[j] != 'T') && (line[j] != '?')) {
          printf("Allele-type invalid for hap%d, snp%d. Exiting...\n", i + 1, j + 1);
          exit(1);
        }
      }

      cond_hap_count ++;
    
    }



    if (((int) floor(i / ploidy)) == ind_val) {

      for (j = 0; j < dat -> nsnps; j++) {
        if (line[j] == '0')
          dat -> ind_chromosomes[ind_hap_count][j] = 0;
        if (line[j] == '1')
          dat -> ind_chromosomes[ind_hap_count][j] = 1;
        if (line[j] == 'A')
          dat -> ind_chromosomes[ind_hap_count][j] = 2;
        if (line[j] == 'C')
          dat -> ind_chromosomes[ind_hap_count][j] = 3;
        if (line[j] == 'G')
          dat -> ind_chromosomes[ind_hap_count][j] = 4;
        if (line[j] == 'T')
          dat -> ind_chromosomes[ind_hap_count][j] = 5;
        if (line[j] == '?')
          dat -> ind_chromosomes[ind_hap_count][j] = 9;
        if ((line[j] != '0') && (line[j] != '1') && (line[j] != 'A') && (line[j] != 'G') && (line[j] != 'C') && (line[j] != 'T') && (line[j] != '?')) {
          printf("Allele-type invalid for hap%d, snp%d. Exiting...\n", i + 1, j + 1);
          exit(1);
        }
      }
      ind_hap_count = ind_hap_count + 1;
    }
  }


  // #################################################### reading in bcf ######################################## //  

  // htsFile *in = bcf_open("test","r");

  // if (in == NULL) {
  //   ERROR("Cannot open input vcf %s.\n", strerror(errno));
  //   //return EXIT_FAILURE;
  // }
  
  // bcf_hdr_t *header = bcf_hdr_read(in);

  // if (header == NULL) {
  //   ERROR("Cannot open input header %s.\n", strerror(errno));
  //   //return EXIT_FAILURE;
  // }

  // hts_idx_t *index = bcf_index_load("test.index");

  // if (index == NULL) {
  //   ERROR("Cannot open input index %s.\n", strerror(errno));
  //   //return EXIT_FAILURE;
  // }

  // bcf1_t* bcf_line = bcf_init();

  // int ngt, nsmpl;

  // nsmpl = bcf_hdr_nsamples(header);
  // //int nsnps = 2;

  // printf("Processing %d samples\n", nsmpl);

  // int32_t *gt_arr = NULL, ngt_arr = 0;

  // // int (*al_arr_2d)[nsmpl] = malloc(sizeof(int[(nsmpl * 2)-1][nsnps-1]));

  // while (bcf_read(in, header, bcf_line)==0) {

  //   ngt = bcf_get_genotypes(header, bcf_line, &gt_arr, &ngt_arr);

  //   int bcf_ploidy = ngt/nsmpl;

  //   if (bcf_ploidy != ploidy) {
  //     ERROR("Difference between CP ploidy and bcf ploidy. Exiting.\n");
  //     //return EXIT_FAILURE;  
  //   }

  //   for (int i=0; i < nsmpl; i++) {

  //     int32_t *ptr = gt_arr + i*2;

  //     for (int j=0; j<bcf_ploidy; j++) {

  //       dat -> cond_chromosomes[cond_hap_count][j] = bcf_gt_allele(ptr[j]);
  //       // al_arr_2d[1][1] = a0;

  //     }

  //   }
  // }

  // printf("Test allele is %d\n", al_arr_2d[1][1]);

  // bcf_hdr_destroy(header);
  // bcf_destroy(bcf_line);
  // bcf_close(in);

  // #################################################### reading in bcf ######################################## //  
  
  if (ind_hap_count < ploidy) {
    printf("Could not find haplotype(s) of individual %d in genotype input file! Row missing?? Exiting...\n", ind_val + 1);
    exit(1);
  }
  
  if (cond_hap_count != dat -> ndonorhaps) {
    printf("Something wrong with genotype input file. Exiting...\n");
    exit(1);
  }

  free(firstline);
  free(line);
  return dat;

}
