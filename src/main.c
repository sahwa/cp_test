#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "zlib.h"
#include "data_t.h"
#include "loglik.h"
#include "destroy.h"
#include "reading.h"
#include "usage.h"

int main(int argc, char * argv[]) {
  struct data_t * Data;
  int i, j, k;
  double bpval;
  double totaldonorprobs;
  int log_lik_check;
  int ndonors = 0;

  int ndonorpops, nrecpops;
  int nind, nsites, cond_nhaps, cond_nind, num_rec_drift;
  int donor_count, rec_count, nind_tot, nind_totGENFILE, rec_ind_count, rec_ind_topaint_count, ind_label_find, ind_pop_find;
  int geno_find, recom_find, donorlist_find, outfile_find, ne_find, mut_find, num_found, copy_prop_em_find, recom_em_find, mutation_em_find, mutationALL_em_find, all_versus_all_ind, haploid_ind, unlinked_ind, prior_donor_probs_ind, mutation_rate_ind, print_file9_ind, drift_calc_ind, idfile_find;
  char * step;
  char line[2047];
  //char * bigline=malloc(1000000*sizeof(char));
  char waste[2047];
  char waste2[2047];
  char templab[15];
  char donorname[2047];
  char indname[2047];
  FILE * fd, * fd2, * fd3, * fd4, * fout, * fout2, * fout3, * fout4, * fout5, * fout6, * fout7, * fout8;
  gzFile fout9;
  char * filename = malloc(1000 * sizeof(char));
  char * filenameGEN = malloc(1000 * sizeof(char));
  char * filenameDONORLIST = malloc(1000 * sizeof(char));
  char * filenameID = malloc(1000 * sizeof(char));
  char * filenameOUT = malloc(1000 * sizeof(char));
  srand((unsigned) time(NULL));

  /***********************************************************/
  // DEFAULT VALUES:

  int EMruns = 0; // number of EMruns 
  int samplesTOT = 10; // number of final hidden-state samples desired after E-M is finished
  double N_e = 400000; // scaling constant for recombination rate
  double GlobalMutRate = -9.0; // global mutation rate per site
  double region_size = 100; // number of chunks per region -- used to look at variability in copied chunks across regions in order to estimate "c" in Dan Lawson's fineSTRUCTURE

  /* 'NUISSANCE' PARAMETER DETAILS: */
  //double theta = 0.0001;
  //double theta = -9.0;
  double small_recom_val = 0.000000000000001; // lower limit for small genmap rates

  int start_val = 0;
  int end_val = 0;

  /************************************************************/

  geno_find = 0;
  recom_find = 0;
  donorlist_find = 0;
  idfile_find = 0;
  outfile_find = 0;
  ne_find = 0;
  mut_find = 0;
  copy_prop_em_find = 0;
  recom_em_find = 0;
  mutation_em_find = 0;
  mutationALL_em_find = 0;
  all_versus_all_ind = 0;
  haploid_ind = 0;
  unlinked_ind = 0;
  prior_donor_probs_ind = 0;
  mutation_rate_ind = 0;
  print_file9_ind = 0;
  drift_calc_ind = 0;
  num_found = 0;

  for (i = 1; i < argc; i++) {
    if ((strcmp(argv[i], "-help") == 0) || (strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-h") == 0)) {
      usage();
      exit(1);
    }
    if (strcmp(argv[i], "-g") == 0) {
      geno_find = 1;
      num_found = num_found + 1;
    }
    if (strcmp(argv[i], "-r") == 0) {
      recom_find = 1;
      num_found = num_found + 1;
    }
    if (strcmp(argv[i], "-f") == 0)
      donorlist_find = 1;
    if (strcmp(argv[i], "-t") == 0) {
      idfile_find = 1;
      num_found = num_found + 1;
    }
    if (strcmp(argv[i], "-n") == 0) {
      ne_find = 1;
      num_found = num_found + 1;
    }
    if (strcmp(argv[i], "-M") == 0) {
      mut_find = 1;
      num_found = num_found + 1;
    }
    if (strcmp(argv[i], "-o") == 0) {
      outfile_find = 1;
      num_found = num_found + 1;
    }
    if (strcmp(argv[i], "-ip") == 0)
      copy_prop_em_find = 1;
    if (strcmp(argv[i], "-in") == 0)
      recom_em_find = 1;
    if (strcmp(argv[i], "-im") == 0)
      mutation_em_find = 1;
    if (strcmp(argv[i], "-iM") == 0)
      mutationALL_em_find = 1;
    //       if (strcmp(argv[i],"-c")==0)
    //   condition_recipient_inds_find=1;
    if (strcmp(argv[i], "-a") == 0)
      all_versus_all_ind = 1;
    if (strcmp(argv[i], "-j") == 0)
      haploid_ind = 1;
    if (strcmp(argv[i], "-u") == 0)
      unlinked_ind = 1;
    if (strcmp(argv[i], "-p") == 0)
      prior_donor_probs_ind = 1;
    if (strcmp(argv[i], "-b") == 0)
      print_file9_ind = 1;
    //if (strcmp(argv[i],"-y")==0)
    //indcount_suppress_ind=1;
    if (strcmp(argv[i], "-m") == 0) {
      mutation_rate_ind = 1;
      //num_found=num_found+1;
    }
  }
  if (argc != (num_found * 2 + copy_prop_em_find + recom_em_find + mutation_em_find + mutationALL_em_find + mutation_rate_ind + 4 * donorlist_find + 3 * all_versus_all_ind + haploid_ind + unlinked_ind + prior_donor_probs_ind + print_file9_ind + 1)) {
    printf("Something wrong with input command line (missing arguments?). Exiting....\n");
    usage();
    exit(1);
  }
  if (donorlist_find == 0)
    strcpy(filenameDONORLIST, "NULL");
  if (recom_find == 0)
    strcpy(filename, "NULL");
  if (idfile_find == 0)
    strcpy(filenameID, "NULL");
  if (((geno_find == 0) || (recom_find == 0)) && (unlinked_ind == 0)) {
    printf("Error with command line (Each of -g and -r MUST be specified if data are linked). Exiting...\n\n");
    usage();
    exit(1);
  }
  if ((geno_find == 0) && (unlinked_ind == 1)) {
    printf("Error with command line (-g MUST be specified). Exiting...\n");
    exit(1);
  }
  if ((recom_find == 1) && (unlinked_ind == 1)) {
    printf("Data specified as containing unlinked sites (-u). Ignoring supplied recombination rate file....\n");
  }
  if (all_versus_all_ind == 0 && (idfile_find == 0 || donorlist_find == 0)) {
    printf("Unless performing all-versus-all ('-a' switch), you MUST specify both the id file ('-t' switch) and file listing donor and recipient populations ('-f' switch). Exiting....\n");
    exit(1);
  }
  if ((mutation_em_find == 1) && (mutationALL_em_find == 1)) {
    printf("You have specified to estimate a global mutation (emission) rate and population-specific mutation (emission) rates. Please choose only one of the '-im' and '-iM' switches. Exiting...\n");
    exit(1);
  }
  if ((mutation_rate_ind == 1) && (mut_find == 1)) {
    printf("You have provided values for both a global mutation (emission) rate ('-M') and population-specific mutation (emission) rates ('-m'). Please choose only one of the '-m' and '-M' switches. Exiting...\n");
    exit(1);
  }
  if ((mutation_rate_ind == 1) && (mutationALL_em_find == 1)) {
    printf("You have specified to estimate a global mutation (emission) rate; will ignore population-specific mutation (emission) rates in %s. If you wish to use donor-specific mutation rates, use the '-im' switch.\n", filenameDONORLIST);
  }

  for (i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-g") == 0) {
      if (argv[(i + 1)][0] == '-') {
        printf("Something wrong with input command line (missing arguments?). Exiting....\n");
        usage();
        exit(1);
      }
      strcpy(filenameGEN, argv[(i + 1)]);
      if (outfile_find == 0) {
        strcpy(filenameOUT, argv[(i + 1)]);
        fout = fopen(strcat(filenameOUT, ".samples.out"), "w");
        strcpy(filenameOUT, argv[(i + 1)]);
        fout2 = fopen(strcat(filenameOUT, ".priorprobs.out"), "w");
        strcpy(filenameOUT, argv[(i + 1)]);
        fout3 = fopen(strcat(filenameOUT, ".EMprobs.out"), "w");
        strcpy(filenameOUT, argv[(i + 1)]);
        fout4 = fopen(strcat(filenameOUT, ".chunkcounts.out"), "w");
        strcpy(filenameOUT, argv[(i + 1)]);
        fout5 = fopen(strcat(filenameOUT, ".chunklengths.out"), "w");
        strcpy(filenameOUT, argv[(i + 1)]);
        fout6 = fopen(strcat(filenameOUT, ".mutationprobs.out"), "w");
        strcpy(filenameOUT, argv[(i + 1)]);
        fout7 = fopen(strcat(filenameOUT, ".regionchunkcounts.out"), "w");
        strcpy(filenameOUT, argv[(i + 1)]);
        fout8 = fopen(strcat(filenameOUT, ".regionsquaredchunkcounts.out"), "w");
        strcpy(filenameOUT, argv[(i + 1)]);
        fout9 = gzopen(strcat(filenameOUT, ".copyprobsperlocus.out.gz"), "w");
        strcpy(filenameOUT, argv[(i + 1)]);
      }
    }
    if (strcmp(argv[i], "-r") == 0) {
      if (argv[(i + 1)][0] == '-') {
        printf("Something wrong with input command line (missing arguments?). Exiting....\n");
        usage();
        exit(1);
      }
      strcpy(filename, argv[(i + 1)]);
    }
    if (strcmp(argv[i], "-f") == 0) {
      if ((argv[(i + 1)][0] == '-') || (argv[(i + 2)][0] == '-') || (argv[(i + 3)][0] == '-')) {
        printf("Something wrong with input command line (missing arguments?). Exiting....\n");
        usage();
        exit(1);
      }
      strcpy(filenameDONORLIST, argv[(i + 1)]);
      start_val = atoi(argv[(i + 2)]);
      end_val = atoi(argv[(i + 3)]);
      if ((end_val < start_val) || (start_val < 0) || (end_val < 0)) {
        printf("Invalid start_ind/stop_ind vals ('-f' switch). If you want to paint each recipient individual using every donor individual, use '-f <donorlist.filein> 0 0'. Exiting...\n");
        exit(1);
      }
      if (start_val > 0) start_val = start_val - 1;
    }
    if (strcmp(argv[i], "-t") == 0) {
      if (argv[(i + 1)][0] == '-') {
        printf("Something wrong with input command line (missing arguments?). Exiting....\n");
        usage();
        exit(1);
      }
      strcpy(filenameID, argv[(i + 1)]);
    }
    if (strcmp(argv[i], "-i") == 0) {
      if (argv[(i + 1)][0] == '-') {
        printf("Something wrong with input command line (missing arguments?). Exiting....\n");
        usage();
        exit(1);
      }
      EMruns = atoi(argv[(i + 1)]);
      if (EMruns < 0) {
        printf("Number of EM runs must be at least 0. Exiting...\n");
        exit(1);
      }
      if ((EMruns > 0) && (copy_prop_em_find == 0) && (recom_em_find == 0) && (mutation_em_find == 0) && (mutationALL_em_find == 0)) {
        printf("You have specified to perform E-M iterations, but have not specified which parameter(s) to maximize. If using '-i' switch, please specify at least one of '-in', '-ip', '-im', and/or '-iM'. Exiting...\n");
        exit(1);
      }
    }
    if (strcmp(argv[i], "-s") == 0) {
      if (argv[(i + 1)][0] == '-') {
        printf("Something wrong with input command line (missing arguments?). Exiting....\n");
        usage();
        exit(1);
      }
      samplesTOT = atoi(argv[(i + 1)]);
      if (samplesTOT < 0) {
        printf("Number of samples must be >= 0. Exiting...\n");
        exit(1);
      }
    }
    if (strcmp(argv[i], "-n") == 0) {
      if (argv[(i + 1)][0] == '-') {
        printf("Something wrong with input command line (missing arguments?). Exiting....\n");
        usage();
        exit(1);
      }
      N_e = atof(argv[(i + 1)]);
      if (N_e <= 0) {
        printf("Recombination scaling parameter N_e must be > 0. Exiting...\n");
        exit(1);
      }
    }
    if (strcmp(argv[i], "-M") == 0) {
      if (argv[(i + 1)][0] == '-') {
        printf("Something wrong with input command line (missing arguments?). Exiting....\n");
        usage();
        exit(1);
      }
      GlobalMutRate = atof(argv[(i + 1)]);
      if (GlobalMutRate == 0) GlobalMutRate = -9;
      if (GlobalMutRate < 0)
        printf("Mutation (emission) parameter must be > 0. Using Li & Stephens (2003) version of Watterson's estimate instead of user-supplied value...\n");
    }
    if (strcmp(argv[i], "-k") == 0) {
      if (argv[(i + 1)][0] == '-') {
        printf("Something wrong with input command line (missing arguments?). Exiting....\n");
        usage();
        exit(1);
      }
      region_size = atof(argv[(i + 1)]);
      if (region_size < 1) {
        printf("Region_size must be >= 1. Exiting...\n");
        exit(1);
      }
    }
    if (strcmp(argv[i], "-a") == 0) {
      if ((argv[(i + 1)][0] == '-') || (argv[(i + 2)][0] == '-')) {
        printf("Something wrong with input command line (missing arguments?). Exiting....\n");
        usage();
        exit(1);
      }
      start_val = atoi(argv[(i + 1)]);
      end_val = atoi(argv[(i + 2)]);
      if ((end_val < start_val) || (start_val < 0) || (end_val < 0)) {
        printf("Invalid start_ind/stop_ind vals ('-a' switch). If you want to condition each individual on every other individual, use '-a 0 0'. Exiting...\n");
        exit(1);
      }
      if (start_val > 0) start_val = start_val - 1;
    }
    if (strcmp(argv[i], "-o") == 0) {
      if (argv[(i + 1)][0] == '-') {
        printf("Something wrong with input command line (missing arguments?). Exiting....\n");
        usage();
        exit(1);
      }
      strcpy(filenameOUT, argv[(i + 1)]);
      fout = fopen(strcat(filenameOUT, ".samples.out"), "w");
      strcpy(filenameOUT, argv[(i + 1)]);
      fout2 = fopen(strcat(filenameOUT, ".priorprobs.out"), "w");
      strcpy(filenameOUT, argv[(i + 1)]);
      fout3 = fopen(strcat(filenameOUT, ".EMprobs.out"), "w");
      strcpy(filenameOUT, argv[(i + 1)]);
      fout4 = fopen(strcat(filenameOUT, ".chunkcounts.out"), "w");
      strcpy(filenameOUT, argv[(i + 1)]);
      fout5 = fopen(strcat(filenameOUT, ".chunklengths.out"), "w");
      strcpy(filenameOUT, argv[(i + 1)]);
      fout6 = fopen(strcat(filenameOUT, ".mutationprobs.out"), "w");
      strcpy(filenameOUT, argv[(i + 1)]);
      fout7 = fopen(strcat(filenameOUT, ".regionchunkcounts.out"), "w");
      strcpy(filenameOUT, argv[(i + 1)]);
      fout8 = fopen(strcat(filenameOUT, ".regionsquaredchunkcounts.out"), "w");
      strcpy(filenameOUT, argv[(i + 1)]);
      fout9 = gzopen(strcat(filenameOUT, ".copyprobsperlocus.out.gz"), "w");
      strcpy(filenameOUT, argv[(i + 1)]);
    }
  }
  if (fout == NULL) {
    printf("error opening closing file1\n");
    exit(1);
  }
  if (fout2 == NULL) {
    printf("error opening closing file2\n");
    exit(1);
  }
  if (fout3 == NULL) {
    printf("error opening closing file3\n");
    exit(1);
  }
  if (fout4 == NULL) {
    printf("error opening closing file4\n");
    exit(1);
  }
  if (fout5 == NULL) {
    printf("error opening closing file5\n");
    exit(1);
  }
  if (fout6 == NULL) {
    printf("error opening closing file6\n");
    exit(1);
  }
  if (fout7 == NULL) {
    printf("error opening closing file7\n");
    exit(1);
  }
  if (fout8 == NULL) {
    printf("error opening closing file8\n");
    exit(1);
  }
  if (fout9 == NULL) {
    printf("error opening closing file9\n");
    exit(1);
  }

  fd = fopen(filenameGEN, "r");
  if (fd == NULL) {
    printf("error opening %s\n", filenameGEN);
    exit(1);
  }

  if (fgets(line, 2047, fd) == NULL) {
    printf("Error!\n");
  }

  sscanf(line, "%d", & nind_totGENFILE);
  nind_totGENFILE = nind_totGENFILE / (2 - haploid_ind);
  fclose(fd);
  if (idfile_find == 0) {
    // get total number of inds by reading in phase file:
    // fd = fopen(filenameGEN,"r");
    //if (fd == NULL) { printf("error opening %s\n",filenameGEN); exit(1);}
    //fgets(line,2047,fd);
    //sscanf(line,"%d",&nind_tot);
    //nind_tot=nind_tot/(2-haploid_ind);
    /* ??VERY UNCLEAR WHY THIS DOESN'T WORK??
      fgets(line,2047,fd);
      fgets(line,2047,fd);
      nind_tot=0;
      while(!feof(fd))
  {
    if (fgets(bigline,2047,fd)!=NULL)
      nind_tot=nind_tot+1;
  }
      fclose(fd);
      printf("%d\n",nind_tot);
      nind_tot=nind_tot/(2-haploid_ind);
      printf("%d\n",nind_tot);
      */
  }
  if (idfile_find == 0) nind_tot = nind_totGENFILE;
  if (idfile_find == 1) {
    // open id file (to get total number of inds)
    fd4 = fopen(filenameID, "r");
    if (fd4 == NULL) {
      printf("error opening %s\n", filenameID);
      exit(1);
    }
    nind_tot = 0;
    while (!feof(fd4)) {
      if (fgets(line, 2047, fd4) != NULL)
        nind_tot = nind_tot + 1;
    }
    fclose(fd4);
    if (nind_tot != nind_totGENFILE) {
      printf("number of inds in %s (%d) does not match number of inds in %s (%d)\n. Exiting....\n", filenameID, nind_tot, filenameGEN, nind_totGENFILE);
      exit(1);
    }
  }
  // find inds to include (and get population labels for each ind):
  int * pop_vec_tot = malloc(((2 - haploid_ind) * nind_tot) * sizeof(int));
  char ** ind_label_vec = malloc(nind_tot * sizeof(char * ));
  char ** pop_label_vec = malloc(nind_tot * sizeof(char * ));
  for (i = 0; i < nind_tot; i++) {
    ind_label_vec[i] = malloc(1000 * sizeof(char));
    pop_label_vec[i] = malloc(1000 * sizeof(char));
  }
  int * include_ind_vec = malloc(nind_tot * sizeof(int));
  if (idfile_find == 0) {
    for (i = 0; i < nind_tot; i++) {
      include_ind_vec[i] = 1;
      sprintf(templab, "%d", i + 1);
      strcpy(indname, "IND");
      strcat(indname, templab);
      strcpy(ind_label_vec[i], indname);
    }
    nind = nind_tot;
  }
  if (idfile_find == 1) {
    fd4 = fopen(filenameID, "r");
    if (fd4 == NULL) {
      printf("error opening %s\n", filenameID);
      exit(1);
    }
    nind = nind_tot;
    for (i = 0; i < nind_tot; i++) {

      if (fgets(line, 2047, fd4) == NULL) {
        printf("Error!\n");
      }

      step = line;
      reading( & step, "%s", waste);
      strcpy(ind_label_vec[i], waste);
      reading( & step, "%s", waste);
      strcpy(pop_label_vec[i], waste);
      reading( & step, "%s", waste);
      include_ind_vec[i] = 1;
      if (strcmp(waste, "0") == 0) {
        include_ind_vec[i] = 0;
        nind = nind - 1;
      }
    }
    fclose(fd4);
  }
  if (all_versus_all_ind == 1) {
    if (donorlist_find == 1 && idfile_find == 1) {
      for (i = 0; i < nind_tot; i++) {
        if (include_ind_vec[i] == 1) {
          // open donor-list file (to get information on number of donor/rec pops)
          fd3 = fopen(filenameDONORLIST, "r");
          if (fd3 == NULL) {
            printf("error opening %s\n", filenameDONORLIST);
            exit(1);
          }
          ind_pop_find = 0;
          while (!feof(fd3)) {
            if (fgets(line, 2047, fd3) != NULL) {
              step = line;
              reading( & step, "%s", waste);
              reading( & step, "%s", waste2);
              if (strcmp(waste2, "D") != 0 && strcmp(waste2, "R") != 0) {
                printf("Second column of %s must contain either a 'D' (for 'donor') or 'R' (for 'recipient'). Exiting....\n", filenameDONORLIST);
                exit(1);
              }
              if (strcmp(pop_label_vec[i], waste) == 0) {
                ind_pop_find = 1;
                break;
              }
            }
          }
          fclose(fd3);
          if (ind_pop_find == 0) {
            printf("Population label for ind %s (%s) not found in %s, so excluding this individual....\n", ind_label_vec[i], pop_label_vec[i], filenameDONORLIST);
            include_ind_vec[i] = 0;
            nind = nind - 1;
          }
        }
      }
    }

    ndonors = (2 - haploid_ind) * nind - 2 + haploid_ind;
    ndonorpops = ndonors / (2 - haploid_ind);
    nrecpops = 1;
    donor_count = 0;
    for (i = 0; i < nind_tot; i++) {
      for (j = 0; j < (2 - haploid_ind); j++)
        pop_vec_tot[(2 - haploid_ind) * i + j] = -9;
      if (include_ind_vec[i] == 1) {
        for (j = 0; j < (2 - haploid_ind); j++)
          pop_vec_tot[(2 - haploid_ind) * i + j] = donor_count;
        donor_count = donor_count + 1;
      }
    }
  }
  if (all_versus_all_ind == 0) {
    // open donor-list file (to get information on number of donor/rec pops)
    fd3 = fopen(filenameDONORLIST, "r");
    if (fd3 == NULL) {
      printf("error opening %s\n", filenameDONORLIST);
      exit(1);
    }
    ndonorpops = 0;
    nrecpops = 0;
    while (!feof(fd3)) {
      if (fgets(line, 2047, fd3) != NULL) {
        step = line;
        reading( & step, "%s", waste);
        reading( & step, "%s", waste2);
        if (strcmp(waste2, "D") != 0 && strcmp(waste2, "R") != 0) {
          printf("Second column of %s must contain either a 'D' (for 'donor') or 'R' (for 'recipient'). Exiting....\n", filenameDONORLIST);
          exit(1);
        }
        if (strcmp(waste2, "D") == 0)
          ndonorpops = ndonorpops + 1;
        if (strcmp(waste2, "R") == 0)
          nrecpops = nrecpops + 1;
      }
    }
    fclose(fd3);
    if (ndonorpops == 0) {
      printf("No donor populations found in %s! Exiting....\n", filenameDONORLIST);
      exit(1);
    }
    if (nrecpops == 0) {
      printf("No recipient populations found in %s! Exiting....\n", filenameDONORLIST);
      exit(1);
    }
  }
  int * ndonorhaps_tot = malloc(ndonorpops * sizeof(int));
  int * nrechaps = malloc(nrecpops * sizeof(int));
  double * ndonorprobs = malloc(ndonorpops * sizeof(double));
  double * ndonormutrates = malloc(ndonorpops * sizeof(double));
  char ** donor_pop_vec = malloc(ndonorpops * sizeof(char * ));
  for (i = 0; i < ndonorpops; i++)
    donor_pop_vec[i] = malloc(1000 * sizeof(char));
  char ** rec_pop_vec = malloc(nrecpops * sizeof(char * ));
  for (i = 0; i < nrecpops; i++)
    rec_pop_vec[i] = malloc(1000 * sizeof(char));
  if (all_versus_all_ind == 1) {
    nrechaps[0] = 2 - haploid_ind;
    strcpy(rec_pop_vec[0], "RECPOP");
    for (i = 0; i < ndonorpops; i++) {
      sprintf(templab, "%d", i + 1);
      strcpy(donorname, "DONOR");
      strcat(donorname, templab);
      strcpy(donor_pop_vec[i], donorname);
      ndonorhaps_tot[i] = 2 - haploid_ind;
    }
    cond_nhaps = (2 - haploid_ind) * nind;
  }
  if (all_versus_all_ind == 0) {
    for (i = 0; i < ndonorpops; i++)
      ndonorhaps_tot[i] = 0;
    for (i = 0; i < nrecpops; i++)
      nrechaps[i] = 0;
    // open donor-list file (to get information on donor/rec pop labels)
    fd3 = fopen(filenameDONORLIST, "r");
    if (fd3 == NULL) {
      printf("error opening %s\n", filenameDONORLIST);
      exit(1);
    }
    donor_count = 0;
    rec_count = 0;
    totaldonorprobs = 0.0;
    while (!feof(fd3)) {
      if (fgets(line, 2047, fd3) != NULL) {
        step = line;
        reading( & step, "%s", waste);
        reading( & step, "%s", waste2);
        if (strcmp(waste2, "R") == 0) {
          strcpy(rec_pop_vec[rec_count], waste);
          rec_count = rec_count + 1;
        }
        if (strcmp(waste2, "D") == 0) {
          strcpy(donor_pop_vec[donor_count], waste);
          if (prior_donor_probs_ind == 1) {
            reading( & step, "%lf", & ndonorprobs[donor_count]);
            if (ndonorprobs[donor_count] <= 0.0000000000000001) {
              printf("Donor copying probabilities in %s must be > 0. Exiting...\n", filenameDONORLIST);
              exit(1);
            }
            totaldonorprobs = totaldonorprobs + ndonorprobs[donor_count];
          }
          if (mutation_rate_ind == 1) {
            if (prior_donor_probs_ind == 0) reading( & step, "%s", waste);
            reading( & step, "%lf", & ndonormutrates[donor_count]);
            if (ndonormutrates[donor_count] > 1) {
              printf("Donor mutation (emission) probabilities must be <= 1 in %s (use a negative number to specify default). Exiting...\n", filenameDONORLIST);
              exit(1);
            }
          }
          donor_count = donor_count + 1;
        }
      }
    }
    fclose(fd3);
    if (prior_donor_probs_ind == 1) {
      if ((totaldonorprobs > 1.00001) || (totaldonorprobs < 0.99999)) {
        printf("Probabilities across all donors in %s does not sum to 1.0 (instead sums to %lf)! Rescaling to sum to 1.0....\n", filenameDONORLIST, totaldonorprobs);
        for (k = 0; k < ndonorpops; k++)
          ndonorprobs[k] = ndonorprobs[k] / totaldonorprobs;
      }
    }

    // open id file (to get information on donor pop hap numbers and do some checks if necessary)
    fd4 = fopen(filenameID, "r");
    if (fd4 == NULL) {
      printf("error opening %s\n", filenameID);
      exit(1);
    }
    nind = nind_tot;
    for (i = 0; i < nind_tot; i++) {
      for (j = 0; j < (2 - haploid_ind); j++)
        pop_vec_tot[(2 - haploid_ind) * i + j] = -9;

      if (fgets(line, 2047, fd4) == NULL) {
        printf("Error!\n");
      }

      step = line;
      reading( & step, "%s", waste);
      reading( & step, "%s", waste);
      reading( & step, "%s", waste);
      if (include_ind_vec[i] == 0) nind = nind - 1;
      if (include_ind_vec[i] == 1) {
        ind_label_find = 0;
        for (k = 0; k < ndonorpops; k++) {
          if (strcmp(pop_label_vec[i], donor_pop_vec[k]) == 0) {
            ind_label_find = 1;
            for (j = 0; j < (2 - haploid_ind); j++)
              pop_vec_tot[(2 - haploid_ind) * i + j] = k;
            ndonorhaps_tot[k] = ndonorhaps_tot[k] + 2 - haploid_ind;
            break;
          }
        }
        for (k = 0; k < nrecpops; k++) {
          if (strcmp(pop_label_vec[i], rec_pop_vec[k]) == 0) {
            ind_label_find = 1;
            nrechaps[k] = nrechaps[k] + 2 - haploid_ind;
            break;
          }
        }
        if (ind_label_find == 0) {
          include_ind_vec[i] = 0;
          nind = nind - 1;
          printf("Individual %d (%s) in %s has a population label (%s) not found in %s, so is being excluded.....\n", i + 1, ind_label_vec[i], filenameID, pop_label_vec[i], filenameDONORLIST);
          //printf("Individual %d in %s has a population label (%s) not found in %s, which is not allowed (unless '-a' switch used)! Exiting....\n",ind_count+1;filenameID,pop_label_vec[indcount],filenameDONORLIST);
          //exit(1);
        }
      }
    }
    fclose(fd4);
    ndonors = 0;
    for (k = 0; k < ndonorpops; k++) {
      if (ndonorhaps_tot[k] == 0) {
        printf("No individuals found in %s with population label %s found in %s, which is not allowed (unless '-a' switch used)! Exiting....\n", filenameID, donor_pop_vec[k], filenameDONORLIST);
        exit(1);
      }
      for (j = 0; j < nrecpops; j++) {
        if (strcmp(donor_pop_vec[k], rec_pop_vec[j]) == 0) {
          printf("WARNING: Population %s is specified as a donor and as a recipient. Will allow 'self-copying' of own population label in recipient individuals....\n", rec_pop_vec[j]);
          //ndonorhaps_tot[k]=ndonorhaps_tot[k]-2+haploid_ind;
          break;
        }
      }
      ndonors = ndonors + ndonorhaps_tot[k];
    }
    cond_nhaps = 0;
    for (k = 0; k < nrecpops; k++) {
      if (nrechaps[k] == 0) {
        printf("No individuals found in %s with population label %s found in %s, which is not allowed (unless '-a' switch used)! Exiting....\n", filenameID, rec_pop_vec[k], filenameDONORLIST);
        exit(1);
      }
      cond_nhaps = cond_nhaps + nrechaps[k];
    }
  }
  cond_nind = cond_nhaps / (2 - haploid_ind);
  // (0) INITIALIZE copy_prob_tot, copy_probSTART_tot, MutProb_vec_tot:
  double * MutProb_vec_tot = malloc(((2 - haploid_ind) * nind_tot) * sizeof(double));
  double * copy_prob_tot = malloc(((2 - haploid_ind) * nind_tot) * sizeof(double));
  double * copy_probSTART_tot = malloc(((2 - haploid_ind) * nind_tot) * sizeof(double));
  for (i = 0; i < nind_tot; i++) {
    for (j = 0; j < (2 - haploid_ind); j++) {
      MutProb_vec_tot[(2 - haploid_ind) * i + j] = 0;
      copy_prob_tot[(2 - haploid_ind) * i + j] = -9;
      copy_probSTART_tot[(2 - haploid_ind) * i + j] = -9;
    }
    if (include_ind_vec[i] == 1 && prior_donor_probs_ind == 0) {
      for (j = 0; j < (2 - haploid_ind); j++) {
        //copy_prob_tot[(2-haploid_ind)*i+j] = 1.0/ndonors;
        //copy_probSTART_tot[(2-haploid_ind)*i+j] = 1.0/ndonors;
        copy_prob_tot[(2 - haploid_ind) * i + j] = 1.0;
        copy_probSTART_tot[(2 - haploid_ind) * i + j] = 1.0;
      }
    }
    if (include_ind_vec[i] == 1 && prior_donor_probs_ind == 1) {
      for (k = 0; k < ndonorpops; k++) {
        if (pop_vec_tot[(2 - haploid_ind) * i] == k) {
          for (j = 0; j < (2 - haploid_ind); j++) {
            copy_prob_tot[(2 - haploid_ind) * i + j] = ndonorprobs[k] / ndonorhaps_tot[k];
            copy_probSTART_tot[(2 - haploid_ind) * i + j] = ndonorprobs[k] / ndonorhaps_tot[k];
          }
          break;
        }
      }
    }
    if (include_ind_vec[i] == 1 && (mutation_rate_ind == 0 || mutationALL_em_find == 1)) {
      for (j = 0; j < (2 - haploid_ind); j++)
        MutProb_vec_tot[(2 - haploid_ind) * i + j] = GlobalMutRate;
    }
    if (include_ind_vec[i] == 1 && mutation_rate_ind == 1 && mutationALL_em_find == 0) {
      for (k = 0; k < ndonorpops; k++) {
        if (pop_vec_tot[(2 - haploid_ind) * i] == k) {
          for (j = 0; j < (2 - haploid_ind); j++)
            MutProb_vec_tot[(2 - haploid_ind) * i + j] = ndonormutrates[k];
          break;
        }
      }
    }
  }

  // open haplotype file (to get information on haplotype numbers)
  fd = fopen(filenameGEN, "r");
  if (fd == NULL) {
    printf("error opening %s\n", filenameGEN);
    exit(1);
  }
  Data = ReadData(fd, 0, 2 - haploid_ind, include_ind_vec, pop_label_vec, ndonorpops, donor_pop_vec, pop_vec_tot, copy_prob_tot, copy_probSTART_tot, MutProb_vec_tot, ndonorhaps_tot, all_versus_all_ind);

  nsites = Data -> nsnps;
  if (idfile_find == 1 && (nind_tot != (((int) Data -> nhaps) / (2 - haploid_ind)))) {
    printf("Total number of individuals in %s does not match that in %s. Exiting....\n", filenameGEN, filenameID);
    exit(1);
  }

  if (all_versus_all_ind == 1) printf("Will condition each individual on every other individual...\n");
  if (all_versus_all_ind == 1 && idfile_find == 1) printf("Will use %s for population labels in output (even though performing all-versus-all; i.e. '-a' switch)\n", filenameID);
  //if (all_versus_all_ind==1 && donorlist_find==1) printf("As you have specified '-a' switch for all-versus-all, will ignore '-f' input, including file %s listing donor/recipient populations.\n",filenameDONORLIST);
  if (all_versus_all_ind == 1 && donorlist_find == 1 && idfile_find == 0) printf("As you have specified '-a' switch for all-versus-all but have not specified label file with '-t', will ignore '-f' input, i.e. file %s listing populations to include in all-versus-all analysis.\n", filenameDONORLIST);
  if (all_versus_all_ind == 1 && donorlist_find == 1 && idfile_find == 1) printf("As you have specified '-a' switch for all-versus-all and supplied '-t' and '-f' files, will only consider (non-excluded) individuals in file %s that are from all populations listed in file %s.\n", filenameID, filenameDONORLIST);
  if (prior_donor_probs_ind == 1)
    printf("Using specified prior donor probs from input file....\n");
  if ((mutation_rate_ind == 1) && (mutationALL_em_find == 0))
    printf("Using specified donor population mutation rates from input file....\n");
  if (copy_prop_em_find == 1)
    printf("Running E-M to estimate copying proportions....\n");
  if (recom_em_find == 1)
    printf("Running E-M to estimate N_e....\n");
  if (mutation_em_find == 1)
    printf("Running E-M to estimate mutation (emission) probabilities....\n");
  if (mutationALL_em_find == 1)
    printf("Running E-M to estimate global mutation (emission) probability....\n");
  printf(" Total number of individuals found: %d\n", nind_tot);
  printf(" Total number of individuals included for analysis: %d\n", nind);
  if (all_versus_all_ind == 0) {
    printf(" Number of donor pops: %d\n", ndonorpops);
    printf(" Number of recipient pops: %d\n", nrecpops);
  }
  if (all_versus_all_ind == 1) {
    printf(" Number of donor pops: %d (all-versus-all)\n", ndonorpops);
    printf(" Number of recipient pops: %d (all-versus-all)\n", nrecpops);
  }
  if (print_file9_ind == 1) {
    gzprintf(fout9, "pos");
    if (all_versus_all_ind == 0) {
      for (k = 0; k < ndonorpops; k++)
        gzprintf(fout9, " %s", donor_pop_vec[k]);
    }
    if (all_versus_all_ind == 1) {
      for (j = 0; j < nind_tot; j++) {
        if (include_ind_vec[j] == 1)
          gzprintf(fout9, " %s", ind_label_vec[j]);
      }
    }
    gzprintf(fout9, "\n");
  }
  //if (ne_find==0) N_e=N_e/ndonors;

  if (start_val >= cond_nind) {
    printf("Fewer recipient individuals than expected -- either your '-a' switch is mis-specified or there is something wrong with your input file %s? Exiting....\n", filenameGEN);
    exit(1);
  }
  if (drift_calc_ind == 0) num_rec_drift = 0;
  if (drift_calc_ind == 1 && num_rec_drift > cond_nind) {
    printf("You have specified to store %d individuals in memory using the '-d' switch, but you have only %d recipient individuals. Will set this value to %d.\n", num_rec_drift, cond_nind, cond_nind);
    num_rec_drift = cond_nind;
  }
  if (drift_calc_ind == 1 && num_rec_drift == 0) num_rec_drift = cond_nind;
  if (end_val == 0 || end_val > cond_nind) end_val = cond_nind;
  int * recipient_ind_vec = malloc(nind_tot * sizeof(int));
  rec_ind_count = 0;
  rec_ind_topaint_count = 0;
  for (i = 0; i < nind_tot; i++) {
    recipient_ind_vec[i] = 0;
    if (idfile_find == 0) {
      if (i >= start_val && i < end_val) {
        recipient_ind_vec[i] = 1;
        rec_ind_topaint_count = rec_ind_topaint_count + 1;
      }
    }
    if (idfile_find == 1) {
      if (all_versus_all_ind == 0 && include_ind_vec[i] == 1) {
        ind_label_find = 0;
        for (k = 0; k < nrecpops; k++) {
          if (strcmp(pop_label_vec[i], rec_pop_vec[k]) == 0) {
            ind_label_find = 1;
            break;
          }
        }
        if (ind_label_find == 1) {
          if (rec_ind_count >= start_val && rec_ind_count < end_val) {
            recipient_ind_vec[i] = 1;
            rec_ind_topaint_count = rec_ind_topaint_count + 1;
          }
          rec_ind_count = rec_ind_count + 1;
        }
      }
      if (all_versus_all_ind == 1 && include_ind_vec[i] == 1) {
        if (rec_ind_count >= start_val && rec_ind_count < end_val) {
          recipient_ind_vec[i] = 1;
          rec_ind_topaint_count = rec_ind_topaint_count + 1;
        }
        rec_ind_count = rec_ind_count + 1;
      }
    }
  }

  if (haploid_ind == 1) printf("Assuming all inds are haploid....\n");
  if (unlinked_ind == 1) {
    printf("Assuming sites are unlinked....\n");
    //EMruns=0;
  }
  //printf(" Number of donor haplotypes = %d\n Number of recipient haplotypes = %d\n",ndonors,(2-haploid_ind)*nind-ndonors);
  printf(" Number of donor haplotypes = %d\n Number of recipient haplotypes = %d\n", ndonors, cond_nhaps);
  if (ne_find == 1) printf(" Number of EM-runs = %d\n Number of samples = %d\n N_e value = %lf\n Region size = %lf\n", EMruns, samplesTOT, N_e, region_size);
  if (ne_find == 0) printf(" Number of EM-runs = %d\n Number of samples = %d\n N_e value = %lf (divided by number of donor haplotypes)\n Region size = %lf\n", EMruns, samplesTOT, N_e, region_size);
  printf(" Global mutation value = %lf\n", GlobalMutRate);
  printf(" Painting %d of %d recipient individuals....\n", rec_ind_topaint_count, cond_nind);

  if (ne_find == 1) fprintf(fout, "EM_iter = %d (N_e = %d / copy_prop = %d / mutation = %d / mutationGLOBAL = %d), nsamples = %d, N_e_start = %lf, region_size = %lf, haplotype dataset = %s, genmap dataset = %s, pop-list dataset = %s, id-label dataset = %s\n", EMruns, recom_em_find, copy_prop_em_find, mutation_em_find, mutationALL_em_find, samplesTOT, N_e, region_size, filenameGEN, filename, filenameDONORLIST, filenameID);
  if (ne_find == 0) fprintf(fout, "EM_iter = %d (N_e = %d / copy_prop = %d / mutation = %d / mutationGLOBAL = %d), nsamples = %d, N_e_start = %lf (divided by number of donor haplotypes), region_size = %lf, haplotype dataset = %s, genmap dataset = %s, pop-list dataset = %s, id-label dataset = %s\n", EMruns, recom_em_find, copy_prop_em_find, mutation_em_find, mutationALL_em_find, samplesTOT, N_e, region_size, filenameGEN, filename, filenameDONORLIST, filenameID);

  // (i) GET SAMPLES, RUN E-M:
  // open recomb map file:
  double * recom_map = malloc((Data -> nsnps - 1) * sizeof(double));
  if ((unlinked_ind == 1) && (recom_find == 0)) {
    for (j = 0; j < (Data -> nsnps - 1); j++) recom_map[j] = -9.0;
  }
  if (recom_find == 1) {
    fd2 = fopen(filename, "r");
    if (fd2 == NULL) {
      printf("error opening recom map input file: %s\n", filename);
      exit(1);
    }

    if (fgets(line, 2047, fd2) == NULL) { // header
      printf("Error!\n");
    }

    for (j = 0; j < (Data -> nsnps - 1); j++) {

      if (fgets(line, 2047, fd2) == NULL) { // header
        printf("Error!\n");
      }

      step = line;

      reading( & step, "%lf", & bpval); // basepair position
      if (bpval != Data -> positions[j]) {
        printf("basepair positions do not match between %s and %s at basepair %d. Exiting....\n", filename, filenameGEN, j + 1);
        exit(1);
      }
      reading( & step, "%lf", & recom_map[j]);
      if (recom_map[j] >= 0 && recom_map[j] <= small_recom_val) {
        printf("recom rate very low at basepair %lf (%lf). Assuming recomb rate between this snp and next one is %lf....\n", Data -> positions[j], recom_map[j], small_recom_val);
        recom_map[j] = small_recom_val;
      }
      if (recom_map[j] < 0) {
        printf("recom rate < 0 at basepair %lf. Assuming transition probability of 1 between this snp and next one....\n", Data -> positions[j]);
        //printf("recom rate must be > 0 (basepair %lf)!! Exiting....\n",Data->positions[j]);
        //exit(1);
      }
    }

    if (fgets(line, 2047, fd2) == NULL) { // header
      printf("Error!\n");
    }

    step = line;
    reading( & step, "%lf", & bpval); // basepair position
    if (bpval != Data -> positions[(Data -> nsnps - 1)]) {
      printf("basepair positions do not match between %s and %s at basepair %d. Exiting....\n", filename, filenameGEN, Data -> nsnps);
      exit(1);
    }
    fclose(fd2);
  }
  // check ordering of SNPs (only allowed to be less than previous position if recom_map<0 at position -- i.e. suggesting new chromosome):
  for (i = 0; i < Data -> nsnps; i++) {
    if (i > 0) {
      if (Data -> positions[i] <= Data -> positions[(i - 1)] && (recom_map[(i - 1)] >= 0)) {
        printf("positions in %s are not listed in increasing order at basepairs %lf and %lf. Exiting....\n", filenameGEN, Data -> positions[(i - 1)], Data -> positions[i]);
        exit(1);
      }
    }
  }
  DestroyData(Data, 2 - haploid_ind);
  fclose(fd);

  EMruns = EMruns + 1;

  /* print-out headers for copy-props, chunk counts, lengths, and differences: */
  fprintf(fout2, "Recipient");
  fprintf(fout4, "Recipient");
  fprintf(fout5, "Recipient");
  fprintf(fout6, "Recipient");
  fprintf(fout7, "Recipient num.regions");
  fprintf(fout8, "Recipient num.regions");
  if (all_versus_all_ind == 0) {
    for (k = 0; k < ndonorpops; k++) {
      fprintf(fout2, " %s", donor_pop_vec[k]);
      fprintf(fout4, " %s", donor_pop_vec[k]);
      fprintf(fout5, " %s", donor_pop_vec[k]);
      fprintf(fout6, " %s", donor_pop_vec[k]);
      fprintf(fout7, " %s", donor_pop_vec[k]);
      fprintf(fout8, " %s", donor_pop_vec[k]);
    }
  }
  if (all_versus_all_ind == 1) {
    for (j = 0; j < nind_tot; j++) {
      if (include_ind_vec[j] == 1) {
        fprintf(fout2, " %s", ind_label_vec[j]);
        fprintf(fout4, " %s", ind_label_vec[j]);
        fprintf(fout5, " %s", ind_label_vec[j]);
        fprintf(fout6, " %s", ind_label_vec[j]);
        fprintf(fout7, " %s", ind_label_vec[j]);
        fprintf(fout8, " %s", ind_label_vec[j]);
      }
    }
  }
  fprintf(fout2, "\n");
  fprintf(fout4, "\n");
  fprintf(fout5, "\n");
  fprintf(fout6, "\n");
  fprintf(fout7, "\n");
  fprintf(fout8, "\n");

  log_lik_check = loglik(nind_tot, ndonors, & nsites, (2 - haploid_ind) * nind, N_e, recom_map, MutProb_vec_tot, samplesTOT, ndonorpops, ndonorhaps_tot, include_ind_vec, ind_label_vec, pop_label_vec, donor_pop_vec, nrecpops, rec_pop_vec, copy_prob_tot, copy_probSTART_tot, pop_vec_tot, region_size, EMruns, copy_prop_em_find, recom_em_find, ne_find, mutation_em_find, mutationALL_em_find, all_versus_all_ind, prior_donor_probs_ind, rec_ind_count, recipient_ind_vec, filenameGEN, donorlist_find, haploid_ind, unlinked_ind, print_file9_ind, drift_calc_ind, num_rec_drift, rec_ind_topaint_count, fout, fout2, fout3, fout4, fout5, fout6, fout7, fout8, fout9);
  if (log_lik_check != 1) {
    printf("Algorithm failed. Check input files and parameters. Exiting....\n");
    exit(1);
  }

  for (i = 0; i < nind_tot; i++) {
    free(ind_label_vec[i]);
    free(pop_label_vec[i]);
  }
  for (i = 0; i < ndonorpops; i++)
    free(donor_pop_vec[i]);
  for (i = 0; i < nrecpops; i++)
    free(rec_pop_vec[i]);
  free(ind_label_vec);
  free(pop_label_vec);
  free(donor_pop_vec);
  free(rec_pop_vec);
  free(include_ind_vec);
  free(recipient_ind_vec);
  free(recom_map);
  free(filename);
  free(filenameGEN);
  free(filenameDONORLIST);
  free(filenameID);
  free(copy_prob_tot);
  free(copy_probSTART_tot);
  free(pop_vec_tot);
  free(ndonorhaps_tot);
  free(nrechaps);
  free(ndonorprobs);
  free(ndonormutrates);
  free(MutProb_vec_tot);
  //free(bigline);

  fclose(fout);
  fclose(fout2);
  fclose(fout3);
  fclose(fout4);
  fclose(fout5);
  fclose(fout6);
  fclose(fout7);
  fclose(fout8);
  gzclose(fout9);

  printf("Finished!\n");

  return 0;
}