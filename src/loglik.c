#include <stdio.h>
#include <stdlib.h>
#include "data_t.h"
#include "destroy.h"
#include "reading.h"
#include "sampler.h"
#include "zlib.h"

int loglik(int nind_tot, int nhaps_startpop, int * p_nloci, int p_nhaps, double N_e_start, double * recom_map, double * MutProb_vec_tot, int nsampTOT, int ndonorpops, int * ndonorhaps_tot, int * include_ind_vec, char ** ind_label_vec, char ** pop_label_vec, char ** donor_pop_vec, int nrecpops, char ** rec_pop_vec, double * copy_prob_tot, double * copy_probSTART_tot, int * pop_vec_tot, double region_size, int EMruns, int estimate_copyprob_ind, int estimate_recom_ind, int ne_find, int estimate_mutation_ind, int estimate_mutationALL_ind, int all_versus_all_ind, int prior_donor_probs_ind, int num_rec_ind, int * recipient_ind_vec, char * filename, int donorlist_ind, int haploid_ind, int unlinked_ind, int print_file9_ind, int drift_calc_ind, int num_rec_drift, int rec_ind_topaint_count, FILE * fout, FILE * fout2, FILE * fout3, FILE * fout4, FILE * fout5, FILE * fout6, FILE * fout7, FILE * fout8, gzFile fout9) {
  double small_copy_val = 0.000000000000001; // (!!!) copy props per hap not allowed to go below this value, even if E-M wants to make them lower (!!!)
  int nhaps_condpop, nind_condpop;

  FILE * fd;
  struct data_t * Data;

  int i, j, m, n, r, h;
  int nhaps, num_regions_tot;
  int included_count, included_count_bigloop;

  double sum_total_diff;
  double * total_back_prob = malloc(ndonorpops * sizeof(double));
  double * total_back_probSTART = malloc(ndonorpops * sizeof(double));
  double * total_counts = malloc(ndonorpops * sizeof(double));
  double * total_lengths = malloc(ndonorpops * sizeof(double));
  double * total_differences = malloc(ndonorpops * sizeof(double));
  double * total_region_counts = malloc(ndonorpops * sizeof(double));
  double * total_squared_region_counts = malloc(ndonorpops * sizeof(double));
  double * snp_info_measure_final = malloc(ndonorpops * sizeof(double));
  double N_e_new;
  double N_e = 0;

  double total_prob, total_probSTART;
  double ** copy_prob_pop = malloc(2 * sizeof(double * ));

  int * newhap = malloc(( * p_nloci) * sizeof(int));
  int drift_malloc_size, drift_malloc_size2, drift_malloc_size3;

  drift_malloc_size = num_rec_drift;
  drift_malloc_size2 = num_rec_ind;
  drift_malloc_size3 = * p_nloci;

  if (drift_malloc_size == 0) {
    drift_malloc_size = 1;
    drift_malloc_size2 = 1;
    drift_malloc_size3 = 1;
  }

  double * correlated_drift_vec = malloc(drift_malloc_size * sizeof(double));

  nhaps_condpop = p_nhaps - nhaps_startpop;
  nind_condpop = nhaps_condpop / (2 - haploid_ind);
  double ** correlated_drift_calc = malloc(drift_malloc_size2 * sizeof(double * ));
  double ** * RecChromProbArray = malloc(drift_malloc_size * sizeof(double ** ));
  for (i = 0; i < drift_malloc_size2; i++)
    correlated_drift_calc[i] = malloc(nind_condpop * sizeof(double));
  for (i = 0; i < drift_malloc_size; i++) {
    RecChromProbArray[i] = malloc(ndonorpops * sizeof(double * ));
    for (j = 0; j < ndonorpops; j++)
      RecChromProbArray[i][j] = malloc(drift_malloc_size3 * sizeof(double));
  }
  for (i = 0; i < 2; i++)
    copy_prob_pop[i] = malloc(ndonorpops * sizeof(double));
  for (i = 0; i < ndonorpops; i++) {
    total_back_prob[i] = 0.0;
    total_back_probSTART[i] = 0.0;
  }

  total_prob = 0.0;
  total_probSTART = 0.0;
  int * allelic_type_count_vec = malloc( * p_nloci * sizeof(int));
  int * found_vec = malloc(6 * sizeof(int));
  included_count_bigloop = 0;
  for (m = 0; m < nind_tot; m++) {
    if (recipient_ind_vec[m] == 1) {
      printf(" .....Painting recipient individual %d of %d......\n", included_count_bigloop + 1, rec_ind_topaint_count);
      fprintf(fout3, "%s\n", ind_label_vec[m]);

      fd = fopen(filename, "r");
      if (fd == NULL) {
        printf("error opening %s\n", filename);
        exit(1);
      }
      Data = ReadData(fd, m, 2 - haploid_ind, include_ind_vec, pop_label_vec, ndonorpops, donor_pop_vec, pop_vec_tot, copy_prob_tot, copy_probSTART_tot, MutProb_vec_tot, ndonorhaps_tot, all_versus_all_ind);

      for (i = 0; i < Data -> ndonorhaps; i++) {
        Data -> copy_prob_new[i] = Data -> copy_prob[i];
        Data -> copy_prob_newSTART[i] = Data -> copy_probSTART[i];
        if (prior_donor_probs_ind == 0) {
          Data -> copy_prob_new[i] = Data -> copy_prob[i] / Data -> ndonorhaps;
          Data -> copy_prob_newSTART[i] = Data -> copy_probSTART[i] / Data -> ndonorhaps;
        }
        Data -> MutProb_vec_new[i] = Data -> MutProb_vec[i];
      }

      // find number of alleles per snp (this is NOT every used, but perhaps should be to get default mutation rate correct):
      for (j = 0; j < Data -> nsnps; j++) {
        allelic_type_count_vec[j] = 0;
        for (i = 0; i < 6; i++)
          found_vec[i] = 0;
        for (i = 0; i < Data -> ndonorhaps; i++) {
          if ((Data -> cond_chromosomes[i][j] == 0) && (found_vec[0] == 0)) {
            allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
            found_vec[0] = 1;
          }
          if ((Data -> cond_chromosomes[i][j] == 1) && (found_vec[1] == 0)) {
            allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
            found_vec[1] = 1;
          }
          if ((Data -> cond_chromosomes[i][j] == 2) && (found_vec[2] == 0)) {
            allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
            found_vec[2] = 1;
          }
          if ((Data -> cond_chromosomes[i][j] == 3) && (found_vec[3] == 0)) {
            allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
            found_vec[3] = 1;
          }
          if ((Data -> cond_chromosomes[i][j] == 4) && (found_vec[4] == 0)) {
            allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
            found_vec[4] = 1;
          }
          if ((Data -> cond_chromosomes[i][j] == 5) && (found_vec[5] == 0)) {
            allelic_type_count_vec[j] = allelic_type_count_vec[j] + 1;
            found_vec[5] = 1;
          }
        }
      }

      for (r = 0; r < EMruns; r++) {
        fprintf(fout3, "%d", r);

        total_prob = 0.0;
        total_probSTART = 0.0;
        for (i = 0; i < ndonorpops; i++) {
          total_back_prob[i] = 0.0;
          total_back_probSTART[i] = 0.0;
          total_counts[i] = 0.0;
          total_lengths[i] = 0.0;
          total_differences[i] = 0.0;
          total_region_counts[i] = 0.0;
          total_squared_region_counts[i] = 0.0;
          snp_info_measure_final[i] = 0.0;
        }

        N_e_new = 0.0;
        num_regions_tot = 0;
        for (h = 0; h < (2 - haploid_ind); h++) {
          if (r == 0) {
            if (ne_find == 0) N_e = N_e_start / Data -> ndonorhaps;
            if (ne_find == 1) N_e = N_e_start;
          }

          for (n = 0; n < Data -> nsnps; n++)
            newhap[n] = Data -> ind_chromosomes[h][n];

          nhaps = Data -> ndonorhaps;

          if (r == (EMruns - 1)) {
            fprintf(fout, "HAP %d %s\n", h + 1, ind_label_vec[m]);
            if (print_file9_ind == 1) gzprintf(fout9, "HAP %d %s\n", h + 1, ind_label_vec[m]);
          }

          /* SAMPLE FROM PAC CONDITIONAL ON COPY-PROBS: */
          Data -> back_prob = sampler(newhap, Data -> cond_chromosomes, p_nloci, & nhaps, & p_nhaps, N_e, Data -> MutProb_vec_new, allelic_type_count_vec, recom_map, Data -> positions, nsampTOT, Data -> copy_prob_new, Data -> copy_prob_newSTART, Data -> pop_vec, ndonorpops, region_size, r, EMruns, all_versus_all_ind, haploid_ind, unlinked_ind, included_count_bigloop, Data -> hap_label_vec, print_file9_ind, fout, fout3, fout9);

          N_e_new = N_e_new + Data -> back_prob[0][Data -> ndonorhaps];
          num_regions_tot = num_regions_tot + Data -> back_prob[1][Data -> ndonorhaps];

          /* GET NEW COPY-PROBS BASED ON PAC SAMPLES: */
          for (i = 0; i < Data -> ndonorhaps; i++) {
            total_back_prob[Data -> pop_vec[i]] = total_back_prob[Data -> pop_vec[i]] + Data -> back_prob[0][i];
            total_prob = total_prob + Data -> back_prob[0][i];
            total_back_probSTART[Data -> pop_vec[i]] = total_back_probSTART[Data -> pop_vec[i]] + Data -> back_prob[1][i];
            total_probSTART = total_probSTART + Data -> back_prob[1][i];
            total_counts[Data -> pop_vec[i]] = total_counts[Data -> pop_vec[i]] + Data -> back_prob[2][i];
            total_lengths[Data -> pop_vec[i]] = total_lengths[Data -> pop_vec[i]] + Data -> back_prob[3][i];
            total_differences[Data -> pop_vec[i]] = total_differences[Data -> pop_vec[i]] + Data -> back_prob[4][i];
          }
          for (i = 0; i < ndonorpops; i++) {
            total_region_counts[i] = total_region_counts[i] + Data -> back_prob[5][i];
            total_squared_region_counts[i] = total_squared_region_counts[i] + Data -> back_prob[6][i];
            snp_info_measure_final[i] = snp_info_measure_final[i] + Data -> back_prob[7][i];
          }
        }
        fprintf(fout3, " %.10lf %.10lf\n", N_e, Data -> MutProb_vec_new[0]);
        if (estimate_recom_ind == 1) N_e = N_e_new / (2.0 - haploid_ind);
        for (i = 0; i < ndonorpops; i++) {
          copy_prob_pop[0][i] = total_back_prob[i] / total_prob;
          copy_prob_pop[1][i] = total_back_probSTART[i] / total_probSTART;
        }
        /* RESET COPY-PROBS and MUTATION-PROBS: */
        // (first check for probabilities of 0:)
        for (i = 0; i < ndonorpops; i++) {
          if (copy_prob_pop[0][i] <= 0)
            copy_prob_pop[0][i] = small_copy_val * Data -> ndonorhaps_vec[i];
          if (copy_prob_pop[1][i] <= 0)
            copy_prob_pop[1][i] = small_copy_val * Data -> ndonorhaps_vec[i];
        }
        total_prob = 0.0;
        total_probSTART = 0.0;
        for (j = 0; j < ndonorpops; j++) {
          total_prob = total_prob + copy_prob_pop[0][j];
          total_probSTART = total_probSTART + copy_prob_pop[1][j];
        }
        for (j = 0; j < ndonorpops; j++) {
          copy_prob_pop[0][j] = copy_prob_pop[0][j] / total_prob;
          copy_prob_pop[1][j] = copy_prob_pop[1][j] / total_probSTART;
        }
        if (estimate_copyprob_ind == 1) {
          for (i = 0; i < Data -> ndonorhaps; i++) {
            for (j = 0; j < ndonorpops; j++) {
              if (Data -> pop_vec[i] == j) {
                if (Data -> ndonorhaps_vec[j] > 0) {
                  Data -> copy_prob_new[i] = copy_prob_pop[0][j] / Data -> ndonorhaps_vec[j];
                  Data -> copy_prob_newSTART[i] = copy_prob_pop[1][j] / Data -> ndonorhaps_vec[j];
                }
                if (Data -> ndonorhaps_vec[j] == 0) {
                  Data -> copy_prob_new[i] = 0.0;
                  Data -> copy_prob_newSTART[i] = 0.0;
                }
                break;
              }
            }
          }
        }
        if (estimate_mutation_ind == 1) {
          for (i = 0; i < Data -> ndonorhaps; i++) {
            for (j = 0; j < ndonorpops; j++) {
              if (Data -> pop_vec[i] == j) {
                Data -> MutProb_vec_new[i] = total_differences[j] / ( * p_nloci * (2 - haploid_ind));
                break;
              }
            }
          }
        }

        if (estimate_mutationALL_ind == 1) {
          sum_total_diff = 0.0;
          for (i = 0; i < ndonorpops; i++)
            sum_total_diff = sum_total_diff + total_differences[i] / ( * p_nloci * (2 - haploid_ind));
          for (i = 0; i < Data -> ndonorhaps; i++)
            Data -> MutProb_vec_new[i] = sum_total_diff;
        }

        if (r == (EMruns - 1)) {
          /* print props, lengths, counts, and differences: */
          fprintf(fout2, "%s", ind_label_vec[m]);
          fprintf(fout4, "%s", ind_label_vec[m]);
          fprintf(fout5, "%s", ind_label_vec[m]);
          fprintf(fout6, "%s", ind_label_vec[m]);
          fprintf(fout7, "%s", ind_label_vec[m]);
          fprintf(fout8, "%s", ind_label_vec[m]);
          fprintf(fout7, " %d", num_regions_tot);
          fprintf(fout8, " %d", num_regions_tot);
          if (all_versus_all_ind == 0) {
            for (j = 0; j < ndonorpops; j++) {
              fprintf(fout2, " %lf", copy_prob_pop[0][j]);
              fprintf(fout4, " %lf", total_counts[j]);
              fprintf(fout5, " %lf", total_lengths[j]);
              fprintf(fout6, " %lf", total_differences[j]);
              fprintf(fout7, " %lf", total_region_counts[j]);
              fprintf(fout8, " %lf", total_squared_region_counts[j]);
            }
          }
          if (all_versus_all_ind == 1) {
            included_count = 0;
            for (j = 0; j < nind_tot; j++) {
              if (j == m) {
                fprintf(fout2, " 0.00");
                fprintf(fout4, " 0.00");
                fprintf(fout5, " 0.00");
                fprintf(fout6, " 0.00");
                fprintf(fout7, " 0.00");
                fprintf(fout8, " 0.00");
              }
              if (include_ind_vec[j] == 1 && j != m) {
                fprintf(fout2, " %lf", copy_prob_pop[0][included_count]);
                fprintf(fout4, " %lf", total_counts[included_count]);
                fprintf(fout5, " %lf", total_lengths[included_count]);
                fprintf(fout6, " %lf", total_differences[included_count]);
                fprintf(fout7, " %lf", total_region_counts[included_count]);
                fprintf(fout8, " %lf", total_squared_region_counts[included_count]);
                included_count = included_count + 1;
              }
            }
          }
          fprintf(fout2, "\n");
          fprintf(fout4, "\n");
          fprintf(fout5, "\n");
          fprintf(fout6, "\n");
          fprintf(fout7, "\n");
          fprintf(fout8, "\n");
        }
      }
      DestroyData(Data, 2 - haploid_ind);
      fclose(fd);
    }
    if (include_ind_vec[m] == 1 && recipient_ind_vec[m] == 1)
      included_count_bigloop = included_count_bigloop + 1;
  }

  free(newhap);
  for (i = 0; i < drift_malloc_size2; i++)
    free(correlated_drift_calc[i]);
  for (i = 0; i < drift_malloc_size; i++) {
    for (j = 0; j < ndonorpops; j++)
      free(RecChromProbArray[i][j]);
    free(RecChromProbArray[i]);
  }
  free(RecChromProbArray);
  free(total_back_prob);
  free(total_back_probSTART);
  free(total_counts);
  free(total_lengths);
  free(total_differences);
  free(total_region_counts);
  free(total_squared_region_counts);
  free(snp_info_measure_final);
  free(allelic_type_count_vec);
  free(found_vec);
  free(correlated_drift_vec);
  free(correlated_drift_calc);

  return (1);
}