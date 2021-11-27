#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "zlib.h"
#include "data_t.h"

double ** sampler(int * newh, int ** existing_h, int * p_Nloci, int * p_Nhaps, int * p_nchr, double p_rhobar, double * MutProb_vec, int * allelic_type_count_vec, double * lambda, double * pos, int nsampTOT, double * copy_prob, double * copy_probSTART, int * pop_vec, int ndonorpops, double region_size, int run_num, int run_thres, int all_versus_all_ind, int haploid_ind, int unlinked_ind, int ind_val, int * hap_label_vec, int print_file9_ind, FILE * fout, FILE * fout3, gzFile fout9) {
  double rounding_val = 1.0 / 10000000.0; // for regional_counts; c is a bit lame
  //double small_missing_val=0.000000000000001;    // prob used when donor data is missing at a SNP (!!! NOT CURRENTLY IMPLEMENTED !!!)

  int i, j, locus;
  double sum, Theta;
  double prob, total_prob, total_gen_dist;
  double total_prob_from_i_to_i, total_prob_to_i_exclude_i, total_prob_from_i_exclude_i, constant_exclude_i, constant_from_i_to_i, constant_exclude_i_both_sides, total_prob_from_any_to_any_exclude_i;
  double total_regional_chunk_count, total_ind_sum;
  int num_regions;
  double * TransProb = malloc((( * p_Nloci) - 1) * sizeof(double));
  double N_e_new = 0;
  int * sample_state = malloc( * p_Nloci * sizeof(int));
  double ObsStateProb, ObsStateProbPREV;
  double delta;
  //correction to PAC-A rho_est
  double random_unif, random_unifSWITCH;
  double no_switch_prob;
  double ** Alphamat = malloc( * p_Nhaps * sizeof(double * ));
  double * BetavecPREV = malloc( * p_Nhaps * sizeof(double));
  double * BetavecCURRENT = malloc( * p_Nhaps * sizeof(double));
  double Alphasum, Alphasumnew, Betasum, Betasumnew;
  double large_num;
  double * copy_prob_new = malloc( * p_Nhaps * sizeof(double));
  double * copy_prob_newSTART = malloc( * p_Nhaps * sizeof(double));
  double * Alphasumvec = malloc( * p_Nloci * sizeof(double));
  double * expected_transition_prob = malloc(( * p_Nloci - 1) * sizeof(double));
  double * corrected_chunk_count = malloc( * p_Nhaps * sizeof(double));
  double * regional_chunk_count = malloc( * p_Nhaps * sizeof(double));
  double * expected_chunk_length = malloc( * p_Nhaps * sizeof(double));
  double * expected_differences = malloc( * p_Nhaps * sizeof(double));
  double * regional_chunk_count_sum = malloc(ndonorpops * sizeof(double));
  double * regional_chunk_count_sum_final = malloc(ndonorpops * sizeof(double));
  double * regional_chunk_count_sum_squared_final = malloc(ndonorpops * sizeof(double));
  double * ind_snp_sum_vec = malloc(ndonorpops * sizeof(double));
  double * snp_info_measure = malloc(ndonorpops * sizeof(double));
  double * exp_copy_pop = malloc(ndonorpops * sizeof(double));
  double expected_chunk_length_sum, sum_prob;

  double ** copy_prob_new_mat = malloc(8 * sizeof(double * ));
  for (i = 0; i < 8; i++)
    copy_prob_new_mat[i] = malloc(( * p_Nhaps + 1) * sizeof(double));

  for (i = 0; i < * p_Nhaps; i++) {
    Alphamat[i] = malloc( * p_Nloci * sizeof(double));
  }
  for (i = 0; i < ndonorpops; i++) {
    regional_chunk_count_sum[i] = 0.0;
    regional_chunk_count_sum_final[i] = 0.0;
    regional_chunk_count_sum_squared_final[i] = 0.0;
    ind_snp_sum_vec[i] = 0.0;
    snp_info_measure[i] = 0.0;
  }

  // Theta as given in Li and Stephens 
  sum = 0;
  for (i = 1; i < * p_nchr; i++) {
    sum = sum + 1.0 / i;
  }
  Theta = 1.0 / sum;

  for (i = 0; i < * p_Nhaps; i++) {
    if (MutProb_vec[i] < 0) MutProb_vec[i] = 0.5 * Theta / ( * p_Nhaps + Theta);
    //if (MutProb_vec[i]<rounding_val) MutProb_vec[i]=rounding_val;
  }

  // TransProb[i] is probability of copying mechanism "jumping" between
  //   loci i and i+1 
  delta = 1.0;

  if (unlinked_ind == 0 && lambda[0] >= 0) TransProb[0] = 1 - exp(-1 * (pos[1] - pos[0]) * delta * p_rhobar * lambda[0]);
  if (unlinked_ind == 1 || lambda[0] < 0) TransProb[0] = 1.0;
  /*
  if (TransProb[0]<rounding_val)
  {
    printf("Transition prob is too low; will likely cause rounding errors. Exiting...\n");
    exit(1);
  }
  */

  for (locus = 1; locus < * p_Nloci - 1; locus++) {
    delta = 1.0;
    if (unlinked_ind == 0 && lambda[locus] >= 0) TransProb[locus] = 1 - exp(-1 * (pos[locus + 1] - pos[locus]) * delta * p_rhobar * lambda[locus]);
    if (unlinked_ind == 1 || lambda[locus] < 0) TransProb[locus] = 1.0;
    /*
      if (TransProb[locus]<rounding_val)
  {
    printf("Transition prob is too low; will likely cause rounding errors. Exiting...\n");
    exit(1);
  }
      */
  }
  /* FORWARDS ALGORITHM: (Rabiner 1989, p.262) */
  /* INITIALIZATION: */
  Alphasum = 0.0;
  for (i = 0; i < * p_Nhaps; i++) {
    //if (newh[0]!=9 && existing_h[i][0] != -9) ObsStateProb = (1-MutProb_vec[i]) * (newh[0] == existing_h[i][0]) + MutProb_vec[i] * (newh[0] != existing_h[i][0]);
    if (newh[0] != 9) ObsStateProb = (1 - MutProb_vec[i]) * (newh[0] == existing_h[i][0]) + MutProb_vec[i] * (newh[0] != existing_h[i][0]);
    if (newh[0] == 9) ObsStateProb = 1.0;
    //if (existing_h[i][0]==9) ObsStateProb = small_missing_val;
    /*
      if (ObsStateProb<rounding_val)
  {
    printf("Mutation (emission) rate is too low; will likely cause rounding errors. Exiting...\n");
    exit(1);
  }
      */

    Alphamat[i][0] = log(copy_probSTART[i] * ObsStateProb);
    Alphasum = Alphasum + exp(Alphamat[i][0]) * TransProb[0];
  }

  /* INDUCTION: */
  Alphasum = log(Alphasum);
  for (locus = 1; locus < * p_Nloci; locus++) {
    Alphasumnew = 0.0;
    large_num = -1.0 * Alphasum;
    for (i = 0; i < * p_Nhaps; i++) {
      if (newh[locus] != 9) ObsStateProb = (1 - MutProb_vec[i]) * (newh[locus] == existing_h[i][locus]) + MutProb_vec[i] * (newh[locus] != existing_h[i][locus]);
      if (newh[locus] == 9) ObsStateProb = 1.0;
      /*
      if (ObsStateProb<rounding_val)
        {
          printf("Mutation (emission) rate is too low; will likely cause rounding errors. Exiting...\n");
          exit(1);
        }
      */

      Alphamat[i][locus] = log(ObsStateProb * copy_prob[i] * exp(Alphasum + large_num) + ObsStateProb * (1 - TransProb[(locus - 1)]) * exp(Alphamat[i][(locus - 1)] + large_num)) - large_num;
      if (locus < ( * p_Nloci - 1)) Alphasumnew = Alphasumnew + exp(Alphamat[i][locus] + large_num) * TransProb[locus];
      if (locus == ( * p_Nloci - 1)) Alphasumnew = Alphasumnew + exp(Alphamat[i][locus] + large_num);
    }
    Alphasum = log(Alphasumnew) - large_num;
  }
  //if (Alphasum == (Alphasum*5))
  if (isnan(Alphasum)) {
    printf("Negative or NaN likelihood. Could be because emission or transition probabilities are too low??...Exiting...\n");
    //for (i=0; i < *p_Nhaps; i++) printf("%d %lf %lf %lf\n",i,copy_prob[i],log(copy_prob[i]),log(MutProb_vec[i]));
    exit(1);
  }
  fprintf(fout3, " %.10lf", Alphasum);

  for (i = 0; i < * p_Nhaps; i++) {
    copy_prob_new[i] = 0.0;
    corrected_chunk_count[i] = 0.0;
    expected_chunk_length[i] = 0.0;
    expected_differences[i] = 0.0;
    regional_chunk_count[i] = 0.0;
  }
  total_regional_chunk_count = 0.0;
  num_regions = 0;
  if (run_num <= (run_thres - 1)) {
    /* BACKWARDS ALGORITHM: (Rabiner 1989, p.263) */
    /* INITIALIZATION: */
    Betasum = 0.0;
    if (run_num == (run_thres - 1)) {
      for (i = 0; i < ndonorpops; i++)
        exp_copy_pop[i] = 0.0;
    }
    for (i = 0; i < * p_Nhaps; i++) {
      if (newh[( * p_Nloci - 1)] != 9) ObsStateProb = (1 - MutProb_vec[i]) * (newh[( * p_Nloci - 1)] == existing_h[i][( * p_Nloci - 1)]) + MutProb_vec[i] * (newh[( * p_Nloci - 1)] != existing_h[i][( * p_Nloci - 1)]);
      if (newh[( * p_Nloci - 1)] == 9) ObsStateProb = 1.0;
      /*
      if (ObsStateProb<rounding_val)
        {
          printf("Mutation (emission) rate is too low; will likely cause rounding errors. Exiting...\n");
          exit(1);
        }
      */
      BetavecPREV[i] = 0.0;
      Betasum = Betasum + TransProb[( * p_Nloci - 2)] * copy_prob[i] * ObsStateProb * exp(BetavecPREV[i]);
      if (run_num == (run_thres - 1)) exp_copy_pop[pop_vec[i]] = exp_copy_pop[pop_vec[i]] + exp(BetavecPREV[i] + Alphamat[i][( * p_Nloci - 1)] - Alphasum);
      // for estimating new mutation rates:
      expected_differences[i] = expected_differences[i] + exp(Alphamat[i][( * p_Nloci - 1)] - Alphasum) * (newh[( * p_Nloci - 1)] != existing_h[i][( * p_Nloci - 1)]);
    }
    if ((run_num == (run_thres - 1)) && (print_file9_ind == 1)) {
      gzprintf(fout9, "%.0lf", pos[( * p_Nloci - 1)]);
      if (all_versus_all_ind == 0) {
        for (i = 0; i < ndonorpops; i++)
          gzprintf(fout9, " %lf", exp_copy_pop[i]);
      }
      if (all_versus_all_ind == 1) {
        for (i = 0; i < ndonorpops; i++) {
          if (i == ind_val) gzprintf(fout9, " 0.0");
          gzprintf(fout9, " %lf", exp_copy_pop[i]);
        }
        if (ind_val == ndonorpops) gzprintf(fout9, " 0.0");
      }
      gzprintf(fout9, "\n");
    }

    /* INDUCTION: */
    Betasum = log(Betasum);
    /* CALCULATE EXPECTED NUMBER OF TIMES OF COPYING TO EACH DONOR POP (Rabiner 1989, p.263,265 or Scheet/Stephens 2006 Appendix C): */
    for (locus = ( * p_Nloci - 2); locus >= 0; locus--) {
      Betasumnew = 0.0;
      large_num = -1.0 * Betasum;
      total_prob = 0.0;

      constant_exclude_i = 0.5;
      constant_from_i_to_i = 1.0;
      constant_exclude_i_both_sides = 0.0;
      expected_chunk_length_sum = 0.0;
      sum_prob = 0.0;
      if (run_num == (run_thres - 1)) {
        for (i = 0; i < ndonorpops; i++)
          exp_copy_pop[i] = 0.0;
      }
      for (i = 0; i < * p_Nhaps; i++) {
        if (newh[locus] != 9) ObsStateProb = (1 - MutProb_vec[i]) * (newh[locus] == existing_h[i][locus]) + MutProb_vec[i] * (newh[locus] != existing_h[i][locus]);
        if (newh[locus] == 9) ObsStateProb = 1.0;
        if (newh[(locus + 1)] != 9) ObsStateProbPREV = (1 - MutProb_vec[i]) * (newh[(locus + 1)] == existing_h[i][(locus + 1)]) + MutProb_vec[i] * (newh[(locus + 1)] != existing_h[i][(locus + 1)]);
        if (newh[(locus + 1)] == 9) ObsStateProbPREV = 1.0;
        /*
        if ((ObsStateProb<rounding_val) || (ObsStateProbPREV<rounding_val))
    {
      printf("Mutation (emission) rate is too low; will likely cause rounding errors. Exiting...\n");
      exit(1);
    }
        */
        BetavecCURRENT[i] = log(exp(Betasum + large_num) + (1 - TransProb[locus]) * ObsStateProbPREV * exp(BetavecPREV[i] + large_num)) - large_num;
        if (locus > 0) Betasumnew = Betasumnew + TransProb[(locus - 1)] * copy_prob[i] * ObsStateProb * exp(BetavecCURRENT[i] + large_num);
        if (locus == 0) copy_prob_newSTART[i] = exp(Alphamat[i][0] + BetavecCURRENT[i] - Alphasum);
        total_prob = total_prob + exp(Alphamat[i][(locus + 1)] + BetavecPREV[i] - Alphasum) - exp(Alphamat[i][locus] + BetavecPREV[i] - Alphasum) * ObsStateProbPREV * (1 - TransProb[locus]);

        copy_prob_new[i] = copy_prob_new[i] + exp(Alphamat[i][(locus + 1)] + BetavecPREV[i] - Alphasum) - exp(Alphamat[i][locus] + BetavecPREV[i] - Alphasum) * ObsStateProbPREV * (1 - TransProb[locus]);

        total_prob_from_i_to_i = exp(Alphamat[i][locus] + BetavecPREV[i] - Alphasum) * ObsStateProbPREV * (1 - TransProb[locus] + TransProb[locus] * copy_prob[i]);
        total_prob_to_i_exclude_i = exp(Alphamat[i][(locus + 1)] + BetavecPREV[i] - Alphasum) - exp(Alphamat[i][locus] + BetavecPREV[i] - Alphasum) * ObsStateProbPREV * (1 - TransProb[locus] + TransProb[locus] * copy_prob[i]);
        total_prob_from_i_exclude_i = exp(Alphamat[i][locus] + BetavecCURRENT[i] - Alphasum) - exp(Alphamat[i][locus] + BetavecPREV[i] - Alphasum) * ObsStateProbPREV * (1 - TransProb[locus] + TransProb[locus] * copy_prob[i]);
        total_prob_from_any_to_any_exclude_i = 1.0 - exp(Alphamat[i][locus] + BetavecCURRENT[i] - Alphasum) - exp(Alphamat[i][(locus + 1)] + BetavecPREV[i] - Alphasum) + exp(Alphamat[i][locus] + BetavecPREV[i] - Alphasum) * ObsStateProbPREV * (1 - TransProb[locus] + TransProb[locus] * copy_prob[i]);

        regional_chunk_count[i] = regional_chunk_count[i] + (exp(Alphamat[i][(locus + 1)] + BetavecPREV[i] - Alphasum) - exp(Alphamat[i][locus] + BetavecPREV[i] - Alphasum) * ObsStateProbPREV * (1 - TransProb[locus]));
        total_regional_chunk_count = total_regional_chunk_count + (exp(Alphamat[i][(locus + 1)] + BetavecPREV[i] - Alphasum) - exp(Alphamat[i][locus] + BetavecPREV[i] - Alphasum) * ObsStateProbPREV * (1 - TransProb[locus]));
        ind_snp_sum_vec[pop_vec[i]] = ind_snp_sum_vec[pop_vec[i]] + (exp(Alphamat[i][(locus + 1)] + BetavecPREV[i] - Alphasum) - exp(Alphamat[i][locus] + BetavecPREV[i] - Alphasum) * ObsStateProbPREV * (1 - TransProb[locus]));

        //corrected_chunk_count[i]=corrected_chunk_count[i]+(exp(Alphamat[i][(locus+1)]+BetavecPREV[i]-Alphasum)-exp(Alphamat[i][locus]+BetavecPREV[i]- Alphasum)*ObsStateProbPREV*(1-TransProb[locus]))*(1.0+(1.0/(*p_Nhaps))*((p_rhobar*(pos[(locus+1)]-pos[locus])*delta*lambda[locus]/(*p_Nhaps))/(1.0-exp(-1.0*p_rhobar*(pos[(locus+1)]-pos[locus])*delta*lambda[locus]/(*p_Nhaps)))-1.0));
        corrected_chunk_count[i] = corrected_chunk_count[i] + (exp(Alphamat[i][(locus + 1)] + BetavecPREV[i] - Alphasum) - exp(Alphamat[i][locus] + BetavecPREV[i] - Alphasum) * ObsStateProbPREV * (1 - TransProb[locus]));
        if (unlinked_ind == 0 && lambda[locus] >= 0) expected_chunk_length[i] = expected_chunk_length[i] + 100 * (pos[locus + 1] - pos[locus]) * delta * lambda[locus] * (constant_from_i_to_i * total_prob_from_i_to_i + constant_exclude_i * (total_prob_to_i_exclude_i + total_prob_from_i_exclude_i) + constant_exclude_i_both_sides * total_prob_from_any_to_any_exclude_i); // multiply by 100 to get cM
        expected_chunk_length_sum = expected_chunk_length_sum + constant_from_i_to_i * total_prob_from_i_to_i + constant_exclude_i * (total_prob_to_i_exclude_i + total_prob_from_i_exclude_i) + constant_exclude_i_both_sides * total_prob_from_any_to_any_exclude_i;
        // for estimating new mutation rates:
        expected_differences[i] = expected_differences[i] + exp(Alphamat[i][locus] + BetavecCURRENT[i] - Alphasum) * (newh[locus] != existing_h[i][locus]);
        BetavecPREV[i] = BetavecCURRENT[i];

        if (run_num == (run_thres - 1)) exp_copy_pop[pop_vec[i]] = exp_copy_pop[pop_vec[i]] + exp(BetavecCURRENT[i] + Alphamat[i][locus] - Alphasum);

        sum_prob = sum_prob + total_prob_from_i_to_i + total_prob_to_i_exclude_i + total_prob_from_i_exclude_i;
      }

      if ((run_num == (run_thres - 1)) && (print_file9_ind == 1)) {
        gzprintf(fout9, "%.0lf", pos[locus]);
        if (all_versus_all_ind == 0) {
          for (i = 0; i < ndonorpops; i++)
            gzprintf(fout9, " %lf", exp_copy_pop[i]);
        }
        if (all_versus_all_ind == 1) {
          for (i = 0; i < ndonorpops; i++) {
            if (i == ind_val) gzprintf(fout9, " 0.0");
            gzprintf(fout9, " %lf", exp_copy_pop[i]);
          }
          if (ind_val == ndonorpops) gzprintf(fout9, " 0.0");
        }
        gzprintf(fout9, "\n");
      }

      expected_transition_prob[locus] = total_prob;
      if (locus > 0) Betasum = log(Betasumnew) - large_num;

      if ((total_regional_chunk_count + rounding_val) >= region_size) {
        for (i = 0; i < * p_Nhaps; i++) {
          regional_chunk_count_sum[pop_vec[i]] = regional_chunk_count_sum[pop_vec[i]] + regional_chunk_count[i];
          regional_chunk_count[i] = 0.0;
        }
        for (i = 0; i < ndonorpops; i++) {
          regional_chunk_count_sum_final[i] = regional_chunk_count_sum_final[i] + regional_chunk_count_sum[i];
          regional_chunk_count_sum_squared_final[i] = regional_chunk_count_sum_squared_final[i] + pow(regional_chunk_count_sum[i], 2.0);
          regional_chunk_count_sum[i] = 0.0;
        }
        total_regional_chunk_count = 0.0;
        num_regions = num_regions + 1;
      }
      total_ind_sum = 0.0;
      for (i = 0; i < ndonorpops; i++)
        total_ind_sum = total_ind_sum + ind_snp_sum_vec[i];
      for (i = 0; i < ndonorpops; i++) {
        snp_info_measure[i] = snp_info_measure[i] + pow((ind_snp_sum_vec[i] / total_ind_sum), 2.0);
        ind_snp_sum_vec[i] = 0.0;
      }
    }
    for (i = 0; i < ndonorpops; i++)
      snp_info_measure[i] = snp_info_measure[i] / ( * p_Nloci);

    /* CALCULATE EXPECTED NUMBER OF TOTAL TRANSITIONS, IN ORDER TO ESTIMATE N_e (Scheet/Stephens 2006 Appendix C (C3)): */
    total_prob = 0.0;
    total_gen_dist = 0.0;
    for (locus = 0; locus < ( * p_Nloci - 1); locus++) {
      if (unlinked_ind == 0 && lambda[locus] >= 0) total_gen_dist = total_gen_dist + (pos[(locus + 1)] - pos[locus]) * delta * lambda[locus];
      if (unlinked_ind == 0 && lambda[locus] >= 0) total_prob = total_prob + ((p_rhobar * (pos[(locus + 1)] - pos[locus]) * delta * lambda[locus]) / (1.0 - exp(-1.0 * p_rhobar * (pos[(locus + 1)] - pos[locus]) * delta * lambda[locus]))) * expected_transition_prob[locus];
    }
    if (unlinked_ind == 0) N_e_new = total_prob / total_gen_dist;
    if (unlinked_ind == 1) N_e_new = 0.0;

    /* CALCULATE SOMETHING ANALAGOUS TO EXPECTED NUMBER OF TIMES EACH HAP i IS VISITED, CONDITIONAL ON THE OBSERVED DATA (I.E  (27) AND PARAGRAPH UNDER (38) IN RABINER 1989, Proceedings of the IEEE 77(2):257-286), BUT -- AS WE'RE ONLY COUNTING CHUNKS -- SUBTRACT OUT TIMES YOU DO NOT SWITCH */
    for (i = 0; i < * p_Nhaps; i++) corrected_chunk_count[i] = corrected_chunk_count[i] + copy_prob_newSTART[i];
  }

  /* print-out samples if we've done enough iterations: */
  if (run_num == (run_thres - 1)) {

    N_e_new = p_rhobar;

    for (i = 0; i < * p_Nhaps; i++) {
      copy_prob_new[i] = copy_prob[i];
      copy_prob_newSTART[i] = copy_probSTART[i];
    }

    /* calculate Alphasums (for efficient sampling): */
    for (locus = 0; locus < * p_Nloci; locus++) {
      Alphasumvec[locus] = 0.0;
      large_num = Alphamat[0][locus];
      for (i = 1; i < * p_Nhaps; i++) {
        if (Alphamat[i][locus] > large_num)
          large_num = Alphamat[i][locus];
      }
      large_num = -1.0 * large_num;
      for (i = 0; i < * p_Nhaps; i++)
        Alphasumvec[locus] = Alphasumvec[locus] + exp(Alphamat[i][locus] + large_num);
      Alphasumvec[locus] = log(Alphasumvec[locus]) - large_num;
    }

    /* SAMPLING ALGORITHM: (from Falush, Stephens, & Pritchard (2003) Genetics 164:1567-1587) */
    for (j = 0; j < nsampTOT; j++) {
      /* sample last position: */
      total_prob = 0.0;
      large_num = Alphamat[0][( * p_Nloci - 1)];
      for (i = 1; i < * p_Nhaps; i++) {
        if (Alphamat[i][( * p_Nloci - 1)] > large_num)
          large_num = Alphamat[i][( * p_Nloci - 1)];
      }
      large_num = -1.0 * large_num;
      random_unif = (double) rand() / RAND_MAX;
      total_prob = Alphasumvec[( * p_Nloci - 1)];
      prob = 0.0;
      for (i = 0; i < * p_Nhaps; i++) {
        prob = prob + exp(Alphamat[i][( * p_Nloci - 1)] + large_num);
        if (random_unif <= exp(log(prob) - large_num - total_prob)) {
          sample_state[( * p_Nloci - 1)] = i;
          break;
        }
      }

      /* sample remaining positions: */
      for (locus = ( * p_Nloci - 2); locus >= 0; locus--) {
        // first sample prob you switch and see if you need to
        // if you do need to switch, you need to go through the below loop to figure out where to switch to
        large_num = -1.0 * Alphasumvec[locus];
        total_prob = log(exp(Alphasumvec[locus] + large_num) * TransProb[locus] * copy_prob[sample_state[(locus + 1)]] + exp(Alphamat[sample_state[(locus + 1)]][locus] + large_num) * (1.0 - TransProb[locus])) - large_num;
        no_switch_prob = exp(log(exp(Alphamat[sample_state[(locus + 1)]][locus] + large_num) * (1.0 - TransProb[locus])) - large_num - total_prob);
        random_unifSWITCH = (double) rand() / RAND_MAX;
        if (random_unifSWITCH <= no_switch_prob) sample_state[locus] = sample_state[(locus + 1)];

        if (random_unifSWITCH > no_switch_prob) {
          total_prob = 0.0;
          large_num = Alphamat[0][locus];
          for (i = 1; i < * p_Nhaps; i++) {
            if (Alphamat[i][locus] > large_num)
              large_num = Alphamat[i][locus];
          }
          large_num = -1.0 * large_num;

          random_unif = (double) rand() / RAND_MAX;
          total_prob = log(exp(Alphasumvec[locus] + large_num) * TransProb[locus] * copy_prob[sample_state[(locus + 1)]]) - large_num;
          prob = 0.0;
          for (i = 0; i < * p_Nhaps; i++) {
            prob = prob + exp(Alphamat[i][locus] + large_num) * TransProb[locus] * copy_prob[sample_state[(locus + 1)]];
            if (random_unif <= exp(log(prob) - large_num - total_prob)) {
              sample_state[locus] = i;
              break;
            }
          }
        }
      }

      fprintf(fout, "%d", j + 1);
      for (i = 0; i < * p_Nloci; i++)
        fprintf(fout, " %d", hap_label_vec[sample_state[i]]);
      fprintf(fout, "\n");
    }
  }

  for (i = 0; i < * p_Nhaps; i++)
    copy_prob_new_mat[0][i] = copy_prob_new[i];
  for (i = 0; i < * p_Nhaps; i++)
    copy_prob_new_mat[1][i] = copy_prob_newSTART[i];
  for (i = 0; i < * p_Nhaps; i++)
    copy_prob_new_mat[2][i] = corrected_chunk_count[i];
  for (i = 0; i < * p_Nhaps; i++)
    copy_prob_new_mat[3][i] = expected_chunk_length[i];
  for (i = 0; i < * p_Nhaps; i++)
    copy_prob_new_mat[4][i] = expected_differences[i];
  for (i = 0; i < ndonorpops; i++)
    copy_prob_new_mat[5][i] = regional_chunk_count_sum_final[i];
  for (i = 0; i < ndonorpops; i++)
    copy_prob_new_mat[6][i] = regional_chunk_count_sum_squared_final[i];
  for (i = 0; i < ndonorpops; i++)
    copy_prob_new_mat[7][i] = snp_info_measure[i];
  copy_prob_new_mat[0][( * p_Nhaps)] = N_e_new;
  copy_prob_new_mat[1][( * p_Nhaps)] = num_regions;

  for (i = 0; i < * p_Nhaps; i++) {
    free(Alphamat[i]);
  }
  free(Alphamat);
  free(BetavecPREV);
  free(BetavecCURRENT);
  free(TransProb);
  free(sample_state);
  free(expected_transition_prob);
  free(Alphasumvec);
  free(copy_prob_new);
  free(copy_prob_newSTART);
  free(corrected_chunk_count);
  free(regional_chunk_count);
  free(expected_chunk_length);
  free(expected_differences);
  free(regional_chunk_count_sum);
  free(regional_chunk_count_sum_final);
  free(regional_chunk_count_sum_squared_final);
  free(ind_snp_sum_vec);
  free(snp_info_measure);
  free(exp_copy_pop);

  return (copy_prob_new_mat);
}