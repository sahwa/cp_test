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

struct data_t * ReadData(FILE * fd, int ind_val, int ploidy, int * include_ind_vec, char ** pop_label_vec, int num_donor_pops, char ** donor_pop_vec, int * pop_vec_tot, double * copy_prob_tot, double * copy_probSTART_tot, double * MutProb_vec_tot, int * ndonorhaps_tot, int all_versus_all_ind);
