#include <stdio.h>
#include <stdlib.h>
#include "data_t.h"

void DestroyData(struct data_t * dat, int ploidy) {

  int n_args = 8;

  for (int i = 0; i < dat -> ndonorhaps; i++) {
    free(dat -> cond_chromosomes[i]);
  }

  for (int i = 0; i < ploidy; i++) {
    free(dat -> ind_chromosomes[i]);
  }

  free(dat -> cond_chromosomes);
  free(dat -> ind_chromosomes);
  free(dat -> positions);
  free(dat -> lambda);
  free(dat -> copy_prob);
  free(dat -> copy_probSTART);
  free(dat -> MutProb_vec);
  free(dat -> pop_vec);
  free(dat -> hap_label_vec);
  free(dat -> copy_prob_new);
  free(dat -> copy_prob_newSTART);
  free(dat -> MutProb_vec_new);

  for (int i = 0; i < n_args; i++) {
    free(dat -> back_prob[i]);
  }

  free(dat -> back_prob);
  free(dat -> ndonorhaps_vec);
}