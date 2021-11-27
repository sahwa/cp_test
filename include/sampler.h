#include "zlib.h"

double ** sampler(int * newh, int ** existing_h, int * p_Nloci, int * p_Nhaps, int * p_nchr, double p_rhobar, double * MutProb_vec, int * allelic_type_count_vec, double * lambda, double * pos, int nsampTOT, double * copy_prob, double * copy_probSTART, int * pop_vec, int ndonorpops, double region_size, int run_num, int run_thres, int all_versus_all_ind, int haploid_ind, int unlinked_ind, int ind_val, int * hap_label_vec, int print_file9_ind, FILE * fout, FILE * fout3, gzFile fout9);
