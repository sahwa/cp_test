#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void usage() {
  printf("to run: use './ChromoPainterv2' with following options:\n");
  printf("       -g <geno.filein>  (REQUIRED; no default)\n");
  printf("       -r <recommap.filein>  (REQUIRED; no default -- unless using -u switch)\n");
  printf("       -t <labels.filein> file listing id and population labels for each individual (REQUIRED; no default -- unless using -a switch)\n");
  printf("       -f <poplist.filein> <f_1> <f_2>  file listing donor and recipient populations (REQUIRED; no default -- unless using -a switch); paint recipient individuals f_1 through f_2 using all donor population haplotypes (use '-f <donorlist.filein> 0 0' to paint all recipient inds)\n");
  printf("       -i <int>  number of EM iterations for estimating parameters (default=0)\n");
  printf("       -in  maximize over average switch rate parameter using E-M\n");
  printf("       -ip  maximize over copying proportions using E-M\n");
  printf("       -im  maximize over donor population mutation (emission) probabilities using E-M\n");
  printf("       -iM  maximize over global mutation (emission) probability using E-M\n");
  printf("       -s <int>  number of samples per recipient haplotype (default=10)\n");
  printf("       -n <double>  average switch rate parameter start-value (default=400000 divided by total number of donor haplotypes included in analysis)\n");
  printf("       -p  specify to use prior copying probabilities in population list (-f) file\n");
  printf("       -m  specify to use donor population mutation (emission) probabilities in population list (-f) file\n");
  printf("       -M <double>  global mutation (emission) probability (default=Li & Stephen's (2003) fixed estimate)\n");
  printf("       -k <double>  specify number of expected chunks to define a 'region' (default=100)\n");
  printf("       -j  specify that individuals are haploid\n");
  printf("       -u  specify that data are unlinked\n");
  printf("       -a <a_1> <a_2>  paint individuals a_1 through a_2 using every other individual (use '-a 0 0' to paint all inds)\n");
  printf("       -b  print-out zipped file with suffix '.copyprobsperlocus.out' containing prob each recipient copies each donor at every SNP (note: file can be quite large)\n");
  printf("       -o <outfile-prefix>  (default = 'geno.filein')\n");
  printf("       --help  print this menu\n");
}
