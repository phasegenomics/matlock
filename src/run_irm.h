#ifndef RUN_IRM_H
#define RUN_IRM_H


#include <assert.h>
#include <getopt.h>

#include "matrix.h"
#include "irm.h"
#include "count_motif.h"

#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"


struct options{
  int  nclusters;
  int    nsweeps;
  char  *  fasta;
  char  * binmat;
  char  *  motif;

}global_opts;




/**
 * Runs the weighted IRM
 * @param  argv - command line options
 * @param  argc - command line index
 * @return       > 0 if problem
 */
int run_irm(char ** argv, int argc );

/**
 * Print the index and value of the vector
 * @param v - the vector
 * @param l - 0-l
 */
void printv(double * v, datum l);

/**
 * Perform one sweep of the Gibb's sampler
 * @param  mm the model
 * @return   int the number of switches for each iteration
 */
int gibbs_sweep(struct model * mm, int nsweep);

/**
 * Parses the command line options into the global_opts
 * @param  argc [description]
 * @param  argv [description]
 * @return      [description]
 */
int parse_command_line(char ** argv, int argc);

/**
 * [print_irm_runner_usage description]
 */
void print_irm_runner_usage(void);


int normalize_link_weights_by_seq_len(struct model * mm, struct sequenceInfo * sInfo);



double modularity(struct model *mm, double matrix_sum);


#endif /* RUN_IRM_H */
