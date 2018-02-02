
#ifndef READ_COUNT_H
#define READ_COUNT_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <getopt.h>

#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"

/**
 * Counts the number of reads per seqid then prints to stdout
 * @param  fn file pointer
 * @return    0 on success 
 */
int count_reads(char * fn);

#endif
