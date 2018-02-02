#include "read_count.h"
#include "matrix.h"

int count_reads(char * fn)
{

  struct read_pair{
      bam1_t * read1;
      bam1_t * read2;
  }rp;


  htsFile * h = hts_open(fn, "r");



  const htsFormat *fmt = hts_get_format(h);

  int hts_close(htsFile *fp);

  fprintf(stderr, "INFO: detected %s filetype\n", hts_format_file_extension(fmt));

  samFile *in = 0;
  bam_hdr_t *header = NULL;

  if((in = sam_open_format(fn, "r", fmt)) == 0){
      printf("FATAL: failed to open \"%s\" for reading\n",fn);
      return 1;
  }

  if ((header = sam_hdr_read(in)) == 0) {
      fprintf(stderr, "FATAL: failed to read the header \"%s\" \n", fn);
      return 1;
  }

    datum * counts = (datum *) malloc(header->n_targets * sizeof(datum));
    memset(counts, 0, header->n_targets*sizeof(datum));

    rp.read1 = bam_init1();
    rp.read2 = bam_init1();

    int r1;
  	int r2;
  	long int c = 0;

    while (1){

      c+=1;

      r1 = sam_read1(in, header, rp.read1);
      r2 = sam_read1(in, header, rp.read2);

      if( r1 < 0 || r2 < 0 ) break;

      counts[rp.read1->core.tid] += 1;
      counts[rp.read2->core.tid] += 1;

    }

    datum i;
    fprintf(stdout, "##seqid\tlength\tcount_reads\n");
    for(i = 0; i < header->n_targets; i++){
      fprintf(stdout, "%s\t%i\t%i\n", header->target_name[i], header->target_len[i], counts[i]);
    }

    bam_destroy1(rp.read1);
    bam_destroy1(rp.read2);


    return 0;

}
