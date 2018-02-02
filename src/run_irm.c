#include "run_irm.h"

#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })


void print_irm_runner_usage(void){
  	fprintf(stderr, "\nusage: matlock irm [options] -n 25 -f your.fasta -b your.binmat \n\n");
    fprintf(stderr, "   + required:\n");
    fprintf(stderr, "      -c number of clusters [int].\n");
    fprintf(stderr, "      -f fasta file.\n");
    fprintf(stderr, "      -b binmat file.\n");
    fprintf(stderr, "   + options: \n");
  	fprintf(stderr, "      -s number of sweeps [int].\n");
}

int parse_command_line(char ** argv, int argc)
{

  global_opts.nclusters  = -1;
  global_opts.fasta      = NULL;
  global_opts.binmat     = NULL;
  global_opts.nsweeps    = 100;
  global_opts.motif      = NULL;

  if(argc < 2){
    print_irm_runner_usage();
    return 1;
  }

  int             c;
   const char    * short_opt = "hf:b:s:c:m:";
   struct option   long_opt[] =
   {
      {"help",          no_argument,       NULL, 'h'},
      {"fasta",         required_argument, NULL, 'f'},
      {"bin",           required_argument, NULL, 'b'},
      {"cluster",       required_argument, NULL, 'c'},
      {"sweep",         optional_argument, NULL, 'w'},
      {"motif",         required_argument, NULL, 'm'},
      {NULL,            0,                 NULL, 0  }
   };

   while((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1)
   {
      switch(c)
        {
          case 'h':
          {
            print_irm_runner_usage();
            exit(1);
          }
          case 'm':
          {
            global_opts.motif = optarg;
            break;
          }
          case 'f':
          {
            global_opts.fasta = optarg;
            break;
          }
          case 'b':
          {
            global_opts.binmat = optarg;
            break;
          }
          case 'c':
          {
            global_opts.nclusters = atoi(optarg);
            break;
          }

          case 's':
            {
              global_opts.nsweeps = atoi(optarg);
              if(global_opts.nsweeps < 1){
                print_irm_runner_usage();
                exit(1);
              }
              break;
            }
          default:
          {
            /* invalid option */
            fprintf(stderr, "%s: option '-%c' is invalid: ignored\n",
            argv[0], optopt);
            break;
          }
        }
    }


  fprintf(stderr, "INFO: RUNNING: ");
  int i = 0;

   for(;i<argc;i++)
   {
       fprintf(stderr, "%s ", argv[i]);
   }
  fprintf(stderr, "\n");
  return 0;
}


void printv(double * v, datum l)
{
  datum i = 0;
  for(; i < l; i++){
      fprintf(stderr, "%i: %.1f ", i, v[i]);
  }
  fprintf(stderr, "\n");
}


int gibbs_sweep(struct model * mm, int nsweep){


  datum i   = 0;
  datum j   = 0;
  datum k   = 0;
  datum l   = 0;
  datum r   = 0;

  datum x   = 0;
  datum y   = 0;
  datum z   = 0;

  double ms =  matrix_sum(mm->adj);

  double new_score = hot_cold(mm);
  double old_score = hot_cold(mm);

  int nswitch  = 0;

  int counter = 0;

  for(; i < mm->nd; i++){

    if(mm->rowsum[i] == 0) continue;

    x = random_next(mm, i );
    y = random_next(mm, x );

    if(i == y){

      fprintf(stderr, "HERE x: %i y: %i z: %i\n", i, x, y);

      if(mm->labels[x] != mm->labels[i]){
        mm->nswitches[x]+=1;
      }

        mm->labels[x] = mm->labels[i];
    }

  }

  double mod = modularity(mm,ms);
  fprintf(stderr, "modularity %f\n", mod );


  return nswitch;
}

int run_irm(char ** argv, int argc)
{

  int parse_flag = parse_command_line(argv, argc);
  if(parse_flag > 0) return 1;

  // load the sequence info
  struct sequenceInfo * sInfo = load_seq_info(global_opts.fasta, global_opts.motif);
  if(sInfo == NULL) return 1;

  // read the weighted undirected graph into memory
  struct matrix * link_counts = thaw_matrix(global_opts.binmat);

  // create the model
  struct model * mm =  model_init(&link_counts, 1, 1, 1);

  // load grop assignments just cut them into even bins
  init_group_assignment_fixed(mm, global_opts.nclusters);

  normalize_link_weights_by_seq_len(mm, sInfo);

  row_sum(mm);


   // count the links
   int r2 =  count_links(mm);
   assert(r2 == 0);

   assert( check_label_count_weight(mm) == 0);
   assert( check_label_count_sum(mm) == 0);

   int sweepi  = 1;

   for(; sweepi <= (int)global_opts.nsweeps; sweepi++){

       int nswitch = gibbs_sweep(mm, sweepi);
       fprintf(stderr, "INFO: sweep %i had %i switches\n", sweepi, nswitch);
  //     printv(mm->label_counts, mm->next_open +1);
       if(nswitch < 0) break;
  //     printv(mm->cluster_weights, mm->next_open +1);

   }

   print_results(mm);


  return 0;

}



int normalize_link_weights_by_seq_len(struct model * mm, struct sequenceInfo * sinfo)
{

  datum i = 0;
  datum j = 0;

  uint64_t x = 0;
  uint64_t y = 0;
  uint64_t b = 0;

  uint64_t l1 = 0;
  uint64_t l2 = 0;

  fprintf(stderr, "%i\n", matrix_sum(mm->adj));

  for(; i < mm->adj->n1; i++){
    j = 0;
    for(; j < mm->adj->n1; j++){
      if(mm->adj->dat[i][j] == 0) continue;


      b = mm->adj->dat[i][j];

      l1 = mm->seq_lens[i];
      l2 = mm->seq_lens[j];

  if(i == j ){ mm->adj->dat[i][j] = 0; continue; }

     mm->adj->dat[i][j] = mm->adj->dat[i][j] *  (l1*l2 / pow(max(l1,l2), 2));
    // if(mm->adj->dat[i][j] > 5) mm->adj->dat[i][j] = 2;


    }
  }


  fprintf(stderr, "INFO: normalized adjacency matrix by intra contig link counts.\n");


  return 0;
}
