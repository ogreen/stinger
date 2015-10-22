/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#define _XOPEN_SOURCE 600
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#if defined(_OPENMP)
#include "omp.h"
#endif

#define MTMETIS_64BIT_VERTICES
#define MTMETIS_64BIT_EDGES
#define MTMETIS_64BIT_WEIGHTS

#include <mtmetis.h>


#include "stinger_core/stinger_atomics.h"
#include "stinger_utils/stinger_utils.h"
#include "stinger_core/stinger.h"
#include "stinger_utils/timer.h"
#include "stinger_core/xmalloc.h"



#define ACTI(k) (action[2*(k)])
#define ACTJ(k) (action[2*(k)+1])

static int64_t nv, ne, naction;
static int64_t * restrict off;
static int64_t * restrict from;
static int64_t * restrict ind;
static int64_t * restrict weight;
static int64_t * restrict action;

/* handles for I/O memory */
static int64_t * restrict graphmem;
static int64_t * restrict actionmem;

static char * initial_graph_name = INITIAL_GRAPH_NAME_DEFAULT;
static char * action_stream_name = ACTION_STREAM_NAME_DEFAULT;

static long batch_size = BATCH_SIZE_DEFAULT;
static long nbatch = NBATCH_DEFAULT;

static struct stinger * S;

static double * update_time_trace;


void partitionStinger(struct stinger* GSting, const  int64_t nv,const  int64_t numPart, int64_t* verPartitionArray){

  int64_t * off;
  int64_t * ind;
  int64_t * weight;
  int64_t nv_=nv;
  int64_t numPart_=numPart;
  stinger_to_unsorted_csr (GSting, nv, (int64_t**)&off, (int64_t**)&ind, (int64_t**)&weight,NULL, NULL, NULL);

  mtmetis_vtx_t *xadj=(mtmetis_vtx_t*)off;
  mtmetis_adj_t *adjncy=(mtmetis_adj_t*)ind; 
  mtmetis_vtx_t ncon=1;
  mtmetis_vtx_t objval;
/*
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);

  options[METIS_OPTION_CTYPE]=METIS_CTYPE_SHEM;

  options[METIS_OPTION_NUMBERING]=0;
  options[METIS_OPTION_NO2HOP]=1;
*/
//  int temp= METIS_PartGraphRecursive(&nv_, &ncon, xadj, adjncy,  NULL, NULL, NULL, 
//    &numPart_, NULL, NULL, options, &objval, verPartitionArray);
//  int temp= METIS_PartGraphKway(&nv_, &ncon, xadj, adjncy,  NULL, NULL, NULL, 
//    &numPart_, NULL, NULL, options, &objval, verPartitionArray);

int temp= mtmetis_partkway(
    nv, 
    xadj, 
    adjncy, 
    NULL, 
    NULL,
    numPart_, 
    (mtmetis_pid_t*)verPartitionArray, 
    NULL);


  free(off);
  free(ind);
  free(weight);
}

int
main (const int argc, char *argv[])
{
  parse_args (argc, argv, &initial_graph_name, &action_stream_name, &batch_size, &nbatch);
  STATS_INIT();

  load_graph_and_action_stream (initial_graph_name, &nv, &ne, (int64_t**)&off,
	      (int64_t**)&ind, (int64_t**)&weight, (int64_t**)&graphmem,
	      action_stream_name, &naction, (int64_t**)&action, (int64_t**)&actionmem);

  print_initial_graph_stats (nv, ne, batch_size, nbatch, naction);
  BATCH_SIZE_CHECK();

#if defined(_OPENMP)
  OMP("omp parallel") {
  OMP("omp master")
  PRINT_STAT_INT64 ("num_threads", (long int) omp_get_num_threads());
  }
#endif


  update_time_trace = xmalloc (nbatch * sizeof(*update_time_trace));

  /* Convert to STINGER */
  tic ();
  S = stinger_new ();
  stinger_set_initial_edges (S, nv, 0, off, ind, weight, NULL, NULL, -2);
  PRINT_STAT_DOUBLE ("time_stinger", toc ());
  fflush(stdout);
  free(graphmem);


  tic ();
  uint32_t errorCode = stinger_consistency_check (S, nv);
  double time_check = toc ();
  PRINT_STAT_HEX64 ("error_code", (long unsigned) errorCode);
  PRINT_STAT_DOUBLE ("time_check", time_check);

  /* Updates */
  int64_t ntrace = 0;

  mtmetis_vtx_t *firstPart= (mtmetis_vtx_t*)malloc(nv*sizeof(mtmetis_vtx_t));

  mtmetis_vtx_t *currPart= (mtmetis_vtx_t*)malloc(nv*sizeof(mtmetis_vtx_t));
  mtmetis_vtx_t *prevPart= (mtmetis_vtx_t*)malloc(nv*sizeof(mtmetis_vtx_t));
  mtmetis_vtx_t* temp;
  mtmetis_vtx_t *lastPart= (mtmetis_vtx_t*)malloc(nv*sizeof(mtmetis_vtx_t));


  for (int64_t actno = 0; actno < nbatch * batch_size; actno += batch_size)
  {
    tic();
    int64_t diffPartitions=0;
    const int64_t endact = (actno + batch_size > naction ? naction : actno + batch_size);
    int64_t *actions = &action[2*actno];
    int64_t numActions = endact - actno;

    MTA("mta assert parallel")
    MTA("mta block dynamic schedule")
    OMP("omp parallel for")
    for(uint64_t k = 0; k < endact - actno; k++) {
      const int64_t i = actions[2 * k];
      const int64_t j = actions[2 * k + 1];
/* // DISABLING Deletions for this experiment 

      if (i != j && i < 0) {
      	stinger_remove_edge(S, 0, ~i, ~j);
      	stinger_remove_edge(S, 0, ~j, ~i);
      }
*/      
        if (i != j && i >= 0) {
          if(prevPart[i]!=prevPart[j]){
            __sync_add_and_fetch(&diffPartitions,1);
          }
      	stinger_insert_edge (S, 0, i, j, 1, actno+2);
      	stinger_insert_edge (S, 0, j, i, 1, actno+2);
      }
    }

    update_time_trace[ntrace] = toc();
    PRINT_STAT_INT64("Cross partition edges", diffPartitions);
    partitionStinger(S,nv, 4, currPart);
    if(ntrace==0)
      memcpy(firstPart,currPart,nv*sizeof(mtmetis_vtx_t));
    

    temp=prevPart;  prevPart=currPart; currPart=temp;
    ntrace++;
  } /* End of batch */

  memcpy(lastPart,prevPart,nv*sizeof(mtmetis_vtx_t));


  free(firstPart);  free(prevPart);  free(currPart);  free(lastPart);

  /* Print the times */
  double time_updates = 0;
  for (int64_t k = 0; k < nbatch; k++) {
    time_updates += update_time_trace[k];
  }
  PRINT_STAT_DOUBLE ("time_updates", time_updates);
  PRINT_STAT_DOUBLE ("updates_per_sec", (nbatch * batch_size) / time_updates); 

  tic ();
  errorCode = stinger_consistency_check (S, nv);
  time_check = toc ();
  PRINT_STAT_HEX64 ("error_code", (long unsigned) errorCode);
  PRINT_STAT_DOUBLE ("time_check", time_check);

  free(update_time_trace);
  stinger_free_all (S);
  free (actionmem);
  STATS_END();



}
