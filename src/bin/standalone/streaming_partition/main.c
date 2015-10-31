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

//int temp= mtmetis_partkway(
    // nv, 
    // xadj, 
    // adjncy, 
    // NULL, 
    // NULL,
    // numPart_, 
    // (mtmetis_pid_t*)verPartitionArray, 
    // NULL);


void partitionStinger(struct stinger* GSting, const  int64_t nv,const  int64_t numPart, int64_t* verPartitionArray){

  int64_t * off;
  int64_t * ind;
  int64_t * weight;
  int64_t nv_=nv;
  int64_t numPart_=numPart;
  stinger_to_unsorted_csr (GSting, nv+1, (int64_t**)&off, (int64_t**)&ind, (int64_t**)&weight,NULL, NULL, NULL);


  mtmetis_vtx_t *xweight= (mtmetis_vtx_t*)malloc((nv+1)*sizeof(mtmetis_vtx_t));
  for(int v=0;v<nv; v++)
    xweight[v]=1;

  mtmetis_vtx_t *xadj=(mtmetis_vtx_t*)off;
  mtmetis_adj_t *adjncy=(mtmetis_adj_t*)ind; 

  printf("\nnumber of vertices and edges: %ld and %ld \n",nv, xadj[nv+1]);

  double* options=mtmetis_init_options();

  options[MTMETIS_OPTION_NPARTS]=numPart;
//  options[MTMETIS_OPTION_CTYPE]=MTMETIS_CTYPE_SHEM;
  options[MTMETIS_OPTION_PTYPE]=MTMETIS_PTYPE_KWAY;
//  options[MTMETIS_OPTION_CONTYPE]=MTMETIS_CONTYPE_SORT;
//  options[MTMETIS_OPTION_RTYPE]=MTMETIS_RTYPE_GREEDY;
//  options[MTMETIS_OPTION_VERBOSITY]=MTMETIS_VERBOSITY_LOW;

//  options[MTMETIS_OPTION_NRUNS]=10;
//  options[MTMETIS_OPTION_NITER]=10;

//  options[MTMETIS_OPTION_IGNORE]=MTMETIS_IGNORE_VERTEXWEIGHTS+MTMETIS_IGNORE_EDGEWEIGHTS;
  
int temp= mtmetis_partition_explicit(
    nv, 
    xadj, 
    adjncy, 
    xweight, 
    weight,
    options, 
    (mtmetis_pid_t*)verPartitionArray, 
    NULL);
printf("Metis errocode : %d && %d",temp, temp==MTMETIS_SUCCESS);

  free(xweight);
  free(options);
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


    for(int i=0; i<naction;i++)
      if(action[2*i]>nv||action[2*i+1]>nv)
      printf("(%ld , %ld)  , ",action[2*i],action[2*i+1]);

  print_initial_graph_stats (nv, ne, batch_size, nbatch, naction);
  BATCH_SIZE_CHECK();

#if defined(_OPENMP)
  OMP("omp parallel") {
  OMP("omp master")
  PRINT_STAT_INT64 ("num_threads", (long int) omp_get_num_threads());
  }
#endif



  mtmetis_vtx_t *xadj=(mtmetis_vtx_t*)off;
  mtmetis_adj_t *adjncy=(mtmetis_adj_t*)ind; 

  double* options=mtmetis_init_options();

  options[MTMETIS_OPTION_NPARTS]=2;
  options[MTMETIS_OPTION_CTYPE]=MTMETIS_CTYPE_SHEM;
  options[MTMETIS_OPTION_PTYPE]=MTMETIS_PTYPE_ESEP;
  options[MTMETIS_OPTION_CONTYPE]=MTMETIS_CONTYPE_SORT;
  options[MTMETIS_OPTION_RTYPE]=MTMETIS_RTYPE_GREEDY;
  options[MTMETIS_OPTION_NRUNS]=10;
  options[MTMETIS_OPTION_NITER]=10;

  mtmetis_vtx_t *csrPart= (mtmetis_vtx_t*)malloc(nv*sizeof(mtmetis_vtx_t));

  options[MTMETIS_OPTION_IGNORE]=MTMETIS_IGNORE_VERTEXWEIGHTS+MTMETIS_IGNORE_EDGEWEIGHTS;
  mtmetis_partition_explicit(
    nv, 
    xadj, 
    adjncy, 
    NULL, 
    weight,
    options, 
    (mtmetis_pid_t*)csrPart, 
    NULL);


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


  int64_t inserted=0;
  for (int64_t actno = 0; actno < nbatch * batch_size; actno += batch_size)
  {
    tic();
    int64_t diffPartitions=0;
    const int64_t endact = (((actno + batch_size) > naction) ? naction : actno + batch_size);
    int64_t *actions = &action[2*actno];
    int64_t numActions = endact - actno;

    MTA("mta assert parallel")
    MTA("mta block dynamic schedule")
    OMP("omp parallel for")
    for(uint64_t k = 0; k < endact - actno; k++) {
      const int64_t i = actions[2 * k];
      const int64_t j = actions[2 * k + 1];
        if (i != j && i >= 0) {
          if(prevPart[i]!=prevPart[j]){
            __sync_add_and_fetch(&diffPartitions,1);
          }
      	int val=stinger_insert_edge (S, 0, i, j, 1, actno+2);
      	val+=stinger_insert_edge (S, 0, j, i, 1, actno+2);
        __sync_add_and_fetch(&inserted,val);
      }
    }
    update_time_trace[ntrace] = toc();
    PRINT_STAT_INT64("Cross partition edges - before :", diffPartitions);
    partitionStinger(S,nv, 2, currPart);
    if(ntrace==0)
      memcpy(firstPart,currPart,nv*sizeof(mtmetis_vtx_t));

    temp=prevPart;  prevPart=currPart; currPart=temp;
    ntrace++;
   diffPartitions=0;
    MTA("mta assert parallel")
    MTA("mta block dynamic schedule")
    OMP("omp parallel for")
    for(uint64_t k = 0; k < endact - actno; k++) {
      const int64_t i = actions[2 * k];
      const int64_t j = actions[2 * k + 1];
        if (i != j && i >= 0) {
          if(prevPart[i]!=prevPart[j]){
            __sync_add_and_fetch(&diffPartitions,1);
          }
      }
    }
    update_time_trace[ntrace] = toc();
    PRINT_STAT_INT64("Cross partition edges -  after :", diffPartitions);
  } /* End of batch */

  memcpy(lastPart,prevPart,nv*sizeof(mtmetis_vtx_t));

  int64_t diffPartitions=0,diffCsrPartitions=0;

  for (int64_t actno = 0; actno < nbatch * batch_size; actno += batch_size)
  {

    const int64_t endact = (actno + batch_size > naction ? naction : actno + batch_size);
    int64_t *actions = &action[2*actno];
    int64_t numActions = endact - actno;

    MTA("mta assert parallel")
    MTA("mta block dynamic schedule")
    OMP("omp parallel for")
    for(uint64_t k = 0; k < endact - actno; k++) {
      const int64_t i = actions[2 * k];
      const int64_t j = actions[2 * k + 1];
        if (i != j && i >= 0) {
          if(lastPart[i]!=lastPart[j]){
            __sync_add_and_fetch(&diffPartitions,1);
          }
          if(csrPart[i]!=csrPart[j]){
            __sync_add_and_fetch(&diffCsrPartitions,1);
          }
      }
    }

  } /* End of batch */

    PRINT_STAT_INT64("Actual number of inserted edges :", inserted);

    PRINT_STAT_INT64("Cross partition for last iteration :", diffPartitions);
    PRINT_STAT_INT64("Cross partition for CSR :", diffCsrPartitions);


  free(csrPart);

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
