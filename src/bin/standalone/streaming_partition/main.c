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

#define METIS 1

#if METIS==1
#include <metis.h>
#else
#include <mtmetis.h>
#endif

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

int64_t compPartitions(struct stinger* GSting, const  int64_t nv, 
  int64_t* verPartitionBeforeArray, int64_t* verPartitionAfterArray);

void partitionStinger(struct stinger* GSting, const  int64_t nv,const  int64_t numPart, int64_t* verPartitionArray){
  int64_t * off;
  int64_t * ind;
  int64_t * weight;
  int64_t nv_=nv;
  int64_t numPart_=numPart;
  stinger_to_unsorted_csr (GSting, nv+1, (int64_t**)&off, (int64_t**)&ind, (int64_t**)&weight,NULL, NULL, NULL);



#if METIS==1
    idx_t options[METIS_NOPTIONS];
    idx_t ncon=1;
    METIS_SetDefaultOptions(options);

    options[METIS_OPTION_CTYPE]=METIS_CTYPE_SHEM;
    options[METIS_OPTION_PTYPE]=METIS_PTYPE_KWAY;
    options[METIS_OPTION_IPTYPE]=METIS_IPTYPE_EDGE;
    options[METIS_OPTION_OBJTYPE]=METIS_OBJTYPE_CUT;
    options[METIS_OPTION_NUMBERING]=0;
    options[METIS_OPTION_NO2HOP]=1;
    idx_t objval;//= (idx_t*)malloc((numPart_)*sizeof(idx_t));
    idx_t* xadj=(idx_t*)off;
    idx_t* adjncy=(idx_t*)ind;

  //  int temp= METIS_PartGraphRecursive(&nv_, &ncon, xadj, adjncy,  NULL, NULL, NULL, 
  //    &numPart_, NULL, NULL, options, &objval, verPartitionArray);
    int temp= METIS_PartGraphKway(&nv_, &ncon, xadj, adjncy,  NULL, NULL, NULL, 
      &numPart_, NULL, NULL, options, &objval, (idx_t*)verPartitionArray);

//    printf("edgecut size: %ld", objval);
#else
/*
    mtmetis_vtx_t *xadj=(mtmetis_vtx_t*)off;
    mtmetis_adj_t *adjncy=(mtmetis_adj_t*)ind; 


    double* options=mtmetis_init_options();
    mtmetis_vtx_t *xweight= (mtmetis_vtx_t*)malloc((nv+1)*sizeof(mtmetis_vtx_t));
    for(int v=0;v<nv; v++)
      xweight[v]=1;

    printf("\nnumber of vertices and edges: %ld and %ld \n",nv, xadj[nv+1]);


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
*/
#endif

  free(off);
  free(ind);
  free(weight);
}



int main (const int argc, char *argv[]){
  parse_args (argc, argv, &initial_graph_name, &action_stream_name, &batch_size, &nbatch);
  STATS_INIT();
  load_graph_and_action_stream (initial_graph_name, &nv, &ne, (int64_t**)&off,
	      (int64_t**)&ind, (int64_t**)&weight, (int64_t**)&graphmem,
	      action_stream_name, &naction, (int64_t**)&action, (int64_t**)&actionmem);

  print_initial_graph_stats (nv, ne, batch_size, nbatch, naction);
  BATCH_SIZE_CHECK();

  nbatch=naction/batch_size;

#if defined(_OPENMP)
  OMP("omp parallel") {
  OMP("omp master")
  PRINT_STAT_INT64 ("num_threads", (long int) omp_get_num_threads());
  }
#endif

  int64_t initialNe=ne, currentNe=ne, lastDoubleNe=ne, doubledFlag=0;

  /* Convert to STINGER */
  tic (); S = stinger_new ();
  stinger_set_initial_edges (S, nv, 0, off, ind, weight, NULL, NULL, -2);
  PRINT_STAT_DOUBLE ("time_stinger", toc ());
  fflush(stdout);
  free(graphmem);


  tic ();
  uint32_t errorCode = stinger_consistency_check (S, nv);
  double time_check = toc ();
  PRINT_STAT_HEX64 ("error_code", (long unsigned) errorCode);
  PRINT_STAT_DOUBLE ("time_check", time_check);


  int64_t *firstPart  = (int64_t*)malloc(nv*sizeof(int64_t));
  int64_t *currPart   = (int64_t*)malloc(nv*sizeof(int64_t));
  int64_t *prevPart   = (int64_t*)malloc(nv*sizeof(int64_t));
  int64_t *lastPart   = (int64_t*)malloc(nv*sizeof(int64_t));
  int64_t *doublePart = (int64_t*)malloc(nv*sizeof(int64_t));

  int64_t inserted=0,val;
  int64_t batchId=0;

  partitionStinger(S,nv, 2, firstPart);

  // Going through all the batches
  for (int64_t actno = 0; actno < nbatch * batch_size; actno += batch_size,batchId++){
    int64_t diffPartitionsBefore=0,diffPartitionsAfter=0;
    const int64_t endact = (((actno + batch_size) > naction) ? naction : actno + batch_size);
    int64_t *actions = &action[2*actno];
    int64_t numActions = endact - actno;

    // Adding the edges into the network 
    OMP("omp parallel for")
    for(uint64_t k = 0; k < endact - actno; k++) {
      const int64_t i = actions[2 * k];
      const int64_t j = actions[2 * k + 1];
        if (i != j && i >= 0) {
          // Looking for the edges that are cross partition
          if(prevPart[i]!=prevPart[j]){
            __sync_add_and_fetch(&diffPartitionsBefore,1);
          }
        // Counting the number of edges inserted into the network for sanity checking
      	val =stinger_insert_edge (S, 0, i, j, 1, actno+2);
      	val+=stinger_insert_edge (S, 0, j, i, 1, actno+2);
        __sync_add_and_fetch(&inserted,val); 
      }
    }

    partitionStinger(S,nv, 2, currPart);
    currentNe+=inserted;
    int64_t samePartitions=compPartitions(S,nv, prevPart,currPart);
    // Swapping partition arrays every iteration;
    int64_t* temp;temp=prevPart;  prevPart=currPart; currPart=temp;
    
    // Graph has doubled the number of edges
    if (currentNe > (2*lastDoubleNe)){
      lastDoubleNe=currentNe;
      doubledFlag=1;
      memcpy(doublePart,prevPart,nv*sizeof(int64_t));      
    }

    OMP("omp parallel for")
    for(uint64_t k = 0; k < endact - actno; k++) {
      const int64_t i = actions[2 * k];
      const int64_t j = actions[2 * k + 1];
        if (i != j && i >= 0) {
          if(prevPart[i]!=prevPart[j]){
            __sync_add_and_fetch(&diffPartitionsAfter,1);
          }
      }
    }
    PRINT_STAT_INT64("Cross partition edges - before :", diffPartitionsBefore);    
    PRINT_STAT_INT64("Cross partition edges -  after :", diffPartitionsAfter);
    PRINT_STAT_INT64("Same partition before and after", samePartitions);    
  } /* End of batch */

  memcpy(lastPart,prevPart,nv*sizeof(int64_t));

  int64_t diffPartitions=0,diffCsrPartitions=0;

  for (int64_t actno = 0; actno < nbatch * batch_size; actno += batch_size)  {
    const int64_t endact = (actno + batch_size > naction ? naction : actno + batch_size);
    int64_t *actions = &action[2*actno];
    int64_t numActions = endact - actno;

    OMP("omp parallel for")
    for(uint64_t k = 0; k < endact - actno; k++) {
      const int64_t i = actions[2 * k];
      const int64_t j = actions[2 * k + 1];
      if (i != j && i >= 0) {
        if(lastPart[i]!=lastPart[j]){
          __sync_add_and_fetch(&diffPartitions,1);
        }
      }
    }
    
  } 
  /* End of batch */
  PRINT_STAT_INT64("Actual number of inserted edges :", inserted);
  PRINT_STAT_INT64("Cross partition for last iteration :", diffPartitions);

  free(firstPart);  free(prevPart);  free(currPart);  free(lastPart); free(doublePart);

  errorCode = stinger_consistency_check (S, nv);
  PRINT_STAT_HEX64 ("error_code", (long unsigned) errorCode);

  stinger_free_all (S);
  free (actionmem);
  STATS_END();

}


// Given two differnt partitions of the graph this function checks if the vertices stay within 
// the same partition. 
int64_t compPartitions(struct stinger* GSting, const  int64_t nv, int64_t* verPartitionBeforeArray,
    int64_t* verPartitionAfterArray){

    int64_t *countSameAfter= (int64_t*)malloc((nv+1)*sizeof(int64_t));
    int64_t *countSameBefore= (int64_t*)malloc((nv+1)*sizeof(int64_t));
    int64_t samePartitions=0;

    OMP("omp parallel for")
    for(int64_t src=0; src<nv; src++){
      countSameAfter[src]=0;
      countSameBefore[src]=0;
      int edgeCounter=0;
      STINGER_FORALL_EDGES_OF_VTX_BEGIN(GSting,src)
        int64_t dest=STINGER_EDGE_DEST;

        if(verPartitionBeforeArray[dest]==verPartitionBeforeArray[src])
          countSameBefore[src]++;
        if(verPartitionAfterArray[dest]==verPartitionAfterArray[src])
          countSameAfter[src]++;
        edgeCounter++;
      STINGER_FORALL_EDGES_OF_VTX_END();
      double percentage=(double)(countSameAfter[src]-countSameBefore[src])/(double)edgeCounter;
      if(fabs(percentage)<0.3)
        __sync_add_and_fetch(&samePartitions,1);
    }
    free(countSameAfter);
    free(countSameBefore);
    return samePartitions;
}

