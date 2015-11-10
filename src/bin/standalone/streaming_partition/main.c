/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#define _XOPEN_SOURCE 600
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <stdlib.h>

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

#define PARTITIONS_DEFAULT 4

static long batch_size = BATCH_SIZE_DEFAULT;
static long nbatch = NBATCH_DEFAULT;
static long npartitions = PARTITIONS_DEFAULT;

static struct stinger * S;

typedef struct {
  int64_t diffPartBefore;
  int64_t diffPartAfter;
  int64_t insertedEdges;
  int64_t edgesInGraph;
  int64_t graphSizeDoubled;
  int64_t* partAfterDoubling;  
  int64_t nvOverlapFirst;
  int64_t nvOverlapPrevious;
  int64_t nvOverlapDouble;
  int64_t nvOverlapLast;

} iterationInfo;


void prettyPrint(iterationInfo* extraInfo, int64_t nbatch, int64_t nv, int64_t ne, int64_t npartitions){
  printf("\n");
  nv--;
  for(int64_t b=0;b<nbatch; b++){
    printf("%11ld, ",nv);
    printf("%11ld, ",npartitions);
    printf("%11ld, ",extraInfo[b].diffPartBefore);
    printf("%11ld, ",extraInfo[b].diffPartAfter);
    printf("%11ld, ",extraInfo[b].insertedEdges);
    printf("%11ld, ",extraInfo[b].edgesInGraph);
    printf("%1.10lf, ",(double)extraInfo[b].edgesInGraph/(double)ne);

    printf("%1.10lf, ",(double)extraInfo[b].nvOverlapFirst/(double)nv);
    printf("%1.10lf, ",(double)extraInfo[b].nvOverlapPrevious/(double)nv);
    if (extraInfo[b].graphSizeDoubled){
      printf("%1.10lf, ",(double)extraInfo[b].nvOverlapDouble/(double)nv);
      printf("%1.10lf, ",(double)extraInfo[b].nvOverlapLast/(double)nv);
    }

/*    printf("%11ld, ",);
    printf("%11ld, ",);
    printf("%11ld, ",);
    printf("%11ld, ",);

    extraInfo[b].nvOverlapLast    =  compPartitions(S,nv, fullSize,lastPart);
    extraInfo[b].nvOverlapDouble  =  compPartitions(S,nv, halfSize,fullSize);


extraInfo[b].nvOverlapFirst
extraInfo[b].nvOverlapPrevious
extraInfo[b].nvOverlapLast
extraInfo[b].nvOverlapDouble
    if (extraInfo[b].graphSizeDoubled){
    }
*/    
    printf("\n");
  }


}


int64_t compPartitions(struct stinger* GSting, const  int64_t nv, 
  int64_t* verPartitionBeforeArray, int64_t* verPartitionAfterArray);

void parse_args_streaming_partitions (const int argc, char *argv[], char ** graph_name,
    char **action_name, int64_t* batch_num, int64_t* batch_size, int64_t* num_partitions);

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
  parse_args_streaming_partitions (argc,argv, &initial_graph_name,&action_stream_name,
    &nbatch,&batch_size, &npartitions);
//  parse_args (argc, argv, &initial_graph_name, &action_stream_name, &batch_size, &nbatch);
  STATS_INIT();
  load_graph_and_action_stream (initial_graph_name, &nv, &ne, (int64_t**)&off,
	      (int64_t**)&ind, (int64_t**)&weight, (int64_t**)&graphmem,
	      action_stream_name, &naction, (int64_t**)&action, (int64_t**)&actionmem);
  
  print_initial_graph_stats (nv, ne, batch_size, nbatch, naction);
  BATCH_SIZE_CHECK();
  PRINT_STAT_INT64("number of partitions", npartitions);

  nbatch=(naction/batch_size)+((naction%batch_size)>1);
  iterationInfo* extraInfo=(iterationInfo*)malloc(sizeof(iterationInfo)*nbatch);

#if defined(_OPENMP)
  OMP("omp parallel") {
  OMP("omp master")
  PRINT_STAT_INT64 ("num_threads", (long int) omp_get_num_threads());
  }
#endif

  int64_t initialNe=ne, currentNe=ne, lastDoubleNe=ne, doubledFlag=0, finalNe=ne+naction;

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

  int64_t inserted=0,totalInserted=0;
  int64_t batchId=0;

  partitionStinger(S,nv, npartitions, firstPart);

  // Going through all the batches
  for (int64_t actno = 0; actno < nbatch * batch_size; actno += batch_size,batchId++){
    int64_t diffPartitionsBefore=0,diffPartitionsAfter=0;
    const int64_t endact = (((actno + batch_size) > naction) ? naction : actno + batch_size);
    int64_t *actions = &action[2*actno];
    int64_t numActions = endact - actno;
    inserted=0;
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
      	int val =stinger_insert_edge (S, 0, i, j, 1, actno+2);
      	val+=stinger_insert_edge (S, 0, j, i, 1, actno+2);
        __sync_add_and_fetch(&inserted,val); 
      }
    }

    partitionStinger(S,nv, npartitions, currPart);
    // Swapping partition arrays every iteration;
    int64_t* temp;temp=prevPart;  prevPart=currPart; currPart=temp;
    
    OMP("omp parallel for")
    for(uint64_t k = 0; k < endact - actno; k++) {
      const int64_t i = actions[2 * k]; const int64_t j = actions[2 * k + 1];
        if (i != j && i >= 0 && prevPart[i]!=prevPart[j]) {
          __sync_add_and_fetch(&diffPartitionsAfter,1);
      }
    }

    currentNe+=inserted;
    extraInfo[batchId].graphSizeDoubled=0;
    if (currentNe > (2*lastDoubleNe) || batchId==(nbatch-1)){
      lastDoubleNe=currentNe;
      extraInfo[batchId].graphSizeDoubled=1;
      extraInfo[batchId].partAfterDoubling=(int64_t*)malloc(nv*sizeof(int64_t));
      memcpy(extraInfo[batchId].partAfterDoubling,prevPart,nv*sizeof(int64_t));
    }

    extraInfo[batchId].diffPartBefore=diffPartitionsBefore;
    extraInfo[batchId].diffPartAfter=diffPartitionsAfter;
    extraInfo[batchId].insertedEdges=inserted;
    extraInfo[batchId].edgesInGraph = currentNe;
    extraInfo[batchId].nvOverlapPrevious    =  compPartitions(S,nv, prevPart,currPart);
    extraInfo[batchId].nvOverlapFirst =  compPartitions(S,nv, prevPart,firstPart);
    totalInserted+=inserted;
  } /* End of batch */

  memcpy(lastPart,prevPart,nv*sizeof(int64_t));

  int64_t diffPartitions=0,diffCsrPartitions=0;
  for (int64_t actno = 0; actno < nbatch * batch_size; actno += batch_size)  {
    const int64_t endact = (actno + batch_size > naction ? naction : actno + batch_size);
    int64_t *actions = &action[2*actno];
    int64_t numActions = endact - actno;

    OMP("omp parallel for")
    for(uint64_t k = 0; k < endact - actno; k++) {
      const int64_t i = actions[2 * k]; const int64_t j = actions[2 * k + 1];
        if (i != j && i >= 0 && lastPart[i]!=lastPart[j]) {
          __sync_add_and_fetch(&diffPartitions,1);
      }
    }
  } 

  int64_t* halfSize=firstPart, *fullSize;
  for(int64_t b=0;b<nbatch; b++){
    if (extraInfo[b].graphSizeDoubled){
      fullSize=extraInfo[b].partAfterDoubling;
      extraInfo[b].nvOverlapDouble  =  compPartitions(S,nv, halfSize,fullSize);
      extraInfo[b].nvOverlapLast    =  compPartitions(S,nv, fullSize,lastPart);
      halfSize=fullSize;
    }
  }
  prettyPrint(extraInfo,nbatch,nv,currentNe,npartitions);

  for(int64_t b=0;b<nbatch; b++){
    if(extraInfo[b].graphSizeDoubled)
      free(extraInfo[b].partAfterDoubling);
  }

  free(firstPart);  free(prevPart);  free(currPart);  free(lastPart); 
  free(extraInfo);

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



void usage (FILE * out, char *progname)
{
  fprintf (out,
           "%s [--batch-size|-b #] [--nbatch|-n #] [--partitions|-p #] [graph-name.bin] [action-stream.bin]\n"
           "\tDefaults:\n"
           "\t   Batch size, default vale = %d\n"
           "\t   Number of batches, default vale = %d\n"
           "\t   Number of partitions, default vale = %d\n"
           "\t   stinger graph name, default vale = \"%s\"\n"
           "\t   stinger actions stream name, default vale = \"%s\"\n",
           progname,
           BATCH_SIZE_DEFAULT,
           NBATCH_DEFAULT,
           PARTITIONS_DEFAULT, 
           INITIAL_GRAPH_NAME_DEFAULT,
           ACTION_STREAM_NAME_DEFAULT);
}

/**
* @brief Parses command line arguments.
*
* Parses the command line input as given by usage().  Batch size, number of
* batches, initial graph filename, and action stream filename are given by
* to the caller if they were specified on the command line.
*
* @param argc The number of arguments
* @param argv[] The array of arguments
* @param initial_graph_name Path/filename of the initial input graph on disk
* @param action_stream_name Path/filename of the action stream on disk
* @param batch_size Number of edge actions to consider in one batch
* @param nbatch Number of batchs to process
* @param num_partition Number of batchs to process
*/
void parse_args_streaming_partitions (const int argc, char *argv[], char ** graph_name,
    char **action_name, int64_t* batch_num, int64_t* batch_size, int64_t* num_partitions) {
  int k = 1;
  int seenPartition=0, seenBatchSize=0, seenNumBatches=0;
  if (k >= argc)
    return;
  while (k < argc && argv[k][0] == '-') {
    if (0 == strcmp (argv[k], "--percentage") || 0 == strcmp (argv[k], "-p")) {
      if (seenPartition)
        goto err;
      seenPartition = 1;
      ++k;
      if (k >= argc)
        goto err;
      *num_partitions = strtol (argv[k], NULL,10) ;
      if (*num_partitions<=1)
        goto err;
      ++k;
    } 
    else if (0 == strcmp (argv[k], "--batch-size") || 0 == strcmp (argv[k], "-b")) {
      if (seenBatchSize)
        goto err;
      seenBatchSize = 1;
      ++k;
      if (k >= argc)
        goto err;
      *batch_size = strtol (argv[k], NULL,10);
      if (*batch_size<1)
        goto err;
      ++k;
    } 
    else if (0 == strcmp (argv[k], "--nbatch") || 0 == strcmp (argv[k], "-n")) {
      if (seenNumBatches)
        goto err;
      seenNumBatches = 1;
      ++k;
      if (k >= argc)
        goto err;
      *batch_num = strtol (argv[k], NULL,10);
      if (*batch_num<0)
        goto err;
      ++k;
    }     
    else if (0 == strcmp (argv[k], "--help")
             || 0 == strcmp (argv[k], "-h") || 0 == strcmp (argv[k], "-?")) {
      usage (stdout, argv[0]);
      exit (EXIT_SUCCESS);
      return;
    } else if (0 == strcmp (argv[k], "--")) {
      ++k;
      break;
    }
  }
  if (k < argc)
    *graph_name = argv[k++];
  if (k < argc)
    *action_name = argv[k++];
  return;
 err:
  usage (stderr, argv[0]);
  exit (EXIT_FAILURE);
  return;
}



