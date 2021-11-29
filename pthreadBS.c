/*
 bitonic.c

 This file contains two different implementations of the bitonic sort
        recursive  version :  recBitonicSort()
        imperative version :  impBitonicSort()


 The bitonic sort is also known as Batcher Sort.
 For a reference of the algorithm, see the article titled
 Sorting networks and their applications by K. E. Batcher in 1968


 The following codes take references to the codes avaiable at

 http://www.cag.lcs.mit.edu/streamit/results/bitonic/code/c/bitonic.c

 http://www.tools-of-computing.com/tc/CS/Sorts/bitonic_sort.htm

 http://www.iti.fh-flensburg.de/lang/algorithmen/sortieren/bitonic/bitonicen.htm
 */

/*
------- ----------------------
   Nikos Pitsianis, Duke CS
-----------------------------
*/


#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>

struct timeval startwtime, endwtime;
double seq_time;

typedef struct {
	int lo;		// starting point
	int cnt;	// length
	int dir;	// direction (asc or desc)
	int depth;	// depth of tree
}data;

int N;          // data array size
int *a;         // data array to be sorted
int d; 			// depth

const int ASCENDING  = 1;
const int DESCENDING = 0;


void init(void);
void print(void);
void sort(void);
void test(void);
void exchange(int i, int j);
void compare(int i, int j, int dir);
void bitonicMerge(int lo, int cnt, int dir);
void recBitonicSort(int lo, int cnt, int dir);
void impBitonicSort(void);

void* pthread_recBitonicSort(void* myData);
void* pthread_bitonicMerge(void* myData);
int Asc(const void * x, const void * y);
int Desc(const void * x, const void * y);


/** the main program **/
int main(int argc, char **argv) {

  if (argc != 3) {
    printf("Usage: %s q\n  where n=2^q is problem size (power of two)\n",
	   argv[0]);
    exit(1);
  }

  d = atoi(argv[2]);           // 2^d threads  -  d: max depth
  N = 1<<atoi(argv[1]);
  a = (int *) malloc(N * sizeof(int));

  init();

  gettimeofday (&startwtime, NULL);
  impBitonicSort();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Imperative wall clock time = %f\n", seq_time);

  test();

  init();
  gettimeofday (&startwtime, NULL);
  sort();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Recursive wall clock time = %f\n", seq_time);

  test();

  // print();

  // My parallel code

  data data1;


  data1.lo = 0;
  data1.cnt = N;
  data1.dir = ASCENDING;
  data1.depth = 0;
  init();

  gettimeofday (&startwtime, NULL);
  pthread_recBitonicSort(& data1);
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Recursive (pthread) wall clock time = %f\n", seq_time);

  test();

  // qsort

   init();

  gettimeofday (&startwtime, NULL);
  qsort(&a[0], N, sizeof(int), Asc);
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("qsort wall clock time = %f\n", seq_time);

  test();

}

/** -------------- SUB-PROCEDURES  ----------------- **/

/** procedure test() : verify sort results **/
void test() {
  int pass = 1;
  int i;
  for (i = 1; i < N; i++) {
    pass &= (a[i-1] <= a[i]);
  }

  printf(" TEST %s\n",(pass) ? "PASSed" : "FAILed");
}


/** procedure init() : initialize array "a" with data **/
void init() {
  int i;
  for (i = 0; i < N; i++) {
    a[i] = rand() % N; // (N - i);
  }
}

/** procedure  print() : print array elements **/
void print() {
  int i;
  for (i = 0; i < N; i++) {
    printf("%d\n", a[i]);
  }
  printf("\n");
}


/** INLINE procedure exchange() : pair swap **/
void exchange(int i, int j) {
  int t;
  t = a[i];
  a[i] = a[j];
  a[j] = t;
}



/** procedure compare()
   The parameter dir indicates the sorting direction, ASCENDING
   or DESCENDING; if (a[i] > a[j]) agrees with the direction,
   then a[i] and a[j] are interchanged.
**/
void compare(int i, int j, int dir) {
  if (dir==(a[i]>a[j]))
    exchange(i,j);
}




/** Procedure bitonicMerge()
   It recursively sorts a bitonic sequence in ascending order,
   if dir = ASCENDING, and in descending order otherwise.
   The sequence to be sorted starts at index position lo,
   the parameter cbt is the number of elements to be sorted.
 **/
void bitonicMerge(int lo, int cnt, int dir) {
  if (cnt>1) {
    int k=cnt/2;
    int i;
    for (i=lo; i<lo+k; i++)
      compare(i, i+k, dir);
    bitonicMerge(lo, k, dir);
    bitonicMerge(lo+k, k, dir);
  }
}



/** function recBitonicSort()
    first produces a bitonic sequence by recursively sorting
    its two halves in opposite sorting orders, and then
    calls bitonicMerge to make them in the same order
 **/
void recBitonicSort(int lo, int cnt, int dir) {
  if (cnt>1) {
    int k=cnt/2;
    recBitonicSort(lo, k, ASCENDING);
    recBitonicSort(lo+k, k, DESCENDING);
    bitonicMerge(lo, cnt, dir);
  }
}


/** function sort()
   Caller of recBitonicSort for sorting the entire array of length N
   in ASCENDING order
 **/
void sort() {
  recBitonicSort(0, N, ASCENDING);
}



/*
  imperative version of bitonic sort
*/
void impBitonicSort() {

  int i,j,k;

  for (k=2; k<=N; k=2*k) {
    for (j=k>>1; j>0; j=j>>1) {
      for (i=0; i<N; i++) {
	int ij=i^j;
	if ((ij)>i) {
	  if ((i&k)==0 && a[i] > a[ij])
	      exchange(i,ij);
	  if ((i&k)!=0 && a[i] < a[ij])
	      exchange(i,ij);
	}
      }
    }
  }
}


//=================================================================================================================================



void* pthread_recBitonicSort(void* myData) {
  data *dt;


  dt = (data*)myData;

  int l, c, dr, dpt;
  l = dt->lo;
  c = dt->cnt;
  dr = dt->dir;
  dpt = dt->depth;

  if (c>1) {

    if (dpt < d){

    	pthread_t t1, t2;
    	data dt1, dt2;



    	dt1.lo = l;
    	dt1.cnt = c/2;
    	dt1.dir = ASCENDING;
    	dt1.depth = dpt+1;

    	dt2.lo = l+c/2;
    	dt2.cnt = c/2;
    	dt2.dir = DESCENDING;
    	dt2.depth = dpt+1;


    	pthread_create(&t1, NULL, pthread_recBitonicSort,  &dt1);
    	pthread_create(&t2, NULL, pthread_recBitonicSort,  &dt2);


		pthread_join(t1, NULL);
		pthread_join(t2, NULL);

		pthread_bitonicMerge((void*)myData);

    }else{
    	qsort(&a[l], c/2, sizeof(int), Asc);
    	qsort(&a[l+c/2], c/2, sizeof(int), Desc);

    	bitonicMerge(l, c, dr);

    	/*
        if (dr == ASCENDING){                       // makes my function faster
            qsort(&a[l], c, sizeof(int), Asc);
        }else{
            qsort(&a[l], c, sizeof(int), Desc);
        }
    	*/
    }
  }
  return 0;
}



int Asc(const void * x, const void * y){
	return (*(int*)x - *(int*)y);
}


int Desc(const void * x, const void * y){
	return -(*(int*)x - *(int*)y);
}



void* pthread_bitonicMerge(void *myData){

	data *dt;


  	dt = (data*)myData;


  	int l, c, dr, dpt;
  	l = dt->lo;
  	c = dt->cnt;
  	dr = dt->dir;
  	dpt = dt->depth;

  	if (c > 1){
  		int i, k;
  		k=c/2;
  		for (i=l; i<l+k; i++){
      		compare(i, i+c/2, dr);
      	}

      	if (dpt < d){

		  	pthread_t t1, t2;
			data dt1, dt2;


			dt1.lo = l;
			dt1.cnt = c/2;
			dt1.dir = dr;
			dt1.depth = dpt+1;

			dt2.lo = l+c/2;
			dt2.cnt = c/2;
			dt2.dir = dr;
			dt2.depth = dpt+1;

			pthread_create(&t1, NULL, pthread_bitonicMerge, &dt1);
			pthread_create(&t2, NULL, pthread_bitonicMerge,  &dt2);


			pthread_join(t1, NULL);
			pthread_join(t2, NULL);
		}else{

			bitonicMerge(l, c/2, dr);
			bitonicMerge(l+c/2, c/2, dr);

		}
  	}
 return 0;
}









