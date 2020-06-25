#include <cstdio>
#include <cstdlib>
#include <vector>

int main() {
  int n = 50;
  int range = 5;
  std::vector<int> key(n);

  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]); 

}
  printf("\n");

<<<<<<< HEAD
  std::vector<int> bucket(range);
#pragma omp parallel for
  for (int i=0; i<range; i++) {
    bucket[i] = 0;
  }
#pragma omp parallel for shared(bucket)
  for (int i=0; i<n; i++) {
#pragma omp atomic update
    bucket[key[i]]++;
 }


#pragma omp parallel for 
  for (int i=0; i<range; i++) {
    int j = 0;
#pragma omp parallel for reduction(+:j)
  for(int k=0; k<i; k++){
     j += bucket[k];
    }

#pragma omp parallel for
  for (int a=bucket[i]; a>0; a--) {    
=======
  std::vector<int> bucket(range,0); 
#pragma omp parallel for
  for (int i=0; i<n; i++)
#pragma omp atomic update
    bucket[key[i]]++;
  std::vector<int> offset(range,0);
  for (int i=1; i<range; i++) 
    offset[i] = offset[i-1] + bucket[i-1];
#pragma omp parallel for
  for (int i=0; i<range; i++) {
    int j = offset[i];
    for (; bucket[i]>0; bucket[i]--) {
>>>>>>> upstream/master
      key[j++] = i;
    }
  }


  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
