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
      key[j++] = i;
    }
  }


  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
