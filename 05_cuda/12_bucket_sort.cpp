#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void bucketsort(std::vector<int> key, int *range){
  std::vector<int> bucket(range);
  for(int i=0); i<range; i++){
    bucket[i] = 0;
  }
  
  int j = blockIdx.x * blockDim.x + threadIdx.x;
}


int main() {
  int n = 50;
  int range = 5;
  std::vector<int> key(n);
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");
  
  const int X = 32;

//  std::vector<int> bucket(range); 
//  for (int i=0; i<range; i++) {
//    bucket[i] = 0;
//  }
  
//  for (int i=0; i<n; i++) {
//    bucket[key[i]]++;
//  }
 
//  for (int i=0, j=0; i<range; i++) {
//    for (; bucket[i]>0; bucket[i]--) {
//     key[j++] = i;
//    }
//  }

  bucketsort<<<((n+X-1)/X,X)>>>(key,range);
  cudaDeviceSynchronize();
  
  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
  
  cudafree(key);
}
