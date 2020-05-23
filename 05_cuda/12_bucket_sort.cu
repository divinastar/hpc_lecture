#include <cstdio>
#include <cstdlib>
#include <vector>


__global__ void bucketsort(int *key, int *bucket, int n){
  
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i>=n) return;

  int a = key[i];
  atomicAdd(&bucket[a],1);
  __syncthreads();
  
  for(int j = 0, k = 0; j <= i; k++){
     key[i] = k;
     j += bucket[k];
  }
}


int main() {
  int n = 50;
  int range = 5;
  int *key;
  cudaMallocManaged(&key, n*sizeof(int));

  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");
  
  const int X = 32;
  int *bucket;

  cudaMallocManaged(&bucket, range*sizeof(int));  
  
  bucketsort<<<(n+X-1)/X,X,range>>>(key,bucket,n);
  cudaDeviceSynchronize();
  
  for (int i=0; i<n; i++) {
    printf("%d ", key[i]);
  }
  printf("\n");

//  for (int i=0; i<range; i++){
//    printf("%d ", bucket[i]);
//  }
//  printf("\n");
  
  cudaFree(key);
  cudaFree(bucket);
}


