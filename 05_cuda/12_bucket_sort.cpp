#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void update_bucket(int *bucket,int *key) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  atomicAdd(&bucket[key[i]], 1);
}

__global__ void update_key(int *position,int *key){
   int i = blockIdx.x * blockDim.x + threadIdx.x; 
   for(int j=position[i];j<position[i+1];j++){
       key[j] = i;
   }    
}

__global__ void init(int *bucket){
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  bucket[i] = 0;
}

int main() {
  int n = 50;
  int range = 5;
  int *key,*bucket,*position;
  cudaMallocManaged(&key, n*sizeof(int));
  cudaMallocManaged(&bucket, range*sizeof(int)); 
  cudaMallocManaged(&position,(range+1)*sizeof(int));
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");
  
  init<<<1,n>>>(bucket);
  cudaDeviceSynchronize();  
  update_bucket<<<1,n>>>(bucket,key);  
  cudaDeviceSynchronize();  

  // get the information of position of bucket's element
  int cur = 0;
  for (int i=1; i<range; i++) {
    cur += bucket[i-1];  
    position[i] = cur;      
  }  
  position[range] = n;  

  
  update_key<<<1,range>>>(position,key); 
  cudaDeviceSynchronize();  
  
  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
