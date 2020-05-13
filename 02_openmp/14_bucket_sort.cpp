#include <cstdio>
#include <cstdlib>
#include <vector>

int main() {
  int n = 50;
  int range = 5;
  std::vector<int> key(n);
#pragma omp parallel for    
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

  // get the information of position of bucket's element
  std::vector<int> position(range+1);
  int cur = 0;
  for (int i=1; i<range; i++) {
    cur += bucket[i-1];  
    position[i] = cur;      
  }  
  position[range] = n;  
//   printf("position values: ");
//   for (int i=0; i<range+1; i++) {
//     printf("%d ",position[i]);
//   }    
//   printf("\n");
//   int start,end;
#pragma omp parallel for    
  for (int i=0; i<range; i++) {
//      start = position[i],end = position[i+1];
     //just copy
#pragma omp parallel for          
     for(int j=position[i];j<position[i+1];j++){
         key[j] = i;
     }
  }  
    
//   for (int i=0, j=0; i<range; i++) {
//     for (; bucket[i]>0; bucket[i]--) {
//       key[j++] = i;
//     }
//   }

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}

