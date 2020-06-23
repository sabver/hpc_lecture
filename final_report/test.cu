#include <cuda_runtime.h>
#include <iostream>
int main()
{
  cudaDeviceProp prop;
 
  int count;
  cudaGetDeviceCount(&count);
  printf("显卡所支持的cuda处理器数量：%d\n", count);
  for (int i = 0; i < count; ++i){
    cudaGetDeviceProperties(&prop , i);
    printf("----第%d个处理器的基本信息----\n" ,i+1 );
    printf("处理器名称：%s \n" , prop.name );
    printf("计算能力：%d.%d\n" ,prop.major , prop.minor);
    printf("设备上全局内存总量：%dMB\n" ,prop.totalGlobalMem/1024/1024 );
    printf("设备上常量内存总量：%dKB\n", prop.totalConstMem/1024);
    printf("一个线程块中可使用的最大共享内存：%dKB\n", prop.sharedMemPerBlock / 1024);
    printf("一个线程束包含的线程数量：%d\n", prop.warpSize);
    printf("一个线程块中可包含的最大线程数量：%d\n", prop.maxThreadsPerBlock);
    printf("多维线程块数组中每一维可包含的最大线程数量：(%d,%d,%d)\n", prop.maxThreadsDim[0],
      prop.maxThreadsDim[1], prop.maxThreadsDim[2] );
    printf("一个线程格中每一维可包含的最大线程块数量：(%d,%d,%d)\n", prop.maxGridSize[0],
      prop.maxGridSize[1], prop.maxGridSize[2]);
  }
  return 0;
}