#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>

using namespace std;

const static int U_FILE = 1;
const static int V_FILE = 2;
const static int P_FILE = 3;

const int row = 41;
const int col = 41;

const int THREAD_NUM = 1024;

const double dx_cpu = 2.0 / (col - 1);
const double dy_cpu = 2.0 / (row - 1);

__device__ int nx = col; 
__device__ int ny = row; 
__device__ int nt = 500;
__device__ int nit = 50;
__device__ double c = 1.0;
// double dx = 2.0 / (nx - 1);
__device__ double dx = dx_cpu;
// double dy = 2.0 / (ny - 1);
__device__ double dy = dy_cpu;

__device__ double rho = 1.0;
__device__ double nu = 0.1;
__device__ double dt = 0.001;

__device__ int TOTAL_GPU = row * col;
__device__ int PER_GPU = 0;


template < typename Type > std::string to_str (const Type & t)
{
  std::ostringstream os;
  os << t;
  return os.str ();
}


//2-d
string array_2d_to_json(double *vector,int row,int col){
  string result = "[";
  for(int i=0;i<row;i++){
    result += "[";
    for(int j=0;j<col;j++){
      result += to_str(vector[i*col+j]);
      if( j!=col-1 ){
        result += ",";
      }
    }
    result += "]";
    if( i!=row-1 ){
      result += ",";
    }
  }  
  return result+"]";  
}

void print_array(double *vector,int row,int col){
    for(int i=0;i<row;i++){
        for(int j=0;j<col;j++){
            printf("%f,",vector[i*col+j]);
        }
        printf("\n");
    }
}

void write_string_to_file(string str,int flag){
  ofstream outfile;
  switch(flag){
    case U_FILE:
      outfile.open("/home/9/19M38171/t3workspace/hpc-lecture/final_report/u_json.txt");  
      break;
    case V_FILE:
      outfile.open("/home/9/19M38171/t3workspace/hpc-lecture/final_report/v_json.txt");  
      break;
    case P_FILE:
      outfile.open("/home/9/19M38171/t3workspace/hpc-lecture/final_report/p_json.txt");  
      break;    
    default:
      break;  
  }
  outfile << str << endl;
  outfile.close();
}


void zeros(double *vector,int row,int col){
  for(int i=0;i<row*col;i++){
    vector[i] = 0.0;
  }
}
//2-d's index to 1-d's index
__device__ int index(int j,int i){
  //nx is col
  return j*nx+i;
}

// __device__ void zeros_gpu(double *vector,int row,int col){
//   // printf("zeros_gpu\n");
//   for(int i=0;i<row*col;i++){
//     vector[i] = 0.0;
//   }
//   // printf("end\n");
// }

//先假设是单线程，后面改为单个复制，这样子就没有调用的必要了
// __device__ void copy(double *copy,double *origin){
//   for(int i=0;i<nx * ny;i++){
//     copy[i] = origin[i];
//   } 
// }


__device__ void build_up_b(double *b,double *u, double *v,int startIndex,int endIndex){
  int row = ny,col = nx;
  for(int i=startIndex;i<=endIndex;i++){
      if( i / col == 0 || i / col == row - 1 || i % col == 0 || i % col == col-1 ){
        continue;
      }  
      b[i] = (rho * (1 / dt * 
                    ((u[i+1] - u[i-1]) / 
                     (2 * dx) + (v[i+col] - v[i-col]) / (2 * dy)) -
                    pow(((u[i+1] - u[i-1]) / (2 * dx)),2) -
                      2 * ((u[i+col] - u[i-col]) / (2 * dy) *
                           (v[i+1] - v[i-1]) / (2 * dx))-
                          pow(((v[i+col] - v[i-col]) / (2 * dy)),2)));         
  }
}


__device__ void pressure_poisson(double *p, double *b, double *pn,int startIndex,int endIndex){
  // copy(pn,p);
  for(int i=startIndex;i<=endIndex;i++){
    pn[i] = p[i];
  }  
  __syncthreads(); 
  int row = ny,col = nx;
  //q-loop have Data dependence
  for(int q=0;q<nit;q++){    
    // copy(pn,p);
    for(int i=startIndex;i<=endIndex;i++){
      pn[i] = p[i];
    }      
    __syncthreads(); 
    for(int i=startIndex;i<=endIndex;i++){
      if( !(i / col == 0 || i / col == row - 1 || i % col == 0 || i % col == col-1) ){
        p[i] = (((pn[i+1] + pn[i-1]) * pow(dy,2) + 
                            (pn[i+col] + pn[i-col]) * pow(dx,2)) /
                            (2 * (pow(dx,2) + pow(dy,2))) -
                            pow(dx,2) * pow(dy,2) / (2 * (pow(dx,2) + pow(dy,2))) * 
                            b[i]); 
      }    
      // p[j][col-1] = p[j][col-2]
      int rowIndex = i / col,colIndex = i % col;
      p[index(rowIndex,col-1)] = p[index(rowIndex,col-2)];
      p[index(0,colIndex)] = p[index(1,colIndex)];
      p[index(rowIndex,0)] = p[index(rowIndex,1)];
      p[index(row-1,colIndex)] = 0.0;
    }
    __syncthreads();     
  }
}

__device__ void cavity_flow(double *u,double *v,double *p,double *b,double *un,double *vn,double *pn,int startIndex,int endIndex){
  // zeros_gpu(b,ny,nx);
  for(int n=0;n<nt;n++){
    // copy(un,u);
    // copy(vn,v);
    for(int i=startIndex;i<=endIndex;i++){
      un[i] = u[i];
      vn[i] = v[i];
    }
    __syncthreads(); 
    // change b
    build_up_b(b, u, v,startIndex,endIndex);
    __syncthreads(); 
    // change p 
    pressure_poisson(p, b, pn,startIndex,endIndex);
    // __syncthreads(); 
    //lock
    int row = ny,col = nx;    
    for(int i=startIndex;i<=endIndex;i++){
      if( i / col == 0 ){
         // u[0][i] = 0;
        u[i] = 0;
        // v[0][i] = 0;
        v[i] = 0;
      }
      if( i % col == 0 ){
        // u[j][0] = 0;
        u[i] = 0;
        // v[j][0] = 0;
        v[i] = 0;
      }      
      if( i % col == col-1 ){
        // u[j][col-1] = 0;
        u[i] = 0;
        // v[j][col-1] = 0;
        v[i] = 0;

      }
      if( i / col == row-1 ){
        // u[row-1][i] = 1;
        u[i] = 1;
        // v[row-1][i] = 0;
        v[i] = 0;
      }      
      if( i / col == 0 || i / col == row - 1 || i % col == 0 || i % col == col-1 ){
        continue;
      }    
      u[i] = (un[i]-
              un[i] * dt / dx *
             (un[i] - un[i-1]) -
              vn[i] * dt / dy *
             (un[i] - un[i-col]) -
              dt / (2 * rho * dx) * (p[i+1] - p[i-1]) +
              nu * (dt / pow(dx,2) *
             (un[i+1] - 2 * un[i] + un[i-1]) +
              dt / pow(dy,2) *
             (un[i+col] - 2 * un[i] + un[i-col]))); 
      v[i] = (vn[i] -
              un[i] * dt / dx *
             (vn[i] - vn[i-1]) -
              vn[i] * dt / dy *
             (vn[i] - vn[i-col]) -
              dt / (2 * rho * dy) * (p[i+col] - p[i-col]) +
              nu * (dt / pow(dx,2) *
             (vn[i+1] - 2 * vn[i] + vn[i-1]) +
              dt / pow(dy,2) *
             (vn[i+col] - 2 * vn[i] + vn[i-col])));       

    }
    __syncthreads();
  }
}

__global__ void kernel(double *u,double *v,double *p,double *b,double *un,double *vn,double *pn){  
  int threadId = blockIdx.x * blockDim.x + threadIdx.x;  
  if( threadId < TOTAL_GPU ){
    if( TOTAL_GPU % THREAD_NUM == 0 ){
      PER_GPU = TOTAL_GPU/THREAD_NUM;
    }else{
      PER_GPU = (TOTAL_GPU/THREAD_NUM) + 1;
    }
    //if PER_GPU = 0 => 
    int startIndex = threadId*PER_GPU,endIndex = startIndex+PER_GPU-1;
    // if( PER_GPU == 0 ){
    //   startIndex = endIndex = threadId;
    // }   
    if( startIndex >= TOTAL_GPU ){
      startIndex = TOTAL_GPU - 1;
    }
    if( endIndex >= TOTAL_GPU ){
      endIndex = TOTAL_GPU - 1;
    }
    // printf("(%d,%d)(%d,%d)\n",startIndex,endIndex,threadId,PER_GPU);
    cavity_flow(u,v,p,b,un,vn,pn,startIndex,endIndex);
  }    
}

int main(){
  //2-d
  double *u,*v,*p,*b,*un,*vn,*pn;
  cudaMallocManaged(&u, row*col*sizeof(double));
  cudaMallocManaged(&v, row*col*sizeof(double));
  cudaMallocManaged(&p, row*col*sizeof(double));
  cudaMallocManaged(&b, row*col*sizeof(double));
  cudaMallocManaged(&un, row*col*sizeof(double));
  cudaMallocManaged(&vn, row*col*sizeof(double));
  cudaMallocManaged(&pn, row*col*sizeof(double));

  zeros(u,row,col);
  zeros(v,row,col);
  zeros(p,row,col);
  zeros(b,row,col);

  kernel<<<1,THREAD_NUM>>>(u,v,p,b,un,vn,pn);
  cudaError_t error = cudaGetLastError();
  printf("CUDA error: %s\n", cudaGetErrorString(error));

  cudaDeviceSynchronize();  
  cudaError_t error2 = cudaGetLastError();
  printf("CUDA error: %s\n", cudaGetErrorString(error2));   
  string u_json = array_2d_to_json(u,row,col),
         v_json = array_2d_to_json(v,row,col),
         p_json = array_2d_to_json(p,row,col);
  write_string_to_file(u_json,U_FILE);
  write_string_to_file(v_json,V_FILE);
  write_string_to_file(p_json,P_FILE);
  cudaFree(u);
  cudaFree(v);
  cudaFree(p);
  cudaFree(b);
  cudaFree(un);
  cudaFree(vn);
  cudaFree(pn);
}