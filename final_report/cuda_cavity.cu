#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <cooperative_groups.h>

using namespace std;

const static int U_FILE = 1;
const static int V_FILE = 2;
const static int P_FILE = 3;

const int row = 41;
const int col = 41;

const int BLOCK_NUM = row;
const int THREAD_NUM = col;

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

__device__ void zeros_gpu(double *vector,int row,int col){
  // printf("zeros_gpu\n");
  for(int i=0;i<row*col;i++){
    vector[i] = 0.0;
  }
  // printf("end\n");
}
//先假设是单线程，后面改为单个复制，这样子就没有调用的必要了
__device__ void copy(double *copy,double *origin){
  for(int i=0;i<nx * ny;i++){
    copy[i] = origin[i];
  } 
}

//先假设是单线程
__device__ void build_up_b(double *b,double *u, double *v){
  int row = ny,col = nx;
  //j-loop and i-loop have no Data dependence,so it can make Parallelization directly
  for(int j=1;j<row-1;j++){
    for(int i=1;i<col-1;i++){ 
      // b[j][i] = (rho * (1 / dt * 
      //               ((u[j][i+1] - u[j][i-1]) / 
      //                (2 * dx) + (v[j+1][i] - v[j-1][i]) / (2 * dy)) -
      //               pow(((u[j][i+1] - u[j][i-1]) / (2 * dx)),2) -
      //                 2 * ((u[j+1][i] - u[j-1][i]) / (2 * dy) *
      //                      (v[j][i+1] - v[j][i-1]) / (2 * dx))-
      //                     pow(((v[j+1][i] - v[j-1][i]) / (2 * dy)),2)));  
      b[index(j,i)] = (rho * (1 / dt * 
                    ((u[index(j,i+1)] - u[index(j,i-1)]) / 
                     (2 * dx) + (v[index(j+1,i)] - v[index(j-1,i)]) / (2 * dy)) -
                    pow(((u[index(j,i+1)] - u[index(j,i-1)]) / (2 * dx)),2) -
                      2 * ((u[index(j+1,i)] - u[index(j-1,i)]) / (2 * dy) *
                           (v[index(j,i+1)] - v[index(j,i-1)]) / (2 * dx))-
                          pow(((v[index(j+1,i)] - v[index(j-1,i)]) / (2 * dy)),2)));                                                                 
    }
  }
}

//先假设是单线程
__device__ void pressure_poisson(double *p,double *b){
  double *pn;
  pn = (double *)malloc(nx * ny * sizeof(double));  
  copy(pn,p);
  int row = ny,col = nx;
  //q-loop have Data dependence
  for(int q=0;q<nit;q++){    
    copy(pn,p);
    //j-loop and i-loop have no Data dependence,so it can make Parallelization directly
    for(int j=1;j<row-1;j++){
      for(int i=1;i<col-1;i++){        
        // p[j][i] = (((pn[j][i+1] + pn[j][i-1]) * pow(dy,2) + 
        //                   (pn[j+1][i] + pn[j-1][i]) * pow(dx,2)) /
        //                   (2 * (pow(dx,2) + pow(dy,2))) -
        //                   pow(dx,2) * pow(dy,2) / (2 * (pow(dx,2) + pow(dy,2))) * 
        //                   b[j][i]);      
        p[index(j,i)] = (((pn[index(j,i+1)] + pn[index(j,i-1)]) * pow(dy,2) + 
                          (pn[index(j+1,i)] + pn[index(j-1,i)]) * pow(dx,2)) /
                          (2 * (pow(dx,2) + pow(dy,2))) -
                          pow(dx,2) * pow(dy,2) / (2 * (pow(dx,2) + pow(dy,2))) * 
                          b[index(j,i)]);                                                                               
      }
    }  
    for(int j=0;j<row;j++){
      //p[:, -1] = p[:, -2] # dp/dx = 0 at x = 2
      // p[j][col-1] = p[j][col-2];
      p[index(j,col-1)] = p[index(j,col-2)];
    }
    for(int i=0;i<col;i++){
      //p[0, :] = p[1, :]   # dp/dy = 0 at y = 0
      // p[0][i] = p[1][i];
      p[index(0,i)] = p[index(1,i)];
    }
    for(int j=0;j<row;j++){
      //p[:, 0] = p[:, 1]   # dp/dx = 0 at x = 0
      // p[j][0] = p[j][1];
      p[index(j,0)] = p[index(j,1)];
    }
    for(int i=0;i<col;i++){
      //p[-1, :] = 0        # p = 0 at y = 2
      // p[row-1][i] = 0.0;
      p[index(row-1,i)] = 0.0;
    }
    //lock
  }
  free(pn);
  pn = NULL;
}

__device__ void cavity_flow(double *u, double *v, double *p,double *b){
  double *un,*vn;
  un = (double *)malloc(nx * ny * sizeof(double));   
  vn = (double *)malloc(nx * ny * sizeof(double));   
  // zeros_gpu(b,ny,nx);
  for(int n=0;n<nt;n++){
    copy(un,u);
    copy(vn,v);
    // change b
    build_up_b(b, u, v);
    // change p 
    pressure_poisson(p, b);
    //lock
    int row = ny,col = nx;
    //j-loop and i-loop have no Data dependence,so it can make Parallelization directly
    for(int j=1;j<row-1;j++){
      for(int i=1;i<col-1;i++){           
        // u[j][i] = (un[j][i]-
        //                  un[j][i] * dt / dx *
        //                 (un[j][i] - un[j][i-1]) -
        //                  vn[j][i] * dt / dy *
        //                 (un[j][i] - un[j-1][i]) -
        //                  dt / (2 * rho * dx) * (p[j][i+1] - p[j][i-1]) +
        //                  nu * (dt / pow(dx,2) *
        //                 (un[j][i+1] - 2 * un[j][i] + un[j][i-1]) +
        //                  dt / pow(dy,2) *
        //                 (un[j+1][i] - 2 * un[j][i] + un[j-1][i])));        
        u[index(j,i)] = (un[index(j,i)]-
                         un[index(j,i)] * dt / dx *
                        (un[index(j,i)] - un[index(j,i-1)]) -
                         vn[index(j,i)] * dt / dy *
                        (un[index(j,i)] - un[index(j-1,i)]) -
                         dt / (2 * rho * dx) * (p[index(j,i+1)] - p[index(j,i-1)]) +
                         nu * (dt / pow(dx,2) *
                        (un[index(j,i+1)] - 2 * un[index(j,i)] + un[index(j,i-1)]) +
                         dt / pow(dy,2) *
                        (un[index(j+1,i)] - 2 * un[index(j,i)] + un[index(j-1,i)])));                                  
      }
    }  
    //j-loop and i-loop have no Data dependence,so it can make Parallelization directly
    for(int j=1;j<row-1;j++){
      for(int i=1;i<col-1;i++){ 
        // v[j][i] = (vn[j][i] -
        //                 un[j][i] * dt / dx *
        //                (vn[j][i] - vn[j][i-1]) -
        //                 vn[j][i] * dt / dy *
        //                (vn[j][i] - vn[j-1][i]) -
        //                 dt / (2 * rho * dy) * (p[j+1][i] - p[j-1][i]) +
        //                 nu * (dt / pow(dx,2) *
        //                (vn[j][i+1] - 2 * vn[j][i] + vn[j][i-1]) +
        //                 dt / pow(dy,2) *
        //                (vn[j+1][i] - 2 * vn[j][i] + vn[j-1][i])));                
        v[index(j,i)] = (vn[index(j,i)] -
                        un[index(j,i)] * dt / dx *
                       (vn[index(j,i)] - vn[index(j,i-1)]) -
                        vn[index(j,i)] * dt / dy *
                       (vn[index(j,i)] - vn[index(j-1,i)]) -
                        dt / (2 * rho * dy) * (p[index(j+1,i)] - p[index(j-1,i)]) +
                        nu * (dt / pow(dx,2) *
                       (vn[index(j,i+1)] - 2 * vn[index(j,i)] + vn[index(j,i-1)]) +
                        dt / pow(dy,2) *
                       (vn[index(j+1,i)] - 2 * vn[index(j,i)] + vn[index(j-1,i)])));        
      }
    }       
    for(int i=0;i<col;i++){
      // u[0, :]  = 0
      // u[0][i] = 0;
      u[index(0,i)] = 0;
    }  
    for(int j=0;j<row;j++){
      // u[:, 0]  = 0
      // u[j][0] = 0;
      u[index(j,0)] = 0;
    }
    for(int j=0;j<row;j++){
      // u[:, -1] = 0
      // u[j][col-1] = 0;
      u[index(j,col-1)] = 0;
    }
    for(int i=0;i<col;i++){
      // u[-1, :] = 1    # set velocity on cavity lid equal to 1
      // u[row-1][i] = 1;
      u[index(row-1,i)] = 1;
    }
    for(int i=0;i<col;i++){
      // v[0, :]  = 0
      // v[0][i] = 0;
      v[index(0,i)] = 0;
    }
    for(int i=0;i<col;i++){
      // v[-1, :] = 0
      // v[row-1][i] = 0;
      v[index(row-1,i)] = 0;
    }
    for(int j=0;j<row;j++){
      // v[:, 0]  = 0
      // v[j][0] = 0;
      v[index(j,0)] = 0;
    }
    for(int j=0;j<row;j++){
      // v[:, -1] = 0  
      // v[j][col-1] = 0;
      v[index(j,col-1)] = 0;
    }  
    //lock  
  }
  free(un);
  free(vn);
  un = NULL;
  vn = NULL;
}

__global__ void kernel(double *u,double *v,double *p){
  int blockId = blockIdx.x;
  int threadId = threadIdx.x;
  if( blockId < 1 && threadId < 1 ){
    double *b;
    b = (double *)malloc(nx * ny * sizeof(double));
    int test = blockId == 0 && threadId == 0;
    printf("blockId:%d,threadId:%d,result:%d\n",blockId,threadId,test);
    zeros_gpu(b,ny,nx);  
    printf("---------\n");
    cavity_flow(u,v,p,b);
    free(b);
    b = NULL;
  }    
}

int main(){
  //2-d
  // double *u,*v,*p,*b,*un,*vn,*pn,*bn;
  double *u,*v,*p;
  cudaMallocManaged(&u, row*col*sizeof(double));
  cudaMallocManaged(&v, row*col*sizeof(double));
  cudaMallocManaged(&p, row*col*sizeof(double));

  zeros(u,row,col);
  zeros(v,row,col);
  zeros(p,row,col);

  void *args[] = {(void *)&u,  (void *)&v, (void *)&p};
  cudaLaunchCooperativeKernel((void*)kernel, BLOCK_NUM, THREAD_NUM, args); 
  
  // cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu);
  cudaDeviceSynchronize();  
  string u_json = array_2d_to_json(u,row,col),
         v_json = array_2d_to_json(v,row,col),
         p_json = array_2d_to_json(p,row,col);
  write_string_to_file(u_json,U_FILE);
  write_string_to_file(v_json,V_FILE);
  write_string_to_file(p_json,P_FILE);
  cudaFree(u);
  cudaFree(v);
  cudaFree(p);
}