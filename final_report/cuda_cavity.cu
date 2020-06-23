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
// 1024 = > CUDA error: too many blocks in cooperative launch
// https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#grid-synchronization-cg
int THREAD_NUM = 512;

// const int BLOCK_NUM = row;
// const int THREAD_NUM = col;

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

// __device__ __managed__ int count;

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
  if( j*nx+i >= nx * ny ){
    printf("over index:%d,%d\n",j,i);
  }
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
//先假设是单线程，后面改为单个复制，这样子就没有调用的必要了,这里会不会有死锁的问题
__device__ void copy(double *copy,double *origin){
  // printf("copy\n");
  for(int i=0;i<nx * ny;i++){
    copy[i] = origin[i];
  } 
}

__device__ void judge_over_index(int index){
  if( index >= nx * ny ){
    printf("over index:%d,%d\n",index);
  }  
}

__device__ void build_up_b_single_thread(cooperative_groups::grid_group grid,double *b,double *u, double *v){
  // printf("build_up_b:%d\n",threadId);
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

//先假设是单线程，这里不考虑边界越界
__device__ void build_up_b(cooperative_groups::grid_group grid,int threadId,double *b,double *u, double *v){
  // if( blockId < 1 && threadId < 1 ){
  //   printf("build_up_b\n");  
  // }   
  // printf("build_up_b:%d\n",threadId);
  int row = ny,col = nx;
  //j-loop and i-loop have no Data dependence,so it can make Parallelization directly
  b[threadId] = (rho * (1 / dt * 
                    ((u[threadId+1] - u[threadId-1]) / 
                     (2 * dx) + (v[threadId+col] - v[threadId-col]) / (2 * dy)) -
                    pow(((u[threadId+1] - u[threadId-1]) / (2 * dx)),2) -
                      2 * ((u[threadId+col] - u[threadId-col]) / (2 * dy) *
                           (v[threadId+1] - v[threadId-1]) / (2 * dx))-
                          pow(((v[threadId+col] - v[threadId-col]) / (2 * dy)),2)));    
  // for(int j=1;j<row-1;j++){
  //   for(int i=1;i<col-1;i++){ 
  //     // b[j][i] = (rho * (1 / dt * 
  //     //               ((u[j][i+1] - u[j][i-1]) / 
  //     //                (2 * dx) + (v[j+1][i] - v[j-1][i]) / (2 * dy)) -
  //     //               pow(((u[j][i+1] - u[j][i-1]) / (2 * dx)),2) -
  //     //                 2 * ((u[j+1][i] - u[j-1][i]) / (2 * dy) *
  //     //                      (v[j][i+1] - v[j][i-1]) / (2 * dx))-
  //     //                     pow(((v[j+1][i] - v[j-1][i]) / (2 * dy)),2)));  
  //     b[index(j,i)] = (rho * (1 / dt * 
  //                   ((u[index(j,i+1)] - u[index(j,i-1)]) / 
  //                    (2 * dx) + (v[index(j+1,i)] - v[index(j-1,i)]) / (2 * dy)) -
  //                   pow(((u[index(j,i+1)] - u[index(j,i-1)]) / (2 * dx)),2) -
  //                     2 * ((u[index(j+1,i)] - u[index(j-1,i)]) / (2 * dy) *
  //                          (v[index(j,i+1)] - v[index(j,i-1)]) / (2 * dx))-
  //                         pow(((v[index(j+1,i)] - v[index(j-1,i)]) / (2 * dy)),2)));                                                                 
  //   }
  // }
}

__device__ void pressure_poisson_single_thread(cooperative_groups::grid_group grid,double *p,double *b,double *pn){ 
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
}

//先假设是单线程
__device__ void pressure_poisson(cooperative_groups::grid_group grid,int threadId,double *p,double *b,double *pn){
  // if( blockId < 1 && threadId < 1 ){
  //   printf("pressure_poisson\n");  
  // }    
  // printf("pressure_poisson:%d\n",threadId);
  // double *pn;
  // pn = (double *)malloc(nx * ny * sizeof(double));  
  // copy(pn,p);
  int row = ny,col = nx;
  //q-loop have Data dependence
  for(int q=0;q<nit;q++){    
    // copy(pn,p);  
    pn[threadId] = p[threadId];   
    grid.sync();
    p[threadId] = (((pn[threadId+1] + pn[threadId-1]) * pow(dy,2) + 
                          (pn[threadId+col] + pn[threadId-col]) * pow(dx,2)) /
                          (2 * (pow(dx,2) + pow(dy,2))) -
                          pow(dx,2) * pow(dy,2) / (2 * (pow(dx,2) + pow(dy,2))) * 
                          b[threadId]);      
    //j-loop and i-loop have no Data dependence,so it can make Parallelization directly
    // for(int j=1;j<row-1;j++){
    //   for(int i=1;i<col-1;i++){        
    //     // p[j][i] = (((pn[j][i+1] + pn[j][i-1]) * pow(dy,2) + 
    //     //                   (pn[j+1][i] + pn[j-1][i]) * pow(dx,2)) /
    //     //                   (2 * (pow(dx,2) + pow(dy,2))) -
    //     //                   pow(dx,2) * pow(dy,2) / (2 * (pow(dx,2) + pow(dy,2))) * 
    //     //                   b[j][i]);      
    //     p[index(j,i)] = (((pn[index(j,i+1)] + pn[index(j,i-1)]) * pow(dy,2) + 
    //                       (pn[index(j+1,i)] + pn[index(j-1,i)]) * pow(dx,2)) /
    //                       (2 * (pow(dx,2) + pow(dy,2))) -
    //                       pow(dx,2) * pow(dy,2) / (2 * (pow(dx,2) + pow(dy,2))) * 
    //                       b[index(j,i)]);                                                                               
    //   }
    // }  
    grid.sync();
    for(int j=0;j<row;j++){
      //p[:, -1] = p[:, -2] # dp/dx = 0 at x = 2
      // p[j][col-1] = p[j][col-2];
      p[index(j,col-1)] = p[index(j,col-2)];
    }
    grid.sync();
    for(int i=0;i<col;i++){
      //p[0, :] = p[1, :]   # dp/dy = 0 at y = 0
      // p[0][i] = p[1][i];
      p[index(0,i)] = p[index(1,i)];
    }
    grid.sync();
    for(int j=0;j<row;j++){
      //p[:, 0] = p[:, 1]   # dp/dx = 0 at x = 0
      // p[j][0] = p[j][1];
      p[index(j,0)] = p[index(j,1)];
    }
    grid.sync();
    for(int i=0;i<col;i++){
      //p[-1, :] = 0        # p = 0 at y = 2
      // p[row-1][i] = 0.0;
      p[index(row-1,i)] = 0.0;
    }
    grid.sync();
    //lock
  }
  // free(pn);
  // pn = NULL;
}
__device__ void cavity_flow_single_thread(cooperative_groups::grid_group grid,int threadId, double *u, double *v, double *p,double *b,double *un,double *vn,double *pn){
  int row = ny,col = nx;  
  // zeros_gpu(b,ny,nx);
  for(int n=0;n<nt;n++){ 
    // copy
    if( threadId == 60 ){      
      copy(un,u);
      copy(vn,v);
      copy(pn,p);
      //边界判断
      if( threadId / col == 0 || threadId / col == row-1 || threadId % col == 0 || threadId % col == col-1 ){
        continue;
      }
      build_up_b_single_thread(grid,b,u,v);
      pressure_poisson_single_thread(grid,p, b,pn);
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
    }                 
  }
  printf("end:%d\n",threadId);
}
__device__ void cavity_flow(cooperative_groups::grid_group grid,int threadId, double *u, double *v, double *p,double *b,double *un,double *vn,double *pn){
  // if( threadId > 1599 ){
  //   printf("cavity_flow:%d,total:%d\n",threadId < nx * ny,nx * ny);
  // }  
  // if( blockId < 1 && threadId < 1 ){
  //   printf("cavity_flow\n");  
  // }  
  // printf("cavity_flow\n");  
  // double *un,*vn;
  // un = (double *)malloc(nx * ny * sizeof(double));   
  // vn = (double *)malloc(nx * ny * sizeof(double)); 
  int row = ny,col = nx;  
  // zeros_gpu(b,ny,nx);
  for(int n=0;n<nt;n++){
    // printf("threadId:%d,n:%d\n",threadId,n);
    // atomicAdd(&count, 1);
    // printf("%d\n",n);
    // if( blockId < 1 && threadId < 1 ){
    //   printf("cavity_flow\n");  
    // }   
    // copy   
    un[threadId] = u[threadId];
    vn[threadId] = v[threadId];   
    pn[threadId] = p[threadId];   
    // // printf("first\n");
    grid.sync();
    // printf("second\n");
    //边界判断
    if( threadId / col == 0 || threadId / col == row-1 || threadId % col == 0 || threadId % col == col-1 ){
      // printf("threadId:%d\n",threadId);
      continue;
    }
    // copy(un,u);
    // copy(vn,v);
    // printf("---\n");
    // change b
    build_up_b(grid,threadId,b, u, v);
    grid.sync();
    // change p 
    pressure_poisson(grid,threadId,p, b,pn);
    grid.sync();    
    //lock
    //j-loop and i-loop have no Data dependence,so it can make Parallelization directly
    // for(int j=1;j<row-1;j++){
    //   for(int i=1;i<col-1;i++){           
    //     // u[j][i] = (un[j][i]-
    //     //                  un[j][i] * dt / dx *
    //     //                 (un[j][i] - un[j][i-1]) -
    //     //                  vn[j][i] * dt / dy *
    //     //                 (un[j][i] - un[j-1][i]) -
    //     //                  dt / (2 * rho * dx) * (p[j][i+1] - p[j][i-1]) +
    //     //                  nu * (dt / pow(dx,2) *
    //     //                 (un[j][i+1] - 2 * un[j][i] + un[j][i-1]) +
    //     //                  dt / pow(dy,2) *
    //     //                 (un[j+1][i] - 2 * un[j][i] + un[j-1][i])));        
    //     u[index(j,i)] = (un[index(j,i)]-
    //                      un[index(j,i)] * dt / dx *
    //                     (un[index(j,i)] - un[index(j,i-1)]) -
    //                      vn[index(j,i)] * dt / dy *
    //                     (un[index(j,i)] - un[index(j-1,i)]) -
    //                      dt / (2 * rho * dx) * (p[index(j,i+1)] - p[index(j,i-1)]) +
    //                      nu * (dt / pow(dx,2) *
    //                     (un[index(j,i+1)] - 2 * un[index(j,i)] + un[index(j,i-1)]) +
    //                      dt / pow(dy,2) *
    //                     (un[index(j+1,i)] - 2 * un[index(j,i)] + un[index(j-1,i)])));                                  
    //   }
    // }  
    u[threadId] = (un[threadId]-
                         un[threadId] * dt / dx *
                        (un[threadId] - un[threadId-1]) -
                         vn[threadId] * dt / dy *
                        (un[threadId] - un[threadId-col]) -
                         dt / (2 * rho * dx) * (p[threadId+1] - p[threadId-1]) +
                         nu * (dt / pow(dx,2) *
                        (un[threadId+1] - 2 * un[threadId] + un[threadId-1]) +
                         dt / pow(dy,2) *
                        (un[threadId+col] - 2 * un[threadId] + un[threadId-col])));  
    grid.sync();   
    //j-loop and i-loop have no Data dependence,so it can make Parallelization directly
    // for(int j=1;j<row-1;j++){
    //   for(int i=1;i<col-1;i++){ 
    //     // v[j][i] = (vn[j][i] -
    //     //                 un[j][i] * dt / dx *
    //     //                (vn[j][i] - vn[j][i-1]) -
    //     //                 vn[j][i] * dt / dy *
    //     //                (vn[j][i] - vn[j-1][i]) -
    //     //                 dt / (2 * rho * dy) * (p[j+1][i] - p[j-1][i]) +
    //     //                 nu * (dt / pow(dx,2) *
    //     //                (vn[j][i+1] - 2 * vn[j][i] + vn[j][i-1]) +
    //     //                 dt / pow(dy,2) *
    //     //                (vn[j+1][i] - 2 * vn[j][i] + vn[j-1][i])));                
    //     v[index(j,i)] = (vn[index(j,i)] -
    //                     un[index(j,i)] * dt / dx *
    //                    (vn[index(j,i)] - vn[index(j,i-1)]) -
    //                     vn[index(j,i)] * dt / dy *
    //                    (vn[index(j,i)] - vn[index(j-1,i)]) -
    //                     dt / (2 * rho * dy) * (p[index(j+1,i)] - p[index(j-1,i)]) +
    //                     nu * (dt / pow(dx,2) *
    //                    (vn[index(j,i+1)] - 2 * vn[index(j,i)] + vn[index(j,i-1)]) +
    //                     dt / pow(dy,2) *
    //                    (vn[index(j+1,i)] - 2 * vn[index(j,i)] + vn[index(j-1,i)])));        
    //   }
    // }   
    v[threadId] = (vn[threadId] -
                        un[threadId] * dt / dx *
                       (vn[threadId] - vn[threadId-1]) -
                        vn[threadId] * dt / dy *
                       (vn[threadId] - vn[threadId-col]) -
                        dt / (2 * rho * dy) * (p[threadId+col] - p[threadId-col]) +
                        nu * (dt / pow(dx,2) *
                       (vn[threadId+1] - 2 * vn[threadId] + vn[threadId-1]) +
                        dt / pow(dy,2) *
                       (vn[threadId+col] - 2 * vn[threadId] + vn[threadId-col])));         
    grid.sync();   
    for(int i=0;i<col;i++){
      // u[0, :]  = 0
      // u[0][i] = 0;
      u[index(0,i)] = 0;
    }  
    grid.sync();
    for(int j=0;j<row;j++){
      // u[:, 0]  = 0
      // u[j][0] = 0;
      u[index(j,0)] = 0;
    }
    grid.sync();
    for(int j=0;j<row;j++){
      // u[:, -1] = 0
      // u[j][col-1] = 0;
      u[index(j,col-1)] = 0;
    }
    grid.sync();
    for(int i=0;i<col;i++){
      // u[-1, :] = 1    # set velocity on cavity lid equal to 1
      // u[row-1][i] = 1;
      u[index(row-1,i)] = 1;
    }
    grid.sync();
    for(int i=0;i<col;i++){
      // v[0, :]  = 0
      // v[0][i] = 0;
      v[index(0,i)] = 0;
    }
    grid.sync();
    for(int i=0;i<col;i++){
      // v[-1, :] = 0
      // v[row-1][i] = 0;
      v[index(row-1,i)] = 0;
    }
    grid.sync();
    for(int j=0;j<row;j++){
      // v[:, 0]  = 0
      // v[j][0] = 0;
      v[index(j,0)] = 0;
    }
    grid.sync();
    for(int j=0;j<row;j++){
      // v[:, -1] = 0  
      // v[j][col-1] = 0;
      v[index(j,col-1)] = 0;
    } 
    grid.sync();
    //lock  
  }
  // free(un);
  // free(vn);
  // un = NULL;
  // vn = NULL;
  printf("end:%d\n",threadId);
}

__global__ void kernel(double *u,double *v,double *p,double *b,double *un,double *vn,double *pn){
  int threadId = blockIdx.x * blockDim.x + threadIdx.x;  
  if( threadId < nx * ny ){
    // printf("threadId:%d\n",threadId);
    cooperative_groups::grid_group grid = cooperative_groups::this_grid();
    // cavity_flow(grid,threadId,u,v,p,b,un,vn,pn);  
    cavity_flow(grid,threadId,u,v,p,b,un,vn,pn);
  } 
}
//nvcc cuda_cavity.cu  -arch=sm_60 -rdc=true
int main(){
  //2-d
  // double *u,*v,*p,*b,*un,*vn,*pn,*bn;
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

  void *args[] = {(void *)&u,  (void *)&v, (void *)&p, (void *)&b,(void *)&u,  (void *)&v, (void *)&p};  
  int total = row * col;
  
  int BLOCK_NUM = total / THREAD_NUM ;
  if( total % THREAD_NUM > 0 ){
    BLOCK_NUM++;
  }
  printf("total:%d\n",total);
  printf("BLOCK_NUM:%d\n",BLOCK_NUM);
  printf("THREAD_NUM:%d\n",THREAD_NUM);
  cudaLaunchCooperativeKernel((void*)kernel, BLOCK_NUM, THREAD_NUM, args); 
  cudaError_t error = cudaGetLastError();
  printf("CUDA error: %s\n", cudaGetErrorString(error));
  
  cudaDeviceSynchronize();  
  cudaError_t error2 = cudaGetLastError();
  printf("CUDA error: %s\n", cudaGetErrorString(error2));  
  // printf("count:%d\n",count);
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