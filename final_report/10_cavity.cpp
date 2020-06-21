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
const static int B_FILE = 4;

int nx = 41; //col
int ny = 41; //row
int nt = 500;
int nit = 50;
double c = 1.0;
double dx = 2.0 / (nx - 1);
double dy = 2.0 / (ny - 1);
// x = numpy.linspace(0, 2, nx)
// y = numpy.linspace(0, 2, ny)
// X, Y = numpy.meshgrid(x, y)

double rho = 1.0;
double nu = 0.1;
double dt = 0.001;

vector< vector<double> > u(ny); //ny * nx
vector< vector<double> > v(ny);
vector< vector<double> > p(ny);
vector< vector<double> > b(ny);
vector< vector<double> > un(ny); //ny * nx
vector< vector<double> > vn(ny);
vector< vector<double> > pn(ny);
vector< vector<double> > bn(ny);

template < typename Type > std::string to_str (const Type & t)
{
  std::ostringstream os;
  os << t;
  return os.str ();
}

/**
 * s -> "s"
 * @param  s [description]
 * @return   [description]
 */
string add_s(string s){
  return "\""+s+"\"";
}


/**
 * for example,{"0":"[0,0,0,0,0]","1":"[1,1,1,1,1]","2":"[2,2,2,2,2]"}
 * in python:
 * import json
 * myarray = json.loads(jsonstr)
 * print(myarray[str(0)]) = > [0,0,0,0,0]
 * @param  vector [description]
 * @return        [description]
 */
string to_json(vector<string> vector,int len){
  string result = "{";
  for(int i=0;i<len;i++){
    result += ( add_s(to_str(i)) + ":" + vector[i]) ;
    if( i!=len-1 ){
      result += ",";
    }
  }
  return result+"}";
}

//1-d
string array_to_json( vector<double> vector,int len){
  string result = "[";
  for(int i=0;i<len;i++){
    result += to_str(vector[i]);
    if( i!=len-1 ){
      result += ",";
    }
  }
  return result+"]";
}


//2-d
string array_2d_to_json(vector< vector<double> > vector,int row,int col){
  string result = "[";
  for(int i=0;i<row;i++){
    result += array_to_json(vector[i],col);
    if( i!=row-1 ){
      result += ",";
    }
  }  
  return result+"]";  
}


void zeros(vector< vector<double> > &vector){
  for(int y = 0;y<ny;y++){
    for(int x = 0;x<nx;x++){
      vector[y][x] = 0.0;
    }
  }
}

void print_array(vector< vector<double> > vector){
    for(int i=0;i<vector.size();i++){
        for(int j=0;j<vector[i].size();j++){
            printf("%f,",vector[i][j]);
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
    case B_FILE:
      outfile.open("/home/9/19M38171/t3workspace/hpc-lecture/final_report/b_json.txt");  
      break;
    default:
      break;  
  }
  outfile << str << endl;
  outfile.close();
}

/**
 * u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2
 */
void set_value(vector< vector<double> > &vector,double dy_start_point,double dy_end_point,double dx_start_point,double dx_end_point, double value){  
  double dy_start = dy_start_point/dy,
        dy_end = dy_end_point/dy + 1,
        dx_start = dx_start_point/dx,        
        dx_end = dx_end_point/dx + 1;

  for(int i=(int)dy_start;i<=(int)dy_end;i++){
    for(int j=(int)dx_start;j<=(int)dx_end;j++){
      u[i][j] = value;
    }
  }
}

void set_value2(vector< vector<double> > &vector , int row_start,int row_end,int col_start,int col_end,double value ){
  for(int i=row_start;i<=row_end;i++){
    for(int j=col_start;j<=col_end;j++){
      vector[i][j] = value;
    }
  }
}

void copy(vector< vector<double> > &copy,vector< vector<double> > origin){
    for(int i=0;i<origin.size();i++){
      for(int j=0;j<origin[i].size();j++){
          copy[i][j] = origin[i][j];
      }
    }  
}

void build_up_b(vector< vector<double> > &b,double rho, double dt, vector< vector<double> > u, vector< vector<double> > v, double dx, double dy){
  // b[1:-1, 1:-1] = (rho * (1 / dt * 
  //                   ((u[1:-1, 2:] - u[1:-1, 0:-2]) / 
  //                    (2 * dx) + (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
  //                   ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 -
  //                     2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
  //                          (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx))-
  //                         ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2));  
  // [j][i] -> [1:-1, 1:-1]
  // [j][i+1] -> [1:-1, 2:]
  // [j][i-1] -> [1:-1, 0:-2]
  // [j+1][i] -> [2:, 1:-1]
  // [j-1][i] -> [0:-2, 1:-1]
  int row = ny,col = nx;
  for(int j=1;j<row-1;j++){
    for(int i=1;i<col-1;i++){ 
      b[j][i] = (rho * (1 / dt * 
                    ((u[j][i+1] - u[j][i-1]) / 
                     (2 * dx) + (v[j+1][i] - v[j-1][i]) / (2 * dy)) -
                    pow(((u[j][i+1] - u[j][i-1]) / (2 * dx)),2) -
                      2 * ((u[j+1][i] - u[j-1][i]) / (2 * dy) *
                           (v[j][i+1] - v[j][i-1]) / (2 * dx))-
                          pow(((v[j+1][i] - v[j-1][i]) / (2 * dy)),2)));                                  
    }
  }
  // printf("print b:\n");
  // print_array(b);
  // printf("\n");
}

void pressure_poisson(vector< vector<double> > &p,double dx,double dy,vector< vector<double> > b){
  copy(pn,p);
  int row = ny,col = nx;
  for(int q=0;q<nit;q++){    
    copy(pn,p);
    for(int j=1;j<row-1;j++){
      for(int i=1;i<col-1;i++){
        // p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 + 
        //                   (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) /
        //                   (2 * (dx**2 + dy**2)) -
        //                   dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * 
        //                   b[1:-1,1:-1])    
        //[j][i] -> [1:-1, 1:-1]  
        //[j][i+1] -> [1:-1, 2:]
        //[j][i-1] -> [1:-1, 0:-2]
        //[j+1][i] -> [2:, 1:-1]
        //[j-1][i] -> [0:-2, 1:-1]         
        p[j][i] = (((pn[j][i+1] + pn[j][i-1]) * pow(dy,2) + 
                          (pn[j+1][i] + pn[j-1][i]) * pow(dx,2)) /
                          (2 * (pow(dx,2) + pow(dy,2))) -
                          pow(dx,2) * pow(dy,2) / (2 * (pow(dx,2) + pow(dy,2))) * 
                          b[j][i]);                               
      }
    }  
    for(int j=0;j<row;j++){
      //p[:, -1] = p[:, -2] # dp/dx = 0 at x = 2
      p[j][col-1] = p[j][col-2];
    }
    for(int i=0;i<col;i++){
      //p[0, :] = p[1, :]   # dp/dy = 0 at y = 0
      p[0][i] = p[1][i];
    }
    for(int j=0;j<row;j++){
      //p[:, 0] = p[:, 1]   # dp/dx = 0 at x = 0
      p[j][0] = p[j][1];
    }
    for(int i=0;i<col;i++){
      //p[-1, :] = 0        # p = 0 at y = 2
      p[row-1][i] = 0.0;
    }
  }
  // printf("print p:\n");
  // print_array(p);
  // printf("\n");  
}

void cavity_flow(int nt, vector< vector<double> > &u, vector< vector<double> > &v, double dt, double dx, double dy, vector< vector<double> > &p, double rho, double nu){
  zeros(b);
  for(int n=0;n<nt;n++){
    copy(un,u);
    copy(vn,v);
    build_up_b(b, rho, dt, u, v, dx, dy);
    pressure_poisson(p, dx, dy, b);
    int row = ny,col = nx;
    for(int j=1;j<row-1;j++){
      for(int i=1;i<col-1;i++){
        // u[1:-1, 1:-1] = (un[1:-1, 1:-1]-
        //                  un[1:-1, 1:-1] * dt / dx *
        //                 (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
        //                  vn[1:-1, 1:-1] * dt / dy *
        //                 (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
        //                  dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, 0:-2]) +
        //                  nu * (dt / dx**2 *
        //                 (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
        //                  dt / dy**2 *
        //                 (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])));        
        //[j][i] -> [1:-1, 1:-1]  
        //[j][i+1] -> [1:-1, 2:]
        //[j][i-1] -> [1:-1, 0:-2]
        //[j+1][i] -> [2:, 1:-1]
        //[j-1][i] -> [0:-2, 1:-1]            
        u[j][i] = (un[j][i]-
                         un[j][i] * dt / dx *
                        (un[j][i] - un[j][i-1]) -
                         vn[j][i] * dt / dy *
                        (un[j][i] - un[j-1][i]) -
                         dt / (2 * rho * dx) * (p[j][i+1] - p[j][i-1]) +
                         nu * (dt / pow(dx,2) *
                        (un[j][i+1] - 2 * un[j][i] + un[j][i-1]) +
                         dt / pow(dy,2) *
                        (un[j+1][i] - 2 * un[j][i] + un[j-1][i])));                                  
      }
    }  
    for(int j=1;j<row-1;j++){
      for(int i=1;i<col-1;i++){
        // v[1:-1,1:-1] = (vn[1:-1, 1:-1] -
        //                 un[1:-1, 1:-1] * dt / dx *
        //                (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
        //                 vn[1:-1, 1:-1] * dt / dy *
        //                (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
        //                 dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) +
        //                 nu * (dt / dx**2 *
        //                (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
        //                 dt / dy**2 *
        //                (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))     
        v[j][i] = (vn[j][i] -
                        un[j][i] * dt / dx *
                       (vn[j][i] - vn[j][i-1]) -
                        vn[j][i] * dt / dy *
                       (vn[j][i] - vn[j-1][i]) -
                        dt / (2 * rho * dy) * (p[j+1][i] - p[j-1][i]) +
                        nu * (dt / pow(dx,2) *
                       (vn[j][i+1] - 2 * vn[j][i] + vn[j][i-1]) +
                        dt / pow(dy,2) *
                       (vn[j+1][i] - 2 * vn[j][i] + vn[j-1][i])));        
      }
    }       
    for(int i=0;i<col;i++){
      // u[0, :]  = 0
      u[0][i] = 0;
    }  
    for(int j=0;j<row;j++){
      // u[:, 0]  = 0
      u[j][0] = 0;
    }
    for(int j=0;j<row;j++){
      // u[:, -1] = 0
      u[j][col-1] = 0;
    }
    for(int i=0;i<col;i++){
      // u[-1, :] = 1    # set velocity on cavity lid equal to 1
      u[row-1][i] = 1;
    }
    for(int i=0;i<col;i++){
      // v[0, :]  = 0
      v[0][i] = 0;
    }
    for(int i=0;i<col;i++){
      // v[-1, :] = 0
      v[row-1][i] = 0;
    }
    for(int j=0;j<row;j++){
      // v[:, 0]  = 0
      v[j][0] = 0;
    }
    for(int j=0;j<row;j++){
      // v[:, -1] = 0  
      v[j][col-1] = 0;
    }    
  }
}

int main(){
  //2-d
  for(int i=0;i<ny;i++){
    u[i].resize(nx);
    v[i].resize(nx);
    p[i].resize(nx);
    b[i].resize(nx);
    un[i].resize(nx);
    vn[i].resize(nx);
    pn[i].resize(nx);
    bn[i].resize(nx);    
  }
  zeros(u);
  zeros(v);
  zeros(p);
  zeros(b);
  // build_up_b(b, rho, dt, u, v, dx, dy);
  cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu);
  string u_json = array_2d_to_json(u,ny,nx),
         v_json = array_2d_to_json(v,ny,nx),
         p_json = array_2d_to_json(p,ny,nx);
  write_string_to_file(u_json,U_FILE);
  write_string_to_file(v_json,V_FILE);
  write_string_to_file(p_json,P_FILE);
}