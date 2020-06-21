#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>

using namespace std;
int ny = 41; //row
int nx = 41; //col
int nt = 50;
double nu = 0.01;
float dx = 2.0 / (nx - 1);
float dy = 2.0 / (ny - 1);
float dt = dx * dy / nu * 0.25;

vector< vector<double> > u(ny); //ny * nx

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

// string all_to_json(vector< vector<double> > vector,int row,int col){
//   string result = "[";
//   for(int i=0;i<row;i++){
//     result += to_json(vector[i],col);
//     if( i!=row-1 ){
//       result += ",";
//     }    
//   }
//   return result+"]";
// }

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


void ones(){
  for(int y = 0;y<ny;y++){
    for(int x = 0;x<nx;x++){
      u[y][x] = 1;
    }
  }
}

void print_array(){
    for(int i=0;i<u.size();i++){
        for(int j=0;j<u[i].size();j++){
            printf("%f,",u[i][j]);
        }
        printf("\n");
    }
}

void write_string_to_file(string str ){
  ofstream outfile;
  outfile.open("/home/9/19M38171/t3workspace/hpc-lecture/final_report/json.txt");
  outfile << str << endl;
  outfile.close();
}

/**
 * u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2
 */
void set_value(float dy_start_point,float dy_end_point,float dx_start_point,float dx_end_point, double value){  
  float dy_start = dy_start_point/dy,
        dy_end = dy_end_point/dy + 1,
        dx_start = dx_start_point/dx,        
        dx_end = dx_end_point/dx + 1;

  for(int i=(int)dy_start;i<=(int)dy_end;i++){
    for(int j=(int)dx_start;j<=(int)dx_end;j++){
      u[i][j] = value;
    }
  }
}

void set_value2( int row_start,int row_end,int col_start,int col_end,double value ){
  for(int i=row_start;i<=row_end;i++){
    for(int j=col_start;j<=col_end;j++){
      u[i][j] = value;
    }
  }
}

void diffusion( vector< vector<double> > &u,vector< vector<double> > &u_old ){
  //u_old = u.copy()
  u_old.assign(u.begin(), u.end());
  int row = ny,col = nx;  
  for(int j=1;j<row-1;j++){
    for(int i=1;i<col-1;i++){
      // u[1:-1, 1:-1] = (u_old[1:-1,1:-1] + 
      //               nu * dt / dx**2 * (u_old[1:-1, 2:] - 2 * u_old[1:-1, 1:-1] + u_old[1:-1, 0:-2]) +
      //               nu * dt / dy**2 * (u_old[2:,1: -1] - 2 * u_old[1:-1, 1:-1] + u_old[0:-2, 1:-1]))
      // [j][i] -> [1:-1, 1:-1]    
      // [j][i+1] -> [1:-1, 2:]          
      // [j][i-1] -> [1:-1, 0:-2]
      // [j+1][i] -> [2:,1:-1]
      // [j-1][i] -> [0:-2, 1:-1]
      u[j][i] = (u_old[j][i] + 
                    nu * dt / pow(dx,2) * (u_old[j][i+1] - 2 * u_old[j][i] + u_old[j][i-1]) +
                    nu * dt / pow(dy,2) * (u_old[j+1][i] - 2 * u_old[j][i] + u_old[j-1][i]));
      //u[0, :] = 1
      set_value2(0,0,0,col-1,1.0);
      //u[-1, :] = 1
      set_value2(row-1,row-1,0,col-1,1.0);
      //u[:, 0] = 1
      set_value2(0,row-1,0,0,1.0);
      //u[:, -1] = 1
      set_value2(0,row-1,col-1,col-1,1.0);
    }
  }  
}

int main(){

  //2-d
  for(int i=0;i<ny;i++){
    u[i].resize(nx);
  }
  ones();
  set_value(0.5,1.0,0.5,1.0,2.0);
  vector< vector<double> > u_old(ny);
  vector<string> array_str(nt);
  for(int i=0;i<nt;i++){
    diffusion(u,u_old);
    //save u's data    
    string result = array_2d_to_json(u,ny,nx);
    array_str[i] = result;
  }  
  // print_array();
  string json = to_json(array_str,nt);
  write_string_to_file(json); 
}