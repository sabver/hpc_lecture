#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;
int ny = 41; //row
int nx = 41; //col
int nt = 50;
int c = 1;
float dx = 2.0 / (nx - 1);
float dy = 2.0 / (ny - 1);
float dt = dx * 0.2;

vector< vector<double> > u(ny); //ny * nx
vector< vector<double> > v(ny);

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


void ones(vector< vector<double> > &vector){
  for(int y = 0;y<ny;y++){
    for(int x = 0;x<nx;x++){
      vector[y][x] = 1;
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

void write_string_to_file(string str,int isU){
  ofstream outfile;
  if( isU == 1 ){
    outfile.open("/home/9/19M38171/t3workspace/hpc-lecture/final_report/u_json.txt");  
  }else{
    outfile.open("/home/9/19M38171/t3workspace/hpc-lecture/final_report/v_json.txt");
  }  
  outfile << str << endl;
  outfile.close();
}

/**
 * u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2
 */
void set_value(vector< vector<double> > &vector,float dy_start_point,float dy_end_point,float dx_start_point,float dx_end_point, double value){  
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

void set_value2(vector< vector<double> > &vector , int row_start,int row_end,int col_start,int col_end,double value ){
  for(int i=row_start;i<=row_end;i++){
    for(int j=col_start;j<=col_end;j++){
      vector[i][j] = value;
    }
  }
}

void convection( vector< vector<double> > &u,vector< vector<double> > &u_old,vector< vector<double> > &v,vector< vector<double> > &v_old){
//   printf("convection\n");
  //u_old = u.copy()
  u_old.assign(u.begin(), u.end());
  v_old.assign(v.begin(), v.end());
  int row = ny,col = nx;
  for(int j=1;j<row;j++){
    for(int i=1;i<col;i++){
      // u[1:, 1:] = (u_old[1:, 1:] - 
      //              (u_old[1:, 1:] * dt / dx * (u_old[1:, 1:] - u_old[1:, :-1])) -
      //               v_old[1:, 1:] * dt / dy * (u_old[1:, 1:] - u_old[:-1, 1:]))      
      //[j][i] -> [1:, 1:] 
      //[j][(i+col-1)%col] -> [1:, :-1]
      //[(j-1+row)%row][i] -> [:-1, 1:]
      u[j][i] = (u_old[j][i] - 
                 (u_old[j][i] * dt / dx * (u_old[j][i] - u_old[j][(i-1+col)%col])) -
                  v_old[j][i] * dt / dy * (u_old[j][i] - u_old[(j-1+row)%row][i]));

      // v[1:, 1:] = (v_old[1:, 1:] -
      //              (u_old[1:, 1:] * dt / dx * (v_old[1:, 1:] - v_old[1:, :-1])) -
      //              v_old[1:, 1:] * dt / dy * (v_old[1:, 1:] - v_old[:-1, 1:]))
      v[j][i] = (v_old[j][i] -
                 (u_old[j][i] * dt / dx * (v_old[j][i] - v_old[j][(i+col-1)%col])) -
                 v_old[j][i] * dt / dy * (v_old[j][i] - v_old[(j-1+row)%row][i]));
    }
  }
  //u[0, :] = 1
  set_value2(u,0,0,0,col-1,1.0);
  //u[-1, :] = 1
  set_value2(u,row-1,row-1,0,col-1,1.0);
  //u[:, 0] = 1
  set_value2(u,0,row-1,0,0,1.0);
  //u[:, -1] = 1
  set_value2(u,0,row-1,col-1,col-1,1.0);  

  //v[0, :] = 1
  set_value2(v,0,0,0,col-1,1.0);
  //v[-1, :] = 1
  set_value2(v,row-1,row-1,0,col-1,1.0);
  //v[:, 0] = 1
  set_value2(v,0,row-1,0,0,1.0);
  //v[:, -1] = 1
  set_value2(v,0,row-1,col-1,col-1,1.0);   
}

int main(){

  //2-d
  for(int i=0;i<ny;i++){
    u[i].resize(nx);
    v[i].resize(nx);
  }
  ones(u);
  ones(v);
  //u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2
  set_value(u,0.5,1.0,0.5,1.0,2.0);
  //v[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2
  set_value(v,0.5,1.0,0.5,1.0,2.0);
  vector< vector<double> > u_old(ny);
  vector< vector<double> > v_old(ny);
  vector<string> array_u_str(nt);
  vector<string> array_v_str(nt);
  for(int i=0;i<nt;i++){
    convection(u,u_old,v,v_old);
    //save u's data    
    // string result = array_2d_to_json(u,ny,nx);
    // array_str[i] = result;
    array_u_str[i] = array_2d_to_json(u,ny,nx);
    array_v_str[i] = array_2d_to_json(v,ny,nx);
  }  
  // print_array();
  // string json = to_json(array_str,nt);
  // write_string_to_file(json); 
  write_string_to_file(to_json(array_u_str,nt),1);
  write_string_to_file(to_json(array_u_str,nt),0);
}