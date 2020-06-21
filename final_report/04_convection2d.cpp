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
float dx = 2 / (nx - 1);
float dy = 2 / (ny - 1);
float dt = dx * 0.2;

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

string array_to_string( vector<float> vector ){
  string result = "[";
  for(int i=0;i<nx;i++){
    result += to_str(vector[i]);
    if( i!=nx-1 ){
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

int main(){

  //2-d
  for(int i=0;i<ny;i++){
    u[i].resize(nx);
  }
  ones();
  set_value(0.5,1.0,0.5,1.0,2.0);
  print_array();

}