#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

int nx = 41;
float dx = 2.0 / (nx - 1);
int nt = 20;    //number of timesteps we want to calculate
float dt = .025;  //amount of time each timestep covers (delta t)

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

void write_string_to_file(string str ){
  ofstream outfile;
  outfile.open("/home/9/19M38171/t3workspace/hpc-lecture/final_report/json.txt");
  outfile << str << endl;
  outfile.close();
}


vector<float> linspace(int start,int end){
  // printf("start:%d,end:%d,nx:%d\n",start,end,nx);
  vector<float> result(nx);
  float interval = (float)( end - start ) / ( nx - 1 );
  // printf("interval:%f\n",interval);
  result[0] = start;
  result[nx-1] = end;
  for(int i=1;i<nx-1;i++){
    result[i] = result[i-1] + interval;
  }
  return result;
}

vector<float> ones(){
  vector<float> result(nx);
  for(int i=0;i<nx;i++){
    result[i] = 1.0;
  }
  return result;
}

/**
 * u[int(.5 / dx) : int(1 / dx + 1)] = 2
 */
void set_value(float dx,vector<float> &arr,float start_point,float end_point, float value){
//   printf("%f\n",start_point/dx);
//   printf("%f\n",end_point/dx);    
  float start = start_point/dx,end = end_point/dx + 1;
//   printf("start:%f\nend:%f\n",start,end);
  for(int i=(int)start;i<=(int)end;i++){
    arr[i] = value;
  }
}

void convection( vector<float> &u,vector<float> &u_old){
  //u_old = u.copy()
  u_old.assign(u.begin(), u.end());
  for(int i = 1;i<nx;i++){
    u[i] = u_old[i] - u_old[i] * dt / dx * (u_old[i] - u_old[i-1]);
  }
}


int main() {
  vector<float> x = linspace(0,2),u = ones();
  set_value( dx, u, 0.5, 1.0, 2 );
  // for(int i=0;i<nx;i++){
  //   printf("%f,",u[i]);
  // }  
  // printf("\n");  
  vector<float> u_old(nx);
  vector<string> array_str(nt);
  for(int i=0;i<nt;i++){
    convection(u,u_old);
    //save u's data    
    string result = array_to_string(u);
    array_str[i] = result;
    // cout << "result:" << add_s(result) << endl;      
  }
  string json = to_json(array_str,nt);
  // cout << "json:" << json << endl;  
  // const string file_string = "/home/9/19M38171/t3workspace/hpc-lecture/final_report/json.txt";
  write_string_to_file(json);    
}