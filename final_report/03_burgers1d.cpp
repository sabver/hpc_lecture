#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>

#define PI acos(-1)

using namespace std;

int nx = 101;
int nt = 100;    
float dx = 2 * PI / (nx - 1);
float nu = 0.07; 
float sigma = 0.2; 
float dt = dx * nu;

template < typename Type > std::string to_str (const Type & t)
{
  std::ostringstream os;
  os << t;
  return os.str ();
}

void print_u(vector<double> u){
  for(int i=0;i<nx;i++){
    cout << to_str(u[i]) << ",";      
  }
  cout << endl;      
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

string array_to_string( vector<double> vector ){
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


vector<double> linspace(double start,double end){
  vector<double> result(nx);
  double interval = (float)( end - start ) / ( nx - 1 );
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

double ufunc( float t,float x,float nu ){
  double result = -2*nu*(-(-8*t + 2*x)*exp(-pow((-4*t + x),2)/(4*nu*(t + 1)))/(4*nu*(t + 1)) - (-8*t + 2*x - 4*PI)*exp(-pow((-4*t + x - 2*PI),2)/(4*nu*(t + 1)))/(4*nu*(t + 1)))/(exp(-pow((-4*t + x - 2*PI),2)/(4*nu*(t + 1))) + exp(-pow((-4*t + x),2)/(4*nu*(t + 1)))) + 4;
  return result;
}

/**
 * u[int(.5 / dx) : int(1 / dx + 1)] = 2
 */
void set_value(float dx,vector<float> &arr,float start_point,float end_point, float value){  
  float start = start_point/dx,end = end_point/dx + 1;
  for(int i=(int)start;i<=(int)end;i++){
    arr[i] = value;
  }
}

void burgers( vector<double> &u,vector<double> &u_old){
  //u_old = u.copy()
  u_old.assign(u.begin(), u.end());
  for(int i = 1;i<nx-1;i++){    
    u[i] = u_old[i] - u_old[i] * dt / dx *(u_old[i] - u_old[i-1]) + nu * dt / pow(dx,2) * (u_old[i+1] - 2 * u_old[i] + u_old[i-1]);    
  }
  u[0] = u_old[0] - u_old[0] * dt / dx * (u_old[0] - u_old[nx-2]) + nu * dt / pow(dx,2) * (u_old[1] - 2 * u_old[0] + u_old[nx-2]);
  //u[-1] = u[0]
  u[nx-1] = u[0];
}

int main() {
  // double result = ufunc(1, 4, 3);
  // printf("result:%f\n",result);
  // cout << result << endl; 
  vector<double> x = linspace(0.0,2*PI);
  // printf("%f\n",x[nx-1]);
  // set_value( dx, u, 0.5, 1.0, 2 );
  double t = 0.0;
  //u = numpy.asarray([ufunc(t, x0, nu) for x0 in x])
  vector<double> u(nx);
  for(int i=0;i<nx;i++){
    u[i] = ufunc(t, x[i], nu);
  }
  // print_u(u);
  
  vector<double> u_old(nx);
  vector<string> array_str(nt);
  for(int i=0;i<nt;i++){
    burgers(u,u_old);
    //save u's data    
    string result = array_to_string(u);
    array_str[i] = result;
  }
  string json = to_json(array_str,nt);
  write_string_to_file(json);    
}