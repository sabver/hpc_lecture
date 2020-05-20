#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N],r[N],fxi[N],fyi[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = r[i] = fxi[i] = fyi[i] = 0;
  }
  __m256 xvec = _mm256_load_ps(x),
         yvec = _mm256_load_ps(y),
         mvec = _mm256_load_ps(m),
         fxvec = _mm256_load_ps(fx),
         fyvec = _mm256_load_ps(fy),    
         rxvec = _mm256_setzero_ps(),
         ryvec = _mm256_setzero_ps(),
         maskvec = _mm256_setzero_ps(),
         xivec = _mm256_setzero_ps(),
         yivec = _mm256_setzero_ps(),
         recrvec = _mm256_setzero_ps(),
         fxresultv = _mm256_setzero_ps(),
         fyresultv = _mm256_setzero_ps(),
         recr3vec = _mm256_setzero_ps(),
         fxivec = _mm256_setzero_ps(),
         fyivec = _mm256_setzero_ps(),
         zerovec = _mm256_setzero_ps(),
         minus1vec = _mm256_set1_ps(-1);
  for(int i=0; i<N; i++) {         
    xivec = _mm256_set1_ps(x[i]);  
    yivec = _mm256_set1_ps(y[i]);  
    
    //double rx = x[i] - x[j];  
    rxvec = _mm256_sub_ps(xivec,xvec);       
    //double ry = y[i] - y[j];
    ryvec = _mm256_sub_ps(yivec,yvec); 
    //1/sqrt(rx * rx + ry * ry); 
    recrvec = _mm256_rsqrt_ps(_mm256_add_ps(_mm256_mul_ps(rxvec,rxvec),_mm256_mul_ps(ryvec,ryvec)));            
    //mask the recrvec  
    maskvec = _mm256_cmp_ps(recrvec, _mm256_set1_ps(INFINITY), _CMP_NEQ_OQ);               
    recrvec = _mm256_blendv_ps(zerovec,recrvec,maskvec);            
    //1/r * r * r
    recr3vec = _mm256_mul_ps(recrvec,_mm256_mul_ps(recrvec,recrvec));     
    //rx * m[j] / (r * r * r)   
    fxresultv = _mm256_mul_ps(_mm256_mul_ps(rxvec,mvec),recr3vec); 
    //ry * m[j] / (r * r * r)  
    fyresultv = _mm256_mul_ps(_mm256_mul_ps(ryvec,mvec),recr3vec); 
    //reduction fxresultv
    fxivec = _mm256_permute2f128_ps(fxresultv,fxresultv,1);
    fxivec = _mm256_add_ps(fxivec,fxresultv);
    fxivec = _mm256_hadd_ps(fxivec,fxivec);
    fxivec = _mm256_hadd_ps(fxivec,fxivec);  
    // multipy -1
    fxivec = _mm256_mul_ps(fxivec,minus1vec);
    //reduction fyresultv      
    fyivec = _mm256_permute2f128_ps(fyresultv,fyresultv,1);
    fyivec = _mm256_add_ps(fyivec,fyresultv);
    fyivec = _mm256_hadd_ps(fyivec,fyivec);
    fyivec = _mm256_hadd_ps(fyivec,fyivec); 
    // multipy -1
    fyivec = _mm256_mul_ps(fyivec,minus1vec);      
    _mm256_store_ps(fxi,fxivec);  
    _mm256_store_ps(fyi,fyivec);       
    
    printf("%d %g %g\n",i,fxi[i],fyi[i]);      
  }
}
