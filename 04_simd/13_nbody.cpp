#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }

  __m256 xvec = _mm256_load_ps(x);
  __m256 yvec = _mm256_load_ps(y);
  __m256 mvec = _mm256_load_ps(m);
  __m256 fxvec = _mm256_setzero_ps();
  __m256 fyvec = _mm256_setzero_ps();
  __m256 xtemp = _mm256_load_ps(x);
  __m256 ytemp = _mm256_load_ps(y);
  
//Permute constant
  int idx[N] = {1,2,3,4,5,6,7,0};
  __m256i ivec = _mm256_load_si256((__m256i*)idx);
	

  for(int i=0; i<N-1; i++) {
//Shifting xtemp and ytemp
      xtemp = _mm256_permutevar8x32_ps(xtemp,ivec);
      ytemp = _mm256_permutevar8x32_ps(ytemp,ivec);    

//Calculate rx, ry by subtracting xvec with xtemp
      __m256 rx = _mm256_sub_ps(xvec,xtemp);
      __m256 ry = _mm256_sub_ps(yvec,ytemp);

//Calculate rx2, ry2, rx2_ry2, and so on...
      __m256 rx2 = _mm256_mul_ps(rx,rx);
      __m256 ry2 = _mm256_mul_ps(ry,ry);
      __m256 rx2_ry2 = _mm256_add_ps(rx2,ry2);
//      __m256 r = _mm256_sqrt_ps(rx2_ry2);
      __m256 rsqrt = _mm256_rsqrt_ps(rx2_ry2);
      __m256 fxd = _mm256_mul_ps(rsqrt,_mm256_mul_ps(rsqrt, _mm256_mul_ps(rsqrt, _mm256_mul_ps(mvec,rx))));
      __m256 fyd = _mm256_mul_ps(rsqrt,_mm256_mul_ps(rsqrt, _mm256_mul_ps(rsqrt, _mm256_mul_ps(mvec,ry))));
//      __m256 fxd = _mm256_div_ps(_mm256_div_ps(_mm256_div_ps(_mm256_mul_ps(rx,mvec),r),r),r);
//      __m256 fyd = _mm256_div_ps(_mm256_div_ps(_mm256_div_ps(_mm256_mul_ps(ry,mvec),r),r),r);
	fxvec = _mm256_sub_ps(fxvec,fxd);
        fyvec = _mm256_sub_ps(fyvec,fyd);
  }

  _mm256_store_ps(fx, fxvec);
  _mm256_store_ps(fy, fyvec);
  
  for(int i=0; i<N; i++){
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}

