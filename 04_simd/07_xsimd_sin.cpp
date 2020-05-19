#include <cstdio>
#include <cmath>
#include "xsimd/xsimd.hpp"

int main() {
  const int N = 8;
  xsimd::batch<float, N> a;
  for(int i=0; i<N; i++)
    a[i] = i * M_PI / (N - 1);
  float b = xsimd::hadd(a);
  prinf("%g\n",b);
  //for(int i=0; i<N; i++)
    //printf("%g %g\n",a[i],b[i]);
}
