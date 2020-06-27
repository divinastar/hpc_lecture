#include <cstdio>
#include <vector>
#include <algorithm>
#include <iterator>

float build_up_b(int rho, float dt, float dx, float dy, float u, float v){
   float b[][];
   for(int i=0;i>n;i++){
      
   }
   #Periodic BC Pressure @x = 2
   
   #Periodic BC Pressure @x = 0
   
   return b;
}

float pressure_poisson_periodic(float p, float nx, float ny, int nit){
   for(int q=0; q<nit; q++){
   }  
   return p;
}

int main() {
   //Variable Declarations
   int nx = 41;
   int ny = 41;
   int nt = 10;
   int nit = 50;
   int c = 1;
   float dx = 2/(nx - 1);
   float dy = 2/(ny - 1);
   float x[nx+1];
   float y[ny+1];
   float X[nx+1][ny+1];
   float Y[nx+1][ny+1];
   

   //Physical Variables
   int rho = 1;
   float nu = .1;
   int F = 1;
   float dt = .01;
   
   //Initial Conditions
   float u[ny][nx];
   float un[ny][nx];
   float v[ny][nx];
   float vn[ny][nx];
   float p[ny][nx];
   float pn[ny][nx];
   float b[ny][nx];
   
   
   for(int i=0;i<=nx;i++){
      x[i] =  (2-0)*i/nx;
   }
   
   for(int i=0;i<=ny;i++){
      y[i] =  (2-0)*i/nx;
   }
   
   for(int i=0;i<=nx;i++){
      for(j=0;j<=ny;j++){
         X[i][j] = x[i];
         Y[i][j] = y[j];
         u[j][i] = 0;
         un[j][i] = 0;
         v[j][i] = 0;
         vn[j][i] = 0;
         p[j][i] = 1;
         pn[j][i] = 1;
         b[j][i] = 0;
      }
   }
   
   int udiff = 1;
   int stepcount = 0;
   
   while(udiff>.001){
      un = std::copy(u);
      vn = std::copy(v);
      
      b = build_up_b(rho, dt, dx, dy, u, v);
      p = pressure_poisson_periodic(p, dx, dy);
      
      
      
      udiff = (sum(u) - sum(un))/ sum(u)
      stepcount += 1
         
   }
   
   
   
   
      
}
