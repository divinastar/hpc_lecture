#include <cstdio>
#include <iostream>
#include <cmath>

const int nx = 41;
const int ny = 41;
const int nt = 10;
const int nit = 50;
const int c = 1;

float[][] build_up_b(int rho, float dt, float dx, float dy, float u[][], float v[][]){
   float b[ny][nx];
   for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){
         b[j][i] = 0;
      }
   }
   
   for(int i=1;i<nx-1;i++){
      for(int j=1;j<ny-1;j++){
         b[j][i] = (rho*(1/dt* ((u[j][i+1] - u[j][i-1])/(2*dx) 
                              + (v[j+1][i] - v[j-1][i])/(2*dy))
                              - pow((u[j][i+1] - u[j][i-1])/(2*dx),2)
                              - 2*((u[j+1][i] - u[j-1][i])/(2*dy) *
                                    (v[j][i+1]-v[j][i-1])/(2*dx)) -
                                    pow((v[j+1][i] - v[j-1][i])/(2*dy),2)));
      
      }
   }
   //Periodic BC Pressure @x = 2
   for(int j=1;j<ny;j++){
      b[j][nx-1] = (rho*(1/dt*((u[j][0] - u[j][nx-2])/(2*dx)
                              +(v[j+1][nx-1] - v[j-1][nx-1])/(2*dy))
                              -pow((u[j][0] - u[j][nx-2])/(2*dx),2)
                               - 2*((u[j+1][nx-1]-u[j-1][nx-1])/(2*dy) *
                                    (v[j][0] - v[j][nx-2])/(2*dx)) -
                                    pow((v[j+1][nx-1] - v[j-1][nx-1])/(2*dy),2)));
   }
                                   
                                   
   //Periodic BC Pressure @x = 0                               
    for(int j=1;j<ny;j++){
      b[j][0] = (rho*(1/dt*((u[j][1] - u[j][nx-1])/(2*dx)
                              +(v[j+1][0] - v[j-1][0])   /(2*dy))
                              - pow((u[j][1] - u[j][nx-1])/(2*dx),2)
                               - 2*((u[j+1][0]-u[j-1][0])/(2*dy) *
                                    (v[j][1] - v[j][nx-1])/(2*dx)) -
                                    pow((v[j+1][0] - v[j-1][0])/(2*dy),2)));
   }                                  
   return b;
}

float[][] pressure_poisson_periodic(float p[][], float dx, float dy){
   float pn[ny][nx];
   
   for(int q=0; q<nit; q++){
      for(int i=0;i<=nx;i++){
         for(int j=0;j<=ny;j++){
             pn[j][i] = p[j][i];
         }
      }
      
      //Periodic BC Pressure @ x = 2
      //Periodic BC Pressure @ x = 0
      //Wall boundary conditions, pressure
       //dp/dy = 0 at y = 2
       //dp/dy = 0 at y = 0
   }
   return p;
}

int main() {
   //Variable Declarations
   float dx = 2/(nx - 1);
   float dy = 2/(ny - 1);
   float x[nx+1];
   float y[ny+1];
   float X[nx+1][ny+1];
   float Y[nx+1][ny+1];
   
   //Physical Variables
   const int rho = 1;
   const float nu = .1;
   const int F = 1;
   const float dt = .01;
   
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
   
  for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){
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
   float sumu = 0;
   float sumun = 0;
   
   while(udiff>.001){
      for(int i=0; i<nx; i++){
         for(int j=0; j<ny; j++){
            un[j][i] = u[j][i];
            vn[j][i] = v[j][i];
         }    
      }
      
      b = build_up_b(rho, dt, dx, dy, u, v);
      p = pressure_poisson_periodic(p, dx, dy);
      
      //Periodic BC u @ x = 2
      //Periodic BC u @ x = 0
      //Periodic BC v @ x = 2
      //Periodic BC v @ x = 0
      
      sumu = 0;
      sumun = 0;
      for(int i=0;i<nx;i++){
         for(int j=0; j<ny;j++){
            sumu += u[j][i];
            sumun += un[j][i];
         }
      }
      
      udiff = (sumu - sumun)/ sumu ;
      stepcount += 1;     
   }
}
