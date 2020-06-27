#include <cstdio>
#include <vector>
#include <iostream>
#include <cmath>

const int nx = 41;
const int ny = 41;
const int nt = 10;
const int nit = 50;
const int c = 1;

float** build_up_b(int rho, float dt, float dx, float dy, float **u , float **v){
   float **b;
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
                              - 2*((u[j+1][i] - u[j-1][i])/(2*dy)*
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

float** pressure_poisson_periodic(float **p, float **b, float dx, float dy){
   float **pn;
   
   for(int q=0; q<nit; q++){
      for(int i=0;i<nx;i++){
         for(int j=0;j<ny;j++){
             pn[j][i] = p[j][i];
         }
      }

      for(int i=1; i<nx-1;i++){
         for(int j=1;j<ny-1;j++){
            p[j][i]=(((pn[j][i+1]+pn[j][i-1])*pow(dy,2)+
                      (pn[j+1][i]+pn[j-1][i])*pow(dx,2))/
                      (2*(pow(dx,2)+pow(dy,2)))-
                      pow(dx,2)*pow(dy,2)/(2*(pow(dx,2)+pow(dy,2)))*b[j][i]);
         }
      }
      
      //Periodic BC Pressure @ x = 2
      for(int j=1;j<ny-1;j++){
         p[j][nx-1]=(((pn[j][0]+pn[j][nx-2])*pow(dy,2)+
                     (pn[j+1][nx-1]+pn[j-1][nx-1])*pow(dx,2))/
                     (2*(pow(dx,2)+pow(dy,2)))-
                     pow(dx,2)*pow(dy,2)/(2*(pow(dx,2)+pow(dy,2)))*b[j][nx-1]);
      }
      //Periodic BC Pressure @ x = 0
      for(int j=1;j<ny-1;j++){
         p[j][0]=(((pn[j][1]+pn[j][nx-1])*pow(dy,2)+
                     (pn[j+1][0]+pn[j-1][0])*pow(dx,2))/
                     (2*(pow(dx,2)+pow(dy,2)))-
                     pow(dx,2)*pow(dy,2)/(2*(pow(dx,2)+pow(dy,2)))*b[j][0]);
      }
      //Wall boundary conditions, pressure
      p[nx-1] = p[nx-2]; //dp/dy = 0 at y = 2
      p[0] = p[1]; //dp/dy = 0 at y = 0
   }
   return p;
}

int main() {
   //Variable Declarations
   float dx = 2/(nx - 1);
   float dy = 2/(ny - 1);
   float *x;
   float *y;
   float **X;
   float **Y;
   
   //Physical Variables
   const int rho = 1;
   const float nu = .1;
   const int F = 1;
   const float dt = .01;
   
   //Initial Conditions
   float **u;
   float **un;
   float **v;
   float **vn;
   float **p;
   float **pn;
   float **b;
   
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
      p = pressure_poisson_periodic(p, b, dx, dy);
      
      for(int i=1;i<nx;i++){
         for(int j=1;j<ny;j++){
            u[j][i] = (un[j][i] -
                       un[j][i] * dt/dx *
                      (un[j][i] - un[j][i-1]) -
                       vn[j][i] * dt/dy *
                      (un[j][i] - un[j-1][i]) -
                       dt/(2*rho*dx) *
                      (p[j][i+1] - p[j][i-1]) +
                       nu * (dt/pow(dx,2)*
                      (un[j][i+1] - 2*un[j][i] + un[j][i-1]) +
                       dt/pow(dy,2) *
                      (un[j+1][i] - 2*un[j][i] + un[j-1][i])) +
                       F* dt); 

            v[j][i] = (vn[j][i] -
                       un[j][i] * dt/dx *
                      (vn[j][i] - vn[j][i-1]) -
                       vn[j][i] * dt/dy *
                      (vn[j][i] - vn[j-1][i]) -
                       dt/(2*rho*dy) *
                      (p[j+1][i] - p[j-1][i]) +
                       nu * (dt/pow(dx,2)*
                      (vn[j][i+1] - 2*vn[j][i] + vn[j][i-1]) +
                       dt/pow(dy,2) *
                      (vn[j+1][i] - 2*vn[j][i] + vn[j-1][i]))); 
         }
      }

      
      for(int j=1;j<ny;j++){
      //Periodic BC u @ x = 2
            u[j][nx-1] = (un[j][nx-1] -  
                       un[j][nx-1] * dt/dx *
                      (un[j][nx-1] - un[j][nx-2]) -
                       vn[j][nx-1] * dt/dy *
                      (un[j][nx-1] - un[j-1][nx-1]) -
                       dt/(2*rho*dx) *
                      (p[j][0] - p[j][nx-2]) +
                       nu * (dt/pow(dx,2)*
                      (un[j][0] - 2*un[j][nx-1] + un[j][nx-2]) +
                       dt/pow(dy,2) *
                      (un[j+1][nx-1] - 2*un[j][nx-1] + un[j-1][nx-1])) +
                       F* dt); 
      //Periodic BC u @ x = 0
            u[j][0] = (un[j][0] -  
                       un[j][0] * dt/dx *
                      (un[j][0] - un[j][nx-1]) -
                       vn[j][0] * dt/dy *
                      (un[j][0] - un[j-1][0]) -
                       dt/(2*rho*dx) *
                      (p[j][1] - p[j][nx-1]) +
                       nu * (dt/pow(dx,2)*
                      (un[j][1] - 2*un[j][0] + un[j][nx-1]) +
                       dt/pow(dy,2) *
                      (un[j+1][0] - 2*un[j][0] + un[j-1][0])) +
                       F* dt); 

      //Periodic BC v @ x = 2
            v[j][nx-1] = (vn[j][nx-1] -
                       un[j][nx-1] * dt/dx *
                      (vn[j][nx-1] - vn[j][nx-2]) -
                       vn[j][nx-1] * dt/dy *
                      (vn[j][nx-1] - vn[j-1][nx-1]) -
                       dt/(2*rho*dy) *
                      (p[j+1][nx-1] - p[j-1][nx-1]) +
                       nu * (dt/pow(dx,2)*
                      (vn[j][0] - 2*vn[j][nx-1] + vn[j][nx-2]) +
                       dt/pow(dy,2) *
                      (vn[j+1][nx-1] - 2*vn[j][nx-1] + vn[j-1][nx-1]))); 
      //Periodic BC v @ x = 0
            v[j][0] = (vn[j][0] -
                       un[j][0] * dt/dx *
                      (vn[j][0] - vn[j][nx-1]) -
                       vn[j][0] * dt/dy *
                      (vn[j][0] - vn[j-1][0]) -
                       dt/(2*rho*dy) *
                      (p[j+1][0] - p[j-1][0]) +
                       nu * (dt/pow(dx,2)*
                      (vn[j][1] - 2*vn[j][0] + vn[j][nx-1]) +
                       dt/pow(dy,2) *
                      (vn[j+1][0] - 2*vn[j][0] + vn[j-1][0]))); 
      }

      
      //Wall BC: u,v = 0 @ y = 0,2
      for(int j=0; j<nx; j++){
         u[0][j] = 0;
         u[ny-1][j] = 0;
         v[0][j] = 0;
         v[ny-1][j] = 0;
      }
      
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

   std::cout<< stepcount <<std::endl;

}
