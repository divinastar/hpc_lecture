#include <cstdio>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;

const int nx = 41;
const int ny = 41;
//const int nt = 10;
const int nit = 50;
//const int c = 1;

__global__ void build_up_b(float *b, int rho, float dt, float dx, float dy, float *u , float *v){
   //m = j * nx + i
   int m = blockIdx.x * blockDim.x + threadIdx.x;
   int bId = blockIdx.x;
   int tId = threadIdx.x;

   b[m] = 0;
   
   if(bId != 0 && bId != (ny-1) && tId != 0 && tId != (nx-1)){
      b[m] = (rho*(1/dt* ((u[m+1]-u[m-1])/(2*dx)
                              + (v[m+nx] - v[m-nx])/(2*dy))
                              - pow((u[m+1] - u[m-1])/(2*dx),2)
                              - 2*((u[m+nx] - u[m-nx])/(2*dy)*
                                    (v[m+1]-v[m-1])/(2*dx)) -
                                    pow((v[m+nx] - v[m-nx])/(2*dy),2)));
    
   }else if(bId != 0 && bId != (ny-1) && tId == (nx-1)){
      //Periodic BC Pressure @x = 2
      b[m] = (rho*(1/dt*((u[m-nx+1] - u[m-1])/(2*dx)
                              +(v[m+nx] - v[m-nx])/(2*dy))
                              -pow((u[m-nx+1] - u[m-1])/(2*dx),2)
                               - 2*((u[m+nx]-u[m-nx])/(2*dy) *
                                    (v[m-nx+1] - v[m-1])/(2*dx)) -
                                    pow((v[m+nx] - v[m-nx])/(2*dy),2)));
   }else if(bId != 0 && bId != (ny-1) && tId == (0)){
      //Periodic BC Pressure @x = 0                               
         b[m] = (rho*(1/dt*((u[m+1] - u[m+nx-1])/(2*dx)
                                 +(v[m+nx] - v[m-nx])   /(2*dy))
                                 - pow((u[m+1] - u[m+nx-1])/(2*dx),2)
                                  - 2*((u[m+nx]-u[m-nx])/(2*dy) *
                                       (v[m+1] - v[m+nx-1])/(2*dx)) -
                                       pow((v[m+nx] - v[m-nx])/(2*dy),2)));
       
   }
}

__global__ void pressure_poisson_periodic(float *p, float *pn, float *b, float dx, float dy){
   //m = j * nx + i
   int m = blockIdx.x * blockDim.x + threadIdx.x;
   int bId = blockIdx.x;
   int tId = threadIdx.x;

   for(int q=0; q<nit; q++){
      __syncthreads();
      pn[m] = p[m];
      __syncthreads();

      if(bId != 0 && bId != (ny-1) && tId != 0 && tId != (nx-1)){
         p[m]=(((pn[m+1]+pn[m-1])*pow(dy,2)+
                      (pn[m+nx]+pn[m-nx])*pow(dx,2))/
                      (2*(pow(dx,2)+pow(dy,2)))-
                      pow(dx,2)*pow(dy,2)/(2*(pow(dx,2)+pow(dy,2)))*b[m]);
      }else if(bId != 0 && bId != (ny-1) && tId == (nx-1)){
         //Periodic BC Pressure @ x = 2
         p[m]=(((pn[m-nx+1]+pn[m-1])*pow(dy,2)+
                     (pn[m+nx]+pn[m-nx])*pow(dx,2))/
                     (2*(pow(dx,2)+pow(dy,2)))-
                     pow(dx,2)*pow(dy,2)/(2*(pow(dx,2)+pow(dy,2)))*b[m]);
      }else if(bId != 0 && bId != (ny-1) && tId == (0)){
         //Periodic BC Pressure @ x = 0
         p[m]=(((pn[m+1]+pn[m+nx-1])*pow(dy,2)+
                     (pn[m+nx]+pn[m-nx])*pow(dx,2))/
                     (2*(pow(dx,2)+pow(dy,2)))-
                     pow(dx,2)*pow(dy,2)/(2*(pow(dx,2)+pow(dy,2)))*b[m]);
      }else if(bId == 0){
         p[m] = p[m+nx];
      }else if(bId == (ny-1)){
         p[m] = p[m-nx];
      }  
   }
}

__global__ void updated_u_v(float *u, float *v, float *un, float *vn, float *p, float dx, float dy, float dt, float rho, float nu, float F){
   //m = j * nx + i
   int m = blockIdx.x * blockDim.x + threadIdx.x;
   int bId = blockIdx.x;
   int tId = threadIdx.x;
   if(bId != 0 && bId != (ny-1) && tId != 0 && tId != (nx-1)){
      u[m] = (un[m] -
         un[m] * dt/dx *
        (un[m] - un[m-1]) -
         vn[m] * dt/dy *
        (un[m] - un[m-nx]) -
         dt/(2*rho*dx) *
        (p[m+1] - p[m-1]) +
         nu * (dt/pow(dx,2)*
        (un[m+1] - 2*un[m] + un[m-1]) +
         dt/pow(dy,2) *
        (un[m+nx] - 2*un[m] + un[m-nx])) +
         F* dt);

      v[m] = (vn[m] -
            un[m] * dt/dx *
           (vn[m] - vn[m-1]) -
            vn[m] * dt/dy *
           (vn[m] - vn[m-nx]) -
            dt/(2*rho*dy) *
           (p[m+nx] - p[m-nx]) +
            nu * (dt/pow(dx,2)*
           (vn[m+1] - 2*vn[m] + vn[m-1]) +
            dt/pow(dy,2) *
           (vn[m+nx] - 2*vn[m] + vn[m-nx])));

   }else if(bId != 0 && bId != (ny-1) && tId == (nx-1)){
         //Periodic BC u @ x = 2 u[j][nx-1]
         u[m] = (un[m] -  
            un[m] * dt/dx *
           (un[m] - un[m-1]) -
            vn[m] * dt/dy *
           (un[m] - un[m-nx]) -
            dt/(2*rho*dx) *
           (p[m-nx+1] - p[m-1]) +
            nu * (dt/pow(dx,2)*
           (un[m-nx+1] - 2*un[m] + un[m-1]) +
            dt/pow(dy,2) *
           (un[m+nx] - 2*un[m] + un[m-nx])) +
            F* dt); 

         //Periodic BC v @ x = 2
         v[m] = (vn[m] -
            un[m] * dt/dx *
           (vn[m] - vn[m-1]) -
            vn[m] * dt/dy *
           (vn[m] - vn[m-nx]) -
            dt/(2*rho*dy) *
           (p[m+nx] - p[m-nx]) +
            nu * (dt/pow(dx,2)*
           (vn[m-nx+1] - 2*vn[m] + vn[m-1]) +
            dt/pow(dy,2) *
           (vn[m+nx] - 2*vn[m] + vn[m-nx])));

   }else if(bId !=0 && bId != (ny-1) && tId == 0){
      //Periodic BC u @ x = 0
      u[m] = (un[m] -  
         un[m] * dt/dx *
        (un[m] - un[m+nx-1]) -
         vn[m] * dt/dy *
        (un[m] - un[m-nx]) -
         dt/(2*rho*dx) *
        (p[m+1] - p[m+nx-1]) +
         nu * (dt/pow(dx,2)*
        (un[m+1] - 2*un[m] + un[m+nx-1]) +
         dt/pow(dy,2) *
        (un[m+nx] - 2*un[m] + un[m-nx])) +
         F* dt);

      //Periodic BC v @ x = 0
      v[m] = (vn[m] -
         un[m] * dt/dx *
        (vn[m] - vn[m+nx-1]) -
         vn[m] * dt/dy *
        (vn[m] - vn[m-nx]) -
         dt/(2*rho*dy) *
        (p[m+nx] - p[m-nx]) +
         nu * (dt/pow(dx,2)*
        (vn[m+1] - 2*vn[m] + vn[m+nx-1]) +
         dt/pow(dy,2) *
        (vn[m+nx] - 2*vn[m] + vn[m-nx]))); 
   }else if(bId == 0 || bId == (ny-1)){
      //Wall BC: u,v = 0 @ y = 0,2
      u[m] = 0;
      v[m] = 0;
   }

}

int main() {
   //Variable Declarations
   float dx = 2/(nx - 1.0);
   float dy = 2/(ny - 1.0);
   int m;
   
   //Physical Variables
   const float rho = 1.0;
   const float nu = .1;
   const float F = 1.0;
   const float dt = .01;
   
   //Initial Conditions
   float udiff = 1.0;
   int stepcount = 0;
   float sumu = 0.0;
   float sumun = 0.0;

   float *b;
   float *p;
   float *u;
   float *v;
   float *un;
   float *vn;
   float *pn;

   cudaMallocManaged(&b,ny*nx*sizeof(float));
   cudaMallocManaged(&p,ny*nx*sizeof(float));
   cudaMallocManaged(&u,ny*nx*sizeof(float));
   cudaMallocManaged(&v,ny*nx*sizeof(float));
   cudaMallocManaged(&un,ny*nx*sizeof(float));
   cudaMallocManaged(&vn,ny*nx*sizeof(float));
   cudaMallocManaged(&pn,ny*nx*sizeof(float));

   for(int i=0; i<ny*nx; i++) {
      u[i]=0.0;
      v[i]=0.0;
      un[i]=0.0;
      vn[i]=0.0;
      pn[i]=1.0;
   }

   auto tic = chrono::steady_clock::now();

   while(udiff>.001){
      for(int i=0; i<nx; i++){
         for(int j=0; j<ny; j++){
            un[j*nx+i] = u[j*nx+i];
            vn[j*nx+i] = v[j*nx+i];
         }    
      }
      
      build_up_b<<<ny,nx>>>(b,rho,dt,dx,dy,u,v);
      cudaDeviceSynchronize();
      
      //b = build_up_b(rho, dt, dx, dy, u, v);
      pressure_poisson_periodic<<<ny,nx>>>(p,pn,b, dx, dy);
      cudaDeviceSynchronize();
      
      //p = pressure_poisson_periodic(p, b, dx, dy);
      updated_u_v<<<ny,nx>>>(u,v,un,vn,p,dx,dy,dt,rho,nu,F);
      cudaDeviceSynchronize();

      sumu = 0.0;
      sumun = 0.0;
      for(int i=0;i<nx;i++){
         for(int j=0; j<ny;j++){
            m = j*nx+i;
            sumu += u[m];
            sumun += un[m];
         }
      }
      udiff = (sumu - sumun)/ sumu ;
      if(stepcount % 50 ==0){
         std::cout<<"Step: "<<stepcount<<": udiff = "<<udiff<<std::endl;
      }
      stepcount += 1;      
   }
   auto toc = chrono::steady_clock::now();
   double time = chrono::duration<double>(toc - tic).count();
   std::cout<<"Step: "<<stepcount<<": udiff = "<<udiff<<std::endl;
   std::cout<<"Step Count = " << stepcount <<std::endl;
   std::cout<<"Time = " << time <<std::endl;

}
