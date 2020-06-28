#include <cstdio>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;

const int nx = 5;
const int ny = 5;
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
   printf("Hello GPUhahah\n");
   if(bId != 0 && bId != (ny-1) && tId != 0 && tId != (nx-1)){
      printf("Hello GPU\n");
      u[m] = 1;/*(un[m] -
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
         F* dt);*/ 

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
   //std::vector<float> x;
   //std::vector<float> y;
   //std::vector<std::vector<float> > X;
   //std::vector<std::vector<float> > Y;
   int m;
   
   //Physical Variables
   const float rho = 1.0;
   const float nu = .1;
   const float F = 1.0;
   const float dt = .01;
   
   //Initial Conditions
   //std::vector<float> u;
   std::vector<float> un;
   //std::vector<float> v;
   std::vector<float> vn;
   std::vector<float> pn;
   //float *u;
   //float *v;
   /*
   for(int i=0;i<nx;i++){
      x.push_back((2-0)*i/(nx-1));
   }
   
   for(int i=0;i<ny;i++){
      y.push_back((2-0)*i/(ny-1));
   }*/
   
  for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){
         //X[i][j] = x[i];
         //Y[i][j] = y[j];   
         //u.push_back(0);
         un.push_back(0);
         //v.push_back(0);
         vn.push_back(0);
         //p.push_back(1);
         pn.push_back(1);
         //b.push_back(0);
      }
   }

   float udiff = 1.0;
   int stepcount = 0;
   float sumu = 0.0;
   float sumun = 0.0;
   float *b;
   float *p;
   float *u;
   float *v;

   cudaMallocManaged(&b,ny*nx*sizeof(float));
   cudaMallocManaged(&p,ny*nx*sizeof(float));
   cuda_Error_t err1=cudaMallocManaged(&u,ny*nx*sizeof(float));
   cuda_Error_t err2=cudaMallocManaged(&v,ny*nx*sizeof(float));
   std::cout<<err1<<" "<<err2<<std::endl;
   for(int i=0; i<ny*nx; i++) {
      u[i]=0.0;
      v[i]=0.0;
   }

   while(udiff>.001){
      /*for(int i=0; i<nx; i++){
         for(int j=0; j<ny; j++){
            un[j*nx+i] = u[j*nx+i];
            vn[j*nx+i] = v[j*nx+i];
         }    
      }*/
      
      build_up_b<<<ny,nx>>>(b,rho,dt,dx,dy,u,v);
      //b = build_up_b(rho, dt, dx, dy, u, v);
      pressure_poisson_periodic<<<ny,nx>>>(p,pn.data(),b, dx, dy);
      //p = pressure_poisson_periodic(p, b, dx, dy);
      updated_u_v<<<ny,nx>>>(u,v,un.data(),vn.data(),p,dx,dy,dt,rho,nu,F);
      cudaDeviceSynchronize();

      for(int i=0;i<nx;i++){
         for(int j=0; j<ny;j++){
            std::cout<<"u "<<u[j*nx+i]<<std::endl;
            std::cout<<"un "<<un[j*nx+i]<<std::endl;
         }
      }

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
      std::cout<<"udiff "<<udiff<<std::endl;
      stepcount += 1;  
       
   }
   
   std::cout<< stepcount <<std::endl;
   //cudaFree(b);

}
