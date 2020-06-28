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

/*__device__ __managed__ std::vector<float> b;
__device__ __managed__ std::vector<float> p;
__device__ __managed__ std::vector<float> u;
__device__ __managed__ std::vector<float> un;
__device__ __managed__ std::vector<float> v;
__device__ __managed__ std::vector<float> vn;*/

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

std::vector<float> pressure_poisson_periodic(std::vector<float> p, std::vector<float> b, float dx, float dy){
   std::vector<float> pn(nx*ny, 1);
   int m;
   
   for(int q=0; q<nit; q++){
      for(int i=0;i<nx;i++){
         for(int j=0;j<ny;j++){
             pn[j*nx+i] = p[j*nx+i];
         }
      }

      for(int i=1;i<nx-1;i++){
         for(int j=1;j<ny-1;j++){
            m = j*nx+i;
            p[m]=(((pn[m+1]+pn[m-1])*pow(dy,2)+
                      (pn[m+nx]+pn[m-nx])*pow(dx,2))/
                      (2*(pow(dx,2)+pow(dy,2)))-
                      pow(dx,2)*pow(dy,2)/(2*(pow(dx,2)+pow(dy,2)))*b[m]);
         }
      }
      
      for(int j=1;j<ny-1;j++){
         //Periodic BC Pressure @ x = 2
         m = j*nx +nx-1;
         p[m]=(((pn[m-nx+1]+pn[m-1])*pow(dy,2)+
                     (pn[m+nx]+pn[m-nx])*pow(dx,2))/
                     (2*(pow(dx,2)+pow(dy,2)))-
                     pow(dx,2)*pow(dy,2)/(2*(pow(dx,2)+pow(dy,2)))*b[m]);
     
         //Periodic BC Pressure @ x = 0
         m = j*nx;
         p[m]=(((pn[m+1]+pn[m+nx-1])*pow(dy,2)+
                     (pn[m+nx]+pn[m-nx])*pow(dx,2))/
                     (2*(pow(dx,2)+pow(dy,2)))-
                     pow(dx,2)*pow(dy,2)/(2*(pow(dx,2)+pow(dy,2)))*b[m]);
      }

      for(int i=0;i<nx;i++){
         //Wall boundary conditions, pressure
         p[nx*(ny-1)+i] = p[nx*(ny-2)+i]; //dp/dy = 0 at y = 2
         p[i] = p[nx+i]; //dp/dy = 0 at y = 0
      }
      
   }
   return p;
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
   const int rho = 1;
   const float nu = .1;
   const int F = 1;
   const float dt = .01;
   
   //Initial Conditions
   std::vector<float> u;
   std::vector<float> un;
   std::vector<float> v;
   std::vector<float> vn;
   std::vector<float> p;
   std::vector<float> pn;
   //float *u;
   //float *v;
   //std::vector<float> b;
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
         u.push_back(0);
         un.push_back(0);
         v.push_back(0);
         vn.push_back(0);
         p.push_back(1);
         pn.push_back(1);
         //b.push_back(0);
      }
   }
   std::cout<<u.size()<<std::endl;
   std::cout<<un.size()<<std::endl;
   float udiff = 1.0;
   int stepcount = 0;
   float sumu = 0.0;
   float sumun = 0.0;
   
   while(udiff>.001){
      for(int i=0; i<nx; i++){
         for(int j=0; j<ny; j++){
            un[j*nx+i] = u[j*nx+i];
            vn[j*nx+i] = v[j*nx+i];
            //std::cout<<"un "<<un[j*nx+i]<<std::endl;
            //std::cout<<"vn "<<vn[j*nx+i]<<std::endl;
         }    
      }

      float *b;
      
      cudaMallocManaged(&b,dy*dx*sizeof(float));
      build_up_b<<<dy,dx>>>(b,rho,dt,dx,dy,&u,&v);
      
      //b = build_up_b(rho, dt, dx, dy, u, v);
      //p = pressure_poisson_periodic(p, b, dx, dy);
      
      /*for(int i=0; i<nx; i++){
         for(int j=0; j<nx; j++){
            std::cout<<"b "<<b[j*nx+i]<<std::endl;
            std::cout<<"p "<<p[j*nx+i]<<std::endl;
         }
      }*/
      
      
      
      for(int i=1;i<nx-1;i++){
         for(int j=1;j<ny-1;j++){
            m = j*nx + i;

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
         }
      }

      
      for(int j=1;j<ny-1;j++){
            m = j*nx + nx-1 ;
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


            m = j*nx; 
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
      }
      
      //Wall BC: u,v = 0 @ y = 0,2
      
      for(int i=0; i<nx; i++){
         u[i] = 0;
         u[nx*(ny-1)+i] = 0;
         v[i] = 0;
         v[nx*(ny-1)+i] = 0;
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
   cudaDeviceSynchronize();
   std::cout<< stepcount <<std::endl;
   //cudaFree(b);

}
