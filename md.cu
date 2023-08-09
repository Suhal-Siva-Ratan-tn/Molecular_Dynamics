#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>

// #include <cuda/std/atomic>

#define n_part 2048
#define lx 20.0
#define ly 20.0
#define lz 20.0

// __device__ double lx = 20.0;
const double kbT=1.0;
__device__  double sigma=1.0, eps =4.0, rc= 3.0, del_r=0.05, del_v=0.05;
const double tf=50.0, dt=0.005, v_max=sqrt(8);



double velinit(double* vx, double* vy, double* vz,int n){
    double avg_vx, avg_vy, avg_vz,ke;
    int i;

    avg_vx = 0.0; 
    avg_vy = 0.0; 
    avg_vz = 0.0;

    srand(time(NULL));
    for(i=0;i<n;i++){
        vx[i]= (double)rand()/RAND_MAX-0.5;
        vy[i]= (double)rand()/RAND_MAX-0.5;
        vz[i]= (double)rand()/RAND_MAX-0.5;
    }

    for(i=0;i<n;i++){
        avg_vx += vx[i];
        avg_vy += vy[i];
        avg_vz += vz[i];
    }

    double n_real = (double) n;
    avg_vx = avg_vx/n_real;
    avg_vy = avg_vy/n_real;
    avg_vz = avg_vz/n_real;
    
    for(i=0;i<n;i++){
        vx[i] = avg_vx - vx[i];
        vy[i] = avg_vy - vy[i];
        vz[i] = avg_vz - vz[i];
    }

    double sumv2=0.0;

    for(i=0;i<n;i++){
        sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }

    double fs = sqrt(12.0);

    for(i=0;i<n;i++){
        vx[i] *= fs;
        vy[i] *= fs;
        vz[i] *= fs;
        ke = ke + vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }

    ke *= 0.5;

    return ke/n_real;
    
}

__global__ void calcforces(double* lj, double* ljforces_x, double* ljforces_y, double* ljforces_z, double* x, double* y, double* z){
    // unsigned int tid = threadIdx.x;


    //Initializing required parameters
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    int idy = threadIdx.y + blockDim.y * blockIdx.y;
    int id2 = idx * n_part + idy;
    double sigma6 = pow(sigma, 6), sigma12 = pow(sigma, 12);
    double fc = eps * ((12.0 * sigma12 / pow(rc, 13)) - (6.0 * sigma6 / pow(rc, 7)));
    double ufc = eps * ((pow(sigma / rc, 12)) - pow(sigma / rc, 6)) + fc * rc;


    if(idx < n_part && idy < n_part){
        
        double dx = x[idx] - x[idy];
        double dy = y[idx] - y[idy];
        double dz = z[idx] - z[idy];
        double ljforces = 0.0;

        //Nearest-Image Convention
        if (fabs(dx) >= lx/2) dx = (lx-fabs(dx))*((-1.0*dx)/fabs(dx));
        if (fabs(dy) >= ly/2) dy = (ly-fabs(dy))*((-1.0*dy)/fabs(dy));
        if (fabs(dz) >= lz/2) dz = (lz-fabs(dz))*((-1.0*dz)/fabs(dz));

        double dr2 = dx*dx + dy*dy + dz*dz;
        double dr  = sqrt(dr2);
        double dr6  = dr2*dr2*dr2;
        double dr12 = dr6*dr6;
        

       
        //Force Calculation based on the cutoff distance rc
        if (dr <= rc && idx != idy) {

            lj[id2] = eps*((sigma12/dr12)-(sigma6/dr6))-ufc+fc*dr;
            ljforces = eps*((12.0*(sigma12/(dr12*dr)))-(6.0*(sigma6/(dr6*dr))))-fc;
            
            ljforces_x[id2] = ljforces*dx/dr;
            ljforces_y[id2] = ljforces*dy/dr;
            ljforces_z[id2] = ljforces*dz/dr;
           
        }

        else{

            lj[id2] = 0.0;
            ljforces_x[id2] = 0.0;
            ljforces_y[id2] = 0.0;
            ljforces_z[id2] = 0.0;

        }
        __syncthreads();
        
    }
}

__global__ void updatepos(double* vsq, double*fx, double* fy, double* fz, double* x, double* y, double* z, double* xp, double* yp, double* zp, double* xnew, double* ynew, double* znew, double* vx,double* vy,double* vz){
 
    const double dt2 = dt*dt;
    unsigned int i = threadIdx.x + blockDim.x * blockIdx.x;

    //Updating Postions and Velocities(Verlet's Algorithm)
    if(i<n_part){
        xnew[i] = 2.0*x[i]-xp[i]+fx[i]*dt2;
        ynew[i] = 2.0*y[i]-yp[i]+fy[i]*dt2;
        znew[i] = 2.0*z[i]-zp[i]+fz[i]*dt2;
        
        
        vx[i] = (xnew[i]-xp[i])/(2.0*dt);
        vy[i] = (ynew[i]-yp[i])/(2.0*dt);
        vz[i] = (znew[i]-zp[i])/(2.0*dt);

        //Vsq is estimated to calculate Kinetic Energy
        vsq[i] = vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
    }
    
}

__global__ void thermostat(double scalef, double* vx, double* vy, double* vz, double* vsq){
    unsigned int i = threadIdx.x + blockDim.x * blockIdx.x;

    if(i<n_part){
        //Velocities are rescaled to maintain temp
        vx[i] = vx[i]*scalef;
        vy[i] = vy[i]*scalef;
        vz[i] = vz[i]*scalef;  
        
        vsq[i] = vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
    }
}

__global__ void PBC(double* x, double* y, double* z, double* xp, double* yp, double* zp, double* xnew, double* ynew, double* znew, double* vx,double* vy,double* vz){
  
    unsigned int i = threadIdx.x + blockDim.x * blockIdx.x;

    if(i<n_part){

        xnew[i] = xp[i]+2.0*vx[i]*dt;
        ynew[i] = yp[i]+2.0*vy[i]*dt;
        znew[i] = zp[i]+2.0*vz[i]*dt;
        xp[i] = x[i];
        yp[i] = y[i];
        zp[i] = z[i];

        x[i] = xnew[i];
        y[i] = ynew[i];
        z[i] = znew[i];

        // PBC
        if (x[i] >= lx) {
            x[i] = x[i]-lx;
            xp[i] = xp[i]-lx;
        }
        else if (x[i] < 0.0) {
            x[i] = x[i]+lx;
            xp[i] = xp[i]+lx;
        }

        if (y[i] >= ly) {
            y[i] = y[i]-ly;
            yp[i] = yp[i]-ly;
        }
        else if (y[i] < 0.0) {
            y[i] = y[i]+ly;
            yp[i] = yp[i]+ly;
        }

        if (z[i] >= lz) {
            z[i] = z[i]-lz;
            zp[i] = zp[i]-lz;
        }
        else if (z[i] < 0.0) {
            z[i] = z[i]+lz;
            zp[i] = zp[i]+lz;
        }
    }
} 


int main(){

    double vx[n_part], vy[n_part], vz[n_part];
    double x[n_part], y[n_part], z[n_part];
    double xp[n_part], yp[n_part], zp[n_part];

    double  ke, pe, theoryke, scalef;
    double tt;
    int i;

    //Getting Initial Positions from a text file
    FILE *f1 = fopen("intial_pos.txt", "r");

    for (i = 0; i < n_part; i++) {
        fscanf(f1, "%lf %lf %lf", &x[i], &y[i], &z[i]);
    }
    fclose(f1);
    
    //Initializing Velocities
    ke = velinit(vx,vy,vz,n_part);
                                                                        
    //Initializing Device arrays
    double *fx_d,*fy_d,*fz_d;
    double *x_d,*y_d,*z_d;
    double *xp_d,*yp_d,*zp_d,*xnew_d,*ynew_d,*znew_d;
    double *vx_d,*vy_d,*vz_d,*vsq;
    double *lj_d, *ljforces_x, *ljforces_y, *ljforces_z;
    double *d_vec;
    int size1 = n_part*sizeof(double);
    int size2 = n_part*n_part*sizeof(double);
    int tpb = 32;
    dim3 blocksize(tpb,tpb);
    dim3 gridsize((n_part+blocksize.x-1)/blocksize.x, (n_part+blocksize.y-1)/blocksize.y);

    dim3 block_pos(8*tpb);
    dim3 gridsize_pos((n_part+block_pos.x-1)/block_pos.x);

    cublasHandle_t handle;
    cublasCreate(&handle);
    double h_vec[n_part];
    const double alpha = 1.0;
    const double beta = 0.0;

    for(int i=0; i< n_part ;++i){
        h_vec[i] = 1.0;
        
    }
    
    
    cudaMalloc((void**)&x_d, size1);
    cudaMalloc((void**)&y_d, size1);
    cudaMalloc((void**)&z_d, size1);
    cudaMalloc((void**)&xp_d, size1);
    cudaMalloc((void**)&yp_d, size1);
    cudaMalloc((void**)&zp_d, size1);
    cudaMalloc((void**)&xnew_d, size1);
    cudaMalloc((void**)&ynew_d, size1);
    cudaMalloc((void**)&znew_d, size1);

    cudaMalloc((void**)&vx_d, size1);
    cudaMalloc((void**)&vy_d, size1);
    cudaMalloc((void**)&vz_d, size1);
    cudaMalloc((void**)&vsq, size1);


    cudaMalloc((void**)&fx_d, size1);
    cudaMalloc((void**)&fy_d, size1);
    cudaMalloc((void**)&fz_d, size1);

    cudaMalloc((void**)&lj_d, size2);
    cudaMalloc((void**)&ljforces_x, size2);
    cudaMalloc((void**)&ljforces_y, size2);
    cudaMalloc((void**)&ljforces_z, size2);
    cudaMalloc(&d_vec, n_part * sizeof(double));


    cudaMemcpy(x_d, x, size1, cudaMemcpyHostToDevice);
    cudaMemcpy(y_d, y, size1, cudaMemcpyHostToDevice);
    cudaMemcpy(z_d, z, size1, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec, h_vec, size1, cudaMemcpyHostToDevice);


    //Initial force calculation
    calcforces<<<gridsize, blocksize>>>(lj_d,ljforces_x,ljforces_y,ljforces_z,x_d,y_d,z_d);

    cublasDgemv(handle, CUBLAS_OP_T, n_part, n_part, &alpha, ljforces_x, n_part, d_vec, 1, &beta, fx_d, 1);
    cublasDgemv(handle, CUBLAS_OP_T, n_part, n_part, &alpha, ljforces_y, n_part, d_vec, 1, &beta, fy_d, 1);
    cublasDgemv(handle, CUBLAS_OP_T, n_part, n_part, &alpha, ljforces_z, n_part, d_vec, 1, &beta, fz_d, 1);
    

    cublasDasum(handle, n_part*n_part, lj_d, 1, &pe);
    pe = -pe*0.5/n_part;
    // printf("%lf\t", pe);

    for (i = 0; i < n_part; i++) {
        xp[i] = x[i] - vx[i] * dt;
        yp[i] = y[i] - vy[i] * dt;
        zp[i] = z[i] - vz[i] * dt;
    }
    
    cudaMemcpy(xp_d, xp, size1, cudaMemcpyHostToDevice);
    cudaMemcpy(yp_d, yp, size1, cudaMemcpyHostToDevice);
    cudaMemcpy(zp_d, zp, size1, cudaMemcpyHostToDevice);

    cudaMemcpy(vx_d, vx, size1, cudaMemcpyHostToDevice);
    cudaMemcpy(vy_d, vy, size1, cudaMemcpyHostToDevice);
    cudaMemcpy(vz_d, vz, size1, cudaMemcpyHostToDevice);

      
    tt = dt;
    
    int thermo_check;
    theoryke = 1.50*kbT;
     
    while(tt <= tf){

        // calculating lj potential and ljforces component matrices
        calcforces<<<gridsize, blocksize>>>(lj_d,ljforces_x,ljforces_y,ljforces_z,x_d,y_d,z_d);
    
        //Force calculation using Cublas
        cublasDgemv(handle, CUBLAS_OP_T, n_part, n_part, &alpha, ljforces_x, n_part, d_vec, 1, &beta, fx_d, 1);
        cublasDgemv(handle, CUBLAS_OP_T, n_part, n_part, &alpha, ljforces_y, n_part, d_vec, 1, &beta, fy_d, 1);
        cublasDgemv(handle, CUBLAS_OP_T, n_part, n_part, &alpha, ljforces_z, n_part, d_vec, 1, &beta, fz_d, 1);
     
        //Calculating Potential energy per particle
        cublasDasum(handle, n_part*n_part, lj_d, 1, &pe);
        pe = -pe*0.5/n_part;

        //Updating Positions
        updatepos<<<gridsize_pos,block_pos>>>(vsq, fx_d, fy_d, fz_d, x_d, y_d, z_d, xp_d, yp_d, zp_d, xnew_d, ynew_d, znew_d, vx_d, vy_d, vz_d);
        cublasDasum(handle, n_part, vsq, 1, &ke);
        ke = 0.5*ke/n_part;
        
        //Thermostat
        thermo_check = tt/dt;
        if (thermo_check %100 == 0){
                scalef = sqrt(theoryke/ke);
                thermostat<<<gridsize_pos,block_pos>>>(scalef,vx_d,vy_d,vz_d,vsq);
                cublasDasum(handle, n_part, vsq, 1, &ke);
                ke = 0.5*ke/n_part;
        }    


        //Checking Periodic Boundary conditions
        PBC<<<gridsize_pos,block_pos>>>(x_d, y_d, z_d, xp_d, yp_d, zp_d, xnew_d, ynew_d, znew_d, vx_d, vy_d, vz_d);
        printf("%lf\t KE: %lf \t PE: %lf\n", tt,ke,pe);
        
        tt = tt+dt;
    }

    cublasDestroy(handle);
    return 0;
}

