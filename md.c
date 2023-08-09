#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

// #include <cuda/std/atomic>

#define n_part 3000

const double T=1.0, m=1.0, kbT=1.0;
double lx=20.0, ly=20.0, lz=20.0, sigma=1.0, eps =4.0, rc= 3.0, del_r=0.05, del_v=0.05;
const double tf=50.0, dt=0.005;
const double dt2 = dt*dt;


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

double calcforces(double* fx, double* fy, double* fz,double* x, double* y, double* z){
    double dx, dy, dz, dr;
    double x1,x2,y1,y2,z1,z2;
    double lj, ljforces, pe;
    double lxby2=lx/2.0, lyby2=ly/2.0, lzby2=lz/2.0;
    double sigma6 = pow(sigma, 6), sigma12 = pow(sigma, 12), sigma2 = sigma * sigma;
    double fc = eps * ((12.0 * sigma12 / pow(rc, 13)) - (6.0 * sigma6 / pow(rc, 7)));
    double ufc = eps * ((pow(sigma / rc, 12)) - pow(sigma / rc, 6)) + fc * rc;

    int i,j;
    pe = 0.0;

    #pragma acc parallel loop
    for(i=0; i<n_part; i++){
        fx[i] = 0.0; 
        fy[i] = 0.0;
        fz[i] = 0.0;
    }

    #pragma acc parallel loop
    for(i=0; i<n_part; i++){
    x1 = x[i];
    y1 = y[i];
    z1 = z[i];
        #pragma acc loop vector
        for (j=i+1;j <n_part;j++) {
            x2 = x[j];
            y2 = y[j];
            z2 = z[j];
            dx = x1-x2;
            dy = y1-y2;
            dz = z1-z2;

            if (fabs(dx) >= lxby2) dx = (lx-fabs(dx))*((-1.0*dx)/fabs(dx));
            if (fabs(dy) >= lyby2) dy = (ly-fabs(dy))*((-1.0*dy)/fabs(dy));
            if (fabs(dz) >= lzby2) dz = (lz-fabs(dz))*((-1.0*dz)/fabs(dz));

            dr = sqrt(dx*dx + dy*dy + dz*dz);
            if (dr <= rc) {
                lj = eps*((pow(sigma/dr,12))-(pow(sigma/dr, 6)))-ufc+fc*dr;
                ljforces = eps*((12.0*(sigma12/pow(dr,13)))-(6.0*(sigma6/pow(dr,7))))-fc;

                #pragma acc atomic update
                fx[i] = fx[i]+ljforces*(dx/dr);
                #pragma acc atomic update
                fy[i] = fy[i]+ljforces*(dy/dr);
                #pragma acc atomic update
                fz[i] = fz[i]+ljforces*(dz/dr);
                #pragma acc atomic update
                fx[j] = fx[j]-ljforces*(dx/dr);
                #pragma acc atomic update
                fy[j] = fy[j]-ljforces*(dy/dr);
                #pragma acc atomic update
                fz[j] = fz[j]-ljforces*(dz/dr);
                #pragma acc atomic update
                pe = pe+lj;
            }
        }
    }
    return pe/n_part;
}

double update(double*fx, double* fy, double* fz, double* x, double* y, double* z, double* xp, double* yp, double* zp, double* xnew, double* ynew, double* znew, double* vx,double* vy,double* vz,double tt){
    // #pragma acc routine gang
    int i;
    double ke;
    ke =0.0;
    int thermo_check= (tt/0.005);
    #pragma acc parallel loop reduction(+:ke)
    for(i=0;i<n_part;i++){
        xnew[i] = 2.0*x[i]-xp[i]+fx[i]*dt2;
        ynew[i] = 2.0*y[i]-yp[i]+fy[i]*dt2;
        znew[i] = 2.0*z[i]-zp[i]+fz[i]*dt2;

        vx[i] = (xnew[i]-xp[i])/(2.0*dt);
        vy[i] = (ynew[i]-yp[i])/(2.0*dt);
        vz[i] = (znew[i]-zp[i])/(2.0*dt);

        ke = ke+0.5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
        
    }

    if (thermo_check %100 == 0) printf("%lf\t", ke/n_part);
    // printf("Thermo check at %lf : %d\n",tt, thermo_check );

    // Thermostat

    if (thermo_check %100 == 0){
        double theoryke = 1.50*n_part*kbT;
        double scalef = sqrt(theoryke)/sqrt(ke);
        
        #pragma acc parallel loop  
        for(int i=0; i<n_part; i++){
            vx[i] = vx[i]*scalef;
            vy[i] = vy[i]*scalef;
            vz[i] = vz[i]*scalef;
        }
        // printf("%lf %lf\n", theoryke/ke,scalef);
        double lol=0.0;
        #pragma acc parallel loop reduction(+:lol)
        for(int i=0; i<n_part; i++){
            
            lol = lol+0.5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
        }
        
        ke=lol;
        // printf("%lf %lf\n", lol,ke);
    }
    
    #pragma acc loop gang
    for(i=0; i<n_part; i++){
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
    // printf("%lf\n", ke/n_part);
    return ke/n_part;
}

int main(){
    
    // const int max_neigh=(n_part/(lx*ly*lz))*(4.0*3.140*rs/3), max_len_pc=lxby2/del_r, v_max_MB= (80);

    double vx[n_part], vy[n_part], vz[n_part];
    double x[n_part], y[n_part], z[n_part];
    double xp[n_part], yp[n_part], zp[n_part];
    double fx[n_part], fy[n_part], fz[n_part];
    double xnew[n_part], ynew[n_part], znew[n_part];

    double ke, pe;
    double tt;
    int i;

    printf("%lf %lf\n", T, m);
    FILE *f1 = fopen("intial_pos.txt", "r");
    FILE *f2 = fopen("energy.dat", "w");
    for (i = 0; i < n_part; i++) {
        fscanf(f1, "%lf %lf %lf", &x[i], &y[i], &z[i]);
        // printf("%lf %lf %lf\n", x[i], y[i], z[i]);
    }
    fclose(f1);
    
    ke = velinit(vx,vy,vz,n_part);
    // for (i = 0; i < n_part; i++) {
    //     printf("Particle Number %d: %lf\t%lf\t%lf\n", i, vx[i],vy[i],vz[i]);

    // }
    // for(i=0; i<n_part; i++){
    //     printf("%d %lf\n",i,fx[i]);
    // }
    pe=calcforces(fx,fy,fz,x,y,z);
    printf("KE: %lf \t PE: %lf\n", ke,pe);
    for (i = 0; i < n_part; i++) {
        xp[i] = x[i] - vx[i] * dt;
        yp[i] = y[i] - vy[i] * dt;
        zp[i] = z[i] - vz[i] * dt;
    }

    ke=update(fx,fy,fz,x,y,z,xp,yp,zp,xnew,ynew,znew,vx,vy,vz, 0.05);
    // for(i=0; i<n_part; i++){
    //     printf("%d %lf\n",i,fx[i]);
    // }
    tt = dt;

    while(tt<=tf){
        pe=calcforces(fx,fy,fz,x,y,z);
        ke=update(fx,fy,fz,x,y,z,xp,yp,zp,xnew,ynew,znew,vx,vy,vz, tt);
        printf("%lf\t KE: %lf \t PE: %lf\n", tt,ke,pe);
        fprintf(f2,"%lf\t %lf \t %lf \t %lf \n", tt,ke,pe,ke+pe);
        tt = tt+dt;
    }
    return 0;
}