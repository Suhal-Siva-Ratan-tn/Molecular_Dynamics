
module para
    implicit none 
    integer, parameter :: n_part=3000
    real*8, parameter:: T = 1.0d0, m=1.0d0, kbT=1.0d0
    real*8, parameter:: sigma = 1.0d0, eps = 4.0d0, rc = 3.0d0, del_r=0.05, del_v=0.05
    real*8, parameter:: tf=50.0d0, dt=0.005d0, rs=2*rc, v_max=(8.0d0**0.50d0), v_len=4.0d0/del_v
    real*8, parameter:: lx=20.0d0, ly=20.0d0, lz=20.0d0
    real*8, parameter:: sigma6=sigma**6, sigma12=sigma**12, sigma2=sigma*sigma
    real*8, parameter:: fc=eps*((12.0d0*sigma12/(rc**13))-(6.0d0*sigma6/(rc**7)))
    real*8, parameter:: dt2=dt**2, ufc = eps*(((sigma/rc)**12)-(sigma/rc)**6)+fc*rc
    real*8, parameter:: lxby2=lx/2.0d0, lyby2=ly/2.0d0, lzby2=lz/2.0d0, drby2=(lxby2*lxby2+lyby2*lyby2+lzby2*lzby2)
    integer, parameter:: max_neigh=(n_part/(lx*ly*lz))*(4.0d0*3.140d0*rs/3), max_len_pc=lxby2/del_r, v_max_MB= (80)
    real*8 :: vx(n_part), vy(n_part), vz(n_part)
    real*8 :: x(n_part), y(n_part), z(n_part)
    real*8 :: xp(n_part), yp(n_part), zp(n_part)
    real*8 :: fx(n_part), fy(n_part), fz(n_part)
    real*8 :: xnew(n_part), ynew(n_part), znew(n_part)
    integer :: num_neigh(n_part), neighlist(n_part,max_neigh), p
    real*8 :: pc_num(n_part,max_len_pc), pc_den(n_part,max_len_pc), pc_avgden(max_len_pc), max_boltz_vel(v_max_MB)
    real*8 :: dx, dy, dz, dr
    real*8 :: avg_vx,avg_vy, avg_vz, sumv2
    real*8 :: randm, fs, force, ke, pe, te, lj, ljforces, theoryke, scalef, den, val
    real*8 :: x1, x2, y1, y2, z1, z2, tt
    integer :: i,j,k 
end module para

program md_verlet
    use para 
    implicit none
    real *8 :: start, stopp
    call cpu_time(start) 
    !Initialization of positions and velocities
    open (12, file='intial_pos.txt', action='read')
    open(13, file='energy2.dat', action='write')
    open(14, file='pc.dat', action='write')
    open(15, file='mbv.dat', action='write')
    open(16, file='mom.dat', action='write')
    tt=0.0d0
    do i = 1, n_part
        read(12,*) x(i), y(i), z(i)
    end do
    write(*,*) n_part, "Particle Positions Initialized"
    call vel_init
    !$acc enter data copyin(vx, vy, vz)
    !$acc enter data copyin(x,y,z)
    !$acc enter data copyin(xp,yp,zp)
    !$acc enter data copyin(fx,fy,fz)
    !$acc enter data copyin(xnew(:), ynew(:), znew(:))
    call calcforces
    write(*,*) n_part, "Particle Velocities Initialized"
    write(*,*) "Potential Energy", pe/dfloat(n_part)
    write(*,*) "Kinetic Energy", ke/dfloat(n_part)
    !$acc parallel loop present(xp,yp,zp,x,y,z,vx,vy,vz)
    do i = 1, n_part
        xp(i) = x(i)-vx(i)*dt
        yp(i) = y(i)-vy(i)*dt
        zp(i) = z(i)-vz(i)*dt
    end do
    
    tt=dt
    
    do while (tt<= tf)
        call calcforces 
        call update
        if (mod(int(tt/0.005d0),100)==0) call pair_correlation
        write(13,*)tt, ke/dfloat(n_part), pe/dfloat(n_part)!, (ke+pe)/dfloat(n_part)
        write(16,*)tt, avg_vx, avg_vy, avg_vz
      !  write(*,*) tt, ke/dfloat(n_part), pe/dfloat(n_part)!, (ke+pe)/dfloat(n_part)
        tt=tt+dt
    end do
    !$acc exit data delete(vx, vy, vz)
    !$acc exit data delete(x,y,z)
    !$acc exit data delete(xp,yp,zp)
    !$acc exit data delete(fx,fy,fz)
    !$acc exit data delete(xnew, ynew, znew)
    !write(*,*) pc_avgden
    !write(*,*) max_boltz_vel
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    call cpu_time(stopp) 
    write(*,'(a,f10.3,a)')  ' completed in ', stopp-start, ' seconds'
end program md_verlet
subroutine vel_init
    use para 
    implicit none
    avg_vx = 0.0d0; avg_vy=0.0d0; avg_vz=0.0d0
    sumv2=0.0d0
    ke=0.0d0
    do i = 1, n_part
        call random_number(randm); vx(i) = (randm-0.5d0)
        call random_number(randm); vy(i) = (randm-0.5d0)
        call random_number(randm); vz(i) = (randm-0.5d0)
    end do
    
    do i = 1, n_part
        avg_vx=avg_vx+vx(i)
        avg_vy=avg_vy+vy(i)
        avg_vz=avg_vz+vz(i)
    end do
    
    avg_vx=avg_vx/dfloat(n_part)
    avg_vy=avg_vy/dfloat(n_part)
    avg_vz=avg_vz/dfloat(n_part)

    !$acc parallel loop
    do i = 1, n_part
        vx(i) = (avg_vx-vx(i))
        vy(i) = (avg_vy-vy(i))
        vz(i) = (avg_vz-vz(i))
    end do 
    !$acc end parallel loop
    
    do i = 1, n_part
        sumv2=sumv2+(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
    end do
    
 
    fs = dsqrt(12.0d0*kbT/m)
    !$acc parallel loop
    do i = 1, n_part
        vx(i) = vx(i)*fs
        vy(i) = vy(i)*fs
        vz(i) = vz(i)*fs
        !$acc atomic
        ke=ke+(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
    end do
    ke = ke*0.5d0

end subroutine vel_init

! Force Calculations
subroutine calcforces 
    use para 
    implicit none
    pe = 0.0d0
   
    !$acc parallel loop present(fx,fy,fz)
    do i = 1, n_part
        fx(i) = 0.0d0; fy(i) = 0.0d0; fz(i) = 0.0d0
    end do
    !$acc end parallel loop

    !$acc parallel loop present(fx,fy,fz,x,y,z)
    do i = 1, n_part-1
        x1=x(i)
        y1=y(i)
        z1=z(i)
        !$acc loop vector
        do j = i+1, n_part     !1, num_neigh(i)   i+1, n_part
               !neighlist(i,j) j
            x2=x(j)
            y2=y(j)
            z2=z(j)
            dx = x1-x2; dy = y1-y2; dz=z1-z2
            if(abs(dx).ge.lxby2) dx=(lx-abs(dx))*((-1.0d0*dx/abs(dx)))
            if(abs(dy).ge.lyby2) dy=(ly-abs(dy))*((-1.0d0*dy/abs(dy)))
            if(abs(dz).ge.lzby2) dz=(lz-abs(dz))*((-1.0d0*dz/abs(dz)))
            
            dr = dsqrt(dx*dx+dy*dy+dz*dz)
            if (dr <= rc) then 
                lj = eps*((sigma/dr)**12-(sigma/dr)**6)-ufc+fc*dr
                
                ljforces = eps*((12.0d0*((sigma12)/(dr)**13))-(6.0d0*((sigma6/(dr)**7))))-fc
                !$acc atomic
                fx(i) = fx(i) + ljforces*(dx/dr)
                !$acc atomic
                fy(i) = fy(i) + ljforces*(dy/dr)
                !$acc atomic
                fz(i) = fz(i) + ljforces*(dz/dr)
                !$acc atomic
                fx(j) = fx(j) - ljforces*(dx/dr)
                !$acc atomic
                fy(j) = fy(j) - ljforces*(dy/dr)
                !$acc atomic
                fz(j) = fz(j) - ljforces*(dz/dr)
                !$acc atomic
                pe = pe + lj
            end if
        end do 
    end do
end subroutine calcforces

subroutine update
    use para 
    implicit none
    avg_vx=0.0d0; avg_vy=0.0d0; avg_vz=0.0d0; ke=0.0d0


! Position and Velocity Update
    !$acc parallel loop present(fx,fy,fz,x,y,z,xp,yp,zp,xnew,ynew, znew,vx,vy,vz)
    do i = 1, n_part
        xnew(i) = 2.0d0*x(i)-xp(i)+fx(i)*dt2
        ynew(i) = 2.0d0*y(i)-yp(i)+fy(i)*dt2
        znew(i) = 2.0d0*z(i)-zp(i)+fz(i)*dt2
        vx(i) = (xnew(i)-xp(i))/(2.0d0*dt)
        vy(i) = (ynew(i)-yp(i))/(2.0d0*dt)
        vz(i) = (znew(i)-zp(i))/(2.0d0*dt)
        !$acc atomic
        ke = ke + 0.5d0*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
        !$acc atomic
        avg_vx=avg_vx+(vx(i)/dfloat(n_part))
        !$acc atomic
        avg_vy=avg_vy+(vy(i)/dfloat(n_part))
        !$acc atomic
        avg_vz=avg_vz+(vz(i)/dfloat(n_part))
        
    end do
    !$acc end parallel loop
! Thermostat
    if (mod(int(tt/dt),100)==0) then
        theoryke=1.50d0*dfloat(n_part)*kbT
        scalef=dsqrt(theoryke/ke)
        vx=vx*scalef
        vy=vy*scalef
        vz=vz*scalef
        ke=0.0d0
        !$acc parallel loop present(fx,fy,fz,x,y,z,xp,yp,zp,xnew,ynew, znew,vx,vy,vz)
        do i=1,n_part
            !$acc atomic
            ke=ke+0.50d0*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
        end do
    end if
    
    !$acc parallel loop present(fx,fy,fz,x,y,z,xp,yp,zp,xnew,ynew, znew,vx,vy,vz)
    do i = 1, n_part
        xnew(i) = xp(i) + 2.0d0*vx(i)*dt
        ynew(i) = yp(i) + 2.0d0*vy(i)*dt
        znew(i) = zp(i) + 2.0d0*vz(i)*dt
        xp(i)=x(i)
        yp(i)=y(i)
        zp(i)=z(i)
        
        x(i) = xnew(i)
        y(i) = ynew(i)
        z(i) = znew(i)
 
!PBC       
        if (x(i)>lx) then 
            x(i) = x(i)-lx
            xp(i) = xp(i)-lx
        elseif(x(i)<0.0d0) then 
            x(i) = x(i)+lx  
            xp(i) = xp(i)+lx
        end if
        
        if (y(i)>ly) then 
            y(i) = y(i)-ly
            yp(i) = yp(i)-ly
        elseif(y(i)<0.0d0) then 
            y(i) = y(i)+ly 
            yp(i) = yp(i)+ly

        end if
        if (z(i)>lz) then 
            z(i) = z(i)-lz
            zp(i) = zp(i)-lz
        elseif(z(i)<0.0d0) then 
            z(i) = z(i)+lz  
            zp(i) = zp(i)+lz

        end if
    end do
    !$acc end parallel loop
end subroutine update


!PCMB
subroutine pair_correlation
    use para 
    implicit none
    max_boltz_vel=0
    do i = 1, n_part
        do j = 1, max_len_pc
            pc_num(i,j)=0
            pc_den(i,j)=0
        end do
    end do
    do i = 1, n_part-1
        x1=x(i); y1=y(i);z1=z(i)
        do j = i+1, n_part
            if(i .ne. j) then
                x2=x(j); y2=y(j); z2=z(j)
                dx = x2-x1; dy = y2-y1; dz = z2-z1
                if(abs(dx).ge.lxby2) dx=(lx-abs(dx))*((-1.0d0*dx/abs(dx)))
                if(abs(dy).ge.lyby2) dy=(ly-abs(dy))*((-1.0d0*dy/abs(dy)))
                if(abs(dz).ge.lzby2) dz=(lz-abs(dz))*((-1.0d0*dz/abs(dz)))
                dr = dsqrt(dx*dx+dy*dy+dz*dz)
                k  = int(dr/del_r)+1
                if(k<=max_len_pc) then
                    den        = 4*(4*ATAN(1.0d0))*dr*dr*del_r
                    !write(*,*) tt,i,k
                    pc_num(i,k)= pc_num(i,k)+1
                    pc_den(i,k)= pc_den(i,k)+1.0d0/den
                    pc_num(j,k)= pc_num(j,k)+1
                    pc_den(j,k)= pc_den(j,k)+1.0d0/den
                end if
            end if
        end do
        val = dsqrt(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
        k   = val/del_v+1
        if(k<=v_max_MB) then
            max_boltz_vel(k) =  max_boltz_vel(k)+1
        end if
    end do
    max_boltz_vel =  max_boltz_vel/dfloat(n_part)
    pc_avgden     =  sum(pc_den,dim=1)/dfloat(n_part)
    pc_avgden     =  pc_avgden*max_len_pc/sum(pc_avgden)
    
    write(14,*) tt
    write(15,*) tt
    do i = 1, max_len_pc
         write(14,*) (i*del_r), pc_avgden(i)
    end do

    write(14,*)
    write(14,*)

    
    do i = 1, v_max_MB
         write(15,*) (i*del_v), max_boltz_vel(i)
    end do

    write(15,*)
    write(15,*)
         
end subroutine pair_correlation
