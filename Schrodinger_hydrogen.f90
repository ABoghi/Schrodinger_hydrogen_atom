!!!**********************************************************************
!!!*						                	                        *
!!!*               Schrodinger Hydrogen Atom     	         	        *
!!!*							        	                            *
!!!*                Author: Dr. Andrea Boghi  		                    *
!!!*							        	                            *
!!!*  d2psi/ds2 + d2psi/ds2 + d2psi/ds2 + (2/r + PI/np**2)*psi = 0.d0   *
!!!*								                                    *
!!!**********************************************************************
    
Program main_Schrodinger_Hydrogen_Atom
    implicit none
    real*8, allocatable :: x(:,:,:),y(:,:,:),z(:,:,:)
    real*8, ALLOCATABLE :: r(:,:,:),phi(:,:,:),theta(:,:,:)
    real*8, ALLOCATABLE :: psi_r(:,:,:),psi_i(:,:,:),Pr(:,:,:)
    integer n,l,m,ns,i,j,k
    real*8 Ls,PI,h,h_bar,q_e,m_e,epsilon_0,ds,ds2,a_1,E_1,res
    real*8 Pa,v_1, c_1
    logical flag
    parameter(PI=4.d0*datan2(1.d0,1.d0))
    parameter(h=6.62607015d-34) !!! [m2 * kg / s]
    parameter(h_bar=h/(2.d0*PI))
    parameter(q_e=1.60217663d-19) !!! [C]
    parameter(m_e=9.1093837d-31) !!! [kg]
    parameter(epsilon_0=8.854187812813d-12) !!! [F / m]
    parameter(a_1=epsilon_0*(h**2.d0)/(PI*m_e*q_e**2.d0))
    parameter(E_1=m_e*(q_e**4.d0)/(8.d0*(epsilon_0*h)**2.d0))
    parameter(v_1=(q_e**2.d0)/(2.d0*epsilon_0*h))
    parameter(c_1 = 299792458.d0) !!! [m/s]
    CHARACTER(len=80)::fname
    CHARACTER(len=80)::fname_res
    real*8 x_mean,y_mean,z_mean,x_std,y_std,z_std
    real*8 px_mean,py_mean,pz_mean,px_std,py_std,pz_std

    open(1,file='imp_hyd.dat')
    read(1,*) n
    read(1,*) l
    read(1,*) m
    read(1,*) ns
    read(1,*) Ls
    close(1)

    print*, ' a_1 = ',a_1, ' [m]; E_1 = ',E_1,' [J]; v_1 = ',v_1,' [m/s]; beta = ',v_1/c_1

    allocate(x(1:ns,1:ns,1:ns),y(1:ns,1:ns,1:ns),z(1:ns,1:ns,1:ns))
    ALLOCATE(psi_r(1:ns,1:ns,1:ns),psi_i(1:ns,1:ns,1:ns),Pr(1:ns,1:ns,1:ns))
    allocate(r(1:ns,1:ns,1:ns),phi(1:ns,1:ns,1:ns),theta(1:ns,1:ns,1:ns))

    call grid(ns,Ls,x,y,z,r,theta,phi,ds,ds2)

    call hydrogen_wave_function(n,l,m,ns,ds,theta,phi,r,psi_r,psi_i,Pr)

    call Probability(x,ds,Pr,ns,x_mean)
    call Probability(y,ds,Pr,ns,y_mean)
    call Probability(z,ds,Pr,ns,z_mean)

    print*, 'x_mean = ',x_mean, '; y_mean = ',y_mean, '; z_mean = ',z_mean

    call Probability((x-x_mean)**2,ds,Pr,ns,x_std)
    call Probability((y-y_mean)**2,ds,Pr,ns,y_std)
    call Probability((z-z_mean)**2,ds,Pr,ns,z_std)

    x_std = x_std**0.5d0
    y_std = y_std**0.5d0
    z_std = z_std**0.5d0
    print*, 'x_std = ',x_std, '; y_std = ',y_std, '; z_std = ',z_std

    call mean_momentum(psi_r,psi_i,ds,ns,px_mean,py_mean,pz_mean)
    print*, 'px_mean = ',px_mean, '; py_mean = ',py_mean, '; pz_mean = ',pz_mean

    call momentum_standard_deviation(psi_r,psi_i,Pr,ds,ns,px_mean,py_mean,pz_mean,px_std,py_std,pz_std)
    print*, 'px_std = ',px_std, '; py_std = ',py_std, '; pz_std = ',pz_std

    write(fname,111)'Schrodinger_Hydrogen_Solution_n=',int(n),'_l=',int(l),'_m=',int(m),'.csv'
    111 format(a32,i2,a3,i2,a3,i2,a4)

    open(14,file=fname,form='formatted')
    write(14,*) '"Point:0","Point:1","Point:2","psi_r","psi_i","Pr"'
    do k =1,ns
        do j =1,ns
            do i =1,ns
                write(14,101) x(i,j,k),',',y(i,j,k),',',z(i,j,k),',', &
                psi_r(i,j,k),',',psi_i(i,j,k),',',Pr(i,j,k)
            enddo
        enddo
    enddo
    close(14)

    101 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10)


    end


!!!*************************************************
!!!*						         	             *
!!!*            grid                       *
!!!*								                 *
!!!*************************************************

subroutine grid(ns,Ls,x,y,z,r,theta,phi,ds,ds2)
    implicit none
    integer, intent(in) :: ns
    real*8, intent(in) :: Ls
    real*8, intent(out) :: x(1:ns,1:ns,1:ns),y(1:ns,1:ns,1:ns),z(1:ns,1:ns,1:ns),ds,ds2
    real*8, intent(out) :: r(1:ns,1:ns,1:ns),phi(1:ns,1:ns,1:ns),theta(1:ns,1:ns,1:ns)
    integer i,j,k
    
    ds = Ls/(ns-1)
    ds2 = ds**2.d0

    do k=1,ns
        do j=1,ns
            do i=1,ns
                x(i,j,k) = ds*(i-1) -Ls/2.d0
                y(i,j,k) = ds*(j-1) -Ls/2.d0
                z(i,j,k) = ds*(k-1) -Ls/2.d0
            enddo
        enddo
    enddo

    do k = 1,ns
        do j = 1,ns
            do i = 1,ns
                call radius(x(i,j,k),y(i,j,k),z(i,j,k),r(i,j,k))
                call inclination(x(i,j,k),y(i,j,k),z(i,j,k),theta(i,j,k))
                call azimuth(x(i,j,k),y(i,j,k),phi(i,j,k))
            enddo
        enddo
    enddo

    end

!!!*************************************************
!!!*						         	             *
!!!*            radius                       *
!!!*								                 *
!!!*************************************************

subroutine radius(x,y,z,r)
    implicit none
    real*8, intent(in) :: x,y,z
    real*8, intent(out) :: r

    r = dsqrt(x**2.d0 + y**2.d0 + z**2.d0)

    end

!!!*************************************************
!!!*						         	             *
!!!*            inclination                      *
!!!*								                 *
!!!*************************************************

subroutine inclination(x,y,z,theta)
    implicit none
    real*8, intent(in) :: x,y,z
    real*8, intent(out) :: theta

    theta = datan2(z,dsqrt( x**2.d0 + y**2.d0 ))

    end

!!!*************************************************
!!!*						         	             *
!!!*            azimuth                       *
!!!*								                 *
!!!*************************************************

subroutine azimuth(x,y,phi)
    implicit none
    real*8, intent(in) :: x,y
    real*8, intent(out) :: phi

    phi = datan2(y,x)

    end

!!!*************************************************
!!!*						         	             *
!!!*                    dds                        *
!!!*								                 *
!!!*************************************************
    
subroutine  dds(ns,A,DA,ds)
    implicit none
    integer, intent(in) :: ns
    real*8, intent(in) :: A(1:ns),ds
    real*8, intent(out) :: DA(1:ns)
    integer j

    DA(1) =  -(3.d0*A(1) -4.d0*A(2) +A(3))/(2.d0*ds)

    do j=2,ns-1
        DA(j) = (A(j+1) -A(j-1))/(2.d0*ds)
    enddo

    DA(ns) = (3.d0*A(ns) -4.d0*A(ns-1) +A(ns-2))/(2.d0*ds)

    end

!!!*************************************************
!!!*						         	             *
!!!*                    d2ds2                        *
!!!*								                 *
!!!*************************************************

subroutine  d2ds2(ns,A,D2A,ds)
    implicit none
    integer, intent(in) :: ns
    real*8, intent(in) :: A(1:ns),ds
    real*8, intent(out) :: D2A(1:ns)
    real*8 ds2
    integer j

    ds2 = ds*ds
    
    D2A(1) =  (12.d0*a(1) -30.d0*a(2) +24.d0*a(3) -6.d0*a(4))/(6.d0*ds2)
    
    do j=2,ns-1
        D2A(j) = (a(j+1) - 2.d0*a(j) + a(j-1))/ds2
    enddo
    
    D2A(ns) = (12.d0*a(ns) -30.d0*a(ns-1) +24.d0*a(ns-2) -6.d0*a(ns-3))/(6.d0*ds2)
    
    end

!!!*************************************************
!!!*						         	             *
!!!*            Probability                        *
!!!*								                 *
!!!*************************************************

subroutine Probability(a,ds,Pr,ns,Pa)
    implicit none
    integer, intent(in) :: ns
    real*8, intent(in) :: ds,Pr(1:ns,1:ns,1:ns),a(1:ns,1:ns,1:ns)
    real*8, INTENT(OUT) :: Pa
    real*8 integrand
    integer i,j,k

    Pa = 0.d0
    do k =1,ns-1
        do j =1,ns-1
            do i =1,ns-1
                integrand = ( a(i+1,j+1,k+1) * Pr(i+1,j+1,k+1) + a(i,j+1,k+1) * Pr(i,j+1,k+1) &
                            + a(i+1,j,k+1) * Pr(i+1,j,k+1) + a(i+1,j+1,k) * Pr(i+1,j+1,k) &
                            + a(i+1,j,k) * Pr(i+1,j,k) + a(i,j+1,k) * Pr(i,j+1,k) &
                            + a(i,j,k+1) * Pr(i,j,k+1) + a(i,j,k) * Pr(i,j,k) ) / 8.d0
                Pa = Pa + integrand * ds * ds * ds
            enddo
        enddo
    enddo

    end

!!!*************************************************
!!!*						         	             *
!!!*            Probability                        *
!!!*								                 *
!!!*************************************************

subroutine probability_density_function(psi_r,psi_i,Pr,ns)
    implicit none
    integer, intent(in) :: ns
    real*8, intent(in) :: psi_r(1:ns,1:ns,1:ns),Psi_i(1:ns,1:ns,1:ns)
    real*8, INTENT(OUT) :: Pr(1:ns,1:ns,1:ns)
    integer i,j,k

    do k =1,ns
        do j =1,ns
            do i =1,ns
                Pr(i,j,k) = psi_r(i,j,k)**2.d0 + psi_i(i,j,k)**2.d0
            enddo
        enddo
    enddo

    end

!!!**********************************************************
!!!*						         	                    *
!!!*     Associated Legendre Polynomials Trig               *
!!!*								                        *
!!!**********************************************************

subroutine associated_Legendre_polynomials_trig(l,m,theta,P)
    implicit none 
    integer, intent(in) :: l,m
    real*8, intent(in) :: theta
    real*8, INTENT(OUT) :: P
    
    if(abs(m)>l) then
        print*, 'Inconsistency. Quit'
        stop
    endif
        
    if(l==0) then
        
        if(m==0) then
            P = 1.d0
        endif
        
    else if(l==1) then
        
        if(m==0) then
            P = dsin(theta)
        else if(m==1) then
            P = dcos(theta)
        else if(m==-1) then
            P = -dcos(theta)
        endif
            
    else if(l==2) then
        
        if(m==0) then
            P = (3.d0*dsin(theta)**2.d0 - 1.d0)/2.d0
        else if(m==1) then
            P = 3.d0*dsin(theta)*dcos(theta)
        else if(m==-1) THEN
            P = -3.d0*dsin(theta)*dcos(theta)
        else if(m==2) THEN
            P = 3.d0*dcos(theta)**2.d0
        else if(m==-2) THEN
            P = 3.d0*dcos(theta)**2.d0
        endif
            
    else if(l==3) then
        
        if(m==0) THEN
            P = (-3.d0*dsin(theta) +5.d0*dsin(theta)**3.d0)/2.d0
        else if(m==1) THEN
            P = dcos(theta)*(-3.d0 +15.d0*dsin(theta)**2.d0)/2.d0
        else if(m==-1) THEN
            P = -dcos(theta)*(-3.d0 +15.d0*dsin(theta)**2.d0)/2.d0
        else if(m==2) THEN
            P = 15.d0*dsin(theta)*dcos(theta)**2.d0
        else if(m==-2) THEN
            P = 15.d0*dsin(theta)*dcos(theta)**2.d0
        else if(m==3) THEN
            P = 15.d0*dcos(theta)**3.d0
        else if(m==-3) THEN
            P = -15.d0*dcos(theta)**3.d0
        endif
            
    else if(l==4) then
        
        if(m==0) THEN
            P = (3.d0 -30.d0*dsin(theta)**2.d0 +35.d0*dsin(theta)**4.d0)/8.d0
        else if(m==1) THEN
            P = ((-15.d0*dsin(theta) +35.d0*dsin(theta)**3.d0)/2.d0)*dcos(theta)
        else if(m==-1) THEN
            P = -((-15.d0*dsin(theta) +35.d0*dsin(theta)**3.d0)/2.d0)*dcos(theta)
        else if(m==2) THEN
            P = ((-15.d0 +105.d0*dsin(theta)**2.d0)/2.d0)*dcos(theta)**2.d0
        else if(m==-2) THEN
            P = ((-15.d0 +105.d0*dsin(theta)**2.d0)/2.d0)*dcos(theta)**2.d0
        else if(m==3) THEN
            P = 105.d0*dsin(theta)*dcos(theta)**3.d0
        else if(m==-3) THEN
            P = -105.d0*dsin(theta)*dcos(theta)**3.d0
        else if(m==4) THEN
            P = 105.d0*dcos(theta)**4.d0
        else if(m==-4) THEN
            P = 105.d0*dcos(theta)**4.d0
        endif
    
    endif
            
    end

!!!**********************************************************
!!!*						         	                    *
!!!*     Associated Legendre Polynomials Trig               *
!!!*								                        *
!!!**********************************************************

subroutine spherical_harmonics(l,m,ns,theta,phi,Y_r,Y_i)
    implicit none
    integer, intent(in) :: l,m,ns
    real*8, intent(in) :: theta(1:ns,1:ns,1:ns),phi(1:ns,1:ns,1:ns)
    real*8, intent(out) :: Y_r(1:ns,1:ns,1:ns),Y_i(1:ns,1:ns,1:ns)
    real*8 P(1:ns,1:ns,1:ns)
    integer i,j,k

    do k =1,ns
        do j =1,ns
            do i =1,ns
                call associated_Legendre_polynomials_trig(l,m,theta(i,j,k),P(i,j,k))
                Y_r(i,j,k) = P(i,j,k) * dcos(m * phi(i,j,k))
                Y_i(i,j,k) = P(i,j,k) * dsin(m * phi(i,j,k))
            enddo
        enddo
    enddo
    
    end

!!!**********************************************************
!!!*						         	                    *
!!!*        Generalized Laguerre Polynomials               *
!!!*								                        *
!!!**********************************************************

subroutine generalized_laguerre_polynomials(n,l,r,Rnl)
    implicit none
    integer, INTENT(IN) :: n,l
    real*8, INTENT(IN) :: r
    real*8, INTENT(OUT) :: Rnl
    
    if(n < l + 1) then
        print*, 'Inconsistency Quit'
        stop
    endif
        
    if(n==1) THEN
        
        if(l==0) THEN
            Rnl = dexp(-r)
        endif
            
    else if(n==2) THEN
        
        if(l==0) THEN
            Rnl = dexp(-r/2)*(-r+2)
        else if(l==1) THEN
            Rnl = r*dexp(-r/2)
        endif
            
    else if(n==3) THEN
        
        if(l==0) THEN
            Rnl = dexp(-r/3)*((2*r/3)**2 -6*(2*r/3) +6)/2
        else if(l==1) THEN
            Rnl = (2*r/3)*dexp(-r/3)*(4 - (2*r/3))
        else if(l==2) THEN
            Rnl = ((2*r/3)**2)*dexp(-r/3)
        endif
            
    else if(n==4) THEN
        
        if(l==0) THEN
            Rnl = dexp(-r/4)*(-(r/2)**3 +12*(r/2)**2 -36*(r/2) +24)/6
        else if(l==1) THEN
            Rnl = (r/2)*dexp(-r/4)*((r/2)**2 -10*(r/2) +20)/2
        else if(l==2) THEN
            Rnl = ((r/2)**2)*dexp(-r/4)*(6 - r/2)
        else if(l==3) THEN
            Rnl = ((r/2)**3)*dexp(-r/4)
        endif

    endif
            
    
    end

!!!**********************************************************
!!!*						         	                    *
!!!*                Hydrogen Wave Function               *
!!!*								                        *
!!!**********************************************************
    
subroutine hydrogen_wave_function(n,l,m,ns,ds,theta,phi,r,psi_r,psi_i,Pr)
    implicit none
    integer, intent(in) :: n,l,m,ns
    real*8, intent(in) :: theta(1:ns,1:ns,1:ns),phi(1:ns,1:ns,1:ns),r(1:ns,1:ns,1:ns),ds
    real*8, intent(out) :: psi_r(1:ns,1:ns,1:ns),psi_i(1:ns,1:ns,1:ns),Pr(1:ns,1:ns,1:ns)
    real*8 Y_r(1:ns,1:ns,1:ns),Y_i(1:ns,1:ns,1:ns),Rnl(1:ns,1:ns,1:ns)
    integer i,j,k
    real*8 Pa,ones(1:ns,1:ns,1:ns)
    
    call spherical_harmonics(l,m,ns,theta,phi,Y_r,Y_i)

    do k =1,ns
        do j =1,ns
            do i =1,ns
                call generalized_laguerre_polynomials(n,l,r(i,j,k),Rnl(i,j,k))
                psi_r(i,j,k) = Rnl(i,j,k) * Y_r(i,j,k)
                psi_i(i,j,k) = Rnl(i,j,k) * Y_i(i,j,k)
                ones(i,j,k) = 1.d0
            enddo
        enddo
    enddo

    call probability_density_function(psi_r,psi_i,Pr,ns)
    
    call Probability(ones,ds,Pr,ns,Pa)

    if(Pa==0) then
        Pa = 1.d0
    endif

    psi_r = psi_r/Pa
    psi_i = psi_i/Pa
    Pr = Pr/Pa
                
    print*, "Pa = ",Pa
    
    end
    
!!!*************************************************
!!!*						         	             *
!!!*       Calculate Mean Momentum                       *
!!!*								                 *
!!!*************************************************

subroutine mean_momentum(psi_r,psi_i,ds,ns,px_mean,py_mean,pz_mean)
    implicit none
    integer, intent(in) :: ns
    real*8, intent(in) :: ds,psi_r(1:ns,1:ns,1:ns),psi_i(1:ns,1:ns,1:ns)
    real*8, INTENT(OUT) :: px_mean,py_mean,pz_mean
    real*8 dpsi_rdx(1:ns,1:ns,1:ns),dpsi_rdy(1:ns,1:ns,1:ns),dpsi_rdz(1:ns,1:ns,1:ns)
    real*8 dpsi_idx(1:ns,1:ns,1:ns),dpsi_idy(1:ns,1:ns,1:ns),dpsi_idz(1:ns,1:ns,1:ns)
    integer i,j,k
    real*8 integrand_x,integrand_y,integrand_z

    do k =1,ns
        do j =1,ns
            call dds(ns,psi_r(:,j,k),dpsi_rdx(:,j,k),ds)
            call dds(ns,psi_i(:,j,k),dpsi_idx(:,j,k),ds)
        enddo
    enddo

    do k =1,ns
        do i =1,ns
            call dds(ns,psi_r(i,:,k),dpsi_rdy(i,:,k),ds)
            call dds(ns,psi_i(i,:,k),dpsi_idy(i,:,k),ds)
        enddo
    enddo

    do j =1,ns
        do i =1,ns
            call dds(ns,psi_r(i,j,:),dpsi_rdz(i,j,:),ds)
            call dds(ns,psi_i(i,j,:),dpsi_idz(i,j,:),ds)
        enddo
    enddo

    px_mean = 0.d0
    py_mean = 0.d0
    pz_mean = 0.d0
    do k =1,ns-1
        do j =1,ns-1
            do i =1,ns-1
                integrand_x = ( psi_r(i+1,j+1,k+1) * dpsi_idx(i+1,j+1,k+1) - psi_i(i+1,j+1,k+1) * dpsi_rdx(i+1,j+1,k+1) &
                            + psi_r(i,j+1,k+1) * dpsi_idx(i,j+1,k+1) - psi_i(i,j+1,k+1) * dpsi_rdx(i,j+1,k+1) &
                            + psi_r(i+1,j,k+1) * dpsi_idx(i+1,j,k+1) - psi_i(i+1,j,k+1) * dpsi_rdx(i+1,j,k+1) &
                            + psi_r(i+1,j+1,k) * dpsi_idx(i+1,j+1,k) - psi_i(i+1,j+1,k) * dpsi_rdx(i+1,j+1,k) &
                            + psi_r(i+1,j,k) * dpsi_idx(i+1,j,k) - psi_i(i+1,j,k) * dpsi_rdx(i+1,j,k) &
                            + psi_r(i,j+1,k) * dpsi_idx(i,j+1,k) - psi_i(i,j+1,k) * dpsi_rdx(i,j+1,k) &
                            + psi_r(i,j,k+1) * dpsi_idx(i,j,k+1) - psi_i(i,j,k+1) * dpsi_rdx(i,j,k+1) &
                            + psi_r(i,j,k) * dpsi_idx(i,j,k) - psi_i(i,j,k) * dpsi_rdx(i,j,k) ) / 8.d0
                px_mean = px_mean + integrand_x * ds * ds * ds
                integrand_y = ( psi_r(i+1,j+1,k+1) * dpsi_idy(i+1,j+1,k+1) - psi_i(i+1,j+1,k+1) * dpsi_rdy(i+1,j+1,k+1) &
                            + psi_r(i,j+1,k+1) * dpsi_idy(i,j+1,k+1) - psi_i(i,j+1,k+1) * dpsi_rdy(i,j+1,k+1) &
                            + psi_r(i+1,j,k+1) * dpsi_idy(i+1,j,k+1) - psi_i(i+1,j,k+1) * dpsi_rdy(i+1,j,k+1) &
                            + psi_r(i+1,j+1,k) * dpsi_idy(i+1,j+1,k) - psi_i(i+1,j+1,k) * dpsi_rdy(i+1,j+1,k) &
                            + psi_r(i+1,j,k) * dpsi_idy(i+1,j,k) - psi_i(i+1,j,k) * dpsi_rdy(i+1,j,k) &
                            + psi_r(i,j+1,k) * dpsi_idy(i,j+1,k) - psi_i(i,j+1,k) * dpsi_rdy(i,j+1,k) &
                            + psi_r(i,j,k+1) * dpsi_idy(i,j,k+1) - psi_i(i,j,k+1) * dpsi_rdy(i,j,k+1) &
                            + psi_r(i,j,k) * dpsi_idy(i,j,k) - psi_i(i,j,k) * dpsi_rdy(i,j,k) ) / 8.d0
                py_mean = py_mean + integrand_y * ds * ds * ds
                integrand_z = ( psi_r(i+1,j+1,k+1) * dpsi_idz(i+1,j+1,k+1) - psi_i(i+1,j+1,k+1) * dpsi_rdz(i+1,j+1,k+1) &
                            + psi_r(i,j+1,k+1) * dpsi_idz(i,j+1,k+1) - psi_i(i,j+1,k+1) * dpsi_rdz(i,j+1,k+1) &
                            + psi_r(i+1,j,k+1) * dpsi_idz(i+1,j,k+1) - psi_i(i+1,j,k+1) * dpsi_rdz(i+1,j,k+1) &
                            + psi_r(i+1,j+1,k) * dpsi_idz(i+1,j+1,k) - psi_i(i+1,j+1,k) * dpsi_rdz(i+1,j+1,k) &
                            + psi_r(i+1,j,k) * dpsi_idz(i+1,j,k) - psi_i(i+1,j,k) * dpsi_rdz(i+1,j,k) &
                            + psi_r(i,j+1,k) * dpsi_idz(i,j+1,k) - psi_i(i,j+1,k) * dpsi_rdz(i,j+1,k) &
                            + psi_r(i,j,k+1) * dpsi_idz(i,j,k+1) - psi_i(i,j,k+1) * dpsi_rdz(i,j,k+1) &
                            + psi_r(i,j,k) * dpsi_idz(i,j,k) - psi_i(i,j,k) * dpsi_rdz(i,j,k) ) / 8.d0
                pz_mean = pz_mean + integrand_z * ds * ds * ds
            enddo
        enddo
    enddo

    end   

!!!*************************************************
!!!*						         	             *
!!!*       Calculate Momentum Standard Deviation                       *
!!!*								                 *
!!!*************************************************

subroutine momentum_standard_deviation(psi_r,psi_i,Pr,ds,ns,px_mean,py_mean,pz_mean,px_std,py_std,pz_std)
    implicit none
    integer, intent(in) :: ns
    real*8, intent(in) :: ds,psi_r(1:ns,1:ns,1:ns),psi_i(1:ns,1:ns,1:ns)
    real*8, intent(in) :: px_mean,py_mean,pz_mean,Pr(1:ns,1:ns,1:ns)
    real*8, INTENT(OUT) :: px_std,py_std,pz_std
    real*8 d2psi_rdx2(1:ns,1:ns,1:ns),d2psi_rdy2(1:ns,1:ns,1:ns),d2psi_rdz2(1:ns,1:ns,1:ns)
    real*8 d2psi_idx2(1:ns,1:ns,1:ns),d2psi_idy2(1:ns,1:ns,1:ns),d2psi_idz2(1:ns,1:ns,1:ns)
    integer i,j,k
    real*8 integrand_x,integrand_y,integrand_z

    do k =1,ns
        do j =1,ns
            call d2ds2(ns,psi_r(:,j,k),d2psi_rdx2(:,j,k),ds)
            call d2ds2(ns,psi_i(:,j,k),d2psi_idx2(:,j,k),ds)
        enddo
    enddo

    do k =1,ns
        do i =1,ns
            call d2ds2(ns,psi_r(i,:,k),d2psi_rdy2(i,:,k),ds)
            call d2ds2(ns,psi_i(i,:,k),d2psi_idy2(i,:,k),ds)
        enddo
    enddo

    do j =1,ns
        do i =1,ns
            call d2ds2(ns,psi_r(i,j,:),d2psi_rdz2(i,j,:),ds)
            call d2ds2(ns,psi_i(i,j,:),d2psi_idz2(i,j,:),ds)
        enddo
    enddo

    px_std = 0.d0
    py_std = 0.d0
    pz_std = 0.d0
    do k =1,ns-1
        do j =1,ns-1
            do i =1,ns-1
                integrand_x = - ( psi_r(i+1,j+1,k+1) * d2psi_rdx2(i+1,j+1,k+1) + psi_i(i+1,j+1,k+1) * d2psi_idx2(i+1,j+1,k+1) &
                + px_mean * px_mean * Pr(i+1,j+1,k+1) & 
                + psi_r(i,j+1,k+1) * d2psi_rdx2(i,j+1,k+1) + psi_i(i,j+1,k+1) * d2psi_idx2(i,j+1,k+1) &
                + px_mean * px_mean * Pr(i,j+1,k+1) &
                + psi_r(i+1,j,k+1) * d2psi_rdx2(i+1,j,k+1) + psi_i(i+1,j,k+1) * d2psi_idx2(i+1,j,k+1) &
                + px_mean * px_mean * Pr(i+1,j,k+1) &
                + psi_r(i+1,j+1,k) * d2psi_rdx2(i+1,j+1,k) + psi_i(i+1,j+1,k) * d2psi_idx2(i+1,j+1,k) &
                + px_mean * px_mean * Pr(i+1,j+1,k) &
                + psi_r(i+1,j,k) * d2psi_rdx2(i+1,j,k) + psi_i(i+1,j,k) * d2psi_idx2(i+1,j,k) &
                + px_mean * px_mean * Pr(i+1,j,k) &
                + psi_r(i,j+1,k) * d2psi_rdx2(i,j+1,k) + psi_i(i,j+1,k) * d2psi_idx2(i,j+1,k) &
                + px_mean * px_mean * Pr(i,j+1,k) &
                + psi_r(i,j,k+1) * d2psi_rdx2(i,j,k+1) + psi_i(i,j,k+1) * d2psi_idx2(i,j,k+1) &
                + px_mean * px_mean * Pr(i,j,k+1) &
                + psi_r(i,j,k) * d2psi_rdx2(i,j,k) + psi_i(i,j,k) * d2psi_idx2(i,j,k) &
                + px_mean * px_mean * Pr(i,j,k) ) / 8.d0
                px_std = px_std  + integrand_x * ds * ds * ds
                integrand_y = - ( psi_r(i+1,j+1,k+1) * d2psi_rdy2(i+1,j+1,k+1) + psi_i(i+1,j+1,k+1) * d2psi_idy2(i+1,j+1,k+1) &
                + py_mean * py_mean * Pr(i+1,j+1,k+1) &
                + psi_r(i,j+1,k+1) * d2psi_rdy2(i,j+1,k+1) + psi_i(i,j+1,k+1) * d2psi_idy2(i,j+1,k+1) &
                + py_mean * py_mean * Pr(i,j+1,k+1) &
                + psi_r(i+1,j,k+1) * d2psi_rdy2(i+1,j,k+1) + psi_i(i+1,j,k+1) * d2psi_idy2(i+1,j,k+1) &
                + py_mean * py_mean * Pr(i+1,j,k+1) &
                + psi_r(i+1,j+1,k) * d2psi_rdy2(i+1,j+1,k) + psi_i(i+1,j+1,k) * d2psi_idy2(i+1,j+1,k) &
                + py_mean * py_mean * Pr(i+1,j+1,k) & 
                + psi_r(i+1,j,k) * d2psi_rdy2(i+1,j,k) + psi_i(i+1,j,k) * d2psi_idy2(i+1,j,k) &
                + py_mean * py_mean * Pr(i+1,j,k) &
                + psi_r(i,j+1,k) * d2psi_rdy2(i,j+1,k) + psi_i(i,j+1,k) * d2psi_idy2(i,j+1,k) &
                + py_mean * py_mean * Pr(i,j+1,k) &
                + psi_r(i,j,k+1) * d2psi_rdy2(i,j,k+1) + psi_i(i,j,k+1) * d2psi_idy2(i,j,k+1) &
                + py_mean * py_mean * Pr(i,j,k+1) &
                + psi_r(i,j,k) * d2psi_rdy2(i,j,k) + psi_i(i,j,k) * d2psi_idy2(i,j,k) &
                + py_mean * py_mean * Pr(i,j,k) ) / 8.d0
                py_std = py_std  + integrand_y * ds * ds * ds
                integrand_z = - ( psi_r(i+1,j+1,k+1) * d2psi_rdz2(i+1,j+1,k+1) + psi_i(i+1,j+1,k+1) * d2psi_idz2(i+1,j+1,k+1) &
                + pz_mean * pz_mean * Pr(i+1,j+1,k+1) &
                + psi_r(i,j+1,k+1) * d2psi_rdz2(i,j+1,k+1) + psi_i(i,j+1,k+1) * d2psi_idz2(i,j+1,k+1) &
                + pz_mean * pz_mean * Pr(i,j+1,k+1) &
                + psi_r(i+1,j,k+1) * d2psi_rdz2(i+1,j,k+1) + psi_i(i+1,j,k+1) * d2psi_idz2(i+1,j,k+1) &
                + pz_mean * pz_mean * Pr(i+1,j,k+1) &
                + psi_r(i+1,j+1,k) * d2psi_rdz2(i+1,j+1,k) + psi_i(i+1,j+1,k) * d2psi_idz2(i+1,j+1,k) &
                + pz_mean * pz_mean * Pr(i+1,j+1,k) &
                + psi_r(i+1,j,k) * d2psi_rdz2(i+1,j,k) + psi_i(i+1,j,k) * d2psi_idz2(i+1,j,k) &
                + pz_mean * pz_mean * Pr(i+1,j,k) &
                + psi_r(i,j+1,k) * d2psi_rdz2(i,j+1,k) + psi_i(i,j+1,k) * d2psi_idz2(i,j+1,k) &
                + pz_mean * pz_mean * Pr(i,j+1,k) &
                + psi_r(i,j,k+1) * d2psi_rdz2(i,j,k+1) + psi_i(i,j,k+1) * d2psi_idz2(i,j,k+1) &
                + pz_mean * pz_mean * Pr(i,j,k+1) &
                + psi_r(i,j,k) * d2psi_rdz2(i,j,k) + psi_i(i,j,k) * d2psi_idz2(i,j,k) &
                + pz_mean * pz_mean * Pr(i,j,k) ) / 8.d0
                pz_std = pz_std + integrand_z * ds * ds * ds
            enddo
        enddo
    enddo

    px_std = px_std**0.5d0
    py_std = py_std**0.5d0
    pz_std = pz_std**0.5d0

    end   