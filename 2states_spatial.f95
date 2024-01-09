Program main
implicit none
Integer N, dim,it, i, j, M,Num_IPR,Num_Open,idum, Nx,Ny,Nc
Parameter (N=25000000,dim=5, Num_IPR = 20, Nx =41,Nc=21, Ny = 41)
Integer state(Num_IPR)
REAL*8 :: time, dt,fun(dim), diffusion_Ca, step, space_grid
REAL*8:: P, scale, kIPR, Vs, Ks, Jleakin, Vrocc, Vsocc,r,rate_open,rate_close,ran2
REAL*8:: Ksocc, Vp, Kp, kleak, gamma1, gamma2, kdiff, np, ns, Cp0
REAL*8:: k24, km24, q26, q62, phi, Lm42, L, H, Lm24, Lh24, q42a
REAL*8 ::V42, k42, km42, mi42, hi42, q24a, V24, mi24, hi24, q42, q24
REAL*8 :: Po, Jipr, Jdiff, Jleak, Js, Jrocc, Jsocc, Jin, Jpm, Lh42
REAL*8 :: D, m42, h42, m24, h24,Vh42,ah42, Ca_total, ca_total1
REAL*8, DIMENSION(:,:), ALLOCATABLE :: Cac, Cb, Cer,ct, Cp, Cac1, Cb1,Ct1
 Allocate(Cac(0:Nx+1,0:Ny+1),Cb(Nx,Ny), Cer(Nx,Ny),ct(Nx,Ny),Cp(Nx,Ny))
 Allocate(Cac1(Nx,Ny),Cb1(Nc,Nc),Ct1(Nx,Ny))
open(unit=1, file='Ca_new250sec_IP3.dat')
!open(unit = 3, file ='Ca_1um_AD250sec_IP3.dat')
!open(unit = 4, file ='Ca_1um_AD250secchannel.dat')
!open(unit=2, file='snap250sec.dat')

idum = 1
state = 0
h42=0.5d0 !h42
time=0.0d0
M=100000

dt=0.00001d0 ! in the original one,dt=0.00002
Diffusion_Ca = 30.0d0
space_grid= 0.05d0
step  = Diffusion_Ca*dt/(space_grid**2)
!print*,step
!Parameters
P=0.02d0
scale=4010.0d0/(4010.0d0+10500.0d0)
!print*, scale
kIPR=0.05d0!/float(Num_IPR)*(1.0d0-scale)
Vs=10.0d0

Ks=0.26d0
Jleakin=0.03115d0
Vrocc=0.2d0
Vsocc=1.6d0
Ksocc=100d0
Vp=0.8d0
Kp=0.5d0
kleak=0.0032d0
gamma1=100.0d0
gamma2=10.0d0

!!!!!!!!!!!
do i = 1, Nx
  do j =1, Ny
Cac(i,j)=0.05d0 !Cac
Cb(i,j)=0.05d0 !Cb
Ct(i,j)=46.0d0 !Ct

Cac1(i,j)=0.05d0 !Cac
Cb1(i,j)=0.05d0 !Cb
Ct1(i,j)=46.0d0 !Ct
Cer(i,j) = gamma2*(Ct(i,j)-Cac(i,j)-Cb(i,j)/gamma1)
enddo
   enddo
!Cac(1:Nx,0) = Cac(1:Nx,1)
!  Cac(0,1:Ny)=Cac(1,1:Ny)
!  Cac(1:Nx,NY+1) = Cac(1:Nx,Ny)
!  Cac(Nx+1,1:Ny) = Cac(Nx,1:Ny)   



kdiff=10.0d0
np=2.0d0
ns=1.75d0
Cp0=120.0d0
Vh42 = 20.0d0
ah42 = 0.5d0

Do it=1,N
    time = time + dt
   Ca_total = 0
   ca_total1 = 0   
!  Cer = gamma2*(Ct-Cac-Cb/gamma1)   
!   q26=2927.0d0   !WT
    q26=10500.0d0
    q62=4010.0d0
    phi=q26/(q26+q62)
!   q62=1361.0d0   !WT
    !q62=620.7d0   !AD
    !phi=q26/(q26+q62)
    q42a=1.8d0*P**2.0d0/(P**2.0d0+0.34d0)         
    V42=110.0d0*P**2.0d0/(P**2.0d0+0.01d0)        
    k42=0.49d0+0.543d0*P**3.0d0/(P**3.0d0+64.0d0)
    km42=0.41d0+25.0d0*P**3.0d0/(P**3.0d0+274.6d0)
    q24a=1.0d0+5.0d0/(P**2.0d0+0.25d0)
    V24=62.0d0+880.0d0/(P**2.0d0+4.0d0)  

!   k24=0.2d0 !WT
    k24=0.35d0 
    km24=80.0d0   
    Lm42=100.0d0
    L=0.5d0
    H=20.0d0
    Lm24=100.0d0
    Lh24=40.0d0
    Num_Open = 0
    Do j=1,Num_IPR
        r=ran2(idum)
!   print*,r     
      if(state(j)==0)then
          Cp(Nc,Nc)=Cb(Nc,Nc)
           mi42=Cp(Nc,Nc)**3.0d0/(Cp(Nc,Nc)**3.0d0+k42**3.0d0) !use Cp or Cb here
           hi42=km42**3.0d0/(Cp(Nc,Nc)**3.0d0+km42**3.0d0) !use Cp or Cb here
           mi24=Cp(Nc,Nc)**3.0d0/(Cp(Nc,Nc)**3.0d0+k24**3.0d0) !use Cp or Cb here
           hi24=km24**2.0d0/(Cp(Nc,Nc)**2.0d0+km24**2.0d0) !use Cp or Cb here
           Lh42 = ah42
           fun(4)=Lh42*(hi42-h42) !h42
           h42=h42+dt*fun(4)	
           q42=q42a+V42*mi42*h42
           q24=q24a+V24*(1.0d0-mi24*hi24)
            rate_open = q42*dt  
!print*,rate_open 
       if(r .lt. rate_open)then
                state(j) = 1
          end if
          else if(state(j)==1)then
          Cp(Nc,Nc) = Cp0*(Cer(Nc,Nc)/100.0d0)
                        mi42=Cp(Nc,Nc)**3.0d0/(Cp(Nc,Nc)**3.0d0+k42**3.0d0) !use Cp or Cb here
                        hi42=km42**3.0d0/(Cp(Nc,Nc)**3.0d0+km42**3.0d0) !use Cp or Cb here
                        mi24=Cp(Nc,Nc)**3.0d0/(Cp(Nc,Nc)**3.0d0+k24**3.0d0) !use Cp or Cb here
                        hi24=km24**2.0d0/(Cp(Nc,Nc)**2.0d0+km24**2.0d0) !use Cp or Cb here
                        Lh42 = Vh42
                        fun(4)=Lh42*(hi42-h42) !h42
                        h42=h42+dt*fun(4)	
                        q42=q42a+V42*mi42*h42
                        q24=q24a+V24*(1.0d0-mi24*hi24)
                        rate_close = q24*q62/(q62+q26)*dt          
               if(r .lt. rate_close)then
                state(j) = 0
            end if
          end if
  
          if(state(j)==1)then
           Num_Open = Num_Open + 1
        end if
  End Do
!   print*,Num_open
         Po=float(Num_Open)/float(Num_IPR)
         Jipr = 1000*kIPR*Po*(Cer(Nc,Nc)-Cb(Nc,Nc))
         Jdiff=kdiff*(Cb(Nc,Nc)-Cac(Nc,Nc))

       do i = 1, Nx
         do j = 1, Ny
         Jleak=kleak*(Cer(i,j)-Cac(i,j))
         Js = Vs*Cac(i,j)**ns/(Ks**ns+Cac(i,j)**ns)
        Jrocc = Vrocc*P
        Jsocc=Vsocc*Ksocc**4.0d0/(Ksocc**4.0d0+Cer(i,j)**4.0d0)
        Jin = Jleakin+Jrocc+Jsocc
        Jpm=Vp*Cac(i,j)**np/(Kp**np+Cac(i,j)**np)
        fun(1)=Jdiff+Jleak-Js+Jin-Jpm !cac    
        fun(2)=gamma1*(Jipr-Jdiff) !cb
        fun(3)=Jin-Jpm !ct
        fun(5) = Jleak-Js+Jin-Jpm
!!!!   NO FLUX BRDY CONDITION IMPOSED
       
cac(1,j)=cac1(2,j)
cac(Nx,j)=cac1(Nx-1,j)

if((i.ge.2).and.(i.le.(Nx-1)))then
cac(i,1)=cac1(i,2)
cac(i,Ny)=cac1(i,Ny-1) 
end if	 
       
 if((i.eq.Nc).and.(j.eq.Nc))then
  Cac1(i,j)=Cac(i,j)+step*(Cac(i+1,j)+Cac(i-1,j)+Cac(i,j+1)+Cac(i,j-1)-4*Cac(i,j))+&
                 dt*fun(1)      
  Cb1(i,j)=Cb(i,j)+dt*fun(2)
  else 
  Cac1(i,j)=Cac(i,j)+step*(Cac(i+1,j)+Cac(i-1,j)+Cac(i,j+1)+Cac(i,j-1)-4*Cac(i,j))+&
                   dt*fun(5)
  endif    
  Ct1(i,j)=Ct(i,j)+dt*fun(3)
  ! endif
enddo 
enddo
!Cac(1:Nx,0) = Cac(1:Nx,1)
! Cac(0,1:Ny)=Cac(1,1:Ny)
! Cac(1:Nx,NY+1) = Cac(1:Nx,Ny)
! Cac(Nx+1,1:Ny) = Cac(Nx,1:Ny)
    

      do i =1, Nx
            do j = 1,Ny
            
        Cac(i,j)=cac1(i,j)
        Cb(i,j)=Cb1(i,j)
        Ct(i,j)=Ct1(i,j)
          
       Ca_total = Ca_total+Cac(i,j)

   !   if((i.ge.10).and.(i.le.30).and.(j.ge.10).and.(j.le.30))then   
   !   Ca_total1 = Ca_total1+Cac(i,j)
    ! endif 
       end do
           end do
!print*,time,Ca_total/float(Nx*Ny),Num_open 
 write(1,*)time,Ca_total/float(Nx*Ny),Cac(Nc,Nc),Num_open
!write(3,*) time, Ca_total1/float(20*20)    
!write(4,*) time,Num_open
!if (mod(it,200).eq.0) then
!write(2,"(10000(F14.8, '  '))") time
!  do i=1,Nx
 !write(2,"(10000(F14.8, '  '))")(Cac(i,j), j = 1,Ny)
 !enddo
 !endif


End Do 

!return
END program main

FUNCTION ran2(idum)
!
INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
double precision ran2,AM,EPS,RNMX
PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,      &
&     IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,          &
&     IQ2=52774,IR1=12211,IR2=3791,NTAB=32,              &
&     NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
INTEGER idum2,j,k,iv(NTAB),iy
SAVE iv,iy,idum2
DATA idum2/123456789/, iv/NTAB*0/, iy/0/
if (idum.le.0)then
idum=max(-idum,1)
idum2=idum
do  j=NTAB+8,1,-1
k=idum/IQ1
idum=IA1*(idum-k*IQ1)-k*IR1
if (idum.lt.0) idum=idum+IM1
if (j .le. NTAB) iv(j)=idum
end do
iy=iv(1)
end if
k=idum/IQ1
idum=IA1*(idum-k*IQ1)-k*IR1
if (idum.lt.0) idum=idum+IM1
k=idum2/IQ2
idum2=IA2*(idum2-k*IQ2)-k*IR2
if (idum2.lt.0) idum2=idum2+IM2
j=1+iy/NDIV
iy=iv(j)-idum2
iv(j)=idum
if(iy.lt.1)iy=iy+IMM1
ran2=min(AM*iy,RNMX)
END




