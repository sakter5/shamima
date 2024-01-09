I!!  coIde for bifurcation diagram for deterministic IP3R model
!!  10 th february,2022

implicit none
Integer N, dim,it, i, j, Nt
Parameter (dim=4,Nt=1000000)
Integer state
REAL time, dt,rate,rate1,t_open,t_close,x, fun(dim)
REAL ip, scale, kIPR, vs, Ks, Jleakin, Vrocc, Vsocc,r,rate_open,rate_close,ran2
REAL Ksocc, Vp, Kp, kleak, gamma1, gamma2, kdiff, np, ns,v_delta,Cp0, Cp
REAL k24, km24, Cer, q26, q62, phi, Lm42, L, H, Lm24, Lh24, q42a
REAL V42, k42, km42, mi42, hi42, q24a, V24, mi24, hi24, q42, q24
REAL open_po, Jipr, Jdiff, Jleak, Js, Jrocc, Jsocc, Jin, Jpm, Lh42
REAL Cac, Cb, Ct, D, m42, h42, m24, h24,Vh42,ah42
REAL C4, C2,O, C400,C200, O100,Cac00,h4200,Cer00,O00,maxca,minca

open(unit=12, file='bf_vs_wt_cer500_20sec_paraold_final1_08_30_22.dat')

	
	!!! declerations
	
	!idum = 1
!state = 0
!idum = 1
!state = 0

Cac=0.05d0 !Cac
!Cb=0.05d0 !Cb
!Ct=46.0d0 !Ct
h42=0.5d0 !h42
Cer=500.0d0
time=0.0d0
dt=0.0001d0

!Parameters
!scale=4010.0d0/(4010.0d0+10500.0d0)
!scale=1745.0d0/(1745.0d0+611.0d0)
!scale=975.70/(975.70+4550.0)
kIPR=1.0d0!/float(Num_IPR)*(1.0d0-scale)
ip=0.05d0
!Vs=0.5d0
Ks=0.26d0
Jleakin=0.003115d0
Vrocc=0.2d0
Vsocc=0.15d0
Ksocc=100d0
Vp=0.8d0
Kp=0.5d0
kleak=0.0032d0
gamma1=100.0d0
gamma2=10.0d0
kdiff=10.0d0
np=2.0d0
ns=1.75d0
!Cp0=120.0d0
Vh42 = 20.0d0
ah42 = 0.5d0
k24 =0.08d0
km24 = 8.05d0
Lm42=30.0d0
L=0.5d0
H=20.0d0
Lm24=30.0d0
Lh24=10.0d0
C4=1.0d0
C2=0.0d0
O=0.0d0
q26=611.0d0!WT 
q62=1745.0d0!WT
do vs = 0.0d0, 20.0d0, 0.001
!print*, ip
    minca=1.00d0
    maxca=0.00d0
   open(18,file="ca_bf_vs.dat")
do it = 1, Nt
time = time+dt
!print*,ip
!C2 = 1.0d0- (O + C4)
q42a=1.8d0*ip**2.0d0/(ip**2.0d0+0.34d0)
V42=310.0d0*ip**2.0d0/(ip**2.0d0+0.01d0)
k42=1.2d0+0.543d0*ip**3.0d0/(ip**3.0d0+64.0d0)
km42=1.0d0+25.0d0*ip**3.0d0/(ip**3.0d0+274.6d0)
q24a=1.0d0+5.0d0/(ip**2.0d0+0.25d0)
V24=62.0d0+880.0d0/(ip**2.0d0+4.0d0)
mi42=Cac**3.0d0/(Cac**3.0d0+k42**3.0d0) 
hi42=km42**3.0d0/(Cac**3.0d0+km42**3.0d0)
mi24 =Cac**3.0d0/(Cac**3.0d0+k24**3.0d0) 
hi24=km24**3.0d0/(Cac**3.0d0+km24**3.0d0) 
Lh42 = ah42
fun(4)=Lh42*(hi42-h42)
fun=0.5*(hi42-h42) !h42
h42=h42+dt*fun(4)
q42=q42a+V42*mi42*h42
q24=q24a+V24*(1.0d0-mi24*hi24)
C4 = C4+dt*(q24*C2-q42*C4)
C2 = C2+dt*(q62*O+q42*C4-q26*C2-q24*C2)
O = O +dt*(q26*C2-q62*O)
!print*,O
fun(1)=Jipr+Jleak-Js+Jin-Jpm !cac
fun(5) =5.405*(Js-Jleak-Jipr)
!fun(6) = beta*(IP3star-p)+J_delta
!fun(2)=gamma1*(Jipr-Jdiff) !cb
!fun(3)=Jin-Jpm !ct
!open_po=q42*(q26+q62)/(q24*q62+q42*(q26+q62))
!Jipr = 2*kIPR*open_po*(Cer-Cac)
Jipr =5*kIPR*O*(Cer-Cac)
!Jdiff=kdiff*(Cb-Cac)
Jleak=kleak*(Cer-Cac)
Js = vs*Cac**ns/(Ks**ns+Cac**ns)
Jrocc =Vrocc*iP
Jsocc=Vsocc*Ksocc**4.0d0/(Ksocc**4.0d0+Cer**4.0d0)
Jin =Jleakin+Jrocc+Jsocc
Jpm=Vp*Cac**np/(Kp**np+Cac**np)
Cac=Cac+dt*fun(1)
Cer=Cer+dt*fun(5)
!p=p+dt*(fun(6)
!Cb=Cb+dt*fun(2)
!Ct=Ct+dt*fun(3)
!!!!!
!print*, P,Cac,open_Po
write(18,*) Cac
end do
close(18)
open(18,file="ca_bf_vs.dat",status="unknown")
do j=1,Nt
if(j.gt.Nt/2)then
read(18,*) x
if(x.gt.maxca)then
maxca=x
end if
if(x.lt.minca)then
minca=x
end if
end if
end do
print*,maxca,minca
write(12,*) vs,maxca,minca	
close(18)

end do
	
close(12)
	
	
stop
end

