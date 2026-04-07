program test
  use krome_main
  use krome_user
  use krome_commons
  use krome_photo
! Variables & definitions
! n: Number of layers
! nsp: Number of chemical species
  implicit none
  integer,parameter::n=10,nsp=krome_nmols
  real*8, parameter::coef_33 = 0.0078974358974359,coef_34 = 0.0443076923076923,coef_36 = 0.00021025641025641
  integer::i,ierr,istep,j,k,z,idx,ii,tout
  real*8::du(n),dd(n+1),dl(n),diff(n+1),r,dx,Tgas(n+1)
  real*8::dt,t,tend,uall(n+1,nsp),x(nsp),tmax,xvar
  real*8::o(n+1),co(n+1),m(n+1),co2(n+1),so2(n+1),so2a(n+1),so2b(n+1),so2c(n+1),so2d(n+1),so(n+1),soa(n+1),sob(n+1),soc(n+1),sod(n+1)
  real*8::s(n+1),sa(n+1),sb(n+1),sc(n+1),sd(n+1),ocs(n+1),ocsa(n+1),ocsb(n+1),ocsc(n+1),ocsd(n+1)
  real*8::s2o2(n+1),s2o2a(n+1),s2o2b(n+1),s2o2c(n+1),s2o2d(n+1),s2o2e(n+1),s2o2f(n+1),s2o2g(n+1),s2o2h(n+1),s2o2i(n+1),s2o2j(n+1)
  real*8::s2o(n+1),s2oa(n+1),s2ob(n+1),s2oc(n+1),s2od(n+1),s2oe(n+1),s2of(n+1),s2og(n+1),s2oh(n+1),s2oi(n+1),s2oj(n+1)
  real*8::s2(n+1),s2a(n+1),s2b(n+1),s2c(n+1),s2d(n+1),s2e(n+1),s2f(n+1),s2g(n+1),s2h(n+1),s2i(n+1),s2j(n+1)
  real*8::s3(n+1),s3a(n+1),s3b(n+1),s3c(n+1),s3d(n+1),s3e(n+1),s3f(n+1),s3g(n+1),s3h(n+1),s3i(n+1),s3j(n+1),s3k(n+1),s3l(n+1),s3m(n+1),s3n(n+1),s3o(n+1),s3p(n+1),s3q(n+1),s3r(n+1),s3s(n+1),s3t(n+1)
  real*8::s4(n+1),s4a(n+1),s4b(n+1),s4c(n+1),s4d(n+1),s4e(n+1),s4f(n+1),s4g(n+1),s4h(n+1),s4i(n+1),s4j(n+1),s4k(n+1),s4l(n+1),s4m(n+1),s4n(n+1),s4o(n+1),s4p(n+1),s4q(n+1)&
  ,s4r(n+1),s4s(n+1),s4t(n+1),s4u(n+1),s4v(n+1),s4w(n+1),s4x(n+1),s4y(n+1),s4z(n+1),s4aa(n+1),s4bb(n+1),s4cc(n+1),s4dd(n+1),s4ee(n+1),s4ff(n+1),s4gg(n+1),s4hh(n+1),s4ii(n+1),s4jj(n+1),s4kk(n+1)
  real*8::s8(n+1),One32S(n+1),Two32S(n+1),Three32S(n+1),Four32S(n+1),Five32S(n+1),Six32S(n+1),Seven32S(n+1),Eight32S(n+1),One33S(n+1),Two33S(n+1),Three33S(n+1),Four33S(n+1),Five33S(n+1),Six33S(n+1),Seven33S(n+1),Eight33S(n+1)&
  ,One34S(n+1),Two34S(n+1),Three34S(n+1),Four34S(n+1),Five34S(n+1),Six34S(n+1),Seven34S(n+1),Eight34S(n+1),One36S(n+1),Two36S(n+1),Three36S(n+1),Four36S(n+1),Five36S(n+1),Six36S(n+1),Seven36S(n+1),Eight36S(n+1)
  real*8::singleta(n+1),singleta32(n+1),singleta33(n+1),singleta34(n+1),singleta36(n+1)
  real*8::singletb(n+1),singletb32(n+1),singletb33(n+1),singletb34(n+1),singletb36(n+1)
  real*8::triplet(n+1),triplet32(n+1),triplet33(n+1),triplet34(n+1),triplet36(n+1)
  real*8::so3(n+1),so3a(n+1),so3b(n+1),so3c(n+1),so3d(n+1),o2(n+1)
  real*8::Jflux(n+1),flux(n+1,krome_nPhotoBins),Xeflux(krome_nPhotoBins,n+1),Xelamp(krome_nPhotoBins)
  real*8::tau(krome_nPhotoBins,n+1),op(krome_nPhotoBins),tauAll(n+1,krome_nPhotoBins)
  real*8::Tflux(krome_nPhotoBins,n+1),datar(n+1,nsp),TRates(nphotoRea,n+1),photoRates(n+1,nphotoRea)
  real*8::so2_d33(n+1),so2_d34(n+1),so2_d36(n+1),so2_Delta33(n+1),so2_Delta36(n+1),so_d33(n+1),so_d34(n+1),so_d36(n+1),so_Delta33(n+1),so_Delta36(n+1)
  real*8::s2o2_d33(n+1),s2o2_d34(n+1),s2o2_d36(n+1),s2o2_Delta33(n+1),s2o2_Delta36(n+1),s2_d33(n+1),s2_d34(n+1),s2_d36(n+1),s2_Delta33(n+1),s2_Delta36(n+1)
  real*8::s2o_d33(n+1),s2o_d34(n+1),s2o_d36(n+1),s2o_Delta33(n+1),s2o_Delta36(n+1),s_d33(n+1),s_d34(n+1),s_d36(n+1),s_Delta33(n+1),s_Delta36(n+1)
  real*8::so3_d33(n+1),so3_d34(n+1),so3_d36(n+1),so3_Delta33(n+1),so3_Delta36(n+1),ocs_d33(n+1),ocs_d34(n+1),ocs_d36(n+1),ocs_Delta33(n+1),ocs_Delta36(n+1)
  real*8::singleta_d33(n+1),singleta_d34(n+1),singleta_d36(n+1),singleta_Delta33(n+1),singleta_Delta36(n+1),singletb_d33(n+1),singletb_d34(n+1),singletb_d36(n+1),singletb_Delta33(n+1),singletb_Delta36(n+1)
  real*8::triplet_d33(n+1),triplet_d34(n+1),triplet_d36(n+1),triplet_Delta33(n+1),triplet_Delta36(n+1),s3_d33(n+1),s3_d34(n+1),s3_d36(n+1),s3_Delta33(n+1),s3_Delta36(n+1)
  real*8::s4_d33(n+1),s4_d34(n+1),s4_d36(n+1),s4_Delta33(n+1),s4_Delta36(n+1)
  real*8::s8_d33(n+1),s8_d34(n+1),s8_d36(n+1),s8_Delta33(n+1),s8_Delta36(n+1)
!Paramaters for diffusion calculation
  real*8::rout(2),n1(n+1),r1(n+1),ntot,h(n+1)
  real*8::nlayer(n+1),p(n+1),dens(n),n_layer(n+1)

!init krome
  call krome_init()

!Spatial & Time variables
  dx = 44d-1 !grid size (cm) 44d0
  t = 0d0 !time (s)
  dt = 1d0 !time-step (s)  changed 17/10/10
  istep = 0 !integration step
  tmax = 3600d0 !total simulation time (s) 
  !diffusion (space)
  diff(:) = 44d-1 !1d8
!***********************************************************************************************

!Define and read initial concentrations (space,species), cm-3
  open(28,file="initial_chamberXe8Pa_10.txt",status="old")
  do i=1,n+1
	  read(28, *) (datar(i,z), z=1,6) 
    uall(i,krome_idx_32SO2)=datar(i,1)
	  uall(i,krome_idx_33SO2)=datar(i,2)
	  uall(i,krome_idx_34SO2)=datar(i,3)
	  uall(i,krome_idx_36SO2)=datar(i,4)
	  uall(i,krome_idx_CO)=datar(i,5)
	  uall(i,krome_idx_M)=datar(i,6)

	  print *,uall(i,krome_idx_32SO2)
	  print *,uall(i,krome_idx_33SO2)
	  print *,uall(i,krome_idx_34SO2)
	  print *,uall(i,krome_idx_36SO2)
	  print *,uall(i,krome_idx_CO)
	  print *,uall(i,krome_idx_M)
  end do	
  close(28)
  
  uall(:,krome_idx_32SO) =0d0
  uall(:,krome_idx_33SO) =0d0
  uall(:,krome_idx_34SO) =0d0
  uall(:,krome_idx_36SO) =0d0
  uall(:,krome_idx_32S) =0d0
  uall(:,krome_idx_33S) =0d0
  uall(:,krome_idx_34S) =0d0
  uall(:,krome_idx_36S) =0d0
  uall(:,krome_idx_O) =0d0
  uall(:,krome_idx_O2) =0d0
  uall(:,krome_idx_OC32S) =0d0
  uall(:,krome_idx_OC33S) =0d0
  uall(:,krome_idx_OC34S) =0d0
  uall(:,krome_idx_OC36S) =0d0
  uall(:,krome_idx_CO2) =0d0
  uall(:,krome_idx_32S32S) =0d0
  uall(:,krome_idx_32S33S) =0d0
  uall(:,krome_idx_32S34S) =0d0
  uall(:,krome_idx_32S36S) =0d0
  uall(:,krome_idx_33S33S) =0d0
  uall(:,krome_idx_33S34S) =0d0
  uall(:,krome_idx_33S36S) =0d0
  uall(:,krome_idx_34S34S) =0d0
  uall(:,krome_idx_34S36S) =0d0
  uall(:,krome_idx_36S36S) =0d0 

  uall(:,krome_idx_32SO3) =0d0
  uall(:,krome_idx_33SO3) =0d0
  uall(:,krome_idx_34SO3) =0d0
  uall(:,krome_idx_36SO3) =0d0
  
  uall(:,krome_idx_32SO2_1A2) =0d0
  uall(:,krome_idx_33SO2_1A2) =0d0
  uall(:,krome_idx_34SO2_1A2) =0d0
  uall(:,krome_idx_36SO2_1A2) =0d0
  uall(:,krome_idx_32SO2_1B1) =0d0
  uall(:,krome_idx_33SO2_1B1) =0d0
  uall(:,krome_idx_34SO2_1B1) =0d0
  uall(:,krome_idx_36SO2_1B1) =0d0

  uall(:,krome_idx_32SO2_3B1) =0d0
  uall(:,krome_idx_33SO2_3B1) =0d0
  uall(:,krome_idx_34SO2_3B1) =0d0
  uall(:,krome_idx_36SO2_3B1) =0d0
  
  uall(:,krome_idx_O32S32SO) =0d0
  uall(:,krome_idx_O32S33SO) =0d0
  uall(:,krome_idx_O32S34SO) =0d0
  uall(:,krome_idx_O32S36SO) =0d0
  uall(:,krome_idx_O33S33SO) =0d0
  uall(:,krome_idx_O33S34SO) =0d0
  uall(:,krome_idx_O33S36SO) =0d0
  uall(:,krome_idx_O34S34SO) =0d0
  uall(:,krome_idx_O34S36SO) =0d0
  uall(:,krome_idx_O36S36SO) =0d0 

  uall(:,krome_idx_32S32SO) =0d0
  uall(:,krome_idx_32S33SO) =0d0
  uall(:,krome_idx_32S34SO) =0d0
  uall(:,krome_idx_32S36SO) =0d0
  uall(:,krome_idx_33S33SO) =0d0
  uall(:,krome_idx_33S34SO) =0d0
  uall(:,krome_idx_33S36SO) =0d0
  uall(:,krome_idx_34S34SO) =0d0
  uall(:,krome_idx_34S36SO) =0d0
  uall(:,krome_idx_36S36SO) =0d0 

  uall(:,krome_idx_32S32S32S) =0d0
  uall(:,krome_idx_33S33S33S) =0d0
  uall(:,krome_idx_34S34S34S) =0d0
  uall(:,krome_idx_36S36S36S) =0d0
  uall(:,krome_idx_32S32S33S) =0d0
  uall(:,krome_idx_32S32S34S) =0d0
  uall(:,krome_idx_32S32S36S) =0d0
  uall(:,krome_idx_32S33S33S) =0d0
  uall(:,krome_idx_32S34S34S) =0d0
  uall(:,krome_idx_32S36S36S) =0d0
  uall(:,krome_idx_32S33S34S) =0d0
  uall(:,krome_idx_32S33S36S) =0d0
  uall(:,krome_idx_32S34S36S) =0d0
  uall(:,krome_idx_33S33S34S) =0d0
  uall(:,krome_idx_33S33S36S) =0d0
  uall(:,krome_idx_33S34S34S) =0d0
  uall(:,krome_idx_33S34S36S) =0d0
  uall(:,krome_idx_34S34S36S) =0d0
  uall(:,krome_idx_33S36S36S) =0d0
  uall(:,krome_idx_34S36S36S) =0d0
  
  uall(:,krome_idx_32S32S32S32S) =0d0
  uall(:,krome_idx_33S33S33S33S) =0d0
  uall(:,krome_idx_34S34S34S34S) =0d0
  uall(:,krome_idx_36S36S36S36S) =0d0
  uall(:,krome_idx_32S32S32S33S) =0d0
  uall(:,krome_idx_32S32S32S34S) =0d0
  uall(:,krome_idx_32S32S32S36S) =0d0
  uall(:,krome_idx_32S32S33S33S) =0d0
  uall(:,krome_idx_32S32S33S34S) =0d0
  uall(:,krome_idx_32S32S33S36S) =0d0
  uall(:,krome_idx_32S32S34S34S) =0d0
  uall(:,krome_idx_32S32S34S36S) =0d0
  uall(:,krome_idx_32S32S36S36S) =0d0
  uall(:,krome_idx_32S33S33S33S) =0d0
  uall(:,krome_idx_32S33S33S34S) =0d0
  uall(:,krome_idx_32S33S33S36S) =0d0
  uall(:,krome_idx_32S33S34S34S) =0d0
  uall(:,krome_idx_32S33S34S36S) =0d0
  uall(:,krome_idx_32S33S36S36S) =0d0
  uall(:,krome_idx_32S34S34S34S) =0d0
  uall(:,krome_idx_32S34S34S36S) =0d0
  uall(:,krome_idx_32S34S36S36S) =0d0
  uall(:,krome_idx_32S36S36S36S) =0d0
  uall(:,krome_idx_33S33S33S34S) =0d0
  uall(:,krome_idx_33S33S33S36S) =0d0
  uall(:,krome_idx_33S33S34S34S) =0d0
  uall(:,krome_idx_33S33S36S36S) =0d0
  uall(:,krome_idx_33S33S34S36S) =0d0
  uall(:,krome_idx_33S34S34S34S) =0d0
  uall(:,krome_idx_33S34S34S36S) =0d0
  uall(:,krome_idx_33S34S36S36S) =0d0
  uall(:,krome_idx_33S36S36S36S) =0d0
  uall(:,krome_idx_34S34S34S36S) =0d0
  uall(:,krome_idx_34S34S36S36S) =0d0
  uall(:,krome_idx_34S36S36S36S) =0d0
  
  uall(:,krome_idx_32S8) = 0d0
  uall(:,krome_idx_32S733S) = 0d0
  uall(:,krome_idx_32S734S) = 0d0
  uall(:,krome_idx_32S736S) = 0d0
  uall(:,krome_idx_32S633S2) = 0d0
  uall(:,krome_idx_32S633S34S) = 0d0
  uall(:,krome_idx_32S633S36S) = 0d0
  uall(:,krome_idx_32S533S3) = 0d0
  uall(:,krome_idx_32S634S2) = 0d0
  uall(:,krome_idx_32S634S36S) = 0d0
  uall(:,krome_idx_32S533S234S) = 0d0
  uall(:,krome_idx_32S533S34S2) = 0d0
  uall(:,krome_idx_32S534S3) = 0d0
  uall(:,krome_idx_32S636S2) = 0d0
  uall(:,krome_idx_32S533S236S) = 0d0
  uall(:,krome_idx_32S533S34S36S) = 0d0
  uall(:,krome_idx_32S533S36S2) = 0d0
  uall(:,krome_idx_32S534S236S) = 0d0
  uall(:,krome_idx_32S534S36S2) = 0d0
  uall(:,krome_idx_32S536S3) = 0d0
  uall(:,krome_idx_32S433S4) = 0d0
  uall(:,krome_idx_32S433S334S) = 0d0
  uall(:,krome_idx_32S433S336S) = 0d0
  uall(:,krome_idx_32S333S5) = 0d0
  uall(:,krome_idx_32S433S234S2) = 0d0
  uall(:,krome_idx_32S433S234S36S) = 0d0
  uall(:,krome_idx_32S433S34S3) = 0d0
  uall(:,krome_idx_32S333S434S) = 0d0
  uall(:,krome_idx_32S433S236S2) = 0d0
  uall(:,krome_idx_32S433S34S236S) = 0d0
  uall(:,krome_idx_32S433S34S36S2) = 0d0
  uall(:,krome_idx_32S433S36S3) = 0d0
  uall(:,krome_idx_32S333S436S) = 0d0
  uall(:,krome_idx_32S434S4) = 0d0
  uall(:,krome_idx_32S434S336S) = 0d0
  uall(:,krome_idx_32S333S334S2) = 0d0
  uall(:,krome_idx_32S333S234S3) = 0d0
  uall(:,krome_idx_32S333S34S4) = 0d0
  uall(:,krome_idx_32S334S5) = 0d0
  uall(:,krome_idx_32S434S236S2) = 0d0
  uall(:,krome_idx_32S434S36S3) = 0d0
  uall(:,krome_idx_32S333S334S36S) = 0d0
  uall(:,krome_idx_32S333S234S236S) = 0d0
  uall(:,krome_idx_32S333S34S336S) = 0d0
  uall(:,krome_idx_32S334S436S) = 0d0
  uall(:,krome_idx_32S436S4) = 0d0
  uall(:,krome_idx_32S333S336S2) = 0d0
  uall(:,krome_idx_32S333S234S36S2) = 0d0
  uall(:,krome_idx_32S333S236S3) = 0d0
  uall(:,krome_idx_32S333S34S236S2) = 0d0
  uall(:,krome_idx_32S333S34S36S3) = 0d0
  uall(:,krome_idx_32S333S36S4) = 0d0
  uall(:,krome_idx_32S334S336S2) = 0d0
  uall(:,krome_idx_32S334S236S3) = 0d0
  uall(:,krome_idx_32S334S36S4) = 0d0
  uall(:,krome_idx_32S336S5) = 0d0
  uall(:,krome_idx_32S233S6) = 0d0
  uall(:,krome_idx_32S233S534S) = 0d0
  uall(:,krome_idx_32S233S536S) = 0d0
  uall(:,krome_idx_32S33S7) = 0d0
  uall(:,krome_idx_32S233S434S2) = 0d0
  uall(:,krome_idx_32S233S434S36S) = 0d0
  uall(:,krome_idx_32S233S334S3) = 0d0
  uall(:,krome_idx_32S33S634S) = 0d0
  uall(:,krome_idx_32S233S436S2) = 0d0
  uall(:,krome_idx_32S233S334S236S) = 0d0
  uall(:,krome_idx_32S233S334S36S2) = 0d0
  uall(:,krome_idx_32S233S336S3) = 0d0
  uall(:,krome_idx_32S33S636S) = 0d0
  uall(:,krome_idx_32S233S234S4) = 0d0
  uall(:,krome_idx_32S233S234S336S) = 0d0
  uall(:,krome_idx_32S233S34S5) = 0d0
  uall(:,krome_idx_32S33S534S2) = 0d0
  uall(:,krome_idx_32S233S234S236S2) = 0d0
  uall(:,krome_idx_32S233S234S36S3) = 0d0
  uall(:,krome_idx_32S233S34S436S) = 0d0
  uall(:,krome_idx_32S33S534S36S) = 0d0
  uall(:,krome_idx_32S233S236S4) = 0d0
  uall(:,krome_idx_32S233S34S336S2) = 0d0
  uall(:,krome_idx_32S233S34S236S3) = 0d0
  uall(:,krome_idx_32S233S34S36S4) = 0d0
  uall(:,krome_idx_32S233S36S5) = 0d0
  uall(:,krome_idx_32S33S536S2) = 0d0
  uall(:,krome_idx_32S33S36S6) = 0d0
  uall(:,krome_idx_32S234S6) = 0d0
  uall(:,krome_idx_32S234S536S) = 0d0
  uall(:,krome_idx_32S33S434S3) = 0d0
  uall(:,krome_idx_32S33S334S4) = 0d0
  uall(:,krome_idx_32S33S234S5) = 0d0
  uall(:,krome_idx_32S33S34S6) = 0d0
  uall(:,krome_idx_32S34S7) = 0d0
  uall(:,krome_idx_32S234S436S2) = 0d0
  uall(:,krome_idx_32S234S336S3) = 0d0
  uall(:,krome_idx_32S33S434S236S) = 0d0
  uall(:,krome_idx_32S33S334S336S) = 0d0
  uall(:,krome_idx_32S33S234S436S) = 0d0
  uall(:,krome_idx_32S33S34S536S) = 0d0
  uall(:,krome_idx_32S34S636S) = 0d0
  uall(:,krome_idx_32S234S236S4) = 0d0
  uall(:,krome_idx_32S234S36S5) = 0d0
  uall(:,krome_idx_32S33S434S36S2) = 0d0
  uall(:,krome_idx_32S33S334S236S2) = 0d0
  uall(:,krome_idx_32S33S234S336S2) = 0d0
  uall(:,krome_idx_32S33S34S436S2) = 0d0
  uall(:,krome_idx_32S34S536S2) = 0d0
  uall(:,krome_idx_32S236S6) = 0d0
  uall(:,krome_idx_32S33S436S3) = 0d0
  uall(:,krome_idx_32S33S334S36S3) = 0d0
  uall(:,krome_idx_32S33S336S4) = 0d0
  uall(:,krome_idx_32S33S234S236S3) = 0d0
  uall(:,krome_idx_32S33S234S36S4) = 0d0
  uall(:,krome_idx_32S33S236S5) = 0d0
  uall(:,krome_idx_32S33S34S336S3) = 0d0
  uall(:,krome_idx_32S33S34S236S4) = 0d0
  uall(:,krome_idx_32S33S34S36S5) = 0d0
  uall(:,krome_idx_32S34S436S3) = 0d0
  uall(:,krome_idx_32S34S336S4) = 0d0
  uall(:,krome_idx_32S34S236S5) = 0d0
  uall(:,krome_idx_32S34S36S6) = 0d0
  uall(:,krome_idx_32S36S7) = 0d0
  uall(:,krome_idx_33S8) = 0d0
  uall(:,krome_idx_33S734S) = 0d0
  uall(:,krome_idx_33S736S) = 0d0
  uall(:,krome_idx_33S634S2) = 0d0
  uall(:,krome_idx_33S634S36S) = 0d0
  uall(:,krome_idx_33S534S3) = 0d0
  uall(:,krome_idx_33S636S2) = 0d0
  uall(:,krome_idx_33S534S236S) = 0d0
  uall(:,krome_idx_33S534S36S2) = 0d0
  uall(:,krome_idx_33S536S3) = 0d0
  uall(:,krome_idx_33S434S4) = 0d0
  uall(:,krome_idx_33S434S336S) = 0d0
  uall(:,krome_idx_33S334S5) = 0d0
  uall(:,krome_idx_33S434S236S2) = 0d0
  uall(:,krome_idx_33S434S36S3) = 0d0
  uall(:,krome_idx_33S334S436S) = 0d0
  uall(:,krome_idx_33S436S4) = 0d0
  uall(:,krome_idx_33S334S336S2) = 0d0
  uall(:,krome_idx_33S334S236S3) = 0d0
  uall(:,krome_idx_33S334S36S4) = 0d0
  uall(:,krome_idx_33S336S5) = 0d0
  uall(:,krome_idx_33S234S6) = 0d0
  uall(:,krome_idx_33S234S536S) = 0d0
  uall(:,krome_idx_33S34S7) = 0d0
  uall(:,krome_idx_33S234S436S2) = 0d0
  uall(:,krome_idx_33S234S336S3) = 0d0
  uall(:,krome_idx_33S34S636S) = 0d0
  uall(:,krome_idx_33S234S236S4) = 0d0
  uall(:,krome_idx_33S234S36S5) = 0d0
  uall(:,krome_idx_33S34S536S2) = 0d0
  uall(:,krome_idx_33S236S6) = 0d0
  uall(:,krome_idx_33S34S436S3) = 0d0
  uall(:,krome_idx_33S34S336S4) = 0d0
  uall(:,krome_idx_33S34S236S5) = 0d0
  uall(:,krome_idx_33S34S36S6) = 0d0
  uall(:,krome_idx_33S36S7) = 0d0
  uall(:,krome_idx_34S8) = 0d0
  uall(:,krome_idx_34S736S) = 0d0
  uall(:,krome_idx_34S636S2) = 0d0
  uall(:,krome_idx_34S536S3) = 0d0
  uall(:,krome_idx_34S436S4) = 0d0
  uall(:,krome_idx_34S336S5) = 0d0
  uall(:,krome_idx_34S236S6) = 0d0
  uall(:,krome_idx_34S36S7) = 0d0
  uall(:,krome_idx_36S8) = 0d0
!***********************************************************************************************

!read initial layers data
open(34,file="layers_chamberXe8Pa_10.txt",status="old")
do i=1,n+1
   read(34,*) nlayer(i),h(i),p(i),Tgas(i)
   h(i)=h(i)  ! Height in cm
   p(i) = 9.8692326671601d-6 * p(i) !Pa->atm
   Tgas(i)= Tgas(i) +273.15  !temperature in K
   print *,nlayer(i),h(i),p(i),Tgas(i)
end do
!****************************************************************************************
!Prepare initial conditions for photochemistry
call krome_load_photoBin_file("Xelamp_400_Patrick.txt")
Tflux(:,1) = krome_get_photoBinJ() 
 
open(71,file="Initial-fluxesPa.txt",status='old',action='write',form='formatted',position="append")
write (71,*) "Xelamp flux"
do j=1,krome_nPhotoBins
	write (71,*) j, (Tflux(j,i), i=1,n+1)
end do

print *,'Step-3' 
!STOP  

!****************************************************************************************

!Here the main computation starts. 
!Begins by starting the time loop.
do
!  istep = istep + 1 !increse integration step  
  
	call krome_photoBin_restore()
  Tflux(:,1) = krome_get_photoBinJ()
	  
	tauAll(:,:) = 0d0
	
  do i=2,n+1
  
    x(:)= uall(i,:)
	  tauAll(i,:) = tauAll(i-1,:)+krome_get_opacity_size(uall(i,:),Tgas(i),dx)
    tauCell(:) = tauAll(i,:)
    
    if(i==6) then ! delete this line 
      if(mod(istep,3600)==0) then ! if(t==tmax) then  xvar = real(i-1,8)
        call krome_explore_flux(x(:), Tgas(i), 80, t) !krome_explore_flux(x(:), Tgas(i), 80, xvar)
      end if
    end if ! delete this line

   
! reaction here	  
    call krome(x(:),Tgas(i),dt)
    uall(i,:) = x(:)
	  photoRates(i,:) = photoBinRates(:)
	  
	 
    open(29,file="initial_chamberXe8Pa_10.txt",status="old")
	  read(29, *)
	  uall(i,krome_idx_M)=datar(i,6)
    close(29)

  end do
   
!****************************************************************************************
!Diffusion 09/25   
!comes from Ficks law: dc/dt = D d2c /dx2
  do j=1,nsp
    ntot = sum(uall(:,j))
    n_layer(:) = 0d0
	  n_layer(1) = 0d0
	  n_layer(2) = uall(2,j) + 1.4 * (uall(3,j) - uall(2,j))*dt / (dx**2)
      
    do i = 3,n 
		  n_layer(i) = uall(i,j) + 1.4 * (uall(i + 1,j) - 2 * uall(i,j) + uall(i - 1,j))*dt / (dx**2)
    end do

      n_layer(n+1) = uall(n+1,j) + 1.4 * (uall(n,j) - uall(n+1,j))*dt / (dx**2)
	   
      uall(:,j) = n_layer(:)
  end do
!****************************************************************************************
   
  tout = tmax/10
  if(mod(istep,tout)==0) then
   

   open(70,file="rates-evolutionPa.txt",status='old',action='write',form='formatted',position="append")
    write (70,*) "istep", istep, "time", t
    do i=1,n+1
      write (70,14) t, i-1, photoRates(i,:)
      14 format(f10.2,I5,999E18.4e5)
    end do
	 
	 
    so2(:) = uall(:,krome_idx_32SO2)+uall(:,krome_idx_33SO2)+uall(:,krome_idx_34SO2)+uall(:,krome_idx_36SO2)
	 so2a(:) = uall(:,krome_idx_32SO2)
	 so2b(:) = uall(:,krome_idx_33SO2)
	 so2c(:) = uall(:,krome_idx_34SO2)
	 so2d(:) = uall(:,krome_idx_36SO2)
    so(:) = uall(:,krome_idx_32SO)+uall(:,krome_idx_33SO)+uall(:,krome_idx_34SO)+uall(:,krome_idx_36SO)
	 soa(:) = uall(:,krome_idx_32SO)
	 sob(:) = uall(:,krome_idx_33SO)
	 soc(:) = uall(:,krome_idx_34SO)
	 sod(:) = uall(:,krome_idx_36SO)
	 s(:) = uall(:,krome_idx_32S)+uall(:,krome_idx_33S)+uall(:,krome_idx_34S)+uall(:,krome_idx_36S)
	 sa(:) = uall(:,krome_idx_32S)
	 sb(:) = uall(:,krome_idx_33S)
	 sc(:) = uall(:,krome_idx_34S)
	 sd(:) = uall(:,krome_idx_36S)
	 ocs(:) = uall(:,krome_idx_OC32S)+uall(:,krome_idx_OC33S)+uall(:,krome_idx_OC34S)+uall(:,krome_idx_OC36S)
	 ocsa(:) = uall(:,krome_idx_OC32S)
	 ocsb(:) = uall(:,krome_idx_OC33S)
	 ocsc(:) = uall(:,krome_idx_OC34S)
	 ocsd(:) = uall(:,krome_idx_OC36S)
    o(:) = uall(:,krome_idx_O)
	 o2(:) = uall(:,krome_idx_O2)
	 co(:) = uall(:,krome_idx_CO)
	 co2(:) = uall(:,krome_idx_CO2)
	 m(:) = uall(:,krome_idx_M)
	 
	 s2(:) = uall(:,krome_idx_32S32S)+uall(:,krome_idx_32S33S)+uall(:,krome_idx_32S34S)+uall(:,krome_idx_32S36S)&
   +uall(:,krome_idx_33S33S)+uall(:,krome_idx_33S34S)+uall(:,krome_idx_33S36S)+uall(:,krome_idx_34S34S)+uall(:,krome_idx_34S36S)+uall(:,krome_idx_36S36S)
	 s2a(:) = uall(:,krome_idx_32S32S)
	 s2b(:) = uall(:,krome_idx_32S33S)
	 s2c(:) = uall(:,krome_idx_32S34S)
	 s2d(:) = uall(:,krome_idx_32S36S)
   s2e(:) = uall(:,krome_idx_33S33S)
   s2f(:) = uall(:,krome_idx_33S34S)
   s2g(:) = uall(:,krome_idx_33S36S)
   s2h(:) = uall(:,krome_idx_34S34S)
   s2i(:) = uall(:,krome_idx_34S36S)
   s2j(:) = uall(:,krome_idx_36S36S)
	 
	 so3(:) = uall(:,krome_idx_32SO3)+uall(:,krome_idx_33SO3)+uall(:,krome_idx_34SO3)+uall(:,krome_idx_36SO3)
	 so3a(:) = uall(:,krome_idx_32SO3)
	 so3b(:) = uall(:,krome_idx_33SO3)
	 so3c(:) = uall(:,krome_idx_34SO3)
	 so3d(:) = uall(:,krome_idx_36SO3)
	 
   s2o2(:) = uall(:,krome_idx_O32S32SO)+uall(:,krome_idx_O32S33SO)+uall(:,krome_idx_O32S34SO)+uall(:,krome_idx_O32S36SO)&
   +uall(:,krome_idx_O33S33SO)+uall(:,krome_idx_O33S34SO)+uall(:,krome_idx_O33S36SO)+uall(:,krome_idx_O34S34SO)+uall(:,krome_idx_O34S36SO)+uall(:,krome_idx_O36S36SO)
	 s2o2a(:) = uall(:,krome_idx_O32S32SO)
	 s2o2b(:) = uall(:,krome_idx_O32S33SO)
	 s2o2c(:) = uall(:,krome_idx_O32S34SO)
	 s2o2d(:) = uall(:,krome_idx_O32S36SO)
   s2o2e(:) = uall(:,krome_idx_O33S33SO)
   s2o2f(:) = uall(:,krome_idx_O33S34SO)
   s2o2g(:) = uall(:,krome_idx_O33S36SO)
   s2o2h(:) = uall(:,krome_idx_O34S34SO)
   s2o2i(:) = uall(:,krome_idx_O34S36SO)
   s2o2j(:) = uall(:,krome_idx_O36S36SO)

   
   s2o(:) = uall(:,krome_idx_32S32SO)+uall(:,krome_idx_32S33SO)+uall(:,krome_idx_32S34SO)+uall(:,krome_idx_32S36SO)&
   +uall(:,krome_idx_33S33SO)+uall(:,krome_idx_33S34SO)+uall(:,krome_idx_33S36SO)+uall(:,krome_idx_34S34SO)+uall(:,krome_idx_34S36SO)+uall(:,krome_idx_36S36SO)
	 s2oa(:) = uall(:,krome_idx_32S32SO)
	 s2ob(:) = uall(:,krome_idx_32S33SO)
	 s2oc(:) = uall(:,krome_idx_32S34SO)
	 s2od(:) = uall(:,krome_idx_32S36SO)
   s2oe(:) = uall(:,krome_idx_33S33SO)
   s2of(:) = uall(:,krome_idx_33S34SO)
   s2og(:) = uall(:,krome_idx_33S36SO)
   s2oh(:) = uall(:,krome_idx_34S34SO)
   s2oi(:) = uall(:,krome_idx_34S36SO)
   s2oj(:) = uall(:,krome_idx_36S36SO)

	 singleta(:) = uall(:,krome_idx_32SO2_1A2)+uall(:,krome_idx_33SO2_1A2)+uall(:,krome_idx_34SO2_1A2)+uall(:,krome_idx_36SO2_1A2)
	 singleta32(:) = uall(:,krome_idx_32SO2_1A2)
	 singleta33(:) = uall(:,krome_idx_33SO2_1A2)
	 singleta34(:) = uall(:,krome_idx_34SO2_1A2)
	 singleta36(:) = uall(:,krome_idx_36SO2_1A2)
   singletb(:) = uall(:,krome_idx_32SO2_1B1)+uall(:,krome_idx_33SO2_1B1)+uall(:,krome_idx_34SO2_1B1)+uall(:,krome_idx_36SO2_1B1)
	 singletb32(:) = uall(:,krome_idx_32SO2_1B1)
	 singletb33(:) = uall(:,krome_idx_33SO2_1B1)
	 singletb34(:) = uall(:,krome_idx_34SO2_1B1)
	 singletb36(:) = uall(:,krome_idx_36SO2_1B1)

	 triplet(:) = uall(:,krome_idx_32SO2_3B1)+uall(:,krome_idx_33SO2_3B1)+uall(:,krome_idx_34SO2_3B1)+uall(:,krome_idx_36SO2_3B1)
	 triplet32(:) = uall(:,krome_idx_32SO2_3B1)
	 triplet33(:) = uall(:,krome_idx_33SO2_3B1)
	 triplet34(:) = uall(:,krome_idx_34SO2_3B1)
	 triplet36(:) = uall(:,krome_idx_36SO2_3B1)
	 
	 s3(:) = uall(:,krome_idx_32S32S32S)+uall(:,krome_idx_33S33S33S)+uall(:,krome_idx_34S34S34S)+uall(:,krome_idx_36S36S36S)&
   +uall(:,krome_idx_32S32S33S)+uall(:,krome_idx_32S32S34S)+uall(:,krome_idx_32S32S36S)+uall(:,krome_idx_32S33S33S)&
   +uall(:,krome_idx_32S34S34S)+uall(:,krome_idx_32S36S36S)+uall(:,krome_idx_32S33S34S)+uall(:,krome_idx_32S33S36S)&
   +uall(:,krome_idx_32S34S36S)+uall(:,krome_idx_33S33S34S)+uall(:,krome_idx_33S33S36S)+uall(:,krome_idx_33S34S36S)&
   +uall(:,krome_idx_34S34S36S)+uall(:,krome_idx_33S36S36S)+uall(:,krome_idx_34S36S36S)+uall(:,krome_idx_33S34S34S)
	 s3a(:) = uall(:,krome_idx_32S32S32S)
	 s3b(:) = uall(:,krome_idx_33S33S33S)
	 s3c(:) = uall(:,krome_idx_34S34S34S)
	 s3d(:) = uall(:,krome_idx_36S36S36S)
   s3e(:) = uall(:,krome_idx_32S32S33S)
	 s3f(:) = uall(:,krome_idx_32S32S34S)
	 s3g(:) = uall(:,krome_idx_32S32S36S)
	 s3h(:) = uall(:,krome_idx_32S33S33S)
	 s3i(:) = uall(:,krome_idx_32S34S34S)
	 s3j(:) = uall(:,krome_idx_32S36S36S)
	 s3k(:) = uall(:,krome_idx_32S33S34S)
	 s3l(:) = uall(:,krome_idx_32S33S36S)
	 s3m(:) = uall(:,krome_idx_32S34S36S)
	 s3n(:) = uall(:,krome_idx_33S33S34S)
	 s3o(:) = uall(:,krome_idx_33S33S36S)
	 s3p(:) = uall(:,krome_idx_33S34S34S)
	 s3q(:) = uall(:,krome_idx_33S34S36S)
	 s3r(:) = uall(:,krome_idx_34S34S36S)
	 s3s(:) = uall(:,krome_idx_33S36S36S)
	 s3t(:) = uall(:,krome_idx_34S36S36S)

   s4(:) = uall(:,krome_idx_32S32S32S32S)+uall(:,krome_idx_33S33S33S33S)+uall(:,krome_idx_34S34S34S34S)+uall(:,krome_idx_36S36S36S36S)&
   +uall(:,krome_idx_32S32S32S33S)+uall(:,krome_idx_32S32S32S34S)+uall(:,krome_idx_32S32S32S36S)+uall(:,krome_idx_32S32S33S33S)&
   +uall(:,krome_idx_32S32S33S34S)+uall(:,krome_idx_32S32S33S36S)+uall(:,krome_idx_32S32S34S34S)+uall(:,krome_idx_32S32S34S36S)&
   +uall(:,krome_idx_32S32S36S36S)+uall(:,krome_idx_32S33S33S33S)+uall(:,krome_idx_32S33S33S34S)+uall(:,krome_idx_32S33S33S36S)&
   +uall(:,krome_idx_32S33S34S34S)+uall(:,krome_idx_32S33S34S36S)+uall(:,krome_idx_32S33S36S36S)+uall(:,krome_idx_32S34S34S34S)&
   +uall(:,krome_idx_32S34S34S36S)+uall(:,krome_idx_32S34S36S36S)+uall(:,krome_idx_32S36S36S36S)+uall(:,krome_idx_33S33S33S34S)&
   +uall(:,krome_idx_33S33S33S36S)+uall(:,krome_idx_33S33S34S34S)+uall(:,krome_idx_33S33S36S36S)+uall(:,krome_idx_33S33S34S36S)&
   +uall(:,krome_idx_33S34S34S34S)+uall(:,krome_idx_33S34S34S36S)+uall(:,krome_idx_33S34S36S36S)+uall(:,krome_idx_33S36S36S36S)&
   +uall(:,krome_idx_34S34S34S36S)+uall(:,krome_idx_34S34S36S36S)+uall(:,krome_idx_34S36S36S36S) 
   s4a(:) = uall(:,krome_idx_32S32S32S32S)
	 s4b(:) = uall(:,krome_idx_33S33S33S33S)
	 s4c(:) = uall(:,krome_idx_34S34S34S34S)
	 s4d(:) = uall(:,krome_idx_36S36S36S36S)
   s4e(:) = uall(:,krome_idx_32S32S32S33S)
	 s4f(:) = uall(:,krome_idx_32S32S32S34S)
	 s4g(:) = uall(:,krome_idx_32S32S32S36S)
   s4h(:) = uall(:,krome_idx_32S32S33S33S)
   s4i(:) = uall(:,krome_idx_32S32S33S34S)
	 s4j(:) = uall(:,krome_idx_32S32S33S36S)
	 s4k(:) = uall(:,krome_idx_32S32S34S34S)
   s4l(:) = uall(:,krome_idx_32S32S34S36S)
   s4m(:) = uall(:,krome_idx_32S32S36S36S)
	 s4n(:) = uall(:,krome_idx_32S33S33S33S)
	 s4o(:) = uall(:,krome_idx_32S33S33S34S)
   s4p(:) = uall(:,krome_idx_32S33S33S36S)
   s4q(:) = uall(:,krome_idx_32S33S34S34S)
	 s4r(:) = uall(:,krome_idx_32S33S34S36S)
	 s4s(:) = uall(:,krome_idx_32S33S36S36S)
   s4t(:) = uall(:,krome_idx_32S34S34S34S)
	 s4u(:) = uall(:,krome_idx_32S34S34S36S)
	 s4v(:) = uall(:,krome_idx_32S34S36S36S)
	 s4w(:) = uall(:,krome_idx_32S36S36S36S)
   s4x(:) = uall(:,krome_idx_33S33S33S34S)
   s4y(:) = uall(:,krome_idx_33S33S33S36S)
	 s4z(:) = uall(:,krome_idx_33S33S34S34S)
   s4aa(:) = uall(:,krome_idx_33S33S36S36S)
   s4bb(:) = uall(:,krome_idx_33S33S34S36S) 
   s4cc(:) = uall(:,krome_idx_33S34S34S34S)
	 s4dd(:) = uall(:,krome_idx_33S34S34S36S)
	 s4ee(:) = uall(:,krome_idx_33S34S36S36S)
   s4ff(:) = uall(:,krome_idx_33S36S36S36S)
   s4gg(:) = uall(:,krome_idx_34S34S34S36S)
	 s4hh(:) = uall(:,krome_idx_34S34S36S36S)
   s4ii(:) = uall(:,krome_idx_34S36S36S36S) 
  
	 s8(:) = uall(:,krome_idx_32S8) +  uall(:,krome_idx_32S733S) + uall(:,krome_idx_32S734S) +  uall(:,krome_idx_32S736S)&
  +uall(:,krome_idx_32S633S2) +  uall(:,krome_idx_32S633S34S) + uall(:,krome_idx_32S634S2) +  uall(:,krome_idx_32S634S36S)&
  +uall(:,krome_idx_32S636S2) +  uall(:,krome_idx_32S533S3) + uall(:,krome_idx_32S533S234S) +  uall(:,krome_idx_32S533S236S)&
  +uall(:,krome_idx_32S533S34S2) +  uall(:,krome_idx_32S533S34S36S) + uall(:,krome_idx_32S533S36S2) +  uall(:,krome_idx_32S534S3)&
  +uall(:,krome_idx_32S534S236S) +  uall(:,krome_idx_32S534S36S2) + uall(:,krome_idx_32S536S3) +  uall(:,krome_idx_32S433S4)&
  +uall(:,krome_idx_32S433S334S) +  uall(:,krome_idx_32S433S336S) + uall(:,krome_idx_32S433S234S2) +  uall(:,krome_idx_32S433S234S36S)&
  +uall(:,krome_idx_32S433S236S2) +  uall(:,krome_idx_32S433S34S3) + uall(:,krome_idx_32S433S34S236S) +  uall(:,krome_idx_32S433S34S36S2)&
  +uall(:,krome_idx_32S433S36S3) +  uall(:,krome_idx_32S434S4) + uall(:,krome_idx_32S434S336S) +  uall(:,krome_idx_32S434S236S2)&
  +uall(:,krome_idx_32S434S36S3) +  uall(:,krome_idx_32S436S4) + uall(:,krome_idx_32S333S5) +  uall(:,krome_idx_32S333S434S)&
  +uall(:,krome_idx_32S333S436S) +  uall(:,krome_idx_32S333S334S2) + uall(:,krome_idx_32S333S334S36S) +  uall(:,krome_idx_32S333S336S2)&
  +uall(:,krome_idx_32S333S234S3) +  uall(:,krome_idx_32S333S234S236S) +uall(:,krome_idx_32S333S234S36S2) +  uall(:,krome_idx_32S333S236S3)&
  +uall(:,krome_idx_32S333S34S4) +  uall(:,krome_idx_32S333S34S336S) + uall(:,krome_idx_32S333S34S236S2) +  uall(:,krome_idx_32S333S34S36S3)&
  +uall(:,krome_idx_32S333S36S4) +  uall(:,krome_idx_32S334S5) + uall(:,krome_idx_32S334S436S) +  uall(:,krome_idx_32S334S336S2)&
  +uall(:,krome_idx_32S334S236S3) +  uall(:,krome_idx_32S334S36S4) + uall(:,krome_idx_32S336S5) +  uall(:,krome_idx_32S233S6)&
  +uall(:,krome_idx_32S233S534S) +  uall(:,krome_idx_32S233S536S) + uall(:,krome_idx_32S233S434S2) +  uall(:,krome_idx_32S233S434S36S)&
  +uall(:,krome_idx_32S233S436S2) +  uall(:,krome_idx_32S233S334S3) + uall(:,krome_idx_32S233S334S236S) +  uall(:,krome_idx_32S233S334S36S2)&
  +uall(:,krome_idx_32S233S336S3) +  uall(:,krome_idx_32S233S234S4) + uall(:,krome_idx_32S233S234S336S) +  uall(:,krome_idx_32S233S234S236S2)&
  +uall(:,krome_idx_32S233S234S36S3) +  uall(:,krome_idx_32S233S236S4) + uall(:,krome_idx_32S233S34S5) +  uall(:,krome_idx_32S233S34S436S)&
  +uall(:,krome_idx_32S233S34S336S2) +  uall(:,krome_idx_32S233S34S236S3) + uall(:,krome_idx_32S233S34S36S4) +  uall(:,krome_idx_32S233S36S5)&
  +uall(:,krome_idx_32S234S6) +  uall(:,krome_idx_32S234S536S) + uall(:,krome_idx_32S234S436S2) +  uall(:,krome_idx_32S234S336S3)&
  +uall(:,krome_idx_32S234S236S4) +  uall(:,krome_idx_32S234S36S5) + uall(:,krome_idx_32S236S6) +  uall(:,krome_idx_32S33S7)&
  +uall(:,krome_idx_32S33S634S) +  uall(:,krome_idx_32S33S636S) + uall(:,krome_idx_32S33S534S2) +  uall(:,krome_idx_32S33S534S36S)&
  +uall(:,krome_idx_32S33S536S2) +  uall(:,krome_idx_32S33S434S3) + uall(:,krome_idx_32S33S434S236S) +  uall(:,krome_idx_32S33S434S36S2)&
  +uall(:,krome_idx_32S33S436S3) +  uall(:,krome_idx_32S33S334S4) + uall(:,krome_idx_32S33S334S336S) +  uall(:,krome_idx_32S33S334S236S2)&
  +uall(:,krome_idx_32S33S334S36S3) +  uall(:,krome_idx_32S33S336S4) + uall(:,krome_idx_32S33S234S5) +  uall(:,krome_idx_32S33S234S436S)&
  +uall(:,krome_idx_32S33S234S336S2) +  uall(:,krome_idx_32S33S234S236S3) + uall(:,krome_idx_32S33S234S36S4) +  uall(:,krome_idx_32S33S236S5)&
  +uall(:,krome_idx_32S33S34S6) +  uall(:,krome_idx_32S33S34S536S) + uall(:,krome_idx_32S33S34S436S2) +  uall(:,krome_idx_32S33S34S336S3)&
  +uall(:,krome_idx_32S33S34S236S4) +  uall(:,krome_idx_32S33S34S36S5) + uall(:,krome_idx_32S33S36S6)&
  +uall(:,krome_idx_32S34S7) +  uall(:,krome_idx_32S34S636S) + uall(:,krome_idx_32S34S536S2) +  uall(:,krome_idx_32S34S436S3)&
  +uall(:,krome_idx_32S34S336S4) +  uall(:,krome_idx_32S34S236S5) + uall(:,krome_idx_32S34S36S6) +  uall(:,krome_idx_32S36S7)&
  +uall(:,krome_idx_33S8) +  uall(:,krome_idx_33S734S) + uall(:,krome_idx_33S736S) +  uall(:,krome_idx_33S634S2)&
  +uall(:,krome_idx_33S634S36S) +  uall(:,krome_idx_33S636S2) + uall(:,krome_idx_33S534S3) +  uall(:,krome_idx_33S534S236S)&
  +uall(:,krome_idx_33S534S36S2) +  uall(:,krome_idx_33S536S3) + uall(:,krome_idx_33S434S4) +  uall(:,krome_idx_33S434S336S)&
  +uall(:,krome_idx_33S434S236S2) +  uall(:,krome_idx_33S434S36S3) + uall(:,krome_idx_33S436S4) +  uall(:,krome_idx_33S334S5)&
  +uall(:,krome_idx_33S334S436S) +  uall(:,krome_idx_33S334S336S2) + uall(:,krome_idx_33S334S236S3) +  uall(:,krome_idx_33S334S36S4)&
  +uall(:,krome_idx_33S336S5) +  uall(:,krome_idx_33S234S6) + uall(:,krome_idx_33S234S536S) +  uall(:,krome_idx_33S234S436S2)&
  +uall(:,krome_idx_33S234S336S3) +  uall(:,krome_idx_33S234S236S4) + uall(:,krome_idx_33S234S36S5) +  uall(:,krome_idx_33S236S6)&
  +uall(:,krome_idx_33S34S7) +  uall(:,krome_idx_33S34S636S) + uall(:,krome_idx_33S34S536S2) +  uall(:,krome_idx_33S34S436S3)&
  +uall(:,krome_idx_33S34S336S4) +  uall(:,krome_idx_33S34S236S5) +uall(:,krome_idx_33S34S36S6) +  uall(:,krome_idx_33S36S7)&
  +uall(:,krome_idx_34S8) +  uall(:,krome_idx_34S736S) + uall(:,krome_idx_34S636S2) +  uall(:,krome_idx_34S536S3)&
  +uall(:,krome_idx_34S436S4) +  uall(:,krome_idx_34S336S5) + uall(:,krome_idx_34S236S6) +  uall(:,krome_idx_34S36S7)&
  +uall(:,krome_idx_36S8) +  uall(:,krome_idx_32S633S36S)

   
	  do i=2,n+1
	    open(54,file="XePa-all.txt",status='old',action='write',form='formatted',position="append")
      write(54,30) t,i-1,so2(i),so(i),s(i),o(i),singleta(i),triplet(i),ocs(i),co2(i),co(i),m(i),s2(i),so3(i) &
        ,S3(i),s4(i),s8(i),o2(i),s2o2(i),s2o(i)
      30 format(f10.2,I6,18e13.5)
    end do

    do i=2,n+1
	    open(65,file="sulfur-isotope-conc.txt",status='old',action='write',form='formatted',position="append")
      write(65,28) t,i-1,so2a(i),so2b(i),so2c(i),so2d(i),soa(i),sob(i),soc(i),sod(i) &
		 ,sa(i),sb(i),sc(i),sd(i),ocsa(i),ocsb(i),ocsc(i),ocsd(i) &
		 ,singleta32(i),singleta33(i),singleta34(i),singleta36(i),triplet32(i),triplet33(i),triplet34(i),triplet36(i) &
     ,so3a(i),so3b(i),so3c(i),so3d(i),s2a(i),s2b(i),s2c(i),s2d(i),s2e(i),s2f(i),s2g(i),s2h(i),s2i(i),s2j(i)&
     ,s3a(i),s3b(i),s3c(i),s3d(i),s3e(i),s3f(i),s3g(i),s3h(i),s3i(i),s3j(i),s3k(i),s3l(i),s3m(i),s3n(i),s3o(i),s3p(i),s3q(i),s3r(i),s3s(i),s3t(i) &
		 ,s4a(i),s4b(i),s4c(i),s4d(i),s4e(i),s4f(i),s4g(i),s4h(i),s4i(i),s4j(i),s4k(i),s4l(i),s4m(i),s4n(i),s4o(i),s4p(i),s4q(i),s4r(i),s4s(i),s4t(i) &
     ,s4u(i),s4v(i),s4w(i),s4x(i),s4y(i),s4z(i),s4aa(i),s4bb(i),s4cc(i),s4dd(i),s4ee(i),s4ff(i),s4gg(i),s4hh(i),s4ii(i) &
     ,s2o2a(i),s2o2b(i),s2o2c(i),s2o2d(i),s2o2e(i),s2o2f(i),s2o2g(i),s2o2h(i),s2o2i(i),s2o2j(i) &
     ,s2oa(i),s2ob(i),s2oc(i),s2od(i),s2oe(i),s2of(i),s2og(i),s2oh(i),s2oi(i),s2oj(i)
     
		  28 format(f10.2,I6,113e18.7)
      end do
	          
  end if
!****************************************************************************************
  if(t==0 .or.t==3600) then  
   so2_d33(:) = 1000d0 * log((so2b(:)/so2a(:))/coef_33)
   so2_d34(:) = 1000d0 * log((so2c(:)/so2a(:))/coef_34)
   so2_d36(:) = 1000d0 * log((so2d(:)/so2a(:))/coef_36)
   so2_Delta33(:) = so2_d33(:) - 0.515d0*so2_d34(:)
   so2_Delta36(:) = so2_d36(:) - 1.9d0*so2_d34(:)

   so_d33(:) = 1000d0 * log((sob(:)/soa(:))/coef_33)
   so_d34(:) = 1000d0 * log((soc(:)/soa(:))/coef_34)
   so_d36(:) = 1000d0 * log((sod(:)/soa(:))/coef_36)
   so_Delta33(:) = so_d33(:) - 0.515d0*so_d34(:)
   so_Delta36(:) = so_d36(:) - 1.9d0*so_d34(:)

   s_d33(:) = 1000d0 * log((sb(:)/sa(:))/coef_33)
   s_d34(:) = 1000d0 * log((sc(:)/sa(:))/coef_34)
   s_d36(:) = 1000d0 * log((sd(:)/sa(:))/coef_36)
   s_Delta33(:) = s_d33(:) - 0.515d0*s_d34(:)
   s_Delta36(:) = s_d36(:) - 1.9d0*s_d34(:)

   ocs_d33(:) = 1000d0 * log((ocsb(:)/ocsa(:))/coef_33)
   ocs_d34(:) = 1000d0 * log((ocsc(:)/ocsa(:))/coef_34)
   ocs_d36(:) = 1000d0 * log((ocsd(:)/ocsa(:))/coef_36)
   ocs_Delta33(:) = ocs_d33(:) - 0.515d0*ocs_d34(:)
   ocs_Delta36(:) = ocs_d36(:) - 1.9d0*ocs_d34(:)

   so3_d33(:) = 1000d0 * log((so3b(:)/so3a(:))/coef_33)
   so3_d34(:) = 1000d0 * log((so3c(:)/so3a(:))/coef_34)
   so3_d36(:) = 1000d0 * log((so3d(:)/so3a(:))/coef_36)
   so3_Delta33(:) = so3_d33(:) - 0.515d0*so3_d34(:)
   so3_Delta36(:) = so3_d36(:) - 1.9d0*so3_d34(:)

   singleta_d33(:) = 1000d0 * log((singleta33(:)/singleta32(:))/coef_33)
   singleta_d34(:) = 1000d0 * log((singleta34(:)/singleta32(:))/coef_34)
   singleta_d36(:) = 1000d0 * log((singleta36(:)/singleta32(:))/coef_36)
   singleta_Delta33(:) = singleta_d33(:) - 0.515d0*singleta_d34(:)
   singleta_Delta36(:) = singleta_d36(:) - 1.9d0*singleta_d34(:)

   singletb_d33(:) = 1000d0 * log((singletb33(:)/singletb32(:))/coef_33)
   singletb_d34(:) = 1000d0 * log((singletb34(:)/singletb32(:))/coef_34)
   singletb_d36(:) = 1000d0 * log((singletb36(:)/singletb32(:))/coef_36)
   singletb_Delta33(:) = singletb_d33(:) - 0.515d0*singletb_d34(:)
   singletb_Delta36(:) = singletb_d36(:) - 1.9d0*singletb_d34(:)

   triplet_d33(:) = 1000d0 * log((triplet33(:)/triplet32(:))/coef_33)
   triplet_d34(:) = 1000d0 * log((triplet34(:)/triplet32(:))/coef_34)
   triplet_d36(:) = 1000d0 * log((triplet36(:)/triplet32(:))/coef_36)
   triplet_Delta33(:) = triplet_d33(:) - 0.515d0*triplet_d34(:)
   triplet_Delta36(:) = triplet_d36(:) - 1.9d0*triplet_d34(:)

   s2_d33(:) = 1000d0 * log( ((s2b(:) + 2*s2e(:) + s2f(:) + s2g(:))/(2*s2a(:) + s2b(:) + s2c(:) + s2d(:)))/ coef_33)
   s2_d34(:) = 1000d0 * log( ((s2c(:) + s2f(:) + 2*s2h(:) + s2i(:))/(2*s2a(:) + s2b(:) + s2c(:) + s2d(:)))/ coef_34)
   s2_d36(:) = 1000d0 * log( ((s2d(:) + s2g(:) + s2i(:) + 2*s2j(:))/(2*s2a(:) + s2b(:) + s2c(:) + s2d(:)))/ coef_36)
   s2_Delta33(:) = s2_d33(:) - 0.515d0*s2_d34(:)
   s2_Delta36(:) = s2_d36(:) - 1.9d0*s2_d34(:)

   
   s2o_d33(:) = 1000d0 * log( ((s2ob(:) + 2*s2oe(:) + s2of(:) + s2og(:))/(2*s2oa(:) + s2ob(:) + s2oc(:) + s2od(:)))/ coef_33)
   s2o_d34(:) = 1000d0 * log( ((s2oc(:) + s2of(:) + 2*s2oh(:) + s2oi(:))/(2*s2oa(:) + s2ob(:) + s2oc(:) + s2od(:)))/ coef_34)
   s2o_d36(:) = 1000d0 * log( ((s2od(:) + s2og(:) + s2oi(:) + 2*s2oj(:))/(2*s2oa(:) + s2ob(:) + s2oc(:) + s2od(:)))/ coef_36)
   s2o_Delta33(:) = s2o_d33(:) - 0.515d0*s2o_d34(:)
   s2o_Delta36(:) = s2o_d36(:) - 1.9d0*s2o_d34(:)

   s2o2_d33(:) = 1000d0 * log( ((s2o2b(:) + 2*s2o2e(:) + s2o2f(:) + s2o2g(:))/(2*s2o2a(:) + s2o2b(:) + s2o2c(:) + s2o2d(:)))/ coef_33)
   s2o2_d34(:) = 1000d0 * log( ((s2o2c(:) + s2o2f(:) + 2*s2o2h(:) + s2o2i(:))/(2*s2o2a(:) + s2o2b(:) + s2o2c(:) + s2o2d(:)))/ coef_34)
   s2o2_d36(:) = 1000d0 * log( ((s2o2d(:) + s2o2g(:) + s2o2i(:) + 2*s2o2j(:))/(2*s2o2a(:) + s2o2b(:) + s2o2c(:) + s2o2d(:)))/ coef_36)
   s2o2_Delta33(:) = s2o2_d33(:) - 0.515d0*s2o2_d34(:)
   s2o2_Delta36(:) = s2o2_d36(:) - 1.9d0*s2o2_d34(:)

   s3_d33(:) = 1000d0 * log( (3*s3b(:) + s3e(:) + 2*(s3h(:)+ s3n(:)+ s3o(:)) + s3k(:) + s3l(:) + s3p(:) + s3q(:) + s3s(:))&
    /(3*s3a(:) + 2*(s3e(:) + s3f(:) + s3g(:)) + s3h(:) + s3i(:) + s3j(:) + s3k(:) + s3l(:) + s3m(:)) / coef_33)
   s3_d34(:) = 1000d0 * log( (3*s3c(:) + s3f(:) + 2*(s3i(:)+ s3p(:)+ s3r(:)) + s3k(:) + s3m(:) + s3n(:) + s3q(:) + s3t(:))&
    /(3*s3a(:) + 2*(s3e(:) + s3f(:) + s3g(:)) + s3h(:) + s3i(:) + s3j(:) + s3k(:) + s3l(:) + s3m(:)) / coef_34)
   s3_d36(:) = 1000d0 * log( (3*s3d(:) + s3g(:) + 2*(s3j(:)+ s3s(:)+ s3t(:)) + s3l(:) + s3m(:) + s3o(:) + s3q(:) + s3r(:))&
    /(3*s3a(:) + 2*(s3e(:) + s3f(:) + s3g(:)) + s3h(:) + s3i(:) + s3j(:) + s3k(:) + s3l(:) + s3m(:)) / coef_36)
   s3_Delta33(:) = s3_d33(:) - 0.515d0*s3_d34(:)
   s3_Delta36(:) = s3_d36(:) - 1.9d0*s3_d34(:)
 
   s4_d33(:) = 1000d0 * log( ( 4*s4b(:) + 3*(s4n(:) + s4x(:) + s4y(:))  + 2*(s4h(:) + s4o(:) + s4p(:) + s4z(:) + s4aa(:) + s4bb(:)) &
    + (s4e(:) + s4i(:) + s4j(:) + s4q(:) + s4r(:) + s4s(:)  + s4cc(:) + s4dd(:) + s4ee(:) + s4ff(:)) )  / ( 4*s4a(:)  + 3*(s4e(:) + s4f(:) + s4g(:)) &
    + 2*(s4h(:) + s4i(:) + s4j(:) + s4k(:) + s4l(:) + s4m(:)) + (s4n(:) + s4o(:) + s4p(:) + s4q(:) + s4r(:) + s4s(:) &
    + s4t(:) + s4u(:) + s4v(:) + s4w(:)) ) / coef_33 )
   s4_d34(:) = 1000d0 * log( ( 4*s4c(:) + 3*(s4t(:) + s4cc(:) + s4gg(:)) + 2*(s4k(:) + s4q(:) + s4u(:) + s4z(:) + s4dd(:) + s4hh(:)) &
    + ( s4f(:) + s4i(:) + s4l(:) + s4o(:) + s4r(:) + s4v(:) + s4x(:) + s4bb(:) + s4ee(:) + s4ii(:) ) ) / ( 4*s4a(:) + 3*(s4e(:) + s4f(:) + s4g(:)) &
    + 2*(s4h(:) + s4i(:) + s4j(:) + s4k(:) + s4l(:) + s4m(:))  + (s4n(:) + s4o(:) + s4p(:) + s4q(:) + s4r(:) + s4s(:) &
    + s4t(:) + s4u(:) + s4v(:) + s4w(:)) ) / coef_34 )
   s4_d36(:) = 1000d0 * log( ( 4*s4d(:) + 3*(s4w(:) + s4ff(:) + s4ii(:)) + 2*(s4m(:) + s4s(:) + s4v(:)+ s4aa(:) + s4ee(:) + s4hh(:)) &
    + (s4g(:) + s4j(:) + s4l(:) + s4p(:) + s4r(:) + s4u(:) +s4y(:) + s4bb(:) + s4dd(:) + s4gg(:)) ) / ( 4*s4a(:) + 3*(s4e(:) + s4f(:) + s4g(:)) &
    + 2*(s4h(:) + s4i(:) + s4j(:) + s4k(:) + s4l(:) + s4m(:)) + (s4n(:) + s4o(:) + s4p(:) + s4q(:) + s4r(:) + s4s(:) &
    + s4t(:) + s4u(:) + s4v(:) + s4w(:)) ) / coef_36 )
   s4_Delta33(:) = s4_d33(:) - 0.515d0*s4_d34(:)
   s4_Delta36(:) = s4_d36(:) - 1.9d0*s4_d34(:)

   Eight32S(:) = 8*(uall(:,krome_idx_32S8))
   Seven32S(:) = 7*(uall(:,krome_idx_32S733S)+ uall(:,krome_idx_32S734S)+  uall(:,krome_idx_32S736S))
   Six32S(:) = 6*(uall(:,krome_idx_32S633S2)+ uall(:,krome_idx_32S633S34S)+ uall(:,krome_idx_32S634S2)&
               + uall(:,krome_idx_32S634S36S)+ uall(:,krome_idx_32S636S2)+ uall(:,krome_idx_32S633S36S))
   Five32S(:) = 5*(uall(:,krome_idx_32S533S3)+ uall(:,krome_idx_32S533S234S)+ uall(:,krome_idx_32S533S236S)&
              + uall(:,krome_idx_32S533S34S2)+ uall(:,krome_idx_32S533S34S36S)+ uall(:,krome_idx_32S533S36S2)&
              + uall(:,krome_idx_32S534S3)+ uall(:,krome_idx_32S534S236S)+ uall(:,krome_idx_32S534S36S2)&
              + uall(:,krome_idx_32S536S3))
   Four32S(:) = 4*(uall(:,krome_idx_32S433S4)+ uall(:,krome_idx_32S433S334S)+ uall(:,krome_idx_32S433S336S)&
              + uall(:,krome_idx_32S433S234S2)+ uall(:,krome_idx_32S433S234S36S)+ uall(:,krome_idx_32S433S236S2)&
              + uall(:,krome_idx_32S433S34S3)+ uall(:,krome_idx_32S433S34S236S)+ uall(:,krome_idx_32S433S34S36S2)&
              + uall(:,krome_idx_32S433S36S3)+ uall(:,krome_idx_32S434S4)+ uall(:,krome_idx_32S434S336S)&
              + uall(:,krome_idx_32S434S236S2)+ uall(:,krome_idx_32S434S36S3)+ uall(:,krome_idx_32S436S4))
   Three32S(:) = 3*(uall(:,krome_idx_32S333S5)+ uall(:,krome_idx_32S333S434S)+ uall(:,krome_idx_32S333S436S)&
               + uall(:,krome_idx_32S333S334S2)+ uall(:,krome_idx_32S333S334S36S)+ uall(:,krome_idx_32S333S336S2)&
               + uall(:,krome_idx_32S333S234S3)+ uall(:,krome_idx_32S333S234S236S)+ uall(:,krome_idx_32S333S234S36S2)&
               + uall(:,krome_idx_32S333S236S3)+ uall(:,krome_idx_32S333S34S4)+ uall(:,krome_idx_32S333S34S336S)&
               + uall(:,krome_idx_32S333S34S236S2)+ uall(:,krome_idx_32S333S34S36S3)+ uall(:,krome_idx_32S333S36S4)&
               + uall(:,krome_idx_32S334S5)+ uall(:,krome_idx_32S334S436S)+ uall(:,krome_idx_32S334S336S2)&
               + uall(:,krome_idx_32S334S236S3)+ uall(:,krome_idx_32S334S36S4)+ uall(:,krome_idx_32S336S5))
   Two32S(:) =  2*(uall(:,krome_idx_32S233S6)+ uall(:,krome_idx_32S233S534S)+ uall(:,krome_idx_32S233S536S)&
             + uall(:,krome_idx_32S233S434S2)+ uall(:,krome_idx_32S233S434S36S)+ uall(:,krome_idx_32S233S436S2)&
             + uall(:,krome_idx_32S233S334S3)+ uall(:,krome_idx_32S233S334S236S)+ uall(:,krome_idx_32S233S334S36S2)&
             + uall(:,krome_idx_32S233S336S3)+ uall(:,krome_idx_32S233S234S4)+ uall(:,krome_idx_32S233S234S336S)&
             + uall(:,krome_idx_32S233S234S236S2)+ uall(:,krome_idx_32S233S234S36S3)+ uall(:,krome_idx_32S233S236S4)&
             + uall(:,krome_idx_32S233S34S5)+ uall(:,krome_idx_32S233S34S436S)+ uall(:,krome_idx_32S233S34S336S2)&
             + uall(:,krome_idx_32S233S34S236S3)+ uall(:,krome_idx_32S233S34S36S4)+ uall(:,krome_idx_32S233S36S5)&
             + uall(:,krome_idx_32S234S6)+ uall(:,krome_idx_32S234S536S)+ uall(:,krome_idx_32S234S436S2)&
             + uall(:,krome_idx_32S234S336S3)+ uall(:,krome_idx_32S234S236S4)+ uall(:,krome_idx_32S234S36S5)&
             + uall(:,krome_idx_32S236S6))
   One32S(:) =  uall(:,krome_idx_32S33S7)+ uall(:,krome_idx_32S33S634S)+ uall(:,krome_idx_32S33S636S)&
             + uall(:,krome_idx_32S33S534S2)+ uall(:,krome_idx_32S33S534S36S)+ uall(:,krome_idx_32S33S536S2)&
             + uall(:,krome_idx_32S33S434S3)+ uall(:,krome_idx_32S33S434S236S)+ uall(:,krome_idx_32S33S434S36S2)&
             + uall(:,krome_idx_32S33S436S3)+ uall(:,krome_idx_32S33S334S4)+ uall(:,krome_idx_32S33S334S336S)&
             + uall(:,krome_idx_32S33S334S236S2)+ uall(:,krome_idx_32S33S334S36S3)+ uall(:,krome_idx_32S33S336S4)&
             + uall(:,krome_idx_32S33S234S5)+ uall(:,krome_idx_32S33S234S436S)+ uall(:,krome_idx_32S33S234S336S2)&
             + uall(:,krome_idx_32S33S234S236S3)+ uall(:,krome_idx_32S33S234S36S4)+ uall(:,krome_idx_32S33S236S5)&
             + uall(:,krome_idx_32S33S34S6)+ uall(:,krome_idx_32S33S34S536S)+ uall(:,krome_idx_32S33S34S436S2)&
             + uall(:,krome_idx_32S33S34S336S3)+ uall(:,krome_idx_32S33S34S236S4)+ uall(:,krome_idx_32S33S34S36S5)&
             + uall(:,krome_idx_32S33S36S6)+ uall(:,krome_idx_32S34S7)+ uall(:,krome_idx_32S34S636S)&
             + uall(:,krome_idx_32S34S536S2)+ uall(:,krome_idx_32S34S436S3)+ uall(:,krome_idx_32S34S336S4)&
             + uall(:,krome_idx_32S34S236S5)+ uall(:,krome_idx_32S34S36S6)+ uall(:,krome_idx_32S36S7)

   Eight33S(:) = 8*(uall(:,krome_idx_33S8))
   Seven33S(:) = 7*(uall(:,krome_idx_32S33S7)+ uall(:,krome_idx_33S734S)+ uall(:,krome_idx_33S736S))
   Six33S(:) = 6*(uall(:,krome_idx_32S233S6)+ uall(:,krome_idx_32S33S634S)+ uall(:,krome_idx_32S33S636S)&
             + uall(:,krome_idx_33S634S2)+ uall(:,krome_idx_33S634S36S)+ uall(:,krome_idx_33S636S2))
   Five33S(:) = 5*(uall(:,krome_idx_32S333S5)+ uall(:,krome_idx_32S233S534S)+ uall(:,krome_idx_32S233S536S)&
              + uall(:,krome_idx_32S33S534S2)+ uall(:,krome_idx_32S33S534S36S)+ uall(:,krome_idx_32S33S536S2)&
              + uall(:,krome_idx_33S534S3)+ uall(:,krome_idx_33S534S236S)+ uall(:,krome_idx_33S534S36S2)&
              + uall(:,krome_idx_33S536S3))
   Four33S(:) = 4*(uall(:,krome_idx_32S433S4)+ uall(:,krome_idx_32S333S434S)+ uall(:,krome_idx_32S333S436S)&
              + uall(:,krome_idx_32S233S434S2)+ uall(:,krome_idx_32S233S434S36S)+ uall(:,krome_idx_32S233S436S2)&
              + uall(:,krome_idx_32S33S434S3)+ uall(:,krome_idx_32S33S434S236S)+ uall(:,krome_idx_32S33S434S36S2)&
              + uall(:,krome_idx_32S33S436S3)+ uall(:,krome_idx_33S434S4)+ uall(:,krome_idx_33S434S336S)&
              + uall(:,krome_idx_33S434S236S2)+  uall(:,krome_idx_33S434S36S3)+ uall(:,krome_idx_33S436S4))
   Three33S(:) = 3*(uall(:,krome_idx_32S533S3)+ uall(:,krome_idx_32S433S334S)+ uall(:,krome_idx_32S433S336S)&
               + uall(:,krome_idx_32S333S334S2)+ uall(:,krome_idx_32S333S334S36S)+ uall(:,krome_idx_32S333S336S2)&
               + uall(:,krome_idx_32S233S334S3)+ uall(:,krome_idx_32S233S334S236S)+ uall(:,krome_idx_32S233S334S36S2)&
               + uall(:,krome_idx_32S33S334S4)+ uall(:,krome_idx_32S33S334S336S)+ uall(:,krome_idx_32S33S334S236S2)&
               + uall(:,krome_idx_32S33S334S36S3)+ uall(:,krome_idx_32S233S336S3)+ uall(:,krome_idx_32S33S336S4)&
               + uall(:,krome_idx_33S334S5)+ uall(:,krome_idx_33S334S436S)+ uall(:,krome_idx_33S334S336S2)&
               + uall(:,krome_idx_33S334S236S3)+ uall(:,krome_idx_33S334S36S4)+ uall(:,krome_idx_33S336S5))
   Two33S(:) =  2*(uall(:,krome_idx_32S633S2)+ uall(:,krome_idx_32S533S234S)+ uall(:,krome_idx_32S533S236S)&
             + uall(:,krome_idx_32S433S234S2)+ uall(:,krome_idx_32S433S234S36S)+ uall(:,krome_idx_32S433S236S2)&
             + uall(:,krome_idx_32S333S234S3)+ uall(:,krome_idx_32S333S234S236S)+ uall(:,krome_idx_32S333S234S36S2)&
             + uall(:,krome_idx_32S333S236S3)+ uall(:,krome_idx_32S233S234S4)+ uall(:,krome_idx_32S233S234S336S)&
             + uall(:,krome_idx_32S233S234S236S2)+ uall(:,krome_idx_32S233S234S36S3)+ uall(:,krome_idx_32S233S236S4)&
             + uall(:,krome_idx_32S33S234S5)+ uall(:,krome_idx_32S33S234S436S)+ uall(:,krome_idx_32S33S234S336S2)&
             + uall(:,krome_idx_32S33S234S236S3)+ uall(:,krome_idx_32S33S234S36S4)+ uall(:,krome_idx_32S33S236S5)&
             + uall(:,krome_idx_33S234S6)+ uall(:,krome_idx_33S234S536S)+ uall(:,krome_idx_33S234S436S2)&
             + uall(:,krome_idx_33S234S336S3)+ uall(:,krome_idx_33S234S236S4)+ uall(:,krome_idx_33S234S36S5)&
             + uall(:,krome_idx_33S236S6))
   One33S(:) = uall(:,krome_idx_32S733S)+ uall(:,krome_idx_32S633S34S)+ uall(:,krome_idx_32S533S34S2)&
             + uall(:,krome_idx_32S533S34S36S)+ uall(:,krome_idx_32S533S36S2)+ uall(:,krome_idx_32S433S34S3)&
             + uall(:,krome_idx_32S433S34S236S)+ uall(:,krome_idx_32S433S34S36S2)+ uall(:,krome_idx_32S433S36S3)&
             + uall(:,krome_idx_32S333S34S4)+ uall(:,krome_idx_32S333S34S336S)+ uall(:,krome_idx_32S333S34S236S2)&
             + uall(:,krome_idx_32S333S34S36S3)+ uall(:,krome_idx_32S333S36S4)+ uall(:,krome_idx_32S233S34S5)&
             + uall(:,krome_idx_32S233S34S436S)+ uall(:,krome_idx_32S233S34S336S2)+ uall(:,krome_idx_32S233S34S236S3)&
             + uall(:,krome_idx_32S233S34S36S4)+ uall(:,krome_idx_32S233S36S5)+ uall(:,krome_idx_32S33S34S6)&
             + uall(:,krome_idx_32S33S34S536S)+ uall(:,krome_idx_32S33S34S436S2)+ uall(:,krome_idx_32S33S34S336S3)&
             + uall(:,krome_idx_32S33S34S236S4)+ uall(:,krome_idx_32S33S34S36S5)+ uall(:,krome_idx_32S33S36S6)&
             + uall(:,krome_idx_33S34S7)+ uall(:,krome_idx_33S34S636S)+ uall(:,krome_idx_33S34S536S2)&
             + uall(:,krome_idx_33S34S436S3)+ uall(:,krome_idx_33S34S336S4)+ uall(:,krome_idx_33S34S236S5)&
             + uall(:,krome_idx_33S34S36S6)+ uall(:,krome_idx_33S36S7)+ uall(:,krome_idx_32S633S36S)

   Eight34S(:) = 8*(uall(:,krome_idx_34S8))
   Seven34S(:) = 7*(uall(:,krome_idx_32S34S7)+ uall(:,krome_idx_34S736S)+ uall(:,krome_idx_33S34S7))
   Six34S(:) = 6*(uall(:,krome_idx_32S234S6)+ uall(:,krome_idx_32S33S34S6)+ uall(:,krome_idx_33S34S636S)&
             + uall(:,krome_idx_32S34S636S)+ uall(:,krome_idx_34S636S2)+ uall(:,krome_idx_33S234S6))
   Five34S(:) =  5*(uall(:,krome_idx_32S334S5)+ uall(:,krome_idx_32S233S34S5)+ uall(:,krome_idx_32S234S536S)&
              + uall(:,krome_idx_32S33S234S5)+ uall(:,krome_idx_32S33S34S536S)+ uall(:,krome_idx_32S34S536S2)&
              + uall(:,krome_idx_33S234S536S)+ uall(:,krome_idx_34S536S3)+ uall(:,krome_idx_33S334S5)&
              + uall(:,krome_idx_33S34S536S2))
   Four34S(:) = 4*(uall(:,krome_idx_32S434S4)+ uall(:,krome_idx_32S333S34S4)+ uall(:,krome_idx_32S334S436S)&
              + uall(:,krome_idx_32S233S234S4)+ uall(:,krome_idx_32S233S34S436S)+ uall(:,krome_idx_32S33S334S4)&
              + uall(:,krome_idx_32S234S436S2)+ uall(:,krome_idx_32S33S234S436S)+ uall(:,krome_idx_32S33S34S436S2)&
              + uall(:,krome_idx_32S34S436S3)+ uall(:,krome_idx_33S334S436S)+ uall(:,krome_idx_33S234S436S2)&
              + uall(:,krome_idx_34S436S4)+ uall(:,krome_idx_33S434S4)+ uall(:,krome_idx_33S34S436S3))
   Three34S(:) = 3*(uall(:,krome_idx_32S534S3)+ uall(:,krome_idx_32S433S34S3)+ uall(:,krome_idx_32S434S336S)&
               + uall(:,krome_idx_32S333S234S3)+ uall(:,krome_idx_32S333S34S336S)+ uall(:,krome_idx_32S334S336S2)&
               + uall(:,krome_idx_32S233S334S3)+ uall(:,krome_idx_32S233S234S336S)+ uall(:,krome_idx_32S233S34S336S2)&
               + uall(:,krome_idx_32S33S434S3)+ uall(:,krome_idx_32S234S336S3)+ uall(:,krome_idx_32S33S334S336S)&
               + uall(:,krome_idx_32S33S234S336S2)+ uall(:,krome_idx_32S33S34S336S3)+ uall(:,krome_idx_32S34S336S4)&
               + uall(:,krome_idx_33S534S3)+ uall(:,krome_idx_33S434S336S)+ uall(:,krome_idx_33S334S336S2)&
               + uall(:,krome_idx_33S234S336S3)+ uall(:,krome_idx_34S336S5)+ uall(:,krome_idx_33S34S336S4))
   Two34S(:) =  2*(uall(:,krome_idx_32S634S2)+ uall(:,krome_idx_32S533S34S2)+ uall(:,krome_idx_32S534S236S)&
                + uall(:,krome_idx_32S433S234S2)+ uall(:,krome_idx_32S433S34S236S)+ uall(:,krome_idx_32S333S334S2)&
                + uall(:,krome_idx_32S434S236S2)+ uall(:,krome_idx_32S333S234S236S)+ uall(:,krome_idx_32S333S34S236S2)&
                + uall(:,krome_idx_32S334S236S3)+ uall(:,krome_idx_32S233S434S2)+ uall(:,krome_idx_32S233S334S236S)&
                + uall(:,krome_idx_32S33S534S2)+ uall(:,krome_idx_32S233S234S236S2)+ uall(:,krome_idx_32S233S34S236S3)&
                + uall(:,krome_idx_32S33S434S236S)+ uall(:,krome_idx_32S234S236S4)+ uall(:,krome_idx_32S33S334S236S2)&
                + uall(:,krome_idx_32S33S234S236S3)+ uall(:,krome_idx_32S33S34S236S4)+ uall(:,krome_idx_32S34S236S5)&
                + uall(:,krome_idx_33S634S2)+ uall(:,krome_idx_33S534S236S)+ uall(:,krome_idx_33S434S236S2)&
                + uall(:,krome_idx_33S334S236S3)+ uall(:,krome_idx_33S234S236S4)+ uall(:,krome_idx_34S236S6)&
                + uall(:,krome_idx_33S34S236S5))
   One34S(:) = uall(:,krome_idx_32S734S)+ uall(:,krome_idx_32S633S34S)+ uall(:,krome_idx_32S634S36S)&
             + uall(:,krome_idx_32S533S234S)+ uall(:,krome_idx_32S533S34S36S)+ uall(:,krome_idx_32S534S36S2)&
             + uall(:,krome_idx_32S433S334S)+ uall(:,krome_idx_32S433S234S36S)+ uall(:,krome_idx_32S333S434S)&
             + uall(:,krome_idx_32S433S34S36S2)+ uall(:,krome_idx_32S434S36S3)+ uall(:,krome_idx_32S333S334S36S)&
             + uall(:,krome_idx_32S333S234S36S2)+ uall(:,krome_idx_32S333S34S36S3)+ uall(:,krome_idx_32S334S36S4)&
             + uall(:,krome_idx_32S233S534S)+ uall(:,krome_idx_32S233S434S36S)+ uall(:,krome_idx_32S33S634S)&
             + uall(:,krome_idx_32S233S334S36S2)+ uall(:,krome_idx_32S233S234S36S3)+ uall(:,krome_idx_32S33S534S36S)&
             + uall(:,krome_idx_32S233S34S36S4)+ uall(:,krome_idx_32S234S36S5)+ uall(:,krome_idx_32S33S434S36S2)&
             + uall(:,krome_idx_32S33S334S36S3)+ uall(:,krome_idx_32S33S234S36S4)+ uall(:,krome_idx_32S33S34S36S5)&
             + uall(:,krome_idx_32S34S36S6)+ uall(:,krome_idx_33S634S36S)+ uall(:,krome_idx_33S534S36S2)&
             + uall(:,krome_idx_33S434S36S3)+ uall(:,krome_idx_33S334S36S4)+ uall(:,krome_idx_33S234S36S5)&
             + uall(:,krome_idx_33S734S)+ uall(:,krome_idx_34S36S7)+ uall(:,krome_idx_33S34S36S6)

   Eight36S(:) = 8*(uall(:,krome_idx_36S8))
   Seven36S(:) = 7*(uall(:,krome_idx_32S36S7)+ uall(:,krome_idx_33S36S7)+ uall(:,krome_idx_34S36S7))
   Six36S(:) = 6*(uall(:,krome_idx_32S236S6)+ uall(:,krome_idx_32S33S36S6)+ uall(:,krome_idx_32S34S36S6)&
             + uall(:,krome_idx_33S236S6)+ uall(:,krome_idx_33S34S36S6)+ uall(:,krome_idx_34S236S6))
   Five36S(:) = 5*(uall(:,krome_idx_32S336S5)+ uall(:,krome_idx_32S233S36S5)+ uall(:,krome_idx_32S234S36S5)&
              + uall(:,krome_idx_32S33S236S5)+ uall(:,krome_idx_32S33S34S36S5)+ uall(:,krome_idx_32S34S236S5)&
              + uall(:,krome_idx_33S336S5)+ uall(:,krome_idx_33S234S36S5)+ uall(:,krome_idx_33S34S236S5)&
              + uall(:,krome_idx_34S336S5))
   Four36S(:) = 4*(uall(:,krome_idx_32S436S4)+ uall(:,krome_idx_32S333S36S4)+ uall(:,krome_idx_32S334S36S4)&
              + uall(:,krome_idx_32S233S236S4)+ uall(:,krome_idx_32S233S34S36S4)+ uall(:,krome_idx_32S234S236S4)&
              + uall(:,krome_idx_32S33S336S4)+ uall(:,krome_idx_32S33S234S36S4)+ uall(:,krome_idx_32S33S34S236S4)&
              + uall(:,krome_idx_32S34S336S4)+ uall(:,krome_idx_33S436S4)+ uall(:,krome_idx_33S334S36S4)&
              + uall(:,krome_idx_33S234S236S4)+ uall(:,krome_idx_33S34S336S4)+ uall(:,krome_idx_34S436S4))
   Three36S(:) = 3*(uall(:,krome_idx_32S536S3)+ uall(:,krome_idx_32S433S36S3)+ uall(:,krome_idx_32S434S36S3)&
               + uall(:,krome_idx_32S333S236S3)+ uall(:,krome_idx_32S333S34S36S3)+ uall(:,krome_idx_32S334S236S3)&
               + uall(:,krome_idx_32S233S336S3)+ uall(:,krome_idx_32S233S234S36S3)+ uall(:,krome_idx_32S233S34S236S3)&
               + uall(:,krome_idx_32S234S336S3)+ uall(:,krome_idx_32S33S436S3)+ uall(:,krome_idx_32S33S334S36S3)&
               + uall(:,krome_idx_32S33S234S236S3)+ uall(:,krome_idx_32S33S34S336S3)+ uall(:,krome_idx_32S34S436S3)&
               + uall(:,krome_idx_33S536S3)+ uall(:,krome_idx_33S434S36S3)+ uall(:,krome_idx_33S334S236S3)&
               + uall(:,krome_idx_33S234S336S3)+ uall(:,krome_idx_33S34S436S3)+ uall(:,krome_idx_34S536S3))
   Two36S(:) = 2*(uall(:,krome_idx_32S636S2)+ uall(:,krome_idx_32S533S36S2)+ uall(:,krome_idx_32S534S36S2)&
             + uall(:,krome_idx_32S433S236S2)+ uall(:,krome_idx_32S433S34S36S2)+ uall(:,krome_idx_32S434S236S2)&
             + uall(:,krome_idx_32S333S336S2)+ uall(:,krome_idx_32S333S234S36S2)+ uall(:,krome_idx_32S333S34S236S2)&
             + uall(:,krome_idx_32S334S336S2)+ uall(:,krome_idx_32S233S436S2)+ uall(:,krome_idx_32S233S334S36S2)&
             + uall(:,krome_idx_32S233S234S236S2)+ uall(:,krome_idx_32S233S34S336S2)+ uall(:,krome_idx_32S234S436S2)&
             + uall(:,krome_idx_32S33S536S2)+ uall(:,krome_idx_32S33S434S36S2)+ uall(:,krome_idx_32S33S334S236S2)&
             + uall(:,krome_idx_32S33S234S336S2)+ uall(:,krome_idx_32S33S34S436S2)+ uall(:,krome_idx_32S34S536S2)&
             + uall(:,krome_idx_33S636S2)+ uall(:,krome_idx_33S534S36S2)+ uall(:,krome_idx_33S434S236S2)&
             + uall(:,krome_idx_33S334S336S2)+ uall(:,krome_idx_33S234S436S2)+ uall(:,krome_idx_33S34S536S2)&
             + uall(:,krome_idx_34S636S2))
   One36S(:) = uall(:,krome_idx_32S736S)+ uall(:,krome_idx_32S633S36S) + uall(:,krome_idx_32S634S36S)&
             + uall(:,krome_idx_32S533S236S)+ uall(:,krome_idx_32S533S34S36S)+ uall(:,krome_idx_32S534S236S)&
             + uall(:,krome_idx_32S433S336S)+ uall(:,krome_idx_32S433S234S36S)+ uall(:,krome_idx_32S433S34S236S)&
             + uall(:,krome_idx_32S434S336S)+ uall(:,krome_idx_32S333S436S)+ uall(:,krome_idx_32S333S334S36S)&
             + uall(:,krome_idx_32S333S234S236S)+ uall(:,krome_idx_32S333S34S336S)+ uall(:,krome_idx_32S334S436S)&
             + uall(:,krome_idx_32S233S536S)+ uall(:,krome_idx_32S233S434S36S)+ uall(:,krome_idx_32S233S334S236S)&
             + uall(:,krome_idx_32S233S234S336S)+ uall(:,krome_idx_32S233S34S436S)+ uall(:,krome_idx_32S234S536S)&
             + uall(:,krome_idx_32S33S636S)+ uall(:,krome_idx_32S33S534S36S)+ uall(:,krome_idx_32S33S434S236S)&
             + uall(:,krome_idx_32S33S334S336S)+ uall(:,krome_idx_32S33S234S436S)+ uall(:,krome_idx_32S33S34S536S)&
             + uall(:,krome_idx_32S34S636S)+ uall(:,krome_idx_33S736S)+ uall(:,krome_idx_33S634S36S)&
             + uall(:,krome_idx_33S534S236S)+ uall(:,krome_idx_33S434S336S)+ uall(:,krome_idx_33S334S436S)&
             + uall(:,krome_idx_33S234S536S)+ uall(:,krome_idx_33S34S636S)+ uall(:,krome_idx_34S736S)
   
   s8_d33(:) = 1000d0 * log( (One33S(:) + Two33S(:) + Three33S(:) + Four33S(:) + Five33S(:) + Six33S(:) + Seven33S(:)+ Eight33S(:) )&
    / ( One32S(:) + Two32S(:) + Three32S(:) + Four32S(:) + Five32S(:) + Six32S(:) + Seven32S(:) + Eight32S(:) ) / coef_33 )
   
   s8_d34(:) = 1000d0 * log( (One34S(:) + Two34S(:) + Three34S(:) + Four34S(:) + Five34S(:) + Six34S(:) + Seven34S(:) + Eight34S(:) )&
    / ( One32S(:) + Two32S(:) + Three32S(:) + Four32S(:) + Five32S(:) + Six32S(:) + Seven32S(:) + Eight32S(:) ) / coef_34 )

   s8_d36(:) = 1000d0 * log( (One36S(:) + Two36S(:) + Three36S(:) + Four36S(:) + Five36S(:) + Six36S(:) + Seven36S(:) + Eight36S(:) )& 
    / ( One32S(:) + Two32S(:) + Three32S(:) + Four32S(:) + Five32S(:) + Six32S(:) + Seven32S(:) + Eight32S(:) ) / coef_36 )
   
   s8_Delta33(:) = s8_d33(:) - 0.515d0*s8_d34(:)
   s8_Delta36(:) = s8_d36(:) - 1.9d0*s8_d34(:)
 !****************************************************************************************  
	  
   do i=2,n+1
	  open(55,file="sulfur-isotope.txt",status='old',action='write',form='formatted',position="append")
    write(55,27) t,i-1,so2_d33(i),so2_d34(i),so2_d36(i),so2_Delta33(i),so2_Delta36(i) &
     ,so_d33(i),so_d34(i),so_d36(i),so_Delta33(i),so_Delta36(i) &
		 ,s_d33(i),s_d34(i),s_d36(i),s_Delta33(i),s_Delta36(i) &
     ,ocs_d33(i),ocs_d34(i),ocs_d36(i),ocs_Delta33(i),ocs_Delta36(i) &
     ,so3_d33(i),so3_d34(i),so3_d36(i),so3_Delta33(i),so3_Delta36(i) &
		 ,singleta_d33(i),singleta_d34(i),singleta_d36(i),singleta_Delta33(i),singleta_Delta36(i) &
     ,singletb_d33(i),singletb_d34(i),singletb_d36(i),singletb_Delta33(i),singletb_Delta36(i) &
     ,triplet_d33(i),triplet_d34(i),triplet_d36(i),triplet_Delta33(i),triplet_Delta36(i) &
     ,s2_d33(i),s2_d34(i),s2_d36(i),s2_Delta33(i),s2_Delta36(i) &
     ,s2o_d33(i),s2o_d34(i),s2o_d36(i),s2o_Delta33(i),s2o_Delta36(i) &
     ,s2o2_d33(i),s2o2_d34(i),s2o2_d36(i),s2o2_Delta33(i),s2o2_Delta36(i) &
     ,s3_d33(i),s3_d34(i),s3_d36(i),s3_Delta33(i),s3_Delta36(i) &
		 ,s4_d33(i),s4_d34(i),s4_d36(i),s4_Delta33(i),s4_Delta36(i) &
     ,s8_d33(i),s8_d34(i),s8_d36(i),s8_Delta33(i),s8_Delta36(i)    
     
		27 format(f10.2,I6,70e18.7)
   end do
  end if

  if(t==3600) then         
    open(73,file="tau-evolutionPa.txt",status='old',action='write',form='formatted',position="append")
    write (73,*) "time", t !"istep", istep, "time", t/36000d0, "days units"
    do j=1,krome_nPhotoBins
        write (73,*) j, (tauAll(i,j), i=1,n+1)
        13 format(I5,999E17.4e3)
    end do
  end if
    !dump to file
  if(mod(istep,10)==0) then
    print '(F11.2,a2)',t/tmax*1d2," %"
  end if
	 

  t = t + dt !increase timestep
     
	 !break when max time overshot
  if(t>tmax) exit
	 
	  istep = istep + 1 !increase timestep 18/06/25

end do

    !say goodbye
print *,"done, bye!!!"

end program test


