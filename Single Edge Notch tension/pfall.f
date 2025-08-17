!---------------------------------------------------------------------------------!
!----------------------UEL for AT-1, AT-2 and PFCZM-------------------------------!
!		Two phase-field model for UD composites using PUCK failure criteria in 2D !
!   The codes are well tested for AT2, but the code for AT1, PFCZm model is also implemented!
!---------------------------------------------------------------------------------!
!  Copyright (c) 2025, Technical University of Vienna, Austria.
! The code is distributed under a BSD license    

!Author: Pavan Kumar Asur Vijaya Kumar 
! Email: pavan.kumar@ilsb.tuwien.ac.at 
!Author: Aamir Dean
! Email: a.dean@isd.uni-hannover.de
! ** Feel free to write me suggestions, improvements, and additions. 
!the material properties 
! props(1)= Angel of the ply
! props(2)= Lc for the fiber failure
! props(3)= Gc for the fiber failure
! props(4)= Lc for the matrix failure
! props(5)= Gc for the matric failure
! props(6)= A small stabilization parameter (Keep it close to zero, or 1e-7)
! props(7)= KflagF !(AT2=0, AT1=1, PFCZM Bilinear=2, PFCZM Expo=3)
! props(8)= KflagD  Split (0: No split, 1: Amor, 2: Miehe) 
!Note that, PUCK doesnt need any split. Hence, only option 0 is currentluy available 
! props(9)= Fiber dimensionless driving parameter 
! props(10)= Matrix dimensionless driving parameter 
! props(11)= Threshold energy for the Fiber failure
! props(12)=  Threshold energy for the matrix failure
! Elastic properties of the composites in the local ply setting
! props(13)= E11 ! of fiber
! props(14)= E22
! props(15)= E12
! props(16)= Nu12
! props(17)= Nu23
! props(18)= G12
! Strength properties of the composites in the local ply setting
! props(19)= R1T  ! of fiber
! props(20)= R2T
! props(21)= R1C
! props(22)= R2C
! props(23)= R12
!Proprties for the PUCK failure theory in the local ply setting
! props(24)=P21_plus
! props(25)=P21_minus
! props(26)=P22_plus
! props(27)=P22_minus
! props(28)=E11F
! props(29)=V12F
! props(30)=magnitude
!props(31)=alpha

      module kvisual
      implicit none
      real*8 UserVar(70000,16,4)
      integer nelem
      !Here 70000 is the maximum number of elements in each layer 
      !* Change this if the maximum number of elements at each layer is greater than 70000
      ! 16 is the number of local state variables. 
      ! 4 is the number of the integration points. 
      save
      end module

      subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     1 props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     2 kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf,
     3 lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njpro,period)

      use kvisual
      include 'aba_param.inc' !implicit real(a-h o-z)

      dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),svars(*),
     1 energy(*),coords(mcrd,nnode),u(ndofel),du(mlvarx,*),v(ndofel),
     2 a(ndofel),time(2),params(*),jdltyp(mdload,*),adlmag(mdload,*),
     3 ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*)

      parameter(ndim=2,ntens=3,ninpt=4,nsvint=16)

      dimension wght(ninpt),dN(nnode,1),dNdz(ndim,nnode),stran(ntens),
     2 dNdx(ndim,nnode),b(ntens,nnode*ndim),ddsdde_local(ntens,ntens),
     3 stress_glo(ntens),dstran(ntens),statevLocal(nsvint), Sparam(ndim,ndim),
     4 Sinv(ndim,ndim),SinvT(ndim,ndim), flags(6),flags_old(6),stress_loc(ntens)
      dimension drot(3,3),stran_p(ntens),stress_DegGlobal(ntens),ddsdde_degLocal(ntens, ntens)
      dimension grot(3), gtot(3), ddsdde_glo(ntens, ntens),drotInv(3,3)
      dimension stress_deg1(ntens),stiff(8,8), force(8,1), u_con(8,1),stiff_deg(8,8)
      dimension  drotInvT(3,3),stress_degLocal(ntens),drot_E(3,3),stran_E(3),drotT(3,3)
      dimension  strain_global(ntens), strain_local(ntens),ddsdde_deglobal(3,3)
      dimension  stru_f(2,2), stru_m(2,2)
      data wght /1.d0, 1.d0, 1.d0, 1.d0/
		
!     initialising
      do k1=1,ndofel
       rhs(k1,1)=0.d0
      end do
      amatrx=0.d0

!     find number of elements
      if (dtime.eq.0.d0) then
       if (jelem.eq.1) then
        nelem=jelem
       else
        if (jelem.gt.nelem) nelem=jelem
       endif
      endif
      if (jelem.eq.1) then
      Ac=0.d0 ! Could be important to track the length of the crack (Not implemented here)
      endif

!     Reading parameters
	  angel=props(1)
	  xlcf=props(2)
          Gcf=props(3) 
	  xlcif= props(4) 
	  Gcif= props(5) 
      xk= props(6)
      kFlagF= props(7) !(AT2=0, AT1=1, PFCZM Bilinear=2, PFCZM Expo=3)
	  ft=props(7) !only for the PFCZM
	   kflagD=props(8)! Split (0: No split)
	Sf=props(9)
	Sm=props(10)
	! Compute the angle of the fiber orientation in Radians
	   theta=angel*3.1415926535897932384626433832d0/180.0d0
		!compute the rotation tensor for the C matrix
		drot(:,:)=0.0d0 !this is for the stress_global and the elastic matrix
	    ! Compute the rotation tensor based on the angel provided
		drot(1,1)=cos(theta)**2
		drot(1,2)=sin(theta)**2
		drot(2,1)=sin(theta)**2
		drot(2,2)=cos(theta)**2		
		drot(1,3)=2.0d0*cos(theta)*sin(theta)
		drot(2,3)=-2.0d0*cos(theta)*sin(theta) 
		drot(3,1)=-cos(theta)*sin(theta)
		drot(3,2)=cos(theta)*sin(theta)		
		drot(3,3)=cos(theta)**2-sin(theta)**2
		! Compute the Structural Tensors 
		!All these computations are done here to save the computational effort
		!compute the structural tensor for the matrix
		stru_m(1,1)=cos(theta)**2
		stru_m(1,2)=cos(theta)*sin(theta)
		stru_m(2,1)=cos(theta)*sin(theta)
		stru_m(2,2)=sin(theta)**2
		!compute the structural tensor for the fiber
		stru_f(1,1)=sin(theta)**2
		stru_f(1,2)=-cos(theta)*sin(theta)
		stru_f(2,1)=-cos(theta)*sin(theta)
		stru_f(2,2)=cos(theta)**2
		
		!transpose of the rotataion matrix
		drotT=transpose(drot)
		!inverse of the rotation matrix
		call kmatinv3(drot,drotInv)
		!Inverse transpose of the rotation matrix
		drotInvT=transpose(drotInv)
		

! The Loop over Integration points

      do kintk=1,ninpt
!     evaluate shape functions and derivatives
       call kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)
       call kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
       dvol=wght(kintk)*djac ! compute the volume of the element

!     form B-matrix
       b=0.d0
       do inod=1,nnode
        b(1,2*inod-1)=dNdx(1,inod)
        b(2,2*inod)=dNdx(2,inod)
        b(3,2*inod-1)=dNdx(2,inod)
        b(3,2*inod)=dNdx(1,inod)
       end do
	   !Import the state variables from previous step
      call kstatevar(kintk,nsvint,svars,statevLocal,1)

       stress_glo=statevLocal(1:ntens) ! number 1,2,3 are global stresses
       stran(1:ntens)=statevLocal((ntens+1):(2*ntens))!4, 5, 6 are the global strains
       phin_f=statevLocal(2*ntens+1)!7 is the phase field in fiber
       Hnif=statevLocal(2*ntens+2)!8 Driving force of Matrix
       Hnf=statevLocal(2*ntens+3)!9 Driving force of Fiber
       phin_if=statevLocal(2*ntens+4)!10 Phase field of matrix
	  flags_old(1:6)=statevLocal(11:16)!11,12,13,14,15, 16 Flags from the PUCK failure theory
!     compute from nodal values
       phi_f=0.d0 !initialize the phase field fiber
       phi_if=0.d0 !initialize the phase field Matrix
       do inod=1,nnode
        phi_f=phi_f+dN(inod,1)*u(8+inod)
        phi_if=phi_if+dN(inod,1)*u(12+inod)
       end do
	
!Enforce the Maximum value (Min value is not enforced here)	   
       if (phi_f.gt.1.d0) phi_f=1.d0
       if (phi_if.gt.1.d0) phi_if=1.d0

	if(phi_if.gt.10d0) call XIT
!Compute the geometric degradation functions corresponding to fiber and inter-fiber seperately
	call kgeometricdeg(kFlagF,phi_f, w_f, dw_f,ddw_f, cw)
	call kgeometricdeg(kFlagF,phi_if, w_if, dw_if,ddw_if, cw)

!Compute the degradation functions corresponding to fiber and inter-fiber seperately
      call KDegFun(kFlagF,phi_f,xk,a1,g_f,dg_f,ddg_f)
      call KDegFun(kFlagF,phi_if,xk,a1,g_if,dg_if,ddg_if)

!     compute the strain.
!The strain computed as B.U always give global strains
       strain_global=matmul(b,u(1:8)) ! this is already giving shear angles.
	   ! compute the strains in the local setting
	   !Local here means the angel of fiber.  Strain_l=T^{-T}*Strain_g
	   strain_local=matmul(drotInvT,strain_global )
      
       if (dtime.eq.0.d0) phin_f=phi_f
       if (dtime.eq.0.d0) phin_if=phi_if

	!start the material subroutine
       call kumat(props,ddsdde_local,ddsdde_glo,stress_glo,stress_loc,strain_global,ntens,statevLocal,drot)

		
		
       Psi=0.d0 ! Initialise the total energy
       do k1=1,ntens
        Psi=Psi+stress_loc(k1)*strain_local(k1)*0.5d0 ! Total Strain Energy
       end do
	   
	   
		
!compute the energy from fiber 
		Psi_f=0.0d0
		Psi_f=Psi_f+ stress_loc(1)*strain_local(1)*0.5d0
!Compute the energy from inter-fiber 
		Psi_if=0.0d0
		Psi_if=Psi_if+ (stress_loc(2)*strain_local(2)+stress_loc(3)*strain_local(3))*0.5d0
		
		!Initialise the variable for the comparison with PUCK 
		fiber_energy=0.0d0
		Psi_drive_f=0.0d0		
		fmatrix_energy=0.0d0
		Psi_drive_if=0.0d0
			
!compute tthe fiber driving energy by input from PUCK			
			
		if (statevLocal(11).ge.1) then 
		  Psi_drive_f= Sf*(Psi_f/flags_old(3) -1.0d0)	
	   fiber_energy=flags_old(3) ! the energy will constant.  Refers to the first instance when Flag(fiber)=1. 
	   !if the Puck fiber failure is greater than 1, then, driving force is computed
		elseif (statevLocal(11).lt.1) then 
		Psi_drive_f=0.0d0
		fiber_energy=Psi_f ! save the fiber energy for future
		! If the Puck fiber failure is less than 1, no driving force is computed
		endif 

!compute the matrix driving energy by input from PUCK	
		if (statevLocal(12).ge.1) then 
		!if the Puck matrix failure is greater than 1, then, driving force is computed
		  Psi_drive_if= Sm*(Psi_if/flags_old(5) -1.0d0)
	    fmatrix_energy=flags_old(5) 
		elseif (statevLocal(12).lt.1.0d0) then 
		Psi_drive_if=0.0d0
		fmatrix_energy=Psi_if
		! If the Puck matrix failure is less than 1, no driving force is computed
		endif 
		
		
!Enforcing KKt condition for Fiber 		
		Hf=max(Psi_drive_f,Hnf)
!Enforcing KKt condition for inter-Fiber 
		 Hif=max(Psi_drive_if, Hnif)
		
! Save the Flags and the driving forces		
		flags(1)=statevLocal(11)
		flags(2)=statevLocal(12)
		flags(3)=fiber_energy
		flags(4)=Psi_drive_f
		flags(5)=fmatrix_energy
		flags(6)=Psi_drive_if

!Assign and update the state variables for next increment 
       statevLocal(1:3)=stress_glo(1:3) !,1,2,3
       statevLocal(4:6)=strain_global(1:ntens) !4,5,6
       statevLocal(7)=phi_f !7
       statevLocal(8)=Hif !8
       statevLocal(9)=Hf !9
       statevLocal(10)=phi_if !10
	statevLocal(11:16)=flags(1:6) !10,11,12,13,14,15
	  	   
! Collect the Degradation functions from the fiber and matrix
	   gtot(1)=g_f
	   gtot(2)=g_if
	   gtot(3)=min(g_f, g_if)

!Compute the local Degrdation of the Elasticity tensor

		ddsdde_degLocal=0.0d0
		!Degrade the stiffness matrix
		ddsdde_degLocal(1,1)=gtot(1)*ddsdde_local(1,1)
		ddsdde_degLocal(1,2)=gtot(2)*ddsdde_local(1,2)
		ddsdde_degLocal(2,1)=gtot(2)*ddsdde_local(2,1)
		ddsdde_degLocal(2,2)=gtot(2)*ddsdde_local(2,2)
		ddsdde_degLocal(3,3)=gtot(3)*ddsdde_local(3,3)		
		ddsdde_deglobal=0.0d0
!Compute the Global Degrdation of the Elasticity tensor by pull back operator
		ddsdde_deglobal= matmul(drotInv, matmul(ddsdde_degLocal, drotInvT))

		
!degrdaded stresses in the local setting
	   		stress_degLocal(1)=stress_loc(1)*gtot(1)
			stress_degLocal(2)=stress_loc(2)*gtot(2)
			stress_degLocal(3)=stress_loc(3)*gtot(3)
	
!Compute the Degraded stressses in the Gloabl Setting using Pull back for stress
		stress_DegGlobal=matmul(drotInv, stress_degLocal)
		

       call kstatevar(kintk,nsvint,svars,statevLocal,0)
!Compute the AMATRX for the displacement using global degraded C matrix
       amatrx(1:8,1:8)=amatrx(1:8,1:8)+!stiff
     1 dvol*(matmul(matmul(transpose(b),ddsdde_deglobal),b)) 
!Compute the RHS for the displacement using global degraded C matrix
       rhs(1:8,1)=rhs(1:8,1)-!force(1:8,1)
     1  dvol*(matmul(transpose(b),stress_DegGlobal))
       
!Compute the AMATRX for the Fiber contribution 
        amatrx(9:12,9:12)=amatrx(9:12,9:12)
     1 +dvol*(matmul(transpose(dNdx),matmul(stru_f, dNdx))*2*Gcf*xlcf/cw 
     2 +matmul(dN,transpose(dN))*(Gcf/(cw*xlcf)*ddw_f+ddg_f*Hf)) !Fiber

!Compute the AMATRX for the Matrix contribution     
        amatrx(13:16,13:16)=amatrx(13:16,13:16)
     1 +dvol*(matmul(transpose(dNdx),matmul(stru_m, dNdx))*2*Gcif*xlcif/cw 
     2 +matmul(dN,transpose(dN))*(Gcif/(cw*xlcif)*ddw_if+ddg_if*Hif))! inter-fiber

!Compute the RHS for the Fiber contribution	 
        rhs(9:12,1)=rhs(9:12,1)
     1 -dvol*(matmul(transpose(dNdx),matmul(matmul(stru_f, dNdx),u(9:12)))*2*Gcf*xlcf/cw 
     2 + dN(:,1)*(dg_f*Hf+ Gcf*dw_f/(cw*xlcf))) !fiber

!Compute the RHS for the Matrix contribution 	 
	     rhs(13:16,1)=rhs(13:16,1)
     1 -dvol*(matmul(transpose(dNdx),matmul(matmul(stru_m, dNdx),u(13:16)))*2*Gcif*xlcif/cw  
     2 + dN(:,1)*(dg_if*Hif+ Gcif*dw_if/(cw*xlcif))) !inter fiber
	 
		
	 
! output for visualisation
       UserVar(jelem,1:3,kintk)=stress_degLocal(1:3) ! degraded Stress
       UserVar(jelem,4:6,kintk)=strain_local(1:3) ! Strains 
       UserVar(jelem,7:16,kintk)=statevLocal(7:16) ! 2 phase fields, 2 driving forces, 2 flags, 2 Puck energies, 2 current energy

		

      end do       ! end loop on material integration points

      RETURN
      END

      subroutine kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)
c
      include 'aba_param.inc'
c
      parameter (gaussCoord=0.577350269d0)
      dimension dN(nnode,1),dNdz(ndim,*),coord24(2,4)

      data  coord24 /-1.d0, -1.d0,
     2                1.d0, -1.d0,
     3               -1.d0,  1.d0,
     4                1.d0,  1.d0/

!     2D 4-nodes

!     determine (g,h,r)
      g=coord24(1,kintk)*gaussCoord
      h=coord24(2,kintk)*gaussCoord

!     shape functions
      dN(1,1)=(1.d0-g)*(1.d0-h)/4.d0
      dN(2,1)=(1.d0+g)*(1.d0-h)/4.d0
      dN(3,1)=(1.d0+g)*(1.d0+h)/4.d0
      dN(4,1)=(1.d0-g)*(1.d0+h)/4.d0

!     derivative d(Ni)/d(g)
      dNdz(1,1)=-(1.d0-h)/4.d0
      dNdz(1,2)=(1.d0-h)/4.d0
      dNdz(1,3)=(1.d0+h)/4.d0
      dNdz(1,4)=-(1.d0+h)/4.d0

!     derivative d(Ni)/d(h)
      dNdz(2,1)=-(1.d0-g)/4.d0
      dNdz(2,2)=-(1.d0+g)/4.d0
      dNdz(2,3)=(1.d0+g)/4.d0
      dNdz(2,4)=(1.d0-g)/4.d0

      return
      end

      subroutine kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
!     Notation: djac - Jac determinant; xjaci - inverse of Jac matrix
!     dNdx - shape functions derivatives w.r.t. global coordinates
      include 'aba_param.inc'

      dimension xjac(ndim,ndim),xjaci(ndim,ndim),coords(mcrd,nnode),
     1 dNdz(ndim,nnode),dNdx(ndim,nnode),dNdxb(ndim,nnode)

      xjac=0.d0

      do inod=1,nnode
       do idim=1,ndim
        do jdim=1,ndim
         xjac(jdim,idim)=xjac(jdim,idim)+
     1        dNdz(jdim,inod)*coords(idim,inod)
        end do
       end do
      end do

      djac=xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1)
      if (djac.gt.0.d0) then ! jacobian is positive - o.k.
       xjaci(1,1)=xjac(2,2)/djac
       xjaci(2,2)=xjac(1,1)/djac
       xjaci(1,2)=-xjac(1,2)/djac
       xjaci(2,1)=-xjac(2,1)/djac
      else ! negative or zero jacobian
       write(7,*)'WARNING: element',jelem,'has neg. Jacobian'
      endif

      dNdx= matmul(xjaci,dNdz)

	

      return
      end

c*****************************************************************
      subroutine kstatevar(npt,nsvint,statev,statev_ip,icopy)
c
c     Transfer data to/from element-level state variable array from/to
c     material-point level state variable array.
c
      include 'aba_param.inc'

      dimension statev(*),statev_ip(*)

      isvinc=(npt-1)*nsvint     ! integration point increment

      if (icopy.eq.1) then ! Prepare arrays for entry into umat
       do i=1,nsvint
        statev_ip(i)=statev(i+isvinc)
       enddo
      else ! Update element state variables upon return from umat
       do i=1,nsvint
        statev(i+isvinc)=statev_ip(i)
       enddo
      end if

      return
      end

c*****************************************************************
      subroutine kumat(props,ddsdde_local, ddsdde_glo, stress_glo, stress_loc,strain_global,ntens,statev,drot)
c
c     Subroutine with the material model
c
      include 'aba_param.inc' !implicit real(a-h o-z)

      dimension props(*),ddsdde_local(ntens,ntens),stress_glo(ntens),statev(*)
      dimension strain_global(ntens),drot(3,3),stress_loc(ntens),ddsdde_glo(ntens, ntens)
      dimension drotInv(3,3), drotInvT(3,3), drotT(3,3)	
	
!     Initialization
      ddsdde_local=0.d0
      
!transpose of the rotataion matrix
		drotT=transpose(drot)
!inverse of the rotation matrix
		call kmatinv3(drot,drotInv)
!Inverse transpose of the rotation matrix
		drotInvT=transpose(drotInv)
		E11=props(13)
		E22=props(14) 
		xnu12=props(16)  !minor
		g12=props(18) 
		xnu21=xnu12*E11/E22 !major
		
		eg=(1.0d0-xnu12*xnu21)

!compute the Local Elasticty C matrix
		
			ddsdde_local(1,1)=E11/eg
			ddsdde_local(2,2)=E22/eg
			ddsdde_local(1,2)=xnu12*E11/eg
			ddsdde_local(2,1)=xnu12*E11/eg		
			ddsdde_local(3,3)=g12
! Compute the Global Elasticty matrix by Push forard operator			
		ddsdde_glo= matmul(drotInv, matmul(ddsdde_local, drotInvT))
!Compute the stresses at the global level as Sigma_g=C_g*Strain_g
		stress_glo=matmul(ddsdde_glo,strain_global ) !Global
!Compute the stress at the local level.  Stress_l=T*Stress_g
		stress_loc=matmul(drot, stress_glo)
				
!Initiate the Puck Failure criteria. Inputs: Properties, Stress in local, State varibale for the Flags	
	   call puck_criteria(props, ntens, stress_loc, statev)
		
      return
      end

c*****************************************************************
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1 drplde,drpldt,stran,dstran,time,dtime,temp2,dtemp,predef,dpred,
     2 cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      use kvisual
      include 'aba_param.inc' !implicit real(a-h o-z)

      character*8 cmname
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1 ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),
     2 time(2),predef(1),dpred(1),props(nprops),coords(3),drot(3,3),
     3 dfgrd0(3,3),dfgrd1(3,3),jstep(4)

      ddsdde=0.0d0
      noffset=noel-nelem
	  !Write(7,*)' After noffset',noffset
      statev(1:16)=UserVar(noffset,1:16,npt)
	  !Write(7,*)' statev statev',statev
      return
      end
c*********************************************************************
      subroutine kgeometricdeg(kFlagF,phi, w, dw,ddw, cw)
 	include 'aba_param.inc' !implicit real(a-h o-z)


	if (kFlagF==0) then !AT2
	w=phi**2
	dw=2.0d0*phi
	ddw=2.0d0
	cw=2.0d0
	elseif (kFlagF==1) then !AT1
	w=phi
	dw=1.0d0
	ddw=0.0d0
	cw=8.0d0/3.0d0
	else  !PF-CZM linear softening
	w=2.d0*phi-phi**2
	dw=2.0d0-2.0d0*phi
	ddw=-2.0d0
	cw=3.1415926535897932384626433832d0
	endif


      end
c*********************************************************************
      subroutine KDegFun(kFlagF,phi,xk,a1,g,dg,ddg)

      include 'aba_param.inc'




      if (kFlagF.eq.0.or.kFlagF.eq.1) then ! AT2 & AT1 model
       g=(1.d0-phi)**2+xk
       dg=-2.d0*(1.d0-phi)
       ddg=2.d0
      else ! PF-CZM model
       if (kFlagF.eq.2) then ! Linear PF-CZM
	  p  =  2.0d0
        a2 =  -0.5d0
        a3 =  0.0d0
       elseif (kFlagF.eq.3) then ! Exponential PF-CZM
	p  =  2.5d0
        a2 =  2.0d0**(5.d0/3.d0) - 3.0d0
        a3 =  0.0d0
       endif
	fac1=(1.d0-phi)**p
        dfac1= -p*(1.d0 - phi)**(p-1.d0)
        ddfac1=  p*(p-1.d0)*(1.d0-phi)**(p-2.d0)

	fac2=fac1+a1*phi+a1*a2*phi**2.d0+a1*a2*a3*phi**3.d0
        dfac2=dfac1+a1+2.d0*a1*a2*phi+3.d0*a1*a2*a3*phi**2.d0
        ddfac2=ddfac1+2.d0*a1*a2+6.d0*a1*a2*a3*phi
	g=fac1/fac2
        dg=(dfac1*fac2-fac1*dfac2)/(fac2**2.d0)
        ddg=((ddfac1*fac2-fac1*ddfac2)*fac2
     1 -2.d0*(dfac1*fac2-fac1*dfac2)*dfac2)/(fac2**3.d0)

      end if

      end
c**********************************************************************************
      subroutine KDriFor(props,ntens,stran,stress,kFlagF,
     1 kflagD,psit,g,H,Hmin)

      include 'aba_param.inc'

      dimension stress(ntens),stran(ntens),props(*),
     1 Edev(ntens),PS(3),AN(3,3)
	
	    E=200000
	    xnu=0.3d0
	   eg=E/(1.d0+xnu)/2.d0
        elam=(E/(1.d0-2.d0*xnu)-2.d0*eg)/3.d0
        bk=E/(1.d0-2.d0*xnu)/3.d0
	   Gc=props(4)
	   xl=props(3)
	   ft=props(7)

      if (kFlagF.eq.0.or.kFlagF.eq.1) then
       psip=0.d0
       
       if (kflagD.eq.1) then ! Amor et al.
        CALL SPRINC(stran,PS,2,3,1)
        Edev=PS
        trE=PS(1)+PS(2)+PS(3)
        Edev(1:3)=PS(1:3)-trE/3.d0
        EdevS=Edev(1)**2+Edev(2)**2+Edev(3)**2
        trEp=0.5d0*(trE+abs(trE))
        trEn=0.5d0*(trE-abs(trE))
        psip=0.5d0*bk*trEp**2+eg*EdevS
        psin=0.5d0*bk*trEn**2

       elseif (kflagD.eq.2) then ! Miehe et al.
        call SPRINC(stran,PS,2,3,1)
        trp1=(PS(1)+PS(2)+PS(3)+abs(PS(1)+PS(2)+PS(3)))/2.d0
        trn1=(PS(1)+PS(2)+PS(3)-abs(PS(1)+PS(2)+PS(3)))/2.d0
        trp2=0.d0
        trn2=0.d0
        do i=1,3
         trp2=trp2+(PS(i)+abs(PS(i)))**2.d0/4.d0
         trn2=trn2+(PS(i)-abs(PS(i)))**2.d0/4.d0
        end do
        psip=xnu*eg/(1d0-2d0*xnu)*trp1**2d0+eg*trp2
        psin=xnu*eg/(1d0-2d0*xnu)*trn1**2d0+eg*trn2

       else ! no split
        psip=0.d0
        do i=1,ntens
         psip=psip+stress(i)*stran(i)*0.5d0
        end do
       endif

       if (kFlagF.eq.1) Hmin=3.d0*Gc/(16.d0*xl)
       H=max(psit,psip,Hmin)
      else
       call SPRINC((stress),PS,1,3,1)
	    savg=0.5d0*(PS(1)+PS(2))
	    sdif=0.5d0*(PS(1)-PS(2))
	   sdev=sqrt(sdif*sdif+PS(3)*PS(3))
	   smax=savg+sdev

       psip=0.5d0*max(smax,0.d0)**2.0d0/E
       Hmin=0.5d0*ft**2.0d0/E
       H=max(psit,psip,Hmin)
      end if

      end
c***********************************************************************
      subroutine puck_criteria(props, ntens, stress,statev)
	  include 'aba_param.inc' !implicit real(a-h o-z)

      dimension props(*),stress(ntens),statev(*),flags_old(6)
	  dimension flags(6)
!Get the old flags. 	  
	   flags_old(1:6)=statev(11:16)!11,12,13,14,15
	   flags(:)=0 !initiate the new flags
!	   Properties to be saved as input. All the properties are in Local Ply settings
!     From the elastic properties	
		E11= props(13) 
		E22=props(14) 
		E12=props(15)
		xnu12= props(16)
		xnu23=props(17)
			
!properties of strength material 
		r1t=props(19) 
		r2t=props(20) 
		r1c=props(21) 
		r2c=props(22) 
		r12=props(23) 
!CFRP Properties for Puck 
		p21_plus= props(24) 
		p21_minus=props(25) 
		p22_plus=props(26) 
		p22_minus=props(27)
		E11f= props(28)
		v12f=props(29)  
		mag=props(30) 
!Compressive fiber failure shear influence
		alpha= props(31) !0.0d0

!damage activation function
			ft=1
			fc=1
			mt=1
			mc=1

		r22a=r12/(2.0d0*p21_minus)*(dsqrt(1.d0+2.d0*p21_minus*r2c/r12)-1.0d0)
		r12c=r12*dsqrt(1.d0+2.d0*p22_minus)

		
!maxwell-bettis law (minor)
		xnu21=xnu12*E22/E11
			
!Initialize mterial effort
			f_tf=0.0d0
			f_cf=0.0d0
			f_tm=0.0d0
			f_cm=0.0d0
!new flags initialisation
			flag_f=0.0d0 !flags_old(1)
			flag_if=0.d0 !flags_old(2)

!Evaluate the PUCK Failure criteria
!for fiber failure 

!if its already activated, then skip the steps

		if (flags_old(1).eq.1) then 
		flags(1)=1 !fiber is activated already
		!Write(7,*)'fiber is already activated'
		goto 10
		endif


! via f_tf directly 
	  if (stress(1).gt.0.0d0) then  !tension
	    f_tf=stress(1)/r1t     
	    flags(1)=f_tf !fiber failure in tension
		
		   
	  elseif(stress(1).lt.0.0d0) then  !compression
	    f_cf=dsqrt((stress(1)/r1c)**2.0d0+ alpha*(stress(3)/r12)**2.0d0)
		flags(1)=f_cf !fiber failure in compression  	
	  endif
		
10          continue

!Inter-fiber failure 
!check if its already active 
		
	    if (flags_old(2).eq.1) then 
		flags(2)=1 !interfiber is activated already
		goto 20
		endif
!Mode-A (Transverse tension)

		if (stress(2).gt.0.0d0) then 
		f_tm=dsqrt((stress(3)/r12)**2.0d0 + 
     #  (1.0d0-p21_plus*r2t/r12)**2.0d0 *(stress(2)/r2t)**2.0d0)+
     #   p21_plus*stress(2)/r12
	    flags(2)=f_tm
		
			if(f_tm.ge.1.0d0) then 
			flags(2)=1 !
			!Write(7,*)'Matrix is already activated'
			endif	 
		endif
		
		if (stress(2).lt.0.0d0) then  
		!mode-B (transverse compression)
		fac=abs(stress(2)/stress(3))
		fac1=abs(r22a/r12c)
		
		if ((fac.ge.0.0d0).and.(fac.le.fac1)) then 
			f_cm=dsqrt((stress(3)/r12)**2.0d0 + (p21_minus/r12*stress(2))**2.0d0)
     #    + p21_minus/r12*stress(2)
	     flags(2)=f_cm
				 
    	    if(f_cm.ge.1.d0) then 
				flags(2)=1
			endif
! Mode C (Transverse compression)
		elseif (((1/fac).ge.0.0d0 ).and.((1/fac).le.(1/fac1))) then  	
			f_cm=((stress(3)/(2.0d0*(1.0d0+p22_minus)*r12))**2.0d0 +
     #      (stress(2)/r2c)**2.0d0)*r2c/(-stress(2)) 
				
		  	
				flags(2)=f_cm
			if(f_cm.ge.1.0d0) then 
			flags(2)=1
			endif		
		endif 
	  endif

20          continue		
		
		statev(11:16)=flags(1:6) ! save the new flags
		end 
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
      subroutine kmatinv3(drot,drotInv) 
	  include 'aba_param.inc' !implicit real(a-h o-z)

      dimension drot(3,3), drotInv(3,3)	
	  			
			detinv=0.0d0
	
    ! Calculate the inverse determinant of the matrix
       det1=drot(1,1)*drot(2,2)*drot(3,3)-drot(1,1)*drot(2,3)*drot(3,2)  
       det2=-drot(1,2)*drot(2,1)*drot(3,3)+drot(1,2)*drot(3,1)*drot(2,3)
       det3=drot(1,3)*drot(2,1)*drot(3,2)-drot(1,3)*drot(3,1)*drot(2,2)
	 
		det=det1+det2+det3
		
			det=1/det
			
    ! Calculate the inverse of the matrix
       drotInv(1,1) = +det * (drot(2,2)*drot(3,3) - drot(2,3)*drot(3,2))
       drotInv(2,1) = -det * (drot(2,1)*drot(3,3) - drot(2,3)*drot(3,1))
       drotInv(3,1) = +det * (drot(2,1)*drot(3,2) - drot(2,2)*drot(3,1))
       drotInv(1,2) = -det * (drot(1,2)*drot(3,3) - drot(1,3)*drot(3,2))
       drotInv(2,2) = +det * (drot(1,1)*drot(3,3) - drot(1,3)*drot(3,1))
       drotInv(3,2) = -det * (drot(1,1)*drot(3,2) - drot(1,2)*drot(3,1))
       drotInv(1,3) = +det * (drot(1,2)*drot(2,3) - drot(1,3)*drot(2,2))
       drotInv(2,3) = -det * (drot(1,1)*drot(2,3) - drot(1,3)*drot(2,1))
       drotInv(3,3) = +det * (drot(1,1)*drot(2,2) - drot(1,2)*drot(2,1))
	   !Write(7,*)'Matrix drotInv inside loop',drotInv
	   
      end 

