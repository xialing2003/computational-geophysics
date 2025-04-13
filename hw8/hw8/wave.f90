  program wave

  implicit none

  include "constants.h"

! number of timesteps
  integer, parameter :: NSTEP = 400
 
! time step in seconds
  double precision, parameter :: DT = 0.25
 
! fixed boundary conditions
  logical, parameter :: FIXED_BC = .true.
 
! model parameters  (SI)
  double precision, parameter :: LENGTH = 100
  double precision, parameter :: DENSITY = 1
  double precision, parameter :: SHEARMODULUS = 1

  integer ispec,i,j,iglob,itime, gamma,jglob

! Gauss-Lobatto-Legendre points of integration
  double precision, dimension(NGLL) :: xigll
 
! weights
  double precision, dimension(NGLL) :: wgll
 
! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLL,NGLL) :: hprime

! anchors
  double precision, dimension(NSPEC) :: x1,x2

! global grid points
  double precision, dimension(NGLOB) :: x

! material properties
  double precision, dimension(NGLL,NSPEC) :: rho,shear

! Jacobian `matrix' and Jacobian
  double precision, dimension(NGLL,NSPEC) :: dxidx,jacobian

! local mass matrix
  double precision mass_local

! global mass matrix
  double precision, dimension(NGLOB) :: mass_global

! local stiffness matrix
  double precision, dimension(NGLL, NGLL) :: stiffness_local
  double precision stiffness_temp

! displacement and displacement time derivative
  double precision, dimension(NGLOB) :: displ,veloc,accel

! global rhs matrix
  double precision, dimension(NGLOB) :: rhs_global

! local rhs matrix
  double precision rhs_local

! local to global numbering
  integer, dimension(NGLL,NSPEC) :: ibool

! time marching
  double precision deltat,deltatover2,deltatsquareover2
  double precision dh,time_step

! end gradients
  double precision grad_1,grad_NGLOB

! end displacements
  double precision displ_1,displ_NGLOB

! movie
  character(len=50) moviefile

!++++++++++++++++++++++++++++++++++++++++++++++++++

  call define_derivative_matrix(xigll,wgll,hprime)

! evenly spaced anchors between 0 and 1
  do ispec = 1,NSPEC
    x1(ispec) = LENGTH*dble(ispec-1)/dble(NSPEC)
    x2(ispec) = LENGTH*dble(ispec)/dble(NSPEC)
  enddo

! set up the mesh properties
  do ispec = 1,NSPEC
    do i = 1,NGLL
      rho(i,ispec) = DENSITY
      shear(i,ispec) = SHEARMODULUS
      dxidx(i,ispec) = 2. / (x2(ispec)-x1(ispec))
      jacobian(i,ispec) = (x2(ispec)-x1(ispec)) / 2.
    enddo
  enddo

! set up local to global numbering
  iglob = 1
  do ispec = 1,NSPEC
    do i = 1,NGLL
      if(i > 1) iglob = iglob+1
      ibool(i,ispec) = iglob
    enddo
  enddo

! get the global grid points
  do ispec = 1,NSPEC
    do i = 1,NGLL
      iglob = ibool(i,ispec)
      x(iglob) = 0.5*(1.-xigll(i))*x1(ispec)+0.5*(1.+xigll(i))*x2(ispec)
    enddo
  enddo

! calculate the global mass matrix 'mass_global'
!------------------------------------------------------------------------------
  mass_global(:) = 0
  do ispec = 1,NSPEC
    do i = 1,NGLL
      ! local to global index
      iglob = ibool(i,ispec)
      ! local mass matrix
      mass_local = wgll(i)*rho(i,ispec)*jacobian(i,ispec)
      ! assemble global mass matrix
      mass_global(iglob) = mass_global(iglob)+mass_local
    end do
  end do
!------------------------------------------------------------------------------

! estimate the time step 'time_step' DT, DT/2, DT^2/2
!------------------------------------------------------------------------------
  deltat = DT
  deltatover2 = DT/2
  deltatsquareover2 = DT*DT/2
!------------------------------------------------------------------------------
  print *,'time step estimate: ',deltat,' seconds'

! set up the initial and boundary conditions
!------------------------------------------------------------------------------
  do ispec = 1,NSPEC
    do i = 1,NGLL
      iglob = ibool(i,ispec)
      displ(iglob) = exp(-(x(iglob) - 50)**2 * 0.1)
      veloc(iglob) = 0
      accel(iglob) = 0
    end do
  end do

  displ_1 = 0
  displ_NGLOB = 0

  write(moviefile,'(a,i5.5)')'snapshot',0
  open(unit=10,file=moviefile,status='unknown')
  do iglob = 1,NGLOB
    write(10,*) sngl(x(iglob)),sngl(displ(iglob))
  enddo
  close(10)
!------------------------------------------------------------------------------

! initialize
!  displ(:) = exp[−(x − 50)^2 ∗ 0.1]
!  veloc(:) = 
!  accel(:) = 

  do itime = 1,NSTEP

! "predictor" displacement, velocity, initialize acceleration
!------------------------------------------------------------------------------
    displ(:) = displ(:) + deltat*veloc(:) + deltatsquareover2*accel(:)
    veloc(:) = veloc(:) + deltatover2*accel(:)
    rhs_global(:) = 0
    if (FIXED_BC) then
      displ(1) = displ_1
      displ(NGLOB) = displ_NGLOB
    endif
!------------------------------------------------------------------------------
    do ispec = 1,NSPEC

      ! calculate the local stiffness matrix 'stiffness_local'
      do i = 1,NGLL
        do j = 1,NGLL
          stiffness_temp = 0.0
          do gamma = 1,NGLL
            stiffness_temp = stiffness_temp + wgll(gamma)*shear(gamma,ispec)/jacobian(gamma,ispec)*hprime(i,gamma)*hprime(j,gamma)
          enddo
          stiffness_local(i,j) = stiffness_temp
        enddo
      enddo

      ! GLL point loop
      do i = 1,NGLL
        iglob = ibool(i,ispec)
        rhs_local = 0.0
        do j = 1,NGLL
          jglob = ibool(j,ispec)
          rhs_local = rhs_local - stiffness_local(i,j)*displ(jglob)
        enddo

        ! assembly
        rhs_global(iglob) = rhs_global(iglob) + rhs_local
      enddo
    enddo
! boundary conditions
!------------------------------------------------------------------------------
    if (FIXED_BC) then
      do j = 1,NGLL
        jglob = ibool(j,1)
        rhs_global(1) = rhs_global(1) - shear(1,1)*displ(jglob)*hprime(j,i)/jacobian(j,1)
      enddo
      do j = 1,NGLL
        jglob = ibool(j,NSPEC)
        rhs_global(NGLOB) = rhs_global(NGLOB) - shear(NGLL,NSPEC)*displ(jglob)*hprime(j,i)/jacobian(j,NSPEC)
      enddo
    endif
!------------------------------------------------------------------------------

! solver
!------------------------------------------------------------------------------
    accel(:) = rhs_global / mass_global
!------------------------------------------------------------------------------

! "corrector" acceleration, velocity, displacement
!------------------------------------------------------------------------------
    veloc(:) = veloc(:) + deltatover2*accel(:)
!------------------------------------------------------------------------------

! write out snapshots
    if(mod(itime,20) == 0) then
      write(moviefile,'(a,i5.5)')'snapshot',itime
      open(unit=10,file=moviefile,status='unknown')
      do iglob = 1,NGLOB
        write(10,*) sngl(x(iglob)),sngl(displ(iglob))
      enddo
      close(10)
    endif

  enddo ! end time loop

  end program wave
