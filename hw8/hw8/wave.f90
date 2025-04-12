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

  integer ispec,i,j,iglob,itime

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

! temperature and temperature time derivative
  double precision, dimension(NGLOB) :: displ,veloc,accel

! local to global numbering
  integer, dimension(NGLL,NSPEC) :: ibool

! time marching
  double precision deltat,deltatover2,deltatsquareover2
  double precision dh,time_step

! end gradients
  double precision grad_1,grad_NGLOB

! end temperatures
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
! put your codes here
!------------------------------------------------------------------------------

! initialize
!  displ(:) = 
!  veloc(:) = 
!  accel(:) = 

  do itime = 1,NSTEP

! "predictor" displacement, velocity, initialize acceleration
!------------------------------------------------------------------------------
! put your codes here
!------------------------------------------------------------------------------

! boundary conditions
!------------------------------------------------------------------------------
! put your codes here
!------------------------------------------------------------------------------

! solver
!------------------------------------------------------------------------------
! put your codes here
!------------------------------------------------------------------------------

! "corrector" acceleration, velocity, displacement
!------------------------------------------------------------------------------
! put your codes here
!------------------------------------------------------------------------------

! write out snapshots
    if(mod(itime-1,25) == 0) then
      write(moviefile,'(a,i5.5)')'snapshot',itime
      open(unit=10,file=moviefile,status='unknown')
      do iglob = 1,NGLOB
        write(10,*) sngl(x(iglob)),sngl(displ(iglob))
      enddo
      close(10)
    endif

  enddo ! end time loop

  end program wave
