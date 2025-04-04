  program diffusion

  implicit none

  include "constants.h"

! number of spectral elements
  integer, parameter :: NSPEC_1 = 3
  integer, parameter :: NSPEC_2 = 9

! number of timesteps
  integer, parameter :: NSTEP = 500000
 
! time step in seconds
  double precision, parameter :: DT = 10000000. ! s
 
! fixed boundary conditions
  logical, parameter :: FIXED_BC = .true.
 
! model parameters  (SI)
  double precision, parameter :: LENGTH = 3.0d+03 ! m
  double precision, parameter :: DENSITY = 2.5d+03 ! kg/m^3
  double precision, parameter :: THERMALCONDUCTIVITY_1 = 10.0d-01 ! cal/m/s/K
  double precision, parameter :: THERMALCONDUCTIVITY_2 = 10.0d-01 ! cal/m/s/K
  double precision, parameter :: HEATCAPACITY = 0.3d+03 ! cal/kg/K
  double precision, parameter :: TEMPERATURE_0 = 0
  double precision, parameter :: TEMPERATURE_L = 10

  integer ispec,i,j,iglob,itime, gamma,jglob

! Gauss-Lobatto-Legendre points of integration
  double precision, dimension(NGLL) :: xigll
 
! weights
  double precision, dimension(NGLL) :: wgll
 
! array with derivatives of Lagrange polynomials
! hprime(i,j) = h'_i(xigll_j) by definition of the derivative matrix
  double precision, dimension(NGLL,NGLL) :: hprime

! anchors
  double precision, dimension(NSPEC) :: x1,x2

! global grid points
  double precision, dimension(NGLOB) :: x

! material properties
  double precision, dimension(NGLL,NSPEC) :: rho,heat_capacity,thermal_conductivity

! Jacobian `matrix' and Jacobian
  double precision, dimension(NGLL,NSPEC) :: dxidx,jacobian

! local mass matrix
  double precision mass_local

! global mass matrix
  double precision, dimension(NGLOB) :: mass_global

! global rhs matrix
  double precision, dimension(NGLOB) :: rhs_global

! local rhs matrix
  double precision rhs_local

! local stiffness matrix
  double precision, dimension(NGLL, NGLL) :: stiffness_local
  double precision stiffness_temp

! temperature and temperature time derivative
  double precision, dimension(NGLOB) :: temperature,dtemperature_dt

! local to global numbering
  integer, dimension(NGLL,NSPEC) :: ibool

! time marching
  double precision deltat,deltatover2
  double precision dh,diffusivity,time_step

! end fluxes
  double precision flux_1,flux_NGLOB

! end temperatures
  double precision temperature_1,temperature_NGLOB

! derivatives
  double precision dtdx,flux,templ,temp(NGLL)

! movie
  character(len=50) moviefile

!++++++++++++++++++++++++++++++++++++++++++++++++++

  call define_derivative_matrix(xigll,wgll,hprime)

! evenly spaced anchors between 0 and 1
  do ispec = 1,NSPEC_1
    x1(ispec) = LENGTH/2*dble(ispec-1)/dble(NSPEC_1)
    x2(ispec) = LENGTH/2*dble(ispec)/dble(NSPEC_1)
  enddo
  do ispec = 1,NSPEC_2
    x1(NSPEC_1+ispec) = LENGTH/2*(1+dble(ispec-1)/dble(NSPEC_2))
    x2(NSPEC_1+ispec) = LENGTH/2*(1+dble(ispec)/dble(NSPEC_2))
  enddo

! set up the mesh properties
  do ispec = 1,NSPEC
    do i = 1,NGLL
      rho(i,ispec) = DENSITY
      thermal_conductivity(i,ispec) = THERMALCONDUCTIVITY_1
      heat_capacity(i,ispec) = HEATCAPACITY
      dxidx(i,ispec) = 2. / (x2(ispec)-x1(ispec))
      jacobian(i,ispec) = (x2(ispec)-x1(ispec)) / 2.
    enddo
  enddo
  do ispec = NSPEC_1+1,NSPEC
    do i = 1,NGLL
      thermal_conductivity(i,ispec) = THERMALCONDUCTIVITY_2
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
! put your codes here
  mass_global(:) = 0
  do ispec = 1,NSPEC
    do i = 1,NGLL
      ! local to global index
      iglob = ibool(i,ispec)
      ! local mass matrix
      mass_local = wgll(i)*rho(i,ispec)*heat_capacity(i,ispec)*jacobian(i,ispec)
      ! assemble global mass matrix
      mass_global(iglob) = mass_global(iglob)+mass_local
    end do
  end do

! estimate the time step 'time_step'
! put your codes here
  deltat = DT
  deltatover2 = DT/2
  print *,'time step estimate: ',deltatover2,' seconds'

! set up the boundary conditions
! put your codes here
  temperature_1 = TEMPERATURE_0
  temperature_NGLOB = TEMPERATURE_L

! initialize
  temperature(:) = 0.
  dtemperature_dt(:) = 0.

  do itime = 1,NSTEP

  ! update temperature
  ! put your codes here
    temperature(:) = temperature(:) + deltatover2*dtemperature_dt

    rhs_global(:) = 0

    !element loop
    do ispec = 1,NSPEC

      ! calculate the local stiffness matrix 'stiffness_local'
      do i = 1,NGLL
        do j = 1,NGLL

        stiffness_temp = 0.0
          do gamma = 1,NGLL
            stiffness_temp = stiffness_temp + wgll(gamma)*thermal_conductivity(gamma,ispec) &
                              /jacobian(gamma,ispec)*hprime(i,gamma)*hprime(j,gamma)
          enddo
        stiffness_local(i,j) = stiffness_temp
        enddo
      enddo

      ! GLL point loop
      do i = 1,NGLL
        iglob = ibool(i,ispec)

        ! local contribution
        rhs_local = 0.0
        do j = 1,NGLL
          jglob = ibool(j,ispec)
          rhs_local = rhs_local - stiffness_local(i,j)*temperature(jglob)
        enddo

        ! boundary conditions
        if (FIXED_BC) then 
          if (ispec == 1 .and. i == 1) then
            do j = 1,NGLL
              jglob = ibool(j,ispec)
              rhs_local = rhs_local - thermal_conductivity(i,ispec)*temperature(jglob)*hprime(j,i)/jacobian(j,ispec)
            enddo
          else if (ispec == NSPEC .and. i == NGLL) then
            do j = 1,NGLL
              jglob = ibool(j,ispec)
              rhs_local = rhs_local + thermal_conductivity(i,ispec)*temperature(jglob)*hprime(j,i)/jacobian(j,ispec)
            enddo
          endif
        endif

        ! assembly
        rhs_global(iglob) = rhs_global(iglob) + rhs_local
      enddo
    enddo

    ! temperature increment
    dtemperature_dt(:) = rhs_global / mass_global

    ! time scheme: corrector term
    temperature(:) = temperature(:) + deltatover2*dtemperature_dt(:)

    temperature(1) = temperature_1
    temperature(NGLOB) = temperature_NGLOB

! write out snapshots
    if(mod(itime-1,20000) == 0) then
      write(moviefile,10) itime
10    format('snapshot',i6.6)
      open(unit=10,file=moviefile,status='unknown')
      do iglob = 1,NGLOB
        write(10,*) sngl(x(iglob)),sngl(temperature(iglob))
      enddo
      close(10)
    endif

  enddo ! end time loop

  end program diffusion
