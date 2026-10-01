! Developing a 3D version of the velocity verlet.
! Hopefully will be parallelized.
! ToDo:
!   1. Update allocateXVF() to 3D.
!   2. Include Lennard-Jones force field.
!   3. Update getForces() to 3D.
!       * Keep in mind that this is the main cycle. So, must be parallelized.
!         Additionally, several nodes can access the same force 
!         element when using third Newton's law, this must be avoided.
!   4. Update initPositions() to 3D.
!   5. Include initSquareLatticePositions(). 
!       * This is fundamental for avoiding the singularities of Lennar-Jones
!         potential.
!       o For maintaining some randomness in the positions' initialization, a
!         small random displacement around the crystalline point can be 
!         included.
!   6. Update initVelocities() to 3D.
!   7. Update verletStep() to 3D.
!   8. Update storePosition() to 3D.
!   9. Update storeVelocity() to 3D.
!   10. Replace length to volume.
!   11. Check which loops, apart of getForces(), could be paralelized.
!   12. Include a Uniform Random Distribution of velocities for checking the
!       convergence to Maxwell-Boltzmann.
!   13. Include sigma (from Lennard-Jones) to simState
!   14. Include epsilon (from Lennard-Jones) to simState
!
!
! Completed Tasks:
!   1. Done! 2. Done! 3. Done! 4. Done! 5. Done!      
!   6. Done! 7. Done! 8. Done! 9. Done! 10. Done!
!   11.      12. Done!  13. Done!   14. Done!
!
!
! IMPORTANT:
!   -> N = n^3, where n is an integer.
!           * For avoiding problem with initSquareLatticePositions
!
!
!
! NEXT STEPS ToDo
!   1. Add Periodic Boundary Conditions. Done!
!   2. Protect code from division by zero in Lennar-Jones force calculation.
!      Done!
!   3. Add functions for extracting and saving energy = kinetic + potential


! =============================================================================
!                            Auxiliary Modules
! =============================================================================
module simulation
    use, intrinsic :: iso_fortran_env, only: dp => real64   ! Double presicion
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)   ! pi definition

    type :: simState
        integer,  private :: numParticles  ! Number of particles in the system
        integer,  private :: numSteps      ! Number of time steps
        integer,  private :: framePeriod   ! Steps between coordinates storage
        real(dp), private :: volume        ! Volume of the system
        real(dp), private :: mass          ! Particles' mass
        real(dp), private :: dt            ! Size of time step
        real(dp), private :: kb            ! Boltzmann's constant
        real(dp), private :: temp          ! System's temperature
        real(dp), private :: sigma         ! From Lennard-Jones
        real(dp), private :: epsilon0      ! From Lennard-Jones
        real(dp), private :: rCut          ! Max radius to compute force


        character(len=:), allocatable, private          :: dataDirectory
        real(dp), allocatable, private, dimension(:, :) :: position, &
                                                           velocity, &
                                                           force

        contains
            ! Alias the force field for being able to change it independently
            ! procedure, nopass :: forceField => forceLogHarmonic
            procedure :: forceField => forceLennardJones

            ! The main subroutines
            procedure :: initPositions
            procedure :: initSquareLatticePositions
            procedure :: initMaxwellBoltzmannVelocities
            procedure :: initUniformVelocities
            procedure :: getForces
            procedure :: verletStep

            procedure :: storePosition
            procedure :: storeVelocity

            procedure :: allocateXVF

            ! Set subroutines for privates
            procedure :: setNumParticles
            procedure :: setNumSteps
            procedure :: setFramePeriod
            procedure :: setVolume
            procedure :: setMass
            procedure :: setDt
            procedure :: setKb
            procedure :: setTemperature
            procedure :: setSigma
            procedure :: setEpsilon
            procedure :: setRCut
            procedure :: setOutputDir

            ! Get subroutines for privates
            procedure :: getNumParticles
            procedure :: getNumSteps
            procedure :: getFramePeriod
            procedure :: getVolume
            procedure :: getMass
            procedure :: getDt
            procedure :: getKb
            procedure :: getTemperature
            procedure :: getSigma
            procedure :: getRCut
            procedure :: getEpsilon
            ! procedure :: getOutputDir
    end type simState

    contains

    ! =====================
    ! Set and Get params
    ! =====================
    subroutine setNumParticles(self, number)
        class(simState), intent(inout) :: self
        integer,         intent(in)    :: number
        self%numParticles = number
    end subroutine setNumParticles

    subroutine setNumSteps(self, number)
        class(simState), intent(inout) :: self
        integer,         intent(in)    :: number
        self%numSteps = number
    end subroutine setNumSteps

    subroutine setFramePeriod(self, number)
        class(simState), intent(inout) :: self
        integer,         intent(in)    :: number
        self%framePeriod = number
    end subroutine setFramePeriod

    subroutine setVolume(self, number)
        class(simState), intent(inout) :: self
        real(dp),        intent(in)    :: number
        self%volume = number
    end subroutine setVolume

    subroutine setMass(self, number)
        class(simState), intent(inout) :: self
        real(dp),        intent(in)    :: number
        self%mass = number
    end subroutine setMass

    subroutine setDt(self, number)
        class(simState), intent(inout) :: self
        real(dp),        intent(in)    :: number
        self%dt = number
    end subroutine setDt

    subroutine setKb(self, number)
        class(simState), intent(inout) :: self
        real(dp),        intent(in)    :: number
        self%kb = number
    end subroutine setKb

    subroutine setTemperature(self, number)
        class(simState), intent(inout) :: self
        real(dp),        intent(in)    :: number
        self%temp = number
    end subroutine setTemperature

    subroutine setSigma(self, number)
        class(simState), intent(inout) :: self
        real(dp),        intent(in)    :: number
        self%sigma = number
    end subroutine setSigma

    subroutine setEpsilon(self, number)
        class(simState), intent(inout) :: self
        real(dp),        intent(in)    :: number
        self%epsilon0 = number
    end subroutine setEpsilon

    subroutine setRCut(self, number)
        class(simState), intent(inout) :: self
        real(dp),        intent(in)    :: number
        self%rCut = number
    end subroutine setRCut

    subroutine setOutputDir(self, dir)
        class(simState),  intent(inout) :: self
        character(len=*), intent(in)    :: dir
        self%dataDirectory = dir
    end subroutine setOutputDir

    
    ! Get function
    integer function getNumParticles(self)
        class(simState), intent(in) :: self
        getNumParticles = self%numParticles 
    end function getNumParticles

    integer function getNumSteps(self)
        class(simState), intent(in) :: self
        getNumSteps = self%numSteps 
    end function getNumSteps

    integer function getFramePeriod(self)
        class(simState), intent(in) :: self
        getFramePeriod = self%framePeriod
    end function getFramePeriod

    real(dp) function getVolume(self)
        class(simState), intent(in) :: self
        getVolume = self%volume
    end function getVolume

    real(dp) function getMass(self)
        class(simState), intent(in) :: self
        getMass = self%mass
    end function getMass

    real(dp) function getDt(self)
        class(simState), intent(in) :: self
        getDt = self%dt
    end function getDt

    real(dp) function getKb(self)
        class(simState), intent(in) :: self
        getKb = self%kb
    end function getKb

    real(dp) function getTemperature(self)
        class(simState), intent(in) :: self

        ! By equipartition theorem T = (m/(3*Kb*N)) * sum_i(v_i^2)
        getTemperature = self%mass * sum(self%velocity**2) / &
                          (3.0_dp * self%kb * self%numParticles)
    end function getTemperature

    real(dp) function getSigma(self)
        class(simState), intent(in) :: self
        getSigma = self%sigma
    end function getSigma

    real(dp) function getEpsilon(self)
        class(simState), intent(in) :: self
        getEpsilon = self%epsilon0
    end function getEpsilon

    real(dp) function getRCut(self)
        class(simState), intent(in) :: self
        getRCut = self%rCut
    end function getRCut

    ! character(len=:), allocatable function getOutputDir(self)
    !     class(simState),  intent(in) :: self
    !     getOutputDir = self%dataDirectory
    ! end function getOutputDir


    ! Allocate position, velocity and force
    subroutine allocateXVF(self)
        class(simState), intent(inout) :: self
        allocate(self%position(3, self%numParticles))
        allocate(self%velocity(3, self%numParticles))
        allocate(self%force(3, self%numParticles))
    end subroutine allocateXVF

    ! =====================
    ! Force Fields
    ! =====================
    subroutine forceLogHarmonic(x1, x2, f)
        real(dp), intent(in)  :: x1, x2
        real(dp), intent(out) :: f
        real(dp), parameter   :: A = 2_dp
        real(dp), parameter   :: B = 1_dp
        f = B / (x2 - x1) - 2 * A * (x2 - x1)
        f = -f
    end subroutine forceLogHarmonic

    pure subroutine forceLennardJones(self, r1, r2, f)
        class(simState), intent(in)         :: self
        real(dp), intent(in), dimension(3)  :: r1, r2
        real(dp), intent(out), dimension(3) :: f
        real(dp)                            :: r2_val, inv_r2, sigma2
        real(dp)                            :: sr2, sr6, sr12, l_box
        real(dp), dimension(3)              :: dr

        f = 0.0_dp
        l_box = self%volume ** (1.0_dp / 3.0_dp)
        dr = r2 - r1
        dr = dr - l_box * anint(dr / l_box)  ! anint() is nint() but returns real(dp)
        r2_val = sum(dr**2)

        if (r2_val < self%rCut**2 .and. r2_val > 1e-12_dp) then
            sigma2 = self%sigma ** 2
            inv_r2 = 1.0_dp / r2_val
            sr2 = sigma2 * inv_r2
            sr6 = sr2 * sr2 * sr2
            sr12 = sr6 * sr6

            f = - 24.0_dp * self%epsilon0 * inv_r2 * (2.0_dp*sr12 - sr6) * dr
        end if
    end subroutine forceLennardJones

    subroutine getForces(self)
        class(simState), intent(inout) :: self
        real(dp), dimension(3)         :: f
        integer                        :: i, j

        self%force = 0.0_dp
        f = 0.0_dp

        ! The following comments is for parallelization
        ! Newton's third law is sacrified in order to parallelize
        ! default(none) is for explicitely declaring every variable intend
        ! shared(self) is for the shared objects across threads
        ! private(i, j, f) is for the objects that threads handle individually
        ! schedule(static) is for equally workload distribution across threads

        !$omp parallel do default(none)    &
        !$omp             shared(self)     &
        !$omp             private(i, j, f) &
        !$omp             schedule(static)
        do i = 1, self%numParticles
            do j = 1, self%numParticles

                if (j == i) then 
                    cycle
                end if

                call self%forceField(self%position(:, i), self%position(:, j), f)
                self%force(:, i) = self%force(:, i) + f
            end do
        end do
        !$omp end parallel do


    end subroutine getForces


    ! =====================
    ! Random Numbers
    ! =====================
    subroutine initRandomSeed()
        ! This is for using a different seed in every simulation
        integer, allocatable, dimension(:) :: seed_array
        integer                            :: seed_size

        call random_seed(size=seed_size)
        allocate(seed_array(seed_size))

        call system_clock(count=seed_array(1))
        seed_array(2:) = seed_array(1) * 13  ! Arbitrary scramble
        call random_seed(put=seed_array)

        ! Sweeping away
        deallocate(seed_array)
    end subroutine initRandomSeed


    subroutine initPositions(self)
        class(simState), intent(inout) :: self
        real(dp)                       :: l_box

        l_box = self%volume**(1.0_dp/3.0_dp)

        call random_number(self%position)
        self%position = (self%position - 0.5_dp) * l_box
    end subroutine initPositions


    subroutine initSquareLatticePositions(self)
        ! For avoiding the potential's singularities,
        ! A small random displacement can be included.
        class(simState), intent(inout) :: self
        real(dp)                       :: gap  ! Avoid particle at boundary
        real(dp)                       :: l, dx
        real(dp), parameter            :: tol = 1e-9_dp
        integer                        :: i, j, k, n, count

        ! Calculate size of the comptational domain
        l = self%volume**(1.0_dp/3.0_dp)
        gap = 0.1_dp * l

        ! Calculate number of nodes of the lattice per side
        n = ceiling(real(self%numParticles, dp)**(1.0_dp/3.0_dp) - tol)

        ! Caluculate lattice parameter and include case numParticles = 1
        if (n > 1) then
            dx = (l - 2.0_dp*gap) / (n - 1.0_dp)
        else 
            dx = 0_dp
        end if

        count = 1
        lattice_loop: do i = 1, n
            do j = 1, n
                do k = 1, n
                    self%position(1, count) =  -0.5_dp*l + gap + dx*(i-1)
                    self%position(2, count) =  -0.5_dp*l + gap + dx*(j-1)
                    self%position(3, count) =  -0.5_dp*l + gap + dx*(k-1)
                    count = count + 1

                    if (count > self%numParticles) exit lattice_loop
                end do
            end do
        end do lattice_loop
    end subroutine initSquareLatticePositions


    subroutine initMaxwellBoltzmannVelocities(self)
        ! Using Maxwell-Boltzmann distribution of velocities
        class(simState), intent(inout)               :: self
        real(dp)                                     :: std, temp1, alpha
        real(dp), dimension(3, self%numParticles, 2) :: aux
        real(dp), dimension(3)       :: v_cm

        std = sqrt(self%kb * self%temp / self%mass)

        call random_number(aux)

        ! Applying Box-Muller transform
        self%velocity = std * sqrt(-2_dp * &
                                   log(1.0_dp - aux(:, :, 1))) * &
                              cos(2_dp * pi * aux(:, :, 2))

        ! Set Center-of-Mass velocity to zero
        v_cm = sum(self%velocity, dim=2) / real(self%numParticles, dp)
        self%velocity = self%velocity - &
                        spread(v_cm, dim=2, ncopies=self%numParticles)

        ! Previous step might cause a temperature reduction. To avoid this:
        ! By equipartition theorem T = (m/(3*Kb*N)) * sum_i(v_i^2)
        temp1 = self%getTemperature()

        if (temp1 /= self%temp) then
            alpha = sqrt(self%temp / temp1)
            self%velocity = self%velocity * alpha
        end if

    end subroutine initMaxwellBoltzmannVelocities


    subroutine initUniformVelocities(self)
        class(simState), intent(inout) :: self
        real(dp)                       :: temp1, alpha
        real(dp), dimension(3)         :: v_cm

        call random_number(self%velocity)

        v_cm = sum(self%velocity, dim=2) / real(self%numParticles, dp)
        self%velocity = self%velocity - &
                        spread(v_cm, dim=2, ncopies=self%numParticles)

        temp1 = self%getTemperature()

        if (temp1 /= self%temp) then
            alpha = sqrt(self%temp / temp1)
            self%velocity = self%velocity * alpha
        end if
    end subroutine initUniformVelocities

    subroutine applyPBS(self)
        class(simState), intent(inout) :: self
        real(dp)                       :: l_box

        l_box = self%volume ** (1.0_dp / 3.0_dp)
        self%position = self%position - l_box * anint(self%position/l_box)
    end subroutine


    ! =====================
    ! Time evolution
    ! =====================
    subroutine verletStep(self)
        class(simState), intent(inout)                    :: self
        real(dp),        dimension(3, self%numParticles)  :: v_aux

        v_aux = self%velocity + self%dt * self%force * 0.5_dp / self%mass
        self%position = self%position + self%dt * v_aux

        call applyPBS(self)
        call getForces(self)

        self%velocity = v_aux + self%dt * self%force * 0.5_dp / self%mass

    end subroutine verletStep


    ! =====================
    ! Storing data
    ! =====================
    subroutine storePosition(self, i)
        class(simState),  intent(in) :: self
        integer,          intent(in) :: i
        integer,          save       :: io_unit
        integer                      :: status_code, j


        if (i == 1) then
            open(newunit  = io_unit,                               &
                 file     = self%dataDirectory // "position.xyz",  &
                 status   = "replace",                             &
                 action   = "write",                               &
                 iostat   = status_code)
                
                if (status_code /= 0) then
                    print *, "ERROR: Could not open file: ", & 
                          self%dataDirectory // "position.xyz"
                    print *, "Does the target directory exist?"
                    stop 1
                end if
        end if

        if (mod(i, self%framePeriod) == 0) then

            write(io_unit, '(I8)') self%numParticles

            write(io_unit, '(A, I8, A, F12.6)') &
                  "Lennard-Jones MD Frame | Step = ", i, &
                  " | Time = ", real(i, dp) * self%dt

            do j = 1, self%numParticles
                write(io_unit, '(A, 3(1X, F12.6))') "H", self%position(1, j), &
                                                         self%position(2, j), &
                                                         self%position(3, j)
            end do
        end if

        if (i == self%numSteps) then
            close(io_unit)
        end if
    end subroutine storePosition

    subroutine storeVelocity(self, i)
        class(simState),  intent(in) :: self
        integer,          intent(in) :: i
        integer,          save       :: io_unit
        integer                      :: status_code, j


        if (i == 1) then
            open(newunit  = io_unit,                               &
                 file     = self%dataDirectory // "velocity.txt",  &
                 status   = "replace",                             &
                 action   = "write",                               &
                 iostat   = status_code)
                
                if (status_code /= 0) then
                    print *, "ERROR: Could not open file: ", self%dataDirectory // "velocity.txt"
                    print *, "Does the target directory exist?"
                    stop 1
                end if
        end if

        if (mod(i, self%framePeriod) == 0) then
            write(io_unit, '(A, I8, A, F12.6)') &
                  "Lennard-Jones MD Frame | Step = ", i, &
                  " | Time = ", real(i, dp) * self%dt

            do j = 1, self%numParticles
                write(io_unit, '(I12, 3(1X, F12.6))') j, self%velocity(1, j), &
                                                         self%velocity(2, j), &
                                                         self%velocity(3, j)
            end do
        end if

        if (i == self%numSteps) then
            close(io_unit)
        end if
    end subroutine storeVelocity
end module


! =============================================================================
!                               Main Program
! =============================================================================
program main
    use simulation
    implicit none

    ! =====================
    ! Declarations
    ! =====================
    type(simState)              :: state
    integer                     :: i
    integer                     :: numParticles, numSteps, framePeriod
    real(dp)                    :: volume, mass, dt, kb, temp, sigma, epsilon0
    real(dp)                    :: rCut
    character(len=*), parameter :: dataDirectory = "../postpro/"


    ! =====================
    ! Setting up parameters
    ! =====================
    numParticles  = 512            ! Total number of particles
    numSteps      = 5000           ! Total number of time steps
    framePeriod   = 1              ! Steps between Coordinates storage
    volume        = 100.0_dp       ! Volume of the 1D system
    mass          = 1.0_dp         ! Particles' mass
    dt            = 1e-5_dp        ! Size of time step
    kb            = 1.0_dp         ! Boltzmann's constant
    temp          = 0.85_dp        ! Temperature of the system
    sigma         = 1.0_dp         ! For Lennard-Jones
    epsilon0      = 1.0_dp         ! For Lennard-Jones
    rCut          = 2.5_dp         ! For truncating far interactions


    call state%setOutputDir(dataDirectory)
    call state%setNumParticles(numParticles)
    call state%setFramePeriod(framePeriod)
    call state%setNumSteps(numSteps)
    call state%setVolume(volume)
    call state%setMass(mass)
    call state%setDt(dt)
    call state%setKb(kb)
    call state%setTemperature(temp)
    call state%setSigma(sigma)
    call state%setEpsilon(epsilon0)
    call state%setRCut(rCut)


    ! =====================
    ! Initializing
    ! =====================
    call state%allocateXVF()
    
    call initRandomSeed()
    call state%initSquareLatticePositions()
    call state%initUniformVelocities()
    call state%getForces()


    ! =====================
    ! Main Loop
    ! =====================
    do i = 1, numSteps
        call state%verletStep()
        call state%storePosition(i)
        call state%storeVelocity(i)
    end do
end program main
