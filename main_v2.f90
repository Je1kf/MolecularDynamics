! =============================================================================
! For compilation
!                   > gfortran -O3 -fopenmp your_files.f90 -o simulation
! For execution
!                   > export OMP_NUM_THREADS=8
!                   > /usr/bin/time -v ./simulation
!
! -O3 in gfortran call is for asking to the compiler for maximum optimization.
! =============================================================================
!
!
!Developing a 3D version of the velocity verlet.
! Hopefully will be parallelized.
!
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
!   11. On progress 12. Done!  13. Done!   14. Done!
!
!
! IMPORTANT:
!   -> N = n^3, where n is an integer.
!           * For avoiding problem with initSquareLatticePositions
!
!
!
! NEXT STEPS: Done!
!   1. Add Periodic Boundary Conditions. Done!
!   2. Protect code from division by zero in Lennar-Jones force calculation.
!      Done!
!   3. Add functions for extracting and saving energy = kinetic + potential
! 
!
!
! Further steps ToDo!
!   -> Extract Statistics:
!           1. Potential Energy.                            STATE: Done!
!           2. Kinetic Energy.                              STATE: Done!
!               2.1. Update getTemperature function.        STATE: Done!
!           3. Internal Energy.                             STATE: Done! 
!           4. Preassure.                                   STATE: Done!
!           5. Correlation Function. (inside forces loop)   STATE:
!           6. Update temp and press in verletStep          STATE: Better, not.
!
!   -> Include statistics outputs in the output files.      STATE: Done!
!   -> Create output file for the correlation function.     STATE: 


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
        integer,  private :: nbins         ! Number of bins in corr. function
        real(dp), private :: volume        ! Volume of the system
        real(dp), private :: mass          ! Particles' mass
        real(dp), private :: dt            ! Size of time step
        real(dp), private :: kb            ! Boltzmann's constant
        real(dp), private :: temp          ! System's temperature
        real(dp), private :: sigma         ! From Lennard-Jones
        real(dp), private :: epsilon0      ! From Lennard-Jones
        real(dp), private :: rCut          ! Max radius to compute force
        real(dp), private :: U_rCut        ! Potential energy at rCut
        real(dp), private :: press         ! System's Preassure
        real(dp), private :: kEnergy       ! Kinetic energy of the system
        real(dp), private :: uEnergy       ! Potential energy of the system


        character(len=:), allocatable, private          :: dataDirectory
        integer,  allocatable, private, dimension(:)    :: rCount
        real(dp), allocatable, private, dimension(:, :) :: position, &
                                                           velocity, &
                                                           force

        contains
            ! Alias the force field for being able to change it independently
            ! procedure, nopass :: forceField => forceLogHarmonic
            procedure :: forceField     => forceLennardJones
            procedure :: potentialField => potentialLennardJones

            ! The main subroutines
            procedure :: initPositions
            procedure :: initSquareLatticePositions
            procedure :: initMaxwellBoltzmannVelocities
            procedure :: initUniformVelocities
            procedure :: getForces
            procedure :: verletStep
            procedure :: getTotalPotentialEnergy
            procedure :: getTotalKineticEnergy
            procedure :: getPreassure

            procedure :: storeData

            procedure :: allocateXVF
            procedure :: deallocateXVF

            ! Set subroutines for privates
            procedure :: setNumParticles
            procedure :: setNumSteps
            procedure :: setFramePeriod
            procedure :: setNbins
            procedure :: setVolume
            procedure :: setMass
            procedure :: setDt
            procedure :: setKb
            procedure :: setTemperature
            procedure :: setSigma
            procedure :: setEpsilon
            procedure :: setRCut
            procedure :: setOutputDir
            procedure :: setKEnergy
            procedure :: setKEnergy0

            ! Get subroutines for privates
            procedure :: getNumParticles
            procedure :: getNumSteps
            procedure :: getFramePeriod
            procedure :: getNbins
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

    subroutine setNbins(self, number)
        class(simState), intent(inout) :: self
        integer,         intent(in)    :: number
        self%nbins = number
    end subroutine setNbins

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
        real(dp)                       :: sr

        self%rCut = number
        sr = self%sigma / number

        self%U_rCut = 4.0_dp * self%epsilon0 * (sr**12 - sr**6)
    end subroutine setRCut

    subroutine setOutputDir(self, dir)
        class(simState),  intent(inout) :: self
        character(len=*), intent(in)    :: dir
        self%dataDirectory = dir
    end subroutine setOutputDir

    subroutine setKEnergy(self)
        class(simState), intent(inout) :: self

        self%kEnergy = self%getTotalKineticEnergy()
    end subroutine setKEnergy

    subroutine setKEnergy0(self)
        class(simState), intent(inout) :: self

        self%kEnergy = 1.5_dp * self%numParticles * self%kb * self%temp
    end subroutine setKEnergy0

    
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

    integer function getNbins(self)
        class(simState), intent(in) :: self
        getNbins = self%nbins
    end function getNbins

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

        ! By equipartition theorem T = 2*K/(3*Kb*N)
        getTemperature = 2.0_dp * self%kEnergy / &
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
        allocate(self%rCount(self%nbins))

        self%rCount = 0
    end subroutine allocateXVF

    subroutine deallocateXVF(self)
        class(simState), intent(inout) :: self
        deallocate(self%position)
        deallocate(self%velocity)
        deallocate(self%force)
        deallocate(self%rCount)
    end subroutine deallocateXVF

    ! =====================
    ! Potentials
    ! =====================
    pure subroutine potentialLennardJones(self, r1, r2, u)
        class(simState), intent(in)        :: self
        real(dp), intent(in), dimension(3) :: r1, r2
        real(dp), intent(out)              :: u
        real(dp), dimension(3)             :: dr
        real(dp)                           :: l_box, r2_val, sigma2, &
                                              inv_r2, sr2, sr6, sr12

        u = 0.0_dp
        l_box = self%volume ** (1.0_dp / 3.0_dp)
        dr = r2 - r1
        ! anint() is nint() but returns real(dp)
        dr = dr - l_box * anint(dr / l_box)
        r2_val = sum(dr**2)

        if (r2_val < self%rCut**2 .and. r2_val > 1e-12_dp) then
            sigma2 = self%sigma ** 2
            inv_r2 = 1.0_dp / r2_val
            sr2 = sigma2 * inv_r2
            sr6 = sr2 * sr2 * sr2
            sr12 = sr6 * sr6

            u = 4.0_dp * self%epsilon0 * (sr12 - sr6) - self%U_rCut
        end if
    end subroutine

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

    pure subroutine forceLennardJones(self, r1, r2, f, r_mag)
        class(simState), intent(in)         :: self
        real(dp), intent(in),  dimension(3) :: r1, r2
        real(dp), intent(out), dimension(3) :: f
        real(dp), intent(out)               :: r_mag
        real(dp)                            :: r2_val, inv_r2, sigma2
        real(dp)                            :: sr2, sr6, sr12, l_box
        real(dp), dimension(3)              :: dr

        f = 0.0_dp
        l_box = self%volume ** (1.0_dp / 3.0_dp)
        dr = r2 - r1
        dr = dr - l_box * anint(dr / l_box)  ! anint() is nint() but returns real(dp)
        r2_val = sum(dr**2)
        r_mag = sqrt(r2_val)

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
        class(simState), intent(inout)  :: self
        real(dp), dimension(3)          :: f
        integer,  dimension(self%nbins) :: count_aux
        integer                         :: i, j, bin
        real(dp)                        :: dr, l_box

        self%force = 0.0_dp
        f = 0.0_dp
        count_aux = 0
        l_box = self%volume ** (1.0_dp/3.0_dp)

        ! The following comments is for parallelization
        ! Newton's third law is sacrified in order to parallelize
        ! default(none) is for explicitely declaring every variable intend
        ! shared(self) is for the shared objects across threads
        ! private(i, j, f) is for the objects that threads handle individually
        ! schedule(static) is for equally workload distribution across threads

        !$omp parallel do default(none)              &
        !$omp             shared(self, l_box)        &
        !$omp             private(i, j, f, dr, bin)  &
        !$omp             reduction(+:count_aux)     &
        !$omp             schedule(static)
        do i = 1, self%numParticles
            do j = 1, self%numParticles

                if (j == i) then 
                    cycle
                end if

                call self%forceField(self%position(:, i), &
                                     self%position(:, j), &
                                     f, dr)
                self%force(:, i) = self%force(:, i) + f

                ! For correlation function
                if (dr < l_box/2) then  ! sphere for easy computing ideal gas
                    bin = int(dr * self%nbins * 2.0_dp / l_box) + 1
                    count_aux(bin) = count_aux(bin)  + 1
                end if
            end do
        end do
        !$omp end parallel do

        self%rCount = self%rCount + count_aux


    end subroutine getForces


    ! =====================
    ! Extracting Statistics
    ! =====================
    real(dp) function getTotalPotentialEnergy(self)
        class(simState), intent(in)                :: self
        integer                                    :: i, j
        real(dp)                                   :: u, energy

        u = 0.0_dp
        energy = 0.0_dp


        ! reduction(+:X) do this: 
        !       1. Creates an independent copy of X for each thread.
        !       2. Initializes each copy to 0.
        !       3. (loop operations).
        !       3. When the loop end adds up each copy of X with the original
        !          X variable of the main node.
        !
        ! schedule(dynamic) evaluates if a thread finished its task and then 
        ! assigns a new task which haven't been done yet by other threads.
        ! It is not as efficient in task assignment as schedule(static) 
        ! because it needs to evaluate the workload of threads periodically.
        ! But, in this case, with triangular calculation, works perfectly.
        !$omp parallel do default(none)       &
        !$omp             shared(self)        &
        !$omp             private(i, j, u)    &
        !$omp             reduction(+:energy) &
        !$omp             schedule(dynamic)
        do i = 1, self%numParticles
            do j = i+1, self%numParticles

                call self%potentialField(self%position(:, i), &
                                    self%position(:, j), u)

                energy = energy + u
            end do
        end do
        !$omp end parallel do

        getTotalPotentialEnergy = energy
    end function getTotalPotentialEnergy

    real(dp) function getTotalKineticEnergy(self)
        class(simState), intent(in) :: self
        integer                     :: i
        real(dp)                    :: energy

        energy = 0.0_dp

        !$omp parallel do default(none)        &
        !$omp             shared(self)         &
        !$omp             private(i)           &
        !$omp             reduction(+:energy)  &
        !$omp             schedule(static)
        do i = 1, self%numParticles
            energy = energy + sum(self%velocity(:, i)**2)
        end do
        !$omp end parallel do

        getTotalKineticEnergy = 0.5_dp * self%mass * energy
    end function getTotalKineticEnergy

    real(dp) function getPreassure(self)
        class(simState), intent(in) :: self
        real(dp), dimension(3)      :: f, r
        integer                     :: i, j
        real(dp)                    :: aux, l_box, dummy


        f = 0.0_dp
        r = 0.0_dp
        aux = 0.0_dp
        l_box = self%volume ** (1.0_dp / 3.0_dp)

       
        !$omp parallel do default(none)              &
        !$omp             shared(self, l_box)        &
        !$omp             private(i, j, f, r, dummy) &
        !$omp             reduction(+:aux)           &
        !$omp             schedule(dynamic)
        do i = 1, self%numParticles
            do j = i+1, self%numParticles
                r = self%position(:, i) - self%position(:, j)
                r = r - l_box * anint(r / l_box) 

                call self%forceField(self%position(:, i), &
                                     self%position(:, j), &
                                     f, dummy)

                aux = aux + dot_product(f, r)
            end do
        end do
        !$omp end parallel do

        getPreassure = self%numParticles*self%kb*self%getTemperature() / &
                       self%volume + aux / (3.0_dp * self%volume)
    end function getPreassure




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

        ! Update press and temp here?
        !       Ans-> Better not, too expensive to compute at each step.

    end subroutine verletStep


    ! =====================
    ! Storing data
    ! =====================
    subroutine storeData(self, i)
        class(simState),  intent(inout) :: self  ! only for modifying kEnergy
        integer,          intent(in)    :: i
        integer, dimension(5), save     :: io_unit
        integer                         :: status_code, j
        character(len=20), dimension(5) :: files
        real(dp), dimension(5)          :: stateFunctions
        real(dp)                        :: dr

        stateFunctions = 0.0_dp
        dr = self%volume**(1.0_dp/3.0_dp) / (2.0_dp * self%nbins)

        files(1) = "position.xyz"
        files(2) = "velocity.xyz"
        files(3) = "stateFunctions.txt"
        files(4) = "parameters.txt"
        files(5) = "corrFunction.txt"


        if (i == 1) then
            do j = 1, 5
                open(newunit  = io_unit(j),                            &
                     file     = self%dataDirectory // files(j),        &
                     status   = "replace",                             &
                     action   = "write",                               &
                     iostat   = status_code)
                
                if (status_code /= 0) then
                    print *, "ERROR: Could not open file: ", & 
                              self%dataDirectory // files(j)
                    print *, "Does the target directory exist?"
                    stop 1
                end if
            end do


            write(io_unit(3), '(6(A, 4X))') "step", "time", "T", "U", "K", "P"

            ! Saving Simulation Parameters
            write(io_unit(4), '(A, I8)') "numParticles", self%numParticles
            write(io_unit(4), '(A, I8)') "numSteps",     self%numSteps
            write(io_unit(4), '(A, I8)') "framePeriod",  self%framePeriod
            write(io_unit(4), '(A, F12.6)') "volume",    self%volume
            write(io_unit(4), '(A, F12.6)') "mass",      self%mass
            write(io_unit(4), '(A, F12.6)') "dt",        self%dt
            write(io_unit(4), '(A, F12.6)') "kb",        self%kb
            write(io_unit(4), '(A, F12.6)') "temp",      self%temp
            write(io_unit(4), '(A, F12.6)') "sigma",     self%sigma
            write(io_unit(4), '(A, F12.6)') "epsilon",   self%epsilon0
            write(io_unit(4), '(A, F12.6)') "rCut",      self%rCut
            close(io_unit(4))
        end if

        if (mod(i, self%framePeriod) == 0) then

            ! Calculate and store the total kinetic energy in self%kEnergy
            call self%setKEnergy()
            stateFunctions(1) = real(i, dp) * self%dt
            stateFunctions(2) = self%getTemperature()
            stateFunctions(3) = self%getTotalPotentialEnergy()
            stateFunctions(4) = self%kEnergy
            stateFunctions(5) = self%getPreassure()

            ! Positions
            write(io_unit(1), '(I8)') self%numParticles
            write(io_unit(1), '(A, I8, 5(A, F12.6))') &
                  "Step = ",           i,                  &
                  ", Time = ",         stateFunctions(1),  &
                  ", Temperature = ",  stateFunctions(2),  &
                  ", Potential = ",    stateFunctions(3),  &
                  ", Kinetic = ",      stateFunctions(4),  &
                  ", Preassure = ",    stateFunctions(5)              
            do j = 1, self%numParticles
                write(io_unit(1), '(A, 3(1X, F12.6))') "H", &
                                                       self%position(1, j), &
                                                       self%position(2, j), &
                                                       self%position(3, j)
            end do

            ! Velocities
            write(io_unit(2), '(I8)') self%numParticles
            write(io_unit(2), '(A, I8, 5(A, F12.6))') &
                  "Step = ",           i,                  &
                  ", Time = ",         stateFunctions(1),  &
                  ", Temperature = ",  stateFunctions(2),  &
                  ", Potential = ",    stateFunctions(3),  &
                  ", Kinetic = ",      stateFunctions(4),  &
                  ", Preassure = ",    stateFunctions(5)              
            do j = 1, self%numParticles
                write(io_unit(2), '(A, 3(1X, F12.6))') "H", &
                                                       self%velocity(1, j), &
                                                       self%velocity(2, j), &
                                                       self%velocity(3, j)
            end do

            ! State functions
            write(io_unit(3), '(I8, 5(F16.6))') i, stateFunctions(1),  &
                                                   stateFunctions(2),  & 
                                                   stateFunctions(3),  & 
                                                   stateFunctions(4),  & 
                                                   stateFunctions(5)
        end if

        if (i == self%numSteps) then

            do j = 1, self%nbins
                write(io_unit(5), '(F12.6, 4X, I16)') dr * (j - 0.5_dp), &
                                                      self%rCount(j)
            end do



            close(io_unit(1))
            close(io_unit(2))
            close(io_unit(3))
            close(io_unit(5))
        end if
    end subroutine storeData
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
    integer                     :: numParticles, numSteps, framePeriod, nbins
    real(dp)                    :: volume, mass, dt, kb, temp, sigma, epsilon0
    real(dp)                    :: rCut
    character(len=*), parameter :: dataDirectory = "../postpro/"
    ! character(len=*), parameter :: dataDirectory = "../postpro/idealgas/"


    ! =====================
    ! Setting up parameters
    ! =====================
    numParticles  = 512             ! Total number of particles
    numSteps      = 10000           ! Total number of time steps
    framePeriod   = 10              ! Steps between Coordinates storage
    nbins         = 1000            ! Set number of bins for corr function
    volume        = 10000.0_dp      ! Volume of the 1D system
    mass          = 1.0_dp          ! Particles' mass
    dt            = 1e-3_dp         ! Size of time step
    kb            = 1.0_dp          ! Boltzmann's constant
    temp          = 0.001_dp        ! Temperature of the system
    sigma         = 1.0_dp          ! For Lennard-Jones
    epsilon0      = 1.0_dp          ! For Lennard-Jones
    ! epsilon0      = 0.0_dp          ! For Lennard-Jones
    rCut          = 2.5_dp * sigma  ! For truncating far interactions


    call state%setOutputDir(dataDirectory)
    call state%setNumParticles(numParticles)
    call state%setFramePeriod(framePeriod)
    call state%setNbins(nbins)
    call state%setNumSteps(numSteps)
    call state%setVolume(volume)
    call state%setMass(mass)
    call state%setDt(dt)
    call state%setKb(kb)
    call state%setTemperature(temp)
    call state%setSigma(sigma)
    call state%setEpsilon(epsilon0)
    call state%setRCut(rCut)

    call state%setKEnergy0()


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
        ! call state%storeVelocity(i)
        ! call state%storePosition(i)
        call state%storeData(i)
        call state%verletStep()
    end do


    call state%deallocateXVF()
end program main
