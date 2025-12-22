program generate_input_disk
  !============================================================
  ! Generate initial conditions for a thin disk around a central mass
  !============================================================

  use iso_fortran_env, only: real64
  implicit none

  ! Precision
  integer, parameter :: dp = real64

  ! Simulation parameters
  integer, parameter :: n_particles = 50
  real(dp), parameter :: disk_radius = 1.0_dp
  character(len=*), parameter :: filename = "input_disk.dat"

  ! Input file parameters
  real(dp), parameter :: dt = 0.001_dp
  real(dp), parameter :: dt_out = 0.01_dp
  real(dp), parameter :: t_end = 10.0_dp

  real(dp), parameter :: particle_mass = 1.0e-3_dp
  real(dp), parameter :: M_central = 1.0_dp
  real(dp), parameter :: G = 1.0_dp

  ! Vertical dispersion
  real(dp), parameter :: z_dispersion = 0.05_dp
  real(dp), parameter :: vz_dispersion = 0.05_dp

  ! Arrays for positions and velocities
  real(dp), allocatable :: pos_x(:), pos_y(:), pos_z(:)
  real(dp), allocatable :: vel_x(:), vel_y(:), vel_z(:)
  real(dp), allocatable :: radius(:), angle(:)
  integer :: i
  real(dp) :: velocity_mag

  ! Random seed
  call random_seed()

  ! Allocate arrays
  allocate(pos_x(n_particles), pos_y(n_particles), pos_z(n_particles))
  allocate(vel_x(n_particles), vel_y(n_particles), vel_z(n_particles))
  allocate(radius(n_particles), angle(n_particles))

  ! Generate positions
  call random_number(radius)
  radius = sqrt(radius) * disk_radius

  call random_number(angle)
  angle = angle * 2.0_dp * acos(-1.0_dp)  ! 2*pi

  call random_number(pos_z)
  pos_z = (pos_z - 0.5_dp) * 2.0_dp * z_dispersion  ! Normalized approx ±dispersion

  pos_x = radius * cos(angle)
  pos_y = radius * sin(angle)

  ! Generate velocities
  velocity_mag = 0.0_dp
  do i = 1, n_particles
     velocity_mag = sqrt(G * M_central / max(radius(i), 0.05_dp))
     vel_x(i) = -velocity_mag * sin(angle(i))
     vel_y(i) =  velocity_mag * cos(angle(i))
     call random_number(vel_z(i))
     vel_z(i) = (vel_z(i) - 0.5_dp) * 2.0_dp * vz_dispersion
  end do

  ! Write file
  open(unit=10, file=filename, status='replace', action='write')

  ! Global parameters
  write(10, '(3(F12.6,1X),I6)') dt, dt_out, t_end, n_particles+1

  ! Central mass at origin
  write(10, '(F12.6,6F12.6)') M_central, 0.0_dp,0.0_dp,0.0_dp, 0.0_dp,0.0_dp,0.0_dp

  ! Orbiting particles
  do i = 1, n_particles
     write(10, '(F12.6,6F12.6)') particle_mass, pos_x(i), pos_y(i), pos_z(i), &
                                  vel_x(i), vel_y(i), vel_z(i)
  end do

  close(10)

  print *, "File '", trim(filename), "' successfully generated with ", &
           n_particles, " particles + central mass."

  ! Deallocate
  deallocate(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, radius, angle)

end program generate_input_disk
