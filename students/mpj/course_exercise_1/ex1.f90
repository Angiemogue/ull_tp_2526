program leapfrog
  use geometry
  use particle
  implicit none

  integer :: n
  real(dp) :: dt, t_end, t, dt_out, t_out

  ! Softing parameter to avoid division by zero
  real(dp), parameter :: eps = 1.0e-5_dp

  ! Dynamic arrays to store particles and their accelerations
  type(particle3d), allocatable :: particles(:)
  type(vector3d), allocatable :: a(:)

  ! Output variables
  character(len=200) :: output_file
  integer :: out_unit

  output_file = "output.dat"

  call get_input(dt, dt_out, t_end, n, particles)
  call open_output_file(output_file, out_unit)

  allocate(a(n))
  
  call compute_accelerations(particles, a, n, eps)


  ! Leapfrog integration
  t_out = 0.0_dp
  do t = 0.0_dp, t_end, dt
     call half_step_vel_pos(particles, a, n, dt)
     call compute_accelerations(particles, a, n, eps)
     call half_step_vel(particles, a, n, dt)
     
     ! Output
     t_out = t_out + dt
     if (t_out >= dt_out) then
      call write_output(out_unit, t, particles, n)
      t_out = 0.0_dp
     end if

  end do

  close(out_unit)

  ! Clean up
  deallocate(particles)
  deallocate(a)


contains

  subroutine get_input(dt, dt_out, t_end, n, particles)
   implicit none
   
   real(dp), intent(out) :: dt, dt_out, t_end
   integer, intent(out)  :: n
   type(particle3d), allocatable, intent(out) :: particles(:)

   ! Input parameters
   character(len=200) :: input_file
   integer :: i, ios, unit
   
   ! Opening the input file
   unit = 10
   input_file = "input.dat"
   open(unit, file=trim(input_file), status='old', action='read', iostat=ios)
   
   if (ios /= 0) then
      print *, "Error: can't open file ", trim(input_file)
      stop
   end if

   ! Reading the input from the input file
   read (unit, *) dt        ! integration step
   read (unit, *) dt_out    ! time between two positions prints
   read (unit, *) t_end     ! total time of the simulation
   read (unit, *) n         ! number of particles

   allocate(particles(n))

   ! Reading masses, positions and velocities
   do i = 1, n
      read (unit, *) particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, &
                     particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
   end do

   close(unit)

  end subroutine get_input


  subroutine open_output_file(output_file, out_unit)
   implicit none
   character(len=200), intent(in) :: output_file
   integer, intent(out) :: out_unit
   integer :: ios
   
   ! Opening the output file
   out_unit = 20
   open(out_unit, file=trim(output_file), status='replace', action='write', iostat=ios)

   if (ios /= 0) then
     print *, "Error: can't open file ", trim(output_file)
     stop
   end if

  end subroutine open_output_file


  subroutine compute_accelerations(particles, a, n, eps)
   integer, intent(in) :: n
   real(dp), intent(in) :: eps
   type(particle3d), intent(in) :: particles(n)
   type(vector3d), intent(inout) :: a(n)
   integer :: i, j
   real(dp) :: r2, r3
   type(vector3d) :: rji

   do i = 1, n
      do j = i+1, n
         rji = particles(j)%p - particles(i)%p
         r2 = rji%x**2 + rji%y**2 + rji%z**2 + eps**2
         r3 = r2 * sqrt(r2)
         a(i) = a(i) + particles(j)%m * rji / r3    ! Newton force
         a(j) = a(j) - particles(i)%m * rji / r3
      end do
   end do

  end subroutine compute_accelerations


  subroutine half_step_vel_pos(particles, a, n, dt)
   implicit none
   integer, intent(in) :: n
   real(dp), intent(in) :: dt
   type(particle3d), intent(inout) :: particles(:)
   type(vector3d), intent(inout) :: a(:)
   integer :: i
   
   do i = 1, n
      particles(i)%v = particles(i)%v + a(i) * (dt / 2.0_dp)
      particles(i)%p = particles(i)%p + particles(i)%v * dt
      a(i) = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
   end do

  end subroutine half_step_vel_pos


  subroutine half_step_vel(particles, a, n, dt)
   implicit none
   integer, intent(in) :: n
   real(dp), intent(in) :: dt
   type(particle3d), intent(inout) :: particles(:)
   type(vector3d), intent(inout) :: a(:)
   integer :: i
   
   do i = 1, n
       particles(i)%v = particles(i)%v + a(i) * (dt / 2.0_dp)
   end do

  end subroutine half_step_vel


  subroutine write_output(out_unit, t, particles, n)
   implicit none
   integer, intent(in) :: out_unit      ! Unidad de archivo
   real(dp), intent(in) :: t            ! Tiempo actual
   type(particle3d), intent(in) :: particles(:)
   integer, intent(in) :: n
   integer :: i
   
   write(out_unit, '(F10.5,1X)', advance='no') t  ! We update the time
      ! We write the particles positions:
      do i = 1,n
         write(out_unit, '(3F16.10,1X)', advance='no') particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
      end do
   write(out_unit,*)  ! line break

  end subroutine write_output
  
  
end program leapfrog
