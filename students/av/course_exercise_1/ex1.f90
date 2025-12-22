program ex1
   ! Import our custom modules # Lecture 4
   ! particle also gives us access to geometry
   use geometry
   use particle
   implicit none



   !--- Variable declarations ---!

   integer :: i, n_particles
   real(kind=dp) :: dt, t_end, t, dt_out, t_out

   ! Temporal storage variables
   real(kind=dp) :: mass_in, px, py, pz, vx, vy, vz

   ! Main particle array. Allocatable # Lecture 5
   type(particle3d), dimension(:), allocatable :: particles

   ! Aceleration array. Also allocatable
   type(vector3d), dimension(:), allocatable :: accelerations

   ! I/O variables for command line flags # Lecture 7
   character(len=256) :: input_filename, output_filename
   integer :: num_args, input_unit, output_unit, iostat

   ! Temporal storage for arguments
   character(len=256) :: arg_string

   ! Checking for input file name
   logical :: input_file_set = .false.

   ! Selecting the chanel for the files
   parameter (input_unit = 10, output_unit = 20)
   output_filename = "output.dat" ! Default output name



   !--- Command line arguments (-i, -o) ---!

   ! Obtaining the arguments
   num_args = command_argument_count()
   i = 1

   ! Looping through all arguments
   do while (i <= num_args)

      call get_command_argument(i, arg_string)

      ! Checking the corresponding flags
      ! by selecting a case
      select case (trim(arg_string))
      
         ! Input file case
         case ("-i")
            if (i + 1 > num_args) then
               print *, "Error: -i needs an input file name."
               stop 1
            end if

            ! Continue to the following argument
            ! to get input file name
            i = i + 1 
            call get_command_argument(i, input_filename)
            input_file_set = .true.
         
         ! Output file case
         case ("-o")
            if (i + 1 > num_args) then
               print *, "Error: -o needs an output file name."
               stop 1
            end if

            ! Continue to the following argument
            ! to get output file name
            i = i + 1 
            call get_command_argument(i, output_filename)
         
         ! Unknown argument
         case default
            print *, "Error: Unknown argument: ", trim(arg_string)
            print *, "Use: ./ex1 -i <inputfile> [-o <outputfile>]"
            stop 1
      
      end select
      
      ! Continue to the following argument
      i = i + 1

   end do

   ! Cheking that input file was parsed
   if (.not. input_file_set) then
      print *, "Error: Unexpecified input file name with -i."
      print *, "Use: ./ex1 -i <inputfile> [-o <outputfile>]"
      stop 2
   end if



   !--- Reading initial simulation set up file (input file) ---!

   ! Opening the input file for reading it # Lecture 7
   open(unit=input_unit, &
         file=trim(input_filename), &
         status='old', &
         action='read', &
         iostat=iostat)

   if (iostat /= 0) then
      print *, "Error: Could not open input file: ", trim(input_filename)
      stop 3
   end if

   ! Read simulation parameters line by line
   ! Using inline if statement
   read(input_unit, *, iostat=iostat) dt
   if (iostat /= 0) call handle_read_error("dt")
   read(input_unit, *, iostat=iostat) dt_out
   if (iostat /= 0) call handle_read_error("dt_out")
   read(input_unit, *, iostat=iostat) t_end
   if (iostat /= 0) call handle_read_error("t_end")
   read(input_unit, *, iostat=iostat) n_particles
   if (iostat /= 0) call handle_read_error("n_particles")

   ! Input values validation
   if (n_particles <= 0) then
      print *, "Error: Number of particles must be positive."
      stop 4
   end if
   if (dt <= 0.0_dp .or. dt_out <= 0.0_dp) then
      print *, "Error: Time steps (dt, dt_out) must be positive."
      stop 5
   end if

   ! Allocate memory for arrays
   allocate(particles(n_particles), stat=iostat)
   if (iostat /= 0) stop "Error allocating particles array."
   allocate(accelerations(n_particles), stat=iostat)
   if (iostat /= 0) stop "Error allocating accelerations array."


   ! Loop to read the initial state for each particle 
   ! based on input.dat format.
   do i = 1, n_particles

      ! Data for each row: mass, pos x, pos y, pos z, vel x, vel y, vel z
      read(input_unit, *, iostat=iostat) mass_in, px, py, pz, vx, vy, vz

      ! Checking which particle causes the error
      ! // for concatenate strings. achar(iachar('0')+i for convert i into a string # Gemini
      if (iostat /= 0) call handle_read_error("particle data for particle "//achar(iachar('0')+i))
      
      ! Use derived type constructors
      particles(i)%m = mass_in
      particles(i)%p = point3d(px, py, pz)
      particles(i)%v = vector3d(vx, vy, vz)
      
      ! Validate particle mass
      if (particles(i)%m <= 0.0_dp) then
         print *, "Error: Particle mass must be positive. Check particle ", i
         stop 6
      end if
   end do

   close(input_unit)



   !--- Preparing putput file ---!

   ! Open the output file for writing
   open(unit=output_unit, &
         file=output_filename, &
         status='replace', &
         action='write', &
         iostat=iostat)

   ! Cheking if errors
   if (iostat /= 0) then
      print *, "Error: Could not open output file: ", trim(output_filename)
      stop 7
   end if



   !--- Adapting Leapfrog (Lf) algorithm ---!
   
   ! Initialize accelerations
   accelerations = vector3d(0.0_dp, 0.0_dp, 0.0_dp) 
   call calculate_accelerations(particles, n_particles, accelerations)

   ! Initializing the time
   t_out = 0.0_dp
   t = 0.0_dp
   
   ! Initial output at t=0 is the initial state
   call write_output(output_unit, t, particles, n_particles)
   t_out = t_out + dt_out


   ! Infinite do loop with an exit condition
   do
      !--- Lf Step 1: Update velocities by a half time step ---!
      do i = 1, n_particles
         particles(i)%v = particles(i)%v + accelerations(i) * (dt / 2.0_dp)
      end do
      
      !--- Lf Step 2: Update positions by a full time step ---!
      do i = 1, n_particles
         particles(i)%p = particles(i)%p + particles(i)%v * dt
      end do
      
      !--- Lf Step 3: Recompute accelerations at the new positions ---!
      call calculate_accelerations(particles, n_particles, accelerations)
      
      !--- Lf Step 4: Update velocities by the second half time step ---!
      do i = 1, n_particles
         particles(i)%v = particles(i)%v + accelerations(i) * (dt / 2.0_dp)
      end do
      
      ! Increment simulation time
      t = t + dt

      ! Exit loop condition
      if (t > t_end) exit

      ! Lf Step 5: Writing particles' state in the correct timestep
      if (t >= t_out - dt*0.01_dp) then 
         call write_output(output_unit, t, particles, n_particles)
         t_out = t_out + dt_out
      end if
      
   end do

   close(output_unit)

   ! Deallocate dynamic arrays to free memory
   deallocate(particles)
   deallocate(accelerations)

   contains

      ! Internal subroutine for error handling during file read.
      subroutine handle_read_error(variable_name)
         character(len=*), intent(in) :: variable_name
         print *, "Error reading ", trim(variable_name), " from input file."
         print *, "Check file format and content."
         stop 8
      end subroutine handle_read_error

      ! Internal subroutine to calculate accelerations
      subroutine calculate_accelerations(p_array, num_p, acc_array)

         ! Defining internal variables
         type(particle3d), dimension(:), intent(in) :: p_array
         integer, intent(in) :: num_p
         type(vector3d), dimension(:), intent(out) :: acc_array
         
         integer :: idx1, idx2
         real(kind=dp) :: d, dist_cubed
         type(vector3d) :: pos_diff

         ! Initiialize acceleration for each timestep
         acc_array = vector3d(0.0_dp, 0.0_dp, 0.0_dp) 
         do idx1 = 1, num_p
            do idx2 = idx1 + 1, num_p

               ! Calculate vector from particle idx1 to idx2 using subpp function
               pos_diff = p_array(idx2)%p - p_array(idx1)%p
               
               ! Computing distance between particle idx1 to idx2
               d = distance(p_array(idx1)%p, p_array(idx2)%p)
               
               ! Avoid division by zero
               ! for very close particles
               if (d > epsilon(0.0_dp)) then 
                  dist_cubed = d**3
                  acc_array(idx1) = acc_array(idx1) + (p_array(idx2)%m / dist_cubed) * pos_diff
                  acc_array(idx2) = acc_array(idx2) - (p_array(idx1)%m / dist_cubed) * pos_diff
               end if
            end do
         end do
      end subroutine calculate_accelerations

      ! Internal subroutine to write output data
      subroutine write_output(unit_num, current_time, p_array, num_p)

         ! Defining internal variables
         integer, intent(in) :: unit_num, num_p
         real(kind=dp), intent(in) :: current_time
         type(particle3d), dimension(:), intent(in) :: p_array
         integer :: idx
         
         ! Write time and all particle positions (x,y,z) using inline loop
         write(unit_num, *) current_time, (p_array(idx)%p%x, p_array(idx)%p%y, &
                                          p_array(idx)%p%z, idx = 1, num_p)
      end subroutine write_output

end program ex1