program ex2

    use geometry
    use particles
    use tree

    implicit none

    integer :: ierr ! This is to check and deallocate
    character(len=100) :: infile   ! Name of the input file
    character(len=200) :: in ! argument in the terminal for the input file


    call get_command_argument(1,in)

    if (len_trim(in)==0) then ! If input is not given 
        print *, "Input example: ./ex2 <inputfile>"
        stop
    endif

    infile = trim(in) ! remove empty space
    print *, "Reading input file: ", trim(infile)


    ! We use the subroutine to read the file
    call read_input(infile,in, n,dt,dt_out,t_end,part,a)
    
    !! Head node initialization
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(head)
    call calculate_ranges(head)
    head%type = 0
    call nullify_pointers(head)
    
    !! Initial tree creation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1,n
        call find_cell(head,temp_cell,part(i)%p)
        call place_cell(temp_cell,part(i)%p,i)
    end do

    call borrar_empty_leaves(head)
    call calculate_masses(head)
    
    !! Compute initial accelerations
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1, n
      a(i) = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
    end do
    
    call calculate_forces(head)

    call main_loop(part, a, n, dt, dt_out, t_end)

    print *, "Simulation complete :) Particles positions save in output.dat"

    ! Finally we check the status of the arrays and deallocate
    if (allocated(part)) deallocate(part, stat=ierr)
    if (allocated(a)) deallocate(a, stat=ierr)
    
    
    contains

    subroutine read_input(infile, in,  n, dt, dt_out, t_end, part, a)
        character(len=100) :: infile   
        character(len=200) :: in 
        integer, intent(out) :: n
        real(dp), intent(out) :: dt, dt_out, t_end
        type(particle3d), allocatable, intent(out) :: part(:)
        type(vector3d), allocatable, intent(out) :: a(:)
        integer :: i, ierr 

    
        open(unit=10, file=infile, status='old', action='read', iostat=ierr)
        if (ierr /= 0) then 
            print *, "Error: The file doesn't exist or it cannot be open: ", trim(infile)
        stop
        endif

  
        read(10,*) dt
        read(10,*) dt_out
        read(10,*) t_end
        read(10,*) n

    ! We are going to associate the components of the particle type for simplicity, and now we read the last lines, n times to allocate n particles

        allocate(part(n))
        allocate(a(n)) 

        do i = 1,n
        read(10,*) part(i)%m, & 
                part(i)%p%x, part(i)%p%y, part(i)%p%z, &
                part(i)%v%x, part(i)%v%y, part(i)%v%z
        end do

        close(10)

        do i = 1, n
            a(i) = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
        end do

   end subroutine read_input



end program ex2