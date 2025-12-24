program ex2_mpi

    use geometry
    use particles
    use tree_mpi

    implicit none

    integer :: ierr ! This is to check and deallocate
    character(len=100) :: infile   ! Name of the input file
    character(len=200) :: in ! argument in the terminal for the input file
    real(dp) :: t1, t2

    !! MPI initialization
    !!!!!!!!!!!!!!!!!!!!!!!!
    call MPI_INIT ( error )
    call MPI_COMM_SIZE ( MPI_COMM_WORLD, p, error )
    call MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, error )

    call get_command_argument(1,in)

    if (len_trim(in)==0) then ! If input is not given 
        print *, "Input example: ./ex2 <inputfile>"
        stop
    endif

    infile = trim(in) ! remove empty space

    if (my_rank == 0) then 
        print *, "Reading input file: ", trim(infile)
    end if

    ! We use the subroutine to read the file, and initialize MPI
    call read_input(infile,in, n,dt,dt_out,t_end,part,a)
    

    !! Initialize head node
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(head)
    call Calculate_ranges(head)
    head%type = 0
    call Nullify_Pointers(head)
    
 
    !! Initial tree creation
    !!
    !! Note: we only consider
    !! for the tree creation the
    !! particles that belong to the range
    !! of the head node (head), given
    !! by the Calculate_ranges subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1, n
        if (Belongs(part(i)%p, head)) then
            call Find_Cell(head, temp_cell, part(i)%p)
            call Place_Cell(temp_cell, part(i)%p, i)
        end if
    end do
    
    call Borrar_empty_leaves(head)
    call Calculate_masses(head)

     !! Compute initial accelerations
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Here we start the time calculation

    t1 = MPI_WTIME()

    do i = 1, n
      a(i) = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
    end do

    call Calculate_forces(head)
    call MPI_ALLREDUCE(a, total_a, n*3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, error)
    

    

    call main_loop(part, a, total_a, n, dt, dt_out, t_end)

    if (my_rank == 0) then
        ! In order to print only once in the terminal
        print *, "Simulation complete :) Particles positions save in output.dat"
    end if

        
    t2 = MPI_WTIME()

    if (my_rank == 0) then
        print *, "Total execution time (s): ", t2 - t1
    end if

    
    ! Finally we check the status of the arrays and deallocate
    if (allocated(part)) deallocate(part, stat=ierr)
    if (allocated(a)) deallocate(a, stat=ierr)
    if (allocated(total_a)) deallocate(total_a, stat=ierr)

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
        if (my_rank == 0) then
            if (ierr /= 0) then 
                print *, "Error: The file doesn't exist or it cannot be open: ", trim(infile)
            stop
            endif
        endif

    
        !! Matrix initialization
        !! The master reads from file and broadcasts
        !! all variables to all slaves
        !!
        if ( my_rank == 0 ) then

            read(10,*) dt
            read(10,*) dt_out
            read(10,*) t_end
            read(10,*) n
        
            call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(dt_out, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(t_end, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
        
            allocate(part(n))
            allocate(a(n)) 
            allocate(total_a(n))

            do i = 1,n
                read(10,*) part(i)%m, & 
                part(i)%p%x, part(i)%p%y, part(i)%p%z, &
                part(i)%v%x, part(i)%v%y, part(i)%v%z
            end do

            close(10)


            call MPI_BCAST(part%m, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(part%p%x, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(part%p%y, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(part%p%z, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(part%v%x, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(part%v%y, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(part%v%z, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
        else
            call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(dt_out, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(t_end, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)

            allocate(part(n))
            allocate(a(n))
            allocate(total_a(n))
        
            call MPI_BCAST(part%m, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(part%p%x, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(part%p%y, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(part%p%z, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(part%v%x, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(part%v%y, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
            call MPI_BCAST(part%v%z, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
        end if

   end subroutine read_input



end program ex2_mpi