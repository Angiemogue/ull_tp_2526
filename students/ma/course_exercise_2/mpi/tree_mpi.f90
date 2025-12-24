
! We choose the second option of the notes (parallelizing all the code)

module tree_mpi

    use geometry
    use particles

    implicit none
    !! We use "include" instead of "USE", since with "USE" we couldn't use MPI_ALLREDUCE
    !! USE MPI


    include 'mpif.h'


    !! MPI Variables (rank, number of processors, etc.)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: my_rank, p, error
    integer, dimension(MPI_STATUS_SIZE) :: status
    
    !! N-body problem variables
    !!
    !! A new variable with respect to Solution 1 is
    !! total_a. In the array "a" we will calculate how the
    !! particles in our tree affect all particles.
    !! With an MPI_ALLREDUCE we will sum all local "a"s in
    !! "total_a", thus calculating the total acceleration for
    !! all particles (affected by the trees of the 8 octants)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer :: i,j,k,n
    real(dp) :: dt, t_end, t, dt_out, t_out
    real(dp), parameter :: theta = 1
    type(particle3d), dimension(:), allocatable :: part
    type(vector3d), dimension(:), allocatable :: a, total_a

    

    type range ! limits of a cell
        real(dp), dimension(3) :: min,max
    end type range
    
    type cptr
        type(cell), pointer :: ptr
    end type cptr
    
    type cell
        type (range) :: range
        type(particle3d) :: part
        integer :: pos
        integer :: type !! 0 = no particule; 1 = particle; 2 = group of particles
        type(point3d) :: c_o_m !centre of mass
        type (cptr), dimension(2,2,2) :: subcell
    end type cell
    
    type (cell), pointer :: head, temp_cell
    
   

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Calculate_Ranges
    !!
    !!
    !!
    !! Similar to the subroutine in the first !!
    !! solution, but we modify the ranges     !!
    !! of the head node so that each of the   !!
    !! 8 processors gets its corresponding    !!
    !! part.                                   !!
    !!                                         !!
    !! Uses: Convert_rank_octant              !!
    !!
    !!
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Calculate_Ranges(goal)
        type(CELL), pointer :: goal
        real(dp), dimension(3) :: mins, maxs, medios, new_range_min, new_range_max
        real(dp) :: span
        integer, dimension(3) :: octant
        real(dp), allocatable :: xs(:), ys(:), zs(:)

        xs = [ (part(i)%p%x, i=1,n) ]
        ys = [ (part(i)%p%y, i=1,n) ]
        zs = [ (part(i)%p%z, i=1,n) ]
        
        mins = [ minval(xs), minval(ys), minval(zs) ]
        maxs = [ maxval(xs), maxval(ys), maxval(zs) ]

        span = maxval(maxs - mins) * 1.1 !! Add 10% so particles don't fall at the edge
        medios = (maxs + mins) / 2.0
        
        goal%range%min = medios - span/2.0
        goal%range%max = medios + span/2.0
        
        octant = Convert_rank_octant(my_rank)
        new_range_min = Calcular_Range(0, goal, octant)
        new_range_max = Calcular_Range(1, goal, octant)
        
        goal%range%min = new_range_min
        goal%range%max = new_range_max

    end subroutine Calculate_Ranges
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Convert_rank_octant
    !!
    !!
    !! Used by Calculate_ranges
    !!
    !!
    !! Simply makes a mapping of rank !!
    !! to octants.                    !!
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function Convert_rank_octant(rank)
        integer :: rank
        integer, dimension(3) :: Convert_rank_octant
        
        select case (rank)
        case (0)
            Convert_rank_octant = (/1,1,1/)
        case (1)
            Convert_rank_octant = (/1,1,2/)
        case (2)
            Convert_rank_octant = (/1,2,1/)
        case (3)
            Convert_rank_octant = (/1,2,2/)
        case (4)
            Convert_rank_octant = (/2,1,1/)
        case (5)
            Convert_rank_octant = (/2,1,2/)
        case (6)
            Convert_rank_octant = (/2,2,1/)
        case (7)
            Convert_rank_octant = (/2,2,2/)
        end select
    end function Convert_rank_octant
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Find_Cell
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine Find_Cell(root, goal, partp)
        type(point3d) :: partp
        type(CELL),pointer :: root,goal,temp
        integer :: i,j,k
        select case (root%type)
        case (2)
            out: do i = 1,2
                do j = 1,2
                    do k = 1,2
                        if (Belongs(partp,root%subcell(i,j,k)%ptr)) then
                            call find_cell(root%subcell(i,j,k)%ptr,temp,partp)
                            goal => temp
                            exit out
                        end if
                    end do
                end do
            end do out
        case default
            goal => root
        end select
    end subroutine find_cell
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Place_Cell
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine Place_Cell(goal, partp, n)
        type(CELL), pointer :: goal, temp
        type(point3d) :: partp
        integer :: n
        select case (goal%type)
        case (0)
            goal%type = 1
            goal%part%p = partp
            goal%pos = n
        case (1)
            call crear_subcells(goal)
            call find_cell(goal,temp,partp)
            call place_cell(temp,partp,n)
        case default
            print*,"SHOULD NOT BE HERE. ERROR!"
        end select
    end subroutine Place_Cell
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Crear_Subcells
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Crear_Subcells(goal)
        type(CELL), pointer :: goal
        type(particle3d) :: part
        integer :: i, j, k, n
        integer, dimension(3) :: octant
        
        part = goal%part
        goal%type = 2
        
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    octant = (/i,j,k/)
                    allocate(goal%subcell(i,j,k)%ptr)
                    goal%subcell(i,j,k)%ptr%range%min = Calcular_Range(0, goal, octant)
                    goal%subcell(i,j,k)%ptr%range%max = Calcular_Range(1, goal, octant)
                    
                    if (Belongs(part%p, goal%subcell(i,j,k)%ptr)) then
                        goal%subcell(i,j,k)%ptr%part%p = part%p
                        goal%subcell(i,j,k)%ptr%type = 1
                        goal%subcell(i,j,k)%ptr%pos = goal%pos
                    else
                        goal%subcell(i,j,k)%ptr%type = 0
                    end if
                    call Nullify_Pointers(goal%subcell(i,j,k)%ptr)
                end do
            end do
        end do
    end subroutine Crear_Subcells
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Nullify_Pointers
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Nullify_Pointers(goal)
        type(CELL), pointer :: goal
        integer :: i, j, k
        
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    nullify(goal%subcell(i,j,k)%ptr)
                end do
            end do
        end do
    end subroutine Nullify_Pointers
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Belongs
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function Belongs (partp, goal)
        type(point3d) :: partp
        type(CELL), pointer :: goal
        logical :: Belongs
        
        if (partp%x >= goal%range%min(1) .and. &
            partp%x <= goal%range%max(1) .and. &
            partp%y >= goal%range%min(2) .and. &
            partp%y <= goal%range%max(2) .and. &
            partp%z >= goal%range%min(3) .and. &
            partp%z <= goal%range%max(3)) then
            Belongs = .true.
        else
            Belongs = .false.
        end if
    end function Belongs
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Calcular_Range
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function Calcular_Range (what, goal, octant)
        integer :: what, n
        type(CELL), pointer :: goal
        integer, dimension(3) :: octant
        real(dp), dimension(3) :: Calcular_Range, valor_medio
        
        valor_medio = (goal%range%min + goal%range%max) / 2.0
        
        select case (what)
        case (0)
            where (octant == 1)
                Calcular_Range = goal%range%min
            elsewhere
                Calcular_Range = valor_medio
            endwhere
        case (1)
            where (octant == 1)
                Calcular_Range = valor_medio
            elsewhere
                Calcular_Range = goal%range%max
            endwhere
        end select
    end function Calcular_Range
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Borrar_empty_leaves
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine Borrar_empty_leaves(goal)
        type(CELL), pointer :: goal
        integer :: i, j, k
        
        if (associated(goal%subcell(1,1,1)%ptr)) then
            do i = 1, 2
                do j = 1, 2
                    do k = 1, 2
                        call Borrar_empty_leaves(goal%subcell(i,j,k)%ptr)
                        if (goal%subcell(i,j,k)%ptr%type == 0) then
                            deallocate (goal%subcell(i,j,k)%ptr)
                        end if
                    end do
                end do
            end do
        end if
    end subroutine Borrar_empty_leaves
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Borrar_tree
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine Borrar_tree(goal)
        type(CELL), pointer :: goal
        integer :: i, j, k
        
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    if (associated(goal%subcell(i,j,k)%ptr)) then
                        call Borrar_tree(goal%subcell(i,j,k)%ptr)
                        deallocate (goal%subcell(i,j,k)%ptr)
                    end if
                end do
            end do
        end do
    end subroutine Borrar_tree
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Calculate_masses
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine Calculate_masses(goal)
        type(CELL), pointer :: goal
        integer :: i, j, k
        real(dp), dimension(3) :: c_o_m
        real(dp) :: m_temp
        
        goal%part%m = 0.0_dp
        goal%c_o_m = point3d(0.0_dp, 0.0_dp, 0.0_dp)

        select case (goal%type)
        case (1)
            goal%part%m = part(goal%pos)%m
            goal%c_o_m = part(goal%pos)%p
        case (2)
            do i = 1, 2
                do j = 1, 2
                    do k = 1, 2
                        if (associated(goal%subcell(i,j,k)%ptr)) then
                            call Calculate_masses(goal%subcell(i,j,k)%ptr)
                            m_temp = goal%part%m
                            goal%part%m = m_temp + goal%subcell(i,j,k)%ptr%part%m
                            goal%c_o_m = (m_temp * goal%c_o_m + &
                                goal%subcell(i,j,k)%ptr%part%m * goal%subcell(i,j,k)%ptr%c_o_m) / goal%part%m
                        end if
                    end do
                end do
            end do
        end select
    end subroutine Calculate_masses
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Calculate_forces
    !!
    !!
    !!
    !! Deviation from Solution 1. In this     !!
    !! case we will not calculate the forces  !!
    !! exerted on a reduced number of         !!
    !! particles, but we will calculate how   !!
    !! our tree (which only represents        !!
    !! particles that fall in our octant)     !!
    !! affects all particles                  !!
    !!
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine Calculate_forces(head)
        type(CELL), pointer :: head
        integer :: i, j, k, start, end
        
        do i = 1, n
            call Calculate_forces_aux(i, head)
        end do
    end subroutine Calculate_forces
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Calculate_forces_aux
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine Calculate_forces_aux(number, tree)
        type(CELL), pointer :: tree
        integer :: i, j, k, number
        real(dp) :: l, rs, r3
        type(vector3d) :: rji
        
        select case (tree%type)
        case (1)
            if (number .ne. tree%pos) then
               rji = tree%c_o_m - part(number)%p
                rs = distance(tree%c_o_m, part(number)%p)
                r3 = rs**3
                a(number) = a(number) + part(tree%pos)%m * rji / r3
            end if
        case (2)
            l = tree%range%max(1) - tree%range%min(1) !! Range has the same span in 3 dimensions
            rji = tree%c_o_m - part(number)%p
            rs = distance(tree%c_o_m, part(number)%p)
            
            if (l/rs < theta) then
                !! If cluster, we need to check if l/D < theta
                r3 = rs**3
                a(number) = a(number) + tree%part%m * rji / r3
            else
                do i = 1, 2
                    do j = 1, 2
                        do k = 1, 2
                            if (associated(tree%subcell(i,j,k)%ptr)) then
                                call Calculate_forces_aux(number, tree%subcell(i,j,k)%ptr)
                            end if
                        end do
                    end do
                end do
            end if
        end select
    end subroutine Calculate_forces_aux

    subroutine main_loop(part, a,total_a, n, dt, dt_out, t_end)
        type(particle3d), intent(inout) :: part(:)
        type(vector3d), intent(inout) :: a(:)
        type(vector3d), intent(inout):: total_a(:)
        integer, intent(in) :: n
        real(dp), intent(in) :: dt, dt_out, t_end
        real(dp) :: t, t_out
        integer :: i
        
        ! We create the output file
        open(unit=12, file='output.dat', status='replace', action='write')

        t_out = 0.0

        ! If we change the number of particles the code could be slower, therefore it is convenient to write messages in order to clarify that all is 
        ! running smoothly
        if (my_rank == 0) then
            print *, "Starting the simulation..."
            print *, "Number of particles:", n
            print *, "Total time:", t_end, "Timestep:", dt
            print *, "---------------------------------------------"
        end if

        do t = 0.0, t_end, dt
            
            if (my_rank == 0) then
            ! With the following we print the time every 1.0 time units
                if (mod(t, 1.0_dp) < dt) then
                    print *, "t = ", t
                end if
            end if
    
            do i= 1, n
                part(i)%v = part(i)%v + total_a(i) * dt/2.0_dp
                part(i)%p = part(i)%p + part(i)%v * dt
            end do 
    
        
            !! Positions have changed, so we need to delete
            !! and reinitialize the tree
            call Borrar_tree(head)
            call Calculate_ranges(head)
            head%type = 0
            call Nullify_Pointers(head)
        
            do i = 1, n
                if (Belongs(part(i)%p, head)) then
                    call Find_Cell(head, temp_cell, part(i)%p)
                    call Place_Cell(temp_cell, part(i)%p, i)
                end if
            end do
        
            call Borrar_empty_leaves(head)
            call Calculate_masses(head)

            do i = 1, n
                a(i) = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
            end do

            call Calculate_forces(head)
            call MPI_ALLREDUCE(a, total_a, n*3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, error)
            
            do i= 1, n
                part(i)%v = part(i)%v + a(i) * dt/2.0_dp
            end do
        
            !! We only print if we are the master
            
            if (my_rank == 0) then
                t_out = t_out + dt
                if (t_out >= dt_out) then
                    write(12,'(ES20.10)', advance ='no') t ! With, "advance=no" we avoid to make a split of line
                    do i = 1,n
                        write(12,'(3(1X,ES20.10))', advance='no') part(i)%p
                    end do
                    write(12,*) 
                    t_out = 0.0
                end if
            end if
        end do
    
        call MPI_Finalize (error)

    end subroutine main_loop

end module tree_mpi