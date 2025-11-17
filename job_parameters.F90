module job_parameters
    use mpi
    implicit none
    private
    public :: initialize_parameters_mpi
    integer, parameter :: dp = 8

contains

    !-----------------------------------------------------------------
    !  Initialize parameters for an MPI job launched via SLURM array
    !  Root process reads:
    !     - SLURM_ARRAY_TASK_ID  → selects bf
    !     - environment variable "Pr" → reads Prandtl number
    !  Both bf and Pr are broadcast to all ranks.
    !-----------------------------------------------------------------
    subroutine initialize_parameters_mpi(bf, Pr, comm)
        implicit none
        !------------------------------
        ! Arguments
        !------------------------------
        real(dp), intent(out) :: bf, Pr
        integer,  intent(in)  :: comm

        !------------------------------
        ! Locals
        !------------------------------
        integer :: rank, ierr, task_id, nBF, istat
        real(dp), allocatable :: bfValues(:), bfValuesLg(:)
        character(len=256) :: env_str, cmd, output_dir
        real(dp) :: param

        !------------------------------
        ! MPI rank
        !------------------------------
        call MPI_Comm_rank(comm, rank, ierr)

        if (rank == 0) then
            !========================================================
            ! 1. Get SLURM_ARRAY_TASK_ID
            !========================================================
            call get_environment_variable("SLURM_ARRAY_TASK_ID", env_str, status=istat)
            if (istat /= 0) then
                print *, "Error: SLURM_ARRAY_TASK_ID not found."
                call MPI_Abort(comm, 1, ierr)
            endif
            read(env_str, *) task_id

            !========================================================
            ! 2. Construct bfValues array
            !========================================================
            allocate(bfValuesLg(14))
            bfValuesLg = (/ 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp, &
                           10.0_dp, 13.0_dp, 15.0_dp, 16.0_dp, 17.0_dp, 18.0_dp /)

            allocate(bfValues(12 + 3 + 14))   ! 12 linspace + 3 fixed + 14 large
            call linspace(bfValues(1:12), 0.001_dp, 0.01_dp)
            bfValues(13:15) = (/ 0.25_dp, 0.5_dp, 1.0_dp /)
            bfValues(16:)   = bfValuesLg
            nBF = size(bfValues)

            ! Check bounds
            if (task_id < 1 .or. task_id > nBF) then
                print *, "Error: task_id out of bounds for bfValues"
                call MPI_Abort(comm, 2, ierr)
            endif

            ! Select parameter
            bf = bfValues(task_id)
            param = bf
            write(*,'(A,I4,A,F8.3)') "Running SLURM task ", task_id, " with parameter bf = ", bf

            !========================================================
            ! 3. Get Pr from environment variable
            !========================================================
            call get_environment_variable("Pr", env_str, status=istat)
            if (istat /= 0) then
                print *, "Error: Environment variable Pr not found."
                call MPI_Abort(comm, 3, ierr)
            endif
            read(env_str, *) Pr
            write(*,'(A,F8.3)') "Running job with Pr = ", Pr

            !========================================================
            ! 4. Create directory results/param_<bf>
            !========================================================
            write(output_dir, '(A,F6.3)') 'results/param_', param
            output_dir = adjustl(trim(output_dir))
            write(cmd, '(A,A)') 'mkdir -p ', trim(output_dir)
            call execute_command_line(trim(cmd))

            ! Optional: change working directory
            write(cmd, '(A,A)') 'cd ', trim(output_dir)
            call execute_command_line(trim(cmd))

            print *, "Output directory: ", trim(output_dir)
        endif

        !============================================================
        ! Broadcast bf and Pr to all ranks
        !============================================================
        call MPI_Bcast(bf, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
        call MPI_Bcast(Pr, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

        ! Cleanup
        if (allocated(bfValues))   deallocate(bfValues)
        if (allocated(bfValuesLg)) deallocate(bfValuesLg)
    end subroutine initialize_parameters_mpi


    !-------------------------------------------------------------
    ! Utility: linspace(a,b,n)
    !-------------------------------------------------------------
    subroutine linspace(arr, a, b)
        real(dp), intent(out) :: arr(:)
        real(dp), intent(in)  :: a, b
        integer :: i, n
        n = size(arr)
        if (n == 1) then
            arr(1) = a
        else
            do i = 1, n
                arr(i) = a + (b - a) * real(i - 1, dp) / real(n - 1, dp)
            end do
        endif
    end subroutine linspace

end module job_parameters
