module timeRoutines
  use EDQNMstratifiedModules
  use mpi
    implicit none
    private
    public :: run_timeloop_mpi

    integer, parameter :: dp = 8

contains

  subroutine run_timeloop_mpi( &
       bf, nu, D, h, &
       kVals, &
       v_1, v_2, v_3, v_4, v_5, v_6, &  ! initial state already loaded
       EX1_1,EX2_1,Q_1,f1_1,f2_1,f3_1, EX1_2,EX2_2,Q_2,f1_2,f2_2,f3_2, &
       EX1_3,EX2_3,Q_3,f1_3,f2_3,f3_3, EX1_4,EX2_4,Q_4,f1_4,f2_4,f3_4, &
       EX1_5,EX2_5,Q_5,f1_5,f2_5,f3_5, EX1_6,EX2_6,Q_6,f1_6,f2_6,f3_6, &
       weight, triadFlag, outsideCutCell, insideCutCell, CxVals, CyVals, Q11, &
       ExtForcing, du, qjMin, qjMax, ySquaredAvg, A1, A3, &
                                ! restart/time bookkeeping (loaded from checkpoint)
       t, tStar, &
                                ! outputs/updated scalars
       nPrint, counter, &
       comm)





    ! ---- MPI
    integer, intent(in) :: comm
    integer :: ierr, rank, nprocs, root


    ! ---- Inputs
    real(dp), intent(in) :: bf, nu, D, h
    real(dp), intent(in) :: kVals(:)
    real(dp), intent(in) :: EX1_1(:),EX2_1(:),Q_1(:),f1_1(:),f2_1(:),f3_1(:)
    real(dp), intent(in) :: EX1_2(:),EX2_2(:),Q_2(:),f1_2(:),f2_2(:),f3_2(:)
    real(dp), intent(in) :: EX1_3(:),EX2_3(:),Q_3(:),f1_3(:),f2_3(:),f3_3(:)
    real(dp), intent(in) :: EX1_4(:),EX2_4(:),Q_4(:),f1_4(:),f2_4(:),f3_4(:)
    real(dp), intent(in) :: EX1_5(:),EX2_5(:),Q_5(:),f1_5(:),f2_5(:),f3_5(:)
    real(dp), intent(in) :: EX1_6(:),EX2_6(:),Q_6(:),f1_6(:),f2_6(:),f3_6(:)
    real(dp), intent(in) :: weight(:,:,:), CxVals(:,:,:), CyVals(:,:,:)
    integer,  intent(in) :: triadFlag(:,:,:), outsideCutCell(:,:,:), insideCutCell(:,:,:), Q11(:,:,:,:)
    integer,  intent(in) :: qjMin, qjMax
    real(dp), intent(in) :: ExtForcing(:)
    real(dp), intent(in) :: du,  ySquaredAvg, A1, A3

    ! ---- Restart/time variables (loaded)
    real(dp), intent(inout) :: t, tStar
    integer,  intent(inout) :: nPrint, counter

    ! ---- State (initialized from v_* below)
    real(dp), intent(inout) :: v_1(:), v_2(:), v_3(:), v_4(:), v_5(:), v_6(:)

    ! ---- Locals
    integer :: kLength, n, stepsNeeded
    real(dp) :: TauL0, t0, tRestart, tStarRestart, t30, storageInterval
    logical :: done

    ! distributed indexing
    integer :: local_kLength, local_kStart, local_kEnd, remainder, kj, kjl, i
    integer, allocatable :: recvcounts(:), displs(:)

    ! fields/scratch
    real(dp), allocatable :: E(:),EHdir(:),EHpol(:),ET(:),ETH(:),F(:)
    real(dp), allocatable :: Enew(:),EHdirNew(:),EHpolNew(:),ETnew(:),ETHnew(:),Fnew(:)
    real(dp), allocatable :: mu(:), mu1(:), mu3(:), HDIR(:), HPOL(:), HT(:)
    real(dp), allocatable :: Nv_1(:),Nv_2(:),Nv_3(:),Nv_4(:),Nv_5(:),Nv_6(:)
    real(dp), allocatable :: Na_1(:),Na_2(:),Na_3(:),Na_4(:),Na_5(:),Na_6(:)
    real(dp), allocatable :: Nb_1(:),Nb_2(:),Nb_3(:),Nb_4(:),Nb_5(:),Nb_6(:)
    real(dp), allocatable :: Nc_1(:),Nc_2(:),Nc_3(:),Nc_4(:),Nc_5(:),Nc_6(:)
    real(dp), allocatable :: a_1(:),a_2(:),a_3(:),a_4(:),a_5(:),a_6(:)
    real(dp), allocatable :: c_1(:),c_2(:),c_3(:),c_4(:),c_5(:),c_6(:)
    real(dp), allocatable :: localN1(:),localN2(:),localN3(:),localN4(:),localN5(:),localN6(:)
    real(dp) :: forcing(6)

    ! diagnostics
    real(dp) :: KE_local, KE, eps_local, epzilon, Froude, Re_l, BuoyancyRe

    !==================== Initialization from MATLAB preamble ====================
    root = 0
    call MPI_Comm_rank(comm, rank, ierr)
    call MPI_Comm_size(comm, nprocs, ierr)
    kLength = size(kVals)

    ! Initialize spectral fields from v_* (your MATLAB preamble)
    allocate(E(kLength),EHdir(kLength),EHpol(kLength),ET(kLength),ETH(kLength),F(kLength))
    E    = v_1; EHdir = v_2; EHpol = v_3
    ET   = v_4; ETH   = v_5; F     = v_6

    allocate(Enew(kLength),EHdirNew(kLength),EHpolNew(kLength),ETnew(kLength),ETHnew(kLength),Fnew(kLength))
    Enew = v_1; EHdirNew = v_2; EHpolNew = v_3; ETnew = v_4; ETHnew = v_5; Fnew = v_6

    ! Time scales & storage cadence (match MATLAB semantics)
    TauL0 = t / tStar                      ! (t/TauL0)/t = TauL0^{-1} -> so TauL0 = t/tStar

    tStarRestart = t / TauL0
    tRestart = t
    t30 = 30.0_dp * TauL0
    print *, 't30=',t30
    storageInterval = 0.2_dp
    n = 0
    stepsNeeded = int( ((t - t30)*bf / h) * 10.0_dp )  ! not used in loop control, for reporting
    print *, 'amount of steps needed for 10 nt is ',stepsNeeded
    ! t0 = 30*TauL0 (constant reference passed to kernels)
    t0 = (30.0 * t) / tStar

    print *, 'TauL0 ',TauL0,' t=',t,' tStar=',tStar,' tRestart=',tRestart,' t0=',t0
    call MPI_Bcast(TauL0, 1, MPI_DOUBLE_PRECISION, root, comm, ierr)
    call MPI_Bcast(t0,     1, MPI_DOUBLE_PRECISION, root, comm, ierr)

    ! Allocate per-step helpers
    allocate(mu(kLength), mu1(kLength), mu3(kLength))
    allocate(HDIR(kLength), HPOL(kLength), HT(kLength))
    allocate(Nv_1(kLength),Nv_2(kLength),Nv_3(kLength),Nv_4(kLength),Nv_5(kLength),Nv_6(kLength))
    allocate(Na_1(kLength),Na_2(kLength),Na_3(kLength),Na_4(kLength),Na_5(kLength),Na_6(kLength))
    allocate(Nb_1(kLength),Nb_2(kLength),Nb_3(kLength),Nb_4(kLength),Nb_5(kLength),Nb_6(kLength))
    allocate(Nc_1(kLength),Nc_2(kLength),Nc_3(kLength),Nc_4(kLength),Nc_5(kLength),Nc_6(kLength))
    allocate(c_1(kLength),c_2(kLength),c_3(kLength),c_4(kLength),c_5(kLength),c_6(kLength))
    allocate(a_1(kLength),a_2(kLength),a_3(kLength),a_4(kLength),a_5(kLength),a_6(kLength))

    ! Uneven block distribution
    local_kLength = kLength / nprocs
    remainder = mod(kLength, nprocs)
    if (rank < remainder) then
       local_kLength = local_kLength + 1
       local_kStart  = rank * local_kLength + 1
    else
       local_kStart  = remainder * (local_kLength + 1) + (rank - remainder) * local_kLength + 1
    end if
    local_kEnd = local_kStart + local_kLength - 1
    allocate(localN1(local_kLength),localN2(local_kLength),localN3(local_kLength), &
         localN4(local_kLength),localN5(local_kLength),localN6(local_kLength))
    allocate(recvcounts(nprocs), displs(nprocs))
    call MPI_Gather(local_kLength,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,root,comm,ierr)
    if (rank==root) then
       displs(1)=0
       do i = 2, nprocs
          displs(i)=displs(i-1)+recvcounts(i-1)
       end do
    end if
    call MPI_Bcast(recvcounts,nprocs,MPI_INTEGER,root,comm,ierr)
    call MPI_Bcast(displs,    nprocs,MPI_INTEGER,root,comm,ierr)

    done = .false.
    print *, 'ExtForcing(16) = ',ExtForcing(16),' on processor rank ',rank
    print *, 'Entering do while time loop, h = ',h
    BuoyancyRe = -987654321.0

    !==================== MAIN TIME LOOP ====================
    do while ( ((t - t30)*bf <= 10.0_dp) .and. (.not. done) )

       call MPI_Barrier(comm, ierr)
       if (rank==root) then
       if (any(E<0))then
          print *, 'E=',E
          stop
          endif
          endif
       ! n-driven time update (NOT counter)
       t = tRestart + n*h

       ! mu, H from current state (root computes; we broadcast via bcast_state)
       call getMu(E, kVals, mu);  mu1 = A1*mu;  mu3 = A3*mu
!!$        print *, 'substep 1 mu1 = ',mu1
!!$        print *, 'substep 1 E(1) = ',E(1),' E(2) = ',E(2)
!!$        print *, 'kVals =(44) = ',kVals(44)
       !print *,'Q11',Q11(1,70,70,:)
       !print *, 'wieghts',weight(70,70,:)
       call getHvals(kLength, E, ET, EHdir, EHpol, ETH, HDIR, HPOL, HT)
       call bcast_state(HPOL,HDIR,HT,E,F,ET,mu1,mu3,t,comm)

       call MPI_Barrier(comm, ierr)
        print *, 'rank ',rank,' made it to ln 177 timeRoutines'
       !---------- Substep 1: Nv at time t ----------
       do kj = 1, local_kLength
          kjl = local_kStart + kj - 1
          call TwoD_MidpointTest(bf,nu,D,kjl,kVals(kjl),kVals,HPOL,HDIR,HT, &
               E,F,ET,mu1,mu3,t,weight, &
               triadFlag,outsideCutCell,Q11,CxVals,CyVals,insideCutCell, &
               ExtForcing,t0,forcing)

          localN1(kj)=forcing(1); localN2(kj)=forcing(2); localN3(kj)=forcing(3)
          localN4(kj)=forcing(4); localN5(kj)=forcing(5); localN6(kj)=forcing(6)
          print *, 'rank ',rank,' made it to ln 188 timeRoutines'

       end do
       call gatherv6(localN1,localN2,localN3,localN4,localN5,localN6, Nv_1,Nv_2,Nv_3,Nv_4,Nv_5,Nv_6, &
            recvcounts, displs, comm)
       call MPI_Barrier(comm, ierr)
       if (any(Nv_1 /= Nv_1)) then
          print *, "NaN detected in Nv_1; first index = ", minloc(Nv_1, MASK=(Nv_1 /= Nv_1))
          stop
       end if
       !======= get a_*
       if (rank==root) then


          a_1 = EX2_1*v_1 + Q_1*Nv_1;  !first substep solution E.
          Enew = a_1;
          a_2 = EX2_2*v_2 + Q_2*Nv_2;  !first substep solution EHdir.
          EHdirNew = a_2;
          a_3 = EX2_3*v_3 + Q_3*Nv_3;  !first substep solution EHpol.
          EHpolNew = a_3;
          a_4 = EX2_4*v_4 + Q_4*Nv_4;  !first substep solution ET.
          ETnew = a_4;
          a_5 = EX2_5*v_5 + Q_5*Nv_5;  !first substep solution ETH.
          ETHnew = a_5;
          a_6 = EX2_6*v_6 + Q_6*Nv_6;  !first substep solution F.
          Fnew = a_6;

          if (any(Enew /= Enew)) then
             print *, 'Enew nan 209 timeRoutines'
             stop

          endif
          if (any(EX2_1 /= EX2_1)) then
             print *, 'EX2_1 nan'
          end if
          if (any(Q_1 /= Q_1)) then
             print *, 'Q_1 nan'
          endif
          if (any(v_1 /= v_1)) then
             print *, 'v_1 nan'
          endif
          if (any(Nv_1 /= Nv_1)) then
             print *, 'Nv_1 nan'
          endif
       end if


       call bcast6(Enew,EHdirNew,EHpolNew,ETnew,ETHnew,Fnew,comm)



!!$               print *, 'Enew 219 = ', Enew
!!$               print *, 'kVals 220 = ', kVals
       !---------- Substep 2: Na at t + h/2 ----------
       call getMu(Enew,kVals,mu); mu1=A1*mu; mu3=A3*mu
!!$        print *, 'substep2 mu1 = ',mu1
       call getHvals(kLength, Enew, ETnew, EHdirNew, EHpolNew, ETHnew, HDIR, HPOL, HT)
       call bcast_state(HPOL,HDIR,HT,Enew,Fnew,ETnew,mu1,mu3,t+0.5*h,comm)

       do kj = 1, local_kLength
          kjl = local_kStart + kj - 1
          call TwoD_MidpointTest(bf,nu,D,kjl,kVals(kjl),kVals,HPOL,HDIR,HT, &
               Enew,Fnew,ETnew,mu1,mu3,t+0.5_dp*h,weight, &
               triadFlag,outsideCutCell,Q11,CxVals,CyVals,insideCutCell, &
               ExtForcing,t0,forcing)
          localN1(kj)=forcing(1); localN2(kj)=forcing(2); localN3(kj)=forcing(3)
          localN4(kj)=forcing(4); localN5(kj)=forcing(5); localN6(kj)=forcing(6)
       end do
       call gatherv6(localN1,localN2,localN3,localN4,localN5,localN6, Na_1,Na_2,Na_3,Na_4,Na_5,Na_6, &
            recvcounts, displs, comm)
       call MPI_Barrier(comm, ierr)
       !========= COMPUTE b_1, b_2...
       Enew     = EX2_1*v_1 + Q_1*Na_1
       EHdirNew = EX2_2*v_2 + Q_2*Na_2
       EHpolNew = EX2_3*v_3 + Q_3*Na_3
       ETnew    = EX2_4*v_4 + Q_4*Na_4
       ETHnew   = EX2_5*v_5 + Q_5*Na_5
       Fnew     = EX2_6*v_6 + Q_6*Na_6
       call bcast6(Enew,EHdirNew,EHpolNew,ETnew,ETHnew,Fnew,comm)

       !---------- Substep 3: Nb at t + h/2 ----------
       call getMu(Enew,kVals,mu); mu1=A1*mu; mu3=A3*mu
       call getHvals(kLength, Enew, ETnew, EHdirNew, EHpolNew, ETHnew, HDIR, HPOL, HT)
       call bcast_state(HPOL,HDIR,HT,Enew,Fnew,ETnew,mu1,mu3,t+0.5_dp*h,comm)
       call MPI_Barrier(comm, ierr)

       do kj = 1, local_kLength
          kjl = local_kStart + kj - 1
          call TwoD_MidpointTest(bf,nu,D,kjl,kVals(kjl),kVals,HPOL,HDIR,HT, &
               Enew,Fnew,ETnew,mu1,mu3,t+0.5_dp*h,weight, &
               triadFlag,outsideCutCell,Q11,CxVals,CyVals,insideCutCell, &
               ExtForcing,t0,forcing)
          localN1(kj)=forcing(1); localN2(kj)=forcing(2); localN3(kj)=forcing(3)
          localN4(kj)=forcing(4); localN5(kj)=forcing(5); localN6(kj)=forcing(6)
       end do
       call gatherv6(localN1,localN2,localN3,localN4,localN5,localN6, Nb_1,Nb_2,Nb_3,Nb_4,Nb_5,Nb_6, &
            recvcounts, displs, comm)
       call MPI_Barrier(comm, ierr)

       !======== COMPUTE c_1...c_6 prestate
       c_1 = EX2_1*a_1 + Q_1*(2.0*Nb_1-Nv_1); !third substep
       Enew = c_1;
       c_2 = EX2_2*a_2 + Q_2*(2.0*Nb_2-Nv_2); !third substep
       EHdirNew  = c_2;
       c_3 = EX2_3*a_3 + Q_3*(2.0*Nb_3-Nv_3); !third substep
       EHpolNew = c_3;
       c_4 = EX2_4*a_4 + Q_4*(2.0*Nb_4-Nv_4); !third substep
       ETnew = c_4;
       c_5 = EX2_5*a_5 + Q_5*(2.0*Nb_5-Nv_5); !third substep
       ETHnew = c_5;
       c_6 = EX2_6*a_6 + Q_6*(2.0*Nb_6-Nv_6); !third substep for F
       Fnew = c_6;
       call bcast6(Enew,EHdirNew,EHpolNew,ETnew,ETHnew,Fnew,comm)

       !---------- Substep 4: Nc at t + h ----------
       call getMu(Enew,kVals,mu); mu1=A1*mu; mu3=A3*mu
       call getHvals(kLength, Enew, ETnew, EHdirNew, EHpolNew, ETHnew, HDIR, HPOL, HT)
       call bcast_state(HPOL,HDIR,HT,Enew,Fnew,ETnew,mu1,mu3,t+h,comm)
       call MPI_Barrier(comm, ierr)

       do kj = 1, local_kLength
          kjl = local_kStart + kj - 1
          call TwoD_MidpointTest(bf,nu,D,kjl,kVals(kjl),kVals,HPOL,HDIR,HT, &
               Enew,Fnew,ETnew,mu1,mu3,t+h,weight, &
               triadFlag,outsideCutCell,Q11,CxVals,CyVals,insideCutCell, &
               ExtForcing,t0,forcing)
          localN1(kj)=forcing(1); localN2(kj)=forcing(2); localN3(kj)=forcing(3)
          localN4(kj)=forcing(4); localN5(kj)=forcing(5); localN6(kj)=forcing(6)
       end do
       call gatherv6(localN1,localN2,localN3,localN4,localN5,localN6, Nc_1,Nc_2,Nc_3,Nc_4,Nc_5,Nc_6, &
            recvcounts, displs, comm)
       call MPI_Barrier(comm, ierr)

       !---------- Full ETDRK update on root ----------
       if (rank==root) then
          v_1 = EX1_1*v_1 + Nv_1*f1_1 + 2.0*(Na_1+Nb_1)*f2_1 + Nc_1*f3_1
          v_2 = EX1_2*v_2 + Nv_2*f1_2 + 2.0*(Na_2+Nb_2)*f2_2 + Nc_2*f3_2
          v_3 = EX1_3*v_3 + Nv_3*f1_3 + 2.0*(Na_3+Nb_3)*f2_3 + Nc_3*f3_3
          v_4 = EX1_4*v_4 + Nv_4*f1_4 + 2.0*(Na_4+Nb_4)*f2_4 + Nc_4*f3_4
          v_5 = EX1_5*v_5 + Nv_5*f1_5 + 2.0*(Na_5+Nb_5)*f2_5 + Nc_5*f3_5
          v_6 = EX1_6*v_6 + Nv_6*f1_6 + 2.0*(Na_6+Nb_6)*f2_6 + Nc_6*f3_6
          E = v_1; EHdir=v_2; EHpol=v_3; ET=v_4; ETH=v_5; F=v_6

       end if
       call bcast6(E,EHdir,EHpol,ET,ETH,F,comm)
       call bcast6(v_1,v_2,v_3,v_4,v_5,v_6,comm)




       !---------- Diagnostics & checkpoint cadence (root) ----------
       if (rank==root) then

          ! KE = 0.5 ∫ E dk ; epzilon = 2 nu ∫ k^2 E dk
          KE  = trapz_full(.5_dp, kVals, E)
          epzilon = 2.0_dp*nu*trapz_full_k2E(kVals, E)
          Froude = KE / max(epzilon,1.0e-300_dp) / bf
          Re_l   = (KE*KE) / max(epzilon*nu, 1.0e-300_dp)
                
                if ( mod(n,200)== 0 ) then
                        print *, 'step n= ',n,' (t-t30)*bf=Nt=',(t-t30)*bf,' bf =',bf,' Pr=',nu/D
                endif

          ! storage cadence using t/TauL0 crossings
          if ( ( t/TauL0 <= tStarRestart + nPrint*storageInterval ) .and. &
               ( (t+h)/TauL0 >= tStarRestart + nPrint*storageInterval ) ) then
             nPrint = nPrint + 1
             tStar  = t / TauL0
             counter = counter + 1
             call write_checkpoint_bin(bf, v_1,v_2,v_3,v_4,v_5,v_6, t, tStar, counter)
!!$             if (mod(counter,10) == 0) then
!!$                call write_checkpoint_bin(bf, v_1,v_2,v_3,v_4,v_5,v_6, t, tStar, counter)
!!$             end if
          end if

          ! steady-state stop test (tunable tolerance)
          if (abs(BuoyancyRe - Froude*Froude*Re_l) < 1.0e-5_dp) then
             done = .true.
             print *, 'Steady state, Re_b = ,',BuoyancyRe,'Froude = ',Froude,'Re_l= ',Re_l
          else
             done = .false.
             BuoyancyRe = Froude*Froude * Re_l
          end if
       end if
       call MPI_Bcast(done, 1, MPI_LOGICAL, root, comm, ierr)

       ! advance n
       n = n + 1
       
       


    end do

    !==================== cleanup ====================
    deallocate(E,EHdir,EHpol,ET,ETH,F)
    deallocate(Enew,EHdirNew,EHpolNew,ETnew,ETHnew,Fnew)
    deallocate(mu,mu1,mu3,HDIR,HPOL,HT)
    deallocate(Nv_1,Nv_2,Nv_3,Nv_4,Nv_5,Nv_6)
    deallocate(Na_1,Na_2,Na_3,Na_4,Na_5,Na_6)
    deallocate(Nb_1,Nb_2,Nb_3,Nb_4,Nb_5,Nb_6)
    deallocate(Nc_1,Nc_2,Nc_3,Nc_4,Nc_5,Nc_6)
    deallocate(localN1,localN2,localN3,localN4,localN5,localN6)
    deallocate(recvcounts, displs)

  contains
    subroutine bcast_state(HPOL,HDIR,HT,E,F,ET,mu1,mu3,tval,comm_)
      real(dp), intent(inout) :: HPOL(:),HDIR(:),HT(:),E(:),F(:),ET(:),mu1(:),mu3(:)
      real(dp), intent(in) :: tval
      integer,  intent(in)    :: comm_
      integer :: ier
      call MPI_Bcast(HPOL, size(HPOL), MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Bcast(HDIR, size(HDIR), MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Bcast(HT,   size(HT),   MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Bcast(E,    size(E),    MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Bcast(F,    size(F),    MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Bcast(ET,   size(ET),   MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Bcast(mu1,  size(mu1),  MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Bcast(mu3,  size(mu3),  MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Bcast(tval, 1,          MPI_DOUBLE_PRECISION, 0, comm_, ier)
    end subroutine bcast_state

    subroutine bcast6(a,b,c,d,e,f,comm_)
      real(dp), intent(inout) :: a(:),b(:),c(:),d(:),e(:),f(:)
      integer,  intent(in)    :: comm_
      integer :: ier
      call MPI_Bcast(a, size(a), MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Bcast(b, size(b), MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Bcast(c, size(c), MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Bcast(d, size(d), MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Bcast(e, size(e), MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Bcast(f, size(f), MPI_DOUBLE_PRECISION, 0, comm_, ier)
    end subroutine bcast6

    subroutine gatherv6(l1,l2,l3,l4,l5,l6, g1,g2,g3,g4,g5,g6, rc, displs, comm_)
      real(dp), intent(in)    :: l1(:),l2(:),l3(:),l4(:),l5(:),l6(:)
      real(dp), intent(inout) :: g1(:),g2(:),g3(:),g4(:),g5(:),g6(:)
      integer,  intent(in)    :: rc(:), displs(:), comm_
      integer :: ier
      call MPI_Gatherv(l1, size(l1), MPI_DOUBLE_PRECISION, g1, rc, displs, MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Gatherv(l2, size(l2), MPI_DOUBLE_PRECISION, g2, rc, displs, MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Gatherv(l3, size(l3), MPI_DOUBLE_PRECISION, g3, rc, displs, MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Gatherv(l4, size(l4), MPI_DOUBLE_PRECISION, g4, rc, displs, MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Gatherv(l5, size(l5), MPI_DOUBLE_PRECISION, g5, rc, displs, MPI_DOUBLE_PRECISION, 0, comm_, ier)
      call MPI_Gatherv(l6, size(l6), MPI_DOUBLE_PRECISION, g6, rc, displs, MPI_DOUBLE_PRECISION, 0, comm_, ier)
    end subroutine gatherv6

    pure function trapz_full(scale, x, y) result(val)
      real(dp), intent(in) :: scale, x(:), y(:)
      real(dp) :: val
      integer :: j, n
      n = size(x)
      val = 0.0_dp
      do j = 1, n-1
         val = val + scale * (x(j+1)-x(j)) * (y(j+1)+y(j))
      end do
    end function trapz_full

    pure function trapz_full_k2E(x, y) result(val)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: val
      integer :: j, n
      n = size(x)
      val = 0.0_dp
      do j = 1, n-1
         val = val + 0.5_dp * (x(j+1)-x(j)) * ( x(j+1)**2*y(j+1) + x(j)**2*y(j) )
      end do
    end function trapz_full_k2E

    subroutine write_checkpoint_bin(bf_, v1, v2, v3, v4, v5, v6, t_, tStar_, counter_)
    real(dp), intent(in) :: bf_, t_, tStar_
    real(dp), intent(in) :: v1(:), v2(:), v3(:), v4(:), v5(:), v6(:)
    integer, intent(in)  :: counter_
    integer :: fid, n, ios
    character(len=128) :: dirname, fname, bfstr

    n = size(v1)

    ! Format the bf parameter
    write(bfstr, '(G0.6)') bf_
    bfstr = adjustl(bfstr)

    ! Create directory path and filename
    write(dirname, '(A,A)') 'results/param_', trim(bfstr)
    fname = trim(dirname) // '/checkpointOut.bin'

    ! Open file
    fid = 127
    open(unit=fid, file=fname, form='unformatted', access='stream', &
         status='replace', iostat=ios)
    if (ios /= 0) then
        print *, 'ERROR: Cannot open file: ', trim(fname)
        stop
    end if

    ! Write kLength as int32
    write(fid) int(n, 4)

    ! Write arrays separately (not concatenated)
    write(fid) v1
    write(fid) v2
    write(fid) v3
    write(fid) v4
    write(fid) v5
    write(fid) v6

    ! Write scalars
    write(fid) t_
    write(fid) tStar_
    write(fid) real(counter_, kind=8)

    close(fid)

    if (counter_ == 1) then
        print '(A,A)', 'Checkpoint written to: ', trim(fname)
    end if

end subroutine write_checkpoint_bin
  end subroutine run_timeloop_mpi

end module timeRoutines
