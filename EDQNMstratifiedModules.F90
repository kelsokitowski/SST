module EDQNMstratifiedModules
  use integrationModules
  use kernel_interp
  implicit none



  contains

    subroutine getHvals(kLength,E,ET,EHdir,EHpol,ETH,HDIR,HPOL,HT)
      double precision, intent(in) :: E(:),ET(:),EHdir(:),EHpol(:),ETH(:)
      integer, intent(in) :: kLength
      integer :: kj
      double precision, intent(out):: HDIR(kLength),HPOL(kLength),HT(kLength)
      do kj = 1,kLength
    if ( (E(kj)) > 0) then
        HDIR(kj) = EHdir(kj)/E(kj);
        HPOL(kj) = EHpol(kj)/E(kj);
    else
        HDIR(kj)=0.0;
        HPOL(kj)=0.0;
    endif  !if Enew(kj)~=0
    if (ET(kj) > 0) then
        HT(kj) = ETH(kj)/ET(kj);
    else
        HT(kj)=0.0;
    endif  ! if ETnew(kj) ~= 0
    end do
  end subroutine getHvals


! TwoD_MidpointTest.f90
subroutine TwoD_MidpointTest(N,nu,D,kj,k,pVals,HPOL,HDIR,HT,E,F,ET,mu1,mu3,t, &
                              weight,triadFlag,outsideCutCell,Q11,centroidX,centroidY, &
                              insideCutCell,ExtForcing,t0,forcing)
    implicit none
    ! Scalar inputs
    double precision, intent(in) :: N, nu, D, k, t, t0
    integer, intent(in) :: kj

    ! Array inputs
    double precision, dimension(:), intent(in) :: pVals, HPOL, HDIR, HT, E, F, ET, mu1, mu3, ExtForcing
    double precision, dimension(:,:,:), intent(in) :: weight, centroidX, centroidY
    integer, dimension(:,:,:), intent(in) :: triadFlag, outsideCutCell
    integer, dimension(:,:,:,:), intent(in) :: Q11
    integer, dimension(:,:,:), intent(in) :: insideCutCell

    ! Output
    double precision, dimension(6), intent(out) :: forcing

    ! Local variables
    double precision :: a
    double precision :: S_NL_ISO, S_NL_DIR, S_NL_POL, ST_NL_ISO, ST_NL_DIR, SF_NL
    double precision :: phi33
    integer :: kLength, pj, qj
    double precision :: p, q, x, y, z
    double precision, dimension(:), allocatable :: E0, E0T, EF

    ! Kernel variables
    double precision :: kernel1, kernel4
    double precision :: kernel21, kernel22, kernel23, kernel24, kernel25
    double precision :: kernel31, kernel32, kernel33, kernel34, kernel35, kernel36, kernel37
    double precision :: kernel51, kernel52, kernel53, kernel54
    double precision :: kernel61, kernel62, kernel63, kernel64, kernel65, kernel66

    ! Function declarations
    double precision :: thetaCheck

    a = -1.0

    S_NL_ISO = 0.0d0
    S_NL_DIR = 0.0d0
    S_NL_POL = 0.0d0
    ST_NL_ISO = 0.0d0
    ST_NL_DIR = 0.0d0
    SF_NL = 0.0d0

    kLength = size(pVals)

    ! Allocate and compute E0, E0T, EF
    allocate(E0(kLength))
    allocate(E0T(kLength))
    allocate(EF(kLength))

    E0 = (E / (pVals**2.)) / (4.0d0 * acos(-1.0d0))
    E0T = (ET / (pVals**2.)) / (4.0d0 * acos(-1.0d0))
    EF = (F / (pVals**2.)) / (4.0d0 * acos(-1.0d0))


    if (any(E0T /= E0T)) then
       print *, "NaN detected in E0T; first index = ", minloc(E0T, MASK=(E0T /= E0T))
       stop
    end if
    if (any(E0 /= E0)) then
       print *, "NaN detected in E0; first index = ", minloc(E0T, MASK=(E0 /= E0))
       stop
    end if
    if (any(EF /= EF)) then
       print *, "NaN detected in Nv_1; first index = ", minloc(E0T, MASK=(EF /= EF))
       stop
    end if

    do pj = 1, size(pVals)
        do qj = 1, kLength
            p = pVals(pj)
            q = pVals(qj)

            if (triadFlag(pj,qj,kj) == 1) then
                if (outsideCutCell(pj,qj,kj) == 0 .and. (insideCutCell(pj,qj,kj)==0)) then
                    ! Cell center inside domain, uncut
                    kernel1 = E0(qj) * (E0(pj) - E0(kj))

                    kernel21 = E0(qj) * (E0(pj) - E0(kj)) * HPOL(qj)
                    kernel22 = E0(qj) * E0(pj) * HPOL(pj)
                    kernel23 = E0(qj) * (E0(pj) - E0(kj)) * HDIR(qj)
                    kernel24 = E0(qj) * E0(pj) * HDIR(pj)
                    kernel25 = E0(qj) * E0(kj) * HDIR(kj)

                    kernel31 = E0(qj) * E0(pj) * HPOL(pj)
                    kernel32 = E0(qj) * E0(kj) * HPOL(kj)
                    kernel33 = E0(qj) * (E0(pj) - E0(kj)) * HPOL(qj)
                    kernel34 = E0(qj) * E0(pj) * HPOL(pj)
                    kernel35 = E0(qj) * E0(kj) * HPOL(qj)
                    kernel36 = E0(qj) * (E0(pj) - E0(kj)) * HDIR(qj)
                    kernel37 = E0(qj) * E0(pj) * HDIR(pj)

                    kernel4 = E(qj) * (k**2. * ET(pj) - p**2. * ET(kj))

                    kernel51 = E0(qj) * (E0T(pj) - E0T(kj)) * HPOL(qj)
                    kernel52 = E0(qj) * (E0T(pj) - E0T(kj)) * HDIR(qj)
                    kernel53 = E0(qj) * E0T(pj) * HT(pj)
                    kernel54 = E0(qj) * 2.0d0 * E0T(kj) * HT(kj)

                    kernel61 = E0(pj) * EF(qj)
                    kernel62 = E0(pj) * EF(kj)
                    kernel63 = E0(kj) * EF(pj)
                    kernel64 = E0(kj) * EF(qj)
                    kernel65 = E0(qj) * EF(pj)
                    kernel66 = E0(qj) * EF(kj)

                    call waveCosines(k,p,q,x,y,z)
                else

                   if (Q11(1,pj,qj,kj) < 1 .or. Q11(1,pj,qj,kj) > kLength) then
                      print *, "ERROR: Q11 p-index out of range at kj=", kj
                      stop
                   end if
                   if (Q11(2,pj,qj,kj) < 1 .or. Q11(2,pj,qj,kj) > kLength) then
                      print *, "ERROR: Q11 q-index out of range at kj=", kj
                      stop
                   end if
                    ! Cell center cut, use centroid and interpolation
                    kernel1 = interpolateKernel(E0,E0T,centroidX(pj,qj,kj),centroidY(pj,qj,kj), &
                                                k,pVals,pj,qj,Q11,1,kj,HPOL)

                    call interpolateKernel2(E0,centroidX(pj,qj,kj),centroidY(pj,qj,kj), &
                                           pVals,pj,qj,Q11,kj,HPOL,HDIR, &
                                           kernel21,kernel22,kernel23,kernel24,kernel25)

                    call interpolateKernel3(E0,centroidX(pj,qj,kj),centroidY(pj,qj,kj), &
                                           pVals,pj,qj,Q11,kj,HPOL,HDIR, &
                                           kernel31,kernel32,kernel33,kernel34,kernel35,kernel36,kernel37)

                    kernel4 = interpolateKernel(E,ET,centroidX(pj,qj,kj),centroidY(pj,qj,kj), &
                                                k,pVals,pj,qj,Q11,4,kj,HPOL)

                    kernel51 = interpolateKernel(E0,E0T,centroidX(pj,qj,kj),centroidY(pj,qj,kj), &
                                                 k,pVals,pj,qj,Q11,51,kj,HPOL)
                    kernel52 = interpolateKernel(E0,E0T,centroidX(pj,qj,kj),centroidY(pj,qj,kj), &
                                                 k,pVals,pj,qj,Q11,51,kj,HDIR)
                    kernel53 = interpolateKernel(E0,E0T,centroidX(pj,qj,kj),centroidY(pj,qj,kj), &
                                                 k,pVals,pj,qj,Q11,51,kj,HT)
                    kernel54 = interpolateKernel(E0,E0T,centroidX(pj,qj,kj),centroidY(pj,qj,kj), &
                                                 k,pVals,pj,qj,Q11,51,kj,HT)

                    call interpolateKernel6(E0,EF,centroidX(pj,qj,kj),centroidY(pj,qj,kj), &
                                           pVals,pj,qj,Q11,kj, &
                                           kernel61,kernel62,kernel63,kernel64,kernel65,kernel66)

                    call waveCosines(k,centroidX(pj,qj,kj),centroidY(pj,qj,kj),x,y,z)
                    p = centroidX(pj,qj,kj)
                    q = centroidY(pj,qj,kj)

       if ( (kj == 40) .and. (pj == 45) .and. (qj == 40) ) then
             print *, kernel1, centroidX(pj,qj,kj),centroidY(pj,qj,kj)
             endif
          if (kernel1/=kernel1) then
             print *, kernel1
             stop
          endif




           if (kernel21/=kernel21) then
             print *, kernel21
             stop
          endif
           if (kernel22/=kernel22) then
             print *, kernel22
             stop
          endif
           if (kernel23/=kernel23) then
             print *, kernel23
             stop
          endif
           if (kernel24/=kernel24) then
             print *, kernel24
             stop
          endif
           if (kernel25/=kernel25) then
             print *, kernel25
             stop
          endif

                end if

                if (weight(pj,qj,kj)/= weight(pj,qj,kj)) then
                   print *, 'weight'
                   stop
                endif

                thetaCheck = theta(nu,k,p,q,kj,pj,qj,t,mu1);



                ! eqn 1
                S_NL_ISO = S_NL_ISO + weight(pj,qj,kj) * &
                           (theta(nu,k,p,q,kj,pj,qj,t,mu1) * 16.0d0 * (acos(-1.0d0)**2.) * &
                            p**2. * k**2. * q * (x*y + z**3.) * kernel1)

                ! eqn 2
                S_NL_DIR = S_NL_DIR + weight(pj,qj,kj) * &
                           (theta(nu,k,p,q,kj,pj,qj,t,mu1) * 4.0d0 * (acos(-1.0d0)**2.) * &
                            p**2. * k**2. * q * ((y**2. - 1.0d0) * (x*y + z**3.) * kernel21 + &
                            z * (1.0d0 - z**2)**2. * kernel22) + &
                            theta(nu,k,p,q,kj,pj,qj,t,mu1) * 8.0d0 * (acos(-1.0d0)**2.) * &
                            p**2. * k**2. * q * (x*y + z**3.) * &
                            ((3.0d0*y**2. - 1.0d0) * kernel23 + (3.0d0*z**2. - 1.0d0) * kernel24 - &
                             2.0d0 * kernel25))

                ! eqn 3
                S_NL_POL = S_NL_POL + weight(pj,qj,kj) * &
                           (theta(nu,k,p,q,kj,pj,qj,t,mu1) * 4.0d0 * (acos(-1.0d0)**2) * &
                            p**2. * k**2. * q * ((x*y + z**3.) * ((1.0d0 + z**2.) * kernel31 - &
                            4.0d0 * kernel32) + z * (z**2. - 1.0d0) * (1.0d0 + y**2.) * kernel33 + &
                            2.0d0 * z * (z**2. - y**2.) * kernel34 + &
                            2.0d0 * y * x * (z**2. - 1.0d0) * kernel35) + &
                            theta(nu,k,p,q,kj,pj,qj,t,mu1) * 24.0d0 * (acos(-1.0d0)**2.) * &
                            p**2. * k**2. * q * z * (z**2. - 1.0d0) * &
                            ((y**2. - 1.0d0) * kernel36 + (z**2. - 1.0d0) * kernel37))

                ! eqn 4
                ST_NL_ISO = ST_NL_ISO + weight(pj,qj,kj) * &
                            (thetaT(nu,D,k,p,q,t-t0,mu3,qj) * k / p / q * (1.0d0 - y**2.) * kernel4)

                ! eqn 5
                ST_NL_DIR = ST_NL_DIR + weight(pj,qj,kj) * &
                            (4.0d0 * thetaT(nu,D,k,p,q,t-t0,mu3,qj) * (acos(-1.0d0)**2.) * &
                             k**2. * p**2. * q * (x*y + z) * (y**2. - 1.0d0) * kernel51 + &
                             8.0d0 * thetaT(nu,D,k,p,q,t-t0,mu3,qj) * (acos(-1.0d0)**2.) * &
                             k**2. * p**2. * q * (x*y + z) * (3.0d0*y**2. - 1.0d0) * kernel52 + &
                             8.0d0 * thetaT(nu,D,k,p,q,t-t0,mu3,qj) * (acos(-1.0d0)**2.) * &
                             k**2. * p**2. * q * (x*y + z) * &
                             ((3.0d0*z**2. - 1.0d0) * kernel53 - kernel54))

                ! eqn 6
                SF_NL = SF_NL + weight(pj,qj,kj) * &
                        (4.0d0 * (acos(-1.0d0)**2.) * thetaF(nu,D,k,p,q,t-t0,mu3,pj,qj) * &
                         k**2. * p * q * (k * kernel61 * &
                         (1.0d0 + y**2. - z**2. - x*y*z - 2.0d0*y**2.*z**2.) - &
                         2.0d0 * q * (y**3. + x*z) * kernel62) + &
                         4.0d0 * (acos(-1.0d0)**2.) * &
                         thetaF(nu,D,p,k,q,t-t0,mu3,kj,qj) * k**2 * p * q * &
                         ((q * z * (2.0d0*x*y**2. + y*z - x) * kernel63 - &
                           p * y * (x + y*z) * kernel64) + &
                          k * ((1.0d0 - y**2. + z**2. - x*y*z - 2.0d0*y**2.*z**2.) * kernel65 - &
                               2.0d0 * (1.0d0 - y**2.) * kernel66)))

                if (ST_NL_ISO/=ST_NL_ISO) then
                   print *, 'isnan ST_NL_ISO TwoD_MidpointTest'
                   print *, "  mu1(kj)=", mu1(kj), " mu3(kj)=", mu3(kj)
                   print *, 'mul1(pj)=',mu1(pj),' mu1(qj)= ',mu1(qj)
                   print *, 'pj,qj,kj =',pj,qj,kj
                   print *, 'x,y,z=',x,y,z
                   print *, 't=',t
                   print *, 'theta = ',thetaCheck
                   print *, 'D =',D,' nu = ',nu,' p = ',p,' q = ',q,' k = ',k
                   print *, 't0 = ',t0,' t= ',t
                   print *, 'mu1=',mu1
                   stop
                end if

                if (kj == 40 .and. pj == 45 .and. qj == 40) then
    print *, 'Debug at kj=40, pj=45, qj=40'
    print *, 'k, p, q =', k, p, q
    print *, 'x, y, z =', x, y, z
    print *, 'kernel4 =', kernel4
    print *, 'kernel51-54 =', kernel51, kernel52, kernel53, kernel54
    print *, 'kernel61-66 =', kernel61, kernel62, kernel63, kernel64, kernel65, kernel66
    print *, 'thetaT(nu,D,k,p,q,t-t0,mu3,qj) =', thetaT(nu,D,k,p,q,t-t0,mu3,qj)
    print *, 'thetaF(nu,D,k,p,q,t-t0,mu3,pj,qj) =', thetaF(nu,D,k,p,q,t-t0,mu3,pj,qj)
    print *, 'thetaF(nu,D,p,k,q,t-t0,mu3,kj,qj) =', thetaF(nu,D,p,k,q,t-t0,mu3,kj,qj)
    print *, 'mu3(qj) =',mu3(qj),' mu3(pj) = ',mu3(pj),' mu3(kj) = ',mu3(kj)
    print *, 't0 = ',t0,' t= ',t,' t-t0= ',(t-t0)
    print *, 'weight = ',weight(40,45,40)
endif


            end if
        end do  ! end qj loop
    end do  ! end pj loop

    phi33 = 2.0d0 * E(kj) * (1.0d0/3.0d0 + HDIR(kj) + HPOL(kj))

    forcing(1) = S_NL_ISO + a * N * F(kj) + ExtForcing(kj)
    forcing(2) = S_NL_DIR + a * N / 15.0d0 * F(kj)
    forcing(3) = S_NL_POL + a * 2.0d0 * N / 5.0d0 * F(kj)
    forcing(4) = ST_NL_ISO + 2.0d0 * N * F(kj)
    forcing(5) = ST_NL_DIR + 2.0d0 / 15.0d0 * N * F(kj)
    forcing(6) = SF_NL + N * phi33 + a * 2.0d0 * N * ET(kj) * (1.0d0/3.0d0 + HT(kj))

    if (kj == 40) then
       print *, 'forcing =',forcing
       print *, 't0 = ',t0
       print *, 'mu3 = ',mu3

       stop
       endif

    deallocate(E0, E0T, EF)
end subroutine TwoD_MidpointTest







    subroutine getMu(E,kvals,muInt)
    double precision, intent(in) :: E(:), kvals(:)
    double precision, intent(out) :: muInt(:)
    double precision, allocatable :: df(:), partialSums(:),xVals(:)
    integer :: n

    n = size(E)

    allocate(df(n+1))
    allocate(xVals(n+1))
    allocate(partialSums(n))

    df(1) = 0.0
    df(2:) = (kvals**2.) * E
    xVals(1) = 0.0
    xVals(2:) = kVals
    call cumTrapz(xVals, df, partialSums)

    muInt = sqrt( partialSums )
end subroutine getMu


    double precision function theta(nu,k,p,q,kj,pj,qj,t,mu1)
      double precision, intent(in) :: nu,k,p,q,t
      double precision, intent(in) :: mu1(:)
      integer, intent(in) :: kj, pj, qj

      theta = (1.-exp(-(nu*(k**2.+p**2.+q**2.) + mu1(kj)+mu1(pj)+mu1(qj) )*t) )/ &
     (nu*(k**2.+p**2.+q**2.) + mu1(kj)+mu1(pj)+mu1(qj) );

      end function theta

    double precision function thetaT(nu,D,k,p,q,t,mu3,qj)
      double precision, intent(in) :: nu,D,k,p,q,t
      double precision, intent(in) :: mu3(:)
      integer, intent(in) ::qj

      thetaT = ( 1.-exp( -(D*(k**2.+p**2.) +nu*q**2. + mu3(qj) )*t ))/ &
    (D*(k**2.+p**2.) + nu*q**2. +mu3(qj) );

    end function thetaT

    double precision function thetaF(nu,D,k,p,q,t,mu3,pj,qj)
       double precision, intent(in) :: nu,D,k,p,q,t
      double precision, intent(in) :: mu3(:)
      integer, intent(in) ::pj,qj

      thetaF = ( 1.-exp( -(D*k**2.+nu*(p**2. +q**2.) + mu3(pj)+ mu3(qj) )*t ) )/ &
    ( D*k**2.+nu*(p**2.+q**2.) +mu3(qj)+ mu3(pj));

    end function thetaF

    subroutine waveCosines(k,p,q,x,y,z)
      double precision, intent(in):: k,p,q
      double precision, intent(out):: x, y, z
      !cosines of interior angles of triads
      x = -(k**2.-p**2.-q**2.)/(2.*p*q);
      y = -(p**2.-k**2.-q**2.)/(2.*k*q);
      z = -(q**2.-k**2.-p**2.)/(2.*k*p);
    end subroutine waveCosines


subroutine makeTheOutput(Enew,EHdirNew,EHpolNew,ETHnew,ETnew,Fnew,t,Froude,tStar)

    double precision, intent(in) :: Enew(:),EHdirNew(:),EHpolNew(:),ETHnew(:),ETnew(:),Fnew(:)
    double precision, intent(in) :: t,Froude,tStar
    integer :: i

    ! E
    open(unit=10, file="E_stored.csv", status="unknown", position="append")
    write(10,'(E18.8E3)') (Enew(i), i=1,size(Enew))
    close(10)

    ! EHdir
    open(unit=11, file="EHdir_stored.csv", status="unknown", position="append")
    write(11,'(E18.8E3)') (EHdirNew(i), i=1,size(EHdirNew))
    close(11)

    ! EHpol
    open(unit=12, file="EHpol_stored.csv", status="unknown", position="append")
    write(12,'(E18.8E3)') (EHpolNew(i), i=1,size(EHpolNew))
    close(12)

    ! ETH
    open(unit=13, file="ETH_stored.csv", status="unknown", position="append")
    write(13,'(E18.8E3)') (ETHnew(i), i=1,size(ETHnew))
    close(13)

    ! ET
    open(unit=14, file="ET_stored.csv", status="unknown", position="append")
    write(14,'(E18.8E3)') (ETnew(i), i=1,size(ETnew))
    close(14)

    ! F
    open(unit=15, file="F_stored.csv", status="unknown", position="append")
    write(15,'(E18.8E3)') (Fnew(i), i=1,size(Fnew))
    close(15)

    ! Time values
    open(unit=16, file="tValuesStored.csv", status="unknown", position="append")
    write(16,'(E18.8E3)') t
    close(16)

    open(unit=17, file="Froude.csv", status="unknown", position="append")
    write(17,'(E18.8E3)') Froude
    close(17)

    open(unit=18, file="tStar.csv", status="unknown", position="append")
    write(18,'(E18.8E3)') tStar
    close(18)

end subroutine makeTheOutput

 subroutine errorOutput(E,EHdir,EHpol,ETH,ET,F,t,Froude,tStar)

    double precision, intent(in) :: E(:),EHdir(:),EHpol(:),ETH(:),ET(:),F(:)
    double precision, intent(in) :: t,Froude,tStar
    integer :: i

    ! E crash
    open(unit=19, file="Ecrash.csv", status="unknown", position="append")
    write(19,'(E18.8E3)') (E(i), i=1,size(E))
    close(19)

    ! EHdir
    open(unit=20, file="EHdirCrash.csv", status="unknown", position="append")
    write(20,'(E18.8E3)') (EHdir(i), i=1,size(EHdir))
    close(20)

    ! EHpol
    open(unit=21, file="EHpolCrash.csv", status="unknown", position="append")
    write(21,'(E18.8E3)') (EHpol(i), i=1,size(EHpol))
    close(21)

    ! ETH
    open(unit=22, file="ETHcrash.csv", status="unknown", position="append")
    write(22,'(E18.8E3)') (ETH(i), i=1,size(ETH))
    close(22)

    ! ET
    open(unit=23, file="ETcrash.csv", status="unknown", position="append")
    write(23,'(E18.8E3)') (ET(i), i=1,size(ET))
    close(23)

    ! F
    open(unit=24, file="Fcrash.csv", status="unknown", position="append")
    write(24,'(E18.8E3)') (F(i), i=1,size(F))
    close(24)

    ! Time diagnostics
    open(unit=25, file="tCrash.csv", status="unknown", position="append")
    write(25,'(E18.8E3)') t
    close(25)

    open(unit=26, file="FroudeCrash.csv", status="unknown", position="append")
    write(26,'(E18.8E3)') Froude
    close(26)

    open(unit=27, file="tStarCrash.csv", status="unknown", position="append")
    write(27,'(E18.8E3)') tStar
    close(27)

end subroutine errorOutput


end module EDQNMstratifiedModules
