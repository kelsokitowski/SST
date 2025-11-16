module kernel_interp
  implicit none


contains

  !=========================
  ! Bilinear interpolation
  !=========================
  function bilinear(centroidX, centroidY, x1, y1, x2, y2, fQ12, fQ22, fQ21, fQ11) result(theResult)
    implicit none
    ! Input variables - all real double precision scalars
    double precision, intent(in) :: centroidX, centroidY, x1, y1, x2, y2
    double precision, intent(in) :: fQ12, fQ22, fQ21, fQ11
    ! Output variable
    double precision :: theResult
    ! Local variables
    double precision :: x, y
    double precision, dimension(2) :: a, c
    double precision, dimension(2,2) :: b
    double precision, dimension(2) :: bc

    ! Assign for clarity (matching MATLAB code)
    x = centroidX
    y = centroidY

    ! Build vector a
    a(1) = x2 - x
    a(2) = x - x1

    ! Build matrix b (column-major order in Fortran)
    b(1,1) = fQ11
    b(2,1) = fQ21
    b(1,2) = fQ12
    b(2,2) = fQ22

    ! Build vector c
    c(1) = y2 - y
    c(2) = y - y1

    ! Compute b*c
    bc = matmul(b, c)
    if (((x2 - x1) * (y2 - y1)) <= 0.0000000000000000001d0) then
       print *, 'bilinear denominator zero, error'
       stop
       endif

    ! Compute final result: a*(b*c) / ((x2-x1)*(y2-y1))
    theResult = dot_product(a, bc) / ((x2 - x1) * (y2 - y1))

end function bilinear

!!$! waveCosines.f90
!!$subroutine waveCosines(k,p,q,x,y,z)
!!$    implicit none
!!$    real(kind=8), intent(in) :: k, p, q
!!$    real(kind=8), intent(out) :: x, y, z
!!$
!!$    ! cosines of interior angles of triads
!!$    x = -(k**2 - p**2 - q**2) / (2.0d0 * p * q)
!!$    y = -(p**2 - k**2 - q**2) / (2.0d0 * k * q)
!!$    z = -(q**2 - k**2 - p**2) / (2.0d0 * k * p)
!!$end subroutine waveCosines
!!$
!!$
!!$! thetaT.f90
!!$function thetaT(nu,D,k,p,q,t,mu3,qj) result(thetaTval)
!!$    implicit none
!!$    real(kind=8), intent(in) :: nu, D, k, p, q, t
!!$    real(kind=8), dimension(:), intent(in) :: mu3
!!$    integer, intent(in) :: qj
!!$    real(kind=8) :: thetaTval
!!$
!!$    thetaTval = (1.0d0 - exp(-(D * (k**2 + p**2) + nu * q**2 + mu3(qj)) * t)) / &
!!$                (D * (k**2 + p**2) + nu * q**2 + mu3(qj))
!!$end function thetaT
!!$
!!$
!!$! thetaF.f90
!!$function thetaF(nu,D,k,p,q,t,mu3,pj,qj) result(output)
!!$    implicit none
!!$    real(kind=8), intent(in) :: nu, D, k, p, q, t
!!$    real(kind=8), dimension(:), intent(in) :: mu3
!!$    integer, intent(in) :: pj, qj
!!$    real(kind=8) :: output
!!$
!!$    output = (1.0d0 - exp(-(D * k**2 + nu * (p**2 + q**2) + mu3(pj) + mu3(qj)) * t)) / &
!!$             (D * k**2 + nu * (p**2 + q**2) + mu3(qj) + mu3(pj))
!!$end function thetaF
!!$
!!$
!!$! theta.f90
!!$function theta(nu,k,p,q,kj,pj,qj,t,mu1) result(thetaVal)
!!$    implicit none
!!$    real(kind=8), intent(in) :: nu, k, p, q, t
!!$    real(kind=8), dimension(:), intent(in) :: mu1
!!$    integer, intent(in) :: kj, pj, qj
!!$    real(kind=8) :: thetaVal
!!$
!!$    thetaVal = (1.0d0 - exp(-(nu * (k**2 + p**2 + q**2) + mu1(kj) + mu1(pj) + mu1(qj)) * t)) / &
!!$               (nu * (k**2 + p**2 + q**2) + mu1(kj) + mu1(pj) + mu1(qj))
!!$end function theta


! kernel61Val.f90
function kernel61Val(E0,EF,pj,qj) result(kernel61)
    implicit none
    double precision, dimension(:), intent(in) :: E0, EF
    integer, intent(in) :: pj, qj
    double precision :: kernel61

    kernel61 = E0(pj) * EF(qj)
end function kernel61Val


! kernel54Val.f90
function kernel54Val(E0,E0T,HT,qj,kj) result(theResult)
    implicit none
    double precision, dimension(:), intent(in) :: E0, E0T, HT
    integer, intent(in) :: qj, kj
    double precision :: theResult

    theResult = E0(qj) * 2.0d0 * E0T(kj) * HT(kj)
end function kernel54Val


! kernel53Val.f90
function kernel53Val(E0,E0T,HT,pj,qj) result(theResult)
    implicit none
    double precision, dimension(:), intent(in) :: E0, E0T, HT
    integer, intent(in) :: pj, qj
    double precision :: theResult

    theResult = E0(qj) * E0T(pj) * HT(pj)
end function kernel53Val


! kernel52Val.f90
function kernel52Val(E0,E0T,HDIR,pj,qj,kj) result(theResult)
    implicit none
    double precision, dimension(:), intent(in) :: E0, E0T, HDIR
    integer, intent(in) :: pj, qj, kj
    double precision :: theResult

    theResult = E0(qj) * (E0T(pj) - E0T(kj)) * HDIR(qj)
end function kernel52Val


! kernel51Val.f90
function kernel51Val(E0,E0T,HPOL,pj,qj,kj) result(theResult)
    implicit none
   double precision, dimension(:), intent(in) :: E0, E0T, HPOL
    integer, intent(in) :: pj, qj, kj
   double precision :: theResult

    theResult = E0(qj) * (E0T(pj) - E0T(kj)) * HPOL(qj)
end function kernel51Val


! kernel4Val.f90
function kernel4Val(E,ET,k,p,pj,qj,kj) result(theResult)
    implicit none
   double precision, dimension(:), intent(in) :: E, ET
   double precision, intent(in) :: k, p
    integer, intent(in) :: pj, qj, kj
    double precision :: theResult

    theResult = E(qj) * (k**2. * ET(pj) - p**2. * ET(kj))
end function kernel4Val


! kernel1Val.f90
function kernel1Val(E0,pj,qj,kj) result(theResult)
    implicit none
    double precision, dimension(:), intent(in) :: E0
    integer, intent(in) :: pj, qj, kj
   double precision :: theResult

    theResult = E0(qj) * (E0(pj) - E0(kj))
end function kernel1Val



! interpolateKernel6.f90
subroutine interpolateKernel6(E0,EF,centroidX,centroidY,kVals,pj,qj,Q11,kj, &
                               kernel61,kernel62,kernel63,kernel64,kernel65,kernel66)
    implicit none
   double precision, dimension(:), intent(in) :: E0, EF, kVals
   double precision, intent(in) :: centroidX, centroidY
    integer, dimension(:,:,:,:), intent(in) :: Q11
    integer, intent(in) :: pj, qj, kj
   double precision, intent(out) :: kernel61, kernel62, kernel63, kernel64, kernel65, kernel66

    integer :: pj1, pj2, qj1, qj2
   double precision :: p1, p2, q1, q2
   double precision :: fQ11, fQ12, fQ21, fQ22


    pj1 = Q11(1,pj,qj,kj)

    if (pj1 > 0) then
        p1 = kVals(pj1)
        pj2 = pj1 + 1
        p2 = kVals(pj2)

        qj1 = Q11(2,pj,qj,kj)
        q1 = kVals(qj1)
        qj2 = qj1 + 1
        q2 = kVals(qj2)

        ! kernel61
        fQ11 = E0(pj1) * EF(qj1)
        fQ12 = E0(pj1) * EF(qj2)
        fQ21 = E0(pj2) * EF(qj1)
        fQ22 = E0(pj2) * EF(qj2)
        kernel61 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)

        ! kernel62
        fQ11 = E0(pj1) * EF(kj)
        fQ12 = E0(pj1) * EF(kj)
        fQ21 = E0(pj2) * EF(kj)
        fQ22 = E0(pj2) * EF(kj)
        kernel62 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)

        ! kernel63
        fQ11 = E0(kj) * EF(pj1)
        fQ12 = E0(kj) * EF(pj1)
        fQ21 = E0(kj) * EF(pj2)
        fQ22 = E0(kj) * EF(pj2)
        kernel63 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)

        ! kernel64
        fQ11 = E0(kj) * EF(qj1)
        fQ12 = E0(kj) * EF(qj2)
        fQ21 = E0(kj) * EF(qj1)
        fQ22 = E0(kj) * EF(qj2)
        kernel64 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)

        ! kernel65
        fQ11 = E0(qj1) * EF(pj1)
        fQ12 = E0(qj2) * EF(pj1)
        fQ21 = E0(qj1) * EF(pj2)
        fQ22 = E0(qj2) * EF(pj2)
        kernel65 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)

        ! kernel66
        fQ11 = E0(qj1) * EF(kj)
        fQ12 = E0(qj2) * EF(kj)
        fQ21 = E0(qj1) * EF(kj)
        fQ22 = E0(qj2) * EF(kj)
        kernel66 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)
    else
        kernel61 = 0.0d0
        kernel62 = 0.0d0
        kernel63 = 0.0d0
        kernel64 = 0.0d0
        kernel65 = 0.0d0
        kernel66 = 0.0d0
    end if
end subroutine interpolateKernel6


! interpolateKernel3.f90
subroutine interpolateKernel3(E0,centroidX,centroidY,kVals,pj,qj,Q11,kj,HPOL,HDIR, &
                               kernel31,kernel32,kernel33,kernel34,kernel35,kernel36,kernel37)
    implicit none
   double precision, dimension(:), intent(in) :: E0, kVals, HPOL, HDIR
double precision, intent(in) :: centroidX, centroidY
    integer, dimension(:,:,:,:), intent(in) :: Q11
    integer, intent(in) :: pj, qj, kj
  double precision, intent(out) :: kernel31, kernel32, kernel33, kernel34, kernel35, kernel36, kernel37

    integer :: pj1, pj2, qj1, qj2
   double precision :: p1, p2, q1, q2
  double precision :: fQ11, fQ12, fQ21, fQ22


    pj1 = Q11(1,pj,qj,kj)

    if (pj1 > 0) then
        p1 = kVals(pj1)
        pj2 = pj1 + 1
        p2 = kVals(pj2)

        qj1 = Q11(2,pj,qj,kj)
        q1 = kVals(qj1)
        qj2 = qj1 + 1
        q2 = kVals(qj2)

        ! kernel31
        fQ11 = E0(qj1) * E0(pj1) * HPOL(pj1)
        fQ12 = E0(qj2) * E0(pj1) * HPOL(pj1)
        fQ21 = E0(qj1) * E0(pj2) * HPOL(pj2)
        fQ22 = E0(qj2) * E0(pj2) * HPOL(pj2)
        kernel31 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)

        ! kernel32
        fQ11 = E0(qj1) * E0(kj) * HPOL(kj)
        fQ12 = E0(qj2) * E0(kj) * HPOL(kj)
        fQ21 = E0(qj1) * E0(kj) * HPOL(kj)
        fQ22 = E0(qj2) * E0(kj) * HPOL(kj)
        kernel32 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)

        ! kernel33
        fQ11 = E0(qj1) * (E0(pj1) - E0(kj)) * HPOL(qj1)
        fQ12 = E0(qj2) * (E0(pj1) - E0(kj)) * HPOL(qj2)
        fQ21 = E0(qj1) * (E0(pj2) - E0(kj)) * HPOL(qj1)
        fQ22 = E0(qj2) * (E0(pj2) - E0(kj)) * HPOL(qj2)
        kernel33 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)

        ! kernel34
        fQ11 = E0(qj1) * E0(pj1) * HPOL(pj1)
        fQ12 = E0(qj2) * E0(pj1) * HPOL(pj1)
        fQ21 = E0(qj1) * E0(pj2) * HPOL(pj2)
        fQ22 = E0(qj2) * E0(pj2) * HPOL(pj2)
        kernel34 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)

        ! kernel35
        fQ11 = E0(qj1) * E0(kj) * HPOL(qj1)
        fQ12 = E0(qj2) * E0(kj) * HPOL(qj2)
        fQ21 = E0(qj1) * E0(kj) * HPOL(qj1)
        fQ22 = E0(qj2) * E0(kj) * HPOL(qj2)
        kernel35 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)

        ! kernel36
        fQ11 = E0(qj1) * (E0(pj1) - E0(kj)) * HDIR(qj1)
        fQ12 = E0(qj2) * (E0(pj1) - E0(kj)) * HDIR(qj2)
        fQ21 = E0(qj1) * (E0(pj2) - E0(kj)) * HDIR(qj1)
        fQ22 = E0(qj2) * (E0(pj2) - E0(kj)) * HDIR(qj2)
        kernel36 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)

        ! kernel37
        fQ11 = E0(qj1) * E0(pj1) * HDIR(pj1)
        fQ12 = E0(qj2) * E0(pj1) * HDIR(pj1)
        fQ21 = E0(qj1) * E0(pj2) * HDIR(pj2)
        fQ22 = E0(qj2) * E0(pj2) * HDIR(pj2)
        kernel37 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)
    else
        kernel31 = 0.0d0
        kernel32 = 0.0d0
        kernel33 = 0.0d0
        kernel34 = 0.0d0
        kernel35 = 0.0d0
        kernel36 = 0.0d0
        kernel37 = 0.0d0
    end if
end subroutine interpolateKernel3


! interpolateKernel2.f90
subroutine interpolateKernel2(E0,centroidX,centroidY,kVals,pj,qj,Q11,kj,HPOL,HDIR, &
                               kernel21,kernel22,kernel23,kernel24,kernel25)
    implicit none
   double precision, dimension(:), intent(in) :: E0, kVals, HPOL, HDIR
   double precision, intent(in) :: centroidX, centroidY
    integer, dimension(:,:,:,:), intent(in) :: Q11
    integer, intent(in) :: pj, qj, kj
 double precision, intent(out) :: kernel21, kernel22, kernel23, kernel24, kernel25

    integer :: pj1, pj2, qj1, qj2
   double precision :: p1, p2, q1, q2
   double precision :: fQ11, fQ12, fQ21, fQ22


    pj1 = Q11(1,pj,qj,kj)

    if (pj1 > 0) then
        p1 = kVals(pj1)
        pj2 = pj1 + 1
        p2 = kVals(pj2)

        qj1 = Q11(2,pj,qj,kj)
        q1 = kVals(qj1)
        qj2 = qj1 + 1
        q2 = kVals(qj2)

        ! kernel21
        fQ11 = E0(qj1) * (E0(pj1) - E0(kj)) * HPOL(qj1)
        fQ12 = E0(qj2) * (E0(pj1) - E0(kj)) * HPOL(qj2)
        fQ21 = E0(qj1) * (E0(pj2) - E0(kj)) * HPOL(qj1)
        fQ22 = E0(qj2) * (E0(pj2) - E0(kj)) * HPOL(qj2)
        kernel21 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)

        ! kernel22
        fQ11 = E0(qj1) * E0(pj1) * HPOL(pj1)
        fQ12 = E0(qj2) * E0(pj1) * HPOL(pj1)
        fQ21 = E0(qj1) * E0(pj2) * HPOL(pj2)
        fQ22 = E0(qj2) * E0(pj2) * HPOL(pj2)
        kernel22 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)

        ! kernel23
        fQ11 = E0(qj1) * (E0(pj1) - E0(kj)) * HDIR(qj1)
        fQ12 = E0(qj2) * (E0(pj1) - E0(kj)) * HDIR(qj2)
        fQ21 = E0(qj1) * (E0(pj2) - E0(kj)) * HDIR(qj1)
        fQ22 = E0(qj2) * (E0(pj2) - E0(kj)) * HDIR(qj2)
        kernel23 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)

        ! kernel24
        fQ11 = E0(qj1) * E0(pj1) * HDIR(pj1)
        fQ12 = E0(qj2) * E0(pj1) * HDIR(pj1)
        fQ21 = E0(qj1) * E0(pj2) * HDIR(pj2)
        fQ22 = E0(qj2) * E0(pj2) * HDIR(pj2)
        kernel24 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)

        ! kernel25
        fQ11 = E0(qj1) * E0(kj) * HDIR(kj)
        fQ12 = E0(qj2) * E0(kj) * HDIR(kj)
        fQ21 = E0(qj1) * E0(kj) * HDIR(kj)
        fQ22 = E0(qj2) * E0(kj) * HDIR(kj)
        kernel25 = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)
    else
        kernel21 = 0.0d0
        kernel22 = 0.0d0
        kernel23 = 0.0d0
        kernel24 = 0.0d0
        kernel25 = 0.0d0
    end if
end subroutine interpolateKernel2


! interpolateKernel.f90
function interpolateKernel(E,ET,centroidX,centroidY,k,kVals,pj,qj,Q11,eqnNumber,kj,H) result(kernel)
    implicit none
   double precision, dimension(:), intent(in) :: E, ET, kVals, H
   double precision, intent(in) :: centroidX, centroidY, k
    integer, dimension(:,:,:,:), intent(in) :: Q11
    integer, intent(in) :: pj, qj, eqnNumber, kj
   double precision :: kernel

    integer :: pj1, pj2, qj1, qj2
   double precision :: p1, p2, q1, q2
   double precision :: fQ11, fQ12, fQ21, fQ22



    pj1 = Q11(1,pj,qj,kj)

    if (pj1 > 0) then
        p1 = kVals(pj1)
        pj2 = pj1 + 1
        p2 = kVals(pj2)

        qj1 = Q11(2,pj,qj,kj)
        q1 = kVals(qj1)
        qj2 = qj1 + 1
        q2 = kVals(qj2)

        if (eqnNumber == 1) then
            fQ11 = kernel1Val(E,pj1,qj1,kj)
            fQ12 = kernel1Val(E,pj1,qj2,kj)
            fQ21 = kernel1Val(E,pj2,qj1,kj)
            fQ22 = kernel1Val(E,pj2,qj2,kj)
        elseif (eqnNumber == 4) then
            fQ11 = kernel4Val(E,ET,k,p1,pj1,qj1,kj)
            fQ12 = kernel4Val(E,ET,k,p1,pj1,qj2,kj)
            fQ21 = kernel4Val(E,ET,k,p2,pj2,qj1,kj)
            fQ22 = kernel4Val(E,ET,k,p2,pj2,qj2,kj)
        elseif (eqnNumber == 51) then
            fQ11 = kernel51Val(E,ET,H,pj1,qj1,kj)
            fQ12 = kernel51Val(E,ET,H,pj1,qj2,kj)
            fQ21 = kernel51Val(E,ET,H,pj2,qj1,kj)
            fQ22 = kernel51Val(E,ET,H,pj2,qj2,kj)
        elseif (eqnNumber == 52) then
            fQ11 = kernel52Val(E,ET,H,pj1,qj1,kj)
            fQ12 = kernel52Val(E,ET,H,pj1,qj2,kj)
            fQ21 = kernel52Val(E,ET,H,pj2,qj1,kj)
            fQ22 = kernel52Val(E,ET,H,pj2,qj2,kj)
        elseif (eqnNumber == 53) then
            fQ11 = kernel53Val(E,ET,H,pj1,qj1)
            fQ12 = kernel53Val(E,ET,H,pj1,qj2)
            fQ21 = kernel53Val(E,ET,H,pj2,qj1)
            fQ22 = kernel53Val(E,ET,H,pj2,qj2)
        elseif (eqnNumber == 54) then
            fQ11 = kernel54Val(E,ET,H,qj1,kj)
            fQ12 = kernel54Val(E,ET,H,qj2,kj)
            fQ21 = kernel54Val(E,ET,H,qj1,kj)
            fQ22 = kernel54Val(E,ET,H,qj2,kj)
        end if

        kernel = bilinear(centroidX,centroidY,p1,q1,p2,q2,fQ12,fQ22,fQ21,fQ11)
    else
        kernel = 0.0d0
    end if
end function interpolateKernel

end module kernel_interp
