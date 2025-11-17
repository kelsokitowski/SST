module getForcing
  implicit none

contains


subroutine getGaussExtForcing(kVals, E0, eta, nu, ExtForcing)
    !! ------------------------------------------------------------------
    !! Purpose:
    !!   Compute Gaussian external forcing spectrum consistent with
    !!   ε = ν³ / η⁴, centered near the spectral peak of E0.
    !!
    !! Inputs:
    !!   kVals(:)   - wavenumber array (monotonic increasing)
    !!   E0(:)      - energy spectrum array (used to locate forcing peak)
    !!   eta, nu    - Kolmogorov length scale and viscosity
    !!
    !! Output:
    !!   ExtForcing(:) - forcing amplitude: Cval*k*exp(-((k-kF)/d)²)
    !!
    !! Notes:
    !!   - Uses Newton iteration for Cval (analogous to MATLAB fzero)
    !!   - Integral ∫ kernel(x) dx verified ≈ ε
    !! ------------------------------------------------------------------
    implicit none
    integer, parameter :: dp = 8
    real(dp), intent(in)  :: kVals(:), E0(:), eta, nu
    real(dp), intent(out) :: ExtForcing(size(kVals))

    integer :: n, kLIndex, shft, iter, maxIter
    real(dp) :: epzilon, kF, d, Cval, Cguess, err, dC, tol
    real(dp) :: erf_arg, int_est, pi
    real(dp), allocatable :: kernel(:)
    integer :: i

    n  = size(kVals)
    pi = acos(-1.0_dp)

    !--- Locate maximum of E0 (forcing near energy peak)
    kLIndex = maxloc(E0, dim=1)
    if (kLIndex <= 7) kLIndex = 8  ! safety if near boundary

    epzilon = nu**3.0_dp / (eta**4.0_dp)
    shft = 7
    kF = kVals(kLIndex - shft)
    Cguess = erf(kVals(kLIndex - shft))
    d = kVals(min(13, n)) - kF   ! use spacing around index 13 (MATLAB analogue)

    !--- Solve for Cval such that integral = ε
    tol = 1.0e-10_dp
    maxIter = 50
    Cval = Cguess + 100.0_dp

    do iter = 1, maxIter
        erf_arg = (kVals(1) - kF) / d
        err = Cval * d * sqrt(pi) / 2.0_dp * (1.0_dp - erf(erf_arg)) - epzilon
        if (abs(err) < tol) exit
        dC = -err / (d * sqrt(pi) / 2.0_dp * (1.0_dp - erf(erf_arg)))
        Cval = Cval + dC
    end do

    if (iter >= maxIter) then
        print *, 'Warning [getGaussExtForcing]: Newton solver did not fully converge.'
    end if

    !--- Construct Gaussian forcing
    ExtForcing = Cval * kVals * exp( -((kVals - kF)/d)**2. )

    !--- Optional verification: ∫ kernel(x) dx ≈ ε
    allocate(kernel(n))
    do i = 1, n
        kernel(i) = Cval * exp( -((kVals(i) - kF)/d)**2. )
    end do
    int_est = trapz(kVals, kernel)
    if (abs(epzilon - int_est) / epzilon > 1.0e-3_dp) then
        print *, 'Warning [getGaussExtForcing]: integral mismatch =', (epzilon - int_est)/epzilon
    end if
    deallocate(kernel)

contains
    pure function trapz(x, y) result(intg)
        !! Simple trapezoidal integrator for verification
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: intg
        integer :: i, n
        n = size(x)
        intg = 0.0_dp
        do i = 1, n-1
            intg = intg + 0.5_dp * (x(i+1)-x(i)) * (y(i+1)+y(i))
        end do
    end function trapz
end subroutine getGaussExtForcing



end module getForcing
