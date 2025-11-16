module integrationModules
  implicit none


  contains

    !
    !A function to calculate trapezoidal rule for integral with non constant dx spacing
    !
    double precision function trapezoidalIntegration(xvals,Fvals)
      double precision, intent(in) :: xvals(:), Fvals(:)

     double precision :: h, theResult
      integer :: j
     theResult = 0.


     do j = 2,size(xvals,dim = 1)
        h = (xvals(j)-xvals(j-1))/2.;
        theResult = theResult + h*(Fvals(j)+Fvals(j-1));
     end do
     trapezoidalIntegration = theResult;

   end function trapezoidalIntegration

!!$subroutine getVolumes(kVals,kLength,weight)
!!$  double precision, intent(in)::kVals(:)
!!$  integer, intent(in) :: kLength
!!$  double precision, intent(out):: weight(kLength,kLength,kLength)
!!$  double precision :: dkM, dkP, dpM, dpP, dqM, dqP, dV
!!$  integer :: i,j,k,m,n,pj,qj,boundaryFlag,kj
!!$  integer, parameter :: nPointsInInterval = 100;
!!$  double precision :: kLocal(100), kM,kP, localMax, localVolume,p,pLocal(100),pM,pP,q,qLocal(100),qM,qP
!!$  double precision :: diff, intervalMin, intervalMax,linspaceResult(nPointsInInterval)
!!$
!!$n = kLength;
!!$weight = 0.0;
!!$!for kj = 1:n
!!$do kj = 1,n
!!$
!!$
!!$    do pj = 1,n
!!$
!!$
!!$        do qj = 1,n
!!$
!!$            p = kVals(pj);
!!$            q = kVals(qj);
!!$            k = kVals(kj);
!!$            if triadCondition(q,p,k) then
!!$                if (kj >1) then
!!$
!!$                    dkM = (k-kVals(kj-1))/2.0;
!!$
!!$                    kM = k-dkM;
!!$
!!$                else
!!$                    kM = k;
!!$                end if
!!$                if (kj <n) then
!!$                    dkP = (kVals(kj+1)-k)/2.0;
!!$                    kP = k+dkP;
!!$                else
!!$                    kP = k;
!!$                end if
!!$                if (pj>1) then
!!$                dpM = (p-kVals(pj-1))/2.0;
!!$                pM = p-dpM;
!!$                else
!!$                    pM = p;
!!$                end if
!!$                if (pj<n) then
!!$                dpP = (kVals(pj+1)-p)/2.0;
!!$                pP = p+dpP;
!!$                else
!!$                    pP = p;
!!$                end if
!!$                if (qj>1) then
!!$
!!$                dqM = (q-kVals(qj-1))/2.0;
!!$                qM = q-dqM;
!!$                else
!!$                    qM = q;
!!$                end if
!!$                if (qj<n) then
!!$                 dqP = (kVals(qj+1)-q)/2.0;
!!$                qP = q+dqP;
!!$                else
!!$                    qP = q;
!!$                end if
!!$                ! Calculate volume or perform operations with kM, kP, pM, pP, qM, qP here
!!$                !define local cube mesh
!!$                intervalMin = kM;
!!$                intervalMax = kP;
!!$                diff =(intervalMax- intervalMin)/(nPointsInInterval-1);
!!$                linspaceResult(1) = intervalMin;
!!$                do i=2,nPointsInInterval
!!$                   linspaceResult(i) = test1(1) + (i-1)*diff;
!!$                end do
!!$                kLocal = linspaceResult;
!!$
!!$                 intervalMin = pM;
!!$                intervalMax = pP;
!!$                diff =(intervalMax- intervalMin)/(nPointsInInterval-1);
!!$                linspaceResult(1) = intervalMin;
!!$                do i=2,nPointsInInterval
!!$                   linspaceResult(i) = test1(1) + (i-1)*diff;
!!$                end do
!!$                pLocal = linspaceResult;
!!$
!!$                 intervalMin =qM;
!!$                intervalMax = qP;
!!$                diff =(intervalMax- intervalMin)/(nPointsInInterval-1);
!!$                linspaceResult(1) = intervalMin;
!!$                do i=2,nPointsInInterval
!!$                   linspaceResult(i) = test1(1) + (i-1)*diff;
!!$                end do
!!$                qLocal = linspaceResult;
!!$
!!$                localLength = nPointsInInterval;
!!$                dV = (kLocal(2)-kLocal(1))*(pLocal(2)-pLocal(1))*(qLocal(2)-qLocal(1));
!!$                LocalMax = (kP-kM)*(qP-qM)*(pP-pM);
!!$                LocalVolume = 0.0;
!!$                boundaryFlag = 0;
!!$                do i = 1,localLength-1
!!$                    do j = 1,localLength-1
!!$                        do m = 1,localLength-1
!!$                            if triadCondition(qLocal(i),pLocal(j),kLocal(m)) then
!!$                                LocalVolume = LocalVolume + dV;
!!$                            else
!!$                                boundaryFlag = 1;
!!$                            end if
!!$
!!$                        end do
!!$                    end do
!!$                end do
!!$                if (boundaryFlag == 1) then
!!$                    weight(pj,qj) = LocalVolume/LocalMax;
!!$                else
!!$                    weight(pj,qj) = 1.0;
!!$                end if
!!$
!!$
!!$
!!$
!!$            end if !triadCondition
!!$        end do !do p
!!$    end do !do q
!!$
!!$
!!$end do !do kj
!!$
!!$end subroutine getVolumes









   !below works just like Matlab cumTrapz

  subroutine cumTrapz(xvals,Fvals,cumTrapzVals)
    double precision, intent(in) :: xvals(:), Fvals(:)
    double precision :: h, theResult
     double precision, intent(out) :: cumTrapzVals(size(xVals))
      integer :: j
      cumTrapzVals = 0.;
     theResult = 0.
     h = xvals(1)
     do j = 2,size(xvals,dim = 1)
        h = (xvals(j)-xvals(j-1))/2.;
        theResult = theResult + h*(Fvals(j)+Fvals(j-1));
        cumTrapzVals(j-1) = theResult;
     end do
   end subroutine cumTrapz




integer function triadCondition(p,q,k)
  double precision, intent(in) :: p,q,k
  integer :: flag
  flag = 0;
if (q >= (k-p) ) then
    if ( q <= (k+p) ) then
        if (q >= (p-k) ) then
            flag = 1;
        end if
    end if
end if
triadCondition = flag;
end function triadCondition


end module integrationModules
