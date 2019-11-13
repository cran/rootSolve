     
c######################################################################
c
c STEADY-STATE SOLVER - sparse jacobian
c
c FINDS THE ROOT OF A SET OF NONLINEAR EQUATIONS               
c implementation: karline Soetaert, NIOZ-yerseke, the Netherlands
c
c uses linear algebra routines from the sparsekit package:
c
c######################################################################

c                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
c                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
c                !       SOLVING STEADY-STATE         !
c                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
c                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!

c**********************************************************************
     
       SUBROUTINE dsparsekit(xmodel,N,nnz,nsp,time,Svar,dSvar,beta,x,          &
     &                   a,ewt,ian,jan,igp,jgp,maxg,jlu,ju,iwork,iperm,        &
     &                   maxiter,TolChange,atol,rtol,itol,PositivityInt,       &
     &                   Pos,ipos,SteadyStateReachedInt,Precis,niter,          &
     &                   dims,out,nout,Type,droptol,permtol,imethod,           &
     &                   lfill,lenplumx,plu,rwork,pres)

c------------------------------------------------------------------------------*
c Solves a system of nonlinear equations using the Newton-Raphson method       *
c assumes a sparse Jacobian                                                    *
c------------------------------------------------------------------------------*  
      IMPLICIT NONE
  
c number of equations, maximal number of nonzero elements in jacobian 
c length of work arrays, max and actual number of independent groups
      INTEGER    N, nnz, nsp, maxg, NGP     

c actual number of nonzeros, max and actual iterations
      INTEGER  nonzero, maxiter, niter, dims(*)

c indices to nonzero elements and to groups of independent state variables 
      INTEGER ian(*), jan(*), igp(*),jgp(*), pres(*)
c locals....      
      integer iao(N+1), jao(nnz)
      double precision ao(nnz)
      
c state variables, 
      DOUBLE PRECISION Svar(*) 

c Beta : the negative of the rate of change,
      DOUBLE PRECISION BETA(*), dSvar(*), x(*)

c transpose of jacobian  KS: check if transpose or Jacobian
      DOUBLE PRECISION  a(*)

c false if failed - true if variables must be positive 
c positivity either enforced at once (positivity=TRUE) or as a vector of elements
      LOGICAL SteadyStateReached, positivity
      INTEGER SteadyStateReachedInt, PositivityInt
      INTEGER  Ipos, Pos(Ipos)

c tolerances, precision
      INTEGER          itol, Type
      DOUBLE PRECISION rtol(*), atol(*),tolChange
      DOUBLE PRECISION ewt(*), precis(maxIter),maxewt,RelativeChange
         
c working arrays for sparsekit solvers
      INTEGER jlu(*)   ! length lenplumx
      INTEGER ju(*)    ! length N
      INTEGER iwork(*) ! length 2*N
      INTEGER iperm(*) ! length 2*N

      DOUBLE PRECISION plu(*)    ! length lenplumx
      DOUBLE PRECISION rwork(*)  ! length N
c input for sparsekit solvers
      DOUBLE PRECISION droptol   ! drop tolerance for ilut, ilutp
      DOUBLE PRECISION permtol   ! tolerance ratio for column pivoting in ilutp
      INTEGER imethod      ! 1= use ilut, 2= use ilutp
      INTEGER lfill        ! the level of fill-in
      INTEGER lenplumx     ! = nnz+lenplufac*N, length of workarray plu

c model and jacobian function
      EXTERNAL xmodel
      DOUBLE PRECISION out(*)      
      INTEGER          nout(*) 
      DOUBLE PRECISION time
c
      INTEGER i, j, k, ierr
c-------------------------------------------------------------------------------
      SteadyStateReached = .FALSE.
      SteadyStateReachedInt = 0
      Positivity = .FALSE.
      IF (positivityInt > 0.1) Positivity = .TRUE.

      CALL errSET (N, ITOL, RTOL, ATOL, SVAR, EWT)

c determine sparse structure: if Type == 2, 3, 4 :
c a 1-D or 2-D or 3-D PDE model;
c in this case the number of components, dimensions and cyclic bnd are in dims

      CALL xSparseStruct(N, nnz, ian, jan, igp, jgp, maxg, ngp,                &
     &    Svar, ewt, dSvar, beta, xmodel, time, out, nout, nonzero,            &
     &    Type, dims, pres)

c initial guess for x
      DO j=1, N
        x (j) = 0.D0
      ENDDO

c Iterations
      DO I = 1, maxiter 
        niter = I

c Create sparse jacobian - KARLINE CHECK IF NOT TRANPSOSED
        CALL xSparseJacob (N, nnz, ian, jan, igp, jgp, ngp,                    &
     &     Svar, ewt, dSvar, beta, xmodel, time, out, nout, a)
c transpose...
        CALL csrcsc (N,1,1,a,jan,ian,ao,jao,iao)
        do k = 1, N +1
          ian(k) = iao(k)
        enddo
        do k = 1, nnz
          a(k) = ao(k)
          jan(k) = jao(k)
        enddo
        
c Check convergence 
        precis(I) = 0.d0
        maxewt    = 0.d0
        DO k = 1, N
          precis(I) =precis(I)+ abs(beta(k))
          maxewt = MAX(maxewt, abs(BETA(k)/ewt(k)))
        ENDDO
        precis(i) = precis(i)/N
        IF (maxewt .LE. 1) THEN
          SteadyStateReached = .TRUE.
          SteadyStateReachedInt = 1
          EXIT
        ENDIF

c==========================================================================
c LU factorisation, using either ILUT (imethod=1) or ILUTP (imethod=2)
c==========================================================================

      IF (imethod .EQ. 1) THEN
C Use incomplete factorization routine ILUT from SparsKit.
         CALL ILUT (N, a, jan, ian, lfill, droptol, plu, jlu,
     &              ju,lenplumx,rwork,iwork,IERR)
         IF (IERR .NE. 0) THEN
            call rwarn ("error return from ILUT")
         ENDIF
      ELSE
C Use incomplete factorization routine ILUTP from SparsKit.
         CALL ILUTP (N, a, jan, ian, lfill, droptol, permtol, N,
     &               plu, jlu, ju, lenplumx, rwork, iwork, iperm, IERR)
         IF (IERR .NE. 0) THEN
            call rwarn ("error return from ILUTP")
         ENDIF
      ENDIF

C an error occurred
      IF (IERR .NE. 0) CALL warnflagkit(ierr)

C Solve the linear system P*x=c, using elements of P loaded into
C arrays a, jan and ian.
       CALL LUSOL (N, beta, x, plu, jlu, ju)

c Test convergence + new value of state variables
        RelativeChange=0.d0

        DO k=1, N
          RelativeChange  = MAX(RelativeChange,ABS(x(k)))
          Svar(k)         = Svar(k)+x(k)
          IF (Positivity) Svar(k)=MAX(0.D0,Svar(k))
        ENDDO
        IF (.not. positivity .and. ipos .GT. 1) THEN
           DO K = 1, ipos
             Svar(Pos(K)) = MAX(0.D0,Svar(Pos(K)))
           ENDDO
        ENDIF

        IF (RelativeChange<=TolChange)THEN
c last precision reached
          IF (i .LT.  maxiter) THEN
            precis(i+1) = 0.d0
            DO j = 1, N
              beta(j) = 0.D0
            ENDDO
            CALL XMODEL(N,time,Svar,beta,out,nout)

            DO j=1, N
              precis(i+1) = precis(i+1)+abs(beta(j))
            ENDDO
            precis(i+1) = precis(i+1)
            niter = I+1
          ENDIF
          SteadyStateReached = .TRUE.
          STeadyStateReachedInt = 1
          EXIT
        ENDIF

        CALL errSET (N, ITOL, RTOL, ATOL, SVAR, EWT)


      ENDDO

      SteadyStateReachedInt = 0
      IF (SteadyStateReached) SteadyStateReachedInt = 1
      
c put some values in dims - 
c       dims(3) = nsp - esp

      END SUBROUTINE dsparsekit

c                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
c                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
c                !             FUNCTIONS              !
c                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
c                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!

c**********************************************************************
c       WRITE ERROR/WARNINGS OF SPARSEKIT SOLVER                      *
c**********************************************************************


      SUBROUTINE warnflagkit(ierr)
      INTEGER Ierr
c-------------------------------------------------------------------------------

      IF (IERR .GT. 0) THEN 
        call intpr("zero pivot encountered at step nr ",35, IERR, 1)
C        write(msg,'(A35,I10)')"zero pivot encountered at step nr ",IERR
C        call rwarn(msg)
      ELSE IF (IERR .EQ. -1) THEN
        call rwarn("input matrix may be wrong; elimination process ")
        call rwarn("generated a row in L or U ")
        call rwarn("with length exceeding N")
        call rexit("stopped")
      ELSE IF (IERR .EQ. -2) THEN
        call rwarn("matrix L overflows")
        call rwarn("increase value of lenplufac or decrease value of")
        call rwarn("lfill if lenplufac cannot be increased")
        call rexit("stopped")
      ELSE IF (IERR .EQ. -3) THEN
        call rwarn("matrix U overflows")
        call rwarn("increase value of lenplufac or decrease value")
        call rwarn("lfill if lenplufac cannot be increased")
        call rexit("stopped")
      ELSE IF (IERR .EQ. -4) THEN
        call rexit("illegal value for lfill")
      ELSE IF (IERR .EQ. -5) THEN
        call rexit("zero row encountered")
      ENDIF
     
      END SUBROUTINE warnflagkit 
