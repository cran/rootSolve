c------------------------------------------ --------------------------*
c SPARSITY of 1-D PDE problems
c--------------------------------------------------------------------*

      SUBROUTINE sparse1d (Ntot, Nspec, dim, cyclic, nnz, ian, jan)

c--------------------------------------------------------------------*
c Determines the sparsity structure of the Jacobian in 1-D reaction- *
c transport problems.                                                *
c                                                                    *
c two arrays describe the sparsity structure of the jacobian:        *
c                                                                    *
C IAN, of size NNtot + 1,                                            *
C JAN, of size NNZ. (to be determined by the user).                  *
C                                                                    *
C JAN contains the row indices of the nonzero locations of           *
C the jacobian, reading in columnwise order, and                     *
C IAN contains the starting locations in JAN of the descriptions of  *
C columns 1,...,NEQ, with IAN(1) = 1 and IAN(NEQ+1) = NNZ + 1.       *
c                                                                    *
C Thus for each j = 1,...,NEQ, the row indices i of the              *
C nonzero locations in column j are:                                 *
C                        i = JAN(k),  IAN(j) .le. k .lt. IAN(j+1).   *
C                                                                    *
c--------------------------------------------------------------------*

c total number of state variables,number of different species
       INTEGER Ntot, Nspec, dim, cyclic
c maximal number of indices, sparsity arrays
       INTEGER nnz, ian(*), jan(*)
c
       INTEGER N, I, J, ij, K, L, isp
 	     character *80 msg

c check input
       IF (INT(Ntot/Nspec)*Nspec .NE. Ntot) THEN
         write(msg,*)
     &("cannot generate sparse jacobian - N and nspec not compatible")
         call rexit(msg)
       ENDIF

c number of boxes
       N = Ntot/Nspec

       ij     = 1
       ian(1) = 1

       DO i = 1, Nspec
         isp = (i-1)*N
         DO j = 1, N
           K = (i-1)*N+j

c interactions with current, upstream and downstream boxes
          
           jan(ij) = K
           ij      = ij +1
           IF (J<N) THEN
             jan(ij) = K+1
             ij      = ij +1
           ELSE IF (cyclic == 1) THEN
              jan(ij) = isp + 1
              ij      = ij +1
           ENDIF
           IF (J >1) THEN
             jan(ij) = K-1
             ij      = ij +1
           ELSE IF (cyclic == 1) THEN
              jan(ij) = isp + N
              ij      = ij +1
           ENDIF
c interactions with other species in the same box
            DO L = 1, Nspec
              IF (L == i) cycle
              jan(ij) = (L-1)*N+j
              ij = ij +1            
            ENDDO
       
            ian(K+1) = ij

         ENDDO
       ENDDO
       nnz = ij -1

      END SUBROUTINE sparse1d

c--------------------------------------------------------------------*
c SPARSITY of 2-D reaction-transport problems
c--------------------------------------------------------------------*

      SUBROUTINE sparse2d (Ntot, Nspec, dimens, cyclic, nnz, ian, jan)

c--------------------------------------------------------------------*
c Determines the sparsity structure of the Jacobian in 2-D reaction- *
c transport problems.                                                *
c                                                                    *
C IAN  and JAN: see comments in sparse1D                             *
C                                                                    *
c--------------------------------------------------------------------*

c total number of state variables, number of different species
c dimensions of the problem and whether cyclic boundaries or not
       INTEGER Ntot, Nspec, dimens(*), cyclic(*)

c maximal number of indices, sparsity arrays
       INTEGER nnz, ian(*), jan(*)
c
       INTEGER N, I, J, ij, K, L, M , isp
 	     character *80 msg

c check input
       IF (INT(Ntot/Nspec)*Nspec .NE. Ntot) THEN
         write(msg,*)
     &("cannot generate sparse jacobian - N and nspec not compatible")
         call rexit(msg)
       ENDIF

c number of boxes
       N = dimens(1)*dimens(2)

       ij     = 1
       ian(1) = 1

       DO i = 1, Nspec
         isp = (i-1)*N
         DO j = 1, dimens(1)
           DO k = 1, dimens(2)
             M = isp +(j-1)*dimens(2)+k

c interactions with current, upstream and downstream boxes
            jan(ij) = M
            ij      = ij +1

            IF (k<dimens(2)) THEN
              jan(ij) = M+1
              ij      = ij +1
            ELSE IF (cyclic(2) == 1) THEN
              jan(ij) = isp + (j-1)*dimens(2) +1
              ij      = ij +1
            ENDIF

            IF (j<dimens(1)) THEN
              jan(ij) = M+dimens(2)
              ij      = ij +1
            ELSE IF (cyclic(1) == 1) THEN
              jan(ij) = isp + K
              ij      = ij +1
            ENDIF

            IF (j >1) THEN
              jan(ij) = M-dimens(2)
              ij      = ij +1
            ELSE IF (cyclic(1) == 1) THEN
              jan(ij) = isp + (dimens(1)-1)*dimens(2)+ K
              ij      = ij +1
            ENDIF

            IF (k >1) THEN
              jan(ij) = M-1
              ij      = ij +1
            ELSE IF (cyclic(2) == 1) THEN
              jan(ij) = isp + j*dimens(2)
              ij      = ij +1
            ENDIF

c interactions with other species in the same box
            DO L = 1, Nspec
              IF (L == i) cycle
              jan(ij) = (L-1)*N+(j-1)*dimens(2)+k
              ij = ij +1            
            ENDDO
       
            ian(M+1) = ij
           ENDDO
         ENDDO
       ENDDO
       nnz = ij -1
c
      END SUBROUTINE sparse2d

c--------------------------------------------------------------------*
c SPARSITY of 3-D PDE problems
c--------------------------------------------------------------------*

      SUBROUTINE sparse3d (Ntot, Nspec, dimens, cyclic, nnz, ian, jan)

c--------------------------------------------------------------------*
c Determines the sparsity structure of the Jacobian in 3-D PDE       *
c problems.                                                          *
c                                                                    *
C IAN  and JAN: see comments in sparse1D                             *
C                                                                    *
c--------------------------------------------------------------------*

c total number of state variables, number of different species
c dimensions of the problem and whether cyclic boundaries or not
       INTEGER Ntot, Nspec, dimens(*), cyclic(*)

c maximal number of indices, sparsity arrays
       INTEGER nnz, ian(*), jan(*)
c
       INTEGER N, I, J, ij, K, L, M , isp
 	     character *80 msg

c check input
       IF (INT(Ntot/Nspec)*Nspec .NE. Ntot) THEN
         write(msg,*)
     &("cannot generate sparse jacobian - N and nspec not compatible")
         call rexit(msg)
       ENDIF

c number of boxes
       N = dimens(1)*dimens(2)

       ij     = 1
       ian(1) = 1

       DO i = 1, Nspec
         isp = (i-1)*N
         DO j = 1, dimens(1)
           DO k = 1, dimens(2)
             DO ll = 1, dimens(3)
               M = isp +(j-1)*dimens(2)*dimens(3)+(k-1)*dimens(3)+ll

c interactions with current, upstream and downstream boxes
              jan(ij) = M
              ij      = ij +1

              IF (ll<dimens(3)) THEN
                jan(ij) = M+1
                ij      = ij +1
              ENDIF
              IF (ll >1) THEN
                jan(ij) = M-1
                ij      = ij +1
              ENDIF
 
              IF (k<dimens(2)) THEN
                jan(ij) = M+dimens(3)
                ij      = ij +1
              ENDIF
            
              IF (k >1) THEN
                jan(ij) = M-dimens(3)
                ij      = ij +1
              ENDIF

              IF (j<dimens(1)) THEN
                jan(ij) = M+dimens(2)*dimens(3)
                ij      = ij +1
              ENDIF

              IF (j >1) THEN
                jan(ij) = M-dimens(2)*dimens(3)
                ij      = ij +1
              ENDIF


c interactions with other species in the same box
              DO L = 1, Nspec
                IF (L == i) cycle
                jan(ij)=(L-1)*N+(j-1)*dimens(2)*dimens(3)+                       &
     &                     (k-1)*dimens(3)+ll
                ij = ij +1
              ENDDO

              ian(M+1) = ij
            ENDDO
	    ENDDO
        ENDDO
      ENDDO
      nnz = ij -1
c
      END SUBROUTINE sparse3d