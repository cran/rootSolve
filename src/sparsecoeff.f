c--------------------------------------------------------------------*
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

c check input
       IF (INT(Ntot/Nspec)*Nspec .NE. Ntot) THEN

         call rexit
     &("cannot generate sparse jacobian - N and nspec not compatible")
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
C A(I,J) in vector: at position (I-1)*dim(2) + J                     *
C                                                                    *
c--------------------------------------------------------------------*

c total number of state variables, number of different species
c dimensions of the problem and whether cyclic boundaries or not
       INTEGER Ntot, Nspec, dimens(*), cyclic(*)

c maximal number of indices, sparsity arrays
       INTEGER nnz, ian(*), jan(*)
c
       INTEGER N, I, J, ij, K, L, M , isp

c check input
       IF (INT(Ntot/Nspec)*Nspec .NE. Ntot) THEN
         call rexit
     &("cannot generate sparse jacobian - N and nspec not compatible")
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
              jan(ij) = isp + j*dimens(1)
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

      SUBROUTINE interact (ij, NNZ, ian, jan, M, ival)
      INTEGER IJ, NNZ, ival, ian(*), jan(*), i, isave


      isave = 1
c check if not yet present for current state
      DO I = ian(M), ij-1
       IF (jan(I) .EQ. ival) THEN
	       isave = 0
	       exit
       ENDIF
	    ENDDO

c save 
	    IF (isave .EQ. 1) THEN
        IF (ij .GT. nnz) THEN
          CALL rexit
     &("cannot generate sparse jacobian - not enough room for nonzeros")        
        ENDIF
        jan(ij) = ival
        ij = ij +1      
      ENDIF
	 
      END SUBROUTINE interact
       
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
C A(I,J,K) in vector: (I-1)*dim(2)*dim(3) + (J-1)*dim(3) + K         *  
c--------------------------------------------------------------------*

c total number of state variables, number of different species
c dimensions of the problem and whether cyclic boundaries or not
       INTEGER Ntot, Nspec, dimens(*), cyclic(*)

c maximal number of indices, sparsity arrays
       INTEGER nnz, ian(*), jan(*)
c
       INTEGER N, I, J, ij, K, L, M, isp, im

c check input
       IF (INT(Ntot/Nspec)*Nspec .NE. Ntot) THEN
         call rexit
     &("cannot generate sparse jacobian - N and nspec not compatible")
       ENDIF

c number of boxes
       N = dimens(1)*dimens(2)*dimens(3)

       ij     = 1
       ian(1) = 1

       DO i = 1, Nspec
         isp = (i-1)*N
         DO j = 1, dimens(1)
           DO k = 1, dimens(2)
             DO ll = 1, dimens(3)
               M = isp +(j-1)*dimens(2)*dimens(3)+(k-1)*dimens(3)+ll

c interactions with current, upstream and downstream boxes (ival = M, M+1, M-1)
               CALL interact(ij, nnz, ian, jan, M, M) 

               IF (ll<dimens(3)) THEN
                 CALL interact(ij, nnz, ian, jan, M, M+1) 

               ELSE IF (cyclic(3) == 1 .AND. dimens(3) > 2) THEN
	           im = isp + (j-1)*dimens(2)*dimens(3)+(k-1)*dimens(3)+1   
                 CALL interact(ij, nnz, ian, jan, M, im) 
               ENDIF
              
			 IF (ll >1) THEN
                 CALL interact(ij, nnz, ian, jan, M, M-1) 
               ELSE IF (cyclic(3) == 1 .AND. dimens(3) > 2) THEN
                 im = isp + (j-1)*dimens(2)*dimens(3)+k*dimens(3) 
                 CALL interact(ij, nnz, ian, jan, M, im) 
               ENDIF
 
               IF (k<dimens(2)) THEN
                 CALL interact(ij, nnz, ian, jan, M, M+dimens(3)) 
               ELSE IF (cyclic(2) == 1 .AND. dimens(2) > 2) THEN
                 im = isp + (j-1)*dimens(2)*dimens(3)+ll 
                 CALL interact(ij, nnz, ian, jan, M, im) 
               ENDIF
            
               IF (k >1) THEN
                 CALL interact(ij, nnz, ian, jan, M, M -dimens(3))
               ELSE IF (cyclic(2) == 1 .AND. dimens(2) > 2) THEN
                 im = isp +j*dimens(2)*dimens(3) -dimens(3)+ll
                 CALL interact(ij, nnz, ian, jan, M, im) 
               ENDIF

               IF (j<dimens(1)) THEN
                 im = M+dimens(2)*dimens(3)
                 CALL interact(ij, nnz, ian, jan, M, im) 
               ELSE IF (cyclic(1) == 1 .AND. dimens(1) > 2) THEN
                 im = isp +(k-1)*dimens(3)+ll
                 CALL interact(ij, nnz, ian, jan, M, im) 
               ENDIF

               IF (j >1) THEN
                 im = M-dimens(2)*dimens(3)
                 CALL interact(ij, nnz, ian, jan, M, im) 
               ELSE IF (cyclic(1) == 1 .AND. dimens(1) > 2) THEN
	           im = isp +(dimens(1)-1)*dimens(2)*dimens(3)+                    &
     &		            (k-1)*dimens(3)+ll
                 CALL interact(ij, nnz, ian, jan, M, im) 
               ENDIF


c interactions with other species in the same box
               DO L = 1, Nspec
                 IF (L == i) cycle
                 im=(L-1)*N+(j-1)*dimens(2)*dimens(3)+                           &
     &                      (k-1)*dimens(3)+ll
                 CALL interact(ij, nnz, ian, jan, M, im) 
               ENDDO

               ian(M+1) = ij
             ENDDO
	     ENDDO
         ENDDO
      ENDDO
      nnz = ij -1

c
      END SUBROUTINE sparse3d
