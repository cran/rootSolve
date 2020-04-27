c--------------------------------------------------------------------*
c SPARSITY of 2-D or 3-D reaction-transport problems with mapping
c--------------------------------------------------------------------*

      SUBROUTINE updatejan (ij, ii, nnz, jan, pres)
      IMPLICIT NONE      
c--------------------------------------------------------------------*
c If ii is present, updates the array jan                            *
c--------------------------------------------------------------------*
       INTEGER ii, ij, nnz, jan(*), pres(*)
c

          
       IF (pres(II) > 0) THEN
           jan(ij) = pres(II)
           ij      = ij +1
           IF (ij .GT. nnz) THEN
            CALL rexit
     &("cannot generate sparse jacobian - not enough room for nonzeros")
           ENDIF
       ENDIF

      END SUBROUTINE updatejan

      SUBROUTINE sparse2dmap (Ntot, Nspec, dimens, cyclic,                     &
     &                     nnz, ian, jan, pres)

      IMPLICIT NONE      
c--------------------------------------------------------------------*
c Determines the sparsity structure of the Jacobian in 2-D reaction- *
c transport problems with a mapping matrix.                          *
c the states that are present are in 'pres'                          *
C IAN  and JAN: see comments in sparse1D                             *
C A(I,J) in vector: at position (I-1)*dim(2) + J                     *
C                                                                    *
c--------------------------------------------------------------------*

c total number of state variables, number of different species
c dimensions of the problem and whether cyclic boundaries or not
       INTEGER Ntot, Nspec, Mnew, dimens(*), cyclic(*)

c maximal number of indices, sparsity arrays
       INTEGER nnz, ian(*), jan(*)

c pres has the new numbers of those present, -1 or 0 if absent
       INTEGER pres(*)
c
       INTEGER N, I, J, ij, K, L, M , isp, NN

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
             Mnew = pres(M)
             IF (Mnew > 0) THEN
c interactions with current, upstream and downstream boxes
               CALL updatejan(ij, M, nnz, jan, pres)

              IF (k < dimens(2)) THEN
                CALL updatejan(ij, M+1, nnz, jan, pres)

              ELSE IF (cyclic(2) == 1) THEN
                NN = isp + (j-1)*dimens(2) +1
                CALL updatejan(ij, NN, nnz, jan, pres)

              ENDIF

              IF (j<dimens(1)) THEN  
                CALL updatejan(ij, M+dimens(2), nnz, jan, pres)
              ELSE IF (cyclic(1) == 1) THEN
                CALL updatejan(ij, isp + K, nnz, jan, pres)
              ENDIF

              IF (j >1) THEN
                CALL updatejan(ij, M-dimens(2), nnz, jan, pres)
              ELSE IF (cyclic(1) == 1) THEN
                NN = isp + (dimens(1)-1)*dimens(2)+ K
                CALL updatejan(ij, NN, nnz, jan, pres)
              ENDIF

              IF (k >1) THEN
                CALL updatejan(ij, M-1, nnz, jan, pres)
              ELSE IF (cyclic(2) == 1) THEN
                NN = isp + j*dimens(1) 
                CALL updatejan(ij, NN, nnz, jan, pres)
              ENDIF

c interactions with other species in the same box
              DO L = 1, Nspec
                IF (L == i) cycle
                NN = (L-1)*N+(j-1)*dimens(2)+k
                CALL updatejan(ij, NN, nnz, jan, pres)
              ENDDO
       
              ian(Mnew + 1) = ij
             ENDIF
           ENDDO
         ENDDO
       ENDDO
       nnz = ij -1
c
      END SUBROUTINE sparse2dmap



      SUBROUTINE sparse3dmap (Ntot, Nspec, dimens, cyclic,                     &
     &      nnz, ian, jan, pres)

c--------------------------------------------------------------------*
c Determines the sparsity structure of the Jacobian in 3-D PDE       *
c problems.                                                          *
c                                                                    *
C IAN  and JAN: see comments in sparse1D                             *
C                                                                    *
C A(I,J,K) in vector: (I-1)*dim(2)*dim(3) + (J-1)*dim(3) + K         *  
c--------------------------------------------------------------------*
       IMPLICIT NONE
c total number of state variables, number of different species
c dimensions of the problem and whether cyclic boundaries or not
       INTEGER Ntot, Nspec, dimens(*), cyclic(*), pres(*)

c maximal number of indices, sparsity arrays
       INTEGER nnz, ian(*), jan(*)
c
       INTEGER N, I, J, ij, K, L, M, ll, isp, im, Mnew

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
               Mnew = pres(M)
               IF (Mnew > 0) THEN

c interactions with current, upstream and downstream boxes (ival = M, M+1, M-1)
                 CALL updatejan(ij, M, nnz, jan, pres) 

                 IF (ll < dimens(3)) THEN
                   CALL updatejan(ij, M+1, nnz, jan, pres) 

                 ELSE IF (cyclic(3) == 1 .AND. dimens(3) > 2) THEN
               im = isp + (j-1)*dimens(2)*dimens(3)+(k-1)*dimens(3)+1   
                   CALL updatejan(ij, im, nnz, jan, pres) 
                 ENDIF
              
                 IF (ll > 1) THEN
                   CALL updatejan(ij, M-1, nnz, jan, pres) 
                 ELSE IF (cyclic(3) == 1 .AND. dimens(3) > 2) THEN
                   im = isp + (j-1)*dimens(2)*dimens(3)+k*dimens(3) 
                   CALL updatejan(ij, im, nnz, jan, pres) 
                 ENDIF
 
                 IF (k < dimens(2)) THEN
                   CALL updatejan(ij, M+dimens(3), nnz, jan, pres) 
                 ELSE IF (cyclic(2) == 1 .AND. dimens(2) > 2) THEN
                   im = isp + (j-1)*dimens(2)*dimens(3)+ll 
                   CALL updatejan(ij, im, nnz, jan, pres) 
                 ENDIF
            
                 IF (k > 1) THEN
                   CALL updatejan(ij, M-dimens(3), nnz, jan, pres)
                 ELSE IF (cyclic(2) == 1 .AND. dimens(2) > 2) THEN
                   im = isp +j*dimens(2)*dimens(3) -dimens(3)+ll
                   CALL updatejan(ij, im, nnz, jan, pres) 
                 ENDIF

                 IF (j < dimens(1)) THEN
                   im = M+dimens(2)*dimens(3)
                   CALL updatejan(ij, im, nnz, jan, pres) 
                 ELSE IF (cyclic(1) == 1 .AND. dimens(1) > 2) THEN
                   im = isp +(k-1)*dimens(3)+ll
                   CALL updatejan(ij, im, nnz, jan, pres) 
                 ENDIF

                 IF (j > 1) THEN
                   im = M-dimens(2)*dimens(3)
                   CALL updatejan(ij, im, nnz, jan, pres) 
                 ELSE IF (cyclic(1) == 1 .AND. dimens(1) > 2) THEN
                im = isp +(dimens(1)-1)*dimens(2)*dimens(3)+                    &
     &                  (k-1)*dimens(3)+ll
                   CALL updatejan(ij, im, nnz, jan, pres) 
                 ENDIF

c interactions with other species in the same box
                 DO L = 1, Nspec
                   IF (L == i) cycle
                   im=(L-1)*N+(j-1)*dimens(2)*dimens(3)+                        &
     &                        (k-1)*dimens(3)+ll
                   CALL updatejan(ij, im, nnz, jan, pres) 
                 ENDDO

                 ian(MNew+1) = ij
               ENDIF
             ENDDO
           ENDDO
         ENDDO
      ENDDO
      nnz = ij -1

c
      END SUBROUTINE sparse3dmap
