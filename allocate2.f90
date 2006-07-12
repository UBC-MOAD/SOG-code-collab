SUBROUTINE allocate2(STAT)      

      USE mean_param
      USE declarations

      IMPLICIT NONE
     
      INTEGER, DIMENSION(20), INTENT(IN OUT)::STAT


            IF (time_step == 1 .AND. count == 1) THEN

               ALLOCATE(G_shape%s(0:h%i),G_shape%m(0:h%i),G_shape%t(0:h%i), STAT = alloc_stat(16))
               ALLOCATE(K%u%ML(0:h%i),K%s%ML(0:h%i),K%t%ML(0:h%i), STAT = alloc_stat(17))

            ELSE
               DEALLOCATE(G_shape%s,G_shape%m,G_shape%t,K%u%ML,K%s%ML,K%t%ML)

               ALLOCATE(G_shape%s(0:h%i),G_shape%m(0:h%i),G_shape%t(0:h%i), STAT = alloc_stat(16))
               ALLOCATE(K%u%ML(0:h%i),K%s%ML(0:h%i),K%t%ML(0:h%i), STAT = alloc_stat(17))
 
            END IF

END SUBROUTINE allocate2
