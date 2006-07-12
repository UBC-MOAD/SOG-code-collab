SUBROUTINE find_upwell(mm,w,wo,S)

  USE mean_param

  IMPLICIT NONE
  
  TYPE(gr_d),INTENT(IN)::mm      !grid  
  DOUBLE PRECISION, INTENT(IN) :: wo
  DOUBLE PRECISION, DIMENSION(1:mm%M+1),INTENT(OUT):: w
  DOUBLE PRECISION, DIMENSION(1:mm%M+1),INTENT(IN):: S
  DOUBLE PRECISION, DIMENSION(1:mm%M+1)::fwc
  REAL:: d
  INTEGER:: index,xx,f75,d75

!-----------------------------------
!susan's old upwelling plan
!   do index=1,mm%M+1
!      w(index)= wo*index/(mm%M+1)
!   enddo
!-----------------------------------
!susan's new upwelling plan, March 29, 2006

!find the freshwater content
fwc(1)=30-S(index)
DO index=2,mm%M+1
   fwc(index)=30-S(index)+fwc(index-1)
ENDDO

!find the depth of 75% freshness
f75=floor(fwc(mm%M+1)*0.67);   
DO index=2,mm%M+1
   IF (floor(fwc(index))==f75) THEN
      d75=index;
   ELSE IF (floor(fwc(index))+1==f75) THEN 
      d75=index;
   ENDIF
ENDDO

!choose one of three choices for d

xx=2

IF (xx==1) THEN
   d=37.5;
ELSE IF (xx==2) THEN 
   d=2.51*d75;
ELSE IF (xx==3) THEN
    d=2.51*d75   
   IF (d > 80) THEN
      d=80   
   ENDIF
ENDIF

DO index=1,mm%M+1
   IF (index<d) THEN
      w(index)=wo*(1-((1-index/d)**2.0))
   ELSE 
      w(index)=wo
   ENDIF
ENDDO


!open(666,file="output/upwell.dat")
!do index=1,mm%M+1
!write(666,*)w(index)
!enddo
!close(666)

end

