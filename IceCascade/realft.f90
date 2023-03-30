! REALFT

! taken from numerical recipes (see book for further information)

! subroutines called:
! NONE

        module m_realft
          contains
          SUBROUTINE REALFT(DATA,N,ISIGN)
          use m_four1
          implicit none
          REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
          real(8) :: DATA,C1,C2, H1I, H1R, H2I, H2R, N2P3, WIS, WRS

          integer(4) :: I, N, I1, I2, I3, I4, ISIGN

          DIMENSION DATA(*)
          THETA=6.28318530717959D0/2.0D0/DBLE(N)
          C1=0.5_8
          IF (ISIGN.EQ.1) THEN
            C2=-0.5_8
            CALL FOUR1(DATA,N,+1)
          ELSE
            C2=0.5_8
            THETA=-THETA
          ENDIF
          WPR=-2.0D0*DSIN(0.5D0*THETA)**2
          WPI=DSIN(THETA)
          WR=1.0D0+WPR
          WI=WPI
          N2P3=dble(2*N+3)
          DO 11 I=2,N/2+1
            I1=2*I-1
            I2=I1+1
            I3=int(N2P3)-I2
            I4=I3+1
            WRS=dble(WR)
            WIS=dble(WI)
            H1R=C1*(DATA(I1)+DATA(I3))
            H1I=C1*(DATA(I2)-DATA(I4))
            H2R=-C2*(DATA(I2)+DATA(I4))
            H2I=C2*(DATA(I1)-DATA(I3))
            DATA(I1)=H1R+WRS*H2R-WIS*H2I
            DATA(I2)=H1I+WRS*H2I+WIS*H2R
            DATA(I3)=H1R-WRS*H2R+WIS*H2I
            DATA(I4)=-H1I+WRS*H2I+WIS*H2R
            WTEMP=WR
            WR=WR*WPR-WI*WPI+WR
            WI=WI*WPR+WTEMP*WPI+WI
11          CONTINUE
          IF (ISIGN.EQ.1) THEN
            H1R=DATA(1)
            DATA(1)=H1R+DATA(2)
            DATA(2)=H1R-DATA(2)
          ELSE
            H1R=DATA(1)
            DATA(1)=C1*(H1R+DATA(2))
            DATA(2)=C1*(H1R-DATA(2))
            CALL FOUR1(DATA,N,-1)
          ENDIF
          RETURN
          END subroutine realft
        end module m_realft



