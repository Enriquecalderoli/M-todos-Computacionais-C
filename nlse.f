
        PROGRAM chinatown
        implicit none

        integer nsit
        parameter (nsit = 301)  ! modify nsit in tcge

        complex*16 w1(nsit),Q1(nsit,2),x1(nsit),y1(nsit),jimag
        real*8 cosi1(nsit),energ1(nsit),occup1(nsit)
        real*8 d1,part1,hop1,msd1,weg1,shan1,dsq1,norma1
        real*8 delta1,hamilt1,cent1,aux1

        complex*16 w2(nsit),Q2(nsit,2),x2(nsit),y2(nsit)
        real*8 cosi2(nsit),energ2(nsit),occup2(nsit)
        real*8 d2,part2,hop2,msd2,weg2,shan2,dsq2,norma2
        real*8 delta2,hamilt2,cent2,aux2

        real*8 campo,epsilo,hopping,U1,U2,pi,pe,crit
        real*8 timax,time,timestep

        integer np,L,iter,other,i,k,caso,const,limit


        open (unit=10,file='datain',status='unknown')
        open (unit=11,file='dataout',status='unknown')
        open (unit=12,file='status',status='unknown')
        open (unit=13,file='anderson',status='unknown')
        open (unit=14,file='msd',status='unknown')
        open (unit=15,file='hamilt',status='unknown')
        open (unit=16,file='profile',status='unknown')
        open (unit=17,file='wegner',status='unknown')
        open (unit=18,file='shannon',status='unknown')
        open (unit=19,file='centroide',status='unknown')
        open (unit=20,file='desordem', status='unknown')
        open (unit=21,file='hopping', status='unknown')
        open (unit=22,file='particip', status='unknown')
c        open (unit=23,file='perfil', status='unknown')


c	 open (unit=30,file='datain2',status='unknown')
        open (unit=31,file='dataout2',status='unknown')
        open (unit=32,file='status2',status='unknown')
        open (unit=33,file='anderson2',status='unknown')
        open (unit=34,file='msd2',status='unknown')
        open (unit=35,file='hamilt2',status='unknown')
        open (unit=36,file='profile2',status='unknown')
        open (unit=37,file='wegner2',status='unknown')
        open (unit=38,file='shannon2',status='unknown')
        open (unit=39,file='centroide2',status='unknown')
        open (unit=40,file='desordem2', status='unknown')
        open (unit=41,file='hopping2', status='unknown')
        open (unit=42,file='particip2', status='unknown')
c        open (unit=43,file='perfil2', status='unknown')


c****** read initial data from datain and close file 10


        read (10,*) caso,other,timestep,timax,limit
        read (10,*) campo,epsilo,hopping,U1,U2
        close (unit=10)


c****** save input data in status for each run


        write (11,900)
        write (11,100) nsit,caso
        write (11,110) other,timestep,timax,limit
        write (11,111) campo,epsilo,hopping,U1,U2


        write (31,900)
        write (31,100) nsit,caso
        write (31,110) other,timestep,timax,limit
        write (31,111) campo,epsilo,hopping,U1,U2


c****** initialize wannier amplitudes, occups, energies, off-diag Qs and norms


        jimag = (0.d0,1.d0)
        np = nsit
        L  = nsit/2 + 1
        pi = 4.d0*datan(1.d0)
        pe = (dsqrt(5.d0)+1.d0)/2.d0
        crit = 1.d-20
        const = 20
C        limit = 50
        write (11,*) L, L-limit, L+limit


c****** start with caso = 1 or 3 or 5


        do i =1,nsit
         w1(i) = 1.d-30
         w2(i) = 1.d-30
        enddo
        w1(L) = 1.d0
        w2(L) = 1.d0

        if (caso.eq.3) then
          w1(L-1) =1.d0/dsqrt(3.d0)
          w1(L  ) =1.d0/dsqrt(3.d0)
          w1(L+1) =1.d0/dsqrt(3.d0)
          w2(L-1) =1.d0/dsqrt(3.d0)
          w2(L  ) =1.d0/dsqrt(3.d0)
          w2(L+1) =1.d0/dsqrt(3.d0)
        endif

        if (caso.eq.5) then
          w1(L-2) = 1.d0/dsqrt(7.5d0)
          w1(L-1) = 1.d0/dsqrt(5.d0)
          w1(L  ) = 1.d0/dsqrt(3.d0)
          w1(L+1) = 1.d0/dsqrt(5.d0)
          w1(L+2) = 1.d0/dsqrt(7.5d0)
          w2(L-2) = 1.d0/dsqrt(7.5d0)
          w2(L-1) = 1.d0/dsqrt(5.d0)
          w2(L  ) = 1.d0/dsqrt(3.d0)
          w2(L+1) = 1.d0/dsqrt(5.d0)
          w2(L+2) = 1.d0/dsqrt(7.5d0)
        endif


c*****  calculate initial values Q1(i,1) and Q2(i,1) off-diagonal terms


        norma1=0.d0
        do i =1,nsit
              occup1(i)=abs(w1(i))**2
              cosi1(i) = dcos(2.d0*pi*pe*(i-L))
              energ1(i)=epsilo*cosi1(i)+campo*(i-L)
              Q1(i,1) =-jimag*timestep*hopping/4.d0
              norma1=norma1+occup1(i)
        enddo

        norma2=0.d0
        do k =1,nsit
              occup2(k)=abs(w2(k))**2
              cosi2(k) = dcos(2.d0*pi*pe*(k-L))
              energ2(k)=epsilo*cosi2(k)+campo*(k-L)
              Q2(k,1) =-jimag*timestep*hopping/4.d0
              norma2=norma2+occup2(k)
        enddo


c****** write down initial configuration at time = 0


        iter=0
        time=0.d0
        delta1=0.d0
        delta2=0.d0

        write (11,112)
        write (11,122) (energ1(i),i=L-3,L+3)
        write (11,123)

        write (31,112)
        write (31,122) (energ2(k),k=L-3,L+3)
        write (31,123)

        weg1=0.d0
        msd1=0.d0
        shan1=0.d0
        hamilt1=0.d0
        cent1=0.d0
        d1=0.d0
        hop1=0.d0

        weg2=0.d0
        msd2=0.d0
        shan2=0.d0
        hamilt2=0.d0
        cent2=0.d0
        d2=0.d0
        hop2=0.d0

c **********************************

        do i =1,nsit-1
            if (occup1(i).gt.crit) then
                     dsq1=dsqrt(occup1(i))
                     weg1=weg1+occup1(i)*occup1(i)
                     msd1=msd1+(i-L)*(i-L)*occup1(i)
                     shan1=shan1-2.d0*occup1(i)*dlog(dsq1)
                     hop1=hop1+dreal(w1(i)*conjg(w1(i+1))+
     &               conjg(w1(i))*w1(i+1))
                     hamilt1=hamilt1+energ1(i)*occup1(i)+
     &               0.5d0*U2*occup1(i)*occup2(i)-hopping*
     &          dreal(w1(i)*conjg(w1(i+1))+conjg(w1(i))*w1(i+1))
            endif
            if (occup2(i).gt.crit) then
                     dsq2=dsqrt(occup2(i))
                     weg2=weg2+occup2(i)*occup2(i)
                     msd2=msd2+(i-L)*(i-L)*occup2(i)
                     shan2=shan2-2.d0*occup2(i)*dlog(dsq2)
                     hop2=hop2+dreal(w2(i)*conjg(w2(i+1))+
     &               conjg(w2(i))*w2(i+1))
                     hamilt2=hamilt2+energ2(i)*occup2(i)+
     &               0.5d0*U2*occup2(i)*occup1(i)-hopping*
     &          dreal(w2(i)*conjg(w2(i+1))+conjg(w2(i))*w2(i+1))
            endif
            cent1=cent1+(i-L)*occup1(i)
            cent2=cent2+(i-L)*occup2(i)
            d1=d1+cosi1(i)*occup1(i)
            d2=d2+cosi2(i)*occup2(i)
        enddo

            if (occup1(nsit).gt.crit) then
                        dsq1=dsqrt(occup1(nsit))
                        weg1=weg1+occup1(nsit)*occup1(nsit)
                        msd1=msd1+(nsit-L)*(nsit-L)*occup1(nsit)
                        shan1=shan1-2.d0*occup1(nsit)*dlog(dsq1)
                        hamilt1=hamilt1+energ1(nsit)*occup1(nsit)+
     &         0.5d0*U2*occup1(nsit)*occup2(nsit)+0.5d0*U1*weg1
            endif
            if (occup2(nsit).gt.crit) then
                     dsq2=dsqrt(occup2(nsit))
                     weg2=weg2+occup2(nsit)*occup2(nsit)
                     msd2=msd2+(nsit-L)*(nsit-L)*occup2(nsit)
                     shan2=shan2-2.d0*occup2(nsit)*dlog(dsq2)
                     hamilt2=hamilt2+energ2(nsit)*occup2(nsit)+
     &         0.5d0*U2*occup2(nsit)*occup1(nsit)+0.5d0*U1*weg2
            endif
            cent1=(cent1+(nsit-L)*occup1(nsit))/L
            cent2=(cent2+(nsit-L)*occup2(nsit))/L
            d1=d1+cosi1(nsit)*occup1(nsit)
            d2=d2+cosi2(nsit)*occup2(nsit)

              part1=weg1
              weg1=1./weg1

              part2=weg2
              weg2=1./weg2


c************************************


        write (11,220) time,(occup1(i),i=L-2,L+2),delta1,norma1
        write (12,310) time,msd1,weg1,shan1,occup1(L)
        write (13,500) time,occup1(L)      ! anderson
        write (14,500) time,msd1           ! msd
        write (15,500) time,hamilt1        ! hamiltoniano
        write (17,500) time,weg1           ! wegner
        write (18,500) time,shan1          ! shannon
        write (19,500) time,cent1          ! centroide
        write (20,500) time,d1             ! desordem
        write (21,500) time,hop1           ! hopping
        write (22,500) time,part1          ! particip
c        write (23,700) time,(occup1(i),i=L-limit,L+limit)

        do i = L-limit , L+limit
           write (16,600) time,i,occup1(i)
        enddo


        write (31,220) time,(occup2(k),k=L-2,L+2),delta2,norma2
        write (32,310) time,msd2,weg2,shan2,occup2(L)
        write (33,500) time,occup2(L)
        write (34,500) time,msd2
        write (35,500) time,hamilt2
        write (37,500) time,weg2
        write (38,500) time,shan2
        write (39,500) time,cent2
        write (40,500) time,d2
        write (41,500) time,hop2
        write (42,500) time,part2
c        write (43,700) time,(occup2(i),i=L-limit,L+limit)

        do k = L-limit , L+limit
           write (36,600) time,k,occup2(k)
        enddo



c****** here starts loop of time evolution with a Crank-Nicolson scheme


        time=timestep


 10     do i=1,nsit
                Q1(i,2)=0.5d0+ jimag*timestep*(energ1(i)+ 
     &          U1*occup1(i) + U2*occup2(i))/4.d0
                y1(i)=w1(i)
        enddo

        do k=1,nsit
                Q2(k,2)=0.5d0+ jimag*timestep*(energ2(k)+ 
     &          U1*occup2(k) + U2*occup1(k))/4.d0
                y2(k)=w2(k)
        enddo


c****** tcge: gaussian elimination for tridiagonal complex matrix

        call tcge (Q1,y1,x1,np)
        call tcge (Q2,y2,x2,np)


c****** up-date wannier amplitudes, normas, deltas, occups

        norma1=0.d0
        delta1=0.d0
        do i=1,nsit
           w1(i)=x1(i)-w1(i)
           aux1=abs(w1(i))**2
           norma1=norma1+aux1
           delta1=delta1+dabs(aux1-occup1(i))
           occup1(i)=aux1
        enddo


        norma2=0.d0
        delta2=0.d0
        do k=1,nsit
           w2(k)=x2(k)-w2(k)
           aux2=abs(w2(k))**2
           norma2=norma2+aux2
           delta2=delta2+dabs(aux2-occup2(k))
            occup2(k)=aux2
        enddo


c****** calculate variables msd, wegner, shannon, centroid, hamiltonian


        iter = iter +1


        weg1=0.d0
        msd1=0.d0
        shan1=0.d0
        hamilt1=0.d0
        cent1=0.d0
        d1=0.d0
        hop1=0.d0

        weg2=0.d0
        msd2=0.d0
        shan2=0.d0
        hamilt2=0.d0
        cent2=0.d0
        d2=0.d0
        hop2=0.d0


c **********************************


        do i =1,nsit-1
            if (occup1(i).gt.crit) then
                     dsq1=dsqrt(occup1(i))
                     weg1=weg1+occup1(i)*occup1(i)
                     msd1=msd1+(i-L)*(i-L)*occup1(i)
                     shan1=shan1-2.d0*occup1(i)*dlog(dsq1)
                     hop1=hop1+dreal(w1(i)*conjg(w1(i+1))+
     &               conjg(w1(i))*w1(i+1))
                     hamilt1=hamilt1+energ1(i)*occup1(i)+
     &               0.5d0*U2*occup1(i)*occup2(i)-hopping*
     &          dreal(w1(i)*conjg(w1(i+1))+conjg(w1(i))*w1(i+1))
            endif
            if (occup2(i).gt.crit) then
                     dsq2=dsqrt(occup2(i))
                     weg2=weg2+occup2(i)*occup2(i)
                     msd2=msd2+(i-L)*(i-L)*occup2(i)
                     shan2=shan2-2.d0*occup2(i)*dlog(dsq2)
                     hop2=hop2+dreal(w2(i)*conjg(w2(i+1))+
     &               conjg(w2(i))*w2(i+1))
                     hamilt2=hamilt2+energ2(i)*occup2(i)+
     &               0.5d0*U2*occup2(i)*occup1(i)-hopping*
     &          dreal(w2(i)*conjg(w2(i+1))+conjg(w2(i))*w2(i+1))
            endif
            cent1=cent1+(i-L)*occup1(i)
            cent2=cent2+(i-L)*occup2(i)
            d1=d1+cosi1(i)*occup1(i)
            d2=d2+cosi2(i)*occup2(i)
        enddo

            if (occup1(nsit).gt.crit) then
                     dsq1=dsqrt(occup1(nsit))
                     weg1=weg1+occup1(nsit)*occup1(nsit)
                     msd1=msd1+(nsit-L)*(nsit-L)*occup1(nsit)
                     shan1=shan1-2.d0*occup1(nsit)*dlog(dsq1)
                     hamilt1=hamilt1+energ1(nsit)*occup1(nsit)+
     &         0.5d0*U2*occup1(nsit)*occup2(nsit)+0.5d0*U1*weg1
            endif
            if (occup2(nsit).gt.crit) then
                     dsq2=dsqrt(occup2(nsit))
                     weg2=weg2+occup2(nsit)*occup2(nsit)
                     msd2=msd2+(nsit-L)*(nsit-L)*occup2(nsit)
                     shan2=shan2-2.d0*occup2(nsit)*dlog(dsq2)
                     hamilt2=hamilt2+energ2(nsit)*occup2(nsit)+
     &         0.5d0*U2*occup2(nsit)*occup1(nsit)+0.5d0*U1*weg2
             endif
             cent1=(cent1+(nsit-L)*occup1(nsit))/L
             cent2=(cent2+(nsit-L)*occup2(nsit))/L
             d1=d1+cosi1(nsit)*occup1(nsit)
             d2=d2+cosi2(nsit)*occup2(nsit)

             part1=weg1
             weg1=1./weg1

             part2=weg2
             weg2=1./weg2



c****** writes-down results every other and go to next time slice


        if (mod(iter,other).eq.0) then

        write (11,220) time,(occup1(i),i=L-2,L+2),delta1,norma1
        write (12,310) time,msd1,weg1,shan1,occup1(L)
        write (13,500) time,occup1(L)      ! anderson
        write (14,500) time,msd1           ! msd
        write (15,500) time,hamilt1        ! hamiltoniano
        write (17,500) time,weg1           ! wegner
        write (18,500) time,shan1          ! shannon
        write (19,500) time,cent1          ! centroide
        write (20,500) time,d1             ! desordem
        write (21,500) time,hop1           ! hopping
        write (22,500) time,part1          ! particip
c        write (23,700) time,(occup1(i),i=L-limit,L+limit)


        write (31,220) time,(occup2(k),k=L-2,L+2),delta2,norma2
        write (32,310) time,msd2,weg2,shan2,occup2(L)
        write (33,500) time,occup2(L)
        write (34,500) time,msd2
        write (35,500) time,hamilt2
        write (37,500) time,weg2
        write (38,500) time,shan2
        write (39,500) time,cent2
        write (40,500) time,d2
        write (41,500) time,hop2
        write (42,500) time,part2
c        write (43,700) time,(occup2(i),i=L-limit,L+limit)

        endif


c****** write-down profile of wavepacket every const*other from L-limit to L+limit


        if (mod(iter,const*other).eq.0) then

        do i = L-limit , L+limit
            write (16,600) time,i,occup1(i)
        enddo

        do k = L-limit , L+limit
            write (36,600) time,k,occup2(k)
        enddo

        endif


c   ready to return to loop or exit  ************************


        time = time + timestep

        if (time.le.timax) goto 10



c****** formats, close and stop


 900    format ('Program chinatown')
 100    format ('Calculo para',I6,' sitios, caso',I3)
 110    format ('other,timestep,timax,limit: ',I5,F9.5,F10.3,I5)
 111    format ('campo,epsilo,hopping,U1,U2:',5F9.4)
 112    format ('energ(i),i=L-3,L+3:')
 122    format (7F12.7)
 123    format ('time,(occup(i),i=L-2,L+2),delta,norma:')
 220    format (F10.4,2x,5F11.7,1x,2E12.4)
 310    format (F10.4,5E16.8)
 500    format (F10.4,E18.8)
 600    format (F10.3,I6,E20.8)
 700    format (F10.3,100E16.6)


        close (unit=11) ! dataout
        close (unit=12) ! status
        close (unit=13) ! anderson
        close (unit=14) ! msd
        close (unit=15) ! hamilt
        close (unit=16) ! profile
        close (unit=17) ! wegner
        close (unit=18) ! shannon
        close (unit=19) ! centroide
        close (unit=20) ! desordem
        close (unit=21) ! hopping
        close (unit=22) ! particip
c        close (unit=23) ! perfil


        close (unit=31) ! dataout2
        close (unit=32) ! status2
        close (unit=33) ! anderson2
        close (unit=34) ! msd2
        close (unit=35) ! hamilt2
        close (unit=36) ! profile2
        close (unit=37) ! wegner2
        close (unit=38) ! shannon2
        close (unit=39) ! centroide2
        close (unit=40) ! desordem2
        close (unit=41) ! hopping2
        close (unit=42) ! particip2
c        close (unit=43) ! perfil2



        STOP 'terminei lchinatown'
        END

c******************************************************************


        SUBROUTINE tcge(a,b,x,n)
        integer i,n,nsit
        parameter (nsit = 301) ! modificar nsit em chinatown
        complex*16 a(n,2),b(n),x(n),coeff
        complex*16 alfa(nsit),beta(nsit)
        if (nsit.ne.n) stop 'nsit differs from n'
        do 10 i=1,n
           alfa(i)=a(i,1)
           beta(i)=a(i,2)
 10     continue
c       forward elimination
        do 20 i=2,n
           coeff=alfa(i-1)/beta(i-1)
           beta(i)=beta(i)-coeff*alfa(i-1)
           b(i)=b(i)-coeff*b(i-1)
 20     continue
c       backward elimination
        x(n)=b(n)/beta(n)
        do 30 i=n-1,1,-1
           x(i)=(b(i)-alfa(i)*x(i+1))/beta(i)
 30     continue
        return
        END
