      subroutine gband(a,d,x,n,m,eps,ierr,ifrst)
!     ********** GBAND **************************************

      implicit real*8(a-h,o-z)
      dimension a(1), d(1), x(1)
      !REAL A(2100),D(100),x(100),eps
      !INTEGER ierr,ifrst,n,m



!
!     PURPOSE: solution of system of equations with band matrix by standard
!              elimination without pivoting.
!
!     From "Petroleum Reservoir Simulation", K. Aziz and A. Settari, 1979.
!
!
!        Input/Output:
!       a : one-dimensional array !ontaining the band of the matrix
!           stored by rows. The required dimension of a is
!
!                    n*(2*m+1)-m*m-m
!
!     Input Only:
!       m   : Number of diagonals above the main diagonal.number of
!             diagonals below the main diagonal is also m, therefore
!             the total bandwidth is 2*m+1
!       n   : Number of equations (unknowns)
!       d   : r.h.s. vector
!       eps : The element on the main diagonal is compared to eps
!             during elimination.if it is smaller, value of the counter
!             ierr is incremented.
!       ifrst: =0, matrix is inverted and the inverse is stored
!                  in place of the original matrix .
!              >0, routine assumes that the matrix has been inverted
!                  and will calculate solution corresponding
!                  to the new r.h.s. vector.
!     Output Only:
!       ierr: number of times the element on the main diagonal was
!             smaller than eps
!       x   : solution vector
!
!     Update: R.C. Wattenbarger, changed impliled goto's to if/thens
!
        ierr = 0
        j    = 1
  	  do 10 i=1,n
           ie   = m
           ieaux= m
           if(m.gt.n-i) ie = n-i
           if(m.gt.i  ) ieaux= i
           ie1 = ie+ieaux
           mbig= ie
           j1  = j+ie1
           j2  = j1
           if(ifrst.eq.0 .and. abs(a(j))-eps.le.0.0)ierr=ierr+1
           if(mbig.gt.0) then
             do 20 j0=1,mbig
                ss = a(j1)/a(j)
                if(ifrst.eq.0) then
                   do 30 k=1,mbig
                     j1k = j1+k
                     jk  = j+k
                     a(j1k)=a(j1k) -a(jk)*ss
   30              continue
                endif
                iaux=j0+i
                d(iaux)=d(iaux)-d(i)*ss
                ie   = m
                ieaux= m
                if(m.gt.n-iaux) ie = n-iaux
                if(m.gt.iaux  ) ieaux= iaux
                ie1 = ie+ieaux
                j1  = j1+ie1
   20        continue
          endif
          j=j2+1
   10   continue
        j   = j-m-1
        np1 = n+1
        do 40 iinv=1,n
           i = np1-iinv
           ie= m
           if(i+m-n.gt.0)ie=n-i
           mbig = ie
           x(i) = d(i)
           if(mbig.gt.0) then
             do 50 k=1,mbig
                ik  = i+k
                jk  = j+k
                x(i)= x(i)-x(ik)*a(jk)
   50        continue
           endif
           x(i)= x(i)/a(j)
           ie  = m
           ieaux=m
           if(m.gt.np1-i) ie   = np1-i
           if(m.gt.i-1  ) ieaux= i-1
           ie1=ie+ieaux
           j  =j-ie1-1
   40   continue
        return
        end
