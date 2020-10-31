        Subroutine SOLVER(NX,NY,NZ,FO,FW,A,B,C,D,E,F,G,GB,RHS,DEL)
        Implicit Real*8(A-H,O-Z)

        Dimension GB(23098),RHS(250),DEL(250)
        Dimension A(NX,NY,NZ,4),B(NX+1,NY,NZ,4),C(0:NX,NY,NZ,4),
     &		  D(NX,NY+1,NZ,4),E(NX,0:NY,NZ,4),F(NX,NY,NZ+1,4),
     &		  G(NX,NY,0:NZ,4)
        Dimension FO(NX,NY,NZ),FW(NX,NY,NZ)

C ---------------------------------------------------------------------
C Set band matrix ............
C ---------------------------------------------------------------------
      n=2*nx*ny*nz
      m1=2*nz
      m2=ny*m1
      m=m2+1
      if(nx.eq.1 .and. ny.eq.1) then
         m=3
      else if(nx.eq.1) then
         m=m1+1
      else if(ny.eq.1) then
         m1=2
      end if
      ngb=n*(2*m+1)-m*(m+1)
      do i=1,ngb
         gb(i)=0.
      enddo
      irhs=0
      do i=1,nx
      do j=1,ny
      do k=1,nz
         irhs=irhs+1
         rhs(irhs)=-FO(i,j,k)
         irhs=irhs+1
         rhs(irhs)=-FW(i,j,k)
         ii=k+nz*(j-1)+nz*ny*(i-1)
         do l=1,2
            imat=2*ii-(2-l)
            igb=(imat-1)*(2*m+1)+m+1
            if(imat.le.m) then
               do im=1,imat
                  igb=igb-(m+1-im)
               enddo
            else
               igb=igb-m*(m+1)/2
            end if
            if(imat.gt.n-m+1) then
               do im=1,imat-(n-m+1)
                  igb=igb-im
               enddo
            end if
            if(l.eq.2) igb=igb-1
            l2=2*l
            gb(igb)=a(i,j,k,l2-1)
            gb(igb+1)=a(i,j,k,l2)
            if(k.ge.2) then
               gb(igb-2)=f(i,j,k,l2-1)
               gb(igb-1)=f(i,j,k,l2)
            end if
            if(k.le.nz-1) then
               gb(igb+2)=g(i,j,k,l2-1)
               gb(igb+3)=g(i,j,k,l2)
            end if
            if(j.ge.2) then
               gb(igb-m1)=d(i,j,k,l2-1)
               gb(igb-m1+1)=d(i,j,k,l2)
            end if
            if(j.le.ny-1) then
               gb(igb+m1)=e(i,j,k,l2-1)
               gb(igb+m1+1)=e(i,j,k,l2)
            end if
            if(i.ge.2) then
               gb(igb-m2)=b(i,j,k,l2-1)
               gb(igb-m2+1)=b(i,j,k,l2)
            end if
            if(i.le.nx-1) then
               gb(igb+m2)=c(i,j,k,l2-1)
               gb(igb+m2+1)=c(i,j,k,l2)
            end if
         enddo
      enddo
      enddo
      enddo
c
      call gband(gb,rhs,del,n,m,eps,ierr,0)
c
      return
c
      end
