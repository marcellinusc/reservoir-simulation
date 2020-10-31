module jacobian
    use import_data
    use functions
    use initial_conditions
    use well
    use potential
    use derivation

    implicit none
    real(kind=8)                                :: delt,delz
    real(kind=8),allocatable,dimension(:,:,:)   :: pn,swn
    real(kind=8),allocatable,dimension(:,:,:)   :: FoRes,FwRes
    real(kind=8),allocatable,dimension(:,:,:,:) :: a,b,c,d,e,f,g,acc,ss

    contains

    subroutine jacob()
        integer                                 :: i,j,k,l
        real(kind=8)                            :: pkk,pnn,skk,snn
        real(kind=8)                            :: skkim1,skkip1,skkjm1,skkjp1,skkkm1,skkkp1
        real(kind=8)                            :: pkkim1,pkkip1,pkkjm1,pkkjp1,pkkkm1,pkkkp1
        real(kind=8)                            :: btxw,btxo,ctxw,ctxo,dtyw,dtyo,etyw,etyo,ftzw,ftzo,gtzw,gtzo
        real(kind=8)                            :: acco,accw,FoTrA,FwTrA
        real(kind=8),allocatable,dimension(:)   :: dT

        if(.not.allocated(a))   allocate(a(nx,ny,nz,4))
        if(.not.allocated(b))   allocate(b(nx+1,ny,nz,4))
        if(.not.allocated(c))   allocate(c(0:nx,ny,nz,4))
        if(.not.allocated(d))   allocate(d(nx,ny+1,nz,4))
        if(.not.allocated(e))   allocate(e(nx,0:ny,nz,4))
        if(.not.allocated(f))   allocate(f(nx,ny,nz+1,4))
        if(.not.allocated(g))   allocate(g(nx,ny,0:nz,4))
        if(.not.allocated(acc)) allocate(acc(nx,ny,nz,4))
        if(.not.allocated(ss))  allocate(ss(nx,ny,nz,4))

        if(.not.allocated(pn))    allocate(pn(nx,ny,nz))
        if(.not.allocated(swn))   allocate(swn(nx,ny,nz))
        if(.not.allocated(FoRes)) allocate(FoRes(nx,ny,nz))
        if(.not.allocated(FwRes)) allocate(FwRes(nx,ny,nz))

        do k=1,nz
        do j=1,ny
        do i=1,nx
            skk=swk(i,j,k)
            snn=swn(i,j,k)
            pkk=pk(i,j,k)
            pnn=pn(i,j,k)
            !
            ! B node
            !
            if(i==1) then
            btxw=0.d0
            btxo=0.d0
                else
                pkkim1=pk(i-1,j,k)
                skkim1=swk(i-1,j,k)
                delz=0.d0
                call deriv(skk,skkim1,pkk,pkkim1,delz,ibw(i,j,k),ibw(i,j,k),pgeox,dT,btxw,btxo)
                do l=1,4
                    b(i,j,k,l)=dT(l)
                    b(nx+1,j,k,l)=0.d0
                end do
            endif
            !
            ! C node
            !
            if(i==nx) then
            ctxw=0.d0
            ctxo=0.d0
                else
                pkkip1=pk(i+1,j,k)
                skkip1=swk(i+1,j,k)
                delz=0.d0
                call deriv(skk,skkip1,pkk,pkkip1,delz,icw(i,j,k),icw(i,j,k),pgeox,dT,ctxw,ctxo)
                do l=1,4
                    c(i,j,k,l)=dT(l)
                    c(0,j,k,l)=0.d0
                end do
            endif
            !
            ! D node
            !
            if(j==1) then
            dtyw=0.d0
            dtyo=0.d0
                else
                pkkjm1=pk(i,j-1,k)
                skkjm1=swk(i,j-1,k)
                delz=0.d0
                call deriv(skk,skkjm1,pkk,pkkjm1,delz,idw(i,j,k),idw(i,j,k),pgeoy,dT,dtyw,dtyo)
                do l=1,4
                    d(i,j,k,l)=dT(l)
                    d(i,ny+1,k,l)=0.d0
                end do
            endif
            !
            ! E node
            !
            if(j==ny) then
            etyw=0.d0
            etyo=0.d0
                else
                pkkjp1=pk(i,j+1,k)
                skkjp1=swk(i,j+1,k)
                delz=0.d0
                call deriv(skk,skkjp1,pkk,pkkjp1,delz,iew(i,j,k),iew(i,j,k),pgeoy,dT,etyw,etyo)
                do l=1,4
                    e(i,j,k,l)=dT(l)
                    e(i,0,k,l)=0.d0
                end do
            endif
            !
            ! F node
            !
            if(k==1) then
            ftzw=0.d0
            ftzo=0.d0
                else
                pkkkm1=pk(i,j,k-1)
                skkkm1=swk(i,j,k-1)
                delz=-tz/nz
                call deriv(skk,skkkm1,pkk,pkkkm1,delz,ifw(i,j,k),ifo(i,j,k),pgeoz,dT,ftzw,ftzo)
                do l=1,4
                    f(i,j,k,l)=dT(l)
                    f(i,j,nz+1,l)=0.d0
                end do
            endif
            !
            ! G node
            !
            if(k==nz) then
            gtzw=0.d0
            gtzo=0.d0
                else
                pkkkp1=pk(i,j,k+1)
                skkkp1=swk(i,j,k+1)
                delz=tz/nz
                call deriv(skk,skkkp1,pkk,pkkkp1,delz,igw(i,j,k),igo(i,j,k),pgeoz,dT,gtzw,gtzo)
                do l=1,4
                    g(i,j,k,l)=dT(l)
                    g(i,j,0,l)=0.d0
                end do
            endif

            FoTrA = btxo+ctxo+dtyo+etyo+ftzo+gtzo
            FwTrA = btxw+ctxw+dtyw+etyw+ftzw+gtzw

            acco = (Vb/delt)*((phi(pkk)*(1-skk)/Bo(pkk))-phi(pnn)*(1-snn)/Bo(pnn))
            accw = Vb/delt*((phi(pkk)*skk/Bw(pkk)-(phi(pnn)*snn/Bw(pnn))))

            acc(i,j,k,1) = -Vb/delt*phi(pkk)/Bo(pkk)
            acc(i,j,k,2) =  Vb/delt*(dphidp(pkk)*(1-skk)/Bo(pkk)-(1-skk)*phi(pkk)*dBodP(pkk)/(Bo(pkk)**2))
            acc(i,j,k,3) =  Vb/delt*phi(pkk)/Bw(pkk)
            acc(i,j,k,4) =  Vb/delt*(dphidP(pkk)*skk/Bw(pkk)-skk*phi(pkk)*dBwdP(pkk)/(Bw(pkk)**2))

            if(i==wloc(1,1).AND.j==wloc(1,2).AND.k==wloc(1,3)) then
                FoRes(i,j,k) = -acco+FoTrA
                FwRes(i,j,k) = -accw+FwTrA-qw(1)
            elseif(i==wloc(2,1).AND.j==wloc(2,2).AND.k==wloc(2,3)) then
                FoRes(i,j,k) = -acco+FoTrA-qo(2)
                FwRes(i,j,k) = -accw+FwTrA-qw(2)
            else
                FoRes(i,j,k) = -acco+FoTrA
                FwRes(i,j,k) = -accw+FwTrA
            end if


                if(i==wloc(1,1).AND.j==wloc(1,2).AND.k==wloc(1,3)) then
                    ss(i,j,k,1) = dqodSw(1)
                    ss(i,j,k,2) = dqodP(1)
                    ss(i,j,k,3) = dqwdSw(1)
                    ss(i,j,k,4) = dqwdP(1)
                elseif(i==wloc(2,1).AND.j==wloc(2,2).AND.k==wloc(2,3)) then
                    ss(i,j,k,1) = dqodSw(2)
                    ss(i,j,k,2) = dqodP(2)
                    ss(i,j,k,3) = dqwdSw(2)
                    ss(i,j,k,4) = dqwdP(2)
                else
                    do l=1,4
                        ss(i,j,k,l)=0.d0
                    end do
                end if
        end do
        end do
        end do

        do i=1,nx
        do j=1,ny
        do k=1,nz
            do l=1,4
                a(i,j,k,l) = -b(i+1,j,k,l)-c(i-1,j,k,l)-d(i,j+1,k,l)-e(i,j-1,k,l)-f(i,j,k+1,l)-g(i,j,k-1,l)-ss(i,j,k,l)-acc(i,j,k,l)
            end do
        end do
        end do
        end do
    end subroutine
end module
