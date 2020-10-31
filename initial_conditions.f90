module initial_conditions
    use functions
    use table_lookup
    use import_data

    implicit none
    real(kind=8)                              :: OOIP,OWIP
    real(kind=8)                              :: PGEOx,PGEOy,PGEOz
    real(kind=8)                              :: Vb
    real(kind=8)                              :: remoil,remwat
    real(kind=8),allocatable,dimension(:,:,:) :: Pk,Sok,Swk
    !k means at exact the interface in between vertical (z-direction) grids

    contains

    subroutine IC()
        implicit none
        integer                    :: i,j,k
        real(kind=8)               :: dx,dy,dz
        real(kind=8)               :: Pave
        real(kind=8),dimension(nz) :: Pia,Pi

        if(.not.allocated(Pk))     allocate(Pk(nx,ny,nz))
        if(.not.allocated(Sok))    allocate(Sok(nx,ny,nz))
        if(.not.allocated(Swk))    allocate(Swk(nx,ny,nz))

        dx = tx/nx
        dy = ty/ny
        dz = tz/nz

        Vb = dx*dy*dz

        PGEOx=6.3283d-3*dy*dz*kx/dx
        PGEOy=6.3283d-3*dx*dz*ky/dy
        PGEOz=6.3283d-3*dx*dy*kz/dz

        do i=1,nx
            do j=1,ny
                do k=1,nz
                    Swk(i,j,k) = Swi
                    if(k==1) then
                        Pi(k)      = Pinit
                        Pk(i,j,k) = Pi(k)
                    else
                        Pia(k)=Pi(k-1)
                        do
                            Pave       = 0.5d0*(Pia(k)+Pi(k-1))
                            Pi(k)      = Pi(k-1)+dno(Pave)*dz
                            Pk(i,j,k)  = Pi(k)
                            if((Pi(k)-Pia(k))/Pi(k)<0.1d-5) then
                                exit
                            else
                                Pia(k) = Pi(k)
                            end if
                        end do
                    end if
                end do
            end do
        end do

        OOIP = 0d0
        OWIP = 0d0

        do i=1,nx
        do j=1,ny
        do k=1,nz
            OOIP = OOIP+(Vb*phi(Pk(i,j,k))*(1-Swk(i,j,k))/Bo(Pk(i,j,k)))
            OWIP = OWIP+(Vb*phi(Pk(i,j,k))*Swk(i,j,k)/Bw(Pk(i,j,k)))
            !OGIP=OGIP+Vb*phi*(Rs*So/Bo+Sg/Bg)
            !since So=1-Sw at initial condition, Sg=0, then Sg*Bg=0d0
        end do
        end do
        end do
    end subroutine

    subroutine remaining()
        implicit none
        integer     :: i,j,k

        remoil = 0d0
        remwat = 0d0

        do i=1,nx
        do j=1,ny
        do k=1,nz
            remoil = remoil+(Vb*phi(Pk(i,j,k))*(1-Swk(i,j,k))/Bo(Pk(i,j,k)))
            remwat = remwat+(Vb*phi(Pk(i,j,k))*Swk(i,j,k)/Bw(Pk(i,j,k)))
            !OGIP=OGIP+Vb*phi*(Rs*So/Bo+Sg/Bg)
            !since So=1-Sw at initial condition, Sg=0, then Sg*Bg=0d0
        end do
        end do
        end do
    end subroutine
end module
