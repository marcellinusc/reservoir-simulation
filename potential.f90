module potential
    use functions
    use import_data
    use initial_conditions

    implicit none

    integer,allocatable,dimension(:,:,:)    :: ibw,icw,idw,iew
    integer,allocatable,dimension(:,:,:)    :: ifw,igw
    integer,allocatable,dimension(:,:,:)    :: ifo,igo

    contains

    subroutine flag()
        implicit none

        integer                             :: i,j,k
        real(kind=8)                        :: Pave,dz

        if(.not.allocated(ibw))allocate(ibw(nx,ny,nz))
        if(.not.allocated(icw))allocate(icw(nx,ny,nz))
        if(.not.allocated(idw))allocate(idw(nx,ny,nz))
        if(.not.allocated(iew))allocate(iew(nx,ny,nz))
        if(.not.allocated(ifw))allocate(ifw(nx,ny,nz))
        if(.not.allocated(igw))allocate(igw(nx,ny,nz))
        if(.not.allocated(ifo))allocate(ifo(nx,ny,nz))
        if(.not.allocated(igo))allocate(igo(nx,ny,nz))


        ibw = 0
        icw = 0
        idw = 0
        iew = 0
        ifw = 0
        igw = 0
        ifo = 0
        igo = 0

        dz  = tz/nz

        do i=1,nx
        do j=1,ny
        do k=1,nz
            !node B
            if(i>1) then
                if(Pk(i-1,j,k)>Pk(i,j,k)) ibw(i,j,k)=1
                !since oil potential flow behaves same with water does as long in x & y direction, ibo=ibw and so on
            end if
            !node C
            if(i<nx) then
                if(Pk(i+1,j,k)>Pk(i,j,k)) icw(i,j,k)=1
            end if
            !node D
            if(j>1) then
                if(Pk(i,j-1,k)>Pk(i,j,k)) idw(i,j,k)=1
            end if
            !node E
            if(j<ny) then
                if(Pk(i,j+1,k)>Pk(i,j,k)) iew(i,j,k)=1
            end if
            !node F
            if(k>1) then
                Pave = (Pk(i,j,k-1)+Pk(i,j,k))/2
                if(Pk(i,j,k-1)+dnw(Pave)*dz>Pk(i,j,k)) ifw(i,j,k)=1
                if(Pk(i,j,k-1)+dno(Pave)*dz>Pk(i,j,k)) ifo(i,j,k)=1
            end if
            !node G
            if(k<nz) then
                Pave = (Pk(i,j,k+1)+Pk(i,j,k))/2
                if(Pk(i,j,k+1)-dnw(Pave)*dz>Pk(i,j,k)) igw(i,j,k)=1
                if(Pk(i,j,k+1)-dno(Pave)*dz>Pk(i,j,k)) igo(i,j,k)=1
            end if
        end do
        end do
        end do
    end subroutine
end module
