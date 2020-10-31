module import_data
    implicit none
    integer                                 :: npvto,npvtg,nrl
    integer                                 :: nx,ny,nz
    real(kind=8)                            :: tx,ty,tz
    real(kind=8)                            :: Pinit,Swi
    real(kind=8)                            :: phii,cr,Pref
    real(kind=8)                            :: kx,ky,kz
    real(kind=8)                            :: rhoosc,rhogsc,rhowsc
    real(kind=8)                            :: Pwref,Bwref,cw,uwref,cuw
    real(kind=8),allocatable,dimension(:)   :: tabRs,tabPo,tabBo,tabuo
    real(kind=8),allocatable,dimension(:)   :: tabSw,tabkrw,tabkro,tabPcow
    real(kind=8),allocatable,dimension(:)   :: tabPg,tabBg,tabug
    integer                                 :: nwell
    integer                                 :: nratetemp
    integer,allocatable,dimension(:)        :: nrate
    integer,allocatable,dimension(:,:)      :: wloc
    real(kind=8),allocatable,dimension(:,:) :: tabt,tabqt
    real(kind=8)                            :: check

    contains

    subroutine fetch()
        implicit none
        integer     :: i
        integer                                 :: iw,ir

        open(1, file='Data.txt')

        call comment(5)
        read(1,*) nx,ny,nz
        call comment(5)
        read(1,*) tx,ty,tz
        call comment(5)
        read(1,*) Pinit,Swi
        call comment(6)
        read(1,*) phii,cr,Pref
        call comment(5)
        read(1,*) kx,ky,kz
        call comment(4)
        read(1,*) nrl
        allocate (tabSw(nrl),tabkrw(nrl),tabkro(nrl),tabPcow(nrl))
        call comment(6)
        do i=1,nrl
            read(1,*) tabSw(i),tabkrw(i),tabkro(i),tabPcow(i)
        end do
        call comment(5)
        read(1,*) rhoosc,rhogsc,rhowsc
        call comment(5)
        read(1,*) npvto,npvtg
        allocate(tabRs(npvto),tabPo(npvto),tabBo(npvto),tabuo(npvto))
        allocate(tabPg(npvtg),tabBg(npvtg),tabug(npvtg))
        call comment(7)
        do i=1,npvto
            read(1,*) tabRs(i),tabPo(i),tabBo(i),tabuo(i)
        end do
        call comment(6)
        do i=1,npvtg
            read(1,*) tabPg(i),tabBg(i),tabug(i)
        end do
        call comment(8)
        read(1,*) Pwref,Bwref,cw,uwref,cuw
        call comment(4)
        read(1,*) nwell
        allocate(nrate(nwell))
        allocate(wloc(nwell,3))
        call comment(4)
        do iw=1,nwell
            read(1,*) wloc(iw,1),wloc(iw,2),wloc(iw,3) !x,y,z coordinate
        end do
        call comment(4)
        read(1,*) nratetemp
        allocate(tabt(nwell,nratetemp),tabqt(nwell,nratetemp))
        do iw=1,nwell
            nrate(iw)=nratetemp
        end do
        do iw=1,nwell
            call comment(5)
            do ir=1,nrate(iw) !1=injector 2=producer
                read(1,*) tabt(iw,ir),tabqt(iw,ir)
            end do
        end do
        !read(1,*) check
        close(1)
    end subroutine

    subroutine comment(n)
        implicit none
        integer     :: i,n
        character   :: temp
        do i=1,n
            read(1,*) temp
        end do
    end subroutine
end module
