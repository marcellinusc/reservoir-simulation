module well
    use functions
    use import_data
    use initial_conditions

    implicit none

    real(kind=8)                            :: t
    real(kind=8),allocatable,dimension(:)   :: qo,qw
    real(kind=8),allocatable,dimension(:)   :: dqodP,dqwdP
    real(kind=8),allocatable,dimension(:)   :: dqodSw,dqwdSw

    contains

    subroutine rate()
        implicit none
        integer         :: iw,ir
        real(kind=8)    :: Pw,Sww
        real(kind=8)    :: qt
        real(kind=8)    :: M,fw

        if(.not.allocated(qo))     allocate(qo(nwell))
        if(.not.allocated(qw))     allocate(qw(nwell))
        if(.not.allocated(dqodP))  allocate(dqodP(nwell))
        if(.not.allocated(dqwdP))  allocate(dqwdP(nwell))
        if(.not.allocated(dqodSw)) allocate(dqodSw(nwell))
        if(.not.allocated(dqwdSw)) allocate(dqwdSw(nwell))

        do iw=1,nwell
            Pw=Pk(wloc(iw,1),wloc(iw,2),wloc(iw,3))
            Sww=Swk(wloc(iw,1),wloc(iw,2),wloc(iw,3))

            ir=1
            do while(t>tabt(iw,ir).and.ir<=nrate(nwell))
                ir=ir+1
            end do
            qt=tabqt(iw,ir)*5.6146d0

            if (kro(Sww) <= 1d-4) then
                qw(iw)      = qt
                qo(iw)      = 0.
                dqwdP(iw)   = 0.
                dqwdSw(iw)  = 0.
                dqodP(iw)   = 0.
                dqodSw(iw)  = 0.
            else if (krw(Sww) <= 1d-4) then
                qw(iw)      = 0.
                qo(iw)      = qt
                dqwdP(iw)   = 0.
                dqwdSw(iw)  = 0.
                dqodP(iw)   = 0.
                dqodSw(iw)  = 0.
            else
                if(iw==1) then
                    fw  = 1
                else
                    M   = uo(Pw)*Bo(Pw)*krw(Sww)/uw(Pw)/Bw(Pw)/kro(Sww)
                    fw  = M/(1+M)
                end if
                qw(iw)      = fw*qt
                qo(iw)      = (1-fw)*qt
                dqwdP(iw)   = qw(iw)*(1-fw)*(duodP(Pw)/uo(Pw)+dBodP(Pw)/Bo(Pw)-duwdP(Pw)/uw(Pw)-dBwdP(Pw)/Bw(Pw))
                dqwdSw(iw)  = qw(iw)*(1-fw)*(dkrwdSw(Sww)/krw(Sww)-dkrodSw(Sww)/kro(Sww))
                dqodP(iw)   = -dqwdP(iw)
                dqodSw(iw)  = -dqwdSw(iw)
            end if
        end do
    end subroutine
end module
