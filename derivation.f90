module derivation
    use import_data
    use functions
    use initial_conditions
    use potential

    contains

    subroutine deriv(Soa,Son,Pa,Pn,delz,ijkw,ijko,PGEO,dT,FwTr,FoTr)
        !Soa    : oil saturation at the block a
        !Son    : oil saturation at the neighboring block
        !Pa     : pressure at the block a
        !Pn     : pressure at the neighboring block
        !delz   : depth difference between a and neighboring blocks
        !ijkw   : up-stream indicator for water phase
        !ijko   : up-stream indicator for oil phase
        !PGEO   : pseudogeometrical factor
        !dT(i)  : derivative of transmissibility term
                  !dT(1)=dFo/dSw,
                  !dT(2)=dFo/dP,
                  !dT(3)=dFw/dSw,
                  !dT(4)=dFw/dP
        !FoTr   : transmissibility term for water
        !FwTr   : transmissibility term for oil
        implicit none
        integer                     :: ijkw,ijko !to substitute potential of water and oil respectively
        real(kind=8)                :: Pave,Soa,Son,Pa,Pn,PGEO,krok,dkrok,krwk,dkrwk,delz
        real(kind=8)                :: Tro,Trw,FoTr,FwTr,dTrodSwn,dTrwdSwn,dTrodPn,dTrwdPn
        real(kind=8), allocatable, dimension(:)   :: dT

        if(.not. allocated(dT)) allocate(dT(4))

        Pave=(Pa+Pn)/2

        if(ijko==1) then
            krok  = kro(Son)
            dkrok = dkrodSw(Son)
            else
                krok  = kro(Soa)
                dkrok = 0.
        end if

        if(ijkw==1) then
            krwk  = krw(Son)
            dkrwk = dkrwdSw(Son)
            else
                krwk  = krw(Soa)
                dkrwk = 0.
        end if

        Tro  = krok/(uo(Pave)*Bo(Pave))*PGEO
        Trw  = krwk/(uw(Pave)*Bw(Pave))*PGEO
        FoTr = Tro*(Pn-Pa)-Tro*dno(Pave)*delz
        FwTr = Trw*(Pn-Pa)-Trw*dnw(Pave)*delz

        dTrodSwn = PGEO/(uo(Pave)*Bo(Pave))*dkrok
        dTrwdSwn = PGEO/(uw(Pave)*Bw(Pave))*dkrwk
        dTrodPn  = -Tro/2*(duodP(Pave)/uo(Pave)+dBodP(Pave)/Bo(Pave))
        dTrwdPn  = -Trw/2*(duwdP(Pave)/uw(Pave)+dBwdP(Pave)/Bw(Pave))

        dT(1) = dTrodSwn*((Pn-Pa)-dno(Pave)*delz)
        dT(2) = dTrodPn*(Pn-Pa)+Tro-(dTrodPn*dno(Pave)+Tro*ddnodP(Pave)/2)*delz
        dT(3) = dTrwdSwn*((Pn-Pa)-dnw(Pave)*delz)
        dT(4) = dTrwdPn*(Pn-Pa)+Trw-(dTrwdPn*dnw(Pave)+Trw*ddnwdP(Pave)/2)*delz
    end subroutine
end module
