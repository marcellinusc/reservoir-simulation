module functions
    use table_lookup
    use import_data

    contains

    !==============================OIL==============================
    function Rs(P) !volume in MSCF unit
        implicit none
        real(kind=8)    :: Rstemp,Rs,P
        call lookup(tabPo,tabRs,npvto,P,Rstemp)
        Rs      = Rstemp/5.6146d-3
        return
    end function

    function Bo(P)
        implicit none
        real(kind=8)    :: Bo,P
        call lookup(tabPo,tabBo,npvto,P,Bo)
        return
    end function

    function uo(P)
        implicit none
        real(kind=8)    :: uo,P
        call lookup(tabPo,tabuo,npvto,P,uo)
        return
    end function

    function rhoo(P)
        implicit none
        real(kind=8)    :: rhoo,P
        rhoo    = (rhoosc+Rs(P)*rhogsc)/Bo(P)
    end function

    function dno(P) !volume in cuft unit
        implicit none
        real(kind=8)    :: dno,P
        dno     = rhoo(P)/144d0
        return
    end function

    function dRsdP(P)
        implicit none
        real(kind=8):: dRsdP,P,xx,y,yy
        y       = Rs(P)
        xx      = P-0.001
        yy      = Rs(xx)
        dRsdP   = (y-yy)/(P-xx)
        return
    end function

    function dBodP(P)
        implicit none
        real(kind=8)    :: dBodP,P,xx,y,yy
        y       = Bo(P)
        xx      = P-0.001
        yy      = Bo(xx)
        dBodP   = (y-yy)/(P-xx)
        return
    end function

    function duodP(P)
        implicit none
        real(kind=8):: duodP,P,xx,y,yy
        y       = uo(P)
        xx      = P-0.001
        yy      = uo(xx)
        duodP   = (y-yy)/(P-xx)
        return
    end function

    function drhoodP(P)
        implicit none
        real(kind=8)  :: drhoodP,P,xx,y,yy
        y       = rhoo(P)
        xx      = P-0.001
        yy      = rhoo(xx)
        drhoodP = (y-yy)/(P-xx)
        return
    end function

    function ddnodP(P)
        implicit none
        real(kind=8)    :: ddnodP,P,xx,y,yy
        y       = dno(P)
        xx      = P-0.001
        yy      = dno(xx)
        ddnodP  = (y-yy)/(P-xx)
        return
    end function

    !==============================GAS==============================
    function Bg(P) !volume in MSCF unit
        implicit none
        real(kind=8)    :: Bgtemp,Bg,P
        call lookup(tabPg,tabBg,npvtg,P,Bgtemp)
        Bg      = Bgtemp*5.6146d-3
        return
    end function

    function ug(P)
        implicit none
        real(kind=8)    :: ug,P
        call lookup(tabPg,tabug,npvtg,P,ug)
        return
    end function

    function rhog(P)
        implicit none
        real(kind=8)    :: rhog,P
        rhog    = rhogsc/Bg(P)
    end function

    function dng(P) !volume in cuft unit
        implicit none
        real(kind=8)    :: dng,P
        dng     = rhog(P)/144d0
        return
    end function

    function dBgdP(P)
        implicit none
        real(kind=8)    :: dBgdP,P,xx,y,yy
        y       = Bg(P)
        xx      = P-0.001
        yy      = Bg(xx)
        dBgdP   = (y-yy)/(P-xx)
        return
    end function

    function dugdP(P)
        implicit none
        real(kind=8):: dugdP,P,xx,y,yy
        y       = ug(P)
        xx      = P-0.001
        yy      = ug(xx)
        dugdP   = (y-yy)/(P-xx)
        return
    end function

    function drhogdP(P)
        implicit none
        real(kind=8)  :: drhogdP,P,xx,y,yy
        y       = rhog(P)
        xx      = P-0.001
        yy      = rhog(xx)
        drhogdP = (y-yy)/(P-xx)
        return
    end function

    !=============================WATER============================
    function Bw(P)
        implicit none
        real(kind=8)   :: X,P,Bw
        X       = cw*(P-Pwref)
        Bw      = Bwref*(1-X)
        return
    end function

    function uw(P)
        implicit none
        real(kind=8)   :: Y,P,uw
        Y       = cuw*(P-Pwref)
        uw      = uwref+Y
        return
    end function

    function rhow(P)
        implicit none
        real(kind=8)    :: rhow,P
        rhow    = rhowsc/Bw(P)
    end function

    function dnw(P) !volume in cubic ft unit
        implicit none
        real(kind=8)    :: dnw,P
        dnw     = rhow(P)/144d0
        return
    end function

    function dBwdP(P)
        implicit none
        real(kind=8)    :: dBwdP,P,xx,y,yy
        y       = Bw(P)
        xx      = P-0.001
        yy      = Bw(xx)
        dBwdP   = (y-yy)/(P-xx)
        return
    end function

    function duwdP(P)
        implicit none
        real(kind=8)    :: duwdP,P,xx,y,yy
        y       = uw(P)
        xx      = P-0.001
        yy      = uw(xx)
        duwdP   = (y-yy)/(P-xx)
        return
    end function

    function drhowdP(P)
        implicit none
        real(kind=8)    :: drhowdP,P,xx,y,yy
        y       = rhow(P)
        xx      = P-0.001
        yy      = rhow(xx)
        drhowdP = (y-yy)/(P-xx)
        return
    end function

    function ddnwdP(P)
        implicit none
        real(kind=8)    :: ddnwdP,P,xx,y,yy
        y       = dnw(P)
        xx      = P-0.001
        yy      = dnw(xx)
        ddnwdP  = (y-yy)/(P-xx)
        return
    end function

    !=============================ROCK=============================
    function phi(P)
        implicit none
        real(kind=8)    :: phi,P
        phi     = phii*exp(cr*(P-Pref))
    end function

    function kro(Sw)
        implicit none
        real(kind=8)    :: kro,Sw
        call lookup(tabSw,tabkro,nrl,Sw,kro)
        return
    end function

    function krw(Sw)
        implicit none
        real(kind=8)    :: krw,Sw
        call lookup(tabSw,tabkrw,nrl,Sw,krw)
        return
    end function

    function Pcow(Sw)
        implicit none
        real(kind=8)    :: Pcow,Sw
        call lookup(tabSw,tabPcow,nrl,Sw,Pcow)
        return
    end function

    function dphidP(P)
        implicit none
        real(kind=8)  :: dphidP,P,xx,y,yy
        y       = phi(P)
        xx      = P-0.001
        yy      = phi(xx)
        dphidP  = (y-yy)/(P-xx)
        return
    end function

    function dkrodSw(Sw)
        implicit none
        real(kind=8)  :: dkrodSw,Sw,xx,y,yy
        y       = kro(Sw)
        xx      = Sw-0.001
        yy      = kro(xx)
        dkrodSw = (y-yy)/(Sw-xx)
        return
    end function

    function dkrwdSw(Sw)
        implicit none
        real(kind=8)  :: dkrwdSw,Sw,xx,y,yy
        y       = krw(Sw)
        xx      = Sw-0.001
        yy      = krw(xx)
        dkrwdSw = (y-yy)/(Sw-xx)
        return
    end function

    function dPcowdSw(Sw)
        implicit none
        real(kind=8)  :: dPcowdSw,Sw,xx,y,yy
        y       = Pcow(Sw)
        xx      = Sw-0.001
        yy      = Pcow(xx)
        dPcowdSw= (y-yy)/(Sw-xx)
        return
    end function
end module
