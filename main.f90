Program res_simulation
	use import_data
	use table_lookup
	use functions
	use initial_conditions
	use potential
	use well
	use derivation
	use jacobian

    implicit none
    integer                                     :: ndn,ndm,ndgb,itime,iw,iter,itermax,i,j,k,l
    real(kind=8)                                :: tmax,tnew,delp_max,delsw_max,fw_max,fo_max,cum_prod_wat,cum_prod_oil,cum_inj_wat
    real(kind=8)                                :: eps_fw,eps_fo,eps_s,eps_p,dplim,dslim,dt_old
    real(kind=8)                                :: qwat,qoil,qinj,WC,WOR,err_wat,err_oil,delpn_max,delsn_max,dt_s,dt_p
    real(kind=8),allocatable,dimension(:)       :: rhs,del,gb
    real(kind=8),allocatable,dimension(:,:,:)   :: delsw,delp

    call fetch()

    ndn  = 2*nx*ny*nz
    ndm  = 2*ny*nz+1
    ndgb = ndn*(2*ndm+1)-ndm*(ndm+1)

    eps_Fw  = 5d0; eps_Fo = 1d0; eps_s = 0.001d0; eps_p = 0.1d0
    dplim   = 50.; dslim = 0.02d0
    itermax = 15

    cum_prod_wat = 0
    cum_prod_oil = 0
    cum_inj_wat  = 0

    t    = 0.d0
    delt = 0.1d0
    tmax = 7300.

    allocate(gb(ndgb),del(ndn),rhs(ndn),delsw(nx,ny,nz),delp(nx,ny,nz))

    call IC()

    open(2,file="Fluid Reserve.txt")
        write(2,*) "Reservoir Fluid In Place"
        write(2,*) ""
        write(2,1) "OOIP = ",OOIP/5.6146d6," MMSTB"
        write(2,1) "OWIP = ",OWIP/5.6146d6," MMSTB"
        1 format(a10,f7.3,a7)
        write(2,*) ""
        write(2,*) "  Pk(3,3,1)                 Pk(3,3,2)                 Pk(3,3,3)                 Pk(3,3,4)                 Pk(3,&
                    &3,5)"
        write(2,*) (pk(3,3,k),k=1,nz)
    close(2)

    open(3,file="Simulation Result.txt")
    write(3,*) "Simulation Result by Fortran 90 Coding"
    write(3,*) ""
    write(3,*) "   i         t      dt       q_inj       q_oil       q_wat          WC         WOR             Wi             Np  &
                &           Wp     Error_w     Error_o     p_inj    p_prod"

    itime=0
    do
        t    = t+delt
        tnew = t
        if(tnew>tmax) then
            t    = tmax
            delt = delt-(tnew-t)
        end if
        pn  = pk
        swn = Swk

        100 continue
        call flag()

        iter = 0
        do
            iter = iter+1
            call rate()
            call jacob()
            call solver(nx,ny,nz,Fores,Fwres,a,b,c,d,e,f,g,gb,rhs,del)

            l=0
            do i=1,nx
            do j=1,ny
            do k=1,nz
                l=l+1
                delsw(i,j,k)=del(l)
                Swk(i,j,k)=Swk(i,j,k)+delsw(i,j,k)

                l=l+1
                delp(i,j,k)=del(l)
                Pk(i,j,k)=Pk(i,j,k)+delp(i,j,k)
            end do
            end do
            end do

            delp_max  = maxval(dabs(delp))
            delsw_max = maxval(dabs(delsw))
            fw_max    = maxval(dabs(FwRes))
            fo_max    = maxval(dabs(FoRes))

            if(fw_max<eps_fw.AND.fo_max<eps_fo.AND.delp_max<eps_p.AND.delsw_max<eps_s)then
                exit
            else if(iter>itermax)then
                delt = 0.5d0*delt
                t    = t-delt
                goto 100
            else
                cycle
            end if
        end do

        do iw=1,nwell
            if(qw(iw)>0.d0)then
                cum_prod_wat = cum_prod_wat+qw(iw)*delt
                cum_prod_oil = cum_prod_oil+qo(iw)*delt
            end if
            if(qw(iw)<0.d0)then
                cum_inj_wat  = cum_inj_wat+dabs(qw(iw))*delt
            end if
        end do

        call remaining()

        err_wat = (OWIP-remwat-cum_prod_wat+cum_inj_wat)/OWIP
        err_oil = (OOIP-remoil-cum_prod_oil)/OOIP

        do iw=1,nwell
            if(qw(iw)>0.d0)then
                qwat = qw(iw)
                qoil = qo(iw)
            end if
            if(qw(iw)<0.d0)then
                qinj = dabs(qw(iw))
            end if
        end do

        WC  = qwat/(qwat+qoil)
        WOR = qwat/qoil

        write (*,*) "Time-loop = ", itime
        write (*,*) "============================================="
        write (*,4) "PwbInj    = ", pk(wloc(1,1),wloc(1,2),wloc(1,3)), "psi    "
        write (*,4) "PwbProd   = ", pk(wloc(2,1),wloc(2,2),wloc(2,3)), "psi    "
        write (*,4) "Time      = ", t, "days   "
        write (*,4) "Time inc. = ", delt, "days   "
        write (*,4) "Wat rate  = ", qwat, "SCF/D  "
        write (*,4) "Oil rate  = ", qoil, "SCF/D  "
        write (*,4) "Inj rate  = ", qinj, "SCF/D  "
        write (*,4) "Cum. wat  = ", cum_prod_wat*1d-3, "MSCF/D "
        write (*,4) "Cum. oil  = ", cum_prod_oil*1d-3, "MSCF/D "
        write (*,4) "Cum. inj  = ", cum_inj_wat*1d-3, "MSCF/D "
        write (*,4) "Watercut  = ", WC*100, "%      "
        write (*,4) "WOR       = ", WOR, " "
        write (*,5) "Err. wat  = ", err_wat, " "
        write (*,5) "Err. oil  = ", err_oil, " "
        write (*,*) "============================================="
        4 format (a14,f10.2,a8)
        5 format (a14,e12.4,a8)

        write (3,2) itime, t, delt, qinj, qoil, qwat, WC, WOR, cum_inj_wat*1d-3, cum_prod_oil*1d-3, cum_prod_wat*1d-3, err_wat, &
        err_oil, pk(wloc(1,1),wloc(1,2),wloc(1,3)), pk(wloc(2,1),wloc(2,2),wloc(2,3))
        2 format (i5, f10.2, f8.1, 3f12.2, 2f12.6, 3e15.5, 2e12.2, 2f10.2)

        if (t>=tmax) then
            write(3,*) ""
            write(3,3) "RF = ",cum_prod_oil/OOIP*100,"%"
            write(*,*) "RF = ",cum_prod_oil/OOIP*100,"%"
            3 format(a5,f10.6,a1)
            close(3)
            stop
        else if (delt<=1d-6) then
            write(*,*) "Warning, not convergen"
            close(3)
            stop
        else
            delpn_max = maxval(dabs(Pk-pn))
            delsn_max = maxval(dabs(Swk-swn))
            dt_old    = delt
            dt_s      = dslim/delsn_max
            dt_p      = dplim/delpn_max
            delt      = delt*min(dt_s,dt_p)
            if(delt/dt_old > 2.d0) delt = dt_old*2.d0
            if(delt>30.d0)         delt = 30.d0
            itime = itime + 1
        end if
    end do
end program res_simulation
