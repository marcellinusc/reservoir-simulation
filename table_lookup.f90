module table_lookup
    contains

    subroutine lookup(tabx, taby, n, x, y)
        implicit none
        integer        :: i,n
        real(kind=8)   :: x,y
        real(kind=8), dimension(n)   :: tabx, taby

        i=1
        do while(x>tabx(i) .AND. i<=n)
            i=i+1
        end do
        y=taby(i-1)+(taby(i)-taby(i-1))*(x-tabx(i-1))/(tabx(i)-tabx(i-1))
    end subroutine
end module
