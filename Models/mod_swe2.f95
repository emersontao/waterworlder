MODULE swe2025 

    use mod_fint
    
    implicit none 

    integer,parameter :: Nwmf=12,Nirr=3,Nmass=20,Nage=17
    INTEGER, parameter :: ndim=4     !Number of dimensions
    INTEGER, parameter :: swe_ng(ndim)=(/Nwmf,Nirr,Nmass,Nage/) !Number of points along each dimension (i.e., [nx, ny, ...])

    double precision,dimension(Nwmf,Nirr,Nmass,Nage) :: swe_radius,swe_errcode
    double precision,dimension(Nwmf+Nirr+Nmass+Nage) :: swe_grid = (/ &
    0.001d0,0.01d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, &
    0.7d0, 0.8d0, 0.9d0, 1.0D0, & ! 0 to 0.9 with step 0.1
    500.0d0, 600.0d0, 700.0d0, & ! 500 to 700 with step 100
    0.2d0       ,  0.254855d0  ,  0.32475535d0,  0.41382762d0,  0.52733018d0,&
    0.67196366d0,  0.85626648d0,  1.09111896d0,  1.39038559d0,  1.77173358d0,&
    2.25767578d0,  2.87689978d0,  3.66596142d0,  4.67144294d0,  5.95270288d0,&
    7.58538038d0,  9.66586048d0, 12.31696422d0, 15.69519941d0, 20.d0,&
    0.001d0,0.0015d0,0.002d0,0.003d0,0.005d0,0.01d0,&
                        0.02d0,0.03d0,0.05d0,&
                        0.1d0,0.2d0,0.5d0,&
                        1.0d0,2.0d0,5.0d0,&
                        10d0,20d0         /) ! Final log-spaced values
    real(8) :: swe_axis1(Nwmf)=(/0.001d0,0.01d0, 0.1d0, 0.2d0, 0.3d0, &
    0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, 0.9d0, 1.0D0/), &
    swe_axis2(Nirr)=(/500.0d0, 600.0d0, 700.0d0/), &
    swe_axis3(Nmass)=(/0.2d0       ,  0.254855d0  ,  0.32475535d0,  &
    0.41382762d0,  0.52733018d0,&
    0.67196366d0,  0.85626648d0,  1.09111896d0,  1.39038559d0,  1.77173358d0,&
    2.25767578d0,  2.87689978d0,  3.66596142d0,  4.67144294d0,  5.95270288d0,&
    7.58538038d0,  9.66586048d0, 12.31696422d0, 15.69519941d0, 20.d0/), &
    swe_axis4(Nage)=(/0.001d0,0.0015d0,0.002d0,0.003d0,0.005d0,0.01d0,&
                        0.02d0,0.03d0,0.05d0,&
                        0.1d0,0.2d0,0.5d0,&
                        1.0d0,2.0d0,5.0d0,&
                        10d0,20d0/)

    contains

    subroutine load_data(fname,top_level)
        integer :: i,j,k,l
        integer,intent(in) :: top_level !0 for 20 mbar, 1 for 1 microbar
        real(8) :: ra,trash
        character(len=4) :: trash_char
        character(len=*), intent(in) :: fname

        open(10,file=fname,status='old',action='read')
        ! Header
        do i=1,40
           read(10,*)
        end do

        ! Data
        do i=1,Nwmf
            do j=1,Nirr
                do k=1,Nmass
                    do l=1,Nage
                        if(top_level.eq.0) then
                            read(10,*) trash_char,trash,trash,trash,trash, &
                                    trash,ra,trash,trash, trash, &
                                    trash,trash,trash,trash,trash, &
                                    trash,trash
                        else if (top_level.eq.1) then
                            read(10,*) trash_char,trash,trash,trash,trash, &
                                    trash,trash,ra,trash,trash, &
                                    trash,trash,trash,trash,trash, &
                                    trash,trash
                        else
                            write(6,*) "Error in SWE2025 interpolator."
                            write(6,*) "Wrong top level:", top_level 
                            stop
                        end if
                        swe_radius(i,j,k,l) = ra
                        ! print*, trash_char,trash,ra
                        ! stop
                        !iop_errcode(i,j,k,l) = errcode
                    end do
                read(10,*)
                read(10,*)
                end do
            end do
        end do
        close(10)

    end subroutine load_data

    function interpol_swe(wmf,tirr,mp,age) result(rp)
        implicit NONE
    
        real(8) :: rp
        !REAL(8) :: x(ndim)  !The interpolation point
        !REAL(8) :: grid(SUM(ng)) !Coordinates, one dimension after the other (i.e., [x(:), y(:), ...])
        !REAL(8),allocatable :: values(:,:,:,:) !Tabulated values, linearized in column major format
    
        real(8),intent(in) :: wmf,tirr,mp,age
        !real(8) :: logg
        !real(8) :: z20,z01,matm,tsurf
    
        !allocate(values(ng(1),ng(2),ng(3),ng(4)))
    
        !open()
        rp = fint(ndim,(/wmf,tirr,mp,age/),swe_ng,sum(swe_ng),product(swe_ng),swe_grid,swe_radius)
    end function interpol_swe

end module swe2025