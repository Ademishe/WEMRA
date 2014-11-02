module geometry_and_data
use array_work
use com_prog
implicit none

integer*4 ialpha
contains

subroutine indat
    integer is
    real*8 z(sk)

    if(ss00.eq.1) then
        call geom1
    else
        call geom2
    end if
    call parameters1
    call gazel1
    call eid
    dt=2.0d0*pi

    write (7,12) w0,v0,rel_factor
!      write (9,12) w0,v0,rel_factor,beam_voltage
12  format(2x,'omega0=',(f8.3,1x),2x,'vo =',f8.5,2x,'rel_factor=',f8.3,2x,'beam_voltage=',f8.5)

    zs(1)=v0*2.0d0*pi/w0     ! координата начала системы

    do is=1,sk        ! во всей программе is- номер регулярного участка
        zs(is+1)=zs(is)+dz(is)  ! координаты начал участков
        write (7,*) zs(is),rt(is)
    end do

    zss=zs(sk+1)              ! координата конца последнего участка
    zs(sk+2)=zss+3.0d0*2.0d0*pi/w0
    mkk=int(zss/(zs(1)/mk0))+1  ! число частиц во всей системе
    return
end subroutine


subroutine parameters1
    ce=(0.0d0,1.0d0)
    ee=(1.0d0,0.0d0)
    pi = 3.14159265359d0
    gam0= w0/3.0d0

    constq = 2.0d0*beam_curr/mk0/(w0**2)*2.0d0*pi  !для уравн.возбуждения

    const1 = -9.0d0/(511.0d0*1000.0d0*w0) !для уравн. движения
    rb(:)=rb0

    if ( kluch_beam.eq.0) then    ! холодная задача
        const1=0.0d0
        constq=0.0d0
    end if

    if ( kluch_beam.eq.2) then     ! заданный ток
        const1=0.0d0
    end if

    rel_factor = 1.0d0+beam_voltage/511.0d0
    v0 = 3.0d0*sqrt(1.0d0-1.0d0/(rel_factor**2))
    return
end subroutine


subroutine read_bessel_zeros
  integer ia,i_r,i_alpha
  open (32,file='bessel_zeros.cfg')
    do i_alpha=0,10
        do i_r=1,50
            read (32,*) ia, mu(i_r,i_alpha)
        end do
    end do
  close(32)
  print *, "Reading of Bessel equation zeros from bessel_zeros.cfg is complete"
  return
end subroutine


subroutine gazel1
! определение постоянных распространения, волновых сопротивлений
! и множителя для поля в  области пучка
  integer is,inr,ina
  do is=1,sk
    do ina=0,nka
      do  inr=1,nkr
        if (mu(inr,ina)/rt(is).le. gam0) nopen(is)=inr
        gam(is,inr,ina)=conjg(sqrt((gam0**2-(mu(inr,ina)/rt(is))**2)*ee))
        zn(is,inr,ina)=376.7d0*(gam(is,inr,ina))/gam0
      end do
    end do
  end do
  return
end subroutine


subroutine eid
    integer inr,is,ina,inbm
    do is=1,sk
        do ina=0,nka
            do  inr=1,nkr
                do  inbm=1,nbeam
                    eznbm(is,inr,ina,inbm)= ce*sqrt(abs(zn(is,inr,ina))/pi)* &
                      (mu(inr,ina)/(abs(gam(is,inr,ina))*(rt(is)**2)))* &
                      bessel_jn(0, mu(inr,ina)*rb(inbm)/rt(is))/bessel_jn(1, mu(inr,ina))
                end do
            end do
        end do
    end do
    return
end subroutine


subroutine geom1
! односекционная система
	integer is,iss
	integer sk01,sk00

! ***********************задание геометрии системы **********************
!           геометрия отражателя

	  do iss=1,sdop01                  ! участок гладкого волновода
	   dz(iss)                                = dz01
	   rt(iss)                                = rt01
        end do

	  do iss=1,sdop02                       ! неоднородность
	   dz(sdop01+iss)                          = dz02
         rt(sdop01+iss)                          = rt02
        end do
	 do iss=1,sdop03                       ! неоднородность
	   dz(sdop01+sdop02+iss)                   = dz03
         rt(sdop01+sdop02+iss)                   = rt03
        end do
         sk00=sdop01+sdop02+sdop03

!         геометрия секции

    if (kluch_structur(1).eq.1) then ! неоднородность в виде прямоугольника

        do is=1,sk1     ! цикл по периодам
            do iss=1,sdop11                  ! участок гладкого волновода
                dz((sdop11+sdop12)*(is-1)+sk00+iss)         = dz11
                rt((sdop11+sdop12)*(is-1)+sk00+iss)         = rt1
            end do

            do iss=1,sdop12                       ! неоднородность
                dz((sdop11+sdop12)*(is-1)+sdop11+sk00+iss)   = dz12
                rt((sdop11+sdop12)*(is-1)+sdop11+sk00+iss)   = rt2
            end do
        end do
    end if

    if (kluch_structur(1).eq.2) then
! неоднородность в виде полутора с пьедесталом
        do is=1,sk1
            do iss=1,sdop11
                dz((sdop11+sdop12)*(is-1)+sk00+iss)         = dz11
                rt((sdop11+sdop12)*(is-1)+sk00+iss)         = rt1
            end do

            do iss=1,sdop12
                dz((sdop11+sdop12)*(is-1)+sdop11+sk00+iss) = dz12
                rt((sdop11+sdop12)*(is-1)+sdop11+sk00+iss) = rt1-0.18d0-sqrt(0.325d0**2*(1.0d0-(5.5d0-1.0d0*iss)**2/25.0d0))
            end do
        end do

    end if
	sk01=sk1*(sdop11+sdop12)+sk00    ! номер последнего участка первой секции

    do iss=1,sdop11    ! участок гладкого волновода в конце системы
        dz(sk01+iss)            = dz11
        rt(sk01+iss)            = rt1
    end do

    periodz1 = sdop11*dz11+sdop12*dz12
    periodz2 = sdop31*dz31+sdop32*dz32
end subroutine


subroutine geom2
! двухсекционная система с трубой дрейфа

    integer is,iss
    integer sk00,sk01,sk02,sk03

! ***********************задание геометрии системы **********************
	  do iss=1,sdop01                  ! участок гладкого волновода
	   dz(iss)                                = dz01
	   rt(iss)                                = rt01
        end do

	  do iss=1,sdop02                       ! неоднородность
	   dz(sdop01+iss)                          = dz02
         rt(sdop01+iss)                          = rt02
        end do

	  do iss=1,sdop03                       ! участок гладкого волновода
	   dz(sdop01+sdop02+iss)                   = dz03
         rt(sdop01+sdop02+iss)                   = rt03
        end do

	   sk00=sdop01+sdop02+sdop03
! ************************  1-я секция **********************************

      if (kluch_structur(1).eq.1) then ! неоднородность в виде прямоугольника

	 do is=1,sk1     ! цикл по периодам
	  do iss=1,sdop11                  ! участок гладкого волновода
	   dz((sdop11+sdop12)*(is-1)+sk00+iss)         = dz11
	   rt((sdop11+sdop12)*(is-1)+sk00+iss)         = rt1
        end do

	  do iss=1,sdop12                       ! неоднородность
	   dz((sdop11+sdop12)*(is-1)+sdop11+sk00+iss)   = dz12
         rt((sdop11+sdop12)*(is-1)+sdop11+sk00+iss)   = rt2
        end do
	 end do

	end if

    if (kluch_structur(1).eq.2) then
! неоднородность в виде полутора с пьедесталом
        do is=1,sk1
            do iss=1,sdop11
                dz((sdop11+sdop12)*(is-1)+sk00+iss)         = dz11
                rt((sdop11+sdop12)*(is-1)+sk00+iss)         = rt1
            end do

            do iss=1,sdop12
                dz((sdop11+sdop12)*(is-1)+sdop11+sk00+iss) = dz12
                rt((sdop11+sdop12)*(is-1)+sdop11+sk00+iss) = rt1-0.28d0-sqrt(0.325d0**2*(1.0d0-(5.5d0-1.0d0*iss)**2/25.0d0))
            end do
        end do
    end if

	sk01=sk1*(sdop11+sdop12)+sk00    ! номер последнего участка первой секции

! ************************ труба дрейфа **********************************

	 do is=1,sk2     ! цикл по периодам
	  do iss=1,sdop21                  ! участок гладкого волновода
	   dz(sk01+(sdop21+sdop22)*(is-1)+iss)         = dz21
	   rt(sk01+(sdop21+sdop22)*(is-1)+iss)         = rt1
        end do

	  do iss=1,sdop22                    ! тоже участок гладкого волновода
	   dz(sk01+(sdop21+sdop22)*(is-1)+sdop21+iss)   = dz21
         rt(sk01+(sdop21+sdop22)*(is-1)+sdop21+iss)   = rt1
        end do
	 end do

	continue
	sk02=sk01+sk2*(sdop21+sdop22)    ! номер последнего участка трубы дрейфа

! ************************  2-я секция ***************************************

      if (kluch_structur(3).eq.1) then ! неоднородность в виде прямоугольника

	 do is=1,sk3     ! цикл по периодам
	  do iss=1,sdop31                  ! участок гладкого волновода
	   dz(sk02+(sdop31+sdop32)*(is-1)+iss)         = dz31
	   rt(sk02+(sdop31+sdop32)*(is-1)+iss)         = rt1
        end do

	  do iss=1,sdop32                       ! неоднородность
	   dz(sk02+(sdop31+sdop32)*(is-1)+sdop31+iss)   = dz32
         rt(sk02+(sdop31+sdop32)*(is-1)+sdop31+iss)   = rt2
        end do
	 end do

	end if

      if (kluch_structur(3).eq.2) then
! неоднородность в виде полутора с пьедесталом
	 do is=1,sk3
        do iss=1,sdop31
	   dz(sk02+(sdop31+sdop32)*(is-1)+iss)         = dz31
	   rt(sk02+(sdop31+sdop32)*(is-1)+iss)         = rt1
        end do

	  do iss=1,sdop32
	   dz(sk02+(sdop31+sdop32)*(is-1)+sdop31+iss)   = dz32
         rt(sk02+(sdop31+sdop32)*(is-1)+sdop31+iss)   = rt1-0.28d0-sqrt(0.325d0**2*(1.0d0-(5.5d0-1.0d0*iss)**2/25.0d0))
        end do
	 end do

	end if

    sk03=sk02+sk3*(sdop31+sdop32)

    do iss=1,sdop31    ! участок гладкого волновода в конце системы
        dz(sk03+iss)            = dz31
        rt(sk03+iss)            = rt1
    end do

    periodz1 =sdop11*dz11+sdop12*dz12
    periodz2 =sdop31*dz31+sdop32*dz32
end subroutine

end module geometry_and_data
