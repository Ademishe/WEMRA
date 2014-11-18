module geometry_and_data
use array_work
use com_prog
use omp_lib
implicit none
contains

subroutine indat
    integer is
    real*8 z(sk)

    if(ss00.eq.1) then
        call geom1
    else
        call geom2
    end if
    ! print *, rt(:)
    call parameters1
    call gazel1
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
  integer ibeam
  real*8 delta_alpha

  ce=(0.0d0,1.0d0)
  ee=(1.0d0,0.0d0)
  pi = 3.14159265359d0
  gam0= w0/3.0d0
  constq = 2.0d0*beam_curr/mk0/(w0**2)*2.0d0*pi  !для уравн.возбуждения; здесь не должно быть 2pi
  const1 = -9.0d0/(511.0d0*1000.0d0*w0) !для уравн. движения

  delta_alpha = 2.0d0*pi / nbeam
  do ibeam = 1, nbeam
    alpha(ibeam) = (ibeam - 1) * delta_alpha + alpha_shift
    rb(ibeam) = rb0 / (1.0d0 + ellips*abs(cos(alpha(ibeam))))
  end do
  print *, rb(:)
  print *, alpha(:)

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
  integer ia,i_r,i_alpha, i, index
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
  integer is,inr,ina
  do is=1,sk
    do  inr = 1, nkr
      do ina = 0, nka
        gam(is,inr,ina)=conjg(sqrt((gam0**2-(mu(inr,ina)/rt(is))**2)*ee))
        zn(is,inr,ina)=376.7d0*(gam(is,inr,ina))/gam0
      end do
    end do
  end do
  return
end subroutine


subroutine eid
  integer inr, is, ina, inbm
  complex*16 E, integral, temp_int, a, b
  integer thread_num, thread_id
  real*8 accuracy, chi
  accuracy = 0.0001d0
!!$omp parallel default(private) shared(eznbm, ernbm, ephinbm, rt, gam, zn, rb, mu, alpha, accuracy, sk, nbeam, nka, nkr, ce, ee)
!!$omp do
  do is = 1, sk
    ! thread_id = omp_get_thread_num()
    ! thread_num = omp_get_num_threads()
    ! print *, "***eid***",  is, sk!, thread_id, thread_num
    do inbm = 1, nbeam
      do ina = 0, nka
        do inr = 1, nkr
          temp_int = 0.0d0
          chi = mu(inr, ina)/rt(is)
          if (ina.ne.0) then
            temp_int = (ina*ina/chi/chi) * integrate(bessel_mult, 0.0d0, rt(is), accuracy) * integrate(der_cossin_mult, 0.0d0, 6.2831853d0, accuracy)
          end if
          a = integrate(der_bessel_mult, 0.0d0, rt(is), accuracy)
          b = integrate(cossin_mult, 0.0d0, 6.2831853d0, accuracy)
          integral = a * b + temp_int
          E = sqrt(zn(is, inr, ina) / abs(zn(is, inr, ina)) * chi * chi * conjg(zn(is, inr, ina)) / gam(is, inr, ina) / conjg(gam(is, inr, ina)) / integral)
          ! print *, E, integral, a, b, temp_int
          eznbm(is,inr,ina,inbm) = bessel_jn(ina, chi*rb(inbm)) * cos(ina*alpha(inbm)) * E
          ernbm(is,inr,ina,inbm) = -ce * gam(is, inr, ina) / chi * 0.5d0 * (bessel_jn(ina-1, chi*rb(inbm)) - bessel_jn(ina+1, chi*rb(inbm))) * cos(ina*alpha(inbm)) * E
          ephinbm(is,inr,ina,inbm) = -ce * gam(is, inr, ina) / chi / chi * ina / rb(inbm) * bessel_jn(ina, chi*rb(inbm)) * sin(ina*alpha(inbm)) * E
        end do
      end do
    end do
  end do
!!$omp end do
!!$omp end parallel
  return
  contains

  real*8 function bessel_mult(r)
    real*8 r
    bessel_mult = bessel_jn(ina, chi*r)**2 / r
    return
  end function

  real*8 function der_bessel_mult(r)
    real*8 r
    der_bessel_mult = (bessel_jn(ina-1, chi*r) - bessel_jn(ina+1, chi*r))**2 / 4.0d0 * r
    return
  end function

  real*8 function cossin_mult(r)
    real*8 r
    cossin_mult = cos(ina*r)**2
    return
  end function

  real*8 function der_cossin_mult(r)
    real*8 r
    der_cossin_mult = sin(ina*r)**2
    return
  end function
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
