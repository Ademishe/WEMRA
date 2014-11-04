!****************************************************************************
!
!  program: nest_multy_mode_nes
!
!  purpose: main program
!
!****************************************************************************

program nest_multy_mode_nes
		use array_work
		use beam_nes
		use field_nes
		use geometry_and_data
		use com_prog

    implicit none

    integer*2  day, month, year
    integer*2  hour, minute, second, hund,skolko
		real*8     betazd
		complex*16 ampl_enter,ampl_exit

    character*11 pname,sname,fname,sumname,dann
    character*9 snamevh
    character*12 powename,fielname,specname
    character*1 e
    character*4 f

		integer iw,stepwrite,iomega,nmodealpha,nmoder


    character*1 string(60)

!    open(36,file='interes.dat')
    open(33,file='vremja.dat')
    open(35,file='phi.dat')
    open (3,file='multi5.cfg')
    open (4,file='amp.dat')
    open (5,file='spectrum.dat')
    open (9,file='iter.dat')
    open (7,file='sum.dat')
    open (21,file='power.dat')
    open (22,file='eta.dat')
    open (23,file='yarray.dat')
    open (24,file='field1.dat')
    open (25,file='fieldez2.dat')
    open (26,file='fielder2.dat')
!  	call gettim(hour, minute, second, hund)
!		write (36,*) hour, minute,second


! ***************  data read  ****************************************
    read (3,1000)                             string
    read (3,*)    ss00
    read (3,1000)                             string
    read (3,*)    kluch_beam,ktimemax,stepwrite
    read (3,1000)                             string
    read (3,*)    wh0,iw,dw,betazd
    read (3,1000)                             string
    read (3,*)    nkr, nmoder
    read (3,1000)                             string
    read (3,*)    nka, nmodealpha
    read (3,1000)                             string
    read (3,*)    ampl_enter,ampl_exit
    read (3,1000)                             string
    read (3,*)    beam_curr,beam_voltage,rb0,nbeam, mk0
    read (3,1000)                             string
    read (3,*)    ampl_mod, step_prmt, prmt_4

    read (3,1000)                             string
    read(3,*)     dz01,dz02,dz03
    read (3,1000)                             string
    read (3,*)    rt01,rt02,rt03
    read (3,1000)                             string
    read (3,*)    sdop01,sdop02, sdop03

    read (3,1000)                             string
    read (3,*)    sk1, kluch_structur(1)
    read (3,1000)                             string
    read(3,*)     rt1,rt2
    read (3,1000)                             string
    read (3,*)    dz11,dz12
    read (3,1000)                             string
    read (3,*)    sdop11,sdop12

    read (3,1000)                             string
    read (3,*)    sk2, kluch_structur(2)
    read (3,1000)                             string
    read (3,*)    dz21,dz22
    read (3,1000)                             string
    read (3,*)    sdop21,sdop22

    read (3,1000)                             string
    read (3,*)    sk3, kluch_structur(3)
    read (3,1000)                             string
    read (3,*)    dz31,dz32
    read (3,1000)                             string
    read (3,*)    sdop31,sdop32
		close(3)

    if (ss00.eq.1) then
        sk=sk1*(sdop11+sdop12)+sdop11+sdop01+sdop02+sdop03
    else
        sk=sdop01+sdop02+sdop03+sk1*(sdop11+sdop12)+sk2*(sdop21+sdop22)+sk3*(sdop31+sdop32)+sdop31
    end if
		print *, "sk = ", sk
    call configarr !выделение памяти для основных массивов
		call init_main_arrays
		call read_bessel_zeros

! амплитуды волн, приходящих в систему
    bplus0(nmoder,nmodealpha)		= ampl_enter
    bminussk(nmoder,nmodealpha)	= ampl_exit


!      sk - число регулярных участков, dz- период
!      nk-число мод, w0- рабочая частота
!      curn- ток потока, beam_voltage- напряжение
!      rb- радиус потока, mkk- число частиц на период
!      ktimemax - число итераций ( максимальное время расчета)
!      wh0,iw,dw -данные для спектрального анализа (начальная частота, число точек, шаг)
!      stepwrite -шаг выдачи траекторий и распределения полей

  do iomega=1,iw
		w0=wh0+(iomega-1)*dw
!			beam_voltage = (1.0d0/sqrt(1-((dz1*sdop1+dz2*sdop2)*w0/betazd/3.0d0)**2)-1.0d0)*511.0d0
		call indat            ! проведение вспомогательных расчетов geometry_and_data

		if (iomega.eq.1) then
			call configfield   ! выделение памяти для вспомогательных матриц
			call config_beam
		end if
		call init_arrays
		call matrix_construct ! определение статических матриц
		call manager(stepwrite,iw)
  end do

	call freeallmem

  close (4)
  close (5)
  close (9)
  close (7)
  close (21)
  close (22)
  close (23)
  close (24)
  close (25)
!    call gettim(hour, minute, second, hund)
!    write (36,*) hour, minute,second

    close(33)
    close(35)
!    close(36)

 277     format(i2,'f',i2)
 278     format(i2,'f',i2)
 1000    format(60(a1))

end program nest_multy_mode_nes


subroutine manager (stepwrite,iw)
!     управляющая процессом во времени

	use array_work
	use beam_nes
	use field_nes

  implicit none
  integer k,i11,is,in,stepwrite,iw,i,im,istwr
  real*8 tau,rel_coord

!     double complex fw(200)
!     double precision fff(200)
  complex fffw(256),fffb(256)


  dxnplus = 0.0d0
  dxnminus= 0.0d0
  xbplus(:,:,:)  = 0.0d0
  xbminus(:,:,:) = 0.0d0
  dxbplus = 0.0d0
  dxbminus= 0.0d0
  etaplus = 0.0d0
  etaminus= 0.0d0
  xbplus(1,:,:) = bplus0
  xnplus(1,:,:) = bplus0
  xnminus(sk,:,:)=bminussk
  xbminus(sk,:,:)=bminussk

	do inbeam=1,nbeam
  	call beam_begin
    k=0
    if(kluch_beam.eq.0)  then
    	do im=1,mkk(inbeam)
      	rel_coord=ddz*(im-1)-yarray(2*im)
        write (23,903) k,im, yarray(2*im),rel_coord,yarray(2*im-1),mkk,mk0
    	end do
    end if
  end do

!            main cycle
  tau=0.0d0
  istwr=0
  do k=1,ktimemax
!     то такое istwr? - это счетчик шагов при итерировании, нужен для записи на диск
  	istwr=istwr+1
    tau=tau+dt
    call matrix_dynamic
  	call progonka
    call field_calc

    i11=1
    if (k.ge.(ktimemax-128)) then
  		i11=1+k-(ktimemax-128)
      fffw(i11)=xbplus(sk,1,0)
      fffb(i11)=xbminus(1,1,0)
    end if
    etaplus = 0.0d0
    etaminus= 0.0d0
   !   xnplus  = xbplus
    !  xnminus = xbminus
  !    dxnplus = dxbplus
!	dxnminus= dxbminus
  	do inbeam=1,nbeam
  		if (kluch_beam.eq.2 .and. k.eq.1) then 	   ! режим заданного тока
    		call beam_calc (tau)
				print* ,'w0=',w0,'k=',k
      end if
      if(kluch_beam.eq.1)  then           ! самосогласованная задача
      	call beam_calc (tau)
      end if
          ! if(kluch_beam.eq.0)  then           ! холодная задача
! empty condition???
          ! end if
    end do
    if (istwr.eq.stepwrite) then
  		if (kluch_beam.eq.1) then
!            write (22,*) 'k=',k
!	  	  write(4,*) k
!		  write(4,*)'xbplus'
!		  write(4,*) xbplus
!            write(4,*)'xbminus'
!		  write(4,*) xbminus
      	do inbeam=1,nbeam
        	do im=1,mkk(inbeam)
          	rel_coord=ddz*(im-1)-yarray(2*im)
            write (23,903) k,im, yarray(2*im),rel_coord,yarray(2*im-1),mkk,mk0
          end do
        end do
      end if
      call field_structure(k,1)
      call field_structure(k,2)
      istwr=0
    end if
    if(kluch_beam.eq.1) then
    	call field_power(k,2)
    end if
  end do ! конец основного цикла

  if(kluch_beam.eq.0) then
  	call field_power(ktimemax,2)
    call field_power(ktimemax,1)
  else
    call field_power(ktimemax,1)
    call field_power(ktimemax,2)
    call field_structure (ktimemax,1)
    call field_structure (ktimemax,2)
  end if
  if (ktimemax.ge.128) then
    call spektr(fffw,fffb,w0,wh)
  end if

903 format(1x,'k=',i5,2x,'im=',i5,2x,e12.5,2x,e12.5,2x,e12.5,2x,'mkk=',i6,2x,'mk0=',i5)

	return
end

subroutine spektr(seq1,seq2,w0, w1)
    use fft_mod
    integer    n,i11,i
    parameter (n=128)

!      integer    i, nout
    real       twopi,w0,w1, fow(128),bcw(128),ww
    complex    c,cexp, coef(128), seq1(128),seq2(128)
    intrinsic  cexp

!
!      c     = (0.,1.)
!      twopi = 2.0*const('pi')
!                                  here we compute (2*pi*i/n)*3.
!      h = (twopi*c/n)*3.
!                                  this loop fills out the data vector
!                                  with a pure exponential signal of
!                                  frequency 3.

!                                  compute the fourier transform of seq
!      call fftcf (n, seq1, coef)
    call fft(coef)
!                                  get output unit number and print
!                                  results
    do 10  i=1, 64
        fow(i)=abs(coef(i+64))
        fow(i+64)=abs(coef(i))
10      continue
!        call fftcf (n, seq2, coef)
        call fft(coef)
        do 11  i=1, 64
            bcw(i) =abs(coef(i+64))
            bcw(i+64) =abs(coef(i))
11          continue
            call norma1(fow,n)
            call norma1(bcw,n)

            do 21 i11=1,128
                ww=real((i11-1)/128.0)+0.5
                write (5,921) ww,fow(i11),bcw(i11)
21              continue
921             format(1x,'w=',f9.4,2x,'fw=',e9.2,2x,'bw',e9.2)
!      call umach (2, nout)
!      write (nout,99998)
!99998 format (9x, 'index', 8x, 'seq', 15x, 'coef')
!      write (nout,99999) (i, seq(i), coef(i), i=1,n)
!99999 format (1x, i11, 5x,'(',f5.2,',',f5.2,')',
!     &                 5x,'(',f5.2,',',f5.2,')')
    return
end

subroutine norma1(a,n)
    real a(256),aa
    integer n,i
    aa=0.0
    do 1 i=1,n
        aa=max1(aa,a(i))
1       continue
        do 2 i=1,n
            a(i)=a(i)/aa
2           continue
    return
end

subroutine configarr
	use array_work
  allocate(nopen(sk))
  allocate(zs(sk+2),dz(sk),rt(sk),mu(50,0:10),nu(20),rb(nbeam), &
					alpha(nbeam),mkk(nbeam))
  allocate(gam(sk,nkr,0:nka),zn(sk,nkr,0:nka), &
					eznbm(sk,nkr,0:nka,nbeam),ernbm(sk,nkr,0:nka,nbeam), &
          xnplus(sk,nkr,0:nka),xnminus(sk,nkr,0:nka), &
          dxnplus(sk,nkr,0:nka),dxnminus(sk,nkr,0:nka), &
          xbplus(sk,nkr,0:nka),xbminus(sk,nkr,0:nka), &
          dxbplus(sk,nkr,0:nka),dxbminus(sk,nkr,0:nka), &
          bplus0(nkr,0:nka),bminussk(nkr,0:nka), &
          etaplus(sk,nkr,0:nka),etaminus(sk,nkr,0:nka))
  allocate(amplxplus(ktimemax),amplxminus(ktimemax))
  return
end subroutine configarr

subroutine init_main_arrays
	use array_work
	nopen(:)						= 0
	zs(:)								= 0.0d0
	dz(:)								= 0.0d0
	rt(:)								= 0.0d0
	mu(:,:)							= 0.0d0
	nu(:)								= 0.0d0
	rb(:)								= 0.0d0
	alpha(:)						= 0.0d0
	mkk(:)							= 0
	gam(:,:,:)					= (0.0d0, 0.0d0)
	zn(:,:,:)						= (0.0d0, 0.0d0)
	eznbm(:,:,:,:)			= (0.0d0, 0.0d0)
	ernbm(:,:,:,:)			= (0.0d0, 0.0d0)
	xnplus(:,:,:)				= (0.0d0, 0.0d0)
	xnminus(:,:,:)			= (0.0d0, 0.0d0)
	dxnplus(:,:,:)			= (0.0d0, 0.0d0)
	dxnminus(:,:,:)			= (0.0d0, 0.0d0)
	xbplus(:,:,:)				= (0.0d0, 0.0d0)
	xbminus(:,:,:)			= (0.0d0, 0.0d0)
	dxbplus(:,:,:)			= (0.0d0, 0.0d0)
	dxbminus(:,:,:)			= (0.0d0, 0.0d0)
	bplus0(:,:)					= (0.0d0, 0.0d0)
	bminussk(:,:)				= (0.0d0, 0.0d0)
	etaplus(:,:,:)			= (0.0d0, 0.0d0)
	etaminus(:,:,:)			= (0.0d0, 0.0d0)
	amplxplus(:)				= 0.0d0
	amplxminus(:)				= 0.0d0
end subroutine init_main_arrays


subroutine config_beam
!     распределение памяти для массивов,
!      используемых при интегрировании уравнений движения
    use array_work
    integer err2,dimyarray1,dimaux1,dimyarray2,dimyarray0

    dimyarray1=2*mkk(1)*nbeam+800   !+ns4+100
    dimyarray2=mkk(1)*nbeam+400
    dimyarray0=2*mkk(1)+200
    dimaux1=8

    allocate(yarray(dimyarray1),yarray1(dimyarray1), &
				yarray0(dimyarray0,nbeam), &
        dery(dimyarray1), &
        aux(dimaux1,dimyarray1), &
        velocity(dimyarray2), velocity1(dimyarray2))
    return
end subroutine


subroutine configfield
	use array_work

  allocate(aa1(sk,nkr,nkr,0:nka),aa2(sk,nkr,nkr,0:nka), &
            aa3(sk,nkr,nkr,0:nka), &
            alfaplus(sk,nkr,nkr,0:nka),alfaminus(sk,nkr,nkr,0:nka), &
            betaplus(sk,nkr,nkr,0:nka),betaminus(sk,nkr,nkr,0:nka), &
            d1plus(sk,nkr,nkr,0:nka),d1minus(sk,nkr,nkr,0:nka), &
            d2plus(sk,nkr,nkr,0:nka),d2minus(sk,nkr,nkr,0:nka), &
            db1plus(sk,nkr,nkr,0:nka),db1minus(sk,nkr,nkr,0:nka), &
            db2plus(sk,nkr,nkr,0:nka),db2minus(sk,nkr,nkr,0:nka), &
            dd1plus(sk,nkr,nkr,0:nka),dd1minus(sk,nkr,nkr,0:nka), &
            dd2plus(sk,nkr,nkr,0:nka),dd2minus(sk,nkr,nkr,0:nka), &
            ddb1plus(sk,nkr,nkr,0:nka),ddb1minus(sk,nkr,nkr,0:nka), &
            ddb2plus(sk,nkr,nkr,0:nka),ddb2minus(sk,nkr,nkr,0:nka), &
            dddb1(sk,nkr,nkr,0:nka),dddb2(sk,nkr,nkr,0:nka), &
            hie1(sk,nkr,nkr,0:nka),hie2(sk,nkr,nkr,0:nka), &
            hih1(sk,nkr,nkr,0:nka),hih2(sk,nkr,nkr,0:nka), &
            psipsi1(sk,nkr,nkr,0:nka),psipsi2(sk,nkr,nkr,0:nka), &
            uste1inv(sk,nkr,nkr,0:nka),uste2inv(sk,nkr,nkr,0:nka), &
            usth1inv(sk,nkr,nkr,0:nka),usth2inv(sk,nkr,nkr,0:nka))

!  массивы, необходимые всегда для вычисления правых частей в ур. прогонки
  allocate(ab(sk,nkr,0:nka), &
            b1plus(sk,nkr,nkr,0:nka), b1minus(sk,nkr,nkr,0:nka), &
            b2plus(sk,nkr,nkr,0:nka), b2minus(sk,nkr,nkr,0:nka), &
            bb1plus(sk,nkr,nkr,0:nka), bb1minus(sk,nkr,nkr,0:nka), &
            bb2plus(sk,nkr,nkr,0:nka), bb2minus(sk,nkr,nkr,0:nka), &
            brne1(sk,nkr,0:nka),brne2(sk,nkr,0:nka), &
            brnh1(sk,nkr,0:nka), &
            brnh2(sk,nkr,0:nka), &
            ddb1plusinv(sk,nkr,nkr,0:nka), &
            ddb1minusinv(sk,nkr,nkr,0:nka), &
            ddb2plusinv(sk,nkr,nkr,0:nka), &
            ddb2minusinv(sk,nkr,nkr,0:nka), &
            dddb1inv(sk,nkr,nkr,0:nka), dddb2inv(sk,nkr,nkr,0:nka), &
            fgammaplus(sk,nkr,nkr,0:nka), &
            fgammaminus(sk,nkr,nkr,0:nka), &
            hie2inv(sk,nkr,nkr,0:nka),hih2inv(sk,nkr,nkr,0:nka), &
            hipsisk(sk,nkr,0:nka), &
            psie1(sk,nkr,nkr,0:nka), psie2(sk,nkr,nkr,0:nka), &
            psih1(sk,nkr,nkr,0:nka), psih2(sk,nkr,nkr,0:nka), &
            psie1inv(sk,nkr,nkr,0:nka),psie2inv(sk,nkr,nkr,0:nka), &
            psih1inv(sk,nkr,nkr,0:nka),psih2inv(sk,nkr,nkr,0:nka), &
            psipsi1inv(sk,nkr,nkr,0:nka), &
            psipsi2inv(sk,nkr,nkr,0:nka), &
            psibrn1(sk,nkr,0:nka),psibrn2(sk,nkr,0:nka), &
            psihi1(sk,nkr,nkr,0:nka), psihi2(sk,nkr,nkr,0:nka), &
            psihi3(sk,nkr,nkr,0:nka), psihi4(sk,nkr,nkr,0:nka), &
            rab1(nkr,nkr), rab2(nkr,nkr), &
            rnplus(sk,nkr,0:nka), rnminus(sk,nkr,0:nka), &
            rnnplus(sk,nkr,0:nka), rnnminus(sk,nkr,0:nka), &
            uste1(sk,nkr,nkr,0:nka), uste2(sk,nkr,nkr,0:nka), &
            usth1(sk,nkr,nkr,0:nka), usth2(sk,nkr,nkr,0:nka), &
            xrab(nkr), &

!    это для прогонки
            alphap(sk,nkr,nkr),betap(sk+1,nkr), &
            rabp(nkr,nkr),rabpinv(nkr,nkr), &
            xrabp(nkr), exit_sum(nkr, 0:nka))
	allocate(temp(nkr,nkr), temp2(nkr,nkr))
  allocate(ipiv(nkr), ipiv2(nkr))
  continue
  return
end subroutine configfield


subroutine field_structure (k,kluch2)

! определение распределения поля по радиусу

  use array_work
	use field_nes
  implicit none

  integer inr,is,nr,k,kluch2,ir,iz,nz, ina
	real*8       radius,koordz
	complex*16 er(30),ez(30),a1,a2,a3,a4,ezenter,ezmid1,ezmid2,ezexit, ezaxis,ezbeam,ez65,ezmid(30)

	if (kluch2.eq.1) then
!     выдача поперечного распределения ez
  	nr=30
    er=0.0d0
    ez=0.0d0
    ezmid=0.0d0
    do ir=1,nr
    	ezenter=0.0d0
      ezmid1 =0.0d0
      ezmid2 =0.0d0
      ezexit =0.0d0
      radius = rt(sk)/30.0*ir
      do  inr=1,nkr
      	do ina=0,nka
        	is=1
          a3 = (xbplus(is,inr,ina)-xbminus(is,inr,ina) )
          a4 = (xbplus(is,inr,ina)+xbminus(is,inr,ina) )
          a1 = sqrt(abs(zn(is,inr,ina))/pi)*((gam(is,inr,ina))/(( abs(gam(is,inr,ina)))*(rt(is)) ) )
          a2 = bessel_jn(0, mu(inr,ina)*radius/rt(is))/bessel_jn(1, mu(inr,ina))*ce*mu(inr,ina)/(gam(is,inr,ina)*rt(is))
          ezenter=ezenter-a1*a2*a3
          a2 = bessel_jn(1, mu(inr,ina)*radius/rt(is))/bessel_jn(1, mu(inr,ina))
          er(ir)=er(ir)-a1*a2*a4

          is=(sk-1)/2
          a3 = ( xbplus(is,inr,ina)-xbminus(is,inr,ina) )
          a4 = ( xbplus(is,inr,ina)+xbminus(is,inr,ina) )
          a1 = sqrt(abs(zn(is,inr,ina))/pi)*((gam(is,inr,ina))/(( abs(gam(is,inr,ina)))*(rt(is)) ) )
          a2 = bessel_jn(0, mu(inr,ina)*radius/rt(is))/bessel_jn(1, mu(inr,ina))*ce*mu(inr,ina)/(gam(is,inr,ina)*rt(is))
          ezmid1=ezmid1-a1*a2*a3

          is=(sk-1)/2+1
          a3 = ( xbplus(is,inr,ina)-xbminus(is,inr,ina) )
          a4 = ( xbplus(is,inr,ina)+xbminus(is,inr,ina) )
          a1 = sqrt(abs(zn(is,inr,ina))/pi)*((gam(is,inr,ina))/(( abs(gam(is,inr,ina)))*(rt(is)) ) )
          a2 = bessel_jn(0, mu(inr,ina)*radius/rt(is))/bessel_jn(1, mu(inr,ina))*ce*mu(inr,ina)/(gam(is,inr,ina)*rt(is))
          ezmid2=ezmid2-a1*a2*a3

          is=sk
          a3 = ( xbplus(is,inr,ina)-xbminus(is,inr,ina) )
          a4 = ( xbplus(is,inr,ina)+xbminus(is,inr,ina) )
          a1 = sqrt(abs(zn(is,inr,ina))/pi)*((gam(is,inr,ina))/(( abs(gam(is,inr,ina)))*(rt(is)) ) )
          a2 = bessel_jn(0, mu(inr,ina)*radius/rt(is))/bessel_jn(1, mu(inr,ina))*ce*mu(inr,ina)/(gam(is,inr,ina)*rt(is))
          ezexit=ezexit-a1*a2*a3

          do iz=1,15
          	is=sk/2+2*iz
            a3 = ( xbplus(is,inr,ina)-xbminus(is,inr,ina) )
            a4 = ( xbplus(is,inr,ina)+xbminus(is,inr,ina) )
            a1 = sqrt(abs(zn(is,inr,ina))/pi)*((gam(is,inr,ina))/(( abs(gam(is,inr,ina)))*(rt(is)) ) )
            a2 = bessel_jn(0, mu(inr,ina)*radius/rt(is))/bessel_jn(1, mu(inr,ina))*ce*mu(inr,ina)/(gam(is,inr,ina)*rt(is))
            ezmid(iz)=ezmid(iz)-a1*a2*a3
          end do
        end do
      end do  ! конец суммирования по модам
			radius = rt(sk)/30.0*ir
			write(24,10) w0,k, radius,abs(ezenter), abs(ezmid1),abs(ezmid2) , abs(ezexit)

!     write(24,12) w0,k, radius,(real(ezmid(iz)),iz=1,15)
    end do
  end if

  if (kluch2.eq.2) then
!      выдача продольного распределения ez

  	do is=1,sk
      ezbeam =0.0d0
      ezaxis =0.0d0
      ez65   =0.0d0
      ez=0.0d0
      er=0.0d0
      koordz=zs(is)+dz(is)/2.0d0
      do  inr=1,nkr
      	do ina=0,nka
        	a3 = ( xbplus(is,inr,ina)-xbminus(is,inr,ina) )
          a4 = ( xbplus(is,inr,ina)+xbminus(is,inr,ina) )
          a1 = sqrt(abs(zn(is,inr,ina))/pi)*((gam(is,inr,ina))/(( abs(gam(is,inr,ina)))*(rt(is)) ) )

          do ir=1,9
						radius=0.125d0*(ir-1)*(rt2)
						a2 = bessel_jn(0, mu(inr,ina)*radius/rt(is))/bessel_jn(1, mu(inr,ina))*(ce*mu(inr,ina))/(gam(is,inr,ina)*rt(is))
						ez(ir)=ez(ir)+a1*a2*a3
						a2 = bessel_jn(1, mu(inr,ina)*radius/rt(is))/bessel_jn(1, mu(inr,ina))
            er(ir)=er(ir)-a1*a2*a4
          end do
        end do
      end do
    end do ! конец суммирования по модам
    write(25,11)  w0, k,is, koordz, rt(is),(real(ez(ir)),aimag(ez(ir)),ir=1,9)
    write(26,11)  w0, k,is, koordz, rt(is),(real(er(ir)),aimag(er(ir)),ir=1,9)
  end if

10  format('w0=',f9.4,1x,'k=',i4,1x,'rad_ez',1x,5(e10.2,1x))
11  format('w0=',f9.4,2x,'k=',i4,2x,i4,2x,'is_z_ez_er=',f6.2,2x,f6.2,1x,18(e10.2,1x))
12  format('w0=',f9.4,1x,'k=',i4,1x,'rad_ez',1x,16(e10.2,1x))

  return
end subroutine


subroutine field_power(k,kluch1)

! определение распределения поля по радиусу

  use array_work
	use field_nes

  integer inr,is,k,kluch1
	real power,powerenter,powerexit,powersum,powerim,powerplus,powerminus, dva_d_na_lambda
!	double complex xrab1(20)
  if (kluch1.eq.1) then
    do is=1,sk
      do ina=0,nka
        do inr=1,nkr
          xrab(inr)= (xbplus(is,inr,ina)+xbminus(is,inr,ina))*zn(is,inr,ina)/abs(zn(is,inr,ina))
        end do
        power=0.5d0*real(dot_product(xrab(:),(xbplus(is,:,ina)-xbminus(is,:,ina))))/(beam_curr*beam_voltage*1000.0d0)
        powerplus=0.5d0*real(dot_product(xbplus(is,:,ina),(xbplus(is,:,ina))))/(beam_curr*beam_voltage*1000.0d0)
        powerminus =0.5d0*real(dot_product(xbminus(is,:,ina),(xbminus(is,:,ina))))/(beam_curr*beam_voltage*1000.0d0)
        powerim=0.5d0*aimag(dot_product(xrab(:),(xbplus(is,:,ina)-xbminus(is,:,ina))))/(beam_curr*beam_voltage*1000.0d0)
        write (9,921) k,0,is,power,powerplus,powerminus,powerim
      end do
    end do
  end if
  if (kluch1.eq.2) then
    do inr=1,nkr
      do ina=0,nka
        exit_sum(inr, ina) = (xbplus(sk,inr,ina)+xbminus(sk,inr,ina))*zn(sk,inr,ina)/abs(zn(sk,inr,ina))
        xrabp(inr)= (xbplus(1,inr,ina)+xbminus(1,inr,ina))*zn(1,inr,ina)/abs(zn(1,inr,ina))
      end do
		end do
		do ina = 0, nka
			powerexit = powerexit + 0.5d0*real(dot_product(exit_sum(:, ina),(xbplus(sk,:,ina)-xbminus(sk,:,ina))))/(beam_curr*beam_voltage*1000.0d0)
		end do
		! print *, "beam_curr =", beam_curr, "beam_voltage", beam_voltage
    powerenter=0.5d0*real(dot_product(xrabp(:),(xbplus(1,:,0)-xbminus(1,:,0))))/(beam_curr*beam_voltage*1000.0d0)
    powersum= -( powerenter-powerexit) + sum
		print *, "POWER_EXIT =", powerexit
    ! print* ,'w0=',w0,'k',k,'power=',powerexit
    dva_d_na_lambda=periodz1*w0/(pi*3.0d0)

    write (21,922) dva_d_na_lambda,w0,beam_voltage,k,0,powerenter,powerexit,sum,powersum
  end if
921 format(1x,'k=',i5,2x,i4,2x,'is_p_p+_p-=',i4,2x,e10.3,2x,e10.3,2x,e10.3,2x,'im(p)=',e9.2)
922 format('2d/lam=',f7.3,2x,'w0=',f7.3,2x,'u0=',f8.2,2x,'k_pent_pex_sum_psum=',i4,2x,i4,2x,4(e10.3,2x))

    return
end subroutine


subroutine freeallmem
	use array_work

	deallocate(temp, temp2)
	deallocate(ipiv, ipiv2)
	deallocate(gam,zn,eznbm,ernbm, &
							xnplus,xnminus, &
							dxnplus,dxnminus, &
							xbplus,xbminus, &
							dxbplus,dxbminus, &
							bplus0,bminussk, &
							etaplus,etaminus, &
							amplxplus,amplxminus)
	deallocate(alpha, &
							yarray0,yarray,yarray1, &
							dery,aux, &
							velocity, velocity1)
	deallocate(aa1, aa2,aa3, &
							alfaplus, alfaminus, &
							betaplus, betaminus, &
							ab, &
				b1plus, b1minus, &
				b2plus, b2minus, &
				bb1plus, bb1minus, &
				bb2plus, bb2minus, &
				brne1,brne2,brnh1,brnh2, &
				ddb1plusinv, ddb1minusinv, &
				ddb2plusinv, ddb2minusinv, &
				dddb1inv, dddb2inv, &
				fgammaplus, fgammaminus, &
				hie2inv,hih2inv, &
				hipsisk, &
				psie1, psie2, &
				psih1, psih2, &
				psie1inv,psie2inv, &
				psih1inv,psih2inv, &
				psipsi1inv, psipsi2inv, &
				psibrn1,psibrn2, &
				psihi1, psihi2, &
				psihi3, psihi4, &
				rab1, rab2, &
				rnplus, rnminus,  &
				rnnplus, rnnminus, &
				uste1, uste2, &
				usth1, usth2, &
				xrab, &
				alphap,betap, &
				rabp,rabpinv,xrabp)
!	deallocate(bm1,bm2,b1m1)
	deallocate(d1plus,d1minus, &
			d2plus,d2minus, &
			db1plus, db1minus, &
			db2plus, db2minus, &
			dd1plus,dd1minus, &
			dd2plus,dd2minus, &
			ddb1plus, ddb1minus, &
			ddb2plus, ddb2minus, &
			dddb1, dddb2, &
			hie1, hie2,hih1, hih2, &
			psipsi2,psipsi1, &
			uste1inv,uste2inv, &
			usth1inv,usth2inv, &
			exit_sum)
	return
end subroutine
