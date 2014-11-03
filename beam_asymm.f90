module beam_nes
    use array_work
    use com_prog

	implicit none

	contains

subroutine beam_begin
    integer i
    mkns=mkk(inbeam)*2
    ddz=zs(1)/mk0

    do i=1,mkk(inbeam)
        yarray0(2*i,inbeam)=ddz*(i-1)+ampl_mod*sin(2.0d0*pi*(i-1)/mk0)
        yarray0(2*i-1,inbeam)=rel_factor*v0
        dery(2*i)=v0/w0
        dery(2*i-1)=0.0d0
        velocity(i)=v0
        yarray(2*i-1)=yarray0(2*i-1,inbeam)
        yarray(2*i)=yarray0(2*i,inbeam)
    end do

    return
end subroutine

subroutine beam_calc(tau)
    real*8 tau
	integer imm, k12, i, mk, mkk1, is,inr, ina
    double precision sum0,prmt(5)
    do i=1,2*mkk(inbeam)
        yarray=yarray0(:,inbeam)
    end do

    prmt(1)=tau
    prmt(2)=tau+dt
    prmt(3)=step_prmt*dt
    prmt(4)=prmt_4
    prmt(5)=0.0d0
    imm=0
    call drkgs(prmt,mkns,imm)
!           do i=1,mkk
!	    velocity(i)=yarray(2*i-1)/sqrt((yarray(2*i-1)/3.0d0)**2+1)
!	        yarray(2*i)=yarray(2*i)+velocity(i)*dt/w0
!           end do

    print *, 'imm=', imm
    sum=0.0d0    ! энергия, приобретенная вылетевшими частицами
    k12=mk0
    mk=0         !счетчик вылетающих частиц

    do i=1,mkk(inbeam)
! **************** сортировка частиц **************************
        if(velocity(i).ge.3.0d0) then
            print *,'error, particle velocity > c, k =',k12
        end if

        if(yarray(2*i).ge. zss .or. yarray(2*i-1).lt.0.3d0 ) then
            sum=sum+(1.0d0/sqrt(1-(velocity(i)/3.0d0)**2)-1.0d0)
            mk=mk+1
!                              эти частицы выходят из рассмотрения
        else
            k12=k12+1
            yarray1(2*k12)  = yarray(2*i)
            yarray1(2*k12-1)= yarray(2*i-1)
            velocity1(k12)  = velocity(i)
!          присвоение оставшимся частицам новых номеров
        end if
!****************конец сортировки*****************************
    end do

	     sum=(sum-mk0*(rel_factor-1.0d0))/(mk0*(rel_factor-1.0d0))

               k12=mk0
               mkk1=mkk(inbeam)-mk+mk0
               mkk(inbeam)=mkk1
               mkns=2*mkk(inbeam)

             do i=1,mk0
                 yarray1(2*i)=ddz*(i-1)+ampl_mod*sin(2.0d0*pi*(i-1)/mk0)
	           yarray1(2*i-1)=rel_factor*v0
                 dery(2*i)=v0/w0
                 dery(2*i-1)=0.0d0
	           velocity1(i)=v0

             end do

             do i=1,mkk(inbeam)
	           velocity(i)=velocity1(i)
	       end do

             do  i=1,mkns
                 yarray(i)=yarray1(i)
             end do

        yarray0(:,inbeam)=yarray(:)

       if (kluch_beam.eq.2) then
	   do is=1,sk
	   do inr=1,nkr
         do ina=0,nka
	    write (22,901) is,inr,abs(etaplus(is,inr,ina)),abs(etaminus(is,inr,ina)),mk0,mkk(inbeam)

	   end do
	   end do
         end do
       end if
       print *,mkk(inbeam),mk0,mk,'  sum=',sum

  901 format(1x,'is=',i5,2x,'in=',i4,2x,'etaplus=',e12.5,2x,'etaminus=',e12.5,2x,'mk0=',i6,2x,'mkk=',i5)
      return
	end subroutine

      end module beam_nes
!******************************************************************

subroutine fct(ttau)
!               вычисление правой части в уравнениях движения
!               вызывается из drkgs
      use array_work
	real*8 ttau
      double precision  deltaz
      integer is, im, inr

      do im=1,mkk(inbeam)
        is=0
    3   is=is+1
      if((zs(is)-yarray(2*im))) 3,3,4       ! на каком рег. участке нах.частица

    4  continue

!
 !     do while (zs(is)-yarray(2*im).lt.0.0d0)
!
!	   is=is+1
!

      if (is.eq.1.or.is.eq.sk+2) then
      dery(2*im-1)=0.0d0
      dery(2*im)=velocity(im)/w0
	 else
      dery(2*im-1)=0.0d0
      deltaz=yarray(2*im)-zs(is-1)-dz(is-1)/2.0d0
!               deltaz - текущая координата в локальной системе координат

         do inr=1,nkr
         do ina=0,nka

!                        уравнения движения
    dery(2*im-1)=dery(2*im-1)+const1*dreal(eznbm(is-1,inr,ina,inbeam)*(cos(ina*alpha(inbeam))* &
        ((xnplus(is-1,inr,ina)+dxnplus(is-1,inr,ina)*deltaz)*exp(ce*(ttau-gam(is-1,inr,ina)*deltaz))- &
        (xnminus(is-1,inr,ina)+dxnminus(is-1,inr,ina)*deltaz)*exp(ce*(ttau+gam(is-1,inr,ina)*deltaz)))+ &
        sin(ina*alpha(inbeam))*((xnplus(is-1,inr,ina)+dxnplus(is-1,inr,ina)*deltaz)* &
        exp(ce*(ttau-gam(is-1,inr,ina)*deltaz))-(xnminus(is-1,inr,ina)+dxnminus(is-1,inr,ina)*deltaz)* &
        exp(ce*(ttau+gam(is-1,inr,ina)*deltaz)))))
    dery(2*im)=velocity(im)/(w0)

            end do
            end do

        end if
    end do        !конец цикла по частицам
    return
end subroutine


subroutine outp(ttau,irec,ndim,prmt,ttau0,itemp, ktemp)

    use array_work
!      вычисление наведенных токов. управляется из drkgs
    double precision  deltaz,ttau0,ttau,prmt(5)
    integer is, im, inr,inbm,irec,ndim,itemp,ktemp

    if (prmt(2)-ttau.lt.prmt(3)/20.0d0) then
        prmt(5)=2.0d0
    end if             ! условие окончания интегрирования


    do im=1,mkk(inbeam)
        velocity(im)=yarray(2*im-1)/sqrt((yarray(2*im-1)/3.0d0)**2.0d0)

        is=0
3       is=is+1
        if((zs(is)-yarray(2*im))) 3,4,4       ! на каком рег. участке нах.частица

4           continue
        if (is.gt.1.and.is.lt.sk+2) then
            deltaz=yarray(2*im)-zs(is-1)-dz(is-1)/2.0d0

!      deltaz - текущая координата в лок.сист.координат


            do inr=1,nkr
                do ina=0,nka

! расчет ета -правой части для уравнений возбуждения

                etaplus(is-1,inr,ina)=(ttau-ttau0)*constq* &
                        sin(ina*alpha(inbeam))* &
                        (dconjg(eznbm(is-1,inr,ina,inbeam))* &
                        exp(ce*((-ttau)+dconjg(gam(is-1,inr,ina))*deltaz)))* &
                        velocity(im)+etaplus(is-1,inr,ina)

                etaminus(is-1,inr,ina)=-(ttau-ttau0)*constq* &
                        sin(ina*alpha(inbeam))* &
                        (dconjg(eznbm(is-1,inr,ina,inbeam))* &
                        (exp(ce*(-ttau-dconjg(gam(is-1,inr,ina))*deltaz))* &
                        velocity(im))+etaminus(is-1,inr,ina))

                etaplus(is-1,inr,ina)=(ttau-ttau0)*constq* &
                        cos(ina*alpha(inbeam))* &
                        (dconjg(eznbm(is-1,inr,ina,inbeam))* &
                        exp(ce*((-ttau)+dconjg(gam(is-1,inr,ina))*deltaz)))* &
                        velocity(im)+etaplus(is-1,inr,ina)

                etaminus(is-1,inr,ina)=-(ttau-ttau0)*constq* &
                        cos(ina*alpha(inbeam))* &
                        (dconjg(eznbm(is-1,inr,ina,inbeam))* &
                        (exp(ce*(-ttau-dconjg(gam(is-1,inr,ina))*deltaz))* &
                        velocity(im))+etaminus(is-1,inr,ina))

                end do     ! конец цикла по модам
            end do

        end if
    end do !конец цикла по частицам

      ttau0=ttau
!      write (9,10) ktemp,  itemp, ttau, irec

      itemp=itemp+1

!	if (itemp.eq.390) then
!	itemp=itemp
!	end if

!	end if
10  format(1x,'ktemp',i5,2x,'itemp',i4,2x,'ttau=  ',f17.14,2x,'irec=',i5)
    return
end subroutine
