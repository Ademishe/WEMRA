module beam_nes
use array_work
use com_prog
implicit none

contains

subroutine beam_begin
  integer i, ibeam

  mkns = mkk*2
  ddz = zs(1)/mk0
  do i = 1, mkk
    do ibeam = 1, nbeam
      all_yarray(2*i, ibeam) = ddz*(i-1)+ampl_mod*sin(2.0d0*pi*(i-1)/mk0)
      all_yarray(2*i-1, ibeam) = rel_factor*v0
      all_dery(2*i, ibeam) = v0/w0
      all_dery(2*i-1, ibeam) = 0.0d0
      velocity(i, ibeam) = v0
    end do
  end do

  return
end subroutine


subroutine beam_calc (tau)
	real*8 tau
	integer imm, k12, i, mk, mkk1, is, in, ina, ibeam
  double precision sum0,prmt(5)
  do ibeam = 1, nbeam
    yarray(:) = all_yarray(:, ibeam)
    dery(:) = all_dery(:, ibeam)
    prmt(1)=tau
    prmt(2)=tau+dt
    prmt(3)=step_prmt*dt
    prmt(4)=prmt_4
    prmt(5)=0.0d0
    imm=0
    call drkgs(prmt, mkns, imm, ibeam)
	  print *, 'imm=', imm
	  sum=0.0d0
    k12=mk0
    mk=0
    do i=1,mkk
      if(velocity(i, ibeam).ge.3.0d0) then
        print *,'error, particle velocity > c, k =',k12
      end if
      if(yarray(2*i).ge.zss .or. yarray(2*i-1).lt.0.3d0 ) then
        sum=sum+(1.0d0/sqrt(1-(velocity(i, ibeam)/3.0d0)**2)-1.0d0)
        mk=mk+1
      else
        k12=k12+1
        yarray1(2*k12)  = yarray(2*i)
        yarray1(2*k12-1)= yarray(2*i-1)
        velocity1(k12)  = velocity(i, ibeam)
      end if
    end do

    sum=(sum-mk0*(rel_factor-1.0d0))/(mk0*(rel_factor-1.0d0))
    k12=mk0
    mkk1=mkk-mk+mk0
    mkk=mkk1
    mkns=2*mkk

    do i=1,mk0
      yarray1(2*i)=ddz*(i-1)+ampl_mod*sin(2.0d0*pi*(i-1)/mk0)
      yarray1(2*i-1)=rel_factor*v0
      dery(2*i)=v0/w0
      dery(2*i-1)=0.0d0
	    velocity1(i)=v0
    end do

    do i=1,mkk
      velocity(i, ibeam)=velocity1(i)
	  end do

    do i=1,mkns
      yarray(i)=yarray1(i)
    end do
    all_dery(:, ibeam) = dery(:)
    all_yarray(:, ibeam) = yarray(:)
  end do ! end do ina

  if (kluch_beam.eq.2 .or. kluch_beam.eq.3) then
    open(unit = 22,file='eta.dat', access = 'APPEND')
  	do is = 1, sk
      do in = 1, nkr
        do ina = 0, nka
          write (22,901) is,in,abs(etaplus(is,in,ina)),abs(etaminus(is,in,ina)),mk0,mkk
        end do
      end do
    end do
    close(22)
  end if
  print *,mkk,mk0,mk,'  sum=',sum

901 format(1x,'is=',i5,2x,'in=',i4,2x,'etaplus=',e12.5,2x,'etaminus=',e12.5,2x,'mk0=',i6,2x,'mkk=',i5)
  return
end subroutine

end module beam_nes


subroutine fct(ttau, ibeam)
	use array_work
	real*8 ttau
  double precision  deltaz
  integer is, im, in, ina, ibeam

  do im=1,mkk
    is=0
3   is=is+1
    if((zs(is) - yarray(2*im))) 3,3,4
4   continue

    if (is.eq.1 .or. is.eq.sk+2) then
      dery(2*im-1)=0.0d0
      dery(2*im)=velocity(im,ibeam)/w0
    else
      dery(2*im-1)=0.0d0
      deltaz=yarray(2*im)-zs(is-1)-dz(is-1)/2.0d0
      do in = 1, nkr
        do ina = 0, nka
          dery(2*im-1) = dery(2*im-1) + const1*dreal(eznbm(is-1, in, ina, ibeam) * ((xnplus(is-1,in,ina)+dxnplus(is-1,in,ina)*deltaz)* &
            exp(ce*(ttau-gam(is-1,in,ina)*deltaz))-(xnminus(is-1,in,ina) + dxnminus(is-1,in,ina)*deltaz)*exp(ce*(ttau+gam(is-1,in,ina)*deltaz))))
          dery(2*im) = velocity(im, ibeam)/w0
        end do
      end do
    end if
  end do
  return
end subroutine

subroutine outp(ttau,irec,ndim,prmt,ttau0,itemp, ktemp, ibeam)
  use array_work
  double precision deltaz,ttau0,ttau,prmt(5)
  integer is, im, in, irec, ndim, itemp, ktemp, ibeam, ina

  if (prmt(2)-ttau.lt.prmt(3)/20.0d0) then
    prmt(5)=2.0d0
	end if

  do 2 im=1,mkk
    velocity(im,ibeam) = yarray(2*im-1)/sqrt((yarray(2*im-1)/3.0d0)**2 + 1.0d0)
    is = 0
3   is = is+1
    if((zs(is)-yarray(2*im))) 3,4,4
4   continue
    if (is.gt.1 .and. is.lt.sk+2) then
      deltaz=yarray(2*im)-zs(is-1)-dz(is-1)/2.0d0

      do in = 1, nkr
        do ina = 0, nka
          etaplus(is-1,in,ina)=(ttau-ttau0)*constq* &
                      (dconjg(eznbm(is-1,in,ina,ibeam))*exp(ce*(-ttau+dconjg(gam(is-1,in,ina))*deltaz))* &
                      velocity(im,ibeam)) + etaplus(is-1,in,ina)

          etaminus(is-1,in,ina)=-(ttau-ttau0)*constq* &
                      (dconjg(eznbm(is-1,in,ina,ibeam))*exp(ce*(-ttau-dconjg(gam(is-1,in,ina))*deltaz))* &
                      velocity(im,ibeam)) + etaminus(is-1,in,ina)
        end do
      end do
    end if
2   continue
    ttau0=ttau
    itemp=itemp+1
10  format(1x,'ktemp',i5,2x,'itemp',i4,2x,'ttau=  ',f17.14,2x,'irec=',i5)
    return
end subroutine
