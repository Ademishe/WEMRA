module field_nes
  use array_work
  use com_prog
  implicit none

  contains

subroutine init_arrays
  aa1(:,:,:,:)				  = (0.0d0, 0.0d0)
  aa2(:,:,:,:)				  = (0.0d0, 0.0d0)
  aa3(:,:,:,:)				  = (0.0d0, 0.0d0)
  fgammaplus(:,:,:,:,:)   = (0.0d0,0.0d0)
  fgammaminus(:,:,:,:,:)  = (0.0d0,0.0d0)
  d1plus(:,:,:,:,:)       = (0.0d0,0.0d0)
  d1minus(:,:,:,:,:)      = (0.0d0,0.0d0)
  d2plus(:,:,:,:,:)       = (0.0d0,0.0d0)
  d2minus(:,:,:,:,:)      = (0.0d0,0.0d0)
  dd1plus(:,:,:,:,:)      = (0.0d0,0.0d0)
  dd1minus(:,:,:,:,:)     = (0.0d0,0.0d0)
  dd2plus(:,:,:,:,:)      = (0.0d0,0.0d0)
  dd2minus(:,:,:,:,:)     = (0.0d0,0.0d0)
  b1plus(:,:,:,:,:)       = (0.0d0,0.0d0)
  b1minus(:,:,:,:,:)      = (0.0d0,0.0d0)
  b2plus(:,:,:,:,:)       = (0.0d0,0.0d0)
  b2minus(:,:,:,:,:)      = (0.0d0,0.0d0)
  bb1plus(:,:,:,:,:)      = (0.0d0,0.0d0)
  bb1minus(:,:,:,:,:)     = (0.0d0,0.0d0)
  bb2plus(:,:,:,:,:)      = (0.0d0,0.0d0)
  bb2minus(:,:,:,:,:)     = (0.0d0,0.0d0)
  uste1(:,:,:,:,:)        = (0.0d0,0.0d0)
  uste2(:,:,:,:,:)        = (0.0d0,0.0d0)
  usth2(:,:,:,:,:)        = (0.0d0,0.0d0)
  usth1(:,:,:,:,:)        = (0.0d0,0.0d0)
  rnplus(:,:,:)         = (0.0d0,0.0d0)
  rnminus(:,:,:)        = (0.0d0,0.0d0)
  rnnplus(:,:,:)        = (0.0d0,0.0d0)
  rnnminus(:,:,:)       = (0.0d0,0.0d0)
end subroutine


subroutine matrix_construct
  integer is,i,j,ina,jna
	complex*16 fdplus,fdminus, fbminus, fbplus
  integer info

  do is=1,sk
    do i=1,nkr
      do ina=0,nka
        do j=1,nkr
          do jna=0,nka
            if(is.gt.1 .and. (rt(is-1)-rt(is)).ge.(0.001d0)) then !уменьшается
              uste1(is,i,j,ina,jna)=(0.0d0,0.0d0)
              uste2(is,i,j,ina,jna)=pfunk(is,i,j,ina,jna,1,0)
              usth1(is,i,j,ina,jna)=conjg(pfunk(is,j,i,jna,ina,1,0))
              usth2(is,i,j,ina,jna)=(0.0d0,0.0d0)
            end if
            if(is.gt.1 .and. (rt(is)-rt(is-1)).gt. (0.001d0)) then !увеличивается
              uste1(is,i,j,ina,jna)=pfunk(is,i,j,ina,jna,0,1)   !менял i и j местами
              uste2(is,i,j,ina,jna)=(0.0d0,0.0d0)
              usth1(is,i,j,ina,jna)=(0.0d0,0.0d0)
              usth2(is,i,j,ina,jna)=conjg(pfunk(is,j,i,jna,ina,0,1))
            end if
          end do !end jna
        end do !end j

        fbplus= (abs(gam(is,i,ina))+gam0*gam0/abs(gam(is,i,ina)))/w0
        fbminus=(abs(gam(is,i,ina))-gam0*gam0/abs(gam(is,i,ina)))/w0
        fdplus=zn(is,i,ina)/abs(zn(is,i,ina))+conjg(zn(is,i,ina))/abs(zn(is,i,ina))
        fdminus=-zn(is,i,ina)/abs(zn(is,i,ina))+conjg(zn(is,i,ina))/abs(zn(is,i,ina))

        d1plus(is,i,i,ina,ina)=(fexp(is,i,-1,-1,ina)-fexp(is,i,1,-1,ina))*fdplus
        d1minus(is,i,i,ina,ina)=(fexp(is,i,-1,1,ina)-fexp(is,i,1,1,ina))*fdminus
        d2plus(is,i,i,ina,ina)=(fexp(is,i,1,1,ina)-fexp(is,i,-1,1,ina))*(-fdminus)
        d2minus(is,i,i,ina,ina)=(fexp(is,i,1,-1,ina)-fexp(is,i,-1,-1,ina))*(-fdplus)
        dd1plus(is,i,i,ina,ina)=(fexp(is,i,-1,-1,ina)+fexp(is,i,1,-1,ina))*fdplus*dz(is)/2.0d0
        dd1minus(is,i,i,ina,ina)=(fexp(is,i,-1,1,ina)+fexp(is,i,1,1,ina))*fdminus*dz(is)/2.0d0
        dd2plus(is,i,i,ina,ina)= (fexp(is,i,1,1,ina)+fexp(is,i,-1,1,ina))*(-fdminus)*dz(is)/2.0d0
        dd2minus(is,i,i,ina,ina)=(fexp(is,i,1,-1,ina)+fexp(is,i,-1,-1,ina))*(-fdplus)*dz(is)/2.0d0

        if (i.le.nopen(is)) then
          b1plus(is,i,i,ina,ina)= -fbplus*dz(is)
          b1minus(is,i,i,ina,ina)=-fbminus*funk1(is,i,ina)
          b2plus(is,i,i,ina,ina)= -fbminus*funk1(is,i,ina)
          b2minus(is,i,i,ina,ina)= -fbplus*dz(is)
          bb1plus(is,i,i,ina,ina)= (0.0d0,0.0d0)
          bb1minus(is,i,i,ina,ina)=-fbminus*(funk1(is,i,ina)/(2.0d0*ce*gam(is,i,ina))-funk2(is,i,ina)*dz(is)/2.0d0)
          bb2plus(is,i,i,ina,ina)= -fbminus*(funk2(is,i,ina)*dz(is)/2.0d0-funk1(is,i,ina)/(2.0d0*ce*gam(is,i,ina)))
          bb2minus(is,i,i,ina,ina)= (0.0d0,0.0d0)
        else
          b1plus(is,i,i,ina,ina)=  -fbplus*funk1(is,i,ina)
          b1minus(is,i,i,ina,ina)= -fbminus*dz(is)
          b2plus(is,i,i,ina,ina)=  -fbminus*dz(is)
          b2minus(is,i,i,ina,ina)= -fbplus*funk1(is,i,ina)
          bb1plus(is,i,i,ina,ina)= -fbplus*(funk1(is,i,ina)/(2.0d0*ce*gam(is,i,ina))-funk2(is,i,ina)*dz(is)/2.0d0)
          bb1minus(is,i,i,ina,ina)= (0.0d0,0.0d0)
          bb2plus(is,i,i,ina,ina)= (0.0d0,0.0d0)
          bb2minus(is,i,i,ina,ina)= -fbplus*(-funk1(is,i,ina)/(2.0d0*ce*gam(is,i,ina))+funk2(is,i,ina)*dz(is)/2.0d0)
        end if
        if(is.gt.1 .and. (rt(is-1)-rt(is)).gt. 0.001d0) then
          uste1(is,i,i,ina,ina)=zn(is-1,i,ina)/abs(zn(is-1,i,ina))
          usth2(is,i,i,ina,ina)=conjg(zn(is,i,ina))/abs(zn(is,i,ina))
        end if
        if(is.gt.1 .and. (rt(is)-rt(is-1)).gt. 0.001d0) then
          uste2(is,i,i,ina,ina)=zn(is,i,ina)/abs(zn(is,i,ina))
          usth1(is,i,i,ina,ina)=conjg(zn(is-1,i,ina))/abs(zn(is-1,i,ina))
        end if
        if(is.gt.1 .and. abs(rt(is)-rt(is-1)).le. 0.001) then
          uste1(is,i,i,ina,ina)= (1.0d0,0.0d0)
          uste2(is,i,i,ina,ina)= (1.0d0,0.0d0)
          usth1(is,i,i,ina,ina)= (1.0d0,0.0d0)
          usth2(is,i,i,ina,ina)= (1.0d0,0.0d0)
        end if
        fgammaplus(is,i,i,ina,ina)=exp(-ce*gam(is,i,ina)*dz(is)/2.0d0)
        fgammaminus(is,i,i,ina,ina)=exp(ce*gam(is,i,ina)*dz(is)/2.0d0)
      end do ! end ina
    end do ! end i

    db1plus(is,:,:,:,:)=d1plus(is,:,:,:,:)*dt/w0-b1plus(is,:,:,:,:)
    db1minus(is,:,:,:,:)=d1minus(is,:,:,:,:)*dt/w0-b1minus(is,:,:,:,:)
    ddb1plus(is,:,:,:,:)=dd1plus(is,:,:,:,:)*dt/w0-bb1plus(is,:,:,:,:)
    ddb1minus(is,:,:,:,:)=dd1minus(is,:,:,:,:)*dt/w0-bb1minus(is,:,:,:,:)
    db2plus(is,:,:,:,:)=d2plus(is,:,:,:,:)*dt/w0-b2plus(is,:,:,:,:)
    db2minus(is,:,:,:,:)=d2minus(is,:,:,:,:)*dt/w0-b2minus(is,:,:,:,:)
    ddb2plus(is,:,:,:,:)=dd2plus(is,:,:,:,:)*dt/w0-bb2plus(is,:,:,:,:)
    ddb2minus(is,:,:,:,:)=dd2minus(is,:,:,:,:)*dt/w0-bb2minus(is,:,:,:,:)

!        call dlincg (nkr, ddb1plus(is,:,:,ina),nkr,ddb1plusinv(is,:,:,ina), nkr)
!        call dlincg (nkr, ddb2plus(is,:,:,ina),nkr,ddb2plusinv(is,:,:,ina), nkr)
!        call dlincg (nkr, ddb1minus(is,:,:,ina),nkr,ddb1minusinv(is,:,:,ina),nkr)
!        call dlincg (nkr, ddb2minus(is,:,:,ina),nkr,ddb2minusinv(is,:,:,ina),nkr)

    do ina = 0, nka
      do jna = 0, nka
        ddb1plusinv(is,:,:,ina,jna) = ddb1plus(is,:,:,ina,jna)
        call zgetrf(nkr,nkr,ddb1plusinv(is,:,:,ina,jna),nkr,ipiv,info)
        call zgetri(nkr, ddb1plusinv(is,:,:,ina,jna),nkr,ipiv,temp, nkr, info)
        ddb2plusinv(is,:,:,ina,jna) = ddb2plus(is,:,:,ina,jna)
        call zgetrf(nkr,nkr,ddb2plusinv(is,:,:,ina,jna),nkr,ipiv,info)
        call zgetri(nkr, ddb2plusinv(is,:,:,ina,jna),nkr,ipiv,temp, nkr, info)
        ddb1minusinv(is,:,:,ina,jna) = ddb1minus(is,:,:,ina,jna)
        call zgetrf(nkr,nkr,ddb1minusinv(is,:,:,ina,jna),nkr,ipiv,info)
        call zgetri(nkr, ddb1minusinv(is,:,:,ina,jna),nkr,ipiv,temp, nkr, info)
        ddb2minusinv(is,:,:,ina,jna) = ddb2minus(is,:,:,ina,jna)
        call zgetrf(nkr,nkr,ddb2minusinv(is,:,:,ina,jna),nkr,ipiv,info)
        call zgetri(nkr, ddb2minusinv(is,:,:,ina,jna),nkr,ipiv,temp, nkr, info)

        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),ddb2plusinv(is,:,:,ina,jna),nkr, ddb1plus(is,:,:,ina,jna),nkr,(0.0d0,0.0d0),dddb1(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),ddb2minusinv(is,:,:,ina,jna),nkr, db1minus(is,:,:,ina,jna),nkr,(1.0d0,0.0d0),dddb1(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),ddb1plusinv(is,:,:,ina,jna),nkr, ddb2plus(is,:,:,ina,jna),nkr,(0.0d0,0.0d0),dddb2(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),ddb1minusinv(is,:,:,ina,jna),nkr, ddb2minus(is,:,:,ina,jna),nkr,(1.0d0,0.0d0),dddb2(is,:,:,ina,jna),nkr)

    !        call dlincg (nkr, dddb1(is,:,:,ina),nkr,dddb1inv(is,:,:,ina), nkr)
    !        call dlincg (nkr, dddb2(is,:,:,ina),nkr,dddb2inv(is,:,:,ina), nkr)

        dddb1inv(is,:,:,ina,jna) = dddb1(is,:,:,ina,jna)
        call zgetrf(nkr,nkr,dddb1inv(is,:,:,ina,jna),nkr,ipiv,info)
        call zgetri(nkr, dddb1inv(is,:,:,ina,jna),nkr,ipiv,temp, nkr, info)
        dddb2inv(is,:,:,ina,jna) = dddb2(is,:,:,ina,jna)
        call zgetrf(nkr,nkr,dddb2inv(is,:,:,ina,jna),nkr,ipiv,info)
        call zgetri(nkr, dddb2inv(is,:,:,ina,jna),nkr,ipiv,temp, nkr, info)

        call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),ddb2plusinv(is,:,:,ina,jna),nkr,db1plus(is,:,:,ina,jna),nkr,(0.0d0,0.0d0),rab1,nkr)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),ddb2minusinv(is,:,:,ina,jna),nkr,db1minus(is,:,:,ina,jna),nkr,(1.0d0,0.0d0),rab1,nkr)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),dddb1inv(is,:,:,ina,jna),nkr,rab1,nkr,(0.0d0,0.0d0),alfaplus(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),ddb1plusinv(is,:,:,ina,jna),nkr,db1plus(is,:,:,ina,jna),nkr,(0.0d0,0.0d0),rab1,nkr)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),ddb1minusinv(is,:,:,ina,jna),nkr,db1minus(is,:,:,ina,jna),nkr,(1.0d0,0.0d0),rab1,nkr)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),dddb2inv(is,:,:,ina,jna),nkr,rab1,nkr,(0.0d0,0.0d0),alfaminus(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),ddb2plusinv(is,:,:,ina,jna),nkr,db2plus(is,:,:,ina,jna),nkr,(0.0d0,0.0d0),rab1,nkr)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),ddb2minusinv(is,:,:,ina,jna),nkr,db2minus(is,:,:,ina,jna),nkr,(1.0d0,0.0d0),rab1,nkr)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),dddb1inv(is,:,:,ina,jna),nkr,rab1,nkr,(0.0d0,0.0d0),betaplus(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),ddb1plusinv(is,:,:,ina,jna),nkr,db2plus(is,:,:,ina,jna),nkr,(0.0d0,0.0d0),rab1,nkr)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),ddb1minusinv(is,:,:,ina,jna),nkr,db2minus(is,:,:,ina,jna),nkr,(1.0d0,0.0d0),rab1,nkr)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),dddb2inv(is,:,:,ina,jna),nkr,rab1,nkr,(0.0d0,0.0d0),betaminus(is,:,:,ina,jna),nkr)
      end do
    end do

  end do	! end is

  do is=2,sk
    do ina = 0, nka
      do jna = 0, nka
        call zgemm('n','n',nkr,nkr,nkr,dz(is-1)/2.0d0*ee,fgammaplus(is-1,:,:,ina,jna),nkr,alfaplus(is-1,:,:,ina,jna),nkr,(0.0d0,0.0d0),rab1,nkr)
        call zgemm('n','n',nkr,nkr,nkr,dz(is-1)/2.0d0*ee,fgammaminus(is-1,:,:,ina,jna),nkr,alfaminus(is-1,:,:,ina,jna),nkr,(1.0d0,0.0d0),rab1,nkr)

        rab1=rab1+fgammaplus(is-1,:,:,ina,jna)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),uste1(is,:,:,ina,jna),nkr,rab1,nkr,(0.0d0,0.0d0),hie1(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,dz(is-1)/2.0d0*ee,fgammaplus(is-1,:,:,ina,jna),nkr,alfaplus(is-1,:,:,ina,jna),nkr,(0.0d0,0.0d0),rab1,nkr)
        call zgemm('n','n',nkr,nkr,nkr,-dz(is-1)/2.0d0*ee,fgammaminus(is-1,:,:,ina,jna),nkr,alfaminus(is-1,:,:,ina,jna),nkr,(1.0d0,0.0d0),rab1,nkr)

        rab1=rab1+fgammaplus(is-1,:,:,ina,jna)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),usth1(is,:,:,ina,jna),nkr,rab1,nkr,(0.0d0,0.0d0),hih1(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,-dz(is)/2.0d0*ee,fgammaminus(is,:,:,ina,jna),nkr,alfaplus(is,:,:,ina,jna),nkr,(0.0d0,0.0d0),rab1,nkr)
        call zgemm('n','n',nkr,nkr,nkr,-dz(is)/2.0d0*ee,fgammaplus(is,:,:,ina,jna),nkr,alfaminus(is,:,:,ina,jna),nkr,(1.0d0,0.0d0),rab1,nkr)

        rab1=rab1+fgammaminus(is,:,:,ina,jna)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),uste2(is,:,:,ina,jna),nkr,rab1,nkr,(0.0d0,0.0d0),hie2(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,-dz(is)/2.0d0*ee,fgammaminus(is,:,:,ina,jna),nkr,alfaplus(is,:,:,ina,jna),nkr,(0.0d0,0.0d0),rab1,nkr)
        call zgemm('n','n',nkr,nkr,nkr,dz(is)/2.0d0*ee,fgammaplus(is,:,:,ina,jna),nkr,alfaminus(is,:,:,ina,jna),nkr,(1.0d0,0.0d0),rab1,nkr)

        rab1=rab1+fgammaminus(is,:,:,ina,jna)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),usth2(is,:,:,ina,jna),nkr,rab1,nkr,(0.0d0,0.0d0),hih2(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,dz(is-1)/2.0d0*ee,fgammaplus(is-1,:,:,ina,jna),nkr,betaplus(is-1,:,:,ina,jna),nkr,(0.0d0,0.0d0),rab1,nkr)
        call zgemm('n','n',nkr,nkr,nkr,dz(is-1)/2.0d0*ee,fgammaminus(is-1,:,:,ina,jna),nkr,betaminus(is-1,:,:,ina,jna),nkr,(1.0d0,0.0d0),rab1,nkr)

        rab1=rab1+fgammaminus(is-1,:,:,ina,jna)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),uste1(is,:,:,ina,jna),nkr,rab1,nkr,(0.0d0,0.0d0),psie1(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,dz(is-1)/2.0d0*ee,fgammaplus(is-1,:,:,ina,jna),nkr,betaplus(is-1,:,:,ina,jna),nkr,(0.0d0,0.0d0),rab1,nkr)
        call zgemm('n','n',nkr,nkr,nkr,-dz(is-1)/2.0d0*ee,fgammaminus(is-1,:,:,ina,jna),nkr,betaminus(is-1,:,:,ina,jna),nkr,(1.0d0,0.0d0),rab1,nkr)

        rab1=rab1-fgammaminus(is-1,:,:,ina,jna)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),usth1(is,:,:,ina,jna),nkr,rab1,nkr,(0.0d0,0.0d0),psih1(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,-dz(is)/2.0d0*ee,fgammaminus(is,:,:,ina,jna),nkr,betaplus(is,:,:,ina,jna),nkr,(0.0d0,0.0d0),rab1,nkr)
        call zgemm('n','n',nkr,nkr,nkr,-dz(is)/2.0d0*ee,fgammaplus(is,:,:,ina,jna),nkr,betaminus(is,:,:,ina,jna),nkr,(1.0d0,0.0d0),rab1,nkr)

        rab1=rab1+fgammaplus(is,:,:,ina,jna)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),uste2(is,:,:,ina,jna),nkr,rab1,nkr,(0.0d0,0.0d0),psie2(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,-dz(is)/2.0d0*ee,fgammaminus(is,:,:,ina,jna),nkr,betaplus(is,:,:,ina,jna),nkr,(0.0d0,0.0d0),rab1,nkr)
        call zgemm('n','n',nkr,nkr,nkr,dz(is)/2.0d0*ee,fgammaplus(is,:,:,ina,jna),nkr,betaminus(is,:,:,ina,jna),nkr,(1.0d0,0.0d0),rab1,nkr)

        rab1=rab1-fgammaplus(is,:,:,ina,jna)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),usth2(is,:,:,ina,jna),nkr,rab1,nkr,(0.0d0,0.0d0),psih2(is,:,:,ina,jna),nkr)


    !        call dlincg (nkr, psie1(is,:,:,ina,jna),nkr,psie1inv(is,:,:,ina,jna), nkr)
    !        call dlincg (nkr, psie2(is,:,:,ina,jna),nkr,psie2inv(is,:,:,ina,jna), nkr)
    !        call dlincg (nkr, psih1(is,:,:,ina,jna),nkr,psih1inv(is,:,:,ina,jna), nkr)
    !        call dlincg (nkr, psih2(is,:,:,ina,jna),nkr,psih2inv(is,:,:,ina,jna), nkr)
        psie1inv(is,:,:,ina,jna) = psie1(is,:,:,ina,jna)
        call zgetrf(nkr,nkr,psie1inv(is,:,:,ina,jna),nkr,ipiv,info)
        call zgetri(nkr, psie1inv(is,:,:,ina,jna),nkr,ipiv,temp, nkr, info)
        psie2inv(is,:,:,ina,jna) = psie2(is,:,:,ina,jna)
        call zgetrf(nkr,nkr,psie2inv(is,:,:,ina,jna),nkr,ipiv,info)
        call zgetri(nkr, psie2inv(is,:,:,ina,jna),nkr,ipiv,temp, nkr, info)
        psih1inv(is,:,:,ina,jna) = psih1(is,:,:,ina,jna)
        call zgetrf(nkr,nkr,psih1inv(is,:,:,ina,jna),nkr,ipiv,info)
        call zgetri(nkr, psih1inv(is,:,:,ina,jna),nkr,ipiv,temp, nkr, info)
        psih2inv(is,:,:,ina,jna) = psih2(is,:,:,ina,jna)
        call zgetrf(nkr,nkr,psih2inv(is,:,:,ina,jna),nkr,ipiv,info)
        call zgetri(nkr, psih2inv(is,:,:,ina,jna),nkr,ipiv,temp, nkr, info)


        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),psie1inv(is,:,:,ina,jna),nkr,psie2(is,:,:,ina,jna),nkr,(0.0d0,0.0d0),psipsi1(is,:,:,ina,jna),nkr)
         call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psih1inv(is,:,:,ina,jna),nkr,psih2(is,:,:,ina,jna),nkr,(1.0d0,0.0d0),psipsi1(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),psie1inv(is,:,:,ina,jna),nkr,hie1(is,:,:,ina,jna),nkr,(0.0d0,0.0d0),psihi1(is,:,:,ina,jna),nkr)
         call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psih1inv(is,:,:,ina,jna),nkr,hih1(is,:,:,ina,jna),nkr,(1.0d0,0.0d0),psihi1(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psie1inv(is,:,:,ina,jna),nkr,hie2(is,:,:,ina,jna),nkr,(0.0d0,0.0d0),psihi2(is,:,:,ina,jna),nkr)
         call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),psih1inv(is,:,:,ina,jna),nkr,hih2(is,:,:,ina,jna),nkr,(1.0d0,0.0d0),psihi2(is,:,:,ina,jna),nkr)

    !        call dlincg(nkr, psipsi1(is,:,:,ina,jna),nkr,psipsi1inv(is,:,:,ina,jna),nkr)
        psipsi1inv(is,:,:,ina,jna) = psipsi1(is,:,:,ina,jna)
        call zgetrf(nkr,nkr,psipsi1inv(is,:,:,ina,jna),nkr,ipiv,info)
        call zgetri(nkr, psipsi1inv(is,:,:,ina,jna),nkr,ipiv,temp, nkr, info)

      end do !end jna
    end do !end ina
  end do !end is

  do is=1,sk-1
    do ina = 0, nka
      do jna = 0, nka
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),psie2inv(is+1,:,:,ina,jna),nkr,psie1(is+1,:,:,ina,jna),nkr,(0.0d0,0.0d0),psipsi2(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psih2inv(is+1,:,:,ina,jna),nkr,psih1(is+1,:,:,ina,jna),nkr,(1.0d0,0.0d0),psipsi2(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),psie2inv(is+1,:,:,ina,jna),nkr,hie2(is+1,:,:,ina,jna),nkr,(0.0d0,0.0d0),psihi3(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psih2inv(is+1,:,:,ina,jna),nkr,hih2(is+1,:,:,ina,jna),nkr,(1.0d0,0.0d0),psihi3(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psie2inv(is+1,:,:,ina,jna),nkr,hie1(is+1,:,:,ina,jna),nkr,(0.0d0,0.0d0),psihi4(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),psih2inv(is+1,:,:,ina,jna),nkr,hih1(is+1,:,:,ina,jna),nkr,(1.0d0,0.0d0),psihi4(is,:,:,ina,jna),nkr)
!        call dlincg (nkr, psipsi2(is,:,:,ina,jna),nkr,psipsi2inv(is,:,:,ina,jna),nkr)
        psipsi2inv(is,:,:,ina,jna) = psipsi2(is,:,:,ina,jna)
        call zgetrf(nkr,nkr,psipsi2inv(is,:,:,ina,jna),nkr,ipiv,info)
        call zgetri(nkr, psipsi2inv(is,:,:,ina,jna),nkr,ipiv,temp, nkr, info)
      end do !end jna
    end do !end ina
  end do !end is
  do is=2,sk-1
    do ina = 0, nka
      do jna = 0, nka
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),psipsi1inv(is,:,:,ina,jna),nkr,psihi1(is,:,:,ina,jna),nkr,(0.0d0,0.0d0),aa1(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),psipsi1inv(is,:,:,ina,jna),nkr,psihi2(is,:,:,ina,jna),nkr,(0.0d0,0.0d0),aa2(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psipsi2inv(is,:,:,ina,jna),nkr,psihi4(is,:,:,ina,jna),nkr,(1.0d0,0.0d0),aa2(is,:,:,ina,jna),nkr)
        call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psipsi2inv(is,:,:,ina,jna),nkr,psihi3(is,:,:,ina,jna),nkr,(0.0d0,0.0d0),aa3(is,:,:,ina,jna),nkr)
      end do !end jna
    end do !end ina
  end do ! end is

!        call dlincg (nkr, hie2(sk,:,:,ina,jna),nkr,hie2inv(sk,:,:,ina,jna),nkr)
!        call dlincg (nkr, hih2(sk,:,:,ina,jna),nkr,hih2inv(sk,:,:,ina,jna),nkr)
  do ina = 0, nka
    do jna = 0, nka
      hie2inv(sk,:,:,ina,jna) = hie2(sk,:,:,ina,jna)
      call zgetrf(nkr,nkr,hie2inv(sk,:,:,ina,jna),nkr,ipiv,info)
      call zgetri(nkr, hie2inv(sk,:,:,ina,jna),nkr,ipiv,temp, nkr, info)
      hih2inv(sk,:,:,ina,jna) = hih2(sk,:,:,ina,jna)
      call zgetrf(nkr,nkr,hih2inv(sk,:,:,ina,jna),nkr,ipiv,info)
      call zgetri(nkr, hih2inv(sk,:,:,ina,jna),nkr,ipiv,temp, nkr, info)

      aa1(sk,:,:,ina,jna)=psihi1(sk,:,:,ina,jna)
      aa2(sk,:,:,ina,jna)=psihi2(sk,:,:,ina,jna)
    end do !end jna
  end do !end ina

	return
end subroutine matrix_construct


subroutine matrix_dynamic
  integer is
  do is=1,sk
!*********************************************************
!      rn   plus    minus
!**********************************************************
    call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),b1plus(is,:,:,ina),nkr,xnplus(is,:,ina),1,(0.0d0,0.0d0),rnplus(is,:,ina),1)
    call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),bb1plus(is,:,:,ina),nkr,dxnplus(is,:,ina),1,(1.0d0,0.0d0),rnplus(is,:,ina),1)
    call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),b2plus(is,:,:,ina),nkr,xnminus(is,:,ina),1,(1.0d0,0.0d0),rnplus(is,:,ina),1)
    call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),bb2plus(is,:,:,ina),nkr,dxnminus(is,:,ina),1,(1.0d0,0.0d0),rnplus(is,:,ina),1)
    rnplus(is,:,ina)= rnplus(is,:,ina)+etaplus(is,:,ina)

    call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),b1minus(is,:,:,ina),nkr,xnplus(is,:,ina),1,(0.0d0,0.0d0),rnminus(is,:,ina),1)
    call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),bb1minus(is,:,:,ina),nkr,dxnplus(is,:,ina),1,(1.0d0,0.0d0),rnminus(is,:,ina),1)
    call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),b2minus(is,:,:,ina),nkr,xnminus(is,:,ina),1,(1.0d0,0.0d0),rnminus(is,:,ina),1)
    call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),bb2minus(is,:,:,ina),nkr,dxnminus(is,:,ina),1,(1.0d0,0.0d0),rnminus(is,:,ina),1)
    rnminus(is,:,ina)=rnminus(is,:,ina)+etaminus(is,:,ina)
!*********************************************************
!     rnn plus   minus
!*********************************************************

    call zgemv('n',nkr,nkr,(1.0d0,0.0d0),ddb2plusinv(is,:,:,ina),nkr,rnplus(is,:,ina),1,(0.0d0,0.0d0),rab1,1)
    call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),ddb2minusinv(is,:,:,ina),nkr,rnminus(is,:,ina),1,(1.0d0,0.0d0),rab1,1)
    call zgemv('n',nkr,nkr,(1.0d0,0.0d0),dddb1inv(is,:,:,ina),nkr,rab1,1,(0.0d0,0.0d0),rnnplus(is,:,ina),1)
    call zgemv('n',nkr,nkr,(1.0d0,0.0d0),ddb1plusinv(is,:,:,ina),nkr,rnplus(is,:,ina),1,(0.0d0,0.0d0),rab1,1)
    call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),ddb1minusinv(is,:,:,ina),nkr,rnminus(is,:,ina),1,(1.0d0,0.0d0),rab1,1)
    call zgemv('n',nkr,nkr,(1.0d0,0.0d0),dddb2inv(is,:,:,ina),nkr,rab1,1,(0.0d0,0.0d0),rnnminus(is,:,ina),1)
  end do

!************************************************************
!    brne1   brne2  brnh1  brnh2
!************************************************************

    do is=2,sk
!       brne1
      call zgemv('n',nkr,nkr,dz(is-1)/2.0d0*ee,fgammaminus(is-1,:,:,ina),nkr,rnnminus(is-1,:,ina),1,(0.0d0,0.0d0),xrab,1)
      call zgemv('n',nkr,nkr,dz(is-1)/2.0d0*ee,fgammaplus(is-1,:,:,ina),nkr,rnnplus(is-1,:,ina),1,(1.0d0,0.0d0),xrab,1)
      call zgemv('n',nkr,nkr,(1.0d0,0.0d0),uste1(is,:,:,ina),nkr,xrab,1,(0.0d0,0.0d0),brne1(is,:,ina),1)

!     brne2
      call zgemv('n',nkr,nkr,-dz(is)/2.0d0*ee,fgammaminus(is,:,:,ina),nkr,rnnplus(is,:,ina),1,(0.0d0,0.0d0),xrab,1)
      call zgemv('n',nkr,nkr,-dz(is)/2.0d0*(1.0d0,0.0d0),fgammaplus(is,:,:,ina),nkr,rnnminus(is,:,ina),1,(1.0d0,0.0d0),xrab,1)
      call zgemv('n',nkr,nkr,(1.0d0,0.0d0),uste2(is,:,:,ina),nkr,xrab,1,(0.0d0,0.0d0),brne2(is,:,ina),1)

!      brnh1
      call zgemv('n',nkr,nkr,-dz(is-1)/2.0d0*(1.0d0,0.0d0),fgammaminus(is-1,:,:,ina),nkr,rnnminus(is-1,:,ina),1,(0.0d0,0.0d0),xrab,1)
      call zgemv('n',nkr,nkr,dz(is-1)/2.0d0*(1.0d0,0.0d0),fgammaplus(is-1,:,:,ina),nkr,rnnplus(is-1,:,ina),1,(1.0d0,0.0d0),xrab,1)
      call zgemv('n',nkr,nkr,(1.0d0,0.0d0),usth1(is,:,:,ina),nkr,xrab,1,(0.0d0,0.0d0),brnh1(is,:,ina),1)

!      brnh2
      call zgemv('n',nkr,nkr,-dz(is)/2.0d0*(1.0d0,0.0d0),fgammaminus(is,:,:,ina),nkr,rnnplus(is,:,ina),1,(0.0d0,0.0d0),xrab,1)
      call zgemv('n',nkr,nkr,dz(is)/2.0*ee,fgammaplus(is,:,:,ina),nkr,rnnminus(is,:,ina),1,(1.0d0,0.0d0),xrab,1)
      call zgemv('n',nkr,nkr,(1.0d0,0.0d0),usth2(is,:,:,ina),nkr,xrab,1,(0.0d0,0.0d0),brnh2(is,:,ina),1)

!************************************************************
!           psibrn1    and  psibrn2
!************************************************************
        xrab=brnh1(is,:,ina)-brnh2(is,:,ina)
        call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),psih1inv(is,:,:,ina),nkr,xrab,1,(0.0d0,0.0d0),psibrn1(is,:,ina),1)

        xrab=brne1(is,:,ina)-brne2(is,:,ina)
        call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psie1inv(is,:,:,ina),nkr,xrab,1,(1.0d0,0.0d0),psibrn1(is,:,ina),1)

    end do

    do is=1,sk-1
      xrab=brnh2(is+1,:,ina)-brnh1(is+1,:,ina)
      call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),psih2inv(is+1,:,:,ina),nkr,xrab,1,(0.0d0,0.0d0),psibrn2(is,:,ina),1)

      xrab=brne2(is+1,:,ina)-brne1(is+1,:,ina)
      call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psie2inv(is+1,:,:,ina),nkr,xrab,1,(1.0d0,0.0d0),psibrn2(is,:,ina),1)
    end do

!*************************************************************
!      right part in regular equation - ab
!*************************************************************
    do is=2,sk-1
      call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),psipsi1inv(is,:,:,ina),nkr,psibrn1(is,:,ina),1,(0.0d0,0.0d0),ab(is,:,ina),1)
      call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psipsi2inv(is,:,:,ina),nkr,psibrn2(is,:,ina),1,(1.0d0,0.0d0),ab(is,:,ina),1)
    end do
!***************************************************************
!       right part in    enter equation
!***************************************************************
    call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psipsi2inv(2,:,:,ina),nkr,psibrn2(2,:,ina),1,(0.0d0,0.0d0),ab(2,:,ina),1)
    call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),psipsi1inv(2,:,:,ina),nkr,psibrn1(2,:,ina),1,(1.0d0,0.0d0),ab(2,:,ina),1)

    call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psipsi1inv(2,:,:,ina),nkr,psihi1(2,:,:,ina),nkr,(0.0d0,0.0d0),rab1,nkr)
    call zgemv('n',nkr,nkr,(1.0d0,0.0d0),rab1,nkr,bplus0(:,ina),1,(1.0d0,0.0d0),ab(2,:,ina),1)

    call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psipsi1(sk,:,:,ina),nkr,bminussk(:,ina),1,(0.0d0,0.0d0),ab(sk,:,ina),1)
    ab(sk,:,ina)=ab(sk,:,ina)-psibrn1(sk,:,ina)
    return
end subroutine matrix_dynamic


subroutine field_calc
    integer is
!**************************************************************
! determination xb minus
!**************************************************************
    do ina=0,nka
        do is=2,sk-1
            call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psihi3(is,:,:,ina),nkr,xbplus(is+1,:,ina),1,(0.0d0,0.0d0),xrab,1)
            call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psihi4(is,:,:,ina),nkr,xbplus(is,:,:),1,(1.0d0,0.0d0),xrab,1)

            xrab=xrab+psibrn2(is,:,ina)
            call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psipsi2inv(is,:,:,ina),nkr,xrab,1,(0.0d0,0.0d0),xbminus(is,:,ina),1)
            call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psihi1(is,:,:,ina),nkr,xbplus(is-1,:,ina),1,(0.0d0,0.0d0),xrab,1)
            call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psihi2(is,:,:,ina),nkr,xbplus(is,:,ina),1,(1.0d0,0.0d0),xrab,1)

            xrab=xrab+psibrn1(is,:,ina)
            call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psipsi1inv(is,:,:,ina),nkr,xrab,1,(0.0d0,0.0d0),xbminus(is,:,ina),1)

        end do

        call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psihi4(1,:,:,ina),nkr,bplus0(:,ina),1,(0.0d0,0.0d0),xrab,1)
        call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psihi3(1,:,:,ina),nkr,xbplus(2,:,ina),1,(1.0d0,0.0d0),xrab,1)

        xrab=xrab+psibrn2(1,:,ina)
        call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psipsi2inv(1,:,:,ina),nkr,xrab,1,(0.0d0,0.0d0),xbminus(1,:,ina),1)

        do is=1,sk
            call zgemv('n',nkr,nkr,(1.0d0,0.0d0),betaplus(is,:,:,ina),nkr,xbminus(is,:,ina),1,(0.0d0,0.0d0),dxbplus(is,:,ina),1)
            call zgemv('n',nkr,nkr,(1.0d0,0.0d0),alfaplus(is,:,:,ina),nkr,xbplus(is,:,ina),1,(1.0d0,0.0d0),dxbplus(is,:,ina),1)

            dxbplus(is,:,ina)=dxbplus(is,:,ina)+rnnplus(is,:,ina)
            call zgemv('n',nkr,nkr,(1.0d0,0.0d0),betaminus(is,:,:,ina),nkr,xbminus(is,:,ina),1,(0.0d0,0.0d0),dxbminus(is,:,ina),1)
            call zgemv('n',nkr,nkr,(1.0d0,0.0d0),alfaminus(is,:,:,ina),nkr,xbplus(is,:,ina),1,(1.0d0,0.0d0),dxbminus(is,:,ina),1)

            dxbminus(is,:,ina)=dxbminus(is,:,ina)+rnnminus(is,:,ina)
        end do
    end do !ina
    return
end subroutine field_calc



double complex function fexp(is,i,k1,k2,ina)
	integer is,i,k1,k2,ina
!	 complex*16 gamma1,gamma2
!	 real*8 zl
  fexp=exp(ce*k1*(gam(is,i,ina)+k2*conjg(gam(is,i,ina)))*(0.5d0*dz(is)))
end function


double complex function funk1(is,inr,ina)
	integer is,inr,ina
!	 complex*16 gamma1,gamma2
!	 real*8 zl
  funk1=(exp(ce*gam(is,inr,ina)*dz(is))-exp(-ce*gam(is,inr,ina)*dz(is)))/(ce*2*gam(is,inr,ina))
end function


double complex function funk2(is,inr,ina)
  integer is,inr,ina
!	 complex*16 gamma1,gamma2
!	 real*8 zl
  funk2=(exp(ce*gam(is,inr,ina)*dz(is))+exp(-ce*gam(is,inr,ina)*dz(is)))/(ce*2*gam(is,inr,ina))
end function


double complex function pfunk(is,l,m,k,n,Ec,Hc)
!   pfunk(номер участка, порядок, порядок, ключ для знака, номер моды)
  integer is, m, l, n, k, Ec, Hc
  complex*16 const_part, part_1, part_2, E_nm, E_kl, integral, temp_int
  real*8 chi_nm, chi_kl
  part_1 = (0.0d0, 0.0d0)
  temp_int = (0.0d0, 0.0d0)

  chi_nm = mu(m, n)/rt(is-Hc)
  chi_kl = mu(m, n)/rt(is-Hc)
  if (n.ne.0 .and. k.ne.0) then
    temp_int = (n*k/chi_kl/chi_nm) * integrate(bessel_mult, 0.0d0, rt(is-Hc), 0.0001d0) * integrate(der_cossin_mult, 0.0d0, 6.2831853d0, 0.0001d0)
  end if
  integral =  integrate(der_bessel_mult, 0.0d0, rt(is-Hc), 0.0001d0) * integrate(cossin_mult, 0.0d0, 6.2831853d0, 0.00001d0) + temp_int
  E_kl = sqrt(zn(is-Hc, m, n) / abs(zn(is-Hc, m, n)) * chi_nm * chi_kl * conjg(zn(is-Hc, m, n)) / gam(is-Hc, m, n) / conjg(gam(is-Hc, m, n)) / integral)

  chi_nm = mu(l, k)/rt(is-Ec)
  chi_kl = mu(l, k)/rt(is-Ec)
  if (n.ne.0 .and. k.ne.0) then
    temp_int = (n*k/chi_kl/chi_nm) * integrate(bessel_mult, 0.0d0, rt(is-Ec), 0.0001d0) * integrate(der_cossin_mult, 0.0d0, 6.2831853d0, 0.0001d0)
  end if
  integral =  integrate(der_bessel_mult, 0.0d0, rt(is-Ec), 0.0001d0) * integrate(cossin_mult, 0.0d0, 6.2831853d0, 0.0001d0) + temp_int
  E_nm = sqrt(zn(is-Ec, l, k) / abs(zn(is-Ec, l, k)) * chi_nm * chi_kl * conjg(zn(is-Ec, l, k)) / gam(is-Ec, l, k) / conjg(gam(is-Ec, l, k)) / integral)

  chi_nm = mu(m, n)/rt(is-Hc)
  chi_kl = mu(l, k)/rt(is-Ec)
  const_part =  gam(is-Hc, m, n)*conjg(gam(is-Ec, l,k)) / conjg(zn(is-Ec, l, k)) / chi_nm / chi_kl
  part_2 = integrate(der_bessel_mult, 0.0d0, rt(is-Hc), 0.0001d0) * integrate(cossin_mult, 0.0d0, 6.2831853d0, 0.0001d0)
  if (n.ne.0 .and. k.ne.0) then
    part_1 = (n*k/chi_kl/chi_nm) * integrate(bessel_mult, 0.0d0, rt(is-Hc), 0.0001d0) * integrate(der_cossin_mult, 0.0d0, 6.2831853d0, 0.0001d0)
  end if

  pfunk = const_part * (part_1 + part_2) * E_nm * conjg(E_kl)
  if ((mod(m+l,2)).ne.0) then
    pfunk = -pfunk
  end if
  print *, n, m, k, l, pfunk
  return
  contains

  real*8 function bessel_mult(r)
    real*8 r
    bessel_mult = bessel_jn(n, chi_nm*r) * bessel_jn(k, chi_kl*r) / r
    return
  end function

  real*8 function der_bessel_mult(r)
    real*8 r
    der_bessel_mult = (bessel_jn(n-1, chi_nm*r) - bessel_jn(n+1, chi_nm*r)) * &
                      (bessel_jn(k-1, chi_kl*r) - bessel_jn(k+1, chi_kl*r)) / 4.0d0 * r
    return
  end function

  real*8 function cossin_mult(r)
    real*8 r
    cossin_mult = cos(n*r)*cos(k*r)
    return
  end function

  real*8 function der_cossin_mult(r)
    real*8 r
    der_cossin_mult = sin(n*r)*sin(k*r)
    return
  end function
end function


real function integrate(func, a, b, eps)
    implicit none
    real*8 :: a, b, eps
    real*8, external :: func
    real*8 :: curr, prev
    integrate = 0.0d0
    prev = a
    curr = a + eps
    do while (curr < b)
        integrate = integrate + (curr - prev) * 0.5d0 * &
        (func((prev+curr)*0.5d0 - (curr-prev)*0.2886751345948129d0) + &
        func((prev+curr)*0.5d0 + (curr-prev)*0.2886751345948129d0))
        prev = curr
        curr = curr + eps
    end do
    return
end function integrate


subroutine progonka
  integer is
  integer info

!    call dlincg (nkr, aa2(2,:,:,ina),nk,rabpinv, nkr)
  rabpinv(:,:,ina) = aa2(2,:,:,ina)
  call zgetrf(nkr,nkr,rabpinv,nkr,ipiv2,info)
  call zgetri(nkr, rabpinv,nkr,ipiv,temp2, nkr, info)
  call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),rabpinv,nkr,aa3(2,:,:,ina),nkr,(0.0d0,0.0d0),alphap(3,:,:),nkr)
  call zgemv('n',nkr,nkr,(1.0d0,0.0d0),rabpinv,nkr,ab(2,:,ina),1,(0.0d0,0.0d0),betap(3,:),1)

  do is=3,sk-1
    call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),aa1(is,:,:,ina),nkr,alphap(is,:,:),nkr,(0.0d0,0.0d0),rabp,nkr)

    rabp(:,:,ina)=-aa2(is,:,:,ina)-rabp(:,:,ina)
!        call dlincg (nkr, rabp,nkr,rabpinv, nkr)
    rabpinv = rabp
    call zgetrf(nkr,nkr,rabpinv,nkr,ipiv2,info)
    call zgetri(nkr, rabpinv,nkr,ipiv,temp2, nkr, info)
    call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),rabpinv,nkr,aa3(is,:,:,ina),nkr,(0.0d0,0.0d0),alphap(is+1,:,:),nkr)
    call zgemv('n',nkr,nkr,(1.0d0,0.0d0),aa1(is,:,:,ina),nkr,betap(is,:),1,(0.0d0,0.0d0),xrabp,1)
    xrabp=xrabp-ab(is,:,ina)
    call zgemv('n',nkr,nkr,(1.0d0,0.0d0),rabpinv,nkr,xrabp,1,(0.0d0,0.0d0),betap(is+1,:),1)
	end do

  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),aa1(sk,:,:,ina),nkr,alphap(sk,:,:),nkr,(0.0d0,0.0d0),rabp,nkr)
  rabp=-aa2(sk,:,:,:)-rabp
!    call dlincg (nkr, rabp,nkr,rabpinv, nkr)
  rabpinv = rabp
  call zgetrf(nkr,nkr,rabpinv,nkr,ipiv2,info)
  call zgetri(nkr, rabpinv,nkr,ipiv2,temp2, nkr, info)
  call zgemv('n',nkr,nkr,(1.0d0,0.0d0),aa1(sk,:,:,ina),nkr,betap(sk,:),1,(0.0d0,0.0d0),xrabp,1)
  xrabp=xrabp-ab(sk,:,ina)
  call zgemv('n',nkr,nkr,(1.0d0,0.0d0),rabpinv,nkr,xrabp,1,(0.0d0,0.0d0),betap(sk+1,:),1)
  xbplus(sk,:,ina)=betap(sk+1,:)
  do is=sk-1,2,-1
    call zgemv('n',nkr,nkr,(1.0d0,0.0d0),alphap(is+1,:,:),nkr,xbplus(is+1,:,ina),1,(0.0d0,0.0d0),xbplus(is,:,ina),1)
    xbplus(is,:,ina)=betap(is+1,:)+xbplus(is,:,ina)
	end do
	return
end subroutine progonka

end  module field_nes
