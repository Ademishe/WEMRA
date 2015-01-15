module field_nes
  use array_work
  use com_prog
  use omp_lib
  implicit none

  contains

subroutine init_arrays
  aa1(:,:,:,:)				  = (0.0d0, 0.0d0)
  aa2(:,:,:,:)				  = (0.0d0, 0.0d0)
  aa3(:,:,:,:)				  = (0.0d0, 0.0d0)
  fgammaplus(:,:,:,:)   = (0.0d0,0.0d0)
  fgammaminus(:,:,:,:)  = (0.0d0,0.0d0)
  d1plus(:,:,:,:)       = (0.0d0,0.0d0)
  d1minus(:,:,:,:)      = (0.0d0,0.0d0)
  d2plus(:,:,:,:)       = (0.0d0,0.0d0)
  d2minus(:,:,:,:)      = (0.0d0,0.0d0)
  dd1plus(:,:,:,:)      = (0.0d0,0.0d0)
  dd1minus(:,:,:,:)     = (0.0d0,0.0d0)
  dd2plus(:,:,:,:)      = (0.0d0,0.0d0)
  dd2minus(:,:,:,:)     = (0.0d0,0.0d0)
  b1plus(:,:,:,:)       = (0.0d0,0.0d0)
  b1minus(:,:,:,:)      = (0.0d0,0.0d0)
  b2plus(:,:,:,:)       = (0.0d0,0.0d0)
  b2minus(:,:,:,:)      = (0.0d0,0.0d0)
  bb1plus(:,:,:,:)      = (0.0d0,0.0d0)
  bb1minus(:,:,:,:)     = (0.0d0,0.0d0)
  bb2plus(:,:,:,:)      = (0.0d0,0.0d0)
  bb2minus(:,:,:,:)     = (0.0d0,0.0d0)
  uste1(:,:,:,:)        = (0.0d0,0.0d0)
  uste2(:,:,:,:)        = (0.0d0,0.0d0)
  usth2(:,:,:,:)        = (0.0d0,0.0d0)
  usth1(:,:,:,:)        = (0.0d0,0.0d0)
  rnplus(:,:,:)         = (0.0d0,0.0d0)
  rnminus(:,:,:)        = (0.0d0,0.0d0)
  rnnplus(:,:,:)        = (0.0d0,0.0d0)
  rnnminus(:,:,:)       = (0.0d0,0.0d0)
end subroutine

subroutine matrix_construct
	integer is,i,j,kluch_shiv, ina
  integer thread_num, thread_id
	integer info
	complex*16 fdplus,fdminus, fbminus, fbplus
	kluch_shiv=1
	aa1(:,:,:,:) = (0.0d0, 0.0d0)
	aa2(:,:,:,:) = (0.0d0, 0.0d0)
	aa3(:,:,:,:) = (0.0d0, 0.0d0)
	ab(:,:,:) = (0.0d0, 0.0d0)

  do ina = 0, nka
!$omp parallel default(private) shared(ee, ce, dt, ina, w0, gam0, kluch_shiv, dz, nkr, sk, nka, uste1, uste2, usth1, usth2, rt, gam, d1plus, d1minus, d2plus, d2minus, dd1plus, dd1minus, dd2plus, dd2minus, mu, b1plus, b1minus, b2plus, b2minus, bb1plus, bb1minus, bb2plus, bb2minus, zn, fgammaplus, fgammaminus, db1plus, db1minus, ddb1plus, ddb1minus, ddb2plus, ddb2minus, db2plus, db2minus, ddb1plusinv, ddb1minusinv, ddb2plusinv, ddb2minusinv, dddb1, dddb2, dddb1inv, dddb2inv, alfaplus, alfaminus, betaplus, betaminus)
!$omp do
    do is = 1, sk
      thread_id = omp_get_thread_num()
      thread_num = omp_get_num_threads()
      ! print *, "***matrix_construct***", is, sk, ina, thread_id, thread_num
      do i = 1, nkr
        do j = 1, nkr
          if(is.gt.1 .and. (rt(is-1)-rt(is)).ge.(0.001d0)) then !уменьшается
            uste1(is,i,j,ina) = (0.0d0,0.0d0)
            uste2(is,i,j,ina) = pfunk(is,i,j,1,0,ina)
            usth1(is,i,j,ina) = conjg(pfunk(is,j,i,1,0,ina))
            usth2(is,i,j,ina) = (0.0d0,0.0d0)
          end if
          if(is.gt.1 .and. (rt(is)-rt(is-1)).gt. (0.001d0)) then !увеличивается
            uste1(is,i,j,ina) = pfunk(is,i,j,0,1,ina)   !менял i и j местами
            uste2(is,i,j,ina) = (0.0d0,0.0d0)
            usth1(is,i,j,ina) = (0.0d0,0.0d0)
            usth2(is,i,j,ina) = conjg(pfunk(is,j,i,0,1,ina))
          end if
        end do !end do j
			  fbplus = (abs(gam(is,i,ina))+gam0*gam0/abs(gam(is,i,ina)))/w0
        fbminus = (abs(gam(is,i,ina))-gam0*gam0/abs(gam(is,i,ina)))/w0
			  fdplus=zn(is,i,ina)/abs(zn(is,i,ina))+conjg(zn(is,i,ina))/abs(zn(is,i,ina))
			  fdminus=-zn(is,i,ina)/abs(zn(is,i,ina))+conjg(zn(is,i,ina))/abs(zn(is,i,ina))
			  d1plus(is,i,i,ina) = (fexp(is,i,-1,-1,ina)-fexp(is,i,1,-1,ina))*  fdplus
			  d1minus(is,i,i,ina) = (fexp(is,i,-1,1,ina) -fexp(is,i,1,1,ina)) *  fdminus
			  d2plus(is,i,i,ina) = (fexp(is,i,1,1,ina)  -fexp(is,i,-1,1,ina))*(-fdminus)
			  d2minus(is,i,i,ina) = (fexp(is,i,1,-1,ina)-fexp(is,i,-1,-1,ina))*(-fdplus)
			  dd1plus(is,i,i,ina) = (fexp(is,i,-1,-1,ina)+fexp(is,i,1,-1,ina))*fdplus*dz(is)/2.0d0
			  dd1minus(is,i,i,ina) =(fexp(is,i,-1,1,ina)+fexp(is,i,1,1,ina))*fdminus*dz(is)/2.0d0
			  dd2plus(is,i,i,ina) = (fexp(is,i,1,1,ina)+fexp(is,i,-1,1,ina))*(-fdminus)*dz(is)/2.0d0
        dd2minus(is,i,i,ina) = (fexp(is,i,1,-1,ina)+fexp(is,i,-1,-1,ina))*(-fdplus)*dz(is)/2.0d0

        if (mu(i,ina)/rt(is) .le. gam0) then
					b1plus(is,i,i,ina)= -fbplus*dz(is)
					b1minus(is,i,i,ina)=-fbminus*funk1(is,i,ina)
					b2plus(is,i,i,ina)= -fbminus*funk1(is,i,ina)
					b2minus(is,i,i,ina)= -fbplus*dz(is)
					bb1plus(is,i,i,ina)= 0.0d0
					bb1minus(is,i,i,ina)=-fbminus*(funk1(is,i,ina)/(2.0d0*ce*gam(is,i,ina))-funk2(is,i,ina)*dz(is)/2.0d0)
					bb2plus(is,i,i,ina)= -fbminus*(funk2(is,i,ina)*dz(is)/2.0d0-funk1(is,i,ina)/(2.0d0*ce*gam(is,i,ina)))
					bb2minus(is,i,i,ina)= 0.0d0
        else
					b1plus(is,i,i,ina)=  -fbplus*funk1(is,i,ina)
					b1minus(is,i,i,ina)= -fbminus*dz(is)
					b2plus(is,i,i,ina)=  -fbminus*dz(is)
					b2minus(is,i,i,ina)= -fbplus*funk1(is,i,ina)
					bb1plus(is,i,i,ina)= -fbplus*(funk1(is,i,ina)/(2.0d0*ce*gam(is,i,ina))-funk2(is,i,ina)*dz(is)/2.0d0)
					bb1minus(is,i,i,ina)= 0.0d0
					bb2plus(is,i,i,ina)=  0.0d0
          bb2minus(is,i,i,ina)= -fbplus*(-funk1(is,i,ina)/(2.0d0*ce*gam(is,i,ina))+funk2(is,i,ina)*dz(is)/2.0d0)
        end if
        if (kluch_shiv.eq.1) then
          if(is.gt.1 .and. (rt(is-1)-rt(is)).gt. 0.001d0) then
						uste1(is,i,i,ina)=zn(is-1,i,ina)/abs(zn(is-1,i,ina))
            usth2(is,i,i,ina)=conjg(zn(is,i,ina))/abs(zn(is,i,ina))
  				end if
          if(is.gt.1 .and. (rt(is)-rt(is-1)).gt. 0.001d0) then
						uste2(is,i,i,ina)=zn(is,i,ina)/abs(zn(is,i,ina))
            usth1(is,i,i,ina)=conjg(zn(is-1,i,ina))/abs(zn(is-1,i,ina))
          end if
			  end if
        if(is.gt.1 .and. abs(rt(is)-rt(is-1)).le. 0.001) then
					uste1(is,i,i,ina)= 1.0d0
					uste2(is,i,i,ina)= 1.0d0
					usth1(is,i,i,ina)= 1.0d0
          usth2(is,i,i,ina)= 1.0d0
        end if
				fgammaplus(is,i,i,ina)=exp(-ce*gam(is,i,ina)*dz(is)/2.0d0)
        fgammaminus(is,i,i,ina)=exp(ce*gam(is,i,ina)*dz(is)/2.0d0)
      end do
		  db1plus(is,:,:,ina)=d1plus(is,:,:,ina)*dt/w0-b1plus(is,:,:,ina)
		  db1minus(is,:,:,ina)=d1minus(is,:,:,ina)*dt/w0-b1minus(is,:,:,ina)
		  ddb1plus(is,:,:,ina)=dd1plus(is,:,:,ina)*dt/w0-bb1plus(is,:,:,ina)
		  ddb1minus(is,:,:,ina)=dd1minus(is,:,:,ina)*dt/w0-bb1minus(is,:,:,ina)
		  db2plus(is,:,:,ina)=d2plus(is,:,:,ina)*dt/w0-b2plus(is,:,:,ina)
		  db2minus(is,:,:,ina)=d2minus(is,:,:,ina)*dt/w0-b2minus(is,:,:,ina)
		  ddb2plus(is,:,:,ina)=dd2plus(is,:,:,ina)*dt/w0-bb2plus(is,:,:,ina)
		  ddb2minus(is,:,:,ina)=dd2minus(is,:,:,ina)*dt/w0-bb2minus(is,:,:,ina)
!			call dlincg(nkr, ddb1plus(is,:,:,ina),nkr,ddb1plusinv(is,:,:,ina), nkr)
!			call dlincg(nkr, ddb2plus(is,:,:,ina),nkr,ddb2plusinv(is,:,:,ina), nkr)
!			call dlincg(nkr, ddb1minus(is,:,:,ina),nkr,ddb1minusinv(is,:,:,ina), nkr)
!			call dlincg(nkr, ddb2minus(is,:,:,ina),nkr,ddb2minusinv(is,:,:,ina), nkr)
		  ddb1plusinv(is,:,:,ina) = ddb1plus(is,:,:,ina)
		  call zgetrf(nkr,nkr,ddb1plusinv(is,:,:,ina),nkr,ipiv,info)
		  call zgetri(nkr, ddb1plusinv(is,:,:,ina),nkr,ipiv,temp, nkr, info)
		  ddb2plusinv(is,:,:,ina) = ddb2plus(is,:,:,ina)
		  call zgetrf(nkr,nkr,ddb2plusinv(is,:,:,ina),nkr,ipiv,info)
		  call zgetri(nkr, ddb2plusinv(is,:,:,ina),nkr,ipiv,temp, nkr, info)
		  ddb1minusinv(is,:,:,ina) = ddb1minus(is,:,:,ina)
		  call zgetrf(nkr,nkr,ddb1minusinv(is,:,:,ina),nkr,ipiv,info)
		  call zgetri(nkr, ddb1minusinv(is,:,:,ina),nkr,ipiv,temp, nkr, info)
		  ddb2minusinv(is,:,:,ina) = ddb2minus(is,:,:,ina)
		  call zgetrf(nkr,nkr,ddb2minusinv(is,:,:,ina),nkr,ipiv,info)
		  call zgetri(nkr, ddb2minusinv(is,:,:,ina),nkr,ipiv,temp, nkr, info)
		  call zgemm('n','n',nkr,nkr,nkr,( 1.0d0, 0.0d0),ddb2plusinv(is,:,:,ina),nkr,ddb1plus(is,:,:,ina),nkr,(0.0d0,0.0d0),dddb1(is,:,:,ina),nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(-1.0d0, 0.0d0),ddb2minusinv(is,:,:,ina),nkr,ddb1minus(is,:,:,ina),nkr,(1.0d0,0.0d0),dddb1(is,:,:,ina),nkr)
		  call zgemm('n','n',nkr,nkr,nkr,( 1.0d0, 0.0d0),ddb1plusinv(is,:,:,ina),nkr,ddb2plus(is,:,:,ina),nkr,(0.0d0,0.0d0),dddb2(is,:,:,ina),nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(-1.0d0, 0.0d0),ddb1minusinv(is,:,:,ina),nkr,ddb2minus(is,:,:,ina),nkr,(1.0d0,0.0d0),dddb2(is,:,:,ina),nkr)
!			call dlincg (nkr, dddb1(is,:,:,ina),nkr,dddb1inv(is,:,:,ina), nkr)
!			call dlincg (nkr, dddb2(is,:,:,ina),nkr,dddb2inv(is,:,:,ina), nkr)
		  dddb1inv(is,:,:,ina) = dddb1(is,:,:,ina)
		  call zgetrf(nkr,nkr,dddb1inv(is,:,:,ina),nkr,ipiv,info)
		  call zgetri(nkr, dddb1inv(is,:,:,ina),nkr,ipiv,temp, nkr, info)
		  dddb2inv(is,:,:,ina) = dddb2(is,:,:,ina)
		  call zgetrf(nkr,nkr,dddb2inv(is,:,:,ina),nkr,ipiv,info)
		  call zgetri(nkr, dddb2inv(is,:,:,ina),nkr,ipiv,temp, nkr, info)

		  call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),ddb2plusinv(is,:,:,ina),nkr,db1plus(is,:,:,ina),nkr,(0.0d0,0.0d0),rab1,nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),ddb2minusinv(is,:,:,ina),nkr,db1minus(is,:,:,ina),nkr,(1.0d0,0.0d0),rab1,nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),dddb1inv(is,:,:,ina),nkr,rab1,nkr,(0.0d0,0.0d0),alfaplus(is,:,:,ina),nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),ddb1plusinv(is,:,:,ina),nkr,db1plus(is,:,:,ina),nkr,(0.0d0,0.0d0),rab1,nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),ddb1minusinv(is,:,:,ina),nkr,db1minus(is,:,:,ina),nkr,(1.0d0,0.0d0),rab1,nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),dddb2inv(is,:,:,ina),nkr,rab1,nkr,(0.0d0,0.0d0),alfaminus(is,:,:,ina),nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),ddb2plusinv(is,:,:,ina),nkr,db2plus(is,:,:,ina),nkr,(0.0d0,0.0d0),rab1,nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),ddb2minusinv(is,:,:,ina),nkr,db2minus(is,:,:,ina),nkr,(1.0d0,0.0d0),rab1,nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),dddb1inv(is,:,:,ina),nkr,rab1,nkr,(0.0d0,0.0d0),betaplus(is,:,:,ina),nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),ddb1plusinv(is,:,:,ina),nkr,db2plus(is,:,:,ina),nkr,(0.0d0,0.0d0),rab1,nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),ddb1minusinv(is,:,:,ina),nkr,db2minus(is,:,:,ina),nkr,(1.0d0,0.0d0),rab1,nkr)
      call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),dddb2inv(is,:,:,ina),nkr,rab1,nkr,(0.0d0,0.0d0),betaminus(is,:,:,ina),nkr)
    end do
!$omp end do
!$omp end parallel

    do is=2,sk
		  call zgemm('n','n',nkr,nkr,nkr,dz(is-1)/2.0d0*ee,fgammaplus(is-1,:,:,ina),nkr,alfaplus(is-1,:,:,ina),nkr,(0.0d0,0.0d0),rab1,nkr)
		  call zgemm('n','n',nkr,nkr,nkr,dz(is-1)/2.0d0*ee,fgammaminus(is-1,:,:,ina),nkr,alfaminus(is-1,:,:,ina),nkr,(1.0d0,0.0d0),rab1,nkr)
		  rab1=rab1+fgammaplus(is-1,:,:,ina)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),uste1(is,:,:,ina),nkr,rab1,nkr,(0.0d0,0.0d0),hie1(is,:,:,ina),nkr)
		  call zgemm('n','n',nkr,nkr,nkr,dz(is-1)/2.0d0*ee,fgammaplus(is-1,:,:,ina),nkr,alfaplus(is-1,:,:,ina),nkr,(0.0d0,0.0d0),rab1,nkr)
		  call zgemm('n','n',nkr,nkr,nkr,-dz(is-1)/2.0d0*ee,fgammaminus(is-1,:,:,ina),nkr,alfaminus(is-1,:,:,ina),nkr,(1.0d0,0.0d0),rab1,nkr)
		  rab1=rab1+fgammaplus(is-1,:,:,ina)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),usth1(is,:,:,ina),nkr,rab1,nkr,(0.0d0,0.0d0),hih1(is,:,:,ina),nkr)
		  call zgemm('n','n',nkr,nkr,nkr,-dz(is)/2.0d0*ee,fgammaminus(is,:,:,ina),nkr,alfaplus(is,:,:,ina),nkr,(0.0d0,0.0d0),rab1,nkr)
		  call zgemm('n','n',nkr,nkr,nkr,-dz(is)/2.0d0*ee,fgammaplus(is,:,:,ina),nkr,alfaminus(is,:,:,ina),nkr,(1.0d0,0.0d0),rab1,nkr)
		  rab1=rab1+fgammaminus(is,:,:,ina)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),uste2(is,:,:,ina),nkr,rab1,nkr,(0.0d0,0.0d0),hie2(is,:,:,ina),nkr)

		  call zgemm('n','n',nkr,nkr,nkr,-dz(is)/2.0d0*ee,fgammaminus(is,:,:,ina),nkr,alfaplus(is,:,:,ina),nkr,(0.0d0,0.0d0),rab1,nkr)
		  call zgemm('n','n',nkr,nkr,nkr,dz(is)/2.0d0*ee,fgammaplus(is,:,:,ina),nkr,alfaminus(is,:,:,ina),nkr,(1.0d0,0.0d0),rab1,nkr)
		  rab1=rab1+fgammaminus(is,:,:,ina)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),usth2(is,:,:,ina),nkr,rab1,nkr,(0.0d0,0.0d0),hih2(is,:,:,ina),nkr)

		  call zgemm('n','n',nkr,nkr,nkr,dz(is-1)/2.0d0*ee,fgammaplus(is-1,:,:,ina),nkr,betaplus(is-1,:,:,ina),nkr,(0.0d0,0.0d0),rab1,nkr)
		  call zgemm('n','n',nkr,nkr,nkr,dz(is-1)/2.0d0*ee,fgammaminus(is-1,:,:,ina),nkr,betaminus(is-1,:,:,ina),nkr,(1.0d0,0.0d0),rab1,nkr)
		  rab1=rab1+fgammaminus(is-1,:,:,ina)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),uste1(is,:,:,ina),nkr,rab1,nkr,(0.0d0,0.0d0),psie1(is,:,:,ina),nkr)
		  call zgemm('n','n',nkr,nkr,nkr,dz(is-1)/2.0d0*ee,fgammaplus(is-1,:,:,ina),nkr,betaplus(is-1,:,:,ina),nkr,(0.0d0,0.0d0),rab1,nkr)
		  call zgemm('n','n',nkr,nkr,nkr,-dz(is-1)/2.0d0*ee,fgammaminus(is-1,:,:,ina),nkr,betaminus(is-1,:,:,ina),nkr,(1.0d0,0.0d0),rab1,nkr)
		  rab1=rab1-fgammaminus(is-1,:,:,ina)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),usth1(is,:,:,ina),nkr,rab1,nkr,(0.0d0,0.0d0),psih1(is,:,:,ina),nkr)
		  call zgemm('n','n',nkr,nkr,nkr,-dz(is)/2.0d0*ee,fgammaminus(is,:,:,ina),nkr,betaplus(is,:,:,ina),nkr,(0.0d0,0.0d0),rab1,nkr)
		  call zgemm('n','n',nkr,nkr,nkr,-dz(is)/2.0d0*ee,fgammaplus(is,:,:,ina),nkr,betaminus(is,:,:,ina),nkr,(1.0d0,0.0d0),rab1,nkr)
		  rab1=rab1+fgammaplus(is,:,:,ina)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),uste2(is,:,:,ina),nkr,rab1,nkr,(0.0d0,0.0d0),psie2(is,:,:,ina),nkr)
		  call zgemm('n','n',nkr,nkr,nkr,-dz(is)/2.0d0*ee,fgammaminus(is,:,:,ina),nkr,betaplus(is,:,:,ina),nkr,(0.0d0,0.0d0),rab1,nkr)
		  call zgemm('n','n',nkr,nkr,nkr,dz(is)/2.0d0*ee,fgammaplus(is,:,:,ina),nkr,betaminus(is,:,:,ina),nkr,(1.0d0,0.0d0),rab1,nkr)
		  rab1=rab1-fgammaplus(is,:,:,ina)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),usth2(is,:,:,ina),nkr,rab1,nkr,(0.0d0,0.0d0),psih2(is,:,:,ina),nkr)

!			call dlincg (nkr, psie1(is,:,:,ina),nkr,psie1inv(is,:,:,ina), nkr)
!			call dlincg (nkr, psie2(is,:,:,ina),nkr,psie2inv(is,:,:,ina), nkr)
!			call dlincg (nkr, psih1(is,:,:,ina),nkr,psih1inv(is,:,:,ina), nkr)
!			call dlincg (nkr, psih2(is,:,:,ina),nkr,psih2inv(is,:,:,ina), nkr)
		  psie1inv(is,:,:,ina) = psie1(is,:,:,ina)
		  call zgetrf(nkr,nkr,psie1inv(is,:,:,ina),nkr,ipiv,info)
		  call zgetri(nkr, psie1inv(is,:,:,ina),nkr,ipiv,temp, nkr, info)
		  psie2inv(is,:,:,ina) = psie2(is,:,:,ina)
		  call zgetrf(nkr,nkr,psie2inv(is,:,:,ina),nkr,ipiv,info)
		  call zgetri(nkr, psie2inv(is,:,:,ina),nkr,ipiv,temp, nkr, info)
		  psih1inv(is,:,:,ina) = psih1(is,:,:,ina)
		  call zgetrf(nkr,nkr,psih1inv(is,:,:,ina),nkr,ipiv,info)
		  call zgetri(nkr, psih1inv(is,:,:,ina),nkr,ipiv,temp, nkr, info)
  		psih2inv(is,:,:,ina) = psih2(is,:,:,ina)
      call zgetrf(nkr,nkr,psih2inv(is,:,:,ina),nkr,ipiv,info)
      call zgetri(nkr, psih2inv(is,:,:,ina),nkr,ipiv,temp, nkr, info)

      call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),psie1inv(is,:,:,ina),nkr,psie2(is,:,:,ina),nkr,(0.0d0,0.0d0),psipsi1(is,:,:,ina),nkr)
      call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psih1inv(is,:,:,ina),nkr,psih2(is,:,:,ina),nkr,(1.0d0,0.0d0),psipsi1(is,:,:,ina),nkr)
      call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),psie1inv(is,:,:,ina),nkr,hie1(is,:,:,ina),nkr,(0.0d0,0.0d0),psihi1(is,:,:,ina),nkr)
      call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psih1inv(is,:,:,ina),nkr,hih1(is,:,:,ina),nkr,(1.0d0,0.0d0),psihi1(is,:,:,ina),nkr)
      call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psie1inv(is,:,:,ina),nkr,hie2(is,:,:,ina),nkr,(0.0d0,0.0d0),psihi2(is,:,:,ina),nkr)
      call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),psih1inv(is,:,:,ina),nkr,hih2(is,:,:,ina),nkr,(1.0d0,0.0d0),psihi2(is,:,:,ina),nkr)
!			call dlincg (nkr, psipsi1(is,:,:,ina),nkr,psipsi1inv(is,:,:,ina), nkr)
      psipsi1inv(is,:,:,ina) = psipsi1(is,:,:,ina)
      call zgetrf(nkr,nkr,psipsi1inv(is,:,:,ina),nkr,ipiv,info)
      call zgetri(nkr, psipsi1inv(is,:,:,ina),nkr,ipiv,temp, nkr, info)
    end do
    do is=1,sk-1
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),psie2inv(is+1,:,:,ina),nkr,psie1(is+1,:,:,ina),nkr,(0.0d0,0.0d0),psipsi2(is,:,:,ina),nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psih2inv(is+1,:,:,ina),nkr,psih1(is+1,:,:,ina),nkr,(1.0d0,0.0d0),psipsi2(is,:,:,ina),nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),psie2inv(is+1,:,:,ina),nkr,hie2(is+1,:,:,ina),nkr,(0.0d0,0.0d0),psihi3(is,:,:,ina),nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psih2inv(is+1,:,:,ina),nkr,hih2(is+1,:,:,ina),nkr,(1.0d0,0.0d0),psihi3(is,:,:,ina),nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psie2inv(is+1,:,:,ina),nkr,hie1(is+1,:,:,ina),nkr,(0.0d0,0.0d0),psihi4(is,:,:,ina),nkr)
		  call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),psih2inv(is+1,:,:,ina),nkr,hih1(is+1,:,:,ina),nkr,(1.0d0,0.0d0),psihi4(is,:,:,ina),nkr)
!			call dlincg(nkr, psipsi2(is,:,:,ina),nkr,psipsi2inv(is,:,:,ina), nkr)
      psipsi2inv(is,:,:,ina) = psipsi2(is,:,:,ina)
      call zgetrf(nkr,nkr,psipsi2inv(is,:,:,ina),nkr,ipiv,info)
      call zgetri(nkr, psipsi2inv(is,:,:,ina),nkr,ipiv,temp, nkr, info)
    end do
    do is=2,sk-1
  		call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),psipsi1inv(is,:,:,ina),nkr,psihi1(is,:,:,ina),nkr,(0.0d0,0.0d0),aa1(is,:,:,ina),nkr)
      call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),psipsi1inv(is,:,:,ina),nkr,psihi2(is,:,:,ina),nkr,(0.0d0,0.0d0),aa2(is,:,:,ina),nkr)
      call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psipsi2inv(is,:,:,ina),nkr,psihi4(is,:,:,ina),nkr,(1.0d0,0.0d0),aa2(is,:,:,ina),nkr)
      call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psipsi2inv(is,:,:,ina),nkr,psihi3(is,:,:,ina),nkr,(0.0d0,0.0d0),aa3(is,:,:,ina),nkr)
    end do

!		call dlincg(nkr, hie2(sk,:,:,ina),nkr,hie2inv(sk,:,:,ina), nkr)
!		call dlincg(nkr, hih2(sk,:,:,ina),nkr,hih2inv(sk,:,:,ina), nkr)
	  hie2inv(sk,:,:,ina) = hie2(sk,:,:,ina)
	  call zgetrf(nkr,nkr,hie2inv(sk,:,:,ina),nkr,ipiv,info)
    call zgetri(nkr, hie2inv(sk,:,:,ina),nkr,ipiv,temp, nkr, info)
    hih2inv(sk,:,:,ina) = hih2(sk,:,:,ina)
    call zgetrf(nkr,nkr,hih2inv(sk,:,:,ina),nkr,ipiv,info)
    call zgetri(nkr, hih2inv(sk,:,:,ina),nkr,ipiv,temp, nkr, info)

    aa1(sk,:,:,ina)=psihi1(sk,:,:,ina)
    aa2(sk,:,:,ina)=psihi2(sk,:,:,ina)
  end do
	return
end subroutine  matrix_construct



subroutine matrix_dynamic
  integer is, ina
  do ina = 0, nka
    do is=1,sk
      call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),b1plus(is,:,:,ina),nkr,xnplus(is,:,ina),1,(0.0d0,0.0d0),rnplus(is,:,ina),1)
		  call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),bb1plus(is,:,:,ina),nkr,dxnplus(is,:,ina),1,(1.0d0,0.0d0),rnplus(is,:,ina),1)
		  call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),b2plus(is,:,:,ina),nkr,xnminus(is,:,ina),1,(1.0d0,0.0d0),rnplus(is,:,ina),1)
		  call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),bb2plus(is,:,:,ina),nkr,dxnminus(is,:,ina),1,(1.0d0,0.0d0),rnplus(is,:,ina),1)
	    rnplus(is,:,ina)=rnplus(is,:,ina)+etaplus(is,:,ina)

		  call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),b1minus(is,:,:,ina),nkr,xnplus(is,:,ina),1,(0.0d0,0.0d0),rnminus(is,:,ina),1)
		  call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),bb1minus(is,:,:,ina),nkr,dxnplus(is,:,ina),1,(1.0d0,0.0d0),rnminus(is,:,ina),1)
		  call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),b2minus(is,:,:,ina),nkr,xnminus(is,:,ina),1,(1.0d0,0.0d0),rnminus(is,:,ina),1)
		  call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),bb2minus(is,:,:,ina),nkr,dxnminus(is,:,ina),1,(1.0d0,0.0d0),rnminus(is,:,ina),1)
		  rnminus(is,:,ina)=rnminus(is,:,ina)+etaminus(is,:,ina)

		  call zgemv('n',nkr,nkr,(1.0d0,0.0d0),ddb2plusinv(is,:,:,ina),nkr,rnplus(is,:,ina),1,(0.0d0,0.0d0),rab1,1)
		  call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),ddb2minusinv(is,:,:,ina),nkr,rnminus(is,:,ina),1,(1.0d0,0.0d0),rab1,1)
		  call zgemv('n',nkr,nkr,(1.0d0,0.0d0),dddb1inv(is,:,:,ina),nkr,rab1,1,(0.0d0,0.0d0),rnnplus(is,:,ina),1)
		  call zgemv('n',nkr,nkr,(1.0d0,0.0d0),ddb1plusinv(is,:,:,ina),nkr,rnplus(is,:,ina),1,(0.0d0,0.0d0),rab1,1)
		  call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),ddb1minusinv(is,:,:,ina),nkr,rnminus(is,:,ina),1,(1.0d0,0.0d0),rab1,1)
		  call zgemv('n',nkr,nkr,(1.0d0,0.0d0),dddb2inv(is,:,:,ina),nkr,rab1,1,(0.0d0,0.0d0),rnnminus(is,:,ina),1)
    end do

	  do is=2,sk
      call zgemv('n',nkr,nkr,dz(is-1)/2.0d0*ee,fgammaminus(is-1,:,:,ina),nkr,rnnminus(is-1,:,ina),1,(0.0d0,0.0d0),xrab,1)
		  call zgemv('n',nkr,nkr,dz(is-1)/2.0d0*ee,fgammaplus(is-1,:,:,ina),nkr,rnnplus(is-1,:,ina),1,(1.0d0,0.0d0),xrab,1)
		  call zgemv('n',nkr,nkr,(1.0d0,0.0d0),uste1(is,:,:,ina),nkr,xrab,1,(0.0d0,0.0d0),brne1(is,:,ina),1)
		  call zgemv('n',nkr,nkr,-dz(is)/2.0d0*ee,fgammaminus(is,:,:,ina),nkr,rnnplus(is,:,ina),1,(0.0d0,0.0d0),xrab,1)
	    call zgemv('n',nkr,nkr,-dz(is)/2.0d0*(1.0d0,0.0d0),fgammaplus(is,:,:,ina),nkr,rnnminus(is,:,ina),1,(1.0d0,0.0d0),xrab,1)
		  call zgemv('n',nkr,nkr,(1.0d0,0.0d0),uste2(is,:,:,ina),nkr,xrab,1,(0.0d0,0.0d0),brne2(is,:,ina),1)
		  call zgemv('n',nkr,nkr,-dz(is-1)/2.0d0*(1.0d0,0.0d0),fgammaminus(is-1,:,:,ina),nkr,rnnminus(is-1,:,ina),1,(0.0d0,0.0d0),xrab,1)
		  call zgemv('n',nkr,nkr,dz(is-1)/2.0d0*(1.0d0,0.0d0),fgammaplus(is-1,:,:,ina),nkr,rnnplus(is-1,:,ina),1,(1.0d0,0.0d0),xrab,1)
		  call zgemv('n',nkr,nkr,(1.0d0,0.0d0),usth1(is,:,:,ina),nkr,xrab,1,(0.0d0,0.0d0),brnh1(is,:,ina),1)
		  call zgemv('n',nkr,nkr,-dz(is)/2.0d0*(1.0d0,0.0d0),fgammaminus(is,:,:,ina),nkr,rnnplus(is,:,ina),1,(0.0d0,0.0d0),xrab,1)
		  call zgemv('n',nkr,nkr,dz(is)/2.0*ee,fgammaplus(is,:,:,ina),nkr,rnnminus(is,:,ina),1,(1.0d0,0.0d0),xrab,1)
		  call zgemv('n',nkr,nkr,(1.0d0,0.0d0),usth2(is,:,:,ina),nkr,xrab,1,(0.0d0,0.0d0),brnh2(is,:,ina),1)

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

	  do is=2,sk-1
      call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),psipsi1inv(is,:,:,ina),nkr,psibrn1(is,:,ina),1,(0.0d0,0.0d0),ab(is,:,ina),1)
		  call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psipsi2inv(is,:,:,ina),nkr,psibrn2(is,:,ina),1,(1.0d0,0.0d0),ab(is,:,ina),1)
	  end do
	  call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psipsi2inv(2,:,:,ina),nkr,psibrn2(2,:,ina),1,(0.0d0,0.0d0),ab(2,:,ina),1)
	  call zgemv('n',nkr,nkr,(-1.0d0,0.0d0),psipsi1inv(2,:,:,ina),nkr,psibrn1(2,:,ina),1,(1.0d0,0.0d0),ab(2,:,ina),1)

	  call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),psipsi1inv(2,:,:,ina),nkr,psihi1(2,:,:,ina),nkr,(0.0d0,0.0d0),rab1,nkr)
	  call zgemv('n',nkr,nkr,(1.0d0,0.0d0),rab1,nkr,bplus0(:,ina),1,(1.0d0,0.0d0),ab(2,:,ina),1)

	  call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psipsi1(sk,:,:,ina),nkr,bminussk(:,ina),1,(0.0d0,0.0d0),ab(sk,:,ina),1)
	  ab(sk,:,ina)=ab(sk,:,ina)-psibrn1(sk,:,ina)
  end do ! end ina

  return
end subroutine matrix_dynamic


subroutine field_calc
  integer is, ina
  do ina = 0, nka
    do is=2,sk-1
			call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psihi3(is,:,:,ina),nkr,xbplus(is+1,:,ina),1,(0.0d0,0.0d0),xrab,1)
			call zgemv('n',nkr,nkr,(1.0d0,0.0d0),psihi4(is,:,:,ina),nkr,xbplus(is,:,ina),1,(1.0d0,0.0d0),xrab,1)
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
  end do !end ina
  return
end subroutine field_calc

subroutine progonka
  integer info, is, ina
  do ina = 0, nka
    betap(:,:) = (0.0d0, 0.0d0)
    alphap(:,:,:) = (0.0d0, 0.0d0)
    rabp(:,:) = (0.0d0, 0.0d0)
!    call dlincg (nkr, aa2(2,:,:,ina),nkr,rabpinv, nkr)
    rabpinv(:,:) = aa2(2,:,:,ina)
    call zgetrf(nkr,nkr,rabpinv(:,:),nkr,ipiv2,info)
    call zgetri(nkr, rabpinv(:,:),nkr,ipiv2,temp2, nkr, info)
    call zgemm('n','n',nkr,nkr,nkr,(-1.0d0,0.0d0),rabpinv(:,:),nkr,aa3(2,:,:,ina),nkr,(0.0d0,0.0d0),alphap(3,:,:),nkr)
    call zgemv('n',nkr,nkr,(1.0d0,0.0d0),rabpinv(:,:),nkr,ab(2,:,ina),1,(0.0d0,0.0d0),betap(3,:),1)

    do is=3,sk-1
      call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),aa1(is,:,:,ina),nkr,alphap(is,:,:),nkr,(0.0d0,0.0d0),rabp(:,:),nkr)
      rabp(:,:)=-aa2(is,:,:,ina)-rabp(:,:)
!        call dlincg (nkr, rabp,nkr,rabpinv, nkr)
      rabpinv(:,:) = rabp(:,:)
      call zgetrf(nkr,nkr,rabpinv(:,:),nkr,ipiv2,info)
      call zgetri(nkr, rabpinv(:,:),nkr,ipiv2,temp2, nkr, info)
      call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),rabpinv(:,:),nkr,aa3(is,:,:,ina),nkr,(0.0d0,0.0d0),alphap(is+1,:,:),nkr)
      call zgemv('n',nkr,nkr,(1.0d0,0.0d0),aa1(is,:,:,ina),nkr,betap(is,:),1,(0.0d0,0.0d0),xrabp(:),1)
      xrabp(:)=xrabp(:)-ab(is,:,ina)
      call zgemv('n',nkr,nkr,(1.0d0,0.0d0),rabpinv(:,:),nkr,xrabp(:),1,(0.0d0,0.0d0),betap(is+1,:),1)
    end do

    call zgemm('n','n',nkr,nkr,nkr,(1.0d0,0.0d0),aa1(sk,:,:,ina),nkr,alphap(sk,:,:),nkr,(0.0d0,0.0d0),rabp(:,:),nkr)
    rabp=-aa2(sk,:,:,ina)-rabp
!    call dlincg (nkr, rabp,nkr,rabpinv, nkr)
    rabpinv(:,:) = rabp(:,:)
    call zgetrf(nkr,nkr,rabpinv(:,:),nkr,ipiv2,info)
    call zgetri(nkr, rabpinv(:,:),nkr,ipiv2,temp2, nkr, info)
    call zgemv('n',nkr,nkr,(1.0d0,0.0d0),aa1(sk,:,:,ina),nkr,betap(sk,:),1,(0.0d0,0.0d0),xrabp(:),1)
    xrabp(:) = xrabp(:) - ab(sk,:,ina)
    call zgemv('n',nkr,nkr,(1.0d0,0.0d0),rabpinv(:,:),nkr,xrabp(:),1,(0.0d0,0.0d0),betap(sk+1,:),1)
    xbplus(sk,:,ina)=betap(sk+1,:)
    do is=sk-1,2,-1
      call zgemv('n',nkr,nkr,(1.0d0,0.0d0),alphap(is+1,:,:),nkr,xbplus(is+1,:,ina),1,(0.0d0,0.0d0),xbplus(is,:,ina),1)
      xbplus(is,:,ina)=betap(is+1,:)+xbplus(is,:,ina)
    end do
  end do ! end ina
  return
end subroutine progonka

double complex function fexp(is,i,k1,k2,ina)
	integer is,i,k1,k2,ina
  fexp=exp(ce*k1*(gam(is,i,ina)+k2*conjg(gam(is,i,ina)))*(0.5d0*dz(is)))
end function fexp


double complex function funk1(is,inr,ina)
	integer is,inr,ina
  funk1=(exp(ce*gam(is,inr,ina)*dz(is))-exp(-ce*gam(is,inr,ina)*dz(is)))/(ce*2*gam(is,inr,ina))
end function funk1


double complex function funk2(is,inr,ina)
  integer is,inr,ina
  funk2=(exp(ce*gam(is,inr,ina)*dz(is))+exp(-ce*gam(is,inr,ina)*dz(is)))/(ce*2*gam(is,inr,ina))
end function funk2

double complex function pfunk(is,l,m,Ec,Hc,nk_index)
  integer is, m, l, n, k, Ec, Hc, nk_index
  complex*16 const_part, part_1, part_2, E_nm, E_kl, integral, temp_int
  real*8 chi_nm, chi_kl, accuracy
  part_1 = (0.0d0, 0.0d0)
  temp_int = (0.0d0, 0.0d0)
  n = nk_index
  k = nk_index
  accuracy = 0.001d0

  chi_nm = mu(m, n)/rt(is-Hc)
  chi_kl = mu(m, n)/rt(is-Hc)
  if (n.ne.0 .and. k.ne.0) then
    temp_int = (n*k/chi_kl/chi_nm) * integrate(bessel_mult, 0.0d0, rt(is-Hc), accuracy) * integrate(der_cossin_mult, 0.0d0, 6.2831853d0, accuracy)
  end if
  integral =  integrate(der_bessel_mult, 0.0d0, rt(is-Hc), accuracy) * integrate(cossin_mult, 0.0d0, 6.2831853d0, accuracy) + temp_int
  E_kl = sqrt(zn(is-Hc, m, n) / abs(zn(is-Hc, m, n)) * chi_nm * chi_kl * conjg(zn(is-Hc, m, n)) / gam(is-Hc, m, n) / conjg(gam(is-Hc, m, n)) / integral)

  chi_nm = mu(l, k)/rt(is-Ec)
  chi_kl = mu(l, k)/rt(is-Ec)
  if (n.ne.0 .and. k.ne.0) then
    temp_int = (n*k/chi_kl/chi_nm) * integrate(bessel_mult, 0.0d0, rt(is-Ec), accuracy) * integrate(der_cossin_mult, 0.0d0, 6.2831853d0, accuracy)
  end if
  integral =  integrate(der_bessel_mult, 0.0d0, rt(is-Ec), accuracy) * integrate(cossin_mult, 0.0d0, 6.2831853d0, accuracy) + temp_int
  E_nm = sqrt(zn(is-Ec, l, k) / abs(zn(is-Ec, l, k)) * chi_nm * chi_kl * conjg(zn(is-Ec, l, k)) / gam(is-Ec, l, k) / conjg(gam(is-Ec, l, k)) / integral)

  chi_nm = mu(m, n)/rt(is-Hc)
  chi_kl = mu(l, k)/rt(is-Ec)
  const_part =  gam(is-Hc, m, n)*conjg(gam(is-Ec, l,k)) / conjg(zn(is-Ec, l, k)) / chi_nm / chi_kl
  part_2 = integrate(der_bessel_mult, 0.0d0, rt(is-Hc), accuracy) * integrate(cossin_mult, 0.0d0, 6.2831853d0, accuracy)
  if (n.ne.0 .and. k.ne.0) then
    part_1 = (n*k/chi_kl/chi_nm) * integrate(bessel_mult, 0.0d0, rt(is-Hc), accuracy) * integrate(der_cossin_mult, 0.0d0, 6.2831853d0, accuracy)
  end if

  pfunk = const_part * (part_1 + part_2) * E_nm * conjg(E_kl)
  if (mod(m+l,2).ne.0) then
    pfunk = -pfunk
  end if
  ! print *, pfunk, n, m, k, l
  return
  contains

  real*8 function bessel_mult(r, arg1)
    real*8 r, arg1
    bessel_mult = bessel_jn(n, chi_nm*r) * bessel_jn(k, chi_kl*r) / r
    return
  end function

  real*8 function der_bessel_mult(r, arg1)
    real*8 r, arg1
    der_bessel_mult = (bessel_jn(n-1, chi_nm*r) - bessel_jn(n+1, chi_nm*r)) * &
                      (bessel_jn(k-1, chi_kl*r) - bessel_jn(k+1, chi_kl*r)) / 4.0d0 * r
    return
  end function

  real*8 function cossin_mult(r, arg1)
    real*8 r, arg1
    cossin_mult = cos(n*r)*cos(k*r)
    return
  end function

  real*8 function der_cossin_mult(r, arg1)
    real*8 r, arg1
    der_cossin_mult = sin(n*r)*sin(k*r)
    return
  end function
end function pfunk

end  module field_nes
