module com_prog
use array_work
contains


subroutine rk4 ( t0, u0, dt, f, u )
  implicit none

  real ( kind = 8 ) dt
  external f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) t0
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) u
  real ( kind = 8 ) u0
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) u3
  call f ( t0, u0, f0 )

  t1 = t0 + dt / 2.0D+00
  u1 = u0 + dt * f0 / 2.0D+00
  call f ( t1, u1, f1 )

  t2 = t0 + dt / 2.0D+00
  u2 = u0 + dt * f1 / 2.0D+00
  call f ( t2, u2, f2 )

  t3 = t0 + dt
  u3 = u0 + dt * f2
  call f ( t3, u3, f3 )

  u = u0 + dt * ( f0 + 2.0D+00 * f1 + 2.0D+00 * f2 + f3 ) / 6.0D+00
  return
end

subroutine drkgs(prmt, ndim, ihlf, ibeam)
  implicit none
  dimension a(4),b(4),c(4),prmt(5)
  double precision prmt,a,b,c,x,xend,h,aj,bj,cj,r1,r2,delt,ttau0

  integer ndim,ihlf,i,irec,istep,iend,itest,j,imod,itemp,ibeam
	ttau0=prmt(1)
	itemp=1

  do 1 i=1,ndim
    aux(8,i)=.066666666666666667d0*dery(i)
1   continue

    x=prmt(1)
    xend=prmt(2)
    h=prmt(3)
    prmt(5)=0.d0

    call fct(x,ibeam)

    if(h*(xend-x)) 38,37,2
2   a(1)=.5d0
    a(2)=.29289321881345248d0
    a(3)=1.7071067811865475d0
    a(4)=.16666666666666667d0
    b(1)=2.d0
    b(2)=1.d0
    b(3)=1.d0
    b(4)=2.d0
    c(1)=.5d0
    c(2)=.29289321881345248d0
    c(3)=1.7071067811865475d0
    c(4)=.5d0

    do 3 i=1,ndim
      aux(1,i)=yarray(i)
      aux(2,i)=dery(i)
      aux(3,i)=0.d0
      aux(6,i)=0.d0
3     continue
      irec=0
      h=h+h
      ihlf=-1
      istep=0
      iend=0
    4 if((x+h-xend)*h)7,6,5
    5 h=xend-x
    6 iend=1
    7 call outp(x,irec,ndim,prmt,ttau0,itemp,1,ibeam)
      if(prmt(5))40,8,40
    8 itest=0
    9 istep=istep+1
      j=1
   10 aj=a(j)
      bj=b(j)
      cj=c(j)
      do 11 i=1,ndim
      r1=h*dery(i)
      r2=aj*(r1-bj*aux(6,i))
      yarray(i)=yarray(i)+r2
      r2=r2+r2+r2
      aux(6,i)=aux(6,i)+r2-cj*r1
   11 continue
      if(j-4)12,15,15
   12 j=j+1
      if(j-3)13,14,13
   13 x=x+.5d0*h
   14 call fct(x,ibeam)
      goto 10
   15 if(itest)16,16,20
   16 do 17 i=1,ndim
       aux(4,i)=yarray(i)
   17 continue
      itest=1
      istep=istep+istep-2
   18 ihlf=ihlf+1
      x=x-h
      h=.5d0*h
      do 19 i=1,ndim
      yarray(i)=aux(1,i)
      dery(i)=aux(2,i)
      aux(6,i)=aux(3,i)
   19 continue
      goto 9
   20 imod=istep/2
      if(istep-imod-imod)21,23,21
   21 call fct(x,ibeam)
      do 22 i=1,ndim
      aux(5,i)=yarray(i)
      aux(7,i)=dery(i)
   22 continue
      goto 9
   23 delt=0.d0
      do 24 i=1,ndim
      delt=delt+aux(8,i)*dabs(aux(4,i)-yarray(i))
   24 continue
      if(delt-prmt(4))28,28,25
   25 if(ihlf-10)26,36,36
   26 do 27 i=1,ndim
      aux(4,i)=aux(5,i)
   27 continue
      istep=istep+istep-4
      x=x-h
      iend=0
      goto 18
   28 call fct(x,ibeam)
      do 29 i=1,ndim
      aux(1,i)=yarray(i)
      aux(2,i)=dery(i)
      aux(3,i)=aux(6,i)
      yarray(i)=aux(5,i)
      dery(i)=aux(7,i)
   29 continue
      call outp(x-h,ihlf,ndim,prmt,ttau0,itemp,2, ibeam)
      if(prmt(5))40,30,40
   30 do 31 i=1,ndim
      yarray(i)=aux(1,i)
      dery(i)=aux(2,i)
   31 continue
      irec=ihlf
      if(iend)32,32,39
   32 ihlf=ihlf-1
      istep=istep/2
      h=h+h
      if(ihlf)4,33,33
   33 imod=istep/2
      if(istep-imod-imod)4,34,4
   34 if(delt-.02d0*prmt(4))35,35,4
   35 ihlf=ihlf-1
      istep=istep/2
      h=h+h
      goto 4
   36 ihlf=11
      call fct(x,ibeam)
      goto 39
   37 ihlf=12
      goto 39
   38 ihlf=13
   39 call outp(x,ihlf,ndim,prmt,ttau0,itemp,3, ibeam)
   40 return
end subroutine


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

end module com_prog
