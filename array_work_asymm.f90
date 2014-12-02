module array_work
    save
!    sk- число регулярных участков в системе
!    nkr- число вариаций по радиусу
!    nka- число вариаций по углу
!    nbeam- количество пучков
!    mkk- число частиц в системе
!    mk0 - число частиц, которое влетает в систему на каждом временном шаге
!    ns2- 2*nk*sk - необходимо для определения длины массива y
!    ktimemax- число шагов по времени
!    ij1- через сколько шагов будут выдаваться поля и другие величины
!    tau - текущее время ( умноженное на w0)
!    tnorm - текущее время, нормированное на период
!    dt =2pi - нормированный шаг по времени
!    w0-  рабочая частота, необходимая для нормировок
!    rel_factor - релятивистский фактор
!    constq- коэффициент,пропорциональный заряду частицы, используется в fct
!    ddz - начальное расстояние между влетающими частицами
!    v0- скорость потока
!    zss- длина всей системы, необходима для расчета числа частиц в системе
!    beam_curr -ток потока
!    beam_voltage - ускоряющее напряжение
!    rb- радиус потока
!    gam0 = w0/c
!    wh0,wh,dw - частоты, используемые при определении спектра
    integer sk,ss00,sk1,sk2,sk3,sdop01,sdop02,sdop03, &
            sdop11,sdop12,sdop21,sdop22,sdop31,sdop32,mk0,ns2,ns4, &
            mkns,ktimemax, kluch_beam, kluch_structur(3),nka,nkr, &
            nbeam, mkk
    real*8  dt,rb0, ellips, alpha_shift
    real*8  w0,rel_factor,constq, const1, ddz,v0, &
            rt01,rt02,rt03,dz01,dz02,dz03, &
            rt1,rt2,dz11,dz12,dz21,dz22,dz31,dz32,zss, &
            beam_curr, beam_voltage, ampl_mod,gam0,sum, &
            wh0,wh,dw, pi, &
            step_prmt, prmt_4, periodz1, periodz2

    complex*16  ce,ee  !мнимая единица и реальная единица

!    nopen(:)- число открытых мод на каждом регулярном участке
!    y1(:),y(:),dery(:) - масивы для запоминания положения частиц,
!      их скоростей, наведенных токов для каждой моды и их производных
!      (dery(:)), длина массивов 	mkns=mkk*2+ns2.
!    gam(:,:,:) - массив постоянных распространения
!    zn(:,:)  - массив волновых сопротивлений
!    eznbm(:,:,:,:),ernbm(:,:,:,:)-  компоненты поля в области пучка
!    zs(:) - координаты скачков радиуса волновода
!    dz(:)-  длины регулярных участков (участков гладкого волновода)
!    rt(:)-   радиусы регулярных участков волновода
!    mu(:,:),nu(:) -  корни функций бесселя

    real*8, allocatable:: zs(:),dz(:),rt(:),rb(:), mu(:,:),alpha(:), &
        yarray(:),yarray1(:),velocity1(:), &
        dery(:), aux(:,:), velocity(:,:), all_yarray(:,:), all_dery(:,:)

    complex*16,allocatable::gam(:,:,:),zn(:,:,:),eznbm(:,:,:,:), ernbm(:,:,:,:), ephinbm(:,:,:,:), &
        xnplus(:,:,:),xnminus(:,:,:), &
        dxnplus(:,:,:),dxnminus(:,:,:),xbplus(:,:,:), &
        xbminus(:,:,:),dxbplus(:,:,:),dxbminus(:,:,:), &
        bplus0(:,:),bminussk(:,:), &
        etaplus(:,:,:),etaminus(:,:,:)

    real*8, allocatable:: amplxplus(:),amplxminus(:)

    integer, allocatable::ipiv(:), ipiv2(:)

    complex*16, allocatable::d1plus(:,:,:,:),d1minus(:,:,:,:), &
        d2plus(:,:,:,:),d2minus(:,:,:,:),dd1plus(:,:,:,:), &
        dd1minus(:,:,:,:), &
        dd2plus(:,:,:,:),dd2minus(:,:,:,:), &
        b1plus(:,:,:,:), b1minus(:,:,:,:),b2plus(:,:,:,:), &
        b2minus(:,:,:,:), &
        bb1plus(:,:,:,:), bb1minus(:,:,:,:),bb2plus(:,:,:,:), &
        bb2minus(:,:,:,:), &
        uste1(:,:,:,:), uste2(:,:,:,:),usth1(:,:,:,:), usth2(:,:,:,:), &
        fgammaplus(:,:,:,:), fgammaminus(:,:,:,:), &
        db1plus(:,:,:,:), db1minus(:,:,:,:), db2plus(:,:,:,:), db2minus(:,:,:,:), &
        ddb1plus(:,:,:,:), ddb1minus(:,:,:,:), &
        ddb2plus(:,:,:,:), ddb2minus(:,:,:,:), &
        ddb1plusinv(:,:,:,:), ddb1minusinv(:,:,:,:), &
        ddb2plusinv(:,:,:,:), ddb2minusinv(:,:,:,:), &
        dddb1(:,:,:,:), dddb2(:,:,:,:), &
        dddb1inv(:,:,:,:), dddb2inv(:,:,:,:), &
        alfaplus(:,:,:,:), alfaminus(:,:,:,:), &
        betaplus(:,:,:,:), betaminus(:,:,:,:), &
        hie1(:,:,:,:), hie2(:,:,:,:),hih1(:,:,:,:), hih2(:,:,:,:), &
        hie2inv(:,:,:,:),hih2inv(:,:,:,:), &
        psie1(:,:,:,:), psie2(:,:,:,:), psih1(:,:,:,:), psih2(:,:,:,:), &
        psie1inv(:,:,:,:), psie2inv(:,:,:,:), psih1inv(:,:,:,:), &
        psih2inv(:,:,:,:), &
        psipsi1(:,:,:,:), psipsi2(:,:,:,:), &
        psihi1(:,:,:,:), psihi2(:,:,:,:),psihi3(:,:,:,:),psihi4(:,:,:,:), &
        psipsi1inv(:,:,:,:), psipsi2inv(:,:,:,:), &
        uste1inv(:,:,:,:),uste2inv(:,:,:,:),usth1inv(:,:,:,:), &
        usth2inv(:,:,:,:), &
        aa1(:,:,:,:), aa2(:,:,:,:),aa3(:,:,:,:), ab(:,:,:), &
        rab1(:,:), rab2(:,:),xrab(:), &
        rnplus(:,:,:), rnminus(:,:,:), rnnplus(:,:,:), rnnminus(:,:,:), &
        brne1(:,:,:),brne2(:,:,:),brnh1(:,:,:),brnh2(:,:,:), &
        psibrn1(:,:,:),psibrn2(:,:,:), &
        hipsisk(:,:,:),temp(:,:),temp2(:,:), exit_sum(:,:), enter_sum(:,:)
!****************************************************
! для прогонки
!*****************************************************
    complex*16,allocatable:: alphap(:,:,:),betap(:,:),rabp(:,:), &
        rabpinv(:,:),xrabp(:)

end module array_work
