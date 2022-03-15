!#####################################################################
!!
!!  File  tables_rcl.f90
!!  is part of RECOLA (REcursive Computation of One Loop Amplitudes)
!!
!!  Copyright (C) 2015-2017   Stefano Actis, Ansgar Denner, 
!!                            Lars Hofer, Jean-Nicolas Lang, 
!!                            Andreas Scharf, Sandro Uccirati
!!
!!  RECOLA is licenced under the GNU GPL version 3, 
!!         see COPYING for details.
!!
!#####################################################################

  module tables_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use input_rcl
  use collier_interface_rcl
  use loop_functions_rcl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none
  complex(dp) :: deltar,sew0,xl2(3),xu2(3),xd2(3),m2
  real(dp)    :: dalZ

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine binaries_tables

  ! This subroutine creates many tables of the binary arithmetic.
  ! Given the number of external legs "lmax", the one-loop processes 
  ! generates binary numbers which are the sum of some of the primary
  ! binaries 1,2,4,...,2^(lmax+1).
  ! Given such a binary "e" the subroutine computes the following 
  ! tables.
  ! - "levelLeg(e)": 
  !   It says how many primary binaries are contained in "e".
  ! - "vectorLeg(e,:)": 
  !   It is a vector with "lmax+2" components. Each entry "i" of the 
  !   vector is 1 or 0, depending on whether the primary binary 
  !   "2^(i-1)" is contained in "e" or not.
  ! - firstNumber(e): 
  !   It is the lowest primary binary contained in "e". 
  !   Ex.: ! if e = 8 + 16 + 32 + 128, then firstNumber(e) = 8.
  ! - firstGap(e): 
  !   It is the lowest missing primary binary contained in "e". 
  !   Ex.: if e = 1 + 2 + 8 + 16 + 32 + 128, then firstGap(e)= 4.
  ! - firstNumbers(e,:lmax): 
  !   It is a vector whose entries are the ordered (from the lowest 
  !   to the biggest) primary binaries contained in "e". If the 
  !   primary binaries are less than "lmax" the corresponding entries 
  !   of the vector are 0. 
  !   Ex: if e = 1 + 2 + 8 + 16 + 32 + 128 and lmax = 10, then 
  !   firstNumbers(e,1)= 1, firstNumbers(e,2)= 2, 
  !   firstNumbers(e,3)= 8, firstNumbers(e,4)= 16, 
  !   firstNumbers(e,5)= 32, firstNumbers(e,6)= 128, 
  !   firstNumbers(e,7:10)= 0.
  ! - firstGaps(e,:lmax): 
  !   It is a vector whose entries are the ordered (from the lowest 
  !   to the biggest) missing primary binaries contained in "e".
  !   Ex: if e = 1 + 2 + 8 + 16 + 64 and lmax = 8, then 
  !   firstGaps(e,1)= 4, firstGaps(e,2)= 32, firstGaps(e,3)= 128, 
  !   firstGaps(e,4)= 256, firstGaps(e,5)= 512, firstGaps(e,6)= 1024, 
  !   firstGaps(e,7)= 2048.

  integer             :: lmax,e,ee,max,n,g,i,count
  integer,allocatable :: l(:),ns(:),gs(:)

  lmax = legsMax + 2

  max = 2**lmax - 1

  allocate (     levelLeg(0:max))
  allocate (    vectorLeg(0:max,lmax))
  allocate (firstNumber  (0:max))
  allocate (firstGap     (0:max))
  allocate (firstGaps    (0:max,lmax))
  allocate (firstNumbers (0:max,lmax))

  allocate ( l(lmax))
  allocate (ns(lmax))
  allocate (gs(lmax))

  do e = 0,max

    ee = e
    do i = lmax,1,-1
      l(i) = ee/2**(i-1)
      ee  = ee - l(i)*2**(i-1)
    enddo
    levelLeg (e)   = sum(l(:))
    vectorLeg(e,:) =     l(:)

    ee = e
    if (ee.eq.0) then
      firstNumber(e) = 0
    else
      n  = 1
      do while ( mod(ee,2) .eq. 0 ) ! as long as ee is even
        n  = 2*n
        ee = ee/2
      enddo
      firstNumber(e) = n
    endif

    ee = e
    g  = 1
    do while ( mod(ee,2) .eq. 1 ) ! as long as ee is odd
      g  = 2*g
      ee = (ee-1)/2
    enddo
    firstGap(e) = g

    ee    = e
    i     = 1
    count = 1
    do while ( i < lmax+1 )
      if (ee.eq.0) then
        ns(i) = 0
        i    = i + 1
      else
        if (mod(ee,2).eq.0) then ! ee is even
          ee   = ee/2
        else                     ! ee is odd
          ns(i) = count
          ee   = (ee-1)/2
          i    = i + 1
        endif
      endif
      count  = 2*count
    enddo
    firstNumbers(e,:) = ns

    ee    = e
    i     = 1
    count = 1
    do while ( i < lmax+1 )
      if (mod(ee,2).eq.0) then ! ee is even
        gs(i) = count
        ee    = ee/2
        i     = i + 1
      else                     ! ee is odd
        ee    = (ee-1)/2
      endif
      count  = 2*count
    enddo
    firstGaps(e,:) = gs

  enddo

  deallocate ( l)
  deallocate (gs)
  deallocate (ns)

  end subroutine binaries_tables

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine tensors_tables

  ! This routine creates tables for the integer arithmetic of tensors 
  ! of rank "r" from 0 to "Rmax = legsMax".
  ! - gg(0:3,0:3):
  !   It is the metric tensor.
  ! - RtoS(r):
  !   It is the size of the symmetrized tensor index at rank "r":
  !     RtoS(0)= 1, RtoS(1)= 4, RtoS(2)= 10, ...
  ! - ri(r,i), riMin(r), riMax(r), RItoR(ri), RItoI(ri):
  !   The symmetrized index "i" of a tensor of rank "r" is combined 
  !   with the rank "r" to produce a single index "ri". Given the 
  !   rank "r" and the index "i":
  !     ri(r,i) = the value of "ri"
  !     riMin(r) = the minimal value of "ri" for that "r"
  !     riMax(r) = the maximal value of "ri" for that "r"
  !   Given the combined index "ri":
  !     RItoR(ri) = rank "r" for that "ri"
  !     RItoI(ri) = index "i" for that "ri"
  ! - incRI(0:3,0:riMax(Rmax-1)):
  !   The combined symmetrized index "ri" of a tensor is combined 
  !   with a lorentz index 0:3 to create a new combined symmetrized 
  !   index.

  integer              :: Rmax,repeat,cntRI,cntI,r,n0,n1,n2,n3, &
                          riIn,nIn(0:3),nOut(0:3),mu,i
  integer, allocatable :: RItoN0(:),RItoN1(:),RItoN2(:),RItoN3(:), &
                          N0123toRI(:,:,:,:)
  logical, allocatable :: computed(:)

  ! metric tensor
  gg      =   0
  gg(0,0) = + 1
  gg(1,1) = - 1
  gg(2,2) = - 1
  gg(3,3) = - 1

  Rmax = legsMax

  allocate (RtoS      (0:Rmax))
  allocate (riMin     (0:Rmax))
  allocate (riMax     (0:Rmax))

  ! The combined symmetrized index "ri" of a tensor is transformed 
  ! in the 4-index notation (n0,n1,n2,n3), where n0 is the number 
  ! of lorentz indices set to 0, n1 is the number of lorentz 
  ! indices set to 1, etc.
  ! Given the 4 indices "n0,n1,n2,n3":
  !   N0123toRI(n0,n1,n2,n3) = the value of the index "ri"
  ! Given the combined index "ri":
  !   RItoN0(ri) = the value of the index n0
  !   RItoN1(ri) = the value of the index n1
  !   RItoN2(ri) = the value of the index n2
  !   RItoN3(ri) = the value of the index n3

  allocate (N0123toRI (0:Rmax,0:Rmax,0:Rmax,0:Rmax))

  do repeat = 1,2

    cntRI = 0

    do r = 0,Rmax

      cntI = 0

      riMin(r) = cntRI

      do n0 = r,0,-1
        do n1 = r-n0,0,-1
          do n2 = r-n0-n1,0,-1
            n3 = r-n0-n1-n2

            cntI  =  cntI + 1

            if (repeat.eq.2) then
              N0123toRI(n0,n1,n2,n3) = cntRI
              RItoN0(cntRI) = n0
              RItoN1(cntRI) = n1
              RItoN2(cntRI) = n2
              RItoN3(cntRI) = n3
              ri(cntI,r) = cntRI
              RItoR(cntRI) = r
              RItoI(cntRI) = cntI
            endif

            cntRI = cntRI + 1

          enddo
        enddo
      enddo

      RtoS(r) = cntI

      riMax(r) = cntRI - 1

    enddo

    if (repeat.eq.1) then
      riTot = cntRI
      allocate (ri    (RtoS(Rmax),0:Rmax))
      allocate (RItoR (0:riTot))
      allocate (RItoI (0:riTot))
      allocate (RItoN0 (0:riTot))
      allocate (RItoN1 (0:riTot))
      allocate (RItoN2 (0:riTot))
      allocate (RItoN3 (0:riTot))
    endif

  enddo

  allocate (incRI (0:3,0:riMax(Rmax-1)))
  allocate (firstRI (0:3,0:riMax(Rmax-1)))
  allocate (computed(riTot))

  firstRI = .false.
  computed = .false.
  do riIn = riMax(Rmax-1),0,-1
    nIn(0) = RItoN0(riIn)
    nIn(1) = RItoN1(riIn)
    nIn(2) = RItoN2(riIn)
    nIn(3) = RItoN3(riIn)
    do mu = 0,3
      nOut = nIn
      nOut(mu) = nIn(mu) + 1
      i = N0123toRI(nOut(0),nOut(1),nOut(2),nOut(3))
      incRI(mu,riIn) = i
      if(.not.computed(i)) firstRI(mu,riIn) = .true.
      computed(i) = .true.
    enddo
  enddo

  deallocate (computed)

  deallocate (RItoN3)
  deallocate (RItoN2)
  deallocate (RItoN1)
  deallocate (RItoN0)
  deallocate (N0123toRI)

  end subroutine tensors_tables

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine particles_tables

  ! This subroutine compute tables for fields and particles

  integer :: i

  cpar( 1) = 'X_g'
  cpar( 2) = 'X_A'
  cpar( 3) = 'X_Z'
  cpar( 4) = 'X+'
  cpar( 5) = 'X-'
  cpar( 6) = 'X_g~'
  cpar( 7) = 'X_A~'
  cpar( 8) = 'X_Z~'
  cpar( 9) = 'X+~'
  cpar(10) = 'X-~'
  cpar(11) = 'H'
  cpar(12) = 'p0'
  cpar(13) = 'p+'
  cpar(14) = 'p-'
  cpar(15) = 'g'
  cpar(16) = 'A'
  cpar(17) = 'Z'
  cpar(18) = 'W+'
  cpar(19) = 'W-'
  cpar(20) = 'nu_e'
  cpar(21) = 'nu_mu'
  cpar(22) = 'nu_tau'
  cpar(23) = 'u'
  cpar(24) = 'c'
  cpar(25) = 't'
  cpar(26) = 'e-'
  cpar(27) = 'mu-'
  cpar(28) = 'tau-'
  cpar(29) = 'd'
  cpar(30) = 's'
  cpar(31) = 'b'
  cpar(32) = 'nu_e~'
  cpar(33) = 'nu_mu~'
  cpar(34) = 'nu_tau~'
  cpar(35) = 'u~'
  cpar(36) = 'c~'
  cpar(37) = 't~'
  cpar(38) = 'e+'
  cpar(39) = 'mu+'
  cpar(40) = 'tau+'
  cpar(41) = 'd~'
  cpar(42) = 's~'
  cpar(43) = 'b~'

  do i = 1,nFs

    select case (i)
    case (11);    parKind (i) = 1
    case (15,16); parKind (i) = 2
    case (17:19); parKind (i) = 3
    case (20:31); parKind (i) = 4
    case (32:43); parKind (i) = 5
    case default; parKind (i) = 0
    end select

    select case (i)
    case ( 1: 5); cftype(i) = 'x'    ! ghost field
    case ( 6:10); cftype(i) = 'x~'   ! anti-ghost field
    case (11:14); cftype(i) = 's'    ! scalar
    case (15:19); cftype(i) = 'v'    ! vector
    case (20:31); cftype(i) = 'f'    ! fermion
    case (32:43); cftype(i) = 'f~'   ! anti-fermion
    case default; cftype(i) = '0'
    end select

    select case (i)
    case ( 1);          cftype2(i) = 'G'    ! coloured ghost field
    case ( 2: 5);       cftype2(i) = 'x'    ! uncoloured ghost field
    case ( 6);          cftype2(i) = 'G~'   ! coloured anti-ghost field
    case ( 7:10);       cftype2(i) = 'x~'   ! uncoloured anti-ghost field
    case (11:14);       cftype2(i) = 's'    ! scalar
    case (15);          cftype2(i) = 'g'    ! gluon
    case (16:19);       cftype2(i) = 'v'    ! vector
    case (20:22,26:28); cftype2(i) = 'l'    ! lepton
    case (23:25,29:31); cftype2(i) = 'q'    ! quark
    case (32:34,38:40); cftype2(i) = 'l~'   ! anti-lepton
    case (35:37,41:43); cftype2(i) = 'q~'   ! anti-quark
    case default;       cftype2(i) = '0'
    end select

    select case (i)
    case(1:3,6:8,11:12,15:17,20:22,32:34); charged(i) = .false.
    case default;                          charged(i) = .true.
    end select

    select case (i)
    case (1:3,6:8,11:12,15:17,20:22,32:34); threeQ(i) =   0
    case (4,10,13,18,38:40);                threeQ(i) = + 3
    case (5, 9,14,19,26:28);                threeQ(i) = - 3
    case (23:25);                           threeQ(i) = + 2
    case (35:37);                           threeQ(i) = - 2
    case (41:43);                           threeQ(i) = + 1
    case (29:31);                           threeQ(i) = - 1
    end select

    select case (i)
    case (1:5);   Qghost(i) = + 1
    case (6:10);  Qghost(i) = - 1
    case default; Qghost(i) =   0
    end select

  enddo

  ! Number of quarks
  Nq = 0
  do i = 1,nFs
    if (cftype2(i).eq.'q') Nq = Nq + 1
  enddo
  Nf = Nq*1d0

  end subroutine particles_tables

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine masses_tables

  ! This table computes tables for masses

  integer     :: nmz,nmw,nmh,nmu,nmc,nmt,nmel,nmmu,nmta,nmd,nms,nmb, &
                 defm,i,j,pr
  real(dp)    :: zero_mass,m2(6)
  complex(dp) :: cmz2,cmw2,cmh2,cmu2,cmc2,cmt2,cmel2,cmmu2,cmta2, &
                 cmd2,cms2,cmb2

  zero_mass = zerocut *                                         &
              max(mass_z,mass_w,mass_h,mass_el,mass_mu,mass_ta, &
                  mass_u,mass_c,mass_t,mass_d,mass_s,mass_b     )

  nmz  = 2
  cmz2 = mass_z**2 - cId0*width_z*mass_z

  nmw  = 3
  cmw2 = mass_w**2 - cId0*width_w*mass_w

  nmh  = 4
  cmh2 = mass_h**2 - cId0*width_h*mass_h

  defm = 4

  if (mass_u.lt.zero_mass) then
    nmu  = 1
    cmu2 = c0d0
  else
    nmu  = defm + 1; defm = nmu
    cmu2 = cmplx(mass_u**2,kind=dp)
  endif
  if (mass_c.lt.zero_mass) then
    nmc  = 1
    cmc2 = c0d0
  else
    nmc  = defm + 1; defm = nmc
    cmc2 = mass_c**2 - cId0*width_c*mass_c
  endif
  if (mass_t.lt.zero_mass) then
    nmt  = 1
    cmt2 = c0d0
  else
    nmt  = defm + 1; defm = nmt
    cmt2 = mass_t**2 - cId0*width_t*mass_t
  endif
  if (mass_el.lt.zero_mass) then
    nmel = 1
    cmel2 = c0d0
  else
    nmel = defm + 1; defm = nmel
    cmel2 = cmplx(mass_el**2,kind=dp)
  endif
  if (mass_mu.lt.zero_mass) then
    nmmu = 1
    cmmu2 = c0d0
  else
    nmmu = defm + 1; defm = nmmu
    cmmu2 = mass_mu**2 - cId0*width_mu*mass_mu;
  endif
  if (mass_ta.lt.zero_mass) then
    nmta = 1
    cmta2 = c0d0
  else
    nmta = defm + 1; defm = nmta
    cmta2 = mass_ta**2 - cId0*width_ta*mass_ta
  endif
  if (mass_d.lt.zero_mass) then
    nmd  = 1
    cmd2 = c0d0
  else
    nmd  = defm + 1; defm = nmd
    cmd2 = cmplx(mass_d**2,kind=dp)
  endif
  if (mass_s.lt.zero_mass) then
    nms  = 1
    cms2 = c0d0
  else
    nms  = defm + 1; defm = nms
    cms2 = cmplx(mass_s**2,kind=dp)
  endif
  if (mass_b.lt.zero_mass) then
    nmb  = 1
    cmb2 = c0d0
  else
    nmb  = defm + 1; defm = nmb
    cmb2 = mass_b**2 - cId0*width_b*mass_b
  endif

  nmasses = defm

  allocate(cm2n(nmasses)); cm2n(1)  = c0d0

  ! massless: cm2n(1) = 0
  ! light:    cm2n(n) not 0 (n>1)
  ! massive:  cm2n(n) not 0 (n>1)
  cm2n(nmz)  = cmz2
  cm2n(nmw)  = cmw2
  cm2n(nmh)  = cmh2
  cm2n(nmel) = cmel2
  cm2n(nmmu) = cmmu2
  cm2n(nmta) = cmta2
  cm2n(nmu)  = cmu2
  cm2n(nmd)  = cmd2
  cm2n(nmc)  = cmc2
  cm2n(nms)  = cms2
  cm2n(nmt)  = cmt2
  cm2n(nmb)  = cmb2

  do i = 1,nFs
    ! massless: nmf(i) = 1; regf(i) = 1; cm2f(i) = 0
    ! light:    nmf(i) > 1; regf(i) = 2; cm2f(i) = 0
    ! massive:  nmf(i) > 1; regf(i) = 3; cm2f(i) not 0
    select case (i)
    case (3,8,12,17)
      nmf(i) = nmz;       regf(i) = 3; cm2f(i) = cmz2
    case (4,5,9,10,13,14,18,19)
      nmf(i) = nmw;       regf(i) = 3; cm2f(i) = cmw2
    case (11)
      nmf(i) = nmh;       regf(i) = 3; cm2f(i) = cmh2
    case (23,35)
      nmf(i) = nmu
      if (light_u)  then; regf(i) = 2; cm2f(i) = c0d0
      else;               regf(i) = 3; cm2f(i) = cmu2 
      endif
    case (24,36)
      nmf(i) = nmc
      if (light_c)  then; regf(i) = 2; cm2f(i) = c0d0
      else;               regf(i) = 3; cm2f(i) = cmc2 
      endif
    case (25,37)
      nmf(i) = nmt
      if (light_t)  then; regf(i) = 2; cm2f(i) = c0d0
      else;               regf(i) = 3; cm2f(i) = cmt2 
      endif
    case (26,38)
      nmf(i) = nmel
      if (light_el) then; regf(i) = 2; cm2f(i) = c0d0
      else;               regf(i) = 3; cm2f(i) = cmel2
      endif
    case (27,39)
      nmf(i) = nmmu
      if (light_mu) then; regf(i) = 2; cm2f(i) = c0d0
      else;               regf(i) = 3; cm2f(i) = cmmu2
      endif
    case (28,40)
      nmf(i) = nmta
      if (light_ta) then; regf(i) = 2; cm2f(i) = c0d0
      else;               regf(i) = 3; cm2f(i) = cmta2
      endif
    case (29,41)
      nmf(i) = nmd
      if (light_d)  then; regf(i) = 2; cm2f(i) = c0d0
      else;               regf(i) = 3; cm2f(i) = cmd2 
      endif
    case (30,42)
      nmf(i) = nms
      if (light_s)  then; regf(i) = 2; cm2f(i) = c0d0
      else;               regf(i) = 3; cm2f(i) = cms2 
      endif
    case (31,43)
      nmf(i) = nmb
      if (light_b)  then; regf(i) = 2; cm2f(i) = c0d0
      else;               regf(i) = 3; cm2f(i) = cmb2 
      endif
    case default
      nmf(i) = 1;         regf(i) = 1; cm2f(i) = c0d0 
    end select
    if (nmf(i).eq.1) regf(i) = 1
!    call openOutput
!    write(nx,*) 'Particle:',i,cpar(i)
!    write(nx,*) 'Masses in loop functions:',nmf(i),cm2n(nmf(i))
!    write(nx,*) 'Masses in coefficients:  ',regf(i),cm2f(i)
!    write(nx,*) 
  enddo

  cm2pf = cm2f
  do i = 1,nFs
    if (resPar(i)) then
      if (abs(real(aimag(cm2pf(i)),kind=dp)).lt.zerocut) then
        if (warnings.le.warning_limit) then
          warnings = warnings + 1
          call openOutput
          write(nx,*)
          write(nx,*) 'ERROR: set_resonant_particle_rcl called '
          write(nx,*) '       for a particle with vanishing decay width'
          write(nx,*)
          call toomanywarnings
        endif
        call istop (ifail,1)
      endif
      cm2f(i) = real(cm2pf(i),kind=dp)*c1d0
      cm2n(nmf(i)) = real(cm2n(nmf(i)),kind=dp)*c1d0
    endif
!    call openOutput
!    write(nx,*) 'Particle:',i,cpar(i)
!    write(nx,*) 'Masses in loop functions:',nmf(i),cm2n(nmf(i))
!    write(nx,*) 'Masses in coefficients:  ',regf(i),cm2f(i)
!    write(nx,*) 'Masses in res propagators',regf(i),cm2pf(i)
!    write(nx,*) 
  enddo

  if (reg_soft.eq.2) then ! mass regularization for soft singularities
    do pr = 1,prTot
      do i = 1,legsIn(pr)+legsOut(pr)
        j = par(i,pr)
        if ((cftype(j).eq.'f'.or.cftype(j).eq.'f~').and. &
            charged(j).and.(nmf(j).eq.1)) then
          if (warnings.le.warning_limit) then
            warnings = warnings + 1
            call openOutput
            write(nx,*)
            write(nx,*) 'ERROR: Mass regularization for soft ', &
                               'singularities not allowed '
            write(nx,*) '       if massless external charged fermions ', &
                               'are present'
            write(nx,*)
            call toomanywarnings
          endif
          call istop (ifail,1)
        endif
      enddo
    enddo
  endif

  ! set mq2(:), a vector with the ordered quark masses:
  ! mq2(1) .le. mq2(2) .le. mq2(3) .le. mq2(4) .le. mq2(5) .le. mq2(6)
  if (complex_mass_scheme.eq.1) then
    m2(1) = abs(cmd2)
    m2(2) = abs(cmu2)
    m2(3) = abs(cms2)
    m2(4) = abs(cmc2)
    m2(5) = abs(cmb2)
    m2(6) = abs(cmt2)
  else
    m2(1) = real(cmd2,kind=dp)
    m2(2) = real(cmu2,kind=dp)
    m2(3) = real(cms2,kind=dp)
    m2(4) = real(cmc2,kind=dp)
    m2(5) = real(cmb2,kind=dp)
    m2(6) = real(cmt2,kind=dp)
  endif
  mq2 = m2
  do j = 1,6
    do i = j+1,6
      if (m2(i).lt.mq2(j)) then
        mq2(i) = m2(j)
        mq2(j) = m2(i)
      endif
      m2 = mq2
    enddo
  enddo
  if (use_active_qmasses) mq2 = mq**2
!  call openOutput
!  write(nx,*)
!  do i = 1,6
!    write(nx,*) i,sqrt(mq2(i))
!  enddo
!  write(nx,*)

  ! set Nlq, the number of active flavours for alpha_s renormalization. 
  select case (Nfren)
  case (-1)
    ! Nlq is the number of quarks lighter than Qren
    Nlq = 0
    do i = 1,6
      if (mq2(i).lt.Qren**2) Nlq = Nlq + 1
    enddo
  case (3,4,5)
    if (mq2(Nfren+1).eq.0d0) then
      if (warnings.le.warning_limit) then
        warnings = warnings + 1
        call openOutput
        write(nx,*)
        write(nx,*) 'ERROR: Wrong number Nfren of light flavours ', &
                           'for alphas renormaliztion'
        write(nx,*) '       (Nfren can not be smaller than the ', &
                           'number of massless quarks)'
        write(nx,*)
        call toomanywarnings
      endif
      call istop (ifail,1)
    endif
    Nlq = Nfren
  case (6)
    Nlq = Nfren
  case default
    if (warnings.le.warning_limit) then
      warnings = warnings + 1
      call openOutput
      write(nx,*)
      write(nx,*) 'ERROR: Wrong value for the number Nfren of light ', &
                         'flavours for alphas renormaliztion'
      write(nx,*) '       (accepted values are Nfren = -1,3,4,5,6)'
      write(nx,*)
      call toomanywarnings
    endif
    call istop (ifail,1)
  end select

  ! Set some useful variables for masses
  if (complex_mass_scheme.eq.1) then
    mh2 = cm2f(11)
    mz2 = cm2f(17) 
    mw2 = cm2f(18) 
  else
    mh2 = real(cm2f(11),kind=dp)*c1d0
    mz2 = real(cm2f(17),kind=dp)*c1d0
    mw2 = real(cm2f(18),kind=dp)*c1d0
  endif
  mh4 = mh2*mh2 
  mz4 = mz2*mz2; mz6 = mz4*mz2
  mw4 = mw2*mw2; mw6 = mw4*mw2; mw1 = sqrt(mw2); mw3 = mw2*mw1

  ct2 = mw2/mz2;   ct = sqrt(ct2)
  st2 = 1d0 - ct2; st = sqrt(st2)
  st3 = st2*st; st4 = st2*st2; st6 = st4*st2; st8 = st6*st2
  ct3 = ct2*ct; ct4 = ct2*ct2; ct6 = ct2*ct4
  stct = st*ct

  do i = 1,3
    if (complex_mass_scheme.eq.1) then
      mu2(i) = cm2f(i+22)
      ml2(i) = cm2f(i+25)
      md2(i) = cm2f(i+28)
    else
      mu2(i) = real(cm2f(i+22),kind=dp)*c1d0
      ml2(i) = real(cm2f(i+25),kind=dp)*c1d0
      md2(i) = real(cm2f(i+28),kind=dp)*c1d0
    endif
    mu1(i) = sqrt(mu2(i)); mu4(i) = mu2(i)*mu2(i)
    ml1(i) = sqrt(ml2(i)); ml4(i) = ml2(i)*ml2(i)
    md1(i) = sqrt(md2(i)); md4(i) = md2(i)*md2(i)
  enddo

  end subroutine masses_tables

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine couplings_tables

  integer  :: i
  real(dp) :: rew2,rez2

  ! pi
  pi2  = pi*pi

  ! regularization scheme: lam=1 HV, lam=0 FDH
  ! set it to zero (Pittau 1111.4965)
  select case (reguScheme)
  case (1); lam = 0d0
  case (2); lam = 1d0
  end select

  CMscheme = complex_mass_scheme

  ! EW couplings
  rez2 = real(mz2,kind=dp)
  rew2 = real(mw2,kind=dp)
  select case (ew_reno_scheme)
  case (1)                                             ! alpha_gf
    if (algf.lt.0d0) then
      alpha = gf*sq2/pi*rew2*( 1d0 - rew2/rez2 )
    else
      alpha = algf
      gf = alpha/sq2*pi/rew2/( 1d0 - rew2/rez2 )
    endif
  case (2); alpha = al0                                ! alpha_0
  case (3); alpha = alZ                                ! alpha_Z
  end select

  ! charge
  Ql = - 1d0;   Ql2 = Ql*Ql; Ql3 = Ql2*Ql; Ql4 = Ql2*Ql2
  Qu = + 2/3d0; Qu2 = Qu*Qu; Qu3 = Qu2*Qu; Qu4 = Qu2*Qu2
  Qd = - 1/3d0; Qd2 = Qd*Qd; Qd3 = Qd2*Qd; Qd4 = Qd2*Qd2

  ! isospin
  I3n =   1/2d0; I3n2 = I3n*I3n; I3n3 = I3n2*I3n; I3n4 = I3n2*I3n2
  I3l = - 1/2d0; I3l2 = I3l*I3l; I3l3 = I3l2*I3l; I3l4 = I3l2*I3l2
  I3u =   1/2d0; I3u2 = I3u*I3u; I3u3 = I3u2*I3u; I3u4 = I3u2*I3u2
  I3d = - 1/2d0; I3d2 = I3d*I3d; I3d3 = I3d2*I3d; I3d4 = I3d2*I3d2

  ! g+ and g-
  gpn = c0d0        ; gmn = ( I3n          )/stct
  gpl = - Ql * st/ct; gml = ( I3l - st2*Ql )/stct
  gpu = - Qu * st/ct; gmu = ( I3u - st2*Qu )/stct
  gpd = - Qd * st/ct; gmd = ( I3d - st2*Qd )/stct
  gpn2 = gpn**2; gmn2 = gmn**2
  gpl2 = gpl**2; gml2 = gml**2
  gpu2 = gpu**2; gmu2 = gmu**2
  gpd2 = gpd**2; gmd2 = gmd**2

  ! QCD coupling
  als0   = als
  Qren0  = Qren
  Nfren0 = Nfren
  Nlq0   = Nlq
  allocate(als0R(0:1,prTot));  als0R  = als0
  allocate(Qren0R(0:1,prTot)); Qren0R = Qren0
  allocate(Nlq0R(0:1,prTot));  Nlq0R  = Nlq0

  ! number of colours
  Nc  = real(nCs,kind=dp)
  Nc2 = Nc*Nc
  Cf = (Nc2-1d0)/(2d0*Nc)
  Ca = Nc

  ! QCD beta function at 1- and 2-loop level
  do i = 1,6
    beta0(i) =  1/4d0 * ( 11d0 - 2/3d0*i )
    beta1(i) = 1/16d0 * ( 102d0 - 38/3d0*i )
    rb1(i)   = beta1(i)/beta0(i)
  enddo

  end subroutine couplings_tables

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine counterterms_tables

  integer     :: i,nml(3),nmu(3),nmd(3),regl(3),regu(3),regd(3)
  real(dp)    :: el,el2,el3,el4,alphas,gs,gs2,gs3,gs4,         &
                 Duv,Dir,muUV2,muIR2,ratioIm,             &
                 rmh2,rmz2,rmw2,rmz4,rmw4,                     &
                 rmu(3),rml(3),rmd(3),rmu2(3),rml2(3),rmd2(3), &
                 rxu(3),rxl(3),rxd(3),rxu2(3),rxl2(3),rxd2(3)
  complex(dp) :: logIr,               &
                 sehh,seaz0,seazz,sezz,seww,                    &
                 seL,seR,seS,m2dseL,m2dseR,m2dseS,                   &
                 seLqed,seRqed,seSqed,m2dseLqed,m2dseRqed,m2dseSqed, &
                 seLqcd,seRqcd,seSqcd,m2dseLqcd,m2dseRqcd,m2dseSqcd, &
                 AW,AZ,AH,Al(3),Au(3),Ad(3),                         &
! p^2 = Mh^2
                 B0HWW,B0HZZ,B0HHH,DB0HHH,DB0HWW,DB0HZZ, &
                 B0Hll(3),B0Huu(3),B0Hdd(3),       &
                 DB0Hll(3),DB0Huu(3),DB0Hdd(3), &
! p^2 = 0
                 B00WW,B00ll(3),B00uu(3),B00dd(3), &
                 DB00ll(3),DB00uu(3),DB00dd(3), &
! p^2 = Mz^2
                 B0ZWW,B0ZHZ,B0Z00,DB0ZWW,DB0ZHZ,DB0Z00, &
                 B0Zll(3),B0Zuu(3),B0Zdd(3),  &
                 DB0Zuu(3),DB0Zll(3),DB0Zdd(3), &
! p^2 = Mw^2
                 B0W0W,B0WWZ,B0WHW,DB0W0W,DB0WWZ,DB0WHW, &
                 B0W0l(3),B0Wud(3),DB0W0l(3),DB0Wud(3),             &
! p^2 = 0
                 B000W,B00WZ,B00HW,DB000W,DB00WZ,DB00HW, &
                 B000l(3),B00ud(3),DB000l(3),DB00ud(3), &
! p^2 = Mn^2 = 0
                 B100Z(3),B10lW(3), &
! p^2 = Ml^2
                 B0llZ(3),B0llH(3),B0l0W(3), &
                 B1llZ(3),B1l0W(3),B1llH(3), &
                 DB0llZ(3),DB0llH(3),DB0l0W(3), &
                 DB1llZ(3),DB1l0W(3),DB1llH(3), &
                 B0ll0(3),B1ll0(3),DB0ll0(3),DB1ll0(3), &
! p^2 = Mu^2
                 B0uuZ(3),B0uuH(3),B0udW(3), &
                 B1uuZ(3),B1udW(3),B1uuH(3), &
                 DB0uuZ(3),DB0uuH(3),DB0udW(3), &
                 DB1uuZ(3),DB1udW(3),DB1uuH(3), &
                 B0uu0(3),B1uu0(3),DB0uu0(3),DB1uu0(3), &
! p^2 = Md^2
                 B0ddZ(3),B0ddH(3),B0duW(3), &
                 B1ddZ(3),B1duW(3),B1ddH(3), &
                 DB0ddZ(3),DB0ddH(3),DB0duW(3), &
                 DB1ddZ(3),DB1duW(3),DB1ddH(3), &
                 B0dd0(3),B1dd0(3),DB0dd0(3),DB1dd0(3)


  el = 2*sqrt(pi*alpha)
  el2  = el*el
  el3  = el2*el
  el4  = el2*el2

  alphas = als0
  gs = 2*sqrt(pi*alphas)
  gs2  = gs*gs
  gs3  = gs2*gs
  gs4  = gs2*gs2

  ! UV and IR single poles
  Duv = DeltaUV
  Dir = DeltaIR

  ! real squared masses for dimensional UV/IR regularization
  muUV2 = muUV**2
  muIR2 = muIR**2

  ! fermion masses for loop functions
  do i = 1,3
    nmu(i) = nmf(i+22)
    nml(i) = nmf(i+25)
    nmd(i) = nmf(i+28)
    regu(i) = regf(i+22)
    regl(i) = regf(i+25)
    regd(i) = regf(i+28)
    if (complex_mass_scheme.eq.1) then
      xu2(i) = cm2n(nmu(i))
      xl2(i) = cm2n(nml(i))
      xd2(i) = cm2n(nmd(i))
    else
      xu2(i) = real(cm2n(nmu(i)),kind=dp)*c1d0
      xl2(i) = real(cm2n(nml(i)),kind=dp)*c1d0
      xd2(i) = real(cm2n(nmd(i)),kind=dp)*c1d0
    endif
  enddo

  ! real squared masses
  rmh2 = real(mh2,kind=dp)
  rmz2 = real(mz2,kind=dp); rmz4 = rmz2*rmz2
  rmw2 = real(mw2,kind=dp); rmw4 = rmw2*rmw2
  rmu2 = real(mu2,kind=dp); rmu = sqrt(rmu2)
  rml2 = real(ml2,kind=dp); rml = sqrt(rml2)
  rmd2 = real(md2,kind=dp); rmd = sqrt(rmd2)
  rxu2 = real(xu2,kind=dp); rxu = sqrt(rxu2)
  rxl2 = real(xl2,kind=dp); rxl = sqrt(rxl2)
  rxd2 = real(xd2,kind=dp); rxd = sqrt(rxd2)

  ! A0 functions
  AH = A0 (mh2,Duv,muUV2)
  AZ = A0 (mz2,Duv,muUV2)
  AW = A0 (mw2,Duv,muUV2)
  do i = 1,3
    Al(i) = A0 (ml2(i),Duv,muUV2)
    Au(i) = A0 (mu2(i),Duv,muUV2)
    Ad(i) = A0 (md2(i),Duv,muUV2)
  enddo

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dt
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! the expression for the tadpole renormalization constant
! has been derived by hand; excluding overall normalization,
! it agrees with hep-ph/0612122 and Degrassi-Sirlin,
! NP B 383 (1992) 73
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Bosonic part
  dt = el*mw1/(8*pi2*st)*(              &
       - 3*mh2/(8*mw2)*AH               &
       - ( mh2/(8*mw2) + 3/(4*ct2) )*AZ &
       - ( mh2/(4*mw2) + 3/2d0 )*AW     &
       + mz2/(2*ct2) + mw2              &
       )

  ! Fermionic part
  do i = 1,3
    dt = dt &
         + el*mw1/(8*pi2*st)*(    &
           + ml2(i)/mw2*Al(i)     &
           + 3*mu2(i)/mw2*Au(i)   &
           + 3*md2(i)/mw2*Ad(i)   &
         )
  enddo

  if (reguScheme.eq.1) then
    dt = dt + el*mw3/(8*pi2*st)*( 1d0 + 1/(2*ct4) ) 
  endif

  if (.not.loopWEAK) dt = c0d0

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZh, dmh2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZh is computed through an expansion around the squared
! real H mass, and neglecting O(alpha^2) terms, see eq.(4.20)
! of hep-ph/0505042
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dmh2 is computed through an expansion around the squared
! real H mass, and neglecting O(alpha^3) terms, see eq.(4.20)
! of hep-ph/0505042.
! Expressions coded from 0709.1075 [hep-ph], eq.(B.5), for
! the Higgs boson self-energy
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Bosonic contributions to sehh and dZh

   B0HWW =  B0 (rmh2,mw2,mw2,Duv,muUV2)
  DB0HWW = DB0 (rmh2,mw2,mw2)

   B0HZZ =  B0 (rmh2,mz2,mz2,Duv,muUV2)
  DB0HZZ = DB0 (rmh2,mz2,mz2)

   B0HHH =  B0 (rmh2,mh2,mh2,Duv,muUV2)
  DB0HHH = DB0 (rmh2,mh2,mh2)

  sehh = - alpha/(8*pi*st2)*(                         &
           - ( 6*mw2 - 2*rmh2 + mh4/(2*mw2) )*B0HWW   &
           - ( 3d0 + mh2/(2*mw2) )*AW + 6*mw2         &
           - 1/(2*ct2)*(                              &
             + ( 6*mz2 - 2*rmh2 + mh4/(2*mz2) )*B0HZZ &
             + ( 3d0 + mh2/(2*mz2) )*AZ - 6*mz2       &
             )                                        &
           - 3/4d0*( 3*mh4/mw2*B0HHH + mh2/mw2*AH )   &
           )
  dZh = - alpha/(32*pi*st2)*(                                    &
          + 2/ct2/mz2*(                                          &
            - 4*mw2*B0HWW + ( mh4 - 4*rmh2*mw2 + 12*mw4 )*DB0HWW &
            )                                                    &
          + 1/mw2*(                                              &
            - 4*mz2*B0HZZ + ( mh4 - 4*rmh2*mz2 + 12*mz4 )*DB0HZZ &
            )                                                    &
          + 9*mh4/mw2*DB0HHH                                     &
          )

  ! Fermionic contributions to sehh and dZh

  do i = 1,3

     B0Hll(i) =  B0 (rmh2,ml2(i),ml2(i),Duv,muUV2)
    DB0Hll(i) = DB0 (rmh2,ml2(i),ml2(i))

     B0Huu(i) =  B0 (rmh2,mu2(i),mu2(i),Duv,muUV2)
    DB0Huu(i) = DB0 (rmh2,mu2(i),mu2(i))

     B0Hdd(i) =  B0 (rmh2,md2(i),md2(i),Duv,muUV2)
    DB0Hdd(i) = DB0 (rmh2,md2(i),md2(i))

    sehh = sehh                                                     &
           - alpha/(8*pi*st2*mw2)*(                                 &
             +    ml2(i)*( 2*Al(i) + ( 4*ml2(i) - rmh2 )*B0Hll(i) ) &
             + Nc*md2(i)*( 2*Ad(i) + ( 4*md2(i) - rmh2 )*B0Hdd(i) ) &
             + Nc*mu2(i)*( 2*Au(i) + ( 4*mu2(i) - rmh2 )*B0Huu(i) ) &
             )

    dZh = dZh                                                        &
          + alpha/(8*pi*mw2*st2)*(                                   &
            -    ml2(i)*( B0Hll(i) + ( rmh2 - 4*ml2(i) )*DB0Hll(i) ) &
            - Nc*mu2(i)*( B0Huu(i) + ( rmh2 - 4*mu2(i) )*DB0Huu(i) ) &
            - Nc*md2(i)*( B0Hdd(i) + ( rmh2 - 4*md2(i) )*DB0Hdd(i) ) &
            )

  enddo

  dmh2 = sehh - ( mh2 - rmh2 )*dZh

  if (reguScheme.eq.1) then
    dmh2 = dmh2 - 3*alpha*mw2/(4*pi*st2)*( 1d0 + 1/(2*ct4) )
  endif

  if (complex_mass_scheme.eq.0) then
    dZh  = real(dZh ,kind=dp)*c1d0
    dmh2 = real(dmh2,kind=dp)*c1d0
  endif

  if (.not.loopWEAK) then
    dZh  = c0d0
    dmh2 = c0d0
  endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZg
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! gluon field
! See hep-ph/0211352
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Gluonic contribution
  if (reg_soft.eq.1) then
    dZg = 5*alphas/(4*pi)*( Duv - Dir + log(muUV2/muIR2) )
  else
    dZg = 5*alphas/(4*pi)*( Duv + log(muUV2/lambda**2) )
  endif

  ! Fermionic contribution.
  ! The fermionic contribution is IR divergent in the massless
  ! limit; we use a fermion mass regulator/dimensional
  ! regularization
  do i = 1,3
    B00uu(i) = B0 (0d0,xu2(i),xu2(i),Duv,muUV2,Dir,muIR2)
    B00dd(i) = B0 (0d0,xd2(i),xd2(i),Duv,muUV2,Dir,muIR2)
    dZg = dZg - alphas/(6*pi)*B00uu(i)
    dZg = dZg - alphas/(6*pi)*B00dd(i)
  enddo

  if (reguScheme.eq.1) then
    dZg = dZg + alphas/(12*pi)*Nc
  endif

   dZg = real(dZg,kind=dp)*c1d0

  if (.not.loopWEAK) dZg = c0d0

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZgs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! The strong coupling constant gs is renormalized introducing the 
! counterterm dZgs:
!   gs^0 = gs * ( 1 + dZgs )
! where gs^0 and gs are respectively the bare and the renormalized 
! strong coupling constants.
! dZgs is computed through the computation of 1-loop gluon 
! self-energies and depends on the flavour scheme, which is set by 
! the user through the variable "Nfren". According to the value of 
! "Nfren", recola computes the number "Nlq" of light quarks for 
! renormalization. Then
! - For loops with gluons and light quarks, the renormalization is 
!   done in MSbar scheme, where however 
!     DeltaUV = 1/ep - EulerConstant + ln(4*pi)
!   is replaced by 
!     DeltaUV -> DeltaUV + log(muUV**2/Qren**2).
!   In this way also the MSbar counterterm gets a muUV dependence 
!   which cancels the one coming from the corresponding pure-loop 
!   contribution. These contributions to dZgs depend then also on the 
!   renormalization scale Qren.
! - For loops with heavy quarks, the renormalization is done at 
!   p^2=0 and the dependence on the quark mass is kept. These 
!   contributions to dZgs depend on muUV and do not depend on the 
!   renormalizationscale Qren.
! Therefore also alphas renormalization is done such that the 
! muUV-independence of the complete 1-loop amplitude can be 
! numerically checked. For alphas, the MSbar running is expressed in 
! terms of Qren (instead of muUV).
! See hep-ph/0211352 for the "pure" MSbar renormalization.
! If Qren=muUV we have alphas in MSbar scheme as in hep-ph/0211352
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Gluonic contribution
  dZgs = - alphas/(8*pi)*11 * ( Duv - log(Qren0**2/muUV2) )

  ! Fermionic contribution
  ! Light flavours:
  dZgs = dZgs + alphas/(12*pi)*Nlq0*( Duv - log(Qren0**2/muUV2) )
  ! Heavy flavours:
  do i = Nlq0+1,Nq
    if (complex_mass_scheme.eq.1) then; m2 = mq2(i)*c1d0
    else;                               m2 = real(mq2(i),kind=dp)*c1d0
    endif
    dZgs = dZgs + alphas/(12*pi)*( Duv - log(m2/muUV2) )
  enddo

  if (reguScheme.eq.1) then
    dZgs = dZgs + alphas/(24*pi)*Nc
  endif

  dZgs = real(dZgs,kind=dp)*c1d0

  if (.not.loopWEAK) dZgs = c0d0

  dZgs0 = dZgs
  if (.not.allocated(dZgs0R)) then
    allocate(dZgs0R(prTot))
    dZgs0R = dZgs0
  endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZaa
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZaa is computed without any expansion, see eq.(4.7)
! of hep-ph/0505042.
! We use DB0(0d0,m**2,m**2,muUV2,muIR2) = 1/(6*m**2) if m not 0
! We use DB0(0d0,m**2,m**2,muUV2,muIR2) = 0          if m  =  0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Bosonic contributions

  B00WW = B0 (0d0,mw2,mw2,Duv,muUV2)

  dZaa = - alpha/(4*pi)*( - 3*B00WW - 2/3d0 )

  ! Fermionic contribution.
  ! The fermionic contribution is IR divergent in the massless
  ! limit; we can use a fermion mass regulator, or dimensional
  ! regularization (note the explicit separation between UV
  ! and IR poles)

  do i = 1,3

    B00ll(i) =  B0 (0d0,xl2(i),xl2(i),Duv,muUV2,Dir,muIR2)
    B00uu(i) =  B0 (0d0,xu2(i),xu2(i),Duv,muUV2,Dir,muIR2)
    B00dd(i) =  B0 (0d0,xd2(i),xd2(i),Duv,muUV2,Dir,muIR2)
    dZaa = dZaa                &
           - alpha/(3*pi)*(    &
             +    Ql2*B00ll(i) &
             + Nc*Qu2*B00uu(i) &
             + Nc*Qd2*B00dd(i) &
             )

    if (check_Pole) then
      DB00ll(i) = DB0 (0d0,xl2(i),xl2(i))
      DB00uu(i) = DB0 (0d0,xu2(i),xu2(i))
      DB00dd(i) = DB0 (0d0,xd2(i),xd2(i))
      dZaa = dZaa                                      &
             - alpha/(3*pi)*(                          &
               +    Ql2*( 2*ml2(i)*DB00ll(i) - 1/3d0 ) &
               + Nc*Qu2*( 2*mu2(i)*DB00uu(i) - 1/3d0 ) &
               + Nc*Qd2*( 2*md2(i)*DB00dd(i) - 1/3d0 ) &
               )
    endif

  enddo

  if (reguScheme.eq.1) then
    dZaa = dZaa + alpha/(6*pi) 
  endif

  if (complex_mass_scheme.eq.0) then
    dZaa = real(dZaa,kind=dp)*c1d0
  endif

  if (.not.loopWEAK) dZaa = c0d0

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZza, dZaz
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZza is computed without any expansion; it gets only bosonic
! contributions (see eq.(B.2) of 0709.1075 for seaz(k^2=0));
! see eq.(4.10) of hep-ph/0505042
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZaz is computed through an expansion around the squared
! real Z mass, and neglecting O(alpha^2) terms, see eq.(4.10)
! of hep-ph/0505042
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Bosonic contributions to seaz0 (A-Z self-energy at p^2=0) and 
  !                          seazz (A-Z self-energy at p^2=Mz^2)

  B0ZWW = B0 (rmz2,mw2,mw2,Duv,muUV2)

  seaz0 = alpha/(2*ct*pi*st)*( AW - mw2 )

  seazz = alpha/(12*pi)*ct/st/mw2*(                               & ! W
          + ( 4*mw2*( mz2 + 3*mw2 ) + rmz2*( 9*mw2 + mz2/2d0 ) )* & ! W
            B0ZWW                                                 & ! W
          + 2*( mz2 - 6*mw2 )*( AW - mw2 )                        & ! W
          + mz2*rmz2/3d0                                          & ! W
          )

  ! Fermionic contributions seazz (A-Z self-energy at p^2=Mz^2)

  do i = 1,3

    B0Zll(i) = B0 (rmz2,ml2(i),ml2(i),Duv,muUV2)
    B0Zuu(i) = B0 (rmz2,mu2(i),mu2(i),Duv,muUV2)
    B0Zdd(i) = B0 (rmz2,md2(i),md2(i),Duv,muUV2)

    seazz = seazz                                &
            - alpha/(18*pi)*(                    &
              + Ql*( gpl + gml )*(               &
                + 3*( rmz2 + 2*ml2(i) )*B0Zll(i) &
                - 6*( Al(i) - ml2(i) )           &
                - rmz2                           &
                )                                &
              + Nc*Qu*( gpu + gmu )*(            &
                + 3*( rmz2 + 2*mu2(i) )*B0Zuu(i) &
                - 6*( Au(i) - mu2(i) )           &
                - rmz2                           &
                )                                &
              + Nc*Qd*( gpd + gmd )*(            &
                + 3*( rmz2 + 2*md2(i) )*B0Zdd(i) &
                - 6*( Ad(i) - md2(i) )           &
                - rmz2                           &
                )                                &
              )

  enddo

  dZza = 2/mz2*seaz0

  dZaz = - 2/rmz2*seazz + ( mz2 - rmz2 )/rmz2*dZza

  if (reguScheme.eq.1) then
    dZaz = dZaz - alpha*ct/(3*pi*st) 
  endif

  if (complex_mass_scheme.eq.0) then
    dZza = real(dZza,kind=dp)*c1d0
    dZaz = real(dZaz,kind=dp)*c1d0
  endif

  if (.not.loopWEAK) then
    dZza = c0d0
    dZaz = c0d0
  endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZzz, dmz2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZzz is computed through an expansion around the squared
! real Z mass, and neglecting O(alpha^2) terms, see eq.(4.10)
! of hep-ph/0505042
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dmz2 is computed through an expansion around the squared
! real Z mass, and neglecting O(alpha^3) terms, see eq.(4.9)
! of hep-ph/0505042
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Bosonic contributions to sezz (Z self-energy at p^2=Mz^2) and dZzz

   B0ZWW =  B0 (rmz2, mw2, mw2,Duv,muUV2)
  DB0ZWW = DB0 (rmz2, mw2, mw2)

   B0ZHZ =  B0 (rmz2, mh2, mz2,Duv,muUV2)
  DB0ZHZ = DB0 (rmz2, mh2, mz2)

  sezz = - alpha/(48*pi*st2)*(                                    &
           + ct2/mw4*(                                            & ! W
             + (                                                  & ! W
               + rmz2*( 36*mw4 + 4*mw2*mz2 - mz4 )                & ! W
               + 4*mw2*( 12*mw4 + 8*mw2*mz2 - 5*mz4 )             & ! W
               )*B0ZWW                                            & ! W
             - 4*( 12*mw4 - 4*mw2*mz2 + mz4 )*( AW - mw2 )        & ! W
             + 2/3d0*mz2*rmz2*( 4*mw2 - mz2 )                     & ! W
             )                                                    & ! W
           + 1/ct2*(                                              & ! H
             + ( 2*mh2 - 10*mz2 - rmz2 - ( mh2 - mz2 )**2/rmz2 )* & ! H
               B0ZHZ                                              & ! H
             - 2*( AH - mh2 ) - 2*( AZ - mz2 )                    & ! H
             + ( mh2 - mz2 )/rmz2*( AH - AZ )                     & ! H
             - 2/3d0*rmz2                                         & ! H
             )                                                    & ! H
           )

  dZzz = - alpha/(48*pi*st2)*(                                    &
           + ct2/mw4*(                                            & ! W
             + ( mz4 - 4*mw2*mz2 - 36*mw4 )*B0ZWW                 & ! W
             + (                                                  & ! W
               + rmz2*( mz4 - 4*mw2*mz2 - 36*mw4 )                & ! W
               + 4*mw2*( 5*mz4 - 8*mw2*mz2 - 12*mw4 )             & ! W
               )*DB0ZWW                                           & ! W
             + 2/3d0*mz2*( mz2 - 4*mw2 )                          & ! W
             )                                                    & ! W
           + 1/ct2*(                                              & ! H
             + ( 1d0 - ( mh2 - mz2 )**2/rmz4 )*B0ZHZ              & ! H
             - ( 2*mh2 - 10*mz2 - rmz2 - ( mh2 - mz2 )**2/rmz2 )* & ! H
               DB0ZHZ                                             & ! H
             + ( mh2 - mz2 )/rmz4*( AH - AZ )                     & ! H
             + 2/3d0                                              & ! H
             )                                                    & ! H
           )

  ! Fermionic contributions to sezz (Z self-energy at p^2=Mz^2) and dZzz

   B0Z00 =  B0 (rmz2,c0d0,c0d0,Duv,muUV2)
  DB0Z00 = DB0 (rmz2,c0d0,c0d0)

  do i = 1,3

     B0Zll(i) =  B0 (rmz2,ml2(i),ml2(i),Duv,muUV2)
    DB0Zll(i) = DB0 (rmz2,ml2(i),ml2(i))

     B0Zuu(i) =  B0 (rmz2,mu2(i),mu2(i),Duv,muUV2)
    DB0Zuu(i) = DB0 (rmz2,mu2(i),mu2(i))

     B0Zdd(i) =  B0 (rmz2,md2(i),md2(i),Duv,muUV2)
    DB0Zdd(i) = DB0 (rmz2,md2(i),md2(i))

    sezz = sezz                                              &
           - alpha/(6*pi)*(                                  &
             + ( gpn2 + gmn2 )*( - rmz2*B0Z00 + 1/3d0*rmz2 ) &
             + ( gpl2 + gml2 )*(                             &
               - ( rmz2 + 2*ml2(i) )*B0Zll(i)                &
               + 2*( Al(i) - ml2(i) ) + 1/3d0*rmz2           &
               )                                             &
             + 3/(4*st2*ct2)*ml2(i)*B0Zll(i)                 &
             + Nc*(                                          &
               + ( gpu2 + gmu2 )*(                           &
                 - ( rmz2 + 2*mu2(i) )*B0Zuu(i)              &
                 + 2*( Au(i) - mu2(i) ) + 1/3d0*rmz2         &
                 )                                           &
               + 3/(4*st2*ct2)*mu2(i)*B0Zuu(i)               &
               + ( gpd2 + gmd2 )*(                           &
                 - ( rmz2 + 2*md2(i) )*B0Zdd(i)              &
                 + 2*( Ad(i) - md2(i) ) + 1/3d0*rmz2         &
                 )                                           &
               + 3/(4*st2*ct2)*md2(i)*B0Zdd(i)               &
               )                                             &
             )

    dZzz  = dZzz &
            - alpha/(6*pi)*(                                         &
              + ( gpn2 + gmn2 )*( rmz2*DB0Z00 + B0Z00 - 1/3d0 )      &
              + ( gpl2 + gml2 )*                                     &
                ( ( rmz2 + 2*ml2(i) )*DB0Zll(i) + B0Zll(i) - 1/3d0 ) &
              - 3/(4*st2*ct2)*ml2(i)*DB0Zll(i)                       &
              + Nc*(                                                 &
                + ( gpu2 + gmu2 )*(                                  &
                  + ( rmz2 + 2*mu2(i) )*DB0Zuu(i) + B0Zuu(i) - 1/3d0 &
                  )                                                  &
                - 3/(4*st2*ct2)*mu2(i)*DB0Zuu(i)                     &
                + ( gpd2 + gmd2 )*(                                  &
                  + ( rmz2 + 2*md2(i) )*DB0Zdd(i) + B0Zdd(i) - 1/3d0 &
                  )                                                  &
                - 3/(4*st2*ct2)*md2(i)*DB0Zdd(i)                     &
                )                                                    &
              )
  enddo

  dmz2 = sezz - ( mz2 - rmz2 ) * dZzz

  if (reguScheme.eq.1) then
    dZzz = dZzz + alpha*ct2/(6*pi*st2)
    dmz2 = dmz2 - alpha*ct2*mz2/(6*pi*st2)
  endif

  if (complex_mass_scheme.eq.0) then
    dZzz = real(dZzz,kind=dp)*c1d0
    dmz2 = real(dmz2,kind=dp)*c1d0
  endif

  if (.not.loopWEAK) then
    dZzz = c0d0
    dmz2 = c0d0
  endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZw, dmw2, dmw
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZw is computed through an expansion about the squared
! real W mass, and neglecting O(alpha^2) terms, see eq.(4.10)
! of hep-ph/0505042
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dmw2 is computed through an expansion around the squared
! real W mass, and neglecting O(alpha^3) terms, see eq.(4.9)
! of hep-ph/0505042.
! The correction of the complex mass scheme according to eq.(4.29)
! of hep-ph/0505042 is implemented.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dmw = dmw2/(2*mw)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Bosonic contributions to seww (W self-energy at p^2=Mw^2) and dZw

  B0W0W  = B0 (rmw2,c0d0,mw2,Duv,muUV2)
  if (reg_soft.eq.1) then
    ! If mw2 not = rmw2, 
    !   DB0W0W is not IR divergent and muIR2 is irrelevant
    ! If mw2   =   rmw2, 
    !   DB0W0W = - 1/2/mw2*( Dir + log(muIR2/mw2) + 2 )
    DB0W0W = DB0 (rmw2,c0d0,mw2,Dir,muIR2)
  else
    ! If mw2 not = rmw2, 
    !   DB0W0W is not IR divergent and lambda is irrelevant
    ! If mw2   =   rmw2, 
    !   DB0W0W = - 1/2/mw2*( log(lambda/mw2) + 2 )
    DB0W0W = DB0 (rmw2,c0d0,mw2,0d0,lambda**2)
  endif

   B0WWZ =  B0 (rmw2,mw2,mz2,Duv,muUV2)
  DB0WWZ = DB0 (rmw2,mw2,mz2)

   B0WHW =  B0 (rmw2,mh2,mw2,Duv,muUV2)
  DB0WHW = DB0 (rmw2,mh2,mw2)

  seww = - alpha/(12*pi)*(                                        &
           + 2*( 2*mw2 + 5*rmw2 - mw4/rmw2 )*B0W0W                & ! A
           + 2*( mw2/rmw2 - 2 )*AW                                & ! A
           + 4*mw2 + 2/3d0*rmw2                                   & ! A
           + ct2/(4*st2)/mw2*(                                    & ! Z
             + (                                                  & ! Z
               + ( 40*mw2 - mz2 )*rmw2 + 16*mw4 + 54*mw2*mz2      & ! Z
               - 10*mz4 - ( 8*mw2 + mz2 )*( mw2 - mz2 )**2/rmw2   & ! Z
               )*B0WWZ                                            & ! Z
             + ( 8*mw2 + mz2 )*( mw2 - mz2 )/rmw2*( AW - AZ )     & ! Z
             - 2*( 8*mw2 + mz2 )*( AW - mw2 + AZ - mz2 )          & ! Z
             + 2/3d0*rmw2*( 4*mw2 - mz2 )                         & ! Z
             )                                                    & ! Z
           + 1/(4*st2)*(                                          & ! H
             + ( 2*mh2 - 10*mw2 - rmw2 - ( mw2 - mh2 )**2/rmw2 )* & ! H
               B0WHW                                              & ! H
             + ( mw2 - mh2 )/rmw2*( AW - AH )                     & ! H
             - 2*( AW - mw2 + AH - mh2 )                          & ! H
             - 2/3d0*rmw2                                         & ! H
             )                                                    & ! H
           )
  dZw = - alpha/(12*pi)*(                                        &
          - 2*( 5d0 + mw4/rmw4 )*B0W0W                           & ! A
          - 2*( 5*rmw2 + 2*mw2 - mw4/rmw2 )*DB0W0W               & ! A
          + 2*mw2/rmw4*AW                                        & ! A
          - 2/3d0                                                & ! A
          + ct2/(4*st2)/mw2*(                                    & ! Z
            - (                                                  & ! Z
              + 40*mw2 - mz2                                     & ! Z
              + ( 8*mw2 + mz2 )*( mw2 - mz2 )**2/rmw4            & ! Z
              )*B0WWZ                                            & ! Z
            + (                                                  & ! Z
              + ( mz2 - 40*mw2 )*rmw2 - 16*mw4 - 54*mw2*mz2      & ! Z
              + 10*mz4 + ( 8*mw2 + mz2 )*( mw2 - mz2 )**2/rmw2   & ! Z
              )*DB0WWZ                                           & ! Z
            + ( mz2 + 8*mw2 )*( mw2 - mz2 )/rmw4*( AW - AZ )     & ! Z
            + 2/3d0*( mz2 - 4*mw2 )                              & ! Z
            )                                                    & ! Z
          + 1/(4*st2)*(                                          & ! H
            + ( 1d0 - ( mw2 - mh2 )**2/rmw4 )*B0WHW              & ! H
            + ( rmw2 + 10*mw2 - 2*mh2 + ( mw2 - mh2 )**2/rmw2 )* & ! H
              DB0WHW                                             & ! H
            + ( mw2 - mh2 )/rmw4*( AW - AH )                     & ! H
            + 2/3d0                                              & ! H
            )                                                    & ! H
          )

  ! Fermionic contributions to seww (W self-energy at p^2=Mw^2) and dZw

  do i = 1,3

     B0W0l(i) =  B0 (rmw2,  c0d0,ml2(i),Duv,muUV2)
    DB0W0l(i) = DB0 (rmw2,  c0d0,ml2(i))

     B0Wud(i) =  B0 (rmw2,mu2(i),md2(i),Duv,muUV2)
    DB0Wud(i) = DB0 (rmw2,mu2(i),md2(i))

    seww = seww &
           - alpha/(12*pi*st2)*(                                &
             - ( rmw2 - ml2(i)/2d0 - ml4(i)/(2*rmw2) )*B0W0l(i) &
             - ml2(i)/(2*rmw2)*Al(i)                            &
             + Al(i) - ml2(i)                                   &
             + 1/3d0*rmw2                                       &
             + Nc*(                                             &
               - (                                              &
                 + rmw2 - ( mu2(i) + md2(i) )/2d0               &
                 - ( mu2(i) - md2(i) )**2/(2*rmw2)              &
                 )*B0Wud(i)                                     &
               - ( mu2(i) - md2(i) )/(2*rmw2)*( Au(i) - Ad(i) ) &
               + Au(i) - mu2(i) + Ad(i) - md2(i)                &
               + 1/3d0*rmw2                                     &
               )                                                &
             ) 

    dZw = dZw &
          - alpha/(12*pi*st2)*(                                    &
            + ( 1d0 + ml4(i)/(2*rmw4) )*B0W0l(i)                   &
            + ( rmw2 - ml2(i)/2d0 - ml4(i)/(2*rmw2) )*DB0W0l(i)    &
            - ml2(i)/(2*rmw4)*Al(i)                                &
            - 1/3d0                                                &
            + Nc*(                                                 &
              + ( 1d0 + ( mu2(i) - md2(i) )**2/(2*rmw4) )*B0Wud(i) &
              + (                                                  &
                + rmw2 - ( mu2(i) + md2(i) )/2d0                   &
                - ( mu2(i) - md2(i) )**2/(2*rmw2)                  &
                )*DB0Wud(i)                                        &
              - ( mu2(i) - md2(i) )/(2*rmw4)*( Au(i) - Ad(i) )     &
              - 1/3d0                                              &
              )                                                    &
            ) 

  enddo

  dmw2 = seww - ( mw2 - rmw2 )*( dZw + alpha/pi)

  if (reguScheme.eq.1) then
    dZw  = dZw  + alpha/(6*pi*st2)
    dmw2 = dmw2 - alpha*mw2/(6*pi*st2)
  endif

  if (complex_mass_scheme.eq.0) then
    dZw  = real(dZw ,kind=dp)*c1d0
    dmw2 = real(dmw2,kind=dp)*c1d0
  endif

  if (.not.loopWEAK) then
    dZw  = c0d0
    dmw2 = c0d0
  endif

  dmw = dmw2/(2*mw1)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZe
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  select case (ew_reno_scheme)

  case (1) ! alpha_gf scheme

    ! Bosonic contributions to sew0 (W self-energy at p^2=0)

     B000W =  B0 ( 0d0,c0d0, mw2,Duv,muUV2)
     B00WZ =  B0 ( 0d0, mw2, mz2,Duv,muUV2)
     B00HW =  B0 ( 0d0, mh2, mw2,Duv,muUV2)
    DB000W = DB0 ( 0d0,c0d0, mw2)
    DB00WZ = DB0 ( 0d0, mw2, mz2)
    DB00HW = DB0 ( 0d0, mh2, mw2)

    sew0 = - alpha/(12*pi)*(                               &
             + 4*mw2*B000W                                 & ! A
             - 4*( AW - mw2 )                              & ! A
             - 2*mw4*DB000W                                & ! A
             + ct2/(4*st2)/mw2*(                           & ! Z
               + 2*( 8*mw4 + 27*mw2*mz2 - 5*mz4 )*B00WZ    & ! Z
               - 2*( 8*mw2 + mz2 )*( AW - mw2 + AZ - mz2 ) & ! Z
               - ( 8*mw2 + mz2 )*( mw2 - mz2 )**2*DB00WZ   & ! Z
               )                                           & ! Z
             + 1/(4*st2)*(                                 & ! H
               + 2*( mh2 - 5*mw2 )*B00HW                   & ! H
               - ( mw2 - mh2 )**2*DB00HW                   & ! H 
               - 2*( AW - mw2 + AH - mh2 )                 & ! H
               )                                           & ! H
             )

    ! Fermionic contributions to sew0 (W self-energy at p^2=0)

    do i = 1,3
      if (regl(i).eq.3) then
         B000l(i) =  B0 ( 0d0,  c0d0,ml2(i),Duv,muUV2)
        DB000l(i) = DB0 ( 0d0,  c0d0,ml2(i))
        sew0 = sew0                       &
               - alpha/(12*pi*st2)*(      &
                 + 1/2d0*ml2(i)*B000l(i)  &
                 + Al(i) - ml2(i)         &
                 + 1/2d0*ml4(i)*DB000l(i) &
                 ) 
      endif
      if (regd(i).eq.3.or.regu(i).eq.3) then
         B00ud(i) =  B0 ( 0d0,mu2(i),md2(i),Duv,muUV2)
        DB00ud(i) = DB0 ( 0d0,mu2(i),md2(i))
        sew0 = sew0 &
               - alpha/(12*pi*st2)*Nc*(                   &
                 + 1/2d0*( mu2(i) + md2(i) )*B00ud(i)     &
                 + Au(i) - mu2(i) + Ad(i) - md2(i)        &
                 + 1/2d0*( mu2(i) - md2(i) )**2*DB00ud(i) &
                 ) 
      endif
    enddo

    deltar = - dZaa + ct/st*dZza + 1/mw2*sew0                   &
             - ct2/(st2*mz2)*dmz2 + ( ct2/st2 - 1d0 )/mw2*dmw2  &
             + alpha/(4*pi*st2)*                                &
               ( 6d0 + 1/(2*st2)*( 7d0 - 4*st2 )*log(ct2) )

    dZe = - dZaa/2d0 - dZza*st/(2*ct)

    if (.not.check_Pole) dZe = dZe - deltar/2d0

    ! strange UV-IR poles cancel

  case (2) ! alpha_0 scheme

    dZe = - dZaa/2d0 - dZza*st/(2*ct) ! light fermion masses

    ! contains IR poles through dZaa

  case (3) ! alpha_Z scheme, dalZ is real

    dalZ = 0d0

    do i = 1,3
      if (rxl2(i).lt.rmz2) then
        B00ll(i) = B0 ( 0d0,xl2(i),xl2(i),Duv,muUV2,Dir,muIR2)
        B0Zll(i) = B0 (rmz2,ml2(i),ml2(i),Duv,muUV2)
        dalZ = dalZ                                                &
               + real(                                             &
                 alpha/(3*pi)*Ql2*(                                &
                 + 1/3d0                                           &
                 + ( 1d0 + 2*ml2(i)/rmz2 )*( B00ll(i) - B0Zll(i) ) &
                 )                                                 &
                 ,kind=dp)
      endif
      if (rxu2(i).lt.rmz2) then
        B00uu(i) = B0 ( 0d0,xu2(i),xu2(i),Duv,muUV2,Dir,muIR2)
        B0Zuu(i) = B0 (rmz2,mu2(i),mu2(i),Duv,muUV2)
        dalZ = dalZ                                                &
               + real(                                             &
                 alpha/(3*pi)*Nc*Qu2*(                             &
                 + 1/3d0                                           &
                 + ( 1d0 + 2*mu2(i)/rmz2 )*( B00uu(i) - B0Zuu(i) ) &
                 )                                                 &
                 ,kind=dp)
      endif
      if (rxd2(i).lt.rmz2) then
        B00dd(i) = B0 ( 0d0,xd2(i),xd2(i),Duv,muUV2,Dir,muIR2)
        B0Zdd(i) = B0 (rmz2,md2(i),md2(i),Duv,muUV2)
        dalZ = dalZ                                                &
               + real(                                             &
                 alpha/(3*pi)*Nc*Qd2*(                             &
                 + 1/3d0                                           &
                 + ( 1d0 + 2*md2(i)/rmz2 )*( B00dd(i) - B0Zdd(i) ) &
                 )                                                 &
                 ,kind=dp)
      endif
    enddo

    dZe = - dZaa/2d0 - dZza*st/(2*ct) - dalZ/2d0

    ! the IR poles of dalZ cancel with the IR poles of dZaa

  end select

  if (complex_mass_scheme.eq.0) then
    dZe  = real(dZe ,kind=dp)*c1d0
  endif

  if (.not.loopWEAK) dZe = c0d0

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dst and dct
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  dct = ct/2d0*( dmw2/mw2 - dmz2/mz2 )

  dst = - ct/st*dct

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dgpn, dgpl, dgpu, dgpd, dgmn, dgmu, dgmu, dgmd
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  dgpn = c0d0
  dgpl = - Ql*st/ct*( dZe + dst/(st*ct2) )
  dgpu = - Qu*st/ct*( dZe + dst/(st*ct2) )
  dgpd = - Qd*st/ct*( dZe + dst/(st*ct2) )

  dgmn = dgpn + I3n/stct*( dZe + dst*(st2-ct2)/(st*ct2) )
  dgml = dgpl + I3l/stct*( dZe + dst*(st2-ct2)/(st*ct2) )
  dgmu = dgpu + I3u/stct*( dZe + dst*(st2-ct2)/(st*ct2) )
  dgmd = dgpd + I3d/stct*( dZe + dst*(st2-ct2)/(st*ct2) )

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZnL, dZnR
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! neutrino field counterterms
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  do i = 1,3

    B100Z(i) = B1 (0d0,  c0d0,mz2,Duv,muUV2)
    B10lW(i) = B1 (0d0,ml2(i),mw2,Duv,muUV2)

    seL    = - alpha/(4*pi)*(                                      &
               + gmn2*( 2*B100Z(i) + 1d0 )                         &
               + 1/(2*st2)*( ( 2d0 + ml2(i)/mw2 )*B10lW(i) + 1d0 ) &
               )

    seR    = - alpha/(4*pi)*gpn2*( 2*B100Z(i) + 1d0 )

    dZnL(i) = - seL
    dZnR(i) = - seR

    if (reguScheme.eq.1) &
      dZnL(i) = dZnL(i) - alpha/(4*pi*st2)*( I3n2/ct2 + 1/2d0 )

    if (complex_mass_scheme.eq.0) then
      dZnL(i) = real(dZnL(i),kind=dp)*c1d0
      dZnR(i) = real(dZnR(i),kind=dp)*c1d0
    endif

  if (.not.loopWEAK) then
    dZnL(i) = c0d0
    dZnR(i) = c0d0
  endif

    dZnLc(i) = dZnL(i)
    dZnRc(i) = dZnR(i)

  enddo

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZlL, dZlR, dml
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! lepton field and mass counterterms
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  do i = 1,3

    B1llZ(i) = B1 (rml2(i),ml2(i), mz2,Duv,muUV2)
    B1l0W(i) = B1 (rml2(i),  c0d0, mw2,Duv,muUV2)
    B1llH(i) = B1 (rml2(i),ml2(i), mh2,Duv,muUV2)
    B1ll0(i) = B1 (rxl2(i),xl2(i),c0d0,Duv,muUV2,Dir,muIR2)

    seL    = - alpha/(4*pi)*(                               &
               + gml2*( 2*B1llZ(i) + 1d0 )                  &
               + ml2(i)/(4*st2*mw2)*( B1llZ(i) + B1llH(i) ) &
               + 1/(2*st2)*( 2*B1l0W(i) + 1d0 )             &
               )
    seLqed = - alpha/(4*pi)*Ql2*( 2*B1ll0(i) + 1d0 )

    seR    = - alpha/(4*pi)*(                               &
               + gpl2*( 2*B1llZ(i) + 1d0 )                  &
               + ml2(i)/(4*st2*mw2)*( B1llZ(i) + B1llH(i) ) &
               + ml2(i)/(2*st2*mw2)*B1l0W(i)                &
               )
    seRqed = - alpha/(4*pi)*Ql2*( 2*B1ll0(i) + 1d0 )

    if (regl(i).eq.1) then ! massless fermion

      m2dseL    = c0d0
      m2dseLqed = c0d0
      m2dseR    = c0d0
      m2dseRqed = c0d0
      m2dseS    = c0d0
      m2dseSqed = c0d0

    elseif (regl(i).eq.2) then ! light fermion

      if (reg_soft.eq.2) then
        logIR = log(lambda**2/xl2(i)) ! lambda**2 << xl2(i)
      else
        logIR = Dir + log(muIR2/xl2(i))
      endif
      m2dseL    = c0d0
      m2dseLqed = - alpha/(4*pi)*Ql2*( logIR + 3d0 )
      m2dseR    = c0d0
      m2dseRqed = - alpha/(4*pi)*Ql2*( logIR + 3d0 )
      m2dseS    = c0d0
      m2dseSqed = + alpha/(2*pi)*Ql2*( logIR + 2d0 )

    else ! massive fermion

       B0llZ(i) =  B0 (rml2(i),ml2(i), mz2,Duv,muUV2)
      DB0llZ(i) = DB0 (rml2(i),ml2(i), mz2)
      DB1llZ(i) = DB1 (rml2(i),ml2(i), mz2)

       B0l0W(i) =  B0 (rml2(i),  c0d0, mw2,Duv,muUV2)
      DB0l0W(i) = DB0 (rml2(i),  c0d0, mw2)
      DB1l0W(i) = DB1 (rml2(i),  c0d0, mw2)

       B0llH(i) =  B0 (rml2(i),ml2(i), mh2,Duv,muUV2)
      DB0llH(i) = DB0 (rml2(i),ml2(i), mh2)
      DB1llH(i) = DB1 (rml2(i),ml2(i), mh2)

       B0ll0(i) =  B0 (rml2(i),ml2(i),c0d0,Duv,muUV2)
      ratioIm = 1d0 - abs(rxl2(i))/abs(xl2(i))
      if (abs(ratioIm).le.zerocut) then ! for a vanishing width, IR poles show up
        if (reg_soft.eq.1) then
          DB0ll0(i) = - 1/2d0/xl2(i)*( Dir + log(muIR2/xl2(i)) + 2d0 )
          DB1ll0(i) = + 1/2d0/xl2(i)*( Dir + log(muIR2/xl2(i)) + 3d0 )
        else
          DB0ll0(i) = - 1/2d0/xl2(i)*( log(lambda**2/xl2(i)) + 2d0 ) ! lambda**2 << xl2(i)
          DB1ll0(i) = + 1/2d0/xl2(i)*( log(lambda**2/xl2(i)) + 3d0 ) ! lambda**2 << xl2(i)
        endif
      else
        DB0ll0(i) = DB0 (rml2(i),ml2(i),c0d0)
        DB1ll0(i) = DB1 (rml2(i),ml2(i),c0d0)
      endif

      seS    = - alpha/(4*pi)*(                               &
                 + 2*gpl*gml*( 2*B0llZ(i) - 1d0 )             &
                 + ml2(i)/(4*st2*mw2)*( B0llZ(i) - B0llH(i) ) &
                 )
      seSqed = - alpha/(2*pi)*Ql2*( 2*B0ll0(i) - 1d0 )

      m2dseL    = - alpha/(4*pi)*rml2(i)*(                         &
                    + 2*gml2*DB1llZ(i)                             &
                    + ml2(i)/(4*st2*mw2)*( DB1llZ(i) + DB1llH(i) ) &
                    + 1/st2*DB1l0W(i)                              &
                    )
      m2dseLqed = - alpha/(2*pi)*Ql2*rxl2(i)*DB1ll0(i)

      m2dseR    = - alpha/(4*pi)*rml2(i)*(                         &
                    + 2*gpl2*DB1llZ(i)                             &
                    + ml2(i)/(4*st2*mw2)*( DB1llZ(i) + DB1llH(i) ) &
                    + ml2(i)/(2*st2*mw2)*DB1l0W(i)                 &
                    )
      m2dseRqed = - alpha/(2*pi)*Ql2*rxl2(i)*DB1ll0(i)

      m2dseS    = - alpha/(4*pi)*rml2(i)*(                         &
                    + 4*gpl*gml*DB0llZ(i)                          &
                    + ml2(i)/(4*st2*mw2)*( DB0llZ(i) - DB0llH(i) ) &
                    )
      m2dseSqed = - alpha/pi*Ql2*rxl2(i)*DB0ll0(i)

    endif

    dZlL(i) = c0d0
    if (loopWEAK) &
    dZlL(i) = dZlL(i) - seL - m2dseL - m2dseR - 2*m2dseS 
    if (loopQED) &
    dZlL(i) = dZlL(i) - seLqed - m2dseLqed - m2dseRqed - 2*m2dseSqed

    dZlR(i) = c0d0
    if (loopWEAK) &
    dZlR(i) = dZlR(i) - seR - m2dseL - m2dseR - 2*m2dseS
    if (loopQED) &
    dZlR(i) = dZlR(i) - seRqed - m2dseLqed - m2dseRqed - 2*m2dseSqed

    if (reguScheme.eq.1) then
      if (loopWEAK) dZlL(i) = dZlL(i) - alpha/(4*pi)*( gml2 + 1/(2*st2) )
      if (loopQED)  dZlL(i) = dZlL(i) - alpha/(4*pi)*Ql2
      if (loopWEAK) dZlR(i) = dZlR(i) - alpha/(4*pi)*gpl2
      if (loopQED)  dZlR(i) = dZlR(i) - alpha/(4*pi)*Ql2
    endif

    if (complex_mass_scheme.eq.0) then
      dZlL(i) = real(dZlL(i),kind=dp)*c1d0
      dZlR(i) = real(dZlR(i),kind=dp)*c1d0
    endif

    dZlLc(i) = dZlL(i)
    dZlRc(i) = dZlR(i)

    if (regl(i).le.2) then ! massless and light fermion

      dml(i)    = c0d0

    else ! massive fermion

      dml(i) = c0d0
      if (loopWEAK) &
      dml(i) = dml(i) &
               + ml1(i)/2d0*(                      &
                 + seL + seR + 2*seS               &
                 + ( ml2(i) - rml2(i) ) / rml2(i)* &
                   ( m2dseL + m2dseR + 2*m2dseS )  &
                 )
      if (loopQED) &
      dml(i) = dml(i) &
               + ml1(i)/2d0*(                              &
                 + seLqed + seRqed + 2*seSqed              &
                 + ( ml2(i) - rml2(i) ) / rml2(i)*         &
                   ( + m2dseLqed + m2dseRqed + 2*m2dseSqed &
                     - alpha/pi*Ql2 )                      &
                 )

      if (reguScheme.eq.1) then
        if (loopWEAK) &
        dml(i) = dml(i) + alpha/(8*pi)*ml1(i)*(                 &
                          + gpl2 + gml2 + 4*gpl*gml + 1/(2*st2) &
                          )
        if (loopQED) &
        dml(i) = dml(i) + alpha/(4*pi)*ml1(i)*3*Ql2
      endif

    endif

    if (complex_mass_scheme.eq.0) then
      dml(i) = real(dml(i),kind=dp)*c1d0
    endif

  enddo

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZuL, dZuR, dmu
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! up-quark field and mass counterterms
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  do i = 1,3

    B1uuZ(i) = B1 (rmu2(i),mu2(i), mz2,Duv,muUV2)
    B1udW(i) = B1 (rmu2(i),md2(i), mw2,Duv,muUV2)
    B1uuH(i) = B1 (rmu2(i),mu2(i), mh2,Duv,muUV2)
    B1uu0(i) = B1 (rxu2(i),xu2(i),c0d0,Duv,muUV2,Dir,muIR2)

    seL    = - alpha/(4*pi)*(                                      &
               + gmu2*( 2*B1uuZ(i) + 1d0 )                         &
               + mu2(i)/(4*st2*mw2)*( B1uuZ(i) + B1uuH(i) )        &
               + 1/(2*st2)*( ( 2d0 + md2(i)/mw2 )*B1udW(i) + 1d0 ) &
               )
    seLqed = - alpha/(4*pi)*Qu2*( 2*B1uu0(i) + 1d0 )
    seLqcd = - alphas/(4*pi)*Cf*( 2*B1uu0(i) + 1d0 )

    seR    = - alpha/(4*pi)*(                               &
               + gpu2*( 2*B1uuZ(i) + 1d0 )                  &
               + mu2(i)/(4*st2*mw2)*( B1uuZ(i) + B1uuH(i) ) &
               + mu2(i)/(2*st2*mw2)*B1udW(i)                &
               )
    seRqed = - alpha/(4*pi)*Qu2*( 2*B1uu0(i) + 1d0 )
    seRqcd = - alphas/(4*pi)*Cf*( 2*B1uu0(i) + 1d0 )

    if (regu(i).eq.1) then ! massless fermion

      m2dseL    = c0d0
      m2dseLqed = c0d0
      m2dseLqcd = c0d0
      m2dseR    = c0d0
      m2dseRqed = c0d0
      m2dseRqcd = c0d0
      m2dseS    = c0d0
      m2dseSqed = c0d0
      m2dseSqcd = c0d0

    elseif (regu(i).eq.2) then ! light fermion

      if (reg_soft.eq.2) then
        logIR = log(lambda**2/xu2(i)) ! lambda**2 << xu2(i)
      else
        logIR = Dir + log(muIR2/xu2(i))
      endif
      m2dseL    = c0d0
      m2dseLqed = - alpha/(4*pi)*Qu2*( logIR + 3d0 )
      m2dseLqcd = - alphas/(4*pi)*Cf*( logIR + 3d0 )
      m2dseR    = c0d0
      m2dseRqed = - alpha/(4*pi)*Qu2*( logIR + 3d0 )
      m2dseRqcd = - alphas/(4*pi)*Cf*( logIR + 3d0 )
      m2dseS    = c0d0
      m2dseSqed = + alpha/(2*pi)*Qu2*( logIR + 2d0 )
      m2dseSqcd = + alphas/(2*pi)*Cf*( logIR + 2d0 )

    else ! massive fermion

       B0uuZ(i) =  B0 (rmu2(i),mu2(i), mz2,Duv,muUV2)
      DB0uuZ(i) = DB0 (rmu2(i),mu2(i), mz2)
      DB1uuZ(i) = DB1 (rmu2(i),mu2(i), mz2)

       B0udW(i) =  B0 (rmu2(i),md2(i), mw2,Duv,muUV2)
      DB0udW(i) = DB0 (rmu2(i),md2(i), mw2)
      DB1udW(i) = DB1 (rmu2(i),md2(i), mw2)

       B0uuH(i) =  B0 (rmu2(i),mu2(i), mh2,Duv,muUV2)
      DB0uuH(i) = DB0 (rmu2(i),mu2(i), mh2)
      DB1uuH(i) = DB1 (rmu2(i),mu2(i), mh2)

       B0uu0(i) =  B0 (rmu2(i),mu2(i),c0d0,Duv,muUV2)
      ratioIm = 1d0 - abs(rxu2(i))/abs(xu2(i))
      if (abs(ratioIm).le.zerocut) then ! for a vanishing width, IR poles show up
        if (reg_soft.eq.1) then
          DB0uu0(i) = - 1/2d0/xu2(i)*( Dir + log(muIR2/xu2(i)) + 2d0 )
          DB1uu0(i) = + 1/2d0/xu2(i)*( Dir + log(muIR2/xu2(i)) + 3d0 )
        else
          DB0uu0(i) = - 1/2d0/xu2(i)*( log(lambda**2/xu2(i)) + 2d0 ) ! lambda**2 << xu2(i)
          DB1uu0(i) = + 1/2d0/xu2(i)*( log(lambda**2/xu2(i)) + 3d0 ) ! lambda**2 << xu2(i)
        endif
      else
        DB0uu0(i) = DB0 (rmu2(i),mu2(i),c0d0)
        DB1uu0(i) = DB1 (rmu2(i),mu2(i),c0d0)
      endif

      seS    = - alpha/(4*pi)*(                               &
                 + 2*gpu*gmu*( 2*B0uuZ(i) - 1d0 )             &
                 + mu2(i)/(4*st2*mw2)*( B0uuZ(i) - B0uuH(i) ) &
                 + md2(i)/(2*st2*mw2)*B0udW(i)                &
                 )
      seSqed = - alpha/(2*pi)*Qu2*( 2*B0uu0(i) - 1d0 )
      seSqcd = - alphas/(2*pi)*Cf*( 2*B0uu0(i) - 1d0 )

      m2dseL    = - alpha/(4*pi)*rmu2(i)*(                           &
                    + 2*gmu2*DB1uuZ(i)                             &
                    + mu2(i)/(4*st2*mw2)*( DB1uuZ(i) + DB1uuH(i) ) &
                    + 1/(2*st2)*( 2d0 + md2(i)/mw2 )*DB1udW(i)     &
                    )
      m2dseLqed = - alpha/(2*pi)*Qu2*rxu2(i)*DB1uu0(i)
      m2dseLqcd = - alphas/(2*pi)*Cf*rxu2(i)*DB1uu0(i)

      m2dseR    = - alpha/(4*pi)*rmu2(i)*(                           &
                    + 2*gpu2*DB1uuZ(i)                             &
                    + mu2(i)/(4*st2*mw2)*( DB1uuZ(i) + DB1uuH(i) ) &
                    + mu2(i)/(2*st2*mw2)*DB1udW(i)                 &
                    )
      m2dseRqed = - alpha/(2*pi)*Qu2*rxu2(i)*DB1uu0(i)
      m2dseRqcd = - alphas/(2*pi)*Cf*rxu2(i)*DB1uu0(i)

      m2dseS    = - alpha/(4*pi)*rmu2(i)*(                           &
                    + 4*gpu*gmu*DB0uuZ(i)                          &
                    + mu2(i)/(4*st2*mw2)*( DB0uuZ(i) - DB0uuH(i) ) &
                    + md2(i)/(2*st2*mw2)*DB0udW(i)                 &
                    )
      m2dseSqed = - alpha/pi*Qu2*rxu2(i)*DB0uu0(i)
      m2dseSqcd = - alphas/pi*Cf*rxu2(i)*DB0uu0(i)
    endif

    dZuL(i) = c0d0
    if (loopWEAK) &
    dZuL(i) = dZuL(i) - seL - m2dseL - m2dseR - 2*m2dseS 
    if (loopQED) &
    dZuL(i) = dZuL(i) - seLqed - m2dseLqed - m2dseRqed - 2*m2dseSqed
    dZuLqcd(i) = - seLqcd - m2dseLqcd - m2dseRqcd - 2*m2dseSqcd

    dZuR(i) = c0d0
    if (loopWEAK) &
    dZuR(i) = dZuR(i) - seR - m2dseL - m2dseR - 2*m2dseS
    if (loopQED) &
    dZuR(i) = dZuR(i) - seRqed - m2dseLqed - m2dseRqed - 2*m2dseSqed
    dZuRqcd(i) = - seRqcd - m2dseLqcd - m2dseRqcd - 2*m2dseSqcd

    if (reguScheme.eq.1) then
      if (loopWEAK) &
      dZuL(i)    = dZuL(i)    - alpha/(4*pi)*( gmu2 + 1/(2*st2) )
      if (loopQED) &
      dZuL(i)    = dZuL(i)    - alpha/(4*pi)*Qu2
      dZuLqcd(i) = dZuLqcd(i) - alphas/(4*pi)*Cf
      if (loopWEAK) &
      dZuR(i)    = dZuR(i)    - alpha/(4*pi)*gpu2
      if (loopQED) &
      dZuR(i)    = dZuR(i)    - alpha/(4*pi)*Qu2
      dZuRqcd(i) = dZuRqcd(i) - alphas/(4*pi)*Cf
    endif

    if (complex_mass_scheme.eq.0) then
      dZuL(i) = real(dZuL(i),kind=dp)*c1d0
      dZuR(i) = real(dZuR(i),kind=dp)*c1d0
      dZuLqcd(i) = real(dZuLqcd(i),kind=dp)*c1d0
      dZuRqcd(i) = real(dZuRqcd(i),kind=dp)*c1d0
    endif

    dZuLc   (i) = dZuL   (i)
    dZuLcqcd(i) = dZuLqcd(i)

    dZuRc   (i) = dZuR   (i)
    dZuRcqcd(i) = dZuRqcd(i)

    if (regu(i).le.2) then ! massless and light fermion

      dmu(i)    = c0d0
      dmuqcd(i) = c0d0

    else ! massive fermion

      dmu(i) = c0d0
      if (loopWEAK) &
      dmu(i) = dmu(i) &
               + mu1(i)/2d0*(                      &
                 + seL + seR + 2*seS               &
                 + ( mu2(i) - rmu2(i) ) / rmu2(i)* &
                   ( m2dseL + m2dseR + 2*m2dseS )  &
                 )
      if (loopQED) &
      dmu(i) = dmu(i) &
               + mu1(i)/2d0*(                              &
                 + seLqed + seRqed + 2*seSqed              &
                 + ( mu2(i) - rmu2(i) ) / rmu2(i)*         &
                   ( + m2dseLqed + m2dseRqed + 2*m2dseSqed &
                     - alpha/pi*Qu2 )                      &
                 )
      dmuqcd(i) = - alphas/(4*pi)*mu1(i)*Cf*              &
                    ( 3*(Duv + log(muUV2/mu2(i))) + 4d0 )

      if (reguScheme.eq.1) then
        if (loopWEAK) &
        dmu(i)    = dmu(i) + alpha/(8*pi)*mu1(i)*(          &
                      + gpu2 + gmu2 + 4*gpu*gmu + 1/(2*st2) &
                      )
        if (loopQED) &
        dmu(i)    = dmu(i)    + alpha/(4*pi)*mu1(i)*3*Qu2
        dmuqcd(i) = dmuqcd(i) + alphas/(4*pi)*mu1(i)*3*Cf
      endif

    endif

    if (complex_mass_scheme.eq.0) then
      dmu(i) = real(dmu(i),kind=dp)*c1d0
      dmuqcd(i) = real(dmuqcd(i),kind=dp)*c1d0
    endif

  enddo

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dZdL, dZdR, dmd
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! down-quark field and mass counterterms
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  do i = 1,3

    B1ddZ(i) = B1 (rmd2(i),md2(i),mz2,Duv,muUV2)
    B1duW(i) = B1 (rmd2(i),mu2(i),mw2,Duv,muUV2)
    B1ddH(i) = B1 (rmd2(i),md2(i),mh2,Duv,muUV2)
    B1dd0(i) = B1 (rxd2(i),xd2(i),c0d0,Duv,muUV2,Dir,muIR2)

    seL    = - alpha/(4*pi)*(                                      &
               + gmd2*( 2*B1ddZ(i) + 1d0 )                         &
               + md2(i)/(4*st2*mw2)*( B1ddZ(i) + B1ddH(i) )        &
               + 1/(2*st2)*( ( 2d0 + mu2(i)/mw2 )*B1duW(i) + 1d0 ) &
               )
    seLqed = - alpha/(4*pi)*Qd2*( 2*B1dd0(i) + 1d0 )
    seLqcd = - alphas/(4*pi)*Cf*( 2*B1dd0(i) + 1d0 )

    seR    = - alpha/(4*pi)*(                               &
               + gpd2*( 2*B1ddZ(i) + 1d0 )                  &
               + md2(i)/(4*st2*mw2)*( B1ddZ(i) + B1ddH(i) ) &
               + md2(i)/(2*st2*mw2)*B1duW(i)                &
               )
    seRqed = - alpha/(4*pi)*Qd2*( 2*B1dd0(i) + 1d0 )
    seRqcd = - alphas/(4*pi)*Cf*( 2*B1dd0(i) + 1d0 )

    if (regd(i).eq.1) then ! massless fermion

      m2dseL    = c0d0
      m2dseLqed = c0d0
      m2dseLqcd = c0d0
      m2dseR    = c0d0
      m2dseRqed = c0d0
      m2dseRqcd = c0d0
      m2dseS    = c0d0
      m2dseSqed = c0d0
      m2dseSqcd = c0d0

    elseif (regd(i).eq.2) then ! light fermion

      if (reg_soft.eq.2) then
        logIR = log(lambda**2/xd2(i)) ! lambda**2 << xd2(i)
      else
        logIR = Dir + log(muIR2/xd2(i))
      endif
      m2dseL    = c0d0
      m2dseLqed = - alpha/(4*pi)*Qd2*( logIR + 3d0 )
      m2dseLqcd = - alphas/(4*pi)*Cf*( logIR + 3d0 )
      m2dseR    = c0d0
      m2dseRqed = - alpha/(4*pi)*Qd2*( logIR + 3d0 )
      m2dseRqcd = - alphas/(4*pi)*Cf*( logIR + 3d0 )
      m2dseS    = c0d0
      m2dseSqed = + alpha/(2*pi)*Qd2*( logIR + 2d0 )
      m2dseSqcd = + alphas/(2*pi)*Cf*( logIR + 2d0 )

    else ! massive fermion

       B0ddZ(i) =  B0 (rmd2(i),md2(i), mz2,Duv,muUV2)
      DB0ddZ(i) = DB0 (rmd2(i),md2(i), mz2)
      DB1ddZ(i) = DB1 (rmd2(i),md2(i), mz2)

       B0duW(i) =  B0 (rmd2(i),mu2(i), mw2,Duv,muUV2)
      DB0duW(i) = DB0 (rmd2(i),mu2(i), mw2)
      DB1duW(i) = DB1 (rmd2(i),mu2(i), mw2)

       B0ddH(i) =  B0 (rmd2(i),md2(i), mh2,Duv,muUV2)
      DB0ddH(i) = DB0 (rmd2(i),md2(i), mh2)
      DB1ddH(i) = DB1 (rmd2(i),md2(i), mh2)

       B0dd0(i) =  B0 (rmd2(i),md2(i),c0d0,Duv,muUV2)
      ratioIm = 1d0 - abs(rxd2(i))/abs(xd2(i))
      if (abs(ratioIm).le.zerocut) then ! for a vanishing width, IR poles show up
        if (reg_soft.eq.1) then
          DB0dd0(i) = - 1/2d0/xd2(i)*( Dir + log(muIR2/xd2(i)) + 2d0 )
          DB1dd0(i) = + 1/2d0/xd2(i)*( Dir + log(muIR2/xd2(i)) + 3d0 )
        else
          DB0dd0(i) = - 1/2d0/xd2(i)*( log(lambda**2/xd2(i)) + 2d0 ) ! lambda**2 << xd2(i)
          DB1dd0(i) = + 1/2d0/xd2(i)*( log(lambda**2/xd2(i)) + 3d0 ) ! lambda**2 << xd2(i)
        endif
      else
        DB0dd0(i) = DB0 (rmd2(i),md2(i),c0d0)
        DB1dd0(i) = DB1 (rmd2(i),md2(i),c0d0)
      endif

      seS    = - alpha/(4*pi)*(                               &
                 + 2*gpd*gmd*( 2*B0ddZ(i) - 1d0 )             &
                 + md2(i)/(4*st2*mw2)*( B0ddZ(i) - B0ddH(i) ) &
                 + mu2(i)/(2*st2*mw2)*B0duW(i)                &
                 )
      seSqed = - alpha/(2*pi)*Qd2*( 2*B0dd0(i) - 1d0 )
      seSqcd = - alphas/(2*pi)*Cf*( 2*B0dd0(i) - 1d0 )

      m2dseL    = - alpha/(4*pi)*rmd2(i)*(                           &
                    + 2*gmd2*DB1ddZ(i)                             &
                    + md2(i)/(4*st2*mw2)*( DB1ddZ(i) + DB1ddH(i) ) &
                    + 1/(2*st2)*( 2d0 + mu2(i)/mw2 )*DB1duW(i)     &
                    )
      m2dseLqed = - alpha/(2*pi)*Qd2*rxd2(i)*DB1dd0(i)
      m2dseLqcd = - alphas/(2*pi)*Cf*rxd2(i)*DB1dd0(i)

      m2dseR    = - alpha/(4*pi)*rmd2(i)*(                           &
                    + 2*gpd2*DB1ddZ(i)                             &
                    + md2(i)/(4*st2*mw2)*( DB1ddZ(i) + DB1ddH(i) ) &
                    + md2(i)/(2*st2*mw2)*DB1duW(i)                 &
                    )
      m2dseRqed = - alpha/(2*pi)*Qd2*rxd2(i)*DB1dd0(i)
      m2dseRqcd = - alphas/(2*pi)*Cf*rxd2(i)*DB1dd0(i)

      m2dseS    = - alpha/(4*pi)*rmd2(i)*(                         &
                    + 4*gpd*gmd*DB0ddZ(i)                          &
                    + md2(i)/(4*st2*mw2)*( DB0ddZ(i) - DB0ddH(i) ) &
                    + mu2(i)/(2*st2*mw2)*DB0duW(i)                 &
                    )
      m2dseSqed = - alpha/pi*Qd2*rxd2(i)*DB0dd0(i)
      m2dseSqcd = - alphas/pi*Cf*rxd2(i)*DB0dd0(i)

    endif

    dZdL(i) = c0d0
    if (loopWEAK) &
    dZdL(i) = dZdL(i) - seL - m2dseL - m2dseR - 2*m2dseS 
    if (loopQED) &
    dZdL(i) = dZdL(i) - seLqed - m2dseLqed - m2dseRqed - 2*m2dseSqed
    dZdLqcd(i) = - seLqcd - m2dseLqcd - m2dseRqcd - 2*m2dseSqcd

    dZdR(i) = c0d0
    if (loopWEAK) &
    dZdR(i) = dZdR(i) - seR - m2dseL - m2dseR - 2*m2dseS
    if (loopQED) &
    dZdR(i) = dZdR(i) - seRqed - m2dseLqed - m2dseRqed - 2*m2dseSqed
    dZdRqcd(i) = - seRqcd - m2dseLqcd - m2dseRqcd - 2*m2dseSqcd

    if (reguScheme.eq.1) then
      if (loopWEAK) &
      dZdL(i)    = dZdL(i)    - alpha/(4*pi)*( gmu2 + 1/(2*st2) )
      if (loopQED) &
      dZdL(i)    = dZdL(i)    - alpha/(4*pi)*Qu2
      dZdLqcd(i) = dZdLqcd(i) - alphas/(4*pi)*Cf
      if (loopWEAK) &
      dZdR(i)    = dZdR(i)    - alpha/(4*pi)*gpu2
      if (loopQED) &
      dZdR(i)    = dZdR(i)    - alpha/(4*pi)*Qu2
      dZdRqcd(i) = dZdRqcd(i) - alphas/(4*pi)*Cf
    endif

    if (complex_mass_scheme.eq.0) then
      dZdL(i) = real(dZdL(i),kind=dp)*c1d0
      dZdR(i) = real(dZdR(i),kind=dp)*c1d0
      dZdLqcd(i) = real(dZdLqcd(i),kind=dp)*c1d0
      dZdRqcd(i) = real(dZdRqcd(i),kind=dp)*c1d0
    endif

    dZdLc   (i) = dZdL   (i)
    dZdLcqcd(i) = dZdLqcd(i)

    dZdRc   (i) = dZdR   (i)
    dZdRcqcd(i) = dZdRqcd(i)

    if (regd(i).le.2) then ! massless and light fermion

      dmd(i)    = c0d0
      dmdqcd(i) = c0d0

    else ! massive fermion

      dmd(i) = c0d0
      if (loopWEAK) &
      dmd(i) = dmd(i) &
               + md1(i)/2d0*(                      &
                 + seL + seR + 2*seS               &
                 + ( md2(i) - rmd2(i) ) / rmd2(i)* &
                   ( m2dseL + m2dseR + 2*m2dseS )  &
                 )
      if (loopQED) &
      dmd(i) = dmd(i) &
               + md1(i)/2d0*(                              &
                 + seLqed + seRqed + 2*seSqed              &
                 + ( md2(i) - rmd2(i) ) / rmd2(i)*         &
                   ( + m2dseLqed + m2dseRqed + 2*m2dseSqed &
                     - alpha/pi*Qd2 )                      &
                 )
      dmdqcd(i) = - alphas/(4*pi)*md1(i)*Cf*               &
                    ( 3*(Duv + log(muUV2/md2(i))) + 4d0 )

      if (reguScheme.eq.1) then
        if (loopWEAK) &
        dmd(i)    = dmd(i) + alpha/(8*pi)*md1(i)*(          &
                      + gpd2 + gmd2 + 4*gpd*gmd + 1/(2*st2) &
                      )
        if (loopQED) &
        dmd(i)    = dmd(i)    + alpha/(4*pi)*md1(i)*3*Qd2
        dmdqcd(i) = dmdqcd(i) + alphas/(4*pi)*md1(i)*3*Cf
      endif

    endif

    if (complex_mass_scheme.eq.0) then
      dmd(i) = real(dmd(i),kind=dp)*c1d0
      dmdqcd(i) = real(dmdqcd(i),kind=dp)*c1d0
    endif

  enddo

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (write_counterterms) then
    call openOutput
    write(nx,*) 
    write(nx,*) '+++++++++++++++ COUNTERTERMS +++++++++++++++'
    write(nx,*) 
    write(nx,'(a,2e19.11)') 'dt   =',dt  
    write(nx,'(a,2e19.11)') 'dZh  =',dZh 
    write(nx,'(a,2e19.11)') 'dmh2 =',dmh2
    write(nx,'(a,2e19.11)') 'dZg  =',dZg 
    write(nx,'(a,2e19.11)') 'dZgs =',dZgs
    write(nx,'(a,2e19.11)') 'dZaa =',dZaa
    write(nx,'(a,2e19.11)') 'dZza =',dZza
    write(nx,'(a,2e19.11)') 'dZaz =',dZaz
    write(nx,'(a,2e19.11)') 'dZe  =',dZe 
    write(nx,'(a,2e19.11)') 'dZzz =',dZzz
    write(nx,'(a,2e19.11)') 'dmz2 =',dmz2
    write(nx,'(a,2e19.11)') 'dZw  =',dZw 
    write(nx,'(a,2e19.11)') 'dmw2 =',dmw2
    write(nx,'(a,2e19.11)') 'dst  =',dst 
    write(nx,'(a,2e19.11)') 'dct  =',dct 
    write(nx,'(a,2e19.11)') 'dgpn =',dgpn
    write(nx,'(a,2e19.11)') 'dgmn =',dgmn
    write(nx,'(a,2e19.11)') 'dgpu =',dgpu
    write(nx,'(a,2e19.11)') 'dgmu =',dgmu
    write(nx,'(a,2e19.11)') 'dgpl =',dgpl
    write(nx,'(a,2e19.11)') 'dgml =',dgml
    write(nx,'(a,2e19.11)') 'dgpd =',dgpd
    write(nx,'(a,2e19.11)') 'dgmd =',dgmd
    do i = 1,3
    write(nx,'(a,i1,a,2e19.11)') 'dZnL(',i,') =',dZnL(i)
    enddo
    do i = 1,3
    write(nx,'(a,i1,a,2e19.11)') 'dZlL(',i,') =',dZlL(i)
    write(nx,'(a,i1,a,2e19.11)') 'dZlR(',i,') =',dZlR(i)
    write(nx,'(a,i1,a,2e19.11)') 'dml (',i,') =',dml (i)
    enddo
    do i = 1,3
    write(nx,'(a,i1,a,2e19.11)') 'dZuL   (',i,') =',dZuL   (i)
    write(nx,'(a,i1,a,2e19.11)') 'dZuLqcd(',i,') =',dZuLqcd(i)
    write(nx,'(a,i1,a,2e19.11)') 'dZuR   (',i,') =',dZuR   (i)
    write(nx,'(a,i1,a,2e19.11)') 'dZuRqcd(',i,') =',dZuRqcd(i)
    write(nx,'(a,i1,a,2e19.11)') 'dmu    (',i,') =',dmu    (i)
    write(nx,'(a,i1,a,2e19.11)') 'dmuqcd (',i,') =',dmuqcd (i)
    enddo
    do i = 1,3
    write(nx,'(a,i1,a,2e19.11)') 'dZdL   (',i,') =',dZdL   (i)
    write(nx,'(a,i1,a,2e19.11)') 'dZdLqcd(',i,') =',dZdLqcd(i)
    write(nx,'(a,i1,a,2e19.11)') 'dZdR   (',i,') =',dZdR   (i)
    write(nx,'(a,i1,a,2e19.11)') 'dZdRqcd(',i,') =',dZdRqcd(i)
    write(nx,'(a,i1,a,2e19.11)') 'dmd    (',i,') =',dmd    (i)
    write(nx,'(a,i1,a,2e19.11)') 'dmdqcd (',i,') =',dmdqcd (i)
    enddo
    write(nx,*) 
    write(nx,*) 
!    stop
  endif

  end subroutine counterterms_tables

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module tables_rcl

!#####################################################################
