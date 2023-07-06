!
!        Copyright (C) 2000-2023 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM,AF
! 
! This file is distributed under the terms of the GNU 
! General Public License. You can redistribute it and/or 
! modify it under the terms of the GNU General Public 
! License as published by the Free Software Foundation; 
! either version 2, or (at your option) any later version.
!
! This program is distributed in the hope that it will 
! be useful, but WITHOUT ANY WARRANTY; without even the 
! implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License 
! for more details.
!
! You should have received a copy of the GNU General Public 
! License along with this program; if not, write to the Free 
! Software Foundation, Inc., 59 Temple Place - Suite 330,Boston, 
! MA 02111-1307, USA or visit http://www.gnu.org/copyleft/gpl.txt.
!
subroutine SERIAL_HERMITIAN_diagonalization(n,M,E)
 !
 use magma
 use pars,           ONLY:SP
 use drivers,        ONLY:l_nl_optics,l_real_time
 use linear_algebra, ONLY:LINEAR_ALGEBRA_error,&
&                         LINEAR_ALGEBRA_WS_reset,LALGEBRA_WS
 use timing_m,       ONLY:timing
#include<memory.h>
 !
 integer     :: n
 complex(SP) :: M(n,n)
 real(SP)    :: E(n)
 !
 character(64)     :: subname="SERIAL_HERMITIAN_diagonalization"
 integer           :: lwork, lrwork, liwork
 type(LALGEBRA_WS) :: WS
 character (len = 70) :: str_n, str_lwork, str_lrwork, str_liwork
 !
 ! MAGMA
 call magmaf_init()
 call magmaf_print_environment()
 call magmaf_finalize()
 ! FIN MAGMA
 !
 lwork=-1
 lrwork=-1
 liwork=-1
 allocate(WS%v_cmplx(1))
 allocate(WS%v_real(1))
 allocate(WS%v_int(1))
 !
 if(.not.(l_nl_optics.or.l_real_time)) then
   call timing('SERIAL_HERMITIAN_diagonalization',OPR='start')
 endif
 !
 ! Query lwork, lrwork, liwork
 !
 print *,"Query lwork, lrwork, liwork"
#if defined _DOUBLE
 call ZHEEV('V','U',n,M,size(M,1),E,WS%v_cmplx,lwork,WS%v_real,WS%i_fail)
#else
 !call CHEEV('V','U',n,M,size(M,1),E,WS%v_cmplx,lwork,WS%v_real,WS%i_fail)
 call magmaf_cheevd('v','u',n,M,size(M,1),E,WS%v_cmplx,lwork,WS%v_real,lrwork,WS%v_int,liwork,WS%i_fail)
#endif
 !
 ! Allocate v_cmplx, v_real, v_int
 !
 lwork=nint(real(WS%v_cmplx(1)))
 lrwork=nint(WS%v_real(1))
 liwork=WS%v_int(1)
 deallocate(WS%v_cmplx)
 deallocate(WS%v_real)
 deallocate(WS%v_int)
 !
 if(.not.(l_nl_optics.or.l_real_time)) then
   YAMBO_ALLOC(WS%v_cmplx,(lwork))
   YAMBO_ALLOC(WS%v_real,(lrwork))
   YAMBO_ALLOC(WS%v_int,(liwork))
 else
   allocate(WS%v_cmplx(lwork))
   allocate(WS%v_real(lrwork))
   allocate(WS%v_int(liwork))
 endif
 !
 ! Print
 !
 write (str_n, "(A3,I20)") "n: ", n
 write (str_lwork, "(A7,I20)") "lwork: ", lwork
 write (str_lrwork, "(A8,I20)") "lrwork: ", lrwork
 write (str_liwork, "(A8,I20)") "liwork: ", liwork
 print *,str_n
 print *,str_lwork
 print *,str_lrwork
 print *,str_liwork
 !
 ! Compute eigenvalues
 !
#if defined _DOUBLE
 call ZHEEV('V','U',n,M,size(M,1),E,WS%v_cmplx,lwork,WS%v_real,WS%i_fail)
 if(WS%i_fail.ne.0) call LINEAR_ALGEBRA_error(subname,'performing ZHEEV')
#else
 call magmaf_cheevd('v','u',n,M,size(M,1),E,WS%v_cmplx,lwork,WS%v_real,lrwork,WS%v_int,liwork,WS%i_fail)
 !call CHEEV('V','U',n,M,size(M,1),E,WS%v_cmplx,lwork,WS%v_real,WS%i_fail)
 if(WS%i_fail.ne.0) call LINEAR_ALGEBRA_error(subname,'performing CHEEV')
#endif
 !
 call LINEAR_ALGEBRA_WS_reset(WS)
 !
 if(.not.(l_nl_optics.or.l_real_time)) then
   call timing('SERIAL_HERMITIAN_diagonalization',OPR='stop')
 endif
 print *,"Fin diagonalización"
 !
end subroutine
