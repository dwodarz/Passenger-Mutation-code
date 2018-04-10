
module vars
implicit none

   integer :: nsites
   parameter (nsites = 2)

   type spot
      integer :: cell                                  
      integer :: genome(nsites)
	  integer :: defect(nsites)
   end type spot
   
end module vars





program space
USE IFCORE
use vars
implicit none


integer :: n
parameter (n=50)
 
type(spot) sp(n,n)
integer   :: dist, spread_rad
integer   :: total, wt, coop, def, complete, compl_single, compl_double, mutants(nsites)
integer   :: min, max, imin, imax, jmin, jmax, count, idummy, site_count, site_count2, site_count3
integer   :: ch(nsites), nzombies
real      :: ran0, xran, xr1, xr2
real      :: divA, divD, div_coop, div_prob, death_prob, div_wt, death_wt, death_coop, div_alone, death_alone, div_av
REAL      :: div_alone2, div_wt2      
real(8)   :: s1, s2, mu, mult                       
integer  :: wt2, mutants2(nsites)    
real(8)   :: xid, ss, xid2   
integer   :: i, j, ii, i_hits, idx, div_count, coop_count, nocoop_count, onecoop, nocoop, cp(nsites), iruns
logical   :: found, check
character*30 arg

   
call get_command_argument(1, arg)

!write(*,*) arg
!read(*,*)   



!open(66, file='rr.dat', form='formatted' )
open(10, file=arg, form='formatted' )
!open(20, file='av.dat', form='formatted' )
!open(30, file='rch.dat', form='formatted' )
!open(40, file='frac.dat', form='formatted' )
 

 
!write(*,*) ran0()
!write(*,*) 'graphical output (y/n)?'
!read(*,*) answer

 

ss = 0.
DO iruns=1,100000
!write(*,*) iruns


dist = 1
spread_rad = 1


mult = 0*500000

mu = 1E-7  !1E-4
 


s1 = 0.99
s2 = 0.99
 
div_wt = 0.15
div_wt2 = s2*div_wt
div_coop = 0.5   
 
death_wt = 0.01    
death_coop = 0.01
 
 
div_alone  = s1*div_wt   
div_alone2 = s2*div_wt     
death_alone = death_wt
   

  
do i = 1, n       
   do j = 1, n
      sp(i,j)%cell = 0
      sp(i,j)%genome = 0
	  sp(i,j)%defect = 0
   enddo
enddo 



min=1
max=n

do ii = 1,1500  
   xr1 = ran0()*0.9999999
   i = int(xr1*(max+1-min))+min 
   xr2 = ran0()*0.9999999             
   j = int(xr2*(max+1-min))+min  
   
   sp(i,j)%cell = 1
enddo 

 
 





xid = 0. 
xid2 = 0
count = 0      
check = .false.
do while( .true. )   
   xid = xid + 1
   count = count + 1
    
    !write(*,*) xid
    
!   if ( count .eq. 1 ) then
!      count=0
!       call sumcells( sp, n, total, wt, coop, def, complete, compl_single, compl_double, ch, onecoop, nocoop, cp)  
    
!	  if (total == 0) exit
       
      !write(10,'(1x,f20.10,3x,8(i10,3x))') xid, wt, coop, def, complete, total, compl_single, compl_double, compl_single+compl_double
	  !if (div_count > 0 ) write(20,'(1x,f20.10,3x,3(f20.10,3x))') xid, div_av/real(div_count), real(coop_count)/real(div_count), real(nocoop_count)/real(div_count)
      !write(30,'(1x,f20.10,3x,5(i10,3x))') xid, ch(1), ch(2), ch(3), ch(4), ch(5)
	  !write(40,'(1x,f20.10,3x,f20.10)') xid, real(complete)/real(total)
!   endif

   div_count = 0
   div_av = 0
   coop_count = 0
   nocoop_count = 0


   ! update the patches
   ! ran between max and min: int(max-min+1)*Rand + min
  
   !if ( count .eq. 100 ) then
   !   count=0
   !   call sumcells( sp, n, total, wt, coop, def, complete, compl_single, compl_double) 
       
        
	!  if (total == 0) exit
       
     !  write(66,'(1x,f40.20,3x,7(i10,3x))') xid, wt, coop, def, complete, total, compl_single, compl_double
   !endif
   
   
    call sumcells( sp, n, total, wt, coop, def, complete, compl_single, compl_double, mutants, mutants2, wt2)     
 
   if (total==0) exit
   !if (wt==0) exit
  if ( check ) then    ! real(complete)/real(total) > 0.9  OR check
     !ss = ss + xid
	 !write(*,'(i10,3x,f20.10)') iruns, ss/real(iruns)
     write(10,'(1x,i20)') 1
	 exit
  endif
  
  if ( xid ==1500000)  then
       write(10,'(1x,i20)') 0
	   exit
   endif
  
   !if (coop == 0) then
   !  write(*,*) xid, wt, coop, def
   !  read(*,*)
   !endif
      
    imin = 1
   imax = n
   jmin = 1
   jmax = n
      
    i_hits = 0
   do while(i_hits < total)
      xr1 = ran0()*0.9999999
      i = int(xr1*(imax+1-imin))+imin 
      xr2 = ran0()*0.9999999             
      j = int(xr2*(jmax+1-jmin))+jmin  
   
      !if (sp(i,j)%cell .ne. 1 ) cycle
	  if ( sp(i,j)%cell == 0 ) cycle
      i_hits = i_hits + 1
   
   
      if ( i<0 .or. j<0 ) then
         write(*,*) 'error'
         read(*,*) 
      endif
       
      site_count  = 0
      site_count2 = 0
      do idx = 1,nsites
         if ( sp(i,j)%genome(idx) == 1 .or. sp(i,j)%genome(idx) == 2 ) site_count = site_count + 1
         if ( sp(i,j)%genome(idx) == 2 ) site_count2 = site_count2 + 1
      enddo  
                  
      
      if ( site_count == 0 ) then
         if ( sp(i,j)%cell  == 1  ) then                     !######### START  ###########
             div_prob = div_wt
	     else if ( sp(i,j)%cell  == 2 ) then
		    div_prob = div_wt2
	     endif                                                    !######### END   ########### 
		 
         death_prob = death_wt
      else if ( site_count == nsites ) then
        div_prob = div_coop
        death_prob = death_coop
		div_prob = div_prob + 0*real(nsites)/70.
      else
       if ( sp(i,j)%cell  == 1  ) then                    !######### START  ###########
           div_prob = div_alone
	    else if ( sp(i,j)%cell  == 2 ) then
            div_prob = div_alone2
        endif	 
        death_prob = death_alone
        div_prob = div_prob + 0*site_count/70.          !######### END   ###########  
     endif
 
       
      if ( sp(i,j)%cell  == 1 .or. sp(i,j)%cell  == 2  ) then    ! ####### LINE  #########
		  xran = ran0()
		  if (xran < div_prob) then 
			 !write(*,*) div_prob
			 !read(*,*) 
			 !call divide2( sp, i, j, n, dist, spread_rad, site_count ) 
			 call dividemass( sp, i, j, n, dist, spread_rad, site_count, check, xid2, xid, mu, mult )
			 !call dividespace( sp, i, j, n, dist, spread_rad, site_count, check, xid2, xid )
		  else if (xran < div_prob+death_prob) then
			 call death( sp, i, j, n )
			 xid2 = xid2 + 1
		  endif
	  endif
	  
	 
      	
		 
   
    enddo  ! spots  
    
       
 
enddo   ! time



ENDDO


end program 









subroutine death (sp, i, j, n)
use vars
implicit none

integer    :: i, j, n
type(spot) :: sp(n,n)



sp(i,j)%cell = 0
sp(i,j)%genome = -1 !0
sp(i,j)%defect = -1 !0

return
end subroutine death







subroutine zombie (sp, i, j, n)
use vars
implicit none

integer    :: i, j, n
type(spot) :: sp(n,n)



sp(i,j)%cell = 2
sp(i,j)%genome = -1 !0
sp(i,j)%defect = -1 !0

return
end subroutine zombie




  
  
  
  
 subroutine divide( sp, i, j, n, dist, spread_rad, site_count, check, xid2 ) 
use vars
implicit none

integer    :: i, j, n, dist, spread_rad, site_count
type(spot) :: sp(n,n)

real    :: ran0, xran, xran2, xr1, xr2
real(8) :: xid2, mu, mu2
integer :: idx, min, max, sum_sites, mutcount
integer :: ipos, jpos, idist, jdist, isign, jsign
logical :: found, check
real(8) :: ns, ng, q


ns = 1
ng = 3000
q=0.9

mu = 1E-6  !1E-6



xran = ran0()
if ( xran < mu*ng*q ) then
    call zombie( sp, i, j, n )
	return
endif



min = 0
max = spread_rad
 
 

xr1 = ran0()*0.9999999
idist = int(xr1*(max+1-min))+min 
xr2 = ran0()*0.9999999             
jdist = int(xr2*(max+1-min))+min   
 
 
if ( ran0() < 0.5 ) then
   isign = -1
else 
   isign = 1
endif
if ( ran0() < 0.5 ) then
   jsign = -1
else 
   jsign = 1
endif


idist = isign*idist
jdist = jsign*jdist
 
ipos = i + idist
jpos = j + jdist


if (ipos < 1) then
   ipos = n+ipos
else if (ipos > n) then
   ipos = ipos-n
endif

if (jpos < 1) then
   jpos = n+jpos
else if (jpos > n) then
   jpos = jpos-n
endif



if (ipos <= 0 .or. jpos <= 0 .or. ipos>n .or. jpos>n ) then
   write(*,*) 'error: ', ipos, jpos
   read(*,*)
endif


 
if (  sp(ipos,jpos)%cell .ne. 0  )  return   

 
 
 
!if ( sp(i,j)%genome(1)==0 .and. sp(i,j)%genome(2)==0 ) then
!     mu = 1E-5
!else
!    mu = 1E-4
!endif
 

! prevent second mutation
!mutated1 = .false.
!do iidx = 1,nsites
!    if ( sp(i,j)%genome(iidx) == 1 ) mutated1 = .true.
!enddo 
 
mutcount = 0
do idx = 1,nsites
   if ( sp(i,j)%genome(idx) == 0 ) then
    
       ! prevent second mutation	   
	   !mutated2 = .false. 
       !do iidx = 1,nsites 
       !    if ( sp(ipos,jpos)%genome(iidx) == 1 ) mutated2 = .true.
       !enddo
	   
	   xran = ran0()
	   if (  (xran < ns*mu)   ) then   ! .and. (.not. mutated1) .and. (.not. mutated2)
	        
		   sp(ipos,jpos)%genome(idx) = 1
		   mutcount = mutcount + 1
	   else
          sp(ipos,jpos)%genome(idx) = 0
       endif
    endif
    
    if ( sp(i,j)%genome(idx) == 1 ) then
	   xran = ran0()
	   if (xran < 0*mu ) then     !0* backmutations
	      sp(ipos,jpos)%genome(idx) = 0
	   else	   
         sp(ipos,jpos)%genome(idx) = 1 
       endif		  
    endif
    
    if ( sp(i,j)%genome(idx) == 2 ) then
       write(*,*) 'crap'
       read(*,*)
       sp(ipos,jpos)%genome(idx) = sp(i,j)%genome(idx)
    endif
        
enddo

!if (mutcount == 2) then
!    write(*,*) 'double'
!	read(*,*)
!endif

 
sp(ipos,jpos)%cell = 1
  
sum_sites = 0
do idx = 1,nsites
    sum_sites = sum_sites + sp(ipos,jpos)%genome(idx)
enddo
if ( sum_sites == nsites ) check = .true.

!if (check) sp(ipos,jpos)%genome(2)=0

xid2 = xid2 + 1

 !do idx = 1,nsites
 !   if ( sp(i,j)%genome(idx) == 1 ) then
 !      write(*,*) sp(ipos,jpos)%genome(idx)
 !      read(*,*)  
 !   endif
 !enddo
 
 
return
end subroutine divide  
  
  
  
  
  
  
  
! #############  THIS ROUTINE - change mutation rates here  ###########################
 subroutine dividemass( sp, i, j, n, dist, spread_rad, site_count, check, xid2, xid, mu, mult ) 
use vars
implicit none

integer    :: i, j, n, dist, spread_rad, site_count
type(spot) :: sp(n,n)

real    :: ran0, xran, xran2, xr1, xr2, xran3, xran4
real(8) :: xid2, xid, mu, mu2
integer :: idx, iidx, min, max, sum_sites, mutcount 
integer :: ipos, jpos, idist, jdist, isign, jsign
logical :: found, check, mutated1, mutated2
real(8) :: ns, ng, q, mult


if ( xid > 0 ) then
   ns = 1
else
   ns = 0
endif

ng = 0*100
q=1

!mult = 0*5000

!mu = 1E-5  !1E-4
 


!if ( site_count > 0 ) then 
!   found = .false.
!   if ( site_count == nsites ) then
!      found = .true.
!   else
!      call search_coop(sp, i, j, n, dist, found )
!   endif
!   
!   
!   if ( .not. found ) return
!   
!endif







min = 1
max = n
 
 
xr1 = ran0()*0.9999999
ipos = int(xr1*(max+1-min))+min 
xr2 = ran0()*0.9999999             
jpos = int(xr2*(max+1-min))+min   
 
  


if (ipos <= 0 .or. jpos <= 0 .or. ipos>n .or. jpos>n ) then
   write(*,*) 'error: ', ipos, jpos
   read(*,*)
endif


 
if (  sp(ipos,jpos)%cell .ne. 0  )  return  !return  ! i.e. if patch is full

xran = ran0()
if ( xran < mu*ng*q ) then
	return
endif
 
 

 

 
 
if ( sp(i,j)%cell == 1 ) then     
	 
	 xran3 = ran0()
     if ( xran3 < mult*mu ) then
	     sp(ipos,jpos)%cell = 2
         do idx = 1,nsites
             xran4 = ran0()  
             if ( xran4 < mu ) then
			     sp(ipos,jpos)%genome(idx) = 1
		     else
		         sp(ipos,jpos)%genome(idx) = sp(i,j)%genome(idx)
		     endif
		 enddo
	  else
         sp(ipos,jpos)%cell = 1
         do idx = 1,nsites
             xran4 = ran0()  
             if ( xran4 < mu ) then
			     sp(ipos,jpos)%genome(idx) = 1
		     else
		         sp(ipos,jpos)%genome(idx) = sp(i,j)%genome(idx)
		     endif
		 enddo	 
      endif
	    
endif 
 


if ( sp(i,j)%cell == 2 ) then
    sp(ipos,jpos)%cell = 2
    do idx = 1,nsites
        xran4 = ran0()  
        if ( xran4 < mu ) then
			sp(ipos,jpos)%genome(idx) = 1
		else
		    sp(ipos,jpos)%genome(idx) = sp(i,j)%genome(idx)
		endif
	enddo 
endif	 
	 
	 
	 
	 
!####################################################	 
 
 

  
sum_sites = 0
do idx = 1,nsites
    sum_sites = sum_sites + sp(ipos,jpos)%genome(idx)
enddo
if ( sum_sites == nsites ) check = .true.



xid2 = xid2 + 1

 
 
 
return
end subroutine dividemass


  
  
 






 subroutine dividespace( sp, i, j, n, dist, spread_rad, site_count, check, xid2, xid ) 
use vars
implicit none

integer    :: i, j, n, dist, spread_rad, site_count
type(spot) :: sp(n,n)

real    :: ran0, xran, xran2, xr1, xr2, mu, mu2, xran3, xran4
real(8) :: xid2, xid
integer :: idx, iidx, min, max, sum_sites, mutcount 
integer :: ipos, jpos, idist, jdist, isign, jsign
logical :: found, check, mutated1, mutated2
real(8) :: ns, ng, q, mult


if ( xid > 0 ) then
   ns = 1
else
   ns = 0
endif

ng = 0*100
q=1

mult = 0*500

mu = 1E-4  
 


   
 
  

spread_rad = 1
dist = 1


min = 0
max = spread_rad
 
 

xr1 = ran0()*0.9999999
idist = int(xr1*(max+1-min))+min 
xr2 = ran0()*0.9999999             
jdist = int(xr2*(max+1-min))+min   
 
 
if ( ran0() < 0.5 ) then
   isign = -1
else 
   isign = 1
endif
if ( ran0() < 0.5 ) then
   jsign = -1
else 
   jsign = 1
endif


idist = isign*idist
jdist = jsign*jdist
 
ipos = i + idist
jpos = j + jdist


if (ipos < 1) then
   ipos = n+ipos
else if (ipos > n) then
   ipos = ipos-n
endif

if (jpos < 1) then
   jpos = n+jpos
else if (jpos > n) then
   jpos = jpos-n
endif






if (ipos <= 0 .or. jpos <= 0 .or. ipos>n .or. jpos>n ) then
   write(*,*) 'error: ', ipos, jpos
   read(*,*)
endif


 
if (  sp(ipos,jpos)%cell .ne. 0  )  return  

xran = ran0()
if ( xran < mu*ng*q ) then
	return
endif
 
 


 
 
if ( sp(i,j)%cell == 1 ) then     
	 
	 xran3 = ran0()
     if ( xran3 < mult*mu ) then
	     sp(ipos,jpos)%cell = 2
         do idx = 1,nsites
             xran4 = ran0()  
             if ( xran4 < mu ) then
			     sp(ipos,jpos)%genome(idx) = 1
		     else
		         sp(ipos,jpos)%genome(idx) = sp(i,j)%genome(idx)
		     endif
		 enddo
	  else
         sp(ipos,jpos)%cell = 1
         do idx = 1,nsites
             xran4 = ran0()  
             if ( xran4 < mu ) then
			     sp(ipos,jpos)%genome(idx) = 1
		     else
		         sp(ipos,jpos)%genome(idx) = sp(i,j)%genome(idx)
		     endif
		 enddo	 
      endif
	    
endif 
 


if ( sp(i,j)%cell == 2 ) then
    sp(ipos,jpos)%cell = 2
    do idx = 1,nsites
        xran4 = ran0()  
        if ( xran4 < mu ) then
			sp(ipos,jpos)%genome(idx) = 1
		else
		    sp(ipos,jpos)%genome(idx) = sp(i,j)%genome(idx)
		endif
	enddo 
endif	 
	 
	 
	 
	 
!####################################################	 
 
 

  
sum_sites = 0
do idx = 1,nsites
    sum_sites = sum_sites + sp(ipos,jpos)%genome(idx)
enddo
if ( sum_sites == nsites ) check = .true.



xid2 = xid2 + 1


 
 
return
end subroutine dividespace




 
   




 
 
 
 
 
 
 subroutine search_coop(sp, i, j, n, dist, res )
use vars
implicit none

integer    :: i, j, n, dist 
logical    :: res
type(spot) :: sp(n,n)

integer :: ii, jj, idx
integer :: imin, imax, jmin, jmax
integer :: imin2, imax2, jmin2, jmax2
logical :: i_1, i_n, j_1, j_n
logical :: found(nsites) 


imin = i-dist
imax = i+dist
jmin = j-dist
jmax = j+dist



i_1=.false.; i_n=.false.; j_1=.false.; j_n=.false.;
imin2=-1; jmin2=-1; imax2=-1; jmax2=-1;

if (imin < 1) then
   i_1 = .true.
   imin  = 1
   imin2 = n - abs(dist-i)
endif
if (jmin < 1) then
   j_1 = .true.
   jmin  = 1
   jmin2 = n - abs(dist-j)
endif


if (imax > n) then
   i_n = .true.
   imax  = n
   imax2 = dist - (n-i)
endif
if (jmax > n) then
   j_n = .true.
   jmax  = n
   jmax2 = dist - (n-j)
endif


if (i_1) then
   if (imin2<1 .or. imin2>n ) then
      write(*,*) 'imin: ', imin2
      read(*,*)
   endif
endif
if (j_1) then
   if (jmin2<1 .or. jmin2>n ) then
      write(*,*) 'jmin: ', jmin2
      read(*,*)
   endif
endif
if (i_n) then
   if (imax2<1 .or. imax2>n ) then
      write(*,*) 'imax: ', imax2
      read(*,*)
   endif
endif
if (j_n) then
   if (jmax2<1 .or. jmax2>n ) then
      write(*,*) 'jmax: ', jmax2
      read(*,*)
   endif
endif 

 

found = .false.

do idx = 1,nsites
   if ( sp(i,j)%genome(idx) == 1 .or. sp(i,j)%genome(idx) == 2 ) then
      found(idx) = .true.
	  cycle
	endif

   do ii = imin, imax
	   do jj = jmin, jmax
		  if ( ii==i .and. jj==j ) cycle
		  if ( sp(ii,jj)%genome(idx) == 1 ) then
			 found(idx) = .true.
			 goto 10
		  endif
	   enddo
	   
	   if (j_1) then
		  do jj = jmin2, n
			 if ( ii==i .and. jj==j ) cycle
			 if ( sp(ii,jj)%genome(idx) == 1 ) then
				found(idx) = .true.
				goto 10
			 endif
		  enddo
	   endif
	   
	   if (j_n) then
		  do jj = 1, jmax2
			 if ( ii==i .and. jj==j ) cycle
			 if ( sp(ii,jj)%genome(idx) == 1 ) then
				found(idx) = .true.
				goto 10
			 endif
		  enddo
	   endif
   enddo
   
   
   
   
	if (i_1) then
	   do ii = imin2, n
		  do jj = jmin, jmax
			 if ( ii==i .and. jj==j ) cycle
			 if ( sp(ii,jj)%genome(idx) == 1 ) then
				found(idx) = .true.
				goto 10
			 endif
		  enddo
		  
		  if (j_1) then
			 do jj = jmin2, n
				if ( ii==i .and. jj==j ) cycle
				if ( sp(ii,jj)%genome(idx) == 1 ) then
				   found(idx) = .true.
				   goto 10
				endif
			 enddo
		  endif
	   
		  if (j_n) then
			 do jj = 1, jmax2
				if ( ii==i .and. jj==j ) cycle
				if ( sp(ii,jj)%genome(idx) == 1 ) then
				   found(idx) = .true.
				   goto 10
				endif
			 enddo
		  endif
			  
	   enddo
	endif


	if (i_n) then
	   do ii =1, imax2
		  do jj = jmin, jmax
			 if ( ii==i .and. jj==j ) cycle
			 if ( sp(ii,jj)%genome(idx) == 1 ) then
				found(idx) = .true.
			    goto 10
			 endif
		  enddo
		  
		  if (j_1) then
			 do jj = jmin2, n
				if ( ii==i .and. jj==j ) cycle
				if ( sp(ii,jj)%genome(idx) == 1 ) then
				   found(idx) = .true.
				   goto 10
				endif
			 enddo
		  endif
	   
		  if (j_n) then
			 do jj = 1, jmax2
				if ( ii==i .and. jj==j ) cycle
				if ( sp(ii,jj)%genome(idx) == 1 ) then
				   found(idx) = .true.
				   goto 10
				endif
			 enddo
		  endif
		  
		  
	   enddo
	endif
   
    10 continue 
   
enddo  !ncells






do idx = 1,nsites
   if ( .not. found(idx) ) then
      res = .false.
      return
   endif
enddo 

res = .true.  
 
 
 return
 end subroutine search_coop






 
 
 
  
 
 
 ! #############  THIS ROUTINE  ###########################
 subroutine sumcells( sp, n, total, wt, coop, def, complete, compl_single, compl_double, mutants, mutants2, wt2) 
use vars
implicit none

integer    :: n
type(spot) :: sp(n,n)
integer    :: total, wt, coop, def, complete, compl_single, compl_double, mutants(nsites), mutants2(nsites), wt2

integer    :: ii, jj, idx, site_count, site_count2, nmuts


total = 0
wt = 0
wt2 = 0
coop = 0
def = 0
complete = 0


compl_single = 0   
compl_double = 0   

mutants = 0
mutants2 = 0
do ii = 1,n
   do jj = 1,n
      if ( sp(ii,jj)%cell == 0 ) cycle
	  
      total = total + 1
      
      site_count  = 0
      site_count2 = 0
	  nmuts = 0
      do idx = 1,nsites
         if ( sp(ii,jj)%genome(idx) == 1 .or. sp(ii,jj)%genome(idx) == 2 ) site_count  = site_count + 1
         if ( sp(ii,jj)%genome(idx) == 2 ) site_count2  = site_count2 + 1
		 
		 nmuts = nmuts + sp(ii,jj)%genome(idx)
      enddo
	  
	  if (nmuts > 0 .and. sp(ii,jj)%cell == 1  ) mutants(nmuts) = mutants(nmuts) + 1
	  if (nmuts > 0 .and. sp(ii,jj)%cell == 2  ) mutants2(nmuts) = mutants2(nmuts) + 1
	  
      
      
      if ( site_count == 0 ) then
         if ( sp(ii,jj)%cell == 1 ) wt = wt + 1 
		 if ( sp(ii,jj)%cell == 2 ) wt2 = wt2 + 1 
      else if ( site_count > 0 ) then
        if ( site_count2 > 0 ) then
           def = def + 1
        else
           coop = coop + 1
        endif
      endif
      if (site_count == nsites) then
         complete = complete + 1
         if (site_count2 == nsites) then
            compl_double = compl_double + 1
         else if ( site_count==nsites .and. site_count2==1) then
            compl_single = compl_single + 1
         endif
      endif
      
   enddo
enddo



return
end subroutine sumcells
 
 
 

 
 











 












real function ran0()
implicit none

integer, save :: idum=0  
real          :: ran3, xran




if (idum ==0) then
   call random_seed() 
   idum = -1 
endif

call random_number(xran)
ran0 = xran

 

if ( ran0<0 .or. ran0>1 ) then
   write(*,*) 'random number problem'
   read(*,*)
endif


return
end function ran0
