 

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
 
 
integer   :: x, y, z, w  !, xx, yy, zz, ww
real(8)   :: ss, mu, mult
real(8)   :: fs, fx, fy, fz, fw  
real(8)       :: div_wt, div_alone, div_alone2, div_wt2  
real(8)   :: s1, s2
real(8)   :: xid
real(8)       :: ran0, xran, xran2, xran3, xran4
integer   :: i, iruns 
logical   :: created
character*30 :: arg
 
 

   
call get_command_argument(1, arg) 
open(10, file=arg, form='formatted' )
 

 
 

ss = 0.
DO iruns=1,100000
!write(*,*) iruns


 
 mu = 1E-7  !1E-4
 mult = 0*500000  !500




s1 = 0.999  !0.9
s2 = 0.999  !0.9

 
div_wt = 0.15
div_wt2 = s2*div_wt
 
 
 
div_alone  = s1*div_wt  !0.9*
if (s2 <= s1) then
   div_alone2 = s2*div_wt    ! ######### LINE   ##########   ! this was s2*div_wt
else
   div_alone2 = s1*div_wt
endif
  
   
 



  

x = n**2
y = 0
z = 0
w = 0

 

xid = 0.       
created = .false.
do while( .true. )   
   xid = xid + 1
 
  
 
  if (x+y+w+z==0) exit
  if ( created ) then    
	 write(10,'(1x,i20)') 1
	 exit
  endif
  
  
  if ( xid > 50000 ) then
      write(10,'(1x,i20)') 0
	  exit
  endif
  
   !xx = x
   !yy = y
   !zz = z
   !ww = w
 
 
   do i = 1,n**2
       
	   ss = x + y + z + w
	   xran = ran0()
	   if ( xran < real(x)/ss )  then
	       x = x - 1
	   else if ( xran < real(x+y)/ss )  then	
           y = y - 1
	   else if ( xran < real(x+y+z)/ss )  then		   
           z = z - 1
	   else if ( xran < real(x+y+z+w)/ss )  then		  
           w = w - 1
	   else
	       write(*,*) 'koo'
		   read(*,*)
	   endif
	   
	   
	   fs = ( div_wt*real(x) + div_wt2 *real(z) + div_alone*real(y) + div_alone2*real(w) )
	   fx = div_wt*real(x) /  fs
	   fy = div_alone*real(y) / fs
	   fz = div_wt2 *real(z) / fs
	   fw = div_alone2*real(w) / fs
	   
	   xran2 = ran0()
	   if ( xran2 < fx ) then
	  	   xran3 = ran0()
		   if ( xran3 < mult*mu   ) then 
		 	  xran4 = ran0()
		   	  if ( xran4 < 2*mu ) then
			      w = w + 1 
			  else
			       z = z + 1
			  endif
		   else
			   xran4 = ran0()
			   if ( xran4 < 2*mu ) then
			       y = y + 1
			   else
			       x = x + 1
			   endif
			endif		 

	   else if ( xran2 < fx+fy ) then
		   xran3 = ran0()
		  if ( xran3 < mult*mu   ) then 
			  xran4 = ran0()
			  if ( xran4 < mu ) then
				  created = .true.
				  exit
			  else
			      w = w + 1
			  endif
		   else
			   xran4 = ran0()
			   if ( xran4 < mu ) then
				   created = .true.
				   exit
			   else
			       y = y + 1
			   endif
			endif		 
		
	   else if ( xran2 < fx+fy+fz ) then
		 if ( ran0() < 2*mu ) then
		      w = w + 1 
		 else
		     z = z + 1 
		 endif
		 
	   else if ( xran2 < fx+fy+fz+fw ) then
		 if ( ran0() < mu ) then
			 created = .true.
			 exit
		 else
		     w = w + 1
		 endif

       else
           write(*,*) 'boo', xran2, fx+fy+fz+fw
           read(*,*)	
	 
	 endif
	   
   enddo    

   
 
enddo    



ENDDO


end program 

















 












real(8) function ran0()
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
