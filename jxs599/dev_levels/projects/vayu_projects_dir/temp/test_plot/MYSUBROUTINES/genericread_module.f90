 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module readcode_vars 
       integer, parameter :: ok=0
       integer :: openstatus,alloctest=0
   end module  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   module genericread_module
      interface readfile
         module procedure readfileX,readfile_Srb,readdatapoints,readuniform
      end interface 
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine readfileX(filename,N,whereto)
            use readcode_vars 
            implicit none
               character (len=*) :: filename
               integer :: N,i
               real, pointer :: whereto(:)
                  nullify(whereto)
                  allocate (whereto(N),stat=alloctest)
                  open(unit=1,file=filename,&
               &  status="old",action="read",iostat=openstatus)
                  if(openstatus==ok) then
                     read (1,*) whereto
                  else
                    print *, filename,'NOT found'
                    stop
                  endif
                  close(1)
            return 
         end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine readfile_Srb(filename,N,whereto,fesc)   !!this is not in use as fesc is calc. in file
            use readcode_vars 
            implicit none
               character (len=*) :: filename
               integer :: N,i
               real, pointer :: whereto(:)
               real ::fesc,temp 
                 allocate (whereto(N),stat=alloctest)
                  open(unit=1,file=filename,&
               &  status="old",action="read",iostat=openstatus)
                  if(openstatus==ok) then
                     do i=1,N
                        read (1,*),temp
                        whereto(i)=fesc*temp
                     enddo
                  else
                    print *, filename,'NOT found'
                    stop
                  endif
                  close(1)
            return 
         end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine readdatapoints(filename,N,whereto,lower,upper,errdirn)
            use readcode_vars 
            implicit none
               character (len=*) :: filename
               integer :: N,errdirn,i
               real, pointer :: whereto(:),upper(:),lower(:)
               real ::temp,trash 
                  allocate (whereto(N),stat=alloctest)
                  allocate (upper(N),stat=alloctest)
                  allocate (lower(N),stat=alloctest)
                  open(unit=1,file=filename,&
               &  status="old",action="read",iostat=openstatus)
                  if(openstatus==ok) then
                     do i=1,N
                        read (1,*),trash,whereto(i),lower(i),upper(i)
!   30                   format(F9.6,F9.6)         
!                        print *,trash,whereto(i),lower(i),upper(i)
                     enddo
                  else
                    print *, filename,'NOT found'
                    stop
                  endif
                  close(1)
            return 
         end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine readuniform(filename,N,whereto,uniform)
            use readcode_vars 
            implicit none
               character (len=*) :: filename
               integer :: N,i
               logical :: uniform
               real, pointer :: whereto(:)
               real ::temp,trash 
                  allocate (whereto(N),stat=alloctest)
                  open(unit=1,file=filename,&
               &  status="old",action="read",iostat=openstatus)
                  if(openstatus==ok) then
                     do i=1,N
                        !read (1,20),trash,temp
                        read (1,*),trash,temp
                        whereto(i)=temp
   20                   format(F3.1XX,F5.3)         
                     enddo
                  else
                    print *, filename,'NOT found'
                    stop
                  endif
                  close(1)
            return 
         end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end module genericread_module




