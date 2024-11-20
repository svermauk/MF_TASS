! To calculate l2 error.
! Shivani Verma
! 20-03-2021

program l2_error
  implicit none
  integer :: i_s1,i_s2,nbin1,nbin2,n
  real*8 :: dum,sum,norm_sum,l2,fes
  real*8, allocatable :: s1(:),s2(:),fes1(:,:),fes2(:,:),diff(:,:)

  nbin1=97
  nbin2=100

  allocate(s1(nbin1))
  allocate(s2(nbin2))
  allocate(fes1(nbin1,nbin2))
  allocate(fes2(nbin1,nbin2))
  allocate(diff(nbin1,nbin2))

  ! Reading free energy
  open(unit=1,file='intpl_fes.dat',status='old')
  open(unit=2,file='fes_mf40.dat',status='old')
  open(unit=10,file='l2_error.out',status='unknown')

  fes1=0.d0
  fes2=0.d0
  fes=0.d0
  diff=0.d0
  sum=0.d0
  n=0

  do i_s1=1,nbin1
    do i_s2=1,nbin2
      read(1,'(3e16.8)') dum,dum,fes
      fes1(i_s1,i_s2)=fes
      read(2,'(3e16.8)') s1(i_s1),s2(i_s2),fes2(i_s1,i_s2)
      ! Use second and third if statement when there is a huge difference in
      ! free energy at nbin1 and nbin2 due to inaccurate free energy at the
      ! boundaries.
      !if(fes1(i_s1,i_s2).le.15.d0.and.i_s1.lt.nbin1.and.i_s2.lt.nbin2) then
        diff(i_s1,i_s2)=fes2(i_s1,i_s2)-fes1(i_s1,i_s2)
        diff(i_s1,i_s2)=diff(i_s1,i_s2)**2
        sum=sum+diff(i_s1,i_s2)
        n=n+1
        diff(i_s1,i_s2)=sqrt(diff(i_s1,i_s2))
        write(10,'(3e16.8)') s1(i_s1),s2(i_s2),diff(i_s1,i_s2)
      !else
      !  diff(i_s1,i_s2)=100.d0
      !  write(10,'(3e16.8)') s1(i_s1),s2(i_s2),diff(i_s1,i_s2)
      !end if
    end do
    read(1,*)
    read(2,*)
    write(10,*)
  end do

  norm_sum=sum/n
  l2=sqrt(norm_sum)
  print*, "No. of grid points =", n
  print*, "l2 error is: ",l2

  deallocate(s1)
  deallocate(s2)
  deallocate(fes1)
  deallocate(fes2)
  deallocate(diff)
  close(1)
  close(2)
  close(10)

end program
