      program cubic_spline
      implicit none
      integer:: ndata,n,ndim,i,j,ibase
      real(8),allocatable:: dat(:,:),vb(:),vx(:),amat(:,:),amati(:,:)

c.....Read data from STDIN
      read(5,*) ndata
      n= ndata-1
      ndim= 4*n
      allocate(dat(2,ndata),vb(ndim),amat(ndim,ndim),vx(ndim)
     &     ,amati(ndim,ndim))
      do i=1,ndata
        read(5,*) dat(1:2,i)
      enddo
      
c.....Create vector b
      vb(1:ndim)= 0d0
      do i=1,n
        vb(2*i-1)= dat(2,i)
        vb(2*i  )= dat(2,i+1)
      enddo
      
c.....Create matrix A
c.....0th order
      do i=1,n
        amat(2*i-1,4*i-3)= 1d0
        amat(2*i-1,4*i-2)= dat(1,i)
        amat(2*i-1,4*i-1)= dat(1,i)**2
        amat(2*i-1,4*i  )= dat(1,i)**3
        amat(2*i  ,4*i-3)= 1d0
        amat(2*i  ,4*i-2)= dat(1,i+1)
        amat(2*i  ,4*i-1)= dat(1,i+1)**2
        amat(2*i  ,4*i  )= dat(1,i+1)**3
      enddo
c.....1st order derivative
      ibase= 2*n
      do i=1,n-1
        amat(ibase+i,4*i-2)= 1d0
        amat(ibase+i,4*i-1)= 2d0*dat(1,i+1)
        amat(ibase+i,4*i  )= 3d0*dat(1,i+1)**2
        amat(ibase+i,4*(i+1)-2)= -1d0
        amat(ibase+i,4*(i+1)-1)= -2d0*dat(1,i+1)
        amat(ibase+i,4*(i+1)  )= -3d0*dat(1,i+1)**2
      enddo
c.....2nd order derivative
      ibase= 2*n +n-1
      do i=1,n-1
        amat(ibase+i,4*i-1)= 2d0
        amat(ibase+i,4*i  )= 6d0*dat(1,i+1)
        amat(ibase+i,4*(i+1)-1)= -2d0
        amat(ibase+i,4*(i+1)  )= -6d0*dat(1,i+1)
      enddo
c.....Edge data, natural boundary condition (second deriv. = 0d0)
      ibase= 2*n +2*(n-1)
      amat(ibase+1,3)= 2d0
      amat(ibase+1,4)= 6d0*dat(1,1)
      amat(ibase+2,ndim-1)= 2d0
      amat(ibase+2,ndim  )= 6d0*dat(1,ndata)

c      do i=1,ndim
c        write(6,'(20f6.1)') amat(i,1:ndim),vb(i)
c      enddo

c.....Get inverse of matrix A
      amati(1:ndim,1:ndim)= 0d0
      call ludc_inv(ndim,amat,amati)

c.....Compute Ai*b to get x
      vx(1:ndim)= 0d0
      do i=1,ndim
        do j=1,ndim
          vx(i)=vx(i) +amati(i,j)*vb(j)
        enddo
      enddo

      
      write(6,'(4es15.7)') vx(1:ndim)


      end program cubic_spline
c-----------------------------------------------------------------------
c     Local Variables:
c     compile-command: "ifort -o spline cubic_spline.F lasubs.F"
c     End:
