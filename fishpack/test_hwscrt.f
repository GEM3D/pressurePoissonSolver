      program test_hwscrt
      implicit none
      integer maxmx, maxmy
      parameter(maxmx = 5000, maxmy = 5000)

      double precision f(maxmx+1,maxmy+1)
      double precision u_exact(maxmx+1,maxmy+1)

c     # wsize : 4*(N+1) + (13 + INT(LOG2(N+1)))*(M+1)
c     # log2(m+1) approx. log10(m+1)/log10(2) approx. 4*log10(m)
      integer wsize, logmax
      parameter(logmax = 4*4)
      parameter(wsize = 4*(maxmy+1) + (13 + logmax)*(maxmx+1))
      double precision  work(wsize)

      double precision exp_error(10), compare
      logical check(20)

      double precision lambda, pertrb, dx, dy, ax, ay, bx, by, pi,
     1      uij, luij, error, t0, t1, dmaxmyp1
      integer i,j, mbdcnd, nbdcnd, ierror, choice, mx, my, idimf
      integer run_option, lstart, lend

      double precision BDC(maxmx+1), BDD(maxmx+1)
      double precision BDA(maxmy+1), BDB(maxmy+1)
      double precision xe(maxmx+1)
      double precision ye(maxmy+1)

      common /const_pi/ pi

      pi = 4.d0*datan(1.d0)

      exp_error(1) = 3.2905176298125838E-04
      exp_error(2) = 1.8030021919912542d-13
      exp_error(3) = 1.6164847238542279d-13
      exp_error(4) = 1.4896070488212843d-14
      exp_error(5) = 5.0620402068206338d-06
      exp_error(6) = 1.3239388303443178d-06
      exp_error(7) = 5.1503176130340831d-05
      exp_error(8) = 1.0904658141409953d+01
      exp_error(9) = 7.5652624150062769E-02
      exp_error(10) = 3.9115852588400557E-04

      !! --------------------------------------------------
      !! check that input dimensions don't exceed static memory

      dmaxmyp1 = maxmy + 1
      if (wsize < 4*(maxmy+1) + (13 +
     &      int(dlog10(dmaxmyp1)/dlog10(2.d0)))*(maxmx+1)) then
         write(6,*) 'wsize is not large enough'
         stop
      endif

      !! -------------------------------------------------
      !! Run test suite, or time hwscrt?
      write(6,*) 'Enter problem choice : '
      write(6,*) '   1. Run test suite'
      write(6,*) '   2. Test accuracy and timing of hwscrt'
      read(5,*) run_option

      if (run_option == 1) then
         mx = 100
         my = 100
         lstart = 1 !! start and stop of
         lend = 10

      else
c        !! --------------------------------------------------
c        !! Get input data
         write(6,*) 'Enter problem choice : '
         write(6,*) 'Test solutions :'
         write(6,*) '    1. u = sin(2*pi*x)*sin(2*pi*y)'
         write(6,*) '    2. u = 1  '
         write(6,*) '    3. u = x+y '
         write(6,*) '    4. u = 0.5*((x-0.5)**2 + (y-0.5))'
         write(6,*) '    5. u = y**2*x*4       '
         write(6,*) '    6. u = exp(x + y - 1)'
         write(6,*) '    7. u = exp(-alpha*r**2)'
         write(6,*) '    8. u = exp(-w*r**2/2) + exp(x + y)'
         write(6,*) '    9. u = delta '
         write(6,*) '   10. u = sin(pi*x)*cos(2*pi*y) '
         write(6,*) 'Enter problem choice : '
         read(5,*) choice
         if (choice .lt. 1 .or. choice .gt. 10) then
            write(6,*) 'Error : choice must be between 1-8'
            stop
         endif

         write(6, *) 'Enter problem size mx : '
         read(5,*) mx
         my = mx

c         write(6, *) 'Enter problem size my : '
c         read(5,*) my

         lstart = choice
         lend = choice

         if (mx .gt. maxmx) then
            write(6,*) 'M is too big; increase size of maxmx.'
            stop
         endif

      endif

c     !! -----------------------------------------------------
c     !! Setup domain and mesh size
c
c      my = mx  !! assume square domain
      ax = 0.d0
      bx = 1.d0
      ay = 0.d0
      by = 1.d0
      dx = (bx - ax)/float(mx)
      dy = (by - ay)/float(my)

      !! ----------------------------------
      !! Print parameters
      write(6,100) 'mx = ', mx
      write(6,100) 'my = ', my
      write(6,110) 'dx = ', dx
      write(6,110) 'dy = ', dy
      write(6,*) ' '

  100 format(A,I5)
  110 format(A,F16.8)

      !! -----------------------------------
      !! Set up mesh
      do i = 1,mx+1
         xe(i) = ax + (i - 1)*dx
      enddo

      do j = 1,my+1
         ye(j) = ay + (j - 1)*dy
      enddo

      !! -----------------------------------------------------
      !! Set boundary condition types
      !! Boundary conditions in x direction  - dirichlet conditions
      !! Boundary conditions in y direction  - dirichlet conditions
      mbdcnd = 1
      nbdcnd = 1

      !! -------------------------------------------------------
      !! Set up right hand side data and compute exact solution
      !! Dirichlet boundary data is set in first/last row/columns of
      !! f matrix.
c     # ------------------------------------
c     # Node solution
c     # ------------------------------------
      do j = 1,my+1
         BDA(j) = 0
         BDB(j) = 0
      enddo

      do i = 1,mx+1
         BDC(i) = 0
         BDD(i) = 0
      enddo
      lambda = 0.d0


      open(15,file='rhs.out');
      do choice = lstart, lend
         write(6,'(A,I2,A)') 'Setting up problem for test ', choice,
     &         '....'
         do i = 1,mx+1
            do j = 1,my+1
               call compute_u_all(choice,xe(i),ye(j),uij,luij)
               if (i == 1 .or. i == mx+1
     &               .or. j == 1 .or. j == my+1) then
c                 # Dirichlet data passed in in f.
                  f(i,j) = uij
                  write(15,'(2I5,E24.16)') i,j, 0.0
               else
                  f(i,j) = luij
                  write(15,'(2I5,E24.16)') i,j, f(i,j)
               endif
               u_exact(i,j) = uij
            enddo
         enddo
         close(15)

         idimf = maxmx+1

         check(choice) = .true.
c        # Call poisson solver.  Solution is returned in f.
         write(6,'(A)') 'Solving Poisson problem ...'
         call cpu_time(t0)
         call HWSCRT(ax,bx,mx,mbdcnd,bda,bdb,ay,by,my,nbdcnd, bdc, bdd,
     1         lambda,f,idimf,pertrb,ierror,work)
         call cpu_time(t1)
         write(6,'(A,F12.4)') 'Elapsed time for node solution : ', t1-t0

         if (ierror /= 0) then
            write(6,'(A,I1)') 'hwscrt error :  ierror = ', ierror
            write(6,*) 'See file hwscrt.f for documentation'
            stop
         endif

         if (work(1) > wsize) then
            write(6,*) 'hwscrt error : wsize is not large enough; ',work(1)
            write(6,*) 'See file hwscrt.f for documentation'
            stop
         endif

         error = 0.d0
         do i = 2,mx
            do j = 2,my
               error = dmax1(error,abs(u_exact(i,j) - f(i,j)))
            enddo
         enddo
         write(6,'(A,1PE30.16)') 'Error in computed solution : ', error
         if (mx .eq. 100) then
            write(6,'(A,1PE30.16)') 'Expected error :             ',
     &            exp_error(choice)

            compare = dlog10(error/exp_error(choice))

            if (compare .lt. 1.d0) then
c                write(6,*) 'HWSCRT appears to be working properly ',
c      &               'on test ', choice
               check(choice) = .true.
            else
c               write(6,*) 'WARNING : There may be a problem with HWSCRT'
               check(choice) = .false.
            endif
         endif
         write(6,*) ' '
      enddo

      if (run_option == 1 .or. mx == 100) then
         do choice = lstart,lend
            if (check(choice)) then
               write(6,'(A,I2,A)') 'Test ', choice, ' passed'
            else
               write(6,'(A,I2,A)') 'Test ',choice, ' failed.'
            endif
         enddo
         write(6,*) ' '
      endif

      end

      subroutine compute_u_all(choice,x,y,u,lu)
      implicit none

      double precision x,y,u,lu, pi, w, r2, alpha,d, R,
     &      dxy, xi, yj, th
      integer choice, N, i

      common /const_pi/ pi

      !! These are various test problems.  A second order Poisson solver,
      !! using the standard five point stencil should be able to solve for
      !! quadratic expressions exactly;  Others may be exact in special cases

      if (choice == 1) then
         !! Exact, at least for square domains
         u = dsin(2.d0*pi*x)*dsin(2.d0*pi*y)
         lu = -2*(2*pi)**2*u
      elseif (choice == 2) then
         !! Exact
         u = 1.d0
         lu = 0
      elseif (choice == 3) then
         !! Exact
         u = x + y
         lu = 0
      elseif (choice == 4) then
         !! Exact
         u = 0.5d0*(x - 0.5d0)**2 + 0.5d0*(y - 0.5d0)**2
         lu = 2.d0
      elseif (choice == 5) then
         !! Second order accurate results
         u = y**2*x**4
         lu = 2.d0*x**2*(6*y**2 + x**2)
      elseif (choice == 6) then
         !! Second order accurate results
         u = dexp(x + y - 1.d0)
         lu = 2.d0*dexp(x + y - 1.d0)
      elseif (choice == 7) then
         alpha = 4.d0
         r2 = (x - 0.1d0)**2 + y**2
         u = dexp(-r2*alpha)
         lu = (-4.d0*alpha + 4.d0*alpha**2*r2)*dexp(-r2*alpha)

         r2 = (x + 1.0d0)**2 + (y - 0.2)**2
         u = u - 2*exp(-r2*alpha)
         lu = lu -
     &         2.d0*(-4.d0*alpha + 4.d0*alpha**2*r2)*dexp(-r2*alpha)

         r2 = (x - 0.7d0)**2 + (y - 0.8d0)**2
         u = u + exp(-r2*alpha)
         lu = lu +
     &         (-4.d0*alpha + 4.d0*alpha**2*r2)*dexp(-r2*alpha)

      elseif (choice == 8) then
         w = 80000.d0
         r2 = (x-0.5001d0)**2 + (y-0.5001d0)**2
         u = dexp(-w*r2/2.d0) + dexp(x+y)
         lu = (w**2*r2 - 2*w)*dexp(-w*r2/2.d0) + 2.d0*dexp(x + y)

      elseif (choice == 9) then
c        Assume circle centered at (0.5, 0.5)
         N = 500
         R = 0.3d0
         w = 0.05d0
         r2 = (x - 0.5d0)**2 + (y - 0.5d0)**2
         dxy = dsqrt(r2) - R
         if (dabs(dxy) < w) then
            lu = (w - dabs(dxy))/(w*w)
         else
            lu = 0.d0
         endif
         if (dxy > 0) then
            u = dxy
         else
            u = 0
         endif
         return
         if (dabs(dxy) < w) then
            d = 0.d0
            do i = 0,N-1
               th = 2.d0*pi*i/N
               xi = R*dcos(th) + 0.5d0
               yj = R*dsin(th) + 0.5d0
               r2 = (xi - x)**2 + (yj - y)**2
               dxy = dsqrt(r2)
               if (dxy < w) then
                  d = d + (w - dxy)/(w*w)
               endif
            enddo
            lu = d*R*2.d0*pi/dfloat(N)
         else
            lu = 0
         endif
      elseif (choice == 10) then
!        ! Exact, at least for square domains
         u = dsin(pi*x)*dcos(2.d0*pi*y)
         lu = -5*(pi)**2*u
      endif

      end
