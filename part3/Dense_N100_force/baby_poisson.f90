module type_defs
  integer, parameter:: sp =kind(1.0),&
    dp = selected_real_kind(2*precision(1.0_sp)),&
    qp = selected_real_kind(2*precision(1.0_dp))
end module type_defs

module problem_setup
  use type_defs
  implicit none
    integer,  parameter :: Nx = 32
    logical, parameter :: dirichlet_bc = .false.
    logical, parameter :: use_direct = .false.
    character (len=40) :: method= 'Jacobi'
    real(dp), parameter :: TOL = 1.0e-12_dp
end module problem_setup

module arrs
  use type_defs
  implicit none
  real(dp), allocatable, dimension(:) :: u,b,x
end module arrs

module afuns
  
contains
  subroutine apply_1D_laplacian_D(au,u,n)
    use type_defs
    implicit none
    integer, intent(in) :: n
    real(dp), intent(out) :: au(n)
    real(dp), intent(in)  ::  u(n)
    integer :: i
    Au(1) = (u(2) - 2.0_dp*u(1)         )
    Au(n) = (     - 2.0_dp*u(n) + u(n-1))
    do i = 2,n-1
     Au(i)= (u(i+1) - 2.0_dp*u(i) + u(i-1))
    end do
    
  end subroutine apply_1D_laplacian_D

  subroutine apply_1D_laplacian_N(au,u,n)
    use type_defs
    implicit none
    integer, intent(in) :: n
    real(dp), intent(out) :: au(n)
    real(dp), intent(in)  ::  u(n)
    integer :: i
    Au(1) =   (u(2)   -        u(1)         )
    Au(n) =   (       -        u(n) + u(n-1))
    do i = 2,n-1
       Au(i)= (u(i+1) - 2.0_dp*u(i) + u(i-1))
    end do
    
  end subroutine apply_1D_laplacian_N

end module afuns

module iterative_solver_D
  use type_defs
  implicit none
  real(dp), allocatable, dimension(:) :: r_isd,x_isd,Ar_isd
  real(dp) :: TOL_ISD
contains
  subroutine set_tol_isd(tol)
    implicit none
    real(dp) :: tol
    TOL_ISD = tol
  end subroutine set_tol_isd
  
  subroutine allocate_isd(n)
    implicit none
    integer, intent(in) :: n
    allocate(r_isd(n),x_isd(n),Ar_isd(n))
  end subroutine allocate_isd

  subroutine deallocate_isd
    implicit none
    deallocate(r_isd,x_isd,Ar_isd)
  end subroutine deallocate_isd
  
  subroutine steep_descent_d(x,b,n)
    use type_defs
    use afuns
    implicit none
    integer, intent(in) :: n
    real(dp), intent(inout)  :: x(n)
    real(dp), intent(inout)  :: b(n)
    real(dp) :: res_norm2,res_norm20,alpha
    integer :: iter    
    x_isd = 0.0_dp
    r_isd = b
    res_norm20 = sum(r_isd**2)
    res_norm2 = res_norm20
    iter = 0
    do while ((res_norm2/res_norm20 .gt. TOL_ISD**2) &
         .and. (iter .lt. 10000)) 
       call apply_1D_laplacian_D(Ar_isd,r_isd,n)
       alpha = res_norm2 / sum(r_isd*Ar_isd)
       x_isd = x_isd + alpha*r_isd
       r_isd = r_isd - alpha*Ar_isd       
       res_norm2 = sum(r_isd**2)
       iter = iter + 1
       write(*,*) iter, sqrt(res_norm2) 
    end do
    x = x_isd
  end subroutine steep_descent_d


subroutine Gauss_d(x,b,A,n)
  use type_defs
  use afuns
  implicit none
  integer, intent(in) :: n
  real(dp), dimension(:,:), intent(in) :: A
  real(dp), intent(inout)  :: x(n)
  real(dp), intent(in)  :: b(n)
  real(dp), allocatable, dimension(:) :: x_temp
  real(dp) :: r
  integer :: iter,i,j
  x_temp = x+1
  iter = 0
  do while (( sum( (x_temp-x)**2) .gt. TOL_ISD**2) &
       .and. (iter .lt. 1e12))
        x_temp = x
        do i = 1,n
            r = b(i)
            do j = 1,n
                r = r - A(i,j)*x(j)
            end do
            x(i) = x(i) + r/A(i,i)
        end do
        iter = iter + 1
end do
end subroutine Gauss_d

 subroutine Jacobi_d(x,b,A,n)
    use type_defs
    use afuns
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), intent(inout)  :: x(n)
    real(dp), intent(in)  :: b(n)
    real(dp), allocatable, dimension(:) :: x_temp
    real(dp) :: r
    integer :: iter,i,j
    x_temp = x+1
    iter = 0

    do while (( sum( (x_temp-x)**2) .gt. TOL_ISD**2) &
         .and. (iter .lt. 1e12))
          x_temp = x
          do i = 1,n
              r = b(i)
              do j = 1,n
                  r = r - A(i,j)*x_temp(j)
              end do
              x(i) = x_temp(i) + r/A(i,i)
          end do
          iter = iter + 1
    end do
end subroutine Jacobi_d
end module iterative_solver_D

module iterative_solver_N
  use type_defs
  implicit none
  real(dp), allocatable, dimension(:) :: r_isn,x_isn,Ar_isn
  real(dp) :: TOL_ISN
contains

  subroutine set_tol_isn(tol)
    implicit none
    real(dp) :: tol
    TOL_ISN = tol
  end subroutine set_tol_isn
  
  subroutine allocate_isn(n)
    implicit none
    integer, intent(in) :: n
    allocate(r_isn(n),x_isn(n),Ar_isn(n))
  end subroutine allocate_isn
  
  subroutine deallocate_isn
    implicit none
    deallocate(r_isn,x_isn,Ar_isn)
  end subroutine deallocate_isn


subroutine steep_descent_N(x,b,n)
  use type_defs
  use afuns
  implicit none
  integer, intent(in) :: n
  real(dp), intent(inout)  :: x(n)
  real(dp), intent(inout)  :: b(n)
  real(dp) :: res_norm2,res_norm20,alpha
  integer :: iter,i
  x_isn = 0.0_dp
  r_isn = b
    call apply_1D_laplacian_N(Ar_isn,r_isn,n)
    do i = 1,n
    write(*,*) Ar_isn(i)
    end do
  res_norm20 = sum(r_isn**2)
  res_norm2 = res_norm20
  iter = 0
  do while ((res_norm2/res_norm20 .gt. TOL_ISN**2) &
       .and. (iter .lt. 10000))
     call apply_1D_laplacian_N(Ar_isn,r_isn,n)
     alpha = res_norm2 / sum(r_isn*Ar_isn)
     x_isn = x_isn + alpha*r_isn
     r_isn = r_isn - alpha*Ar_isn
     res_norm2 = sum(r_isn**2)
     iter = iter + 1
     write(*,*) iter, sqrt(res_norm2)
  end do
  x = x_isn
end subroutine steep_descent_N

subroutine Gauss_general(x,b,A,n)
  use type_defs
  use afuns
  implicit none
  integer, intent(in) :: n
  real(dp), dimension(:,:), intent(in) :: A
  real(dp), intent(inout)  :: x(n)
  real(dp), intent(in)  :: b(n)
  real(dp), allocatable, dimension(:) :: x_temp
  real(dp) :: r
  integer :: iter,i,j
  x_temp = x+1
  iter = 0
  do while (( sum( (x_temp-x)**2) .gt. TOL_ISN**2) &
       .and. (iter .lt. 1e12))
        x_temp = x
        do i = 1,n
            r = b(i)
            do j = 1,n
                r = r - A(i,j)*x(j)
            end do
            x(i) = x(i) + r/A(i,i)
        end do
        iter = iter + 1
end do
end subroutine Gauss_general
subroutine Gauss_N(x,b,A,n)
use type_defs
use afuns
implicit none
integer, intent(in) :: n
real(dp), dimension(:,:), intent(in) :: A
real(dp), intent(inout)  :: x(n)
real(dp), intent(in)  :: b(n)
real(dp), allocatable, dimension(:) :: x_temp,d
real(dp), allocatable, dimension(:,:) :: A_res,A_d,v,U,inv_C,A_inv_u,U_T,U_T_A_inv_U,inv_U_T_A_inv_U
real(dp) :: r
integer :: iter,i,j,rot_num,it_num
allocate(A_res(n,n))
allocate(A_d(n,n))
allocate(v(n,n))
allocate(U(n,3))
allocate(A_inv_u(n,3))
allocate(inv_C(3,3))
allocate(d(n))

do i = 1,n
    do j = 1,n
        A_res(i,j) = 1
        A_d(i,j) = A(i,j)
    end do
end do

A_res(1,1) = 2
A_d(1,1) = -2
A_res(n,n) = 2
A_d(n,n) = -2

Call jacobi_eigenvalue ( n, A_res, 10000, v, d, it_num, rot_num )
do i = 1,3
    inv_C(i,i) = 1/d(n+1-i)
end do
do i = 1,n
    do j = 1,3
        U(i,j) = V(n+1-i,n+1-j)
    end do
end do
do i = 1,3
call Gauss_general(A_inv_u(:,i),U(:,i),A_d,n)
end do
call Gauss_general(x,b,A_d,n)
U_T = transpose(U)
U_T_A_inv_U = matmul(U_T,A_inv_u)
inv_U_T_A_inv_U = inv(inv_C + U_T_A_inv_U)
x = x - matmul(A_inv_U,matmul(inv_U_T_A_inv_U ,matmul(U_T,x)))
end subroutine Gauss_N

function inv(A) result(Ainv)
  real(dp), dimension(:,:), intent(in) :: A
  real(dp), dimension(size(A,1),size(A,2)) :: Ainv

  real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function inv

subroutine Jacobi_general(x,b,A,n)
  use type_defs
  use afuns
  implicit none
  integer, intent(in) :: n
  real(dp), dimension(:,:), intent(in) :: A
  real(dp), intent(inout)  :: x(n)
  real(dp), intent(in)  :: b(n)
  real(dp), allocatable, dimension(:) :: x_temp
  real(dp) :: r
  integer :: iter,i,j
  x_temp = x+1
  iter = 0
  do while (( sum( (x_temp-x)**2) .gt. TOL_ISN**2) &
       .and. (iter .lt. 1e12))
        x_temp = x
        do i = 1,n
            r = b(i)
            do j = 1,n
                r = r - A(i,j)*x_temp(j)
            end do
            x(i) = x_temp(i) + r/A(i,i)
        end do
        iter = iter + 1
end do
end subroutine Jacobi_general

subroutine Jacobi_N(x,b,A,n)
use type_defs
use afuns
implicit none
integer, intent(in) :: n
real(dp), dimension(:,:), intent(in) :: A
real(dp), intent(inout)  :: x(n)
real(dp), intent(in)  :: b(n)
real(dp), allocatable, dimension(:) :: x_temp,d
real(dp), allocatable, dimension(:,:) :: A_res,A_d,v,U,inv_C,A_inv_u,U_T,U_T_A_inv_U,inv_U_T_A_inv_U
real(dp) :: r
integer :: iter,i,j,rot_num,it_num
allocate(A_res(n,n))
allocate(A_d(n,n))
allocate(v(n,n))
allocate(U(n,3))
allocate(A_inv_u(n,3))
allocate(inv_C(3,3))
allocate(d(n))

do i = 1,n
    do j = 1,n
        A_res(i,j) = 1
        A_d(i,j) = A(i,j)
    end do
end do

A_res(1,1) = 2
A_d(1,1) = -2
A_res(n,n) = 2
A_d(n,n) = -2

Call jacobi_eigenvalue ( n, A_res, 10000, v, d, it_num, rot_num )
do i = 1,3
    inv_C(i,i) = 1/d(n+1-i)
end do
do i = 1,n
    do j = 1,3
        U(i,j) = V(n+1-i,n+1-j)
    end do
end do
do i = 1,3
call Jacobi_general(A_inv_u(:,i),U(:,i),A_d,n)
end do
call Jacobi_general(x,b,A_d,n)
U_T = transpose(U)
U_T_A_inv_U = matmul(U_T,A_inv_u)
inv_U_T_A_inv_U = inv(inv_C + U_T_A_inv_U)
x = x - matmul(A_inv_U,matmul(inv_U_T_A_inv_U ,matmul(U_T,x)))
end subroutine Jacobi_N

subroutine jacobi_eigenvalue ( n, a, it_max, v, d, it_num, rot_num )

!*****************************************************************************80
!
!! JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
!
!  Discussion:
!
!    This function computes the eigenvalues and eigenvectors of a
!    real symmetric matrix, using Rutishauser's modfications of the classical
!    Jacobi rotation method with threshold pivoting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix, which must be square, real,
!    and symmetric.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, real ( kind = 8 ) V(N,N), the matrix of eigenvectors.
!
!    Output, real ( kind = 8 ) D(N), the eigenvalues, in descending order.
!
!    Output, integer ( kind = 4 ) IT_NUM, the total number of iterations.
!
!    Output, integer ( kind = 4 ) ROT_NUM, the total number of rotations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) bw(n)
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) g
  real ( kind = 8 ) gapq
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) rot_num
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) tau
  real ( kind = 8 ) term
  real ( kind = 8 ) termp
  real ( kind = 8 ) termq
  real ( kind = 8 ) theta
  real ( kind = 8 ) thresh
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) zw(n)

  do j = 1, n
    do i = 1, n
      v(i,j) = 0.0D+00
    end do
    v(j,j) = 1.0D+00
  end do

  do i = 1, n
    d(i) = a(i,i)
  end do

  bw(1:n) = d(1:n)
  zw(1:n) = 0.0D+00
  it_num = 0
  rot_num = 0

  do while ( it_num < it_max )

    it_num = it_num + 1
!
!  The convergence threshold is based on the size of the elements in
!  the strict upper triangle of the matrix.
!
    thresh = 0.0D+00
    do j = 1, n
      do i = 1, j - 1
        thresh = thresh + a(i,j) ** 2
      end do
    end do

    thresh = sqrt ( thresh ) / real ( 4 * n, kind = 8 )

    if ( thresh == 0.0D+00 ) then
      exit
    end if

    do p = 1, n
      do q = p + 1, n

        gapq = 10.0D+00 * abs ( a(p,q) )
        termp = gapq + abs ( d(p) )
        termq = gapq + abs ( d(q) )
!
!  Annihilate tiny offdiagonal elements.
!
        if ( 4 < it_num .and. &
             termp == abs ( d(p) ) .and. &
             termq == abs ( d(q) ) ) then

          a(p,q) = 0.0D+00
!
!  Otherwise, apply a rotation.
!
        else if ( thresh <= abs ( a(p,q) ) ) then

          h = d(q) - d(p)
          term = abs ( h ) + gapq

          if ( term == abs ( h ) ) then
            t = a(p,q) / h
          else
            theta = 0.5D+00 * h / a(p,q)
            t = 1.0D+00 / ( abs ( theta ) + sqrt ( 1.0D+00 + theta * theta ) )
            if ( theta < 0.0D+00 ) then
              t = - t
            end if
          end if

          c = 1.0D+00 / sqrt ( 1.0D+00 + t * t )
          s = t * c
          tau = s / ( 1.0D+00 + c )
          h = t * a(p,q)
!
!  Accumulate corrections to diagonal elements.
!
          zw(p) = zw(p) - h
          zw(q) = zw(q) + h
          d(p) = d(p) - h
          d(q) = d(q) + h

          a(p,q) = 0.0D+00
!
!  Rotate, using information from the upper triangle of A only.
!
          do j = 1, p - 1
            g = a(j,p)
            h = a(j,q)
            a(j,p) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = p + 1, q - 1
            g = a(p,j)
            h = a(j,q)
            a(p,j) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = q + 1, n
            g = a(p,j)
            h = a(q,j)
            a(p,j) = g - s * ( h + g * tau )
            a(q,j) = h + s * ( g - h * tau )
          end do
!
!  Accumulate information in the eigenvector matrix.
!
          do j = 1, n
            g = v(j,p)
            h = v(j,q)
            v(j,p) = g - s * ( h + g * tau )
            v(j,q) = h + s * ( g - h * tau )
          end do

          rot_num = rot_num + 1

        end if

      end do
    end do

    bw(1:n) = bw(1:n) + zw(1:n)
    d(1:n) = bw(1:n)
    zw(1:n) = 0.0D+00

  end do
!
!  Restore upper triangle of input matrix.
!
  do j = 1, n
    do i = 1, j - 1
      a(i,j) = a(j,i)
    end do
  end do
!
!  Ascending sort the eigenvalues and eigenvectors.
!
  do k = 1, n - 1

    m = k

    do l = k + 1, n
      if ( d(l) < d(m) ) then
        m = l
      end if
    end do

    if ( m /= k ) then

      t    = d(m)
      d(m) = d(k)
      d(k) = t

      w(1:n)   = v(1:n,m)
      v(1:n,m) = v(1:n,k)
      v(1:n,k) = w(1:n)

    end if

  end do

  return
end subroutine jacobi_eigenvalue

end module iterative_solver_N



program ins
  use type_defs
  use problem_setup
  use arrs
  use afuns
  use iterative_solver_N
  implicit none
  ! This program solves u_xx = f 
  ! on the domain [x] \in [0,1] with either Dirichlet or Neumann BC
  ! hx = 1/Nx
  real(dp) :: hx
  integer :: i,n_iter,N_sys,info  
  real(dp), allocatable, dimension(:,:) :: A
  real(dp), allocatable, dimension(:) :: action_of_A,u_inner
  integer, allocatable, dimension(:) ::  ipiv
  ! Set up the grid
  hx = 1.0_dp/real(Nx,dp)
  allocate(x(0:nx))
  allocate(u(0:nx))
  do i = 0,nx
   x(i) = real(i,dp)*hx
  end do

  if(dirichlet_bc) then
     ! For Dirichlet BC we don't solve for the boundary conditions. 
     N_sys = nx-1
  else
     ! For Neumann BC the solution at the first and last gridpoint is
     ! part of the solve.
     N_sys = nx+2
     ! Depending on what you want to do you may also need a nonsingular system...
  end if
  
  allocate(b(N_sys))

allocate(A(N_sys,N_sys),&
    action_of_A(N_sys),&
    u_inner(N_sys),&
    ipiv(N_sys))
if(dirichlet_bc) then
    do i = 1,N_sys
        u_inner = 0.0
        u_inner(i) = 1.0d0
        call apply_1D_laplacian_D(action_of_A,u_inner,N_sys)
        A(:,i) = action_of_A
    end do
else
    do i = 1,N_sys-1
        u_inner = 0.0
        u_inner(i) = 1.0d0
        call apply_1D_laplacian_N(action_of_A,u_inner,N_sys-1)
        A(:,i) = action_of_A
    end do
    do i = 1,N_sys-1
        A(i,N_sys) = 1
        A(N_sys,i) = 1
    end do
    A(N_sys,N_sys)=0
end if

  ! Set up the problem.
  ! We cook up a problem with a known solution.
  ! Say that u(x) = exp(-x), then we must have that
  ! f = u_xx  = exp(-x),
  ! u(x = 0) = exp(0), u(x = 1) = exp(-1)
  ! or
  ! u_x(x = 0) = -exp(0), u_x(x = 1) = -exp(-1)
  !
  ! We move the hx^2 from the denominator over to the right hand side of the
  ! system of equations.

  if (dirichlet_bc) then
     do i = 1,nx-1
        b(i) = hx*hx*exp(-x(i))
     end do
     ! We must also account for the boundary conditions
     b(1) = b(1) - exp(-x(0))
     b(nx-1) = b(nx-1) - exp(-x(nx))

     if (use_direct) then
        CALL DGETRF(N_sys,N_sys,A,N_sys,ipiv,INFO)
        CALL DGETRS('N',N_sys,1,A,N_sys,IPIV,b,N_sys,INFO)
        u(1:nx-1) = b 
    else
       if (method == "SD") then
            call allocate_isd(N_sys)
            call set_tol_isd(tol)
            call steep_descent_d(u(1:nx),b,N_sys)
            call deallocate_isd
       elseif  (method == "Gauss") then
           call allocate_isd(N_sys)
           call set_tol_isd(tol)
           call Gauss_d(u(1:nx-1),b,A,N_sys)
           call deallocate_isd
       elseif (method == "Jacobi") then
           call allocate_isd(N_sys)
           call set_tol_isd(tol)
           call Jacobi_d(u(1:nx-1),b,A,N_sys)
           call deallocate_isd
       end if
    end if
     u(0) = exp(-x(0))
     u(nx) = exp(-x(nx))
     write(*,*) maxval(abs(u - exp(-x)))
     
  else
     do i = 1,nx+1
        b(i) = hx*hx*exp(-x(i-1))
     end do
     ! We must also account for the boundary conditions
     ! Here we have that u(1) - 2*u(0) + u(-1) = h*h*exp(-x(0))
     ! And we approximate u_x(x = 0) = -exp(0) by
     ! u(1) - u(-1) = 2*h * (-exp(0))
     ! Plug in
     ! 2*u(1) - 2*u(0) = h*h*exp(-x(0)) + 2*h*(-exp(0)) 
     ! We scale this by 1/2 to make the matrix symmetric
     b(1) = 0.5_dp*b(1) - hx*exp(-x(0))
     b(nx+1) = 0.5_dp*b(nx+1) + hx*exp(-x(nx))
     b(nx+2) = 0
    if (use_direct) then
        CALL DGETRF(N_sys,N_sys,A,N_sys,ipiv,INFO)
        CALL DGETRS('N',N_sys,1,A,N_sys,IPIV,b,N_sys,INFO)
        u(0:nx) = b(1:nx+1)
     else
        if (method == "SD") then
             call allocate_isn(N_sys)
             call set_tol_isn(tol)
             call steep_descent_N(u(0:nx),b(1:nx+1),N_sys-1)
             call deallocate_isn
        elseif  (method == "Gauss") then
            call allocate_isn(N_sys)
            call set_tol_isn(tol)
            call Gauss_N(u(0:nx),b,A,N_sys-1)
            call deallocate_isn
        elseif (method == "Jacobi") then
            call allocate_isn(N_sys)
            call set_tol_isn(tol)
            call Jacobi_N(u(0:nx),b,A,N_sys-1)
            call deallocate_isn
        end if
     end if
    write(*,*) maxval(abs(u - (exp(-x) -1.0d0+exp(-1.0d0) )))
  end if



end program ins
