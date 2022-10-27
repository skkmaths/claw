subroutine init_cond(co1)
   use comvar
   implicit none

   real    :: co1(-ng+1:nx+ng, -ng+1:ny+ng)

   real    :: x, y
   integer :: i, j

   print*,'Setting initial condition'

   final_time = 2.0

   ! periodicity conditions
   xperiod = yes
   yperiod = yes

   xmin = -1.0
   xmax =  1.0
   ymin = -1.0
   ymax =  1.0

   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny

   do j=1,ny
      do i=1,nx
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy
         co1(i,j) = sin(2*M_PI*x) * sin(2*M_PI*y)
      enddo
   enddo

end subroutine init_cond
