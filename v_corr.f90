!
!
! If anyone is reading this:
! subroutine spectral_density is adapted from 
! http://www.physics.unlv.edu/~pang/comp/code51.f  , 
! which comes with the heading:
! Please Note:                                                         
!                                                                      
! (1) This computer program is part of the book, "An Introduction to   
!     Computational Physics," written by Tao Pang and published and    
!     copyrighted by Cambridge University Press in 1997.               
!
!

program velocity_correlation
  implicit none
      
    integer,     parameter :: n_frames_total  = 1500
    integer,     parameter :: n_water         = 512
    integer,     parameter :: n_ox            = n_water 
    integer,     parameter :: n_hy            = n_ox * 2
    integer,     parameter :: n_atoms         = n_ox + n_hy
    integer,     parameter :: dp              = kind(0.0d0)
    integer,     parameter :: frames_msd      = 1400
    real(dp),    parameter :: delta_time      = 0.050d0
    real(dp),    parameter :: length_a        = 24.832132 ! Angstrom
    real(dp),    parameter :: length_b        = 24.832132 ! Angstrom
    real(dp),    parameter :: length_c        = 24.832132 ! Angstrom 
    real(dp),    parameter :: box_volume      = length_a * length_b * length_c
    real(dp),    parameter :: o_mass          = 15.995 ! amu
    real(dp),    parameter :: h_mass          = 1.008  ! amu
    real(dp),    parameter :: h2o_mass        = o_mass  + h_mass  + h_mass
    real(dp),    parameter :: system_mass     = n_water * h2o_mass
    real(dp),    parameter :: pi              = acos(-1.0d0)
    real(dp)               :: o_vel (n_frames_total, n_ox, 3)
    real(dp)               :: h_vel (n_frames_total, n_hy, 3)
    real(dp)               :: v_corr (n_frames_total)


    call read_water_xyz (o_vel, h_vel)
    call center_mass_removal (o_vel, h_vel)
    call v_correlation (o_vel, h_vel, v_corr)
    call spectral_density(v_corr) 
    call diffusion_vcorr(v_corr)



  contains

    subroutine read_water_xyz(o_vel,h_vel)

      implicit none
      integer  :: t_step, o, h1, h2
      real(dp) :: o_vel(:,:,:), h_vel(:,:,:)

      open(unit=10, file='water.vxyz', status='old')

      do t_step = 1, n_frames_total
        do o = 1, n_water
          h2 = (2*o)
          h1 = h2-1
          read(10,*) o_vel(t_step, o,  1), o_vel(t_step, o,  2), o_vel(t_step, o,  3)
          read(10,*) h_vel(t_step, h1, 1), h_vel(t_step, h1, 2), h_vel(t_step, h1, 3)
          read(10,*) h_vel(t_step, h2, 1), h_vel(t_step, h2, 2), h_vel(t_step, h2, 3)
        end do 
      end do 

      close(unit=10)

    end subroutine read_water_xyz

    subroutine center_mass_removal(o_vel,h_vel)

      implicit none
      integer  :: t_step, o1, o2, h, h1, h2
      real(dp) :: o_vel(:,:,:), h_vel(:,:,:)
      real(dp) :: x_1, x_2, y_1, y_2, z_1, z_2
      real(dp) :: x_mid, y_mid, z_mid

      outer: do t_step = 1, n_frames_total

        x_mid = 0
        y_mid = 0
        z_mid = 0

        do o1 = 1, n_water
          x_mid = x_mid + o_mass*o_vel(t_step, o1, 1)
          y_mid = y_mid + o_mass*o_vel(t_step, o1, 2)
          z_mid = z_mid + o_mass*o_vel(t_step, o1, 3)
        end do 

        do h1 = 1, n_hy
          x_mid = x_mid + h_mass*h_vel(t_step, h1, 1)
          y_mid = y_mid + h_mass*h_vel(t_step, h1, 2)
          z_mid = z_mid + h_mass*h_vel(t_step, h1, 3)
        end do 

          x_mid = x_mid / system_mass
          y_mid = y_mid / system_mass
          z_mid = z_mid / system_mass

        o_correct: do o1 = 1, n_water
          o_vel(t_step, o1, 1) = o_vel(t_step, o1, 1) - x_mid
          o_vel(t_step, o1, 2) = o_vel(t_step, o1, 2) - y_mid
          o_vel(t_step, o1, 3) = o_vel(t_step, o1, 3) - z_mid
        end do o_correct    

        h_correct: do h1 = 1, n_hy
          h_vel(t_step, h1, 1) = h_vel(t_step, h1, 1) - x_mid
          h_vel(t_step, h1, 2) = h_vel(t_step, h1, 2) - y_mid
          h_vel(t_step, h1, 3) = h_vel(t_step, h1, 3) - z_mid
        end do h_correct

      end do outer

    end subroutine center_mass_removal

    subroutine v_correlation(o_vel,h_vel,v_corr)

      implicit none
      integer  :: t_step, o1, o2, h, h1, h2 
      integer  :: delta_frames, frame1, frame2
      real(dp) :: o_vel(:,:,:), h_vel(:,:,:)
      real(dp) :: v_corr(:)
      real(dp) :: x_1, x_2, y_1, y_2, z_1, z_2
      real(dp) :: dot_a, dot_b, dot_c, dot_sum

      do delta_frames = 1, frames_msd
        v_corr (delta_frames) = 0.0d0
      end do

      do delta_frames = 1, frames_msd
        dot_sum = 0.0d0
        do frame1 = 1, n_frames_total - delta_frames
          frame2 = frame1 + delta_frames
          do o1 = 1, n_ox
            x_1     = o_vel (frame1, o1, 1)
            y_1     = o_vel (frame1, o1, 2)
            z_1     = o_vel (frame1, o1, 3)            
            x_2     = o_vel (frame2, o1, 1)
            y_2     = o_vel (frame2, o1, 2)
            z_2     = o_vel (frame2, o1, 3)
            dot_a   = x_1 * x_2
            dot_b   = y_1 * y_2
            dot_c   = z_1 * z_2
            dot_sum = dot_sum + dot_a + dot_b + dot_c
          end do
          do h1 = 1, n_hy
            x_1     = h_vel (frame1, h1, 1)
            y_1     = h_vel (frame1, h1, 2)
            z_1     = h_vel (frame1, h1, 3)            
            x_2     = h_vel (frame2, h1, 1)
            y_2     = h_vel (frame2, h1, 2)
            z_2     = h_vel (frame2, h1, 3)
            dot_a   = x_1 * x_2
            dot_b   = y_1 * y_2
            dot_c   = z_1 * z_2
            dot_sum = dot_sum + dot_a + dot_b + dot_c
          end do          
        end do
        dot_sum = dot_sum / ((n_frames_total-(delta_frames+1))*n_atoms)
        v_corr (delta_frames) = dot_sum
      end do

      open(unit=11, file='v_corr.dat', status='replace')
      do delta_frames = 1, frames_msd
        write(11,*) delta_frames*delta_time, v_corr(delta_frames)
      end do

    end subroutine v_correlation


    subroutine diffusion_vcorr(v_corr)
      implicit none
      integer  :: delta_frames
      real(dp) :: v_corr(:)
      real(dp) :: trap

      trap = 0.0d0
      do delta_frames = 1, frames_msd
        trap = trap + 0.50d0*(v_corr(delta_frames)+v_corr(delta_frames+1))
      end do
      trap = trap*delta_time/3.0d0

      write(6,*) "diffusion from velocity correlation:", trap


    end subroutine diffusion_vcorr

    subroutine spectral_density(v_corr)
      implicit none
      integer, parameter :: grid = frames_msd
      integer            :: t_step, i_grid, j_grid
      real(dp)           :: v_corr(:)
      real(dp)           :: input_real(grid), input_imag(grid)
      real(dp)           :: output_real(grid), output_imag(grid)
      real(dp)           :: fourier_x = 2.0d0*pi/grid
      real(dp)           :: fourier_q

      do i_grid = 1, grid
        input_real (i_grid) = v_corr(i_grid)
        input_imag (i_grid) = 0.0d0
      end do

      do i_grid = 1, grid
        output_real (i_grid) = 0.0d0
        output_imag (i_grid) = 0.0d0
        do j_grid = 1, grid
          fourier_q = fourier_x*(i_grid-1)*(j_grid-1)
          output_real (i_grid) = output_real(i_grid) + &
                   input_real(j_grid)*cos(fourier_q) + &
                   input_imag(j_grid)*sin(fourier_q)         
          output_imag (i_grid) = output_real(i_grid) + &
                   input_imag(j_grid)*cos(fourier_q) + &
                   input_real(j_grid)*sin(fourier_q)
        end do
      end do

      open(unit=11, file='trans_real.dat', status='replace')
      open(unit=12, file='trans_imag.dat', status='replace')
      do i_grid = 1, grid
        write(11,*) i_grid, output_real(i_grid)/(2.0d0*pi)
        write(12,*) i_grid, output_imag(i_grid)/(2.0d0*pi)
      end do
      close(unit=11)
      close(unit=12)




    end subroutine spectral_density






end program velocity_correlation