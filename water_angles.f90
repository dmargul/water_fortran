program water_hbond_angle_distributions
  implicit none
      
    integer,  parameter :: n_frames_total  = 400
    integer,  parameter :: n_water         = 512
    integer,  parameter :: n_ox            = n_water 
    integer,  parameter :: n_hy            = n_ox * 2
    integer,  parameter :: n_atoms         = n_ox + n_hy
    integer,  parameter :: n_bins          = 100
    integer,  parameter :: n_snaps         = 40
    integer             :: i_snap
    integer             :: n_frames
    integer             :: frames_here
    integer,  parameter :: dp              = SELECTED_REAL_KIND(15)
    real(dp), parameter :: delta_time      = 0.150d0! picoseconds between frames
    real(dp), parameter :: dt              = 0.020d0 ! picoseconds between snapshots
    real(dp), parameter :: pi              = dacos(-1.0d0)
    real(dp), parameter :: angle_min       =   0.0d0 ! Degrees
    real(dp), parameter :: angle_max       = 180.0d0 ! Degrees
    real(dp), parameter :: angle_delta     = (angle_max - angle_min) / real(n_bins)
    real(dp), parameter :: length_a        = 24.832132d0 ! Angstrom
    real(dp), parameter :: length_b        = 24.832132d0 ! Angstrom
    real(dp), parameter :: length_c        = 24.832132d0 ! Angstrom 
    real(dp), parameter :: box_volume      = length_a * length_b * length_c
    real(dp), parameter :: o_mass          = 15.995d0 ! amu
    real(dp), parameter :: h_mass          =  1.008d0 ! amu
    real(dp), parameter :: h2o_mass        = o_mass  + h_mass  + h_mass
    real(dp), parameter :: system_mass     = n_water * h2o_mass
    real(dp)            :: o_pos (n_frames_total, n_ox, 3)
    real(dp)            :: h_pos (n_frames_total, n_hy, 3)
    real(dp)            :: alpha_dist (n_bins)
    real(dp)            ::  beta_dist (n_bins)
    real(dp)            :: theta_dist (n_bins)
    real(dp)            :: angle_axis (n_bins)
    real(dp)            :: alpha_dist_time (n_snaps, n_bins)
    real(dp)            ::  beta_dist_time (n_snaps, n_bins)
    real(dp)            :: theta_dist_time (n_snaps, n_bins)
    real(dp)            :: alpha_conv (n_snaps)
    real(dp)            ::  beta_conv (n_snaps)
    real(dp)            :: theta_conv (n_snaps)

    call read_water_xyz (o_pos, h_pos)
    call center_mass_removal (o_pos, h_pos)
    do i_snap = 1, n_snaps
      frames_here = nint(real(n_frames_total) / (real(n_snaps)/real(i_snap)) )
      call make_distributions (o_pos, h_pos, alpha_dist, beta_dist, theta_dist, angle_axis)
      call norm_distributions (alpha_dist, beta_dist, theta_dist)
      call save_distributions (alpha_dist, beta_dist, theta_dist, alpha_dist_time, beta_dist_time, theta_dist_time)
    end do
    call distribution_convergence (alpha_dist_time, beta_dist_time, theta_dist_time, alpha_conv, beta_conv, theta_conv)

  contains

    subroutine read_water_xyz(o_pos, h_pos)

      implicit none
      integer  :: t_step, o, h1, h2
      real(dp) :: o_pos(:,:,:), h_pos(:,:,:)

      open(unit=10, file='water.xyz', status='old')

      outer: do t_step = 1, n_frames_total
        inner_water: do o = 1, n_water
          h2 = (2*o)
          h1 = h2-1
          read(10,*) o_pos(t_step, o,  1), o_pos(t_step, o,  2), o_pos(t_step, o,  3)
          read(10,*) h_pos(t_step, h1, 1), h_pos(t_step, h1, 2), h_pos(t_step, h1, 3)
          read(10,*) h_pos(t_step, h2, 1), h_pos(t_step, h2, 2), h_pos(t_step, h2, 3)
        end do inner_water
      end do outer

      close(unit=10)

    end subroutine read_water_xyz

    subroutine center_mass_removal(o_pos, h_pos)

      implicit none
      integer :: t_step, o1, o2, h, h1, h2
      real(dp):: o_pos(:,:,:), h_pos(:,:,:)
      real(dp):: x_1, x_2, y_1, y_2, z_1, z_2
      real(dp):: x_mid, y_mid, z_mid

      outer: do t_step = 1, n_frames_total

        x_mid = 0
        y_mid = 0
        z_mid = 0

        inner_o: do o1 = 1, n_water
          x_mid = x_mid + o_mass*o_pos(t_step, o1, 1)
          y_mid = y_mid + o_mass*o_pos(t_step, o1, 2)
          z_mid = z_mid + o_mass*o_pos(t_step, o1, 3)
        end do inner_o

        inner_h1: do h1 = 1, n_hy
          x_mid = x_mid + h_mass*h_pos(t_step, h1, 1)
          y_mid = y_mid + h_mass*h_pos(t_step, h1, 2)
          z_mid = z_mid + h_mass*h_pos(t_step, h1, 3)
        end do inner_h1

          x_mid = x_mid / system_mass
          y_mid = y_mid / system_mass
          z_mid = z_mid / system_mass

        inner_o_correct: do o1 = 1, n_water
          o_pos(t_step, o1, 1) = o_pos(t_step, o1, 1) - x_mid
          o_pos(t_step, o1, 2) = o_pos(t_step, o1, 2) - y_mid
          o_pos(t_step, o1, 3) = o_pos(t_step, o1, 3) - z_mid
        end do inner_o_correct    

        inner_h_correct: do h1 = 1, n_hy
          h_pos(t_step, h1, 1) = h_pos(t_step, h1, 1) - x_mid
          h_pos(t_step, h1, 2) = h_pos(t_step, h1, 2) - y_mid
          h_pos(t_step, h1, 3) = h_pos(t_step, h1, 3) - z_mid
        end do inner_h_correct

      end do outer

    end subroutine center_mass_removal

    subroutine make_distributions(o_pos, h_pos, alpha_dist, beta_dist, theta_dist, angle_axis)

      implicit none
      integer  :: t_step, o1, o2, h, h2
      integer  :: i_OD, i_OA, i_HD, i_HA1, i_HA2, i_bin
      real(dp) :: o_pos(:,:,:), h_pos(:,:,:)
      real(dp) :: alpha_dist(:), beta_dist(:), theta_dist(:), angle_axis(:)
      real(dp) :: x_2, y_2, z_2
      real(dp) :: x_OD, y_OD, z_OD
      real(dp) :: x_OA, y_OA, z_OA
      real(dp) :: x_HD, y_HD, z_HD
      real(dp) :: x_HA, y_HA, z_HA
      real(dp) :: x_HA1, y_HA1, z_HA1
      real(dp) :: x_HA2, y_HA2, z_HA2
      real(dp) :: disp_x2, disp_y2, disp_z2, disp_r2
      real(dp) :: disp_x,  disp_y,  disp_z,  disp_r
      real(dp) :: dist, distance_oh ! for finding donor O for each H
      real(dp) :: dist_OA_HD, dist_OD_HD ! for alpha angle vectors
      real(dp) :: dist_OA_OD, dist_HD_OD ! for beta  angle vectors
      real(dp) :: dist_HD_OA, dist_HA_OA ! for theta angle vectors
      real(dp) :: dot_a, dot_b, dot_c
      real(dp) :: dot_alpha, dot_beta, dot_theta
      real(dp) :: alpha_now, beta_now, theta_now
      real(dp) :: a1_length, a2_length
      real(dp) :: b1_length, b2_length
      real(dp) :: t1_length, t2_length
      real(dp) :: alpha_1(3), alpha_2(3)
      real(dp) ::  beta_1(3),  beta_2(3)
      real(dp) :: theta_1(3), theta_2(3)

      init: do i_bin = 1, n_bins
        dist = angle_min + i_bin * angle_delta  
        angle_axis (i_bin) = dist
        alpha_dist (i_bin) = 0
         beta_dist (i_bin) = 0
        theta_dist (i_bin) = 0
      end do init

      do t_step = 1, frames_here
        do i_HD = 1, n_hy
          x_HD = h_pos(t_step, i_HD, 1)
          y_HD = h_pos(t_step, i_HD, 2)
          z_HD = h_pos(t_step, i_HD, 3)
          if (mod(i_HD,2).eq.0) then
            i_OD = floor(real(i_HD)/2)
          else
            i_OD = floor(real(i_HD)/2) + 1
          end if
          x_OD = o_pos(t_step, i_OD, 1)
          y_OD = o_pos(t_step, i_OD, 2)
          z_OD = o_pos(t_step, i_OD, 3)
          distance_oh = 5.0d0
          do o2 = 1, n_water
            if (o2.ne.i_OD) then
              x_2 = o_pos(t_step, o2, 1)
              y_2 = o_pos(t_step, o2, 2)
              z_2 = o_pos(t_step, o2, 3)
              disp_x   = x_2 - x_HD
              disp_x   = disp_x - length_a*nint(disp_x/length_a)
              disp_x2  = disp_x * disp_x              
              disp_y   = y_2 - y_HD
              disp_y   = disp_y - length_a*nint(disp_y/length_b)
              disp_y2  = disp_y * disp_y              
              disp_z   = z_2 - z_HD
              disp_z   = disp_z - length_a*nint(disp_z/length_c)
              disp_z2  = disp_z * disp_z
              dist     = sqrt(disp_x2 + disp_y2 + disp_z2)
            end if
            if (dist.le.distance_oh) then
              distance_oh = dist
              x_OA = x_2
              y_OA = y_2
              z_OA = z_2
              i_OA = o2
            end if
          end do
          i_HA1 = 2*i_OA
          x_HA = h_pos(t_step, i_HA1, 1)
          y_HA = h_pos(t_step, i_HA1, 2)
          z_HA = h_pos(t_step, i_HA1, 3)

          alpha_1(1) = x_OD - x_HD
          alpha_1(2) = y_OD - y_HD 
          alpha_1(3) = z_OD - z_HD  
          alpha_1(1) = alpha_1(1) - length_a*nint(alpha_1(1)/length_a)    
          alpha_1(2) = alpha_1(2) - length_b*nint(alpha_1(2)/length_b)    
          alpha_1(3) = alpha_1(3) - length_c*nint(alpha_1(3)/length_c)    

          alpha_2(1) = x_OA - x_HD 
          alpha_2(2) = y_OA - y_HD
          alpha_2(3) = z_OA - z_HD  
          alpha_2(1) = alpha_2(1) - length_a*nint(alpha_2(1)/length_a)    
          alpha_2(2) = alpha_2(2) - length_b*nint(alpha_2(2)/length_b)    
          alpha_2(3) = alpha_2(3) - length_c*nint(alpha_2(3)/length_c)  

          beta_1 (1) = x_OA - x_OD 
          beta_1 (2) = y_OA - y_OD  
          beta_1 (3) = z_OA - z_OD  
          beta_1 (1) = beta_1 (1) - length_a*nint(beta_1 (1)/length_a)
          beta_1 (2) = beta_1 (2) - length_b*nint(beta_1 (2)/length_b)
          beta_1 (3) = beta_1 (3) - length_c*nint(beta_1 (3)/length_c)

          beta_2 (1) = x_HD - x_OD 
          beta_2 (2) = y_HD - y_OD  
          beta_2 (3) = z_HD - z_OD   
          beta_2 (1) = beta_2 (1) - length_a*nint(beta_2 (1)/length_a) 
          beta_2 (2) = beta_2 (2) - length_b*nint(beta_2 (2)/length_b)
          beta_2 (3) = beta_2 (3) - length_c*nint(beta_2 (3)/length_c)

          theta_1(1) = x_HD - x_OA  
          theta_1(2) = y_HD - y_OA  
          theta_1(3) = z_HD - z_OA  
          theta_1(1) = theta_1(1) - length_a*nint(theta_1(1)/length_a)
          theta_1(2) = theta_1(2) - length_b*nint(theta_1(2)/length_b)
          theta_1(3) = theta_1(3) - length_c*nint(theta_1(3)/length_c)

          theta_2(1) = x_HA - x_OA  
          theta_2(2) = y_HA - y_OA  
          theta_2(3) = z_HA - z_OA  
          theta_2(1) = theta_2(1) - length_a*nint(theta_2(1)/length_a)
          theta_2(2) = theta_2(2) - length_b*nint(theta_2(2)/length_b)
          theta_2(3) = theta_2(3) - length_c*nint(theta_2(3)/length_c)

          a1_length = sqrt(alpha_1(1)**2 + alpha_1(2)**2 + alpha_1(3)**2)
          a2_length = sqrt(alpha_2(1)**2 + alpha_2(2)**2 + alpha_2(3)**2)

          b1_length = sqrt(beta_1 (1)**2 + beta_1 (2)**2 + beta_1 (3)**2)
          b2_length = sqrt(beta_2 (1)**2 + beta_2 (2)**2 + beta_2 (3)**2)

          t1_length = sqrt(theta_1(1)**2 + theta_1(2)**2 + theta_1(3)**2)
          t2_length = sqrt(theta_2(1)**2 + theta_2(2)**2 + theta_2(3)**2)

          dot_a = alpha_1(1) * alpha_2(1)
          dot_b = alpha_1(2) * alpha_2(2)
          dot_c = alpha_1(3) * alpha_2(3)
          dot_alpha = dot_a + dot_b + dot_c  
          
          dot_a = beta_1 (1) * beta_2 (1)
          dot_b = beta_1 (2) * beta_2 (2)
          dot_c = beta_1 (3) * beta_2 (3)
          dot_beta  = dot_a + dot_b + dot_c  
          
          dot_a = theta_1(1) * theta_2(1)
          dot_b = theta_1(2) * theta_2(2)
          dot_c = theta_1(3) * theta_2(3)
          dot_theta = dot_a + dot_b + dot_c  

          alpha_now =  dacos ( dot_alpha / (a1_length * a2_length) )
          alpha_now = (alpha_now * (180.0d0 / pi))
          i_bin = nint (alpha_now / angle_delta) + 1
          if (i_bin .le. n_bins) then
            alpha_dist(i_bin) = alpha_dist(i_bin) + 1 
          end if  

          beta_now  =  dacos ( dot_beta  / (b1_length * b2_length) )
          beta_now  = (beta_now * (180.0d0 / pi))
          i_bin = nint ( beta_now / angle_delta) + 1
          if (i_bin .le. n_bins) then
            beta_dist(i_bin) =  beta_dist(i_bin) + 1 
          end if 

          theta_now =  dacos ( dot_theta / (t1_length * t2_length) )
          theta_now = (theta_now * (180.0d0 / pi))
          i_bin = nint (theta_now / angle_delta) + 1
          if (i_bin .le. n_bins) then
            theta_dist(i_bin) = theta_dist(i_bin) + 1 
          end if 


        end do
        
      end do 


    end subroutine make_distributions


    subroutine norm_distributions (alpha_dist, beta_dist, theta_dist)

      implicit none
      real(dp) :: dist_area 
      real(dp) :: alpha_dist(:), beta_dist(:), theta_dist(:)
      integer  :: i_bin

      dist_area = 0.0d0
      do i_bin = 1, n_bins
        dist_area = dist_area + alpha_dist(i_bin)
      end do
      dist_area = dist_area * angle_delta
      do i_bin = 1, n_bins
        alpha_dist(i_bin) = alpha_dist(i_bin) / dist_area
      end do      

      dist_area = 0.0d0
      do i_bin = 1, n_bins
        dist_area = dist_area +  beta_dist(i_bin)
      end do
      dist_area = dist_area * angle_delta
      do i_bin = 1, n_bins
         beta_dist(i_bin) =  beta_dist(i_bin) / dist_area
      end do

      dist_area = 0.0d0
      do i_bin = 1, n_bins
        dist_area = dist_area + theta_dist(i_bin)
      end do
      dist_area = dist_area * angle_delta
      do i_bin = 1, n_bins
        theta_dist(i_bin) = theta_dist(i_bin) / dist_area
      end do

    end subroutine norm_distributions

    subroutine save_distributions (alpha_dist, beta_dist, theta_dist, alpha_dist_time, beta_dist_time, theta_dist_time)

      implicit none
      integer  :: i_bin
      real(dp) :: alpha_dist(:), beta_dist(:), theta_dist(:)
      real(dp) :: alpha_dist_time(:,:), beta_dist_time(:,:), theta_dist_time(:,:)

      do i_bin = 1, n_bins
        alpha_dist_time(i_snap, i_bin) = alpha_dist(i_bin)
         beta_dist_time(i_snap, i_bin) =  beta_dist(i_bin)
        theta_dist_time(i_snap, i_bin) = theta_dist(i_bin)
      end do

    end subroutine save_distributions

    subroutine distribution_convergence(alpha_conv_time, beta_conv_time, theta_conv_time, alpha_conv, beta_conv, theta_conv)

      implicit none
      integer  :: i_snap, i_bin
      real(dp) :: alpha_conv_time(:,:), beta_conv_time(:,:), theta_conv_time(:,:)
      real(dp) :: alpha_conv(:), beta_conv(:), theta_conv(:)
      real(dp) :: squ_1, squ_2

      open(unit=11, file='alpha_dist.dat', status='replace')
      open(unit=12, file= 'beta_dist.dat', status='replace')
      open(unit=13, file='theta_dist.dat', status='replace')

      open(unit=21, file='alpha_conv.dat', status='replace')
      open(unit=22, file= 'beta_conv.dat', status='replace')
      open(unit=23, file='theta_conv.dat', status='replace')      

      do i_snap = 1, n_snaps
        alpha_conv(i_snap) = 0.0d0
         beta_conv(i_snap) = 0.0d0
        theta_conv(i_snap) = 0.0d0
        do i_bin = 1, n_bins
          squ_1 = alpha_conv_time(n_snaps, i_bin) - alpha_conv_time(i_snap, i_bin)
          alpha_conv(i_snap) = alpha_conv(i_snap) + abs((squ_1))

          squ_1 =  beta_conv_time(n_snaps, i_bin) - beta_conv_time (i_snap, i_bin)
           beta_conv(i_snap) =  beta_conv(i_snap) + abs((squ_1))

          squ_1 = theta_conv_time(n_snaps, i_bin) - theta_conv_time(i_snap, i_bin)
          theta_conv(i_snap) = theta_conv(i_snap) + abs((squ_1))

        end do
          alpha_conv(i_snap) = alpha_conv(i_snap) / real(n_bins)
           beta_conv(i_snap) =  beta_conv(i_snap) / real(n_bins)
          theta_conv(i_snap) = theta_conv(i_snap) / real(n_bins)

      end do

      do i_snap = 1, n_snaps
        write(21,*) i_snap, alpha_conv(i_snap)
        write(22,*) i_snap,  beta_conv(i_snap)
        write(23,*) i_snap, theta_conv(i_snap)
      end do


      do i_bin = 1, n_bins
        write(11,*) angle_axis(i_bin), alpha_dist(i_bin)   
        write(12,*) angle_axis(i_bin),  beta_dist(i_bin)
        write(13,*) angle_axis(i_bin), theta_dist(i_bin)
      end do 

      close(unit=11)
      close(unit=12)
      close(unit=13)

      close(unit=21)
      close(unit=22)
      close(unit=23)


    end subroutine distribution_convergence

end program water_hbond_angle_distributions