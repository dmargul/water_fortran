program water_polar_distributions
  implicit none
      
    integer,  parameter :: n_frames_total  = 5000
    integer,  parameter :: n_water         = 512
    integer,  parameter :: n_ox            = n_water 
    integer,  parameter :: n_hy            = n_ox * 2
    integer,  parameter :: n_atoms         = n_ox + n_hy
    integer,  parameter :: n_bins          = 200
    integer,  parameter :: n_snaps         = 100
    integer             :: i_snap
    integer             :: n_frames
    integer             :: frames_here
    integer,  parameter :: dp              = kind(0.0d0)
    real(dp), parameter :: q_H             =  0.25983   ! units
    real(dp), parameter :: q_O             = -2.0d0*q_H ! units
    real(dp), parameter :: convert         = 4.80321 ! from TINKER output to Debye
    real(dp), parameter :: delta_time      = 0.050d0 ! picoseconds between frames
    real(dp), parameter :: pi              = acos(-1.0d0)
    real(dp), parameter :: polar_min       = 0.0d0 ! Debye
    real(dp), parameter :: polar_max       = 5.0d0 ! Debye
    real(dp), parameter :: polar_delta     = (polar_max - polar_min ) / real(n_bins)
    real(dp), parameter :: length_a        = 24.832132 ! Angstrom
    real(dp), parameter :: length_b        = 24.832132 ! Angstrom
    real(dp), parameter :: length_c        = 24.832132 ! Angstrom 
    real(dp), parameter :: box_volume      = length_a * length_b * length_c
    real(dp), parameter :: o_mass          = 15.995 ! amu
    real(dp), parameter :: h_mass          = 1.008  ! amu
    real(dp), parameter :: h2o_mass        = o_mass  + h_mass  + h_mass
    real(dp), parameter :: system_mass     = n_water * h2o_mass
    real(dp)            :: o_pos (n_frames_total, n_ox, 3)
    real(dp)            :: h_pos (n_frames_total, n_hy, 3)
    real(dp)            :: o_stc (n_frames_total, n_ox, 3)
    real(dp)            :: h_stc (n_frames_total, n_hy, 3)
    real(dp)            :: o_ind (n_frames_total, n_ox, 3)
    real(dp)            :: h_ind (n_frames_total, n_hy, 3)
    real(dp)            :: distance_oo 
    real(dp)            :: distance_oh 
    real(dp)            :: distance_hh 
    real(dp)            ::    polar_axis (n_bins)
    real(dp)            ::   static_dist (n_bins)
    real(dp)            ::   induce_dist (n_bins)
    real(dp)            ::   monopl_dist (n_bins)
    real(dp)            ::   plrztn_dist (n_bins)
    real(dp)            ::   static_time (n_snaps, n_bins)
    real(dp)            ::   induce_time (n_snaps, n_bins)
    real(dp)            ::   monopl_time (n_snaps, n_bins)
    real(dp)            ::   plrztn_time (n_snaps, n_bins)
    real(dp)            ::   static_conv (n_snaps)
    real(dp)            ::   induce_conv (n_snaps)
    real(dp)            ::   monopl_conv (n_snaps)
    real(dp)            ::   plrztn_conv (n_snaps)
    real(dp)            ::  stc_avg_time (n_snaps)
    real(dp)            ::  ind_avg_time (n_snaps)
    real(dp)            :: mono_avg_time (n_snaps)
    real(dp)            :: pole_avg_time (n_snaps)    
    real(dp)            ::  stc_avg_conv (n_snaps)
    real(dp)            ::  ind_avg_conv (n_snaps)
    real(dp)            :: mono_avg_conv (n_snaps)
    real(dp)            :: pole_avg_conv (n_snaps)
    real(dp)            :: cell_dplsq_time(n_snaps)
    real(dp)            :: cell_dplsq_conv(n_snaps)

    call read_water_xyz_polar (o_pos, h_pos, o_stc, h_stc, o_ind, h_ind)
    call center_mass_removal (o_pos, h_pos)
    do i_snap = 1, n_snaps
      frames_here = nint(real(n_frames_total) / (real(n_snaps)/real(i_snap) ) )

      call molecular_dipoles(o_pos, h_pos, o_stc, h_stc, o_ind, h_ind, static_dist, induce_dist, monopl_dist, plrztn_dist, &
                             polar_axis, stc_avg_time, ind_avg_time, mono_avg_time, pole_avg_time, cell_dplsq_time)

      call norm_distributions (static_dist, induce_dist, monopl_dist, plrztn_dist)

      call save_distributions (static_dist, induce_dist, monopl_dist, plrztn_dist, &
                               static_time, induce_time, monopl_time, plrztn_time)
    end do

    call distribution_convergence (static_time, induce_time, monopl_time, plrztn_time, stc_avg_conv,           &
               ind_avg_conv, mono_avg_conv, pole_avg_conv, static_conv, induce_conv, monopl_conv, plrztn_conv, &
               cell_dplsq_time, cell_dplsq_conv, stc_avg_time, ind_avg_time, mono_avg_time, pole_avg_time)

  contains

    subroutine read_water_xyz_polar(o_pos, h_pos, o_stc, h_stc, o_ind, h_ind)

      implicit none
      integer  :: t_step, o, h1, h2
      real(dp) :: o_pos(:,:,:), h_pos(:,:,:)
      real(dp) :: o_stc(:,:,:), h_stc(:,:,:)
      real(dp) :: o_ind(:,:,:), h_ind(:,:,:)

      open(unit=10, file='water.xyz_polar', status='old')

      outer: do t_step = 1, n_frames_total
        inner_water: do o = 1, n_water
          h2 = (2*o)
          h1 = h2-1
          read(10,*) o_pos(t_step, o,  1), o_pos(t_step, o,  2), o_pos(t_step, o,  3), &
                     o_stc(t_step, o,  1), o_stc(t_step, o,  2), o_stc(t_step, o,  3), & 
                     o_ind(t_step, o,  1), o_ind(t_step, o,  2), o_ind(t_step, o,  3) 

          read(10,*) h_pos(t_step, h1, 1), h_pos(t_step, h1, 2), h_pos(t_step, h1, 3), &
                     h_stc(t_step, h1, 1), h_stc(t_step, h1, 2), h_stc(t_step, h1, 3), & 
                     h_ind(t_step, h1, 1), h_ind(t_step, h1, 2), h_ind(t_step, h1, 3)  

          read(10,*) h_pos(t_step, h2, 1), h_pos(t_step, h2, 2), h_pos(t_step, h2, 3), &
                     h_stc(t_step, h2, 1), h_stc(t_step, h2, 2), h_stc(t_step, h2, 3), & 
                     h_ind(t_step, h2, 1), h_ind(t_step, h2, 2), h_ind(t_step, h2, 3)            
        end do inner_water
      end do outer

      close(unit=10)

    end subroutine read_water_xyz_polar

    subroutine center_mass_removal(o_pos,h_pos)

      implicit none
      integer  :: t_step, o1, o2, h, h1, h2
      real(dp) :: o_pos(:,:,:), h_pos(:,:,:)
      real(dp) :: x_1, x_2, y_1, y_2, z_1, z_2
      real(dp) :: x_mid, y_mid, z_mid

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

    subroutine molecular_dipoles(o_pos, h_pos, o_stc, h_stc, o_ind, h_ind, static_dist, induce_dist, monopl_dist, plrztn_dist, &
                             polar_axis, stc_avg_time, ind_avg_time, mono_avg_time, pole_avg_time, cell_dplsq_time)

      implicit none
      integer  :: t_step, o1, o2, h, h1, h2, i_bin
      real(dp) :: o_pos(:,:,:), h_pos(:,:,:)
      real(dp) :: o_stc(:,:,:), h_stc(:,:,:)
      real(dp) :: o_ind(:,:,:), h_ind(:,:,:)
      real(dp) :: x_1, x_2, y_1, y_2, z_1, z_2
      real(dp) :: x_O,  y_O,  z_O
      real(dp) :: x_H1, y_H1, z_H1
      real(dp) :: x_H2, y_h2, z_H2
      real(dp) :: x_OH1, y_OH1, z_OH1
      real(dp) :: x_OH2, y_OH2, z_OH2
      real(dp) :: disp_x2, disp_y2, disp_z2, disp_r2
      real(dp) :: disp_x, disp_y, disp_z, disp_r
      real(dp) :: static_dist(:), induce_dist(:), monopl_dist(:), plrztn_dist(:), polar_axis(:)
      real(dp) :: stc_avg_time(:), ind_avg_time(:), mono_avg_time(:), pole_avg_time(:) 
      real(dp) :: cell_dplsq_time(:)
      real(dp) :: cell_dpl_x, cell_dpl_y, cell_dpl_z
      real(dp) :: mono_avg, stc_avg, ind_avg, pole_avg
      real(dp) :: pole_ind, pole_mono, pole_stc, pole_total
      real(dp) :: pole_x, pole_y, pole_z
      real(dp) :: dist


      init: do i_bin = 1, n_bins
         dist = polar_min  + i_bin*polar_delta  
         polar_axis(i_bin) = dist
        static_dist(i_bin) = 0
        induce_dist(i_bin) = 0
        monopl_dist(i_bin) = 0
        plrztn_dist(i_bin) = 0
      end do init
      mono_avg = 0.0d0
       stc_avg = 0.0d0
       ind_avg = 0.0d0
      pole_avg = 0.0d0

      cell_dplsq_time(i_snap) = 0.0d0
      outer: do t_step = 1, frames_here
        cell_dpl_x = 0.0d0
        cell_dpl_y = 0.0d0
        cell_dpl_z = 0.0d0
        inner_water: do o1 = 1, n_water
          h2   = (2*o1)
          h1   = h2-1
          x_O  = o_pos(t_step, o1, 1)
          y_O  = o_pos(t_step, o1, 2)
          z_O  = o_pos(t_step, o1, 3)
          x_H1 = h_pos(t_step, h1, 1)
          y_H1 = h_pos(t_step, h1, 2)
          z_H1 = h_pos(t_step, h1, 3)
          x_H2 = h_pos(t_step, h2, 1)
          y_H2 = h_pos(t_step, h2, 2)
          z_H2 = h_pos(t_step, h2, 3)

          cell_dpl_x = cell_dpl_x + x_O*q_O + x_H1*q_H + x_H2*q_H
          cell_dpl_y = cell_dpl_y + y_O*q_O + y_H1*q_H + y_H2*q_H
          cell_dpl_z = cell_dpl_z + z_O*q_O + z_H1*q_H + z_H2*q_H

          x_OH1 = x_H1 - x_O
          y_OH1 = y_H1 - y_O
          z_OH1 = z_H1 - z_O
          x_OH1 = x_OH1 - length_a*nint(x_OH1/length_a)          
          y_OH1 = y_OH1 - length_b*nint(y_OH1/length_b)
          z_OH1 = z_OH1 - length_c*nint(z_OH1/length_c)

          x_OH2 = x_H2 - x_O
          y_OH2 = y_H2 - y_O
          z_OH2 = z_H2 - z_O
          x_OH2 = x_OH2 - length_a*nint(x_OH2/length_a)          
          y_OH2 = y_OH2 - length_b*nint(y_OH2/length_b)
          z_OH2 = z_OH2 - length_c*nint(z_OH2/length_c)

          pole_x = x_OH1 + x_OH2
          pole_y = y_OH1 + y_OH2
          pole_z = z_OH1 + z_OH2
          pole_mono = sqrt(pole_x**2 + pole_y**2 + pole_z**2)
          pole_mono = q_H * convert * pole_mono

          pole_x = o_stc(t_step,o1,1) + h_stc(t_step,h1,1) + h_stc(t_step,h2,1)
          pole_y = o_stc(t_step,o1,2) + h_stc(t_step,h1,2) + h_stc(t_step,h2,2)
          pole_z = o_stc(t_step,o1,3) + h_stc(t_step,h1,3) + h_stc(t_step,h2,3)
          pole_stc = sqrt(pole_x**2 + pole_y**2 + pole_z**2)
          pole_stc = convert * pole_stc 

          cell_dpl_x = cell_dpl_x + pole_x
          cell_dpl_y = cell_dpl_y + pole_y
          cell_dpl_z = cell_dpl_z + pole_z

          pole_x = o_ind(t_step,o1,1) + h_ind(t_step,h1,1) + h_ind(t_step,h2,1)
          pole_y = o_ind(t_step,o1,2) + h_ind(t_step,h1,2) + h_ind(t_step,h2,2)
          pole_z = o_ind(t_step,o1,3) + h_ind(t_step,h1,3) + h_ind(t_step,h2,3)
          pole_ind = sqrt(pole_x**2 + pole_y**2 + pole_z**2)
          pole_ind = convert * pole_ind 

          cell_dpl_x = cell_dpl_x + pole_x
          cell_dpl_y = cell_dpl_y + pole_y
          cell_dpl_z = cell_dpl_z + pole_z

          pole_total = pole_mono + pole_stc + pole_ind    

          mono_avg = mono_avg + pole_mono
           stc_avg =  stc_avg + pole_stc
           ind_avg =  ind_avg + pole_ind
          pole_avg = pole_avg + pole_total

          i_bin = nint (pole_mono / polar_delta) + 1
          if (i_bin .le. n_bins) then
            monopl_dist(i_bin) = monopl_dist(i_bin) + 1 
          end if  

          i_bin = nint (pole_stc / polar_delta) + 1
          if (i_bin .le. n_bins) then
            static_dist(i_bin) = static_dist(i_bin) + 1 
          end if    

          i_bin = nint (pole_ind / polar_delta) + 1
          if (i_bin .le. n_bins) then
            induce_dist(i_bin) = induce_dist(i_bin) + 1 
          end if         

          i_bin = nint (pole_total / polar_delta) + 1
          if (i_bin .le. n_bins) then
            plrztn_dist(i_bin) = plrztn_dist(i_bin) + 1 
          end if

        end do inner_water
        cell_dplsq_time(i_snap) = cell_dplsq_time(i_snap) + &
                       cell_dpl_x**2 + cell_dpl_y**2 + cell_dpl_z**2
      end do outer
        cell_dplsq_time(i_snap) = cell_dplsq_time(i_snap) / frames_here
        mono_avg_time(i_snap) = mono_avg / (frames_here * n_water)
         stc_avg_time(i_snap) =  stc_avg / (frames_here * n_water)
         ind_avg_time(i_snap) =  ind_avg / (frames_here * n_water)
        pole_avg_time(i_snap) = pole_avg / (frames_here * n_water)
    end subroutine molecular_dipoles


    subroutine norm_distributions (static_dist, induce_dist, monopl_dist, plrztn_dist)

      implicit none
      real(dp) :: dist_area 
      real(dp) :: static_dist(:), induce_dist(:), monopl_dist(:), plrztn_dist(:)
      integer  :: i_bin

      dist_area = 0.0d0
      do i_bin = 1, n_bins
        dist_area = dist_area + static_dist(i_bin)
      end do
      dist_area = dist_area * polar_delta
      do i_bin = 1, n_bins
        static_dist(i_bin) = static_dist(i_bin) / dist_area
      end do      

      dist_area = 0.0d0
      do i_bin = 1, n_bins
        dist_area = dist_area +  induce_dist(i_bin)
      end do
      dist_area = dist_area * polar_delta
      do i_bin = 1, n_bins
         induce_dist(i_bin) =  induce_dist(i_bin) / dist_area
      end do

      dist_area = 0.0d0
      do i_bin = 1, n_bins
        dist_area = dist_area + monopl_dist(i_bin)
      end do
      dist_area = dist_area * polar_delta
      do i_bin = 1, n_bins
        monopl_dist(i_bin) = monopl_dist(i_bin) / dist_area
      end do      

      dist_area = 0.0d0
      do i_bin = 1, n_bins
        dist_area = dist_area + plrztn_dist(i_bin)
      end do
      dist_area = dist_area * polar_delta
      do i_bin = 1, n_bins
        plrztn_dist(i_bin) = plrztn_dist(i_bin) / dist_area
      end do

    end subroutine norm_distributions

    subroutine save_distributions (static_dist, induce_dist, monopl_dist, &
               plrztn_dist, static_time, induce_time, monopl_time, plrztn_time)

      implicit none
      integer  :: i_bin
      real(dp) :: static_dist(:), induce_dist(:), monopl_dist(:), plrztn_dist(:)
      real(dp) :: static_time(:,:), induce_time(:,:), monopl_time(:,:), plrztn_time(:,:)

      do i_bin = 1, n_bins
        static_time(i_snap, i_bin) = static_dist(i_bin)
        induce_time(i_snap, i_bin) = induce_dist(i_bin)
        monopl_time(i_snap, i_bin) = monopl_dist(i_bin)
        plrztn_time(i_snap, i_bin) = plrztn_dist(i_bin)
      end do

    end subroutine save_distributions

    subroutine distribution_convergence(static_time, induce_time, monopl_time, plrztn_time, stc_avg_conv,      &
               ind_avg_conv, mono_avg_conv, pole_avg_conv, static_conv, induce_conv, monopl_conv, plrztn_conv, &
               cell_dplsq_time, cell_dplsq_conv, stc_avg_time, ind_avg_time, mono_avg_time, pole_avg_time)

      implicit none
      integer  :: i_snap, i_bin
      real(dp) :: static_time(:,:), induce_time(:,:), monopl_time(:,:), plrztn_time(:,:)
      real(dp) :: static_conv(:), induce_conv(:), monopl_conv(:), plrztn_conv(:)
      real(dp) :: stc_avg_conv(:), ind_avg_conv(:), mono_avg_conv(:), pole_avg_conv(:)
      real(dp) :: stc_avg_time(:), ind_avg_time(:), mono_avg_time(:), pole_avg_time(:)
      real(dp) :: cell_dplsq_time(:), cell_dplsq_conv(:)
      real(dp) :: squ_1, squ_2
      real(dp) :: time_now


      open(unit=11, file='static_dist.dat', status='replace')
      open(unit=12, file='induce_dist.dat', status='replace')
      open(unit=13, file='monopl_dist.dat', status='replace')
      open(unit=14, file='plrztn_dist.dat', status='replace')

      open(unit=21, file='static_dist_conv.dat', status='replace')
      open(unit=22, file='induce_dist_conv.dat', status='replace')
      open(unit=23, file='monopl_dist_conv.dat', status='replace')      
      open(unit=24, file='plrztn_dist_conv.dat', status='replace') 

      open(unit=31, file='static_avg_conv.dat', status='replace')
      open(unit=32, file='induce_avg_conv.dat', status='replace')
      open(unit=33, file='monopl_avg_conv.dat', status='replace')      
      open(unit=34, file='plrztn_avg_conv.dat', status='replace')  
      open(unit=35, file='static_avg_time.dat', status='replace')
      open(unit=36, file='induce_avg_time.dat', status='replace')
      open(unit=37, file='monopl_avg_time.dat', status='replace')      
      open(unit=38, file='plrztn_avg_time.dat', status='replace')  


      open(unit=40, file='celldpl_sq_avg_conv.dat', status='replace')      
      open(unit=41, file='celldpl_sq_avg_time.dat', status='replace')      

      do i_snap = 1, n_snaps
        static_conv(i_snap) = 0.0d0
        induce_conv(i_snap) = 0.0d0
        monopl_conv(i_snap) = 0.0d0
        plrztn_conv(i_snap) = 0.0d0
        do i_bin = 1, n_bins
          squ_1 = static_time(n_snaps, i_bin) - static_time(i_snap, i_bin)
          static_conv(i_snap) = static_conv(i_snap) + abs((squ_1))

          squ_1 = induce_time(n_snaps, i_bin) - induce_time(i_snap, i_bin)
          induce_conv(i_snap) = induce_conv(i_snap) + abs((squ_1))

          squ_1 = monopl_time(n_snaps, i_bin) - monopl_time(i_snap, i_bin)
          monopl_conv(i_snap) = monopl_conv(i_snap) + abs((squ_1))

          squ_1 = plrztn_time(n_snaps, i_bin) - plrztn_time(i_snap, i_bin)
          plrztn_conv(i_snap) = plrztn_conv(i_snap) + abs((squ_1))
        end do
        static_conv(i_snap) = static_conv(i_snap) / real(n_bins)
        induce_conv(i_snap) = induce_conv(i_snap) / real(n_bins)
        monopl_conv(i_snap) = monopl_conv(i_snap) / real(n_bins)
        plrztn_conv(i_snap) = plrztn_conv(i_snap) / real(n_bins)
        
        squ_1 = stc_avg_time   (n_snaps) - stc_avg_time   (i_snap)
        stc_avg_conv   (i_snap) = abs (squ_1)

        squ_1 = ind_avg_time   (n_snaps) - ind_avg_time   (i_snap)
        ind_avg_conv   (i_snap) = abs (squ_1)

        squ_1 = mono_avg_time  (n_snaps) - mono_avg_time  (i_snap)
        mono_avg_conv  (i_snap) = abs (squ_1)

        squ_1 = pole_avg_time  (n_snaps) - pole_avg_time  (i_snap)
        pole_avg_conv  (i_snap) = abs (squ_1)

        squ_1 = cell_dplsq_time(n_snaps) - cell_dplsq_time(i_snap)
        cell_dplsq_conv(i_snap) = abs (squ_1)

      end do

      do i_bin = 1, n_bins
        write(11,*) polar_axis(i_bin), static_dist(i_bin)   
        write(12,*) polar_axis(i_bin), induce_dist(i_bin)
        write(13,*) polar_axis(i_bin), monopl_dist(i_bin)
        write(14,*) polar_axis(i_bin), plrztn_dist(i_bin)
      end do 

      do i_snap = 1, n_snaps
        time_now = delta_time*i_snap*&
              (real(n_frames_total) / (real(n_snaps)/real(i_snap) ) )

        write(21,*) time_now, static_conv(i_snap)
        write(22,*) time_now, induce_conv(i_snap)
        write(23,*) time_now, monopl_conv(i_snap)
        write(24,*) time_now, plrztn_conv(i_snap)        
      end do

      do i_snap = 1, n_snaps
        time_now = delta_time*(real(n_frames_total) / (real(n_snaps)/real(i_snap) ) )

        write(31,*) time_now,  stc_avg_conv(i_snap)
        write(32,*) time_now,  ind_avg_conv(i_snap)
        write(33,*) time_now, mono_avg_conv(i_snap)
        write(34,*) time_now, pole_avg_conv(i_snap)

        write(35,*) time_now,  stc_avg_time(i_snap)
        write(36,*) time_now,  ind_avg_time(i_snap)
        write(37,*) time_now, mono_avg_time(i_snap)
        write(38,*) time_now, pole_avg_time(i_snap)

        write(40,*) time_now, cell_dplsq_conv(i_snap)
        write(41,*) time_now, cell_dplsq_time(i_snap)
      end do

      close(unit=11)
      close(unit=12)
      close(unit=13)
      close(unit=14)

      close(unit=21)
      close(unit=22)
      close(unit=23)
      close(unit=24)

      close(unit=31)
      close(unit=32)
      close(unit=33)
      close(unit=34)
      close(unit=35)
      close(unit=36)

    end subroutine distribution_convergence

end program water_polar_distributions