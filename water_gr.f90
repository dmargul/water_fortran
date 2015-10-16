program water_radial_distributions
  implicit none
      
    integer,     parameter :: n_frames_total  = 500
    integer,     parameter :: n_water         = 512
    integer,     parameter :: n_ox            = n_water 
    integer,     parameter :: n_hy            = n_ox * 2
    integer,     parameter :: n_atoms         = n_ox + n_hy
    integer,     parameter :: oo_pairs        = n_ox * (n_ox - 1) / 2
    integer,     parameter :: oh_pairs        = n_ox * n_hy
    integer,     parameter :: hh_pairs        = n_hy * (n_hy - 1) / 2
    integer,     parameter :: total_oo_pairs  = n_frames_total * oo_pairs
    integer,     parameter :: total_oh_pairs  = n_frames_total * oh_pairs
    integer,     parameter :: total_hh_pairs  = n_frames_total * hh_pairs
    integer,     parameter :: n_bins_gr       = 200
    integer,     parameter :: n_snaps         = 5
    integer,     parameter :: dp              = kind(0.0d0)
    integer                :: i_snap
    integer                :: n_frames
    integer                :: frames_here
    real(dp),    parameter :: delta_time      = 0.150d0 ! picoseconds between frames
    real(dp),    parameter :: pi              = acos(-1.0d0)
    real(dp),    parameter :: distance_min_gr = 0.0d0  ! Angstrom
    real(dp),    parameter :: distance_max_gr = 12.4d0 ! Angstrom
    real(dp),    parameter :: gr_dist_delta   = (distance_max_gr - distance_min_gr) / real(n_bins_gr)
    real(dp),    parameter :: length_a        = 24.832132 ! Angstrom
    real(dp),    parameter :: length_b        = 24.832132 ! Angstrom
    real(dp),    parameter :: length_c        = 24.832132 ! Angstrom 
    real(dp),    parameter :: box_volume      = length_a * length_b * length_c
    real(dp),    parameter :: o_mass          = 15.995 ! amu
    real(dp),    parameter :: h_mass          = 1.008  ! amu
    real(dp),    parameter :: h2o_mass        = o_mass  + h_mass  + h_mass
    real(dp),    parameter :: system_mass     = n_water * h2o_mass
    real(dp)               :: o_pos (n_frames_total, n_ox, 3)
    real(dp)               :: h_pos (n_frames_total, n_hy, 3)
    real(dp)               :: distance_oo 
    real(dp)               :: distance_oh 
    real(dp)               :: distance_hh 
    real(dp)               :: oo_counted (n_bins_gr)
    real(dp)               :: oh_counted (n_bins_gr)
    real(dp)               :: hh_counted (n_bins_gr)
    real(dp)               :: distance_gr (n_bins_gr)
    real(dp)               :: oo_gr (n_bins_gr)
    real(dp)               :: oh_gr (n_bins_gr)
    real(dp)               :: hh_gr (n_bins_gr)
    real(dp)               :: oo_gr_time (n_snaps, n_bins_gr)
    real(dp)               :: oh_gr_time (n_snaps, n_bins_gr)
    real(dp)               :: hh_gr_time (n_snaps, n_bins_gr)
    real(dp)               :: oo_conv (n_snaps)
    real(dp)               :: oh_conv (n_snaps)
    real(dp)               :: hh_conv (n_snaps)

    call read_water_xyz (o_pos, h_pos)
    call center_mass_removal (o_pos, h_pos)
    do i_snap = 1, n_snaps
      frames_here = nint(real(n_frames_total) / (real(n_snaps)/real(i_snap)) )
      call make_distributions (o_pos, h_pos, oo_counted, oh_counted, hh_counted)
      call norm_distributions (oo_counted, oh_counted, hh_counted, oo_gr, oh_gr, hh_gr, distance_gr)
      call save_distributions (oo_gr, oh_gr, hh_gr, oo_gr_time, oh_gr_time, hh_gr_time)
    end do

    call distribution_convergence (oo_gr_time, oh_gr_time, hh_gr_time, oo_conv, oh_conv, hh_conv)

  contains

    subroutine read_water_xyz(o_pos,h_pos)

      implicit none
      integer :: t_step, o, h1, h2
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

    subroutine center_mass_removal(o_pos,h_pos)

      implicit none
      integer :: t_step, o1, o2, h, h1, h2
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

    subroutine make_distributions(o_pos,h_pos,oo_counted,oh_counted,hh_counted)

      implicit none
      integer  :: t_step, o1, o2, h, h1, h2, i_bin
      real(dp) :: o_pos(:,:,:), h_pos(:,:,:)
      real(dp) :: x_1, x_2, y_1, y_2, z_1, z_2
      real(dp) :: disp_x2, disp_y2, disp_z2, disp_r2
      real(dp) :: disp_x, disp_y, disp_z, disp_r
      real(dp) :: oo_counted(:), oh_counted(:), hh_counted(:)
      real(dp) :: distance_oo, distance_oh, distance_hh
      real(dp) :: dist


      init: do i_bin = 1, n_bins_gr
        dist = distance_min_gr + i_bin*gr_dist_delta  
        distance_gr(i_bin) = dist
        oo_counted(i_bin) = 0
        oh_counted(i_bin) = 0
        hh_counted(i_bin) = 0
      end do init

      outer: do t_step = 1, frames_here

        inner_o: do o1 = 1, n_water
          x_1 = o_pos(t_step, o1, 1)
          y_1 = o_pos(t_step, o1, 2)
          z_1 = o_pos(t_step, o1, 3)

          inner_o2: do o2 = o1+1, n_water
            x_2 = o_pos(t_step, o2, 1)
            y_2 = o_pos(t_step, o2, 2)
            z_2 = o_pos(t_step, o2, 3)
            disp_x2  = (x_2 - x_1)**2
            disp_x   = sqrt(disp_x2) 
            disp_x   = disp_x - length_a*nint(disp_x/length_a)
            disp_y2  = (y_2 - y_1)**2    
            disp_y   = sqrt(disp_y2)        
            disp_y   = disp_y - length_b*nint(disp_y/length_b)
            disp_z2  = (z_2 - z_1)**2 
            disp_z   = sqrt(disp_z2)           
            disp_z   = disp_z - length_c*nint(disp_z/length_c)
            disp_r2  = (disp_x**2 + disp_y**2 + disp_z**2)
            distance_oo = sqrt(disp_r2)
            i_bin = nint(distance_oo / gr_dist_delta) + 1
            if (i_bin .le. n_bins_gr) then
              oo_counted(i_bin) = oo_counted(i_bin) + 1
            end if 
          end do inner_o2

          inner_h: do h = 1, n_hy
            x_2 = h_pos(t_step, h, 1)
            y_2 = h_pos(t_step, h, 2)
            z_2 = h_pos(t_step, h, 3)
            disp_x2  = (x_2 - x_1)**2
            disp_x   = sqrt(disp_x2) 
            disp_x   = disp_x - length_a*nint(disp_x/length_a)
            disp_y2  = (y_2 - y_1)**2    
            disp_y   = sqrt(disp_y2)        
            disp_y   = disp_y - length_b*nint(disp_y/length_b)
            disp_z2  = (z_2 - z_1)**2 
            disp_z   = sqrt(disp_z2)           
            disp_z   = disp_z - length_c*nint(disp_z/length_c)
            disp_r2  = (disp_x**2 + disp_y**2 + disp_z**2)
            distance_oh = sqrt(disp_r2)
            i_bin = nint(distance_oh / gr_dist_delta) + 1
            if (i_bin .le. n_bins_gr) then
              oh_counted(i_bin) = oh_counted(i_bin) + 1 
            end if           
          end do inner_h

        end do inner_o

        inner_h1: do h1 = 1, n_hy
          x_1 = h_pos(t_step, h1, 1)
          y_1 = h_pos(t_step, h1, 2)
          z_1 = h_pos(t_step, h1, 3)
          inner_h2: do h2 = h1+1, n_hy
            x_2 = h_pos(t_step, h2, 1)
            y_2 = h_pos(t_step, h2, 2)
            z_2 = h_pos(t_step, h2, 3)
            disp_x2  = (x_2 - x_1)**2
            disp_x   = sqrt(disp_x2) 
            disp_x   = disp_x - length_a*nint(disp_x/length_a)
            disp_y2  = (y_2 - y_1)**2    
            disp_y   = sqrt(disp_y2)        
            disp_y   = disp_y - length_b*nint(disp_y/length_b)
            disp_z2  = (z_2 - z_1)**2 
            disp_z   = sqrt(disp_z2)           
            disp_z   = disp_z - length_c*nint(disp_z/length_c)
            disp_r2  = (disp_x**2 + disp_y**2 + disp_z**2)
            distance_hh = sqrt(disp_r2)
            i_bin = nint(distance_hh / gr_dist_delta) + 1
            if (i_bin .le. n_bins_gr) then
              hh_counted(i_bin) = hh_counted(i_bin) + 1 
            end if           
          end do inner_h2
        end do inner_h1

      end do outer

    end subroutine make_distributions


    subroutine norm_distributions(oo_counted, oh_counted, hh_counted, oo_gr, oh_gr, hh_gr, distance_gr)

      implicit none
      real(dp)     :: gr_sphere_const 
      real(dp)     :: oo_counted(:), oh_counted(:), hh_counted(:)
      real(dp)     :: distance_gr(:), oo_gr(:), oh_gr(:), hh_gr(:)
      real(dp)     :: r_high, r_low
      real(dp)     :: n_theory, n_theory_oo, n_theory_oh, n_theory_hh
      integer  :: i_bin

      gr_sphere_const = 4.0d0 * pi * frames_here  / 3.0d0 

      do i_bin = 1, n_bins_gr
        r_high = real(i_bin) * gr_dist_delta
        r_low = r_high - gr_dist_delta 
        n_theory = gr_sphere_const * (r_high**3 - r_low**3)
        n_theory_oo = n_theory * oo_pairs / box_volume
        n_theory_oh = n_theory * oh_pairs / box_volume
        n_theory_hh = n_theory * hh_pairs / box_volume
        oo_gr(i_bin) = oo_counted(i_bin) / n_theory_oo
        oh_gr(i_bin) = oh_counted(i_bin) / n_theory_oh
        hh_gr(i_bin) = hh_counted(i_bin) / n_theory_hh
      end do

    end subroutine norm_distributions

    subroutine save_distributions(oo_gr, oh_gr, hh_gr, oo_gr_time, oh_gr_time, hh_gr_time)

      implicit none
      integer :: i_bin
      real(dp)    :: oo_gr(:), oh_gr(:), hh_gr(:)
      real(dp)    :: oo_gr_time(:,:), oh_gr_time(:,:), hh_gr_time(:,:)

      do i_bin = 1, n_bins_gr
        oo_gr_time(i_snap, i_bin) = oo_gr(i_bin)
        oh_gr_time(i_snap, i_bin) = oh_gr(i_bin)
        hh_gr_time(i_snap, i_bin) = hh_gr(i_bin)
      end do

    end subroutine save_distributions

    subroutine distribution_convergence(oo_gr_time, oh_gr_time, hh_gr_time, oo_conv, oh_conv, hh_conv)

      implicit none
      integer :: i_snap, i_bin
      real(dp)    :: oo_gr_time(:,:), oh_gr_time(:,:), hh_gr_time(:,:)
      real(dp)    :: oo_conv(:), oh_conv(:), hh_conv(:)
      real(dp)    :: squ_1, squ_2

      open(unit=11, file='oo_gr.dat', status='replace')
      open(unit=12, file='oh_gr.dat', status='replace')
      open(unit=13, file='hh_gr.dat', status='replace')

      open(unit=21, file='oo_gr_conv.dat', status='replace')
      open(unit=22, file='oh_gr_conv.dat', status='replace')
      open(unit=23, file='hh_gr_conv.dat', status='replace')      

      do i_snap = 1, n_snaps
        oo_conv(i_snap) = 0.0d0
        oh_conv(i_snap) = 0.0d0
        hh_conv(i_snap) = 0.0d0
        do i_bin = 1, n_bins_gr
          squ_1 = oo_gr_time(n_snaps, i_bin) - oo_gr_time(i_snap, i_bin)
          oo_conv(i_snap) = oo_conv(i_snap) + abs((squ_1))

          squ_1 = oh_gr_time(n_snaps, i_bin) - oh_gr_time(i_snap, i_bin)
          oo_conv(i_snap) = oo_conv(i_snap) + abs((squ_1))

          squ_1 = hh_gr_time(n_snaps, i_bin) - hh_gr_time(i_snap, i_bin)
          hh_conv(i_snap) = hh_conv(i_snap) + abs((squ_1))
        end do
          oo_conv(i_snap) = oo_conv(i_snap) / real(n_bins_gr)
          oh_conv(i_snap) = oh_conv(i_snap) / real(n_bins_gr)
          hh_conv(i_snap) = hh_conv(i_snap) / real(n_bins_gr)
      end do

      do i_bin = 1, n_bins_gr
        write(11,*) distance_gr(i_bin), oo_gr_time(n_snaps, i_bin)
        write(12,*) distance_gr(i_bin), oh_gr_time(n_snaps, i_bin)
        write(13,*) distance_gr(i_bin), hh_gr_time(n_snaps, i_bin)
      end do 

      do i_snap = 1, n_snaps
        write(21,*) i_snap*(real(n_frames_total) / (real(n_snaps)/real(i_snap) ) ), oo_conv(i_snap)
        write(22,*) i_snap*(real(n_frames_total) / (real(n_snaps)/real(i_snap) ) ), oh_conv(i_snap)
        write(23,*) i_snap*(real(n_frames_total) / (real(n_snaps)/real(i_snap) ) ), hh_conv(i_snap)
      end do

      close(unit=11)
      close(unit=12)
      close(unit=13)   

      close(unit=21)
      close(unit=22)
      close(unit=23)


    end subroutine distribution_convergence

end program water_radial_distributions