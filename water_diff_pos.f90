program water_einstein_diffusion
  implicit none
      
    integer,    parameter :: n_frames_total  = 2200
    integer,    parameter :: frames_msd      = 500 ! maximum time between configurations
    integer,    parameter :: n_water         = 512
    integer,    parameter :: n_ox            = n_water 
    integer,    parameter :: n_hy            = n_ox * 2
    integer,    parameter :: n_atoms         = n_ox + n_hy
    integer,    parameter :: n_snaps         = 10
    integer,    parameter :: dp              = kind(0.0d0)
    integer               :: i_snap
    integer               :: frames_here
    real(dp),   parameter :: delta_time      = 0.050d0 ! picoseconds between frames
    real(dp),   parameter :: pi              = acos(-1.0d0)
    real(dp),   parameter :: length_a        = 24.832132 ! Angstrom
    real(dp),   parameter :: length_b        = 24.832132 ! Angstrom
    real(dp),   parameter :: length_c        = 24.832132 ! Angstrom 
    real(dp),   parameter :: box_volume      = length_a * length_b * length_c
    real(dp),   parameter :: o_mass          = 15.995 ! amu
    real(dp),   parameter :: h_mass          = 1.008  ! amu
    real(dp),   parameter :: h2o_mass        = o_mass  + h_mass  + h_mass
    real(dp),   parameter :: system_mass     = n_water * h2o_mass
    real(dp)              :: o_pos (n_frames_total, n_ox, 3)
    real(dp)              :: h_pos (n_frames_total, n_hy, 3)
    real(dp)              :: time_series(frames_msd)
    real(dp)              :: x_msd (frames_msd)
    real(dp)              :: y_msd (frames_msd)
    real(dp)              :: z_msd (frames_msd)
    real(dp)              ::   msd (frames_msd)
    real(dp)              :: x_msd_time (n_snaps, frames_msd)
    real(dp)              :: y_msd_time (n_snaps, frames_msd)
    real(dp)              :: z_msd_time (n_snaps, frames_msd)
    real(dp)              ::   msd_time (n_snaps, frames_msd)
    real(dp)              :: x_msd_conv (n_snaps)
    real(dp)              :: y_msd_conv (n_snaps)
    real(dp)              :: z_msd_conv (n_snaps)
    real(dp)              ::   msd_conv (n_snaps)


    call read_water_xyz (o_pos, h_pos)
    call center_mass_removal (o_pos, h_pos) 
    do i_snap = 1, n_snaps
      frames_here = nint(real(n_frames_total) / (real(n_snaps)/real(i_snap)) )
      if (frames_here.ge.frames_msd) then
        call msd_calculation (o_pos, h_pos, x_msd, y_msd, z_msd, msd)
        call save_msd (x_msd, y_msd, z_msd, msd, x_msd_time, y_msd_time, z_msd_time, msd_time)
      end if
    end do
    call convergence_msd (x_msd_time, y_msd_time, z_msd_time, msd_time, x_msd_conv, y_msd_conv, z_msd_conv, msd_conv, time_series)

  contains
!-------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------!
    subroutine read_water_xyz(o_pos,h_pos)

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
!-------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------!
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
!-------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------!
    subroutine msd_calculation (o_pos, h_pos, x_msd, y_msd, z_msd, msd)

      implicit none
      integer  :: t_step, o1, o2, h, h1, h2, i_bin
      integer  :: delta_frames, frame1, frame2
      real(dp) :: o_pos(:,:,:), h_pos(:,:,:)
      real(dp) :: x_msd(:), y_msd(:), z_msd(:), msd(:)
      real(dp) :: x_now, x_old, x_tmp
      real(dp) :: y_now, y_old, y_tmp
      real(dp) :: z_now, z_old, z_tmp
      real(dp) :: dist

      do delta_frames = 1, frames_msd
          x_tmp = 0.0d0
          y_tmp = 0.0d0
          z_tmp = 0.0d0
        do frame1 = 1, frames_here - delta_frames
          frame2 = frame1 + delta_frames
          do o1 = 1, n_ox
            x_now = o_pos(frame1, o1, 1)
            y_now = o_pos(frame1, o1, 2)
            z_now = o_pos(frame1, o1, 3)            
            x_old = o_pos(frame2, o1, 1)
            y_old = o_pos(frame2, o1, 2)
            z_old = o_pos(frame2, o1, 3)
            x_tmp = x_tmp + (x_now-x_old)**2
            y_tmp = y_tmp + (y_now-y_old)**2
            z_tmp = z_tmp + (z_now-z_old)**2
          end do
!          do h1 = 1, n_hy
!            x_now = h_pos(frame1, h1, 1)
!            y_now = h_pos(frame1, h1, 2)
!            z_now = h_pos(frame1, h1, 3)            
!            x_old = h_pos(frame2, h1, 1)
!            y_old = h_pos(frame2, h1, 2)
!            z_old = h_pos(frame2, h1, 3)
!            x_tmp = x_tmp + (x_now-x_old)**2
!            y_tmp = y_tmp + (y_now-y_old)**2
!            z_tmp = z_tmp + (z_now-z_old)**2
!          end do
        end do
        x_tmp = x_tmp / ((frames_here-(delta_frames+1))*n_ox)
        y_tmp = y_tmp / ((frames_here-(delta_frames+1))*n_ox)
        z_tmp = z_tmp / ((frames_here-(delta_frames+1))*n_ox)
        x_msd(delta_frames) = x_tmp
        y_msd(delta_frames) = y_tmp
        z_msd(delta_frames) = z_tmp
          msd(delta_frames) = x_tmp + y_tmp + z_tmp
      end do

    end subroutine msd_calculation
!-------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------!
    subroutine save_msd(x_msd, y_msd, z_msd, msd, x_msd_time, y_msd_time, z_msd_time, msd_time)

      implicit none
      integer  :: i_bin
      real(dp) :: x_msd(:), y_msd(:), z_msd(:), msd(:)
      real(dp) :: x_msd_time(:,:), y_msd_time(:,:), z_msd_time(:,:), msd_time(:,:)

      do i_bin = 1, frames_msd
        x_msd_time(i_snap, i_bin) = x_msd(i_bin)
        y_msd_time(i_snap, i_bin) = y_msd(i_bin)
        z_msd_time(i_snap, i_bin) = z_msd(i_bin)
          msd_time(i_snap, i_bin) =   msd(i_bin)
      end do

    end subroutine save_msd
!-------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------!
    subroutine convergence_msd(x_msd_time, y_msd_time, z_msd_time, msd_time, x_msd_conv, &
                               y_msd_conv, z_msd_conv, msd_conv, time_series)

      implicit none
      integer  :: i_snap, i_bin, delta_frames
      real(dp) :: x_msd_time(:,:), y_msd_time(:,:), z_msd_time(:,:), msd_time(:,:)
      real(dp) :: x_msd_conv(:), y_msd_conv(:), z_msd_conv(:), msd_conv(:), time_series(:)
      real(dp) :: squ_1, squ_2
      real(dp) :: time_now

      open(unit=11, file='x_msd.dat', status='replace')
      open(unit=12, file='y_msd.dat', status='replace')
      open(unit=13, file='z_msd.dat', status='replace')
      open(unit=14, file=  'msd.dat', status='replace')

      open(unit=21, file='x_msd_conv.dat', status='replace')
      open(unit=22, file='y_msd_conv.dat', status='replace')
      open(unit=23, file='z_msd_conv.dat', status='replace')      
      open(unit=24, file=  'msd_conv.dat', status='replace')      

      do i_bin = 1, frames_msd
        time_series(i_bin) = i_bin*delta_time
      end do
 
      do i_snap = 1, n_snaps
        frames_here = nint(real(n_frames_total) / (real(n_snaps)/real(i_snap)) )
        if (frames_here.ge.frames_msd) then
          x_msd_conv(i_snap) = 0.0d0
          y_msd_conv(i_snap) = 0.0d0
          z_msd_conv(i_snap) = 0.0d0
            msd_conv(i_snap) = 0.0d0
          do i_bin = 1, frames_msd
            squ_1 = x_msd_time(n_snaps, i_bin) - x_msd_time(i_snap, i_bin)
            x_msd_conv(i_snap) = x_msd_conv(i_snap) + abs((squ_1))

            squ_1 = y_msd_time(n_snaps, i_bin) - y_msd_time(i_snap, i_bin)
            x_msd_conv(i_snap) = x_msd_conv(i_snap) + abs((squ_1))

            squ_1 = z_msd_time(n_snaps, i_bin) - z_msd_time(i_snap, i_bin)
            z_msd_conv(i_snap) = z_msd_conv(i_snap) + abs((squ_1))

            squ_1 =   msd_time(n_snaps, i_bin) -   msd_time(i_snap, i_bin)
              msd_conv(i_snap) =   msd_conv(i_snap) + abs((squ_1))
          end do
            x_msd_conv(i_snap) = x_msd_conv(i_snap) / real(frames_msd)
            y_msd_conv(i_snap) = y_msd_conv(i_snap) / real(frames_msd)
            z_msd_conv(i_snap) = z_msd_conv(i_snap) / real(frames_msd)
              msd_conv(i_snap) =   msd_conv(i_snap) / real(frames_msd)
        else
          x_msd_conv(i_snap) = 0.0d0
          y_msd_conv(i_snap) = 0.0d0
          z_msd_conv(i_snap) = 0.0d0
            msd_conv(i_snap) = 0.0d0
        end if
      end do

      do i_bin = 1, frames_msd
        write(11,*) time_series(i_bin), x_msd_time(n_snaps, i_bin)
        write(12,*) time_series(i_bin), y_msd_time(n_snaps, i_bin)
        write(13,*) time_series(i_bin), z_msd_time(n_snaps, i_bin)
        write(14,*) time_series(i_bin),   msd_time(n_snaps, i_bin)
      end do 

      do i_snap = 1, n_snaps
        time_now = delta_time*(real(n_frames_total) / (real(n_snaps)/real(i_snap) ) )

        write(21,*) time_now, x_msd_conv(i_snap)
        write(22,*) time_now, y_msd_conv(i_snap)
        write(23,*) time_now, z_msd_conv(i_snap)
        write(24,*) time_now,   msd_conv(i_snap)
      end do

      close(unit=11)
      close(unit=12)
      close(unit=13)   
      close(unit=14)

      close(unit=21)
      close(unit=22)
      close(unit=23)
      close(unit=24)

    end subroutine convergence_msd



end program water_einstein_diffusion