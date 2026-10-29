! ============================================================================
! 1D-3V electromagnetic PIC: educational approximation to Klimontovich-Maxwell
! Normalized units: electron |q|=1, m=1, epsilon0=mu0=1, c=1.
! x is periodic; particles have (vx,vy,vz); fields have Ex,Ey,Ez,Bx,By,Bz.
! Longitudinal Ex is reconstructed from Gauss's law; transverse fields use a
! staggered-grid (Yee-style) FDTD update. This is a teaching baseline, not a
! production/research PIC code. See README notes on limitations.
! ============================================================================
program klimontovich_maxwell_1d3v_pic
  use iso_fortran_env, only: real64
  implicit none
  integer, parameter :: dp=real64
  integer, parameter :: nx=128, np=20000, nsteps=2000, ndiag=10
  real(dp), parameter :: pi=acos(-1.0_dp), L=2.0_dp*pi
  real(dp), parameter :: dx=L/real(nx,dp), c=1.0_dp, eps0=1.0_dp, mu0=1.0_dp
  real(dp), parameter :: q=-1.0_dp, mass=1.0_dp, n0=1.0_dp
  real(dp), parameter :: dt=0.005_dp, vth=0.05_dp, amp=1.0e-3_dp
  real(dp), parameter :: macro_weight=n0*L/real(np,dp)
  real(dp) :: x(np), vx(np), vy(np), vz(np)
  real(dp) :: rho(nx), ne(nx), jx(nx), jy(nx), jz(nx)
  real(dp) :: ex(nx), ey(nx), ez(nx)
  real(dp) :: by(nx), bz(nx), by_new(nx), bz_new(nx)
  real(dp) :: t, particle_energy, field_energy, total_energy
  real(dp) :: gauss_rms, xmid
  integer :: i, p, step, unit_diag, unit_field, unit_part

  call random_seed()
  call initialize_particles()
  call deposit_charge()
  call solve_longitudinal_e()

  ! A small right-going transverse electromagnetic wave is the initial perturbation.
  do i=1,nx
     xmid=(real(i,dp)-0.5_dp)*dx
     ey(i)=amp*sin(xmid)
     ez(i)=0.0_dp
     ex(i)=ex(i) ! Ex came from Gauss's law
     by(i)=0.0_dp
     bz(i)=amp*sin((real(i,dp)-0.5_dp)*dx)/c
  end do
  ! B arrays are interpreted at cell faces x=(i-1/2)dx; the above phase is consistent.

  open(newunit=unit_diag,file='diagnostics.csv',status='replace',action='write')
  write(unit_diag,'(a)') 'step,time,particle_energy,electric_energy,magnetic_energy,total_energy,gauss_rms'
  open(newunit=unit_field,file='fields_final.dat',status='replace',action='write')
  open(newunit=unit_part,file='particles_final.dat',status='replace',action='write')

  do step=0,nsteps
     t=real(step,dp)*dt
     call deposit_charge_and_current()
     call solve_longitudinal_e()
     call calculate_diagnostics(particle_energy,field_energy,gauss_rms)
     total_energy=particle_energy+field_energy
     if (mod(step,ndiag)==0 .or. step==nsteps) then
        write(unit_diag,'(i0,6(a,es22.13))') step,',',t,',',particle_energy,',', &
             electric_energy_value(),',',magnetic_energy_value(),',',total_energy,',',gauss_rms
     end if
     if (step==nsteps) exit

     ! Transverse Maxwell update: Faraday, fields staggered in x.
     do i=1,nx
        p=mod(i,nx)+1
        by_new(i)=by(i)+dt*(ez(p)-ez(i))/dx
        bz_new(i)=bz(i)-dt*(ey(p)-ey(i))/dx
     end do
     ! Deposit currents at integer grid points. Longitudinal Ex is constrained by Gauss.
     ! Transverse Ampere update uses the old transverse fields and deposited Jy,Jz.
     do i=1,nx
        p=mod(i-2+nx,nx)+1
        ey(i)=ey(i)-dt*(c*c)*(bz_new(i)-bz_new(p))/dx-dt*jy(i)/eps0
        ez(i)=ez(i)+dt*(c*c)*(by_new(i)-by_new(p))/dx-dt*jz(i)/eps0
     end do
     by=by_new
     bz=bz_new

     ! Non-relativistic Boris pusher (valid while |v| << c).
     call push_particles()
  end do

  do i=1,nx
     write(unit_field,'(7(es20.11,1x))') (real(i,dp)-0.5_dp)*dx,ex(i),ey(i),ez(i),0.0_dp,by(i),bz(i)
  end do
  do p=1,np
     write(unit_part,'(4(es20.11,1x))') x(p),vx(p),vy(p),vz(p)
  end do
  close(unit_diag); close(unit_field); close(unit_part)
  print '(a,i0)', 'Completed steps: ',nsteps
  print '(a,i0)', 'Particles: ',np
  print '(a)', 'Wrote diagnostics.csv, fields_final.dat, particles_final.dat'

contains

  subroutine initialize_particles()
    integer :: k
    real(dp) :: u1,u2,u3,u4
    do k=1,np
       call random_number(x(k)); x(k)=x(k)*L
       call random_number(u1); call random_number(u2)
       call random_number(u3); call random_number(u4)
       u1=max(u1,1.0e-14_dp); u3=max(u3,1.0e-14_dp)
       vx(k)=vth*sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
       vy(k)=vth*sqrt(-2.0_dp*log(u1))*sin(2.0_dp*pi*u2)
       vz(k)=vth*sqrt(-2.0_dp*log(u3))*cos(2.0_dp*pi*u4)
    end do
  end subroutine initialize_particles

  subroutine deposit_charge()
    integer :: k, il, ir
    real(dp) :: s, frac
    ne=0.0_dp
    do k=1,np
       s=x(k)/dx
       il=int(floor(s))+1
       frac=s-floor(s)
       ir=mod(il,nx)+1
       ne(il)=ne(il)+macro_weight*(1.0_dp-frac)/dx
       ne(ir)=ne(ir)+macro_weight*frac/dx
    end do
    rho=n0+q*ne
  end subroutine deposit_charge

  subroutine deposit_charge_and_current()
    integer :: k, il, ir
    real(dp) :: s, frac
    ne=0.0_dp; jx=0.0_dp; jy=0.0_dp; jz=0.0_dp
    do k=1,np
       s=x(k)/dx
       il=int(floor(s))+1
       frac=s-floor(s)
       ir=mod(il,nx)+1
       ne(il)=ne(il)+macro_weight*(1.0_dp-frac)/dx
       ne(ir)=ne(ir)+macro_weight*frac/dx
       jx(il)=jx(il)+q*macro_weight*vx(k)*(1.0_dp-frac)/dx
       jx(ir)=jx(ir)+q*macro_weight*vx(k)*frac/dx
       jy(il)=jy(il)+q*macro_weight*vy(k)*(1.0_dp-frac)/dx
       jy(ir)=jy(ir)+q*macro_weight*vy(k)*frac/dx
       jz(il)=jz(il)+q*macro_weight*vz(k)*(1.0_dp-frac)/dx
       jz(ir)=jz(ir)+q*macro_weight*vz(k)*frac/dx
    end do
    rho=n0+q*ne
  end subroutine deposit_charge_and_current

  subroutine solve_longitudinal_e()
    integer :: k
    real(dp) :: mean_rho
    ! Integrate dEx/dx=rho/eps0 on the periodic domain and remove the mean Ex.
    ! This enforces a discrete Gauss constraint approximately; it is not a
    ! charge-conserving current-deposition scheme.
    mean_rho=sum(rho)/real(nx,dp)
    ex(1)=0.0_dp
    do k=2,nx
       ex(k)=ex(k-1)+(rho(k-1)-mean_rho)*dx/eps0
    end do
    ex=ex-sum(ex)/real(nx,dp)
  end subroutine solve_longitudinal_e

  subroutine push_particles()
    integer :: k, il, ir, left
    real(dp) :: s, frac, epx, epy, epz, bpx, bpy, bpz
    real(dp) :: vminus(3), tvec(3), svec(3), vprime(3), vplus(3), vnew(3)
    real(dp) :: qm, tx,ty,tz
    qm=q/mass
    do k=1,np
       s=x(k)/dx
       il=int(floor(s))+1; frac=s-floor(s); ir=mod(il,nx)+1
       left=mod(il-2+nx,nx)+1
       epx=ex(il)*(1.0_dp-frac)+ex(ir)*frac
       epy=ey(il)*(1.0_dp-frac)+ey(ir)*frac
       epz=ez(il)*(1.0_dp-frac)+ez(ir)*frac
       ! B is face-centered: average adjacent faces to the particle location.
       bpx=0.0_dp
       bpy=by(left)*(1.0_dp-frac)+by(il)*frac
       bpz=bz(left)*(1.0_dp-frac)+bz(il)*frac
       vminus=[vx(k),vy(k),vz(k)]+0.5_dp*qm*dt*[epx,epy,epz]
       tvec=0.5_dp*qm*dt*[bpx,bpy,bpz]
       tx=tvec(1); ty=tvec(2); tz=tvec(3)
       svec=2.0_dp*tvec/(1.0_dp+sum(tvec*tvec))
       vprime=vminus+[vminus(2)*tz-vminus(3)*ty, &
                      vminus(3)*tx-vminus(1)*tz, &
                      vminus(1)*ty-vminus(2)*tx]
       vplus=vminus+[vprime(2)*svec(3)-vprime(3)*svec(2), &
                     vprime(3)*svec(1)-vprime(1)*svec(3), &
                     vprime(1)*svec(2)-vprime(2)*svec(1)]
       vnew=vplus+0.5_dp*qm*dt*[epx,epy,epz]
       vx(k)=vnew(1); vy(k)=vnew(2); vz(k)=vnew(3)
       x(k)=modulo(x(k)+dt*vx(k),L)
    end do
  end subroutine push_particles

  subroutine calculate_diagnostics(pe,fe,grms)
    real(dp), intent(out) :: pe,fe,grms
    real(dp) :: epart, eg, eb, div
    integer :: k, j, jp
    epart=0.0_dp
    do k=1,np
       epart=epart+0.5_dp*macro_weight*mass*(vx(k)**2+vy(k)**2+vz(k)**2)
    end do
    eg=0.5_dp*eps0*dx*sum(ex*ex+ey*ey+ez*ez)
    eb=0.5_dp/mu0*dx*sum(by*by+bz*bz)
    pe=epart; fe=eg+eb
    grms=0.0_dp
    do j=1,nx
       jp=mod(j,nx)+1
       div=(ex(jp)-ex(j))/dx-rho(j)/eps0
       grms=grms+div*div
    end do
    grms=sqrt(grms/real(nx,dp))
  end subroutine calculate_diagnostics

  function electric_energy_value() result(val)
    real(dp) :: val
    val=0.5_dp*eps0*dx*sum(ex*ex+ey*ey+ez*ez)
  end function electric_energy_value

  function magnetic_energy_value() result(val)
    real(dp) :: val
    val=0.5_dp/mu0*dx*sum(by*by+bz*bz)
  end function magnetic_energy_value

end program klimontovich_maxwell_1d3v_pic
