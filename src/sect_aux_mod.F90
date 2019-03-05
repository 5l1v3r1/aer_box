module sect_aux_mod

  use precision_mod

  implicit none
  private

  public write_state

  contains

  subroutine write_state(n_bins,t_now,T_K,p_hPa,ndens,vvH2O,vvH2SO4,vvSO4,out_id,write_header,write_data,box_id)
  ! Write the current state to an output file (ASCII for now)
  integer,            intent(in) :: out_id
  integer,  optional, intent(in) :: t_now
  real(fp), optional, intent(in) :: T_K
  real(fp), optional, intent(in) :: p_hPa
  real(fp), optional, intent(in) :: ndens
  real(fp), optional, intent(in) :: vvH2O
  real(fp), optional, intent(in) :: vvH2So4
  real(fp), optional, intent(in) :: vvSO4(:)
  logical,  optional, intent(in) :: write_header
  logical,  optional, intent(in) :: write_data
  integer,  optional, intent(in) :: box_id
  integer,            intent(in) :: n_bins

  integer :: k

  ! Write out the header line
  if (present(write_header)) then
    if (write_header) then
      write(out_id,'(a16,6(",",a16))',advance='no') 'Time (s)','Box','Temp (K)',&
             'Pres (hPa)','# (molec/cm3)','H2O (v/v)','H2SO4 gas (v/v)'
      do k=1,n_bins
        write(out_id,'(",",a4,I3,a6,3x)',advance='no') 'Bin ',k,' (v/v)'
      end do
      ! Advance to new line
      write(out_id,*)
    end if
  end if

  if (present(write_data).and.write_data) then
    ! Better hope all the necessary data is present
    write(out_id,'(I16,",",I16,5(",",E16.5E4))',advance='no') &
        t_now, box_id, T_K, p_hPa, ndens, vvH2O, vvH2SO4
    ! VMR of each aerosol bin
    do k=1,n_bins
      write(out_id,'(",",E16.5E4)',advance='no') vvSO4(k)
    end do
    ! Advance to new line
    write(out_id,*)
  end if

  end subroutine write_state

end module sect_aux_mod
