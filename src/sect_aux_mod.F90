module sect_aux_mod

  use precision_mod
  use netcdf

  implicit none
  private

  PUBLIC  :: write_state
  PUBLIC  :: close_output
  PUBLIC  :: error_stop
  PUBLIC  :: debug_msg

  contains

  subroutine write_state(n_bins,t_now,T_K,p_hPa,ndens,vvH2O,vvH2SO4,vvSO4,out_id,out_file,write_header,write_data,box_id,rc)
  ! Write the current state to an output file (ASCII for now)
  integer,                      intent(inout) :: out_id
  integer,            optional, intent(in   ) :: t_now
  character(len=*),   optional, intent(in   ) :: out_file
  real(fp),           optional, intent(in   ) :: T_K
  real(fp),           optional, intent(in   ) :: p_hPa
  real(fp),           optional, intent(in   ) :: ndens
  real(fp),           optional, intent(in   ) :: vvH2O
  real(fp),           optional, intent(in   ) :: vvH2So4
  real(fp),           optional, intent(in   ) :: vvSO4(:)
  logical,            optional, intent(in   ) :: write_header
  logical,            optional, intent(in   ) :: write_data
  integer,            optional, intent(in   ) :: box_id
  integer,                      intent(in   ) :: n_bins
  integer,            optional, intent(out  ) :: rc

  ! For NetCDF
  logical                      :: file_found
  integer                      :: idx_dim_time
  integer                      :: idx_dim_expt_id
  integer                      :: idx_var_time
  integer                      :: idx_var_expt_id
  integer                      :: idx_var_spc(n_bins+2)
  integer                      :: idx_var, idx_dim
  character(len=NF90_MAX_NAME) :: dim_names_long
  character(len=NF90_MAX_NAME) :: dim_units
  character(len=NF90_MAX_NAME) :: curr_var_name
  character(len=NF90_MAX_NAME) :: curr_att_name
  character(len=255)           :: err_msg
  character(len=255)           :: date_str, time_str, date_time_str
  character(len=80)            :: att_str
  character(len=NF90_MAX_NAME) :: var_long_name
  character(len=NF90_MAX_NAME) :: var_short_name

  character(len=255),parameter :: srt_name='write_state in sect_aux_mod.F90'

  integer :: k, i_bin, i_spc

  ! Assume OK
  If (present(RC)) RC = 0

  ! Create the output file
  if (present(write_header) .and. write_header) then
    if (.not.present(out_file)) then
      rc = 10
      call debug_msg('Cannot generate file ID without an output file')
      Return
    end if

    ! Delete the file if it already exists
    INQUIRE(FILE=Trim(out_file),Exist=file_found)
    if (file_found) then
      call debug_msg('Output file already exists - overwriting')
      open(newunit=out_id,file=trim(out_file),status='old')
      close(out_id,status='delete')
    end if

    ! Generate the NetCDF file
    RC = NF90_CREATE(path=Trim(out_file),cmode=NF90_NOCLOBBER,&
      ncid=out_id)
    If (RC.ne.0) Then
      Call Debug_Msg('Failed to create NetCDF file')
      Return
    End If

    ! Write the dimensions
    curr_var_name = 'time'
    RC = NF90_DEF_DIM(ncid=out_id,name=trim(curr_var_name),len=0,dimid=idx_dim_time)
    If (RC.ne.0) Then
      Write(err_msg,'(3a,I4)') 'Failed to create dimension: ', trim(curr_var_name), &
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If

    curr_var_name = 'expt_id'
    RC = NF90_DEF_DIM(ncid=out_id,name=trim(curr_var_name),len=n_bins+2,dimid=idx_dim_expt_id)
    If (RC.ne.0) Then
      Write(err_msg,'(3a,I4)') 'Failed to create dimension: ', trim(curr_var_name), &
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If

    ! Write the variables
    curr_var_name = 'time'
    RC = NF90_DEF_VAR(ncid=out_id,name=trim(curr_var_name),xtype=NF90_FLOAT,&
        dimids=(/idx_dim_time/),varid=idx_var_time)
    If (RC.ne.0) Then
      Write(err_msg,'(3a,I4)') 'Failed to create variable: ', trim(curr_var_name), &
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If

    curr_var_name = 'expt_id'
    RC = NF90_DEF_VAR(ncid=out_id,name=trim(curr_var_name),xtype=NF90_FLOAT,&
        dimids=(/idx_dim_expt_id/),varid=idx_var_expt_id)
    If (RC.ne.0) Then
      Write(err_msg,'(3a,I4)') 'Failed to create variable: ', trim(curr_var_name), &
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If

    ! Species concentrations
    curr_var_name = 'spc_SO2'
    RC = NF90_DEF_VAR(ncid=out_id,name=trim(curr_var_name),xtype=NF90_FLOAT,&
        dimids=(/idx_dim_expt_id,idx_dim_time/),varid=idx_var_spc(1))
    If (RC.ne.0) Then
      Write(err_msg,'(3a,I4)') 'Failed to create variable: ', trim(curr_var_name), &
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If

    curr_var_name = 'spc_H2SO4'
    RC = NF90_DEF_VAR(ncid=out_id,name=trim(curr_var_name),xtype=NF90_FLOAT,&
        dimids=(/idx_dim_expt_id,idx_dim_time/),varid=idx_var_spc(2))
    If (RC.ne.0) Then
      Write(err_msg,'(3a,I4)') 'Failed to create variable: ', trim(curr_var_name), &
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If

    do i_bin=1,n_bins
      write(curr_var_name,'(a,I0.3)') 'spc_bin', i_bin
      RC = NF90_DEF_VAR(ncid=out_id,name=trim(curr_var_name),xtype=NF90_FLOAT,&
          dimids=(/idx_dim_expt_id,idx_dim_time/),varid=idx_var_spc(2+i_bin))
      If (RC.ne.0) Then
        Write(err_msg,'(3a,I4)') 'Failed to create variable: ', trim(curr_var_name), &
          '. Error code: ', RC
        Call Debug_Msg(trim(err_msg))
        Return
      End If
    End Do

    ! Define attributes
    curr_att_name='history'
    Call Date_And_Time(DATE=date_str,TIME=time_str)
    Write(date_time_str,'(a,x,a)') Trim(date_str), Trim(time_str)
    Write(att_str,'(2a)') Trim(date_time_str),': box model'
    RC = NF90_PUT_ATT(ncid=out_id, varid=NF90_GLOBAL, &
      name=trim(curr_att_name),values=Trim(att_str))
    If (RC.ne.0) Then
      Write(err_msg,'(3a,I4)') 'Failed to write attribute: ', &
        trim(curr_att_name), '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If
    
    curr_att_name='title'
    att_str = 'Sectional aerosol model output data'
    RC = NF90_PUT_ATT(ncid=out_id, varid=NF90_GLOBAL, &
      name=trim(curr_att_name),values=Trim(att_str))
    If (RC.ne.0) Then
      Write(err_msg,'(3a,I4)') 'Failed to write attribute: ', &
        trim(curr_att_name), '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If

    curr_att_name='format'
    att_str = 'netCDF-4'
    RC = NF90_PUT_ATT(ncid=out_id, varid=NF90_GLOBAL, &
      name=trim(curr_att_name),values=Trim(att_str))
    If (RC.ne.0) Then
      Write(err_msg,'(3a,I4)') 'Failed to write attribute: ', &
        trim(curr_att_name), '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If

    ! Dimension attributes
    idx_var = idx_var_time
    curr_var_name = 'time'

    curr_att_name='standard_name'
    att_str = 'time'
    RC = NF90_PUT_ATT(ncid=out_id, varid=idx_var, &
      name=trim(curr_att_name),values=Trim(att_str))
    If (RC.ne.0) Then
      Write(err_msg,'(5a,I4)') 'Failed to write attribute: ', &
        trim(curr_att_name), ' for variable ', Trim(curr_var_name),&
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If
    curr_att_name='long_name'
    att_str = 'time'
    RC = NF90_PUT_ATT(ncid=out_id, varid=idx_var, &
      name=trim(curr_att_name),values=Trim(att_str))
    If (RC.ne.0) Then
      Write(err_msg,'(5a,I4)') 'Failed to write attribute: ', &
        trim(curr_att_name), ' for variable ', Trim(curr_var_name),&
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If
    curr_att_name='units'
    att_str = 'seconds since 2000-01-01 00:00:00'
    RC = NF90_PUT_ATT(ncid=out_id, varid=idx_var, &
      name=trim(curr_att_name),values=Trim(att_str))
    If (RC.ne.0) Then
      Write(err_msg,'(5a,I4)') 'Failed to write attribute: ', &
        trim(curr_att_name), ' for variable ', Trim(curr_var_name),&
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If
    curr_att_name='axis'
    att_str = 'T'
    RC = NF90_PUT_ATT(ncid=out_id, varid=idx_var, &
      name=trim(curr_att_name),values=Trim(att_str))
    If (RC.ne.0) Then
      Write(err_msg,'(5a,I4)') 'Failed to write attribute: ', &
        trim(curr_att_name), ' for variable ', Trim(curr_var_name),&
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If

    idx_var = idx_var_expt_id
    curr_var_name = 'expt_id'

    curr_att_name='standard_name'
    att_str = 'expt_id'
    RC = NF90_PUT_ATT(ncid=out_id, varid=idx_var, &
      name=trim(curr_att_name),values=Trim(att_str))
    If (RC.ne.0) Then
      Write(err_msg,'(5a,I4)') 'Failed to write attribute: ', &
        trim(curr_att_name), ' for variable ', Trim(curr_var_name),&
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If
    curr_att_name='long_name'
    att_str = 'Experiment ID'
    RC = NF90_PUT_ATT(ncid=out_id, varid=idx_var, &
      name=trim(curr_att_name),values=Trim(att_str))
    If (RC.ne.0) Then
      Write(err_msg,'(5a,I4)') 'Failed to write attribute: ', &
        trim(curr_att_name), ' for variable ', Trim(curr_var_name),&
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If
    curr_att_name='units'
    att_str = 'number'
    RC = NF90_PUT_ATT(ncid=out_id, varid=idx_var, &
      name=trim(curr_att_name),values=Trim(att_str))
    If (RC.ne.0) Then
      Write(err_msg,'(5a,I4)') 'Failed to write attribute: ', &
        trim(curr_att_name), ' for variable ', Trim(curr_var_name),&
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If
    curr_att_name='axis'
    att_str = 'X'
    RC = NF90_PUT_ATT(ncid=out_id, varid=idx_var, &
      name=trim(curr_att_name),values=Trim(att_str))
    If (RC.ne.0) Then
      Write(err_msg,'(5a,I4)') 'Failed to write attribute: ', &
        trim(curr_att_name), ' for variable ', Trim(curr_var_name),&
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If

    ! Loop over the chemical species
    Do i_spc = 1, n_bins + 2
      idx_var = idx_var_spc(i_spc)
      If (i_spc == 1) Then
        curr_var_name  = 'spc_SO2'
        var_long_name  = 'SO2 mixing ratio'
        var_short_name = 'spc_SO2'
      Else If (i_spc == 2) Then
        curr_var_name  = 'spc_H2SO4'
        var_long_name  = 'H2SO4(g) mixing ratio'
        var_short_name = 'spc_H2SO4'
      Else
        write(curr_var_name,'(a,I0.3)') 'spc_bin', i_bin
        write(var_long_name,'(a,I3,a,I3)') 'VMR for aerosol bin ',i_bin,' of ',n_bins
        write(var_short_name,'(a,I0.3)') 'spc_bin',i_bin
      End If
      curr_att_name='standard_name'
      att_str = var_short_name
      RC = NF90_PUT_ATT(ncid=out_id, varid=idx_var, &
        name=trim(curr_att_name),values=Trim(att_str))
      If (RC.ne.0) Then
        Write(err_msg,'(5a,I4)') 'Failed to write attribute: ', &
          trim(curr_att_name), ' for variable ', Trim(curr_var_name),&
          '. Error code: ', RC
        Call Debug_Msg(trim(err_msg))
        Return
      End If
      curr_att_name='long_name'
      att_str = var_long_name
      RC = NF90_PUT_ATT(ncid=out_id, varid=idx_var, &
        name=trim(curr_att_name),values=Trim(att_str))
      If (RC.ne.0) Then
        Write(err_msg,'(5a,I4)') 'Failed to write attribute: ', &
          trim(curr_att_name), ' for variable ', Trim(curr_var_name),&
          '. Error code: ', RC
        Call Debug_Msg(trim(err_msg))
        Return
      End If
      curr_att_name='units'
      att_str = 'v/v'
      RC = NF90_PUT_ATT(ncid=out_id, varid=idx_var, &
        name=trim(curr_att_name),values=Trim(att_str))
      If (RC.ne.0) Then
        Write(err_msg,'(5a,I4)') 'Failed to write attribute: ', &
          trim(curr_att_name), ' for variable ', Trim(curr_var_name),&
          '. Error code: ', RC
        Call Debug_Msg(trim(err_msg))
        Return
      End If
      curr_att_name='_FillValue'
      RC = NF90_PUT_ATT(ncid=out_id, varid=idx_var, &
        name=trim(curr_att_name),values=-1.0e+10)
      If (RC.ne.0) Then
        Write(err_msg,'(5a,I4)') 'Failed to write attribute: ', &
          trim(curr_att_name), ' for variable ', Trim(curr_var_name),&
          '. Error code: ', RC
        Call Debug_Msg(trim(err_msg))
        Return
      End If
      curr_att_name='missing_value'
      RC = NF90_PUT_ATT(ncid=out_id, varid=idx_var, &
        name=trim(curr_att_name),values=-1.0e+10)
      If (RC.ne.0) Then
        Write(err_msg,'(5a,I4)') 'Failed to write attribute: ', &
          trim(curr_att_name), ' for variable ', Trim(curr_var_name),&
          '. Error code: ', RC
        Call Debug_Msg(trim(err_msg))
        Return
      End If
    End Do

    ! Leave define mode
    RC = NF90_ENDDEF(ncid=out_id)
    If (RC.ne.0) Then
      Write(err_msg,'(a,I4)') 'Received error while ending define mode: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If

    !write(out_id,'(a16,6(",",a16))',advance='no') 'Time (s)','Box','Temp (K)',&
    !       'Pres (hPa)','# (molec/cm3)','H2O (v/v)','H2SO4 gas (v/v)'
    !do k=1,n_bins
    !  write(out_id,'(",",a4,I3,a6,3x)',advance='no') 'Bin ',k,' (v/v)'
    !end do
    !! Advance to new line
    !write(out_id,*)
  end if

  if (present(write_data).and.write_data) then
    !! Better hope all the necessary data is present
    !write(out_id,'(I16,",",I16,5(",",E16.5E4))',advance='no') &
    !    t_now, box_id, T_K, p_hPa, ndens, vvH2O, vvH2SO4
    !! VMR of each aerosol bin
    !do k=1,n_bins
    !  write(out_id,'(",",E16.5E4)',advance='no') vvSO4(k)
    !end do
    !! Advance to new line
    !write(out_id,*)
  end if

  end subroutine write_state

  subroutine debug_msg(out_msg)
    character(len=*) :: out_msg
    write(*,*) trim(out_msg)
  end subroutine debug_msg

  subroutine error_stop(out_msg,out_loc,rc)
    character(len=*) :: out_msg, out_loc
    integer,optional,intent(in) :: rc
    write(*,'(a,a)') 'SIMULATION FAILED IN ', trim(out_loc)
    write(*,'(a,a)') 'ERROR MSG:  ',trim(out_msg)
    if (present(rc)) then
      write(*,'(a,I4)') 'ERROR CODE: ', RC
    end if
    stop 10
  end subroutine error_stop

  subroutine close_output(out_id, rc)
    integer, intent(in)  :: out_id
    integer, intent(out) :: rc
   
    rc = nf90_close(ncid=out_id)
 
  end subroutine

end module sect_aux_mod
