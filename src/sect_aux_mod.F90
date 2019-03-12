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

  subroutine write_state(n_bins,n_expt,t_now,i_time,T_K,p_hPa,ndens,vvH2O,vvH2SO4,vvSO4,out_id,out_file,write_header,write_data,box_id,rc)
  ! Write the current state to an output file (ASCII for now)
  integer,                              intent(inout) :: out_id
  integer,                              intent(in   ) :: n_bins
  integer,                              intent(in   ) :: n_expt
  integer,                              intent(out  ) :: rc
  integer,            optional,         intent(in   ) :: t_now
  integer,            optional,         intent(in   ) :: i_time
  character(len=*),   optional,         intent(in   ) :: out_file
  real(fp),           optional, target, intent(in   ) :: T_K(n_expt)
  real(fp),           optional, target, intent(in   ) :: p_hPa(n_expt)
  real(fp),           optional, target, intent(in   ) :: ndens(n_expt)
  real(fp),           optional, target, intent(in   ) :: vvH2O(n_expt)
  real(fp),           optional, target, intent(in   ) :: vvH2SO4(n_expt)
  real(fp),           optional, target, intent(in   ) :: vvSO4(n_expt,n_bins)
  logical,            optional,         intent(in   ) :: write_header
  logical,            optional,         intent(in   ) :: write_data
  integer,            optional,         intent(in   ) :: box_id(n_expt)

  integer, parameter           :: n_gas = 3
  integer, parameter           :: n_env = 3

  ! For NetCDF
  logical                      :: file_found
  integer                      :: idx_dim_time
  integer                      :: idx_dim_expt_id
  integer                      :: idx_var_time
  integer                      :: idx_var_expt_id
  integer                      :: idx_var_spc(n_bins+n_gas)
  integer                      :: idx_var_env(n_env)
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
  character(len=NF90_MAX_NAME) :: var_units

  character(len=255),parameter :: srt_name='write_state in sect_aux_mod.F90'

  integer, allocatable :: var_1D_int(:)
  real(f4),allocatable :: var_1D_R4(:)
  real(f4),allocatable :: var_2D_R4(:,:)
  real(fp),pointer     :: ptr_1D_R8(:)
  integer :: j, k, i_bin, i_spc

  ! Assume OK
  RC = 0

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
    RC = NF90_DEF_DIM(ncid=out_id,name=trim(curr_var_name),len=n_expt,dimid=idx_dim_expt_id)
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
    curr_var_name = 'spc_H2O'
    RC = NF90_DEF_VAR(ncid=out_id,name=trim(curr_var_name),xtype=NF90_FLOAT,&
        dimids=(/idx_dim_expt_id,idx_dim_time/),varid=idx_var_spc(1))
    If (RC.ne.0) Then
      Write(err_msg,'(3a,I4)') 'Failed to create variable: ', trim(curr_var_name), &
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If

    curr_var_name = 'spc_SO2'
    RC = NF90_DEF_VAR(ncid=out_id,name=trim(curr_var_name),xtype=NF90_FLOAT,&
        dimids=(/idx_dim_expt_id,idx_dim_time/),varid=idx_var_spc(2))
    If (RC.ne.0) Then
      Write(err_msg,'(3a,I4)') 'Failed to create variable: ', trim(curr_var_name), &
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If

    curr_var_name = 'spc_H2SO4'
    RC = NF90_DEF_VAR(ncid=out_id,name=trim(curr_var_name),xtype=NF90_FLOAT,&
        dimids=(/idx_dim_expt_id,idx_dim_time/),varid=idx_var_spc(3))
    If (RC.ne.0) Then
      Write(err_msg,'(3a,I4)') 'Failed to create variable: ', trim(curr_var_name), &
        '. Error code: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If

    ! Aerosol bins
    do i_bin=1,n_bins
      write(curr_var_name,'(a,I0.3)') 'spc_bin', i_bin
      RC = NF90_DEF_VAR(ncid=out_id,name=trim(curr_var_name),xtype=NF90_FLOAT,&
          dimids=(/idx_dim_expt_id,idx_dim_time/),varid=idx_var_spc(n_gas+i_bin))
      If (RC.ne.0) Then
        Write(err_msg,'(3a,I4)') 'Failed to create variable: ', trim(curr_var_name), &
          '. Error code: ', RC
        Call Debug_Msg(trim(err_msg))
        Return
      End If
    End Do

    ! Environmental variables
    do k=1,3
       if (k==1) curr_var_name = 'T'
       if (k==2) curr_var_name = 'p'
       if (k==3) curr_var_name = 'dens_air'
       RC = NF90_DEF_VAR(ncid=out_id,name=trim(curr_var_name),xtype=NF90_FLOAT,&
           dimids=(/idx_dim_expt_id,idx_dim_time/),varid=idx_var_env(k))
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
    Do i_spc = 1, n_bins + n_gas
      idx_var = idx_var_spc(i_spc)
      If (i_spc == 1) Then
        curr_var_name  = 'spc_H2O'
        var_long_name  = 'H2O mixing ratio'
        var_short_name = 'spc_H2O'
      Else If (i_spc == 2) Then
        curr_var_name  = 'spc_SO2'
        var_long_name  = 'SO2 mixing ratio'
        var_short_name = 'spc_SO2'
      Else If (i_spc == 3) Then
        curr_var_name  = 'spc_H2SO4'
        var_long_name  = 'H2SO4(g) mixing ratio'
        var_short_name = 'spc_H2SO4'
      Else
        i_bin = i_spc - 3
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

    ! Loop over the environmental variables
    Do k = 1, 3
      idx_var = idx_var_env(k)
      If (k == 1) Then
        curr_var_name  = 'T'
        var_long_name  = 'Temperature'
        var_short_name = 'T'
        var_units      = 'K'
      Else If (k == 2) Then
        curr_var_name  = 'p'
        var_long_name  = 'Pressure'
        var_short_name = 'p'
        var_units      = 'hPa'
      Else If (k == 3) Then
        curr_var_name  = 'dens_air'
        var_long_name  = 'Number density of air'
        var_short_name = 'dens_air'
        var_units      = '#/cm3'
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
      att_str = var_units
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

    ! Finally, write the dimensions (where possible)
    Allocate(var_1D_int(n_expt),stat=rc)
    RC = NF90_PUT_VAR(ncid=out_id,varid=idx_var_expt_id,&
       values=var_1D_Int)
    If (RC.ne.0) Then
      Write(err_msg,'(a,I4)') 'Received error while writing expt ID var: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If
    Deallocate(var_1D_int)

  end if

  if (present(write_data).and.write_data) then
    ! Write new time entry
    RC = NF90_INQ_VARID(ncid=out_id,name='time',varid=idx_var)
    If (RC.ne.0) Then
      Write(err_msg,'(a,I4)') 'Received error while requesting time index: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If
    RC = NF90_PUT_VAR(ncid=out_id,varid=idx_var,&
       values=(/ real(t_now) /),start=(/i_time/))
    If (RC.ne.0) Then
      Write(err_msg,'(a,I4)') 'Received error while writing time var: ', RC
      Call Debug_Msg(trim(err_msg))
      Return
    End If
    ! Assume that the calling routine knows the output time
    Do k=1,(n_bins + n_gas + n_env)
      If (k <= (n_bins + n_gas)) Then
        if (k == 1) then
          ptr_1D_r8 => vvH2O
          curr_var_name = 'spc_H2O'
        else if (k == 2) then
          ! Not actually present..
          cycle
          !ptr_1D_r8 => vvSO2
          !curr_var_name = 'spc_SO2'
        else if (k == 3) then
          ptr_1D_r8 => vvH2SO4
          curr_var_name = 'spc_H2SO4'
        else
          ptr_1D_r8 => vvSO4(:,k-n_gas)
          write(curr_var_name,'(a,I0.3)') 'spc_bin', k-n_gas
        end if
      Else
        j = k - (n_bins + n_gas)
        if (j == 1) then
           ptr_1d_r8 => T_K
           curr_var_name = 'T'
        else if (j == 2) then
           ptr_1d_r8 => p_hPa
           curr_var_name = 'p'
        else if (j == 3) then
           ptr_1d_r8 => ndens
           curr_var_name = 'dens_air'
        end if
      End If
      ! Get the variable index
      RC = NF90_INQ_VARID(ncid=out_id,name=Trim(curr_var_name),varid=idx_var)
      If (RC.ne.0) Then
        Write(err_msg,'(2a)') 'Could not find ID for variable ',Trim(curr_var_name)
        Call Debug_Msg(trim(err_msg))
        Return
      End If
      ! Copy to R4
      allocate(var_2d_r4(n_expt,1),stat=rc)
      if (rc.ne.0) then
        call debug_msg('Could not allocated Var_2D_R4 in Write_State')
        return
      end if
      !var_2d_r4(:,1) = Real(ptr_1d_r8(:),Kind=f4)
      Do J=1,n_expt
        var_2d_r4(J,1) = Real(ptr_1d_r8(j),Kind=f4)
      End Do
      Nullify(ptr_1d_r8)
      ! Write data
      RC = NF90_PUT_VAR(ncid=out_id,varid=idx_var,&
        values=var_2d_r4,start=(/1,i_time/))
      deallocate(var_2d_r4)
    End Do
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
