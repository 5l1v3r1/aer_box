  program run_aerosol
     use precision_mod
     use sect_aer_mod
     use sect_aer_data_mod, only : aer_dry_rad, aer_Vrat
     use sect_aux_mod
     use physconstants, only : rstarg, avo
     !use sect_aer_data_mod

     implicit none

     integer               :: n_boxes
    
     integer               :: t_start, t_sim, t_stop
     integer               :: dt_main, dt_coag
     integer               :: n_bins, rc, as, k
     integer               :: dt_output, t_next_output
     integer               :: output_idx
     logical               :: LDebug, LImplicit_Coag
     logical               :: LNuc, LGrow, LCoag

     ! Simulation variables
     real(fp), allocatable  :: T_K_Vec(:), p_hPa_Vec(:), nDens_Vec(:)
     real(fp), allocatable  :: vvH2O_Vec(:), vvH2SO4_Vec(:)
     real(fp), allocatable  :: Sfc_Ten_Arr(:,:)
     real(fp), allocatable  :: aDen_Arr(:,:)
     real(fp), allocatable  :: aWP_Arr(:,:)
     real(fp), allocatable  :: vvSO4_Arr(:,:)
     real(fp), allocatable  :: rWet_Arr(:,:)

     ! Output settings
     integer :: output_fID, idx_output
     character(len=255) :: output_file

     ! Initial conditions
     real(fp)               :: T_K_min,     T_K_max
     real(fp)               :: p_hPa_min,   p_hPa_max
     real(fp)               :: vvH2SO4_min, vvH2SO4_max
     real(fp)               :: vvH2O_Init,  vvSO2_Init
     
     ! H2SO4 formation rates
     real(fp)               :: H2SO4_per_day, H2SO4_per_step

     ! Initial population
     real(fp), allocatable  :: vvSO4_Init(:)
     real(fp)               :: vvAero_Init, rAero_Init, gsdAero_Init
     real(fp)               :: r_min, r_max, ln_gsd, aero_factor
     real(fp)               :: r_median, erf_next, erf_last, rv_median

     write(*,*) 'Initializing simulation.'

     call read_input('input.box',dt_main,dt_output,&
        dt_coag,t_start,t_stop,n_bins,n_boxes,&
        output_file,T_K_Min,T_K_Max,p_hPa_min,p_hPa_max,&
        vvH2SO4_min,vvH2SO4_max,vvH2O_Init,vvSO2_Init,&
        LImplicit_Coag,LDebug,LNuc,LGrow,LCoag,&
        H2SO4_per_day,vvAero_Init,rAero_Init,gsdAero_Init,&
        rc)
     if (rc.ne.0) Then
       Call error_stop('Failed to read input file','main',rc)
     End If

     if (dt_output < dt_main) dt_output = dt_main

     ! Calculate H2SO4 formation rate
     H2SO4_per_step = H2SO4_per_day * real(dt_main,fp) / (24.0e+0_fp*3600.0e+0_fp)

     call init_sect_aer(n_boxes,n_bins,rc)
     if (rc.ne.0) then
        write(*,*) 'BAD INITIALIZATION'
        stop 20
     end if

     allocate(T_K_Vec(n_boxes),stat=as)
     if (as.ne.0) then
        write(*,*) 'Failed to allocate T_K_Vec'
        stop
     end if
     T_K_Vec(:) = 0.0e+0_fp

     allocate(p_hPa_Vec(n_boxes),stat=as)
     if (as.ne.0) then
        write(*,*) 'Failed to allocate p_hPa_Vec'
        stop
     end if
     p_hPa_Vec(:) = 0.0e+0_fp

     allocate(ndens_Vec(n_boxes),stat=as)
     if (as.ne.0) then
        write(*,*) 'Failed to allocate ndens_Vec'
        stop
     end if
     ndens_Vec(:) = 0.0e+0_fp

     allocate(vvH2O_Vec(n_boxes),stat=as)
     if (as.ne.0) then
        write(*,*) 'Failed to allocate vvH2O_Vec'
        stop
     end if
     vvH2O_Vec(:) = 0.0e+0_fp

     allocate(vvH2SO4_Vec(n_boxes),stat=as)
     if (as.ne.0) then
        write(*,*) 'Failed to allocate vvH2SO4_Vec'
        stop
     end if
     vvH2SO4_Vec(:) = 0.0e+0_fp

     allocate(vvSO4_Arr(n_boxes,n_bins),stat=as)
     if (as.ne.0) then
        write(*,*) 'Failed to allocate vvSO4_Arr'
        stop
     end if
     vvSO4_Arr(:,:) = 0.0e+0_fp

     allocate(aDen_Arr(n_boxes,n_bins),stat=as)
     if (as.ne.0) then
        write(*,*) 'Failed to allocate aDen_Arr'
        stop
     end if
     aDen_Arr(:,:) = 0.0e+0_fp

     allocate(aWP_Arr(n_boxes,n_bins),stat=as)
     if (as.ne.0) then
        write(*,*) 'Failed to allocate aWP_Arr'
        stop
     end if
     aWP_Arr(:,:) = 0.0e+0_fp

     allocate(Sfc_Ten_Arr(n_boxes,n_bins),stat=as)
     if (as.ne.0) then
        write(*,*) 'Failed to allocate Sfc_Ten_Arr'
        stop
     end if
     Sfc_Ten_Arr(:,:) = 0.0e+0_fp

     allocate(rWet_Arr(n_boxes,n_bins),stat=as)
     if (as.ne.0) then
        write(*,*) 'Failed to allocate rWet_Arr'
        stop
     end if
     rWet_Arr(:,:) = 0.0e+0_fp

     ! Set initial conditions
     if (n_boxes.eq.1) then
       T_K_vec(1)     = (T_K_min     + T_K_max    )/2.0
       p_hPa_vec(1)   = (p_hPa_min   + p_hPa_max  )/2.0
       vvH2SO4_vec(1) = (vvH2SO4_min + vvH2SO4_max)/2.0
     else
       do k=1,n_boxes
         p_hPa_vec(K)   = ((dble(k-1)/dble(n_boxes-1)) * &
            (p_hPa_max   - p_hPa_min  )) + p_hPa_min
         T_K_Vec(k)     = ((dble(k-1)/dble(n_boxes-1)) * &
            (T_K_max     - T_K_min    )) + T_K_min
         vvH2SO4_Vec(k) = ((dble(k-1)/dble(n_boxes-1)) * &
            (vvH2SO4_max - vvH2SO4_min)) + vvH2SO4_min
       end do
     end if

     ! Molecules/cm3 ((m3/cm3) * (molec/mol) * (p/RT), where p/RT = n/V = mol/m3)
     ndens_Vec(:) = 1.0e-6 * AVO * p_hPa_Vec * 100.0e+0_fp / (RStarG * T_K_Vec)

     ! All start with 50 ppbv H2O
     vvH2O_Vec(:) = vvH2O_Init

     ! All start with no aerosol
     vvSO4_Arr(:,:) = 0.0e-9_fp

     ! Initialize wet radii = dry radii
     do k = 1, n_boxes
        rWet_Arr(k,:) = aer_dry_rad(:)
     end do

     If (vvAero_Init > 0.0e+0_fp) Then
        ! Calculate the volume in each bin then just scale
        allocate(vvSO4_Init(n_bins),stat=as)
        if (as.ne.0) then
           write(*,*) 'Failed to allocate vvSO4_Init'
           stop
        end if
        vvSO4_Init(:) = 0.0e+0_fp
        ! Constant factor
        ln_gsd = log(gsdAero_Init)
        aero_factor = 1.0 / (sqrt(2.0) * ln_gsd)
        ! Assume we were given the number-weighted mode radius
        r_median = rAero_Init * exp(0.5 * ln_gsd * ln_gsd)
        ! Need volume-weighted median
        rv_median = r_median * exp(3.0 * ln_gsd * ln_gsd)
        ! Get the lower bound of bin 1 based on a hypothetical bin 0
        !r_max = (aer_dry_rad(1)/(aer_Vrat ** (1.0/3.0))) * (2.0 * aer_Vrat / (1.0 + aer_Vrat))**(1.0/3.0)
        r_max = aer_dry_rad(1) * (2.0 / (1.0 + aer_Vrat))**(1.0/3.0)
        erf_next = erf(log(r_max/rv_median)*aero_factor)
        do k=1,n_bins
           erf_last = erf_next
           r_max = aer_dry_rad(k) * (2.0 * aer_Vrat / (1.0 + aer_Vrat))**(1.0/3.0)
           erf_next = erf(log(r_max/rv_median)*aero_factor)
           vvSO4_Init(k) = 0.5 * (erf_next - erf_last)
        end do
        ! Scale up the total
        If (sum(vvSO4_Init) > 0.0) Then
           vvSO4_Init(:) = vvAero_Init * vvSO4_Init(:) / sum(vvSO4_Init)
        Else
           write(*,*) 'Non-positive initial aerosol distribution?'
           stop 20
        End If
        do k=1,n_boxes
           vvSO4_Arr(k,:) = vvSO4_Init(:)
        end do
     End If

     ! Recalculate surface tension, weight pcg, and density on first step
     Sfc_Ten_Arr(:,:) = 0.0e+0_fp
     aWP_Arr    (:,:) = 0.0e+0_fp
     aDen_Arr   (:,:) = 0.0e+0_fp

     ! Prepare output file
     call write_state(n_bins=n_bins,n_expt=n_boxes,out_id=output_fID,&
                      out_file=trim(output_file),&
                      write_header=.True.,rc=rc)
     if (rc.ne.0) then
       call error_stop('Could not prepare output file','main_sect')
     end if

     ! Write the initial state
     idx_output = 1
     call write_state(write_data=.True.,t_now=t_start,T_K=T_K_Vec,&
       p_hPa=p_hPa_Vec,ndens=ndens_Vec,vvH2O=vvH2O_Vec,rWet=rWet_Arr,&
       vvH2SO4=vvH2SO4_Vec,vvSO4=vvSO4_Arr,i_time=idx_output,&
       n_expt=n_boxes,n_bins=n_bins,out_id=output_fID,rc=rc)
     if (rc.ne.0) then
       call error_stop('Could not perform initial write','main_sect')
     end if

     ! State radii
     If (LDebug) Then
       Write(*,*) 'Aerosol radius data:'
       Do k=1,size(aer_dry_rad)
         Write(*,'(" => Bin ",I3,": ",F10.6," um")') k, aer_dry_rad(k)
       End Do
     end if
     
     write(*,*) 'Beginning main time stepping loop.'

     t_sim = t_start
     t_next_output = t_start + dt_output
     do while ( t_sim < t_stop )

        ! Add any H2SO4 being produced
        If (H2SO4_per_day > 0.0e+0_fp) &
        vvH2SO4_Vec(:) = vvH2SO4_Vec(:) + H2SO4_per_step

        ! Simulate dt_main
        call do_sect_aer(n_boxes,aWP_Arr,aDen_Arr,&
                         vvSO4_Arr,Sfc_Ten_Arr,vvH2O_Vec,&
                         vvH2SO4_Vec,rWet_Arr,T_K_Vec,p_hPa_Vec,&
                         ndens_Vec,dt_main,dt_coag,&
                         lnuc, lgrow, lcoag,&
                         limplicit_coag,RC)
        ! Advance time
        t_sim = t_sim + dt_main

        ! Perform output
        if ((t_sim >= t_stop).or.(t_sim >= t_next_output)) then
           if (ldebug) then
              write(*,'(a,I10)') 'Writing data for t = ', t_sim
           end if
           idx_output = idx_output + 1
           call write_state(write_data=.True.,t_now=t_sim,T_K=T_K_Vec,&
             p_hPa=p_hPa_Vec,ndens=ndens_Vec,vvH2O=vvH2O_Vec,rWet=rWet_Arr,&
             vvH2SO4=vvH2SO4_Vec,vvSO4=vvSO4_Arr,i_time=idx_output,&
             n_expt=n_boxes,n_bins=n_bins,out_id=output_fID,rc=rc)
           if (rc.ne.0) then
             call error_stop('Could not perform mid-simulation write','main_sect')
           end if
           if (t_sim >= t_next_output) t_next_output = t_next_output + dt_output
        end if
     end do

     ! Close the output file
     !Close(output_fID)
     Call close_output(output_fID,rc)
     If (rc.ne.0) Then
       write(*,*) 'Could not close output file'
       stop
     End If

     ! Deallocate and clean up
     call cleanup_sect_aer( .False. )

     if (allocated(T_K_Vec    )) deallocate(T_K_Vec)     
     if (allocated(p_hPa_Vec  )) deallocate(p_hPa_Vec)
     if (allocated(ndens_vec  )) deallocate(ndens_vec)
     if (allocated(vvh2o_vec  )) deallocate(vvh2o_vec)
     if (allocated(vvh2so4_vec)) deallocate(vvh2so4_vec)
     if (allocated(sfc_ten_arr)) deallocate(sfc_ten_arr)
     if (allocated(awp_arr    )) deallocate(awp_arr)
     if (allocated(aden_arr   )) deallocate(aden_arr)
     if (allocated(vvso4_arr  )) deallocate(vvso4_arr)
     if (allocated(rWet_Arr   )) deallocate(rWet_Arr )

     write(*,*) 'Simulation complete.'

  end program run_aerosol
