  program run_aerosol
     use precision_mod
     use sect_aer_mod
     use sect_aux_mod
     use physconstants, only : rstarg, avo
     !use sect_aer_data_mod

     implicit none

     integer               :: n_boxes
    
     integer               :: t_start, t_sim, t_delta, t_stop, t_coag
     integer               :: n_bins, rc, as, k

     ! Simulation variables
     real(fp), allocatable  :: T_K_Vec(:), p_hPa_Vec(:), nDens_Vec(:)
     real(fp), allocatable  :: vvH2O_Vec(:), vvH2SO4_Vec(:)
     real(fp), allocatable  :: Sfc_Ten_Arr(:,:)
     real(fp), allocatable  :: aDen_Arr(:,:)
     real(fp), allocatable  :: aWP_Arr(:,:)
     real(fp), allocatable  :: vvSO4_Arr(:,:)

     ! Output settings
     integer :: output_fID
     character(len=255) :: output_file

     write(*,*) 'Initializing simulation.'

     ! Default settings - these could be read from a file
     n_boxes = 10
     t_delta = 600 ! Seconds
     t_coag  = 300 ! 5-minute coagulation step
     t_start = 0
     t_stop  = t_start + (6*60*60) ! 6 hours
     n_bins  = 40
     output_file = 'output.dat'

     call init_sect_aer(n_bins,rc)
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

     ! Simple initial conditions
     T_K_Vec(:) = 280.0e+0_fp
     p_hPa_Vec(:) = 90.0e+0_fp

     ! Molecules/cm3 ((m3/cm3) * (molec/mol) * (p/RT), where p/RT = n/V = mol/m3)
     ndens_Vec(:) = 1.0e-6 * AVO * p_hPa_Vec * 100.0e+0_fp / (RStarG * T_K_Vec)

     write(*,'(2(x,E16.5E4))') ndens_vec(1), ndens_vec(n_boxes)

     ! Start with a variety of different initial H2SO4 concs (in ppbv)
     do k=1,n_boxes
        vvh2so4_Vec(k) = dble(k-1) * 1.0e-9_fp
     end do

     ! All start with 50 ppbv H2O
     vvH2O_Vec(:) = 50.0e-9_fp

     ! All start with no aerosol
     vvSO4_Arr(:,:) = 0.0e+0_fp

     ! Recalculate surface tension, weight pcg, and density on first step
     Sfc_Ten_Arr(:,:) = 0.0e+0_fp
     aWP_Arr    (:,:) = 0.0e+0_fp
     aDen_Arr   (:,:) = 0.0e+0_fp

     ! Prepare output file
     OPEN(UNIT=output_fID,FILE=output_file,ACCESS='SEQUENTIAL',&
          FORM='FORMATTED',STATUS='UNKNOWN')

     ! Write the header
     call write_state(n_bins=n_bins,out_id=output_fID,write_header=.True.)

     ! Write the initial state
     
     write(*,*) 'Beginning main time stepping loop.'

     t_sim = t_start
     do while ( t_sim < t_stop )
        ! Simulate t_delta
        !call do_sect_aer(n_boxes,aWP_Arr,aDen_Arr,&
        !                 vvSO4_Arr,Sfc_Ten_Arr,vvH2O_Vec,&
        !                 vvH2SO4_Vec,T_K_Vec,p_hPa_Vec,&
        !                 ndens_Vec,t_delta,t_coag,RC)
        !! Perform output
        !do k=1,n_boxes
        !  call write_state(write_data=.True.,t_now=t_sim,T_K=T_K_Vec(k),&
        !    p_hPa=p_hPa_Vec(k),ndens=ndens_Vec(k),vvH2O=vvH2O_Vec(k),&
        !    vvH2SO4=vvH2SO4_Vec(k),vvSO4=vvSO4_Arr(k,:),&
        !    n_bins=n_bins,out_id=output_fID)
        !end do
        ! Advance time
        t_sim = t_sim + t_delta
     end do

     ! Close the output file
     Close(output_fID)
     If (AS.ne.0) Then
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

     write(*,*) 'Simulation complete.'

  end program run_aerosol
