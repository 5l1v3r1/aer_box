  program run_aerosol
     use precision_mod
     use sect_aer_mod
     !use sect_aer_data_mod

     implicit none

     integer               :: n_boxes
     real(fp),allocatable  :: T_K_Vec(:), p_hPa_Vec(:)
    
     integer               :: t_start, t_sim, t_delta, t_stop, t_coag
     integer               :: n_bins, rc

     write(*,*) 'Initializing simulation.'
     n_boxes = 10
     t_delta = 600 ! Seconds
     t_coag  = 300 ! 5-minute coagulation step
     t_start = 0
     t_stop  = t_start + (6*60*60) ! 6 hours
     n_bins  = 40

     call init_sect_aer(n_bins,rc)
     if (rc.ne.0) then
        write(*,*) 'BAD INITIALIZATION'
        stop 20
     end if

     allocate(T_K_Vec(n_boxes))
     allocate(p_hPa_Vec(n_boxes))

     ! Simple initial conditions
     T_K_Vec(:) = 290.0e+0_fp
     p_hPa_Vec(:) = 90.0e+0_fp

     write(*,*) 'Beginning main time stepping loop.'

     t_sim = t_start
     do while ( t_sim < t_stop )
        ! Simulate t_delta
        !call do_sect_aer(n_boxes,aWP_Arr,aDen_Arr,&
        !                 vvSO4_Arr,Sfc_Ten_Arr,vvH2O_Vec,&
        !                 vvH2SO4_Vec,T_K_Vec,p_hPa_Vec,&
        !                 ndens_Vec,t_delta,t_coag,RC)
        t_sim = t_sim + t_delta
     end do

     ! Deallocate
     call cleanup_sect_aer( .False. )

     if (allocated(T_K_Vec)) deallocate(T_K_Vec)     
     if (allocated(p_hPa_Vec)) deallocate(p_hPa_Vec)

     write(*,*) 'Simulation complete.'

  end program run_aerosol
