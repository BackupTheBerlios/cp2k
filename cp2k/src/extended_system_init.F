!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      cjm, FEB 20 2001:  added subroutine initialize_extended_parameters
!>      cjm, MAY 03 2001:  reorganized and added separtate routines for
!>                         nhc_part, nhc_baro, nhc_ao, npt
!> \author CJM
! *****************************************************************************
MODULE extended_system_init

  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE extended_system_mapping,         ONLY: nhc_to_barostat_mapping,&
                                             nhc_to_particle_mapping,&
                                             nhc_to_shell_mapping
  USE extended_system_types,           ONLY: debug_isotropic_limit,&
                                             lnhc_parameters_type,&
                                             map_info_type,&
                                             npt_info_type
  USE f77_blas
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: do_thermo_only_master,&
                                             npe_f_ensemble,&
                                             npe_i_ensemble,&
                                             nph_uniaxial_damped_ensemble,&
                                             nph_uniaxial_ensemble,&
                                             npt_f_ensemble,&
                                             npt_i_ensemble
  USE input_cp2k_binary_restarts,      ONLY: read_binary_thermostats_nose
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_remove_values,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: global_constraint_type,&
                                             molecule_type
  USE parallel_rng_types,              ONLY: next_random_number
  USE simpar_types,                    ONLY: simpar_type
  USE termination,                     ONLY: stop_program
  USE thermostat_types,                ONLY: thermostat_info_type
  USE thermostat_utils,                ONLY: get_nhc_energies
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'extended_system_init'

  PUBLIC :: initialize_nhc_part, initialize_nhc_baro, initialize_npt, &
            initialize_nhc_shell

CONTAINS

! *****************************************************************************
!> \author CJM
! *****************************************************************************
  SUBROUTINE initialize_npt ( simpar, globenv, npt_info, cell, work_section, error)

    TYPE(simpar_type), POINTER               :: simpar
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(npt_info_type), DIMENSION(:, :), &
      POINTER                                :: npt_info
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: work_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_npt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ind, j, stat
    LOGICAL                                  :: explicit, failure, restart
    REAL(KIND=dp)                            :: temp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: buffer
    TYPE(section_vals_type), POINTER         :: work_section2

    CALL timeset(routineN,handle)

    NULLIFY (work_section2)

    explicit = .FALSE.
    failure=.FALSE.
    restart = .FALSE.

    CPPrecondition(.NOT.ASSOCIATED(npt_info),cp_fatal_level,routineP,error,failure)

    ! first allocating the npt_info_type if requested
    SELECT CASE ( simpar % ensemble )
    CASE ( npt_i_ensemble, npe_i_ensemble )
       ALLOCATE ( npt_info ( 1, 1 ), STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       npt_info  ( :, : ) % eps = LOG ( cell % deth ) / 3.0_dp
       temp = simpar % temp_baro_ext

    CASE ( npt_f_ensemble, npe_f_ensemble )
       ALLOCATE ( npt_info ( 3, 3 ), STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       temp = simpar % temp_baro_ext

    CASE ( nph_uniaxial_ensemble )
       ALLOCATE ( npt_info ( 1, 1 ), STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       temp = simpar% temp_baro_ext

    CASE (  nph_uniaxial_damped_ensemble )
       ALLOCATE ( npt_info ( 1, 1 ), STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       temp = simpar % temp_baro_ext

    CASE  DEFAULT
       ! Do nothing..
       NULLIFY(npt_info)
    END SELECT

    IF (ASSOCIATED(npt_info)) THEN
       IF (ASSOCIATED(work_section)) THEN
          work_section2 => section_vals_get_subs_vals(work_section,"VELOCITY",error=error)
          CALL section_vals_get(work_section2, explicit=explicit, error=error)
          restart=explicit
          work_section2 => section_vals_get_subs_vals(work_section,"MASS",error=error)
          CALL section_vals_get(work_section2, explicit=explicit, error=error)
          CALL cp_assert(restart.EQV.explicit,cp_failure_level,cp_assertion_failed,routineP,&
               "You need to define both VELOCITY and MASS section (or none) in the BAROSTAT section",&
               error=error, failure=failure)
          restart=explicit.AND.restart
       END IF
       
       IF ( restart ) THEN
          CALL section_vals_val_get(work_section,"VELOCITY%_DEFAULT_KEYWORD_",r_vals=buffer,error=error)
          ind = 0
          DO i = 1, SIZE(npt_info,1)
             DO j = 1, SIZE(npt_info,2)
                ind = ind + 1
                npt_info ( i, j ) % v = buffer(ind)
             END DO
          END DO
          CALL section_vals_val_get(work_section,"MASS%_DEFAULT_KEYWORD_",r_vals=buffer,error=error)
          ind = 0
          DO i = 1, SIZE(npt_info,1)
             DO j = 1, SIZE(npt_info,2)
                ind = ind + 1
                npt_info ( i, j ) % mass = buffer(ind)
             END DO
          END DO
       ELSE
          CALL init_barostat_variables ( npt_info, simpar % tau_cell, temp,  &
               simpar % nfree, simpar % ensemble, simpar % cmass,&
               globenv ,error=error)
       END IF

    END IF

    CALL timestop(handle)

  END SUBROUTINE initialize_npt

! *****************************************************************************
!> \brief fire up the thermostats, if NPT
!> \author CJM
! *****************************************************************************
  SUBROUTINE initialize_nhc_baro ( simpar, para_env, globenv, nhc, nose_section, save_mem, error )

    TYPE(simpar_type), POINTER               :: simpar
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(section_vals_type), POINTER         :: nose_section
    LOGICAL, INTENT(IN)                      :: save_mem
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_nhc_baro', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure, restart
    REAL(KIND=dp)                            :: temp

    CALL timeset(routineN,handle)

    failure = .FALSE.
    restart = .FALSE.

    CALL nhc_to_barostat_mapping ( simpar, nhc, error)

    ! Set up the Yoshida weights
    IF ( nhc % nyosh > 0 ) THEN
       ALLOCATE ( nhc % dt_yosh ( 1 : nhc%nyosh ), STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       CALL set_yoshida_coef ( nhc, simpar % dt )
    END IF

    CALL restart_nose(nhc,nose_section,save_mem,restart,"",para_env,error)

    IF (.NOT.restart) THEN
       ! Initializing thermostat forces and velocities for the Nose-Hoover
       ! Chain variables
       SELECT CASE (simpar % ensemble )
       CASE DEFAULT
          temp =  simpar%temp_baro_ext
       END SELECT
       IF ( nhc % nhc_len /= 0 ) THEN
          CALL init_nhc_variables(nhc,temp,para_env,globenv,error=error)
       END IF
    END IF

    CALL init_nhc_forces ( nhc, error=error)

    CALL timestop(handle)

  END SUBROUTINE initialize_nhc_baro

! *****************************************************************************
!> \author CJM
! *****************************************************************************
  SUBROUTINE initialize_nhc_part ( thermostat_info, simpar, local_molecules,&
       molecule, molecule_kind_set, para_env, globenv, nhc, nose_section,&
       gci, save_mem, binary_restart_file_name, error)

    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(section_vals_type), POINTER         :: nose_section
    TYPE(global_constraint_type), POINTER    :: gci
    LOGICAL, INTENT(IN)                      :: save_mem
    CHARACTER(LEN=*), INTENT(IN)             :: binary_restart_file_name
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_nhc_part', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure, restart

    CALL timeset(routineN,handle)

    failure = .FALSE.
    restart = .FALSE.
    ! fire up the thermostats, if not NVE

    CALL nhc_to_particle_mapping ( thermostat_info, simpar, local_molecules,&
         molecule, molecule_kind_set, nhc, para_env, gci, error )

    ! Set up the Yoshida weights
    IF ( nhc % nyosh > 0 ) THEN
       ALLOCATE ( nhc % dt_yosh ( 1 : nhc%nyosh ), STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       CALL set_yoshida_coef ( nhc, simpar % dt)
    END IF

    CALL restart_nose(nhc,nose_section,save_mem,restart,binary_restart_file_name,&
                      para_env,error)

    IF (.NOT.restart) THEN
       ! Initializing thermostat forces and velocities for the Nose-Hoover
       ! Chain variables
       IF ( nhc % nhc_len /= 0 ) THEN
          CALL init_nhc_variables(nhc,simpar%temp_ext,para_env,globenv,error=error)
       END IF
    END IF

    CALL init_nhc_forces ( nhc, error=error)

    CALL timestop(handle)

  END SUBROUTINE initialize_nhc_part

! *****************************************************************************
!> \author MI
! *****************************************************************************
  SUBROUTINE initialize_nhc_shell( thermostat_info, simpar, local_molecules,&
       molecule, molecule_kind_set, para_env, globenv, nhc, nose_section,&
       gci, save_mem, binary_restart_file_name, error)
 
    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(section_vals_type), POINTER         :: nose_section
    TYPE(global_constraint_type), POINTER    :: gci
    LOGICAL, INTENT(IN)                      :: save_mem
    CHARACTER(LEN=*), INTENT(IN)             :: binary_restart_file_name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_nhc_shell', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure, restart

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CALL nhc_to_shell_mapping(thermostat_info, simpar, local_molecules,&
       molecule, molecule_kind_set, nhc, para_env, gci, error)

    restart = .FALSE.
    ! Set up the Yoshida weights
    IF ( nhc % nyosh > 0 ) THEN
       ALLOCATE ( nhc % dt_yosh ( 1 : nhc%nyosh ), STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       CALL set_yoshida_coef ( nhc, simpar % dt )
    END IF

    CALL restart_nose(nhc,nose_section,save_mem,restart,binary_restart_file_name,&
                      para_env,error)
    
    IF (.NOT.restart) THEN
       ! Initialize thermostat forces and velocities
       ! Chain variables
       IF ( nhc % nhc_len /= 0 ) THEN
          CALL init_nhc_variables (nhc, simpar%temp_sh_ext, para_env, globenv, error=error)
       END IF
    END IF

    CALL init_nhc_forces ( nhc, error=error)

    CALL timestop(handle)

  END SUBROUTINE  initialize_nhc_shell
  
! *****************************************************************************
!> \brief This lists the coefficients for the Yoshida method (higher
!>      order integrator used in NVT)
!> \par History
!>      none
!> \date 14-NOV-2000
! *****************************************************************************
  SUBROUTINE set_yoshida_coef ( nhc, dt )

    TYPE(lnhc_parameters_type), POINTER      :: nhc
    REAL(KIND=dp), INTENT(IN)                :: dt

    CHARACTER(len=*), PARAMETER :: routineN = 'set_yoshida_coef', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp), DIMENSION(nhc%nyosh)      :: yosh_wt

    SELECT CASE (nhc % nyosh)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,'Value not available.')
    CASE (1)
       yosh_wt(1) = 1.0_dp
    CASE (3)
       yosh_wt(1) = 1.0_dp/(2.0_dp-(2.0_dp)**(1.0_dp/3.0_dp))
       yosh_wt(2) = 1.0_dp - 2.0_dp*yosh_wt(1)
       yosh_wt(3) = yosh_wt(1)
    CASE (5)
       yosh_wt(1) = 1.0_dp/(4.0_dp-(4.0_dp)**(1.0_dp/3.0_dp))
       yosh_wt(2) = yosh_wt(1)
       yosh_wt(4) = yosh_wt(1)
       yosh_wt(5) = yosh_wt(1)
       yosh_wt(3) = 1.0_dp - 4.0_dp*yosh_wt(1)
    CASE (7)
       yosh_wt(1) = .78451361047756_dp
       yosh_wt(2) = .235573213359357_dp
       yosh_wt(3) = -1.17767998417887_dp
       yosh_wt(4) = 1.0_dp - 2.0_dp*(yosh_wt(1)+yosh_wt(2)+yosh_wt(3))
       yosh_wt(5) = yosh_wt(3)
       yosh_wt(6) = yosh_wt(2)
       yosh_wt(7) = yosh_wt(1)
    CASE (9)
       yosh_wt(1) = 0.192_dp
       yosh_wt(2) = 0.554910818409783619692725006662999_dp
       yosh_wt(3) = 0.124659619941888644216504240951585_dp
       yosh_wt(4) = -0.843182063596933505315033808282941_dp
       yosh_wt(5) = 1.0_dp - 2.0_dp*(yosh_wt(1)+yosh_wt(2)+&
            yosh_wt(3)+yosh_wt(4))
       yosh_wt(6) = yosh_wt(4)
       yosh_wt(7) = yosh_wt(3)
       yosh_wt(8) = yosh_wt(2)
       yosh_wt(9) = yosh_wt(1)
    CASE (15)
       yosh_wt(1) = 0.102799849391985_dp
       yosh_wt(2) = -0.196061023297549e1_dp
       yosh_wt(3) = 0.193813913762276e1_dp
       yosh_wt(4) = -0.158240635368243_dp
       yosh_wt(5) = -0.144485223686048e1_dp
       yosh_wt(6) = 0.253693336566229_dp
       yosh_wt(7) = 0.914844246229740_dp
       yosh_wt(8) = 1.0_dp - 2.0_dp*(yosh_wt(1)+yosh_wt(2)+&
            yosh_wt(3)+yosh_wt(4)+yosh_wt(5)+yosh_wt(6)+yosh_wt(7))
       yosh_wt(9) = yosh_wt(7)
       yosh_wt(10) = yosh_wt(6)
       yosh_wt(11) = yosh_wt(5)
       yosh_wt(12) = yosh_wt(4)
       yosh_wt(13) = yosh_wt(3)
       yosh_wt(14) = yosh_wt(2)
       yosh_wt(15) = yosh_wt(1)
    END SELECT
    nhc % dt_yosh = dt * yosh_wt / REAL ( nhc % nc,KIND=dp)

  END SUBROUTINE set_yoshida_coef

! *****************************************************************************
!> \brief read coordinate, velocities, forces and masses of the
!>      thermostat from restart file
!> \par History
!>     24-07-07 created 
!> \author MI
! *****************************************************************************
  SUBROUTINE restart_nose(nhc,nose_section,save_mem,restart,&
                          binary_restart_file_name,para_env,error)

    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(section_vals_type), POINTER         :: nose_section
    LOGICAL, INTENT(IN)                      :: save_mem
    LOGICAL, INTENT(OUT)                     :: restart
    CHARACTER(LEN=*), INTENT(IN)             :: binary_restart_file_name
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'restart_nose', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ind, j
    LOGICAL                                  :: explicit
    REAL(KIND=dp), DIMENSION(:), POINTER     :: buffer
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(section_vals_type), POINTER         :: work_section

    CALL timeset(routineN,handle)

    NULLIFY (buffer)
    NULLIFY (work_section)

    IF (LEN_TRIM(binary_restart_file_name) > 0) THEN

       ! Read binary restart file, if available

       CALL read_binary_thermostats_nose(nhc,binary_restart_file_name,para_env,error)
       restart = .TRUE.

    ELSE

       ! Read the default restart file in ASCII format

       explicit = .FALSE.
       restart = .FALSE.

       IF (ASSOCIATED(nose_section)) THEN
          work_section => section_vals_get_subs_vals(nose_section,"VELOCITY",error=error)
          CALL section_vals_get(work_section, explicit=explicit, error=error)
          restart=explicit
          work_section => section_vals_get_subs_vals(nose_section,"COORD",error=error)
          CALL section_vals_get(work_section, explicit=explicit, error=error)
          CALL cp_assert(restart.or..not.explicit,cp_failure_level,cp_assertion_failed,routineP,&
               "You need to define both VELOCITY and COORD and MASS and FORCE section (or none) in the NOSE section",&
               error=error)
          restart=explicit.and.restart
          work_section => section_vals_get_subs_vals(nose_section,"MASS",error=error)
          CALL section_vals_get(work_section, explicit=explicit, error=error)
          CALL cp_assert(restart.or..not.explicit,cp_failure_level,cp_assertion_failed,routineP,&
               "You need to define both VELOCITY and COORD and MASS and FORCE section (or none) in the NOSE section",&
               error=error)
          restart=explicit.and.restart
          work_section => section_vals_get_subs_vals(nose_section,"FORCE",error=error)
          CALL section_vals_get(work_section, explicit=explicit, error=error)
          CALL cp_assert(restart.or..not.explicit,cp_failure_level,cp_assertion_failed,routineP,&
               "You need to define both VELOCITY and COORD and MASS and FORCE section (or none) in the NOSE section",&
               error=error)
          restart=explicit.and.restart
       END IF

       IF (restart) THEN
          map_info => nhc%map_info
          CALL section_vals_val_get(nose_section,"COORD%_DEFAULT_KEYWORD_",r_vals=buffer,error=error)
          DO i = 1, SIZE ( nhc % nvt, 2)
             ind = map_info%index(i)
             ind = (ind-1) * nhc % nhc_len
             DO j = 1, SIZE(nhc % nvt, 1)
                ind = ind + 1
                nhc % nvt(j,i) % eta = buffer ( ind )
             END DO
          END DO
          CALL section_vals_val_get(nose_section,"VELOCITY%_DEFAULT_KEYWORD_",r_vals=buffer,error=error)
          DO i = 1, SIZE ( nhc % nvt, 2)
             ind = map_info%index(i)
             ind = (ind-1) * nhc % nhc_len
             DO j = 1, SIZE(nhc % nvt, 1)
                ind = ind + 1
                nhc % nvt(j,i) % v = buffer ( ind )
             END DO
          END DO
          CALL section_vals_val_get(nose_section,"MASS%_DEFAULT_KEYWORD_",r_vals=buffer,error=error)
          DO i = 1, SIZE ( nhc % nvt, 2)
             ind = map_info%index(i)
             ind = (ind-1) * nhc % nhc_len
             DO j = 1, SIZE(nhc % nvt, 1)
                ind = ind + 1
                nhc % nvt(j,i) % mass = buffer ( ind )
             END DO
          END DO
          CALL section_vals_val_get(nose_section,"FORCE%_DEFAULT_KEYWORD_",r_vals=buffer,error=error)
          DO i = 1, SIZE ( nhc % nvt, 2)
             ind = map_info%index(i)
             ind = (ind-1) * nhc % nhc_len
             DO j = 1, SIZE(nhc % nvt, 1)
                ind = ind + 1
                nhc % nvt(j,i) % f = buffer ( ind )
             END DO
          END DO
       END IF

       IF (save_mem) THEN
         NULLIFY(work_section)
         work_section => section_vals_get_subs_vals(nose_section,"COORD",error=error)
         CALL section_vals_remove_values(work_section, error)
         NULLIFY(work_section)
         work_section => section_vals_get_subs_vals(nose_section,"VELOCITY",error=error)
         CALL section_vals_remove_values(work_section, error)
         NULLIFY(work_section)
         work_section => section_vals_get_subs_vals(nose_section,"FORCE",error=error)
         CALL section_vals_remove_values(work_section, error)
         NULLIFY(work_section)
         work_section => section_vals_get_subs_vals(nose_section,"MASS",error=error)
         CALL section_vals_remove_values(work_section, error)
       END IF

     END IF

     CALL timestop(handle)

  END SUBROUTINE restart_nose

! *****************************************************************************
!> \brief Initializes the NHC velocities to the Maxwellian distribution
!> \par History
!>      none
!> \date 14-NOV-2000
! *****************************************************************************
  SUBROUTINE init_nhc_variables ( nhc,temp_ext, para_env, globenv, error )
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    REAL(KIND=dp), INTENT(IN)                :: temp_ext
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_nhc_variables', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, icount, j, number, &
                                                stat, tot_rn
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: akin, dum, temp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: array_of_rn
    TYPE(map_info_type), POINTER             :: map_info

    CALL timeset(routineN,handle)

    failure = .FALSE.

    ! first initializing the mass of the nhc variables
    nhc % nvt ( :, : ) % mass = nhc % nvt ( :, : ) % nkt * nhc%tau_nhc**2
    nhc % nvt ( :, : ) % eta = 0._dp
    nhc % nvt ( :, : ) % v = 0._dp
    nhc % nvt ( :, : ) % f = 0._dp

    map_info => nhc%map_info
    SELECT CASE (map_info%dis_type)
    CASE ( do_thermo_only_master )  ! for NPT
    CASE DEFAULT
       tot_rn = nhc%glob_num_nhc*nhc%nhc_len

       ALLOCATE(array_of_rn( tot_rn), STAT = stat)
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       array_of_rn(:) = 0.0_dp
    END SELECT

    SELECT CASE (map_info%dis_type )
    CASE ( do_thermo_only_master ) ! for NPT
       ! Map deterministically determined random number to nhc % v
       DO i = 1, nhc%loc_num_nhc
          DO j = 1, nhc % nhc_len
             nhc%nvt(j,i)%v = next_random_number(globenv%gaussian_rng_stream,error=error)
          END DO
       END DO
       
       akin = 0.0_dp
       DO i = 1, nhc%loc_num_nhc
          DO j = 1, nhc % nhc_len
             akin = akin + 0.5_dp * ( nhc % nvt ( j , i ) % mass * &
                  nhc % nvt ( j , i ) % v * &
                  nhc % nvt ( j , i ) % v )
          END DO
       END DO
       number = nhc%loc_num_nhc
       
       ! scale velocities to get the correct initial temperature
       temp = 2.0_dp*akin/REAL(number)
       temp = SQRT(temp_ext/temp)
       DO i = 1, nhc%loc_num_nhc
          DO j = 1, nhc % nhc_len
             nhc % nvt(j,i) % v = temp * nhc % nvt(j,i) % v
             nhc % nvt(j,i) % eta = 0.0_dp
          END DO
       END DO
       
       ! initializing all of the forces on the thermostats
       DO i = 1, nhc%loc_num_nhc
          DO j = 2, nhc % nhc_len
             nhc % nvt(j,i) % f = nhc % nvt(j-1,i) % mass*nhc % nvt(j-1,i) % v* &
                  nhc % nvt(j-1,i) % v - nhc % nvt(j,i) % nkt
             nhc % nvt(j,i) % f = nhc % nvt(j,i) % f/nhc % nvt(j,i) % mass
          END DO
       END DO
       
    CASE DEFAULT
       DO i=1,tot_rn
          array_of_rn(i) = next_random_number(globenv%gaussian_rng_stream,error=error)
       END DO
       ! Map deterministically determined random number to nhc % v
       DO i = 1, nhc%loc_num_nhc
          icount = map_info%index(i)
          icount = (icount-1) * nhc % nhc_len
          DO j = 1, nhc % nhc_len
             icount = icount + 1
             nhc % nvt(j,i) % v   = array_of_rn ( icount )
             nhc % nvt(j,i) % eta = 0.0_dp
          END DO
       END DO
       DEALLOCATE ( array_of_rn, STAT = stat )
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)

       number = nhc%glob_num_nhc
       CALL get_nhc_energies(nhc, dum, akin, para_env, error=error)

       ! scale velocities to get the correct initial temperature
       temp = 2.0_dp*akin/REAL(number)
       temp = SQRT(temp_ext/temp)
       DO i = 1, nhc%loc_num_nhc
          DO j = 1, nhc % nhc_len
             nhc % nvt(j,i) % v = temp * nhc % nvt(j,i) % v
          END DO
       END DO
       
       ! initializing all of the forces on the thermostats
       DO i = 1, nhc%loc_num_nhc
          DO j = 2, nhc % nhc_len
             nhc % nvt(j,i) % f = nhc % nvt(j-1,i) % mass*nhc % nvt(j-1,i) % v* &
                                  nhc % nvt(j-1,i) % v - nhc % nvt(j,i) % nkt
             nhc % nvt(j,i) % f = nhc % nvt(j,i) % f/nhc % nvt(j,i) % mass
          END DO
       END DO

    END SELECT

    CALL timestop(handle)

  END SUBROUTINE init_nhc_variables

! *****************************************************************************
!> \brief Initializes the barostat velocities to the Maxwellian distribution
!> \par History
!>      none
!> \date 14-NOV-2000
! *****************************************************************************
  SUBROUTINE init_barostat_variables ( npt, tau_cell, temp_ext, nfree, ensemble, &
       cmass, globenv, error )
    
    TYPE(npt_info_type), DIMENSION(:, :), &
      INTENT(INOUT)                          :: npt
    REAL(KIND=dp), INTENT(IN)                :: tau_cell, temp_ext
    INTEGER, INTENT(IN)                      :: nfree, ensemble
    REAL(KIND=dp), INTENT(IN)                :: cmass
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_barostat_variables', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, number
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: akin, temp, v

    CALL timeset(routineN,handle)

    failure = .FALSE.

    ! first initializing the mass of the nhc variables

    npt(:,:)%eps = 0.0_dp 
    npt(:,:)%v = 0.0_dp 
    npt(:,:)%f = 0.0_dp 
    SELECT CASE ( ensemble )
    CASE ( npt_i_ensemble )
       npt ( :, : ) % mass = REAL ( nfree + 3,KIND=dp) * temp_ext * tau_cell ** 2
    CASE ( npt_f_ensemble)
       npt ( :, : ) % mass = REAL ( nfree + 3,KIND=dp) * temp_ext * tau_cell ** 2 / 3.0_dp
    CASE ( nph_uniaxial_ensemble,nph_uniaxial_damped_ensemble )
       npt ( :, : ) % mass = cmass
    CASE ( npe_f_ensemble)
       npt ( :, : ) % mass =  REAL ( nfree + 3,KIND=dp) * temp_ext * tau_cell ** 2 / 3.0_dp
    CASE ( npe_i_ensemble )
       npt ( :, : ) % mass = REAL ( nfree + 3,KIND=dp) * temp_ext * tau_cell ** 2
    END SELECT
    ! initializing velocities
    DO i = 1, SIZE ( npt,1)
       DO j = i, SIZE ( npt,2)
          v = next_random_number(globenv%gaussian_rng_stream,error=error)
          ! Symmetrizing the initial barostat velocities to ensure
          ! no rotation of the cell under NPT_F
          npt(j,i) % v = v
          npt(i,j) % v = v
       END DO
    END DO
    
    akin = 0.0_dp
    DO i = 1, SIZE( npt,1 )
       DO j = 1, SIZE( npt,2 )
          akin = akin + 0.5_dp*( npt(j,i) % mass* npt(j,i)%v * npt(j,i)%v)
       END DO
    END DO
    
    number = SIZE ( npt, 1 ) * SIZE ( npt, 2 )
    
    ! scale velocities to get the correct initial temperature
    IF ( number /= 0 ) THEN
       temp = 2.0_dp * akin / REAL ( number,KIND=dp)
       temp = SQRT(temp_ext/temp)
    ENDIF
    DO i = 1, SIZE(npt,1)
       DO j = i, SIZE(npt,2)
          npt(j,i) % v = temp * npt(j,i) % v
          npt(i,j) % v = npt(j,i) % v
          IF (debug_isotropic_limit) THEN
             npt(j,i) % v = 0.0_dp
             npt(i,j) % v = 0.0_dp
             WRITE (*, *) 'DEBUG ISOTROPIC LIMIT| INTIAL v_eps', npt(j,i) % v
          END IF
       END DO
    END DO
    
    ! Zero barostat velocities for nph_uniaxial
    SELECT CASE ( ensemble )
    ! Zero barostat velocities for nph_uniaxial
    CASE ( nph_uniaxial_ensemble, nph_uniaxial_damped_ensemble )
       npt(:,:) % v = 0.0_dp
    END SELECT
    
    CALL timestop(handle)

  END SUBROUTINE init_barostat_variables

! *****************************************************************************
!> \brief Assigns extended parameters from the restart file.
!> \author CJM
! *****************************************************************************
  SUBROUTINE init_nhc_forces ( nhc, error)

    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_nhc_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(nhc),cp_fatal_level,routineP,error,failure)
    ! assign the forces
    DO i = 1, SIZE ( nhc % nvt, 2 )
       DO j = 2, SIZE ( nhc % nvt, 1 )
          nhc % nvt ( j, i ) % f = nhc % nvt ( j-1, i ) % mass * &
                                   nhc % nvt ( j-1, i ) % v **2 -  &
                                   nhc % nvt ( j, i ) % nkt
          nhc % nvt ( j, i ) % f = nhc % nvt ( j, i ) % f / nhc % nvt ( j, i ) % mass
       END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE init_nhc_forces

END MODULE extended_system_init
