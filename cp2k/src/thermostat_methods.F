!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Methods for Thermostats 
!> \author teo [tlaino] - University of Zurich - 10.2007
! *****************************************************************************
MODULE thermostat_methods
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE csvr_system_dynamics,            ONLY: csvr_barostat,&
                                             csvr_particles,&
                                             csvr_shells
  USE csvr_system_init,                ONLY: initialize_csvr_baro,&
                                             initialize_csvr_part,&
                                             initialize_csvr_shell
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE extended_system_dynamics,        ONLY: lnhc_barostat,&
                                             lnhc_particles,&
                                             lnhc_shells
  USE extended_system_init,            ONLY: initialize_nhc_baro,&
                                             initialize_nhc_part,&
                                             initialize_nhc_shell
  USE extended_system_types,           ONLY: npt_info_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE gle_system_dynamics,             ONLY: gle_particles,&
                                             initialize_gle_part
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       do_region_global, do_thermo_csvr, do_thermo_gle, do_thermo_nose, &
       do_thermo_same_as_part, npe_f_ensemble, npe_i_ensemble, &
       npt_f_ensemble, npt_i_ensemble, nve_ensemble, nvt_ensemble
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_remove_values,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: global_constraint_type,&
                                             molecule_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE qmmm_types,                      ONLY: qmmm_env_qm_type
  USE simpar_types,                    ONLY: simpar_type
  USE thermostat_types,                ONLY: allocate_thermostats,&
                                             create_thermostat_type,&
                                             release_thermostat_info,&
                                             release_thermostat_type,&
                                             release_thermostats,&
                                             thermostat_type,&
                                             thermostats_type
  USE thermostat_utils,                ONLY: compute_degrees_of_freedom,&
                                             get_thermostat_energies,&
                                             setup_thermostat_info
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: create_thermostats,&
            apply_thermostat_baro,&
            apply_thermostat_particles,&
            apply_thermostat_shells

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermostat_methods'

CONTAINS

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE create_thermostats(thermostats, md_section, force_env, simpar, &
       para_env, globenv, global_section, error )
    TYPE(thermostats_type), POINTER          :: thermostats
    TYPE(section_vals_type), POINTER         :: md_section
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: global_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_thermostats', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: binary_restart_file_name
    INTEGER                                  :: n_rep, region, thermostat_type
    LOGICAL :: apply_general_thermo, apply_thermo_baro, apply_thermo_shell, &
      explicit_baro, explicit_barostat_section, explicit_part, &
      explicit_shell, failure, save_mem, shell_adiabatic, shell_present
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(section_vals_type), POINTER :: barostat_section, print_section, &
      region_sections, thermo_baro_section, thermo_part_section, &
      thermo_shell_section, work_section

    NULLIFY(qmmm_env)
    failure =.FALSE.
    CALL allocate_thermostats(thermostats, error)
    thermo_part_section  => section_vals_get_subs_vals(md_section,"THERMOSTAT",error=error)
    thermo_shell_section => section_vals_get_subs_vals(md_section,"SHELL%THERMOSTAT",error=error)
    thermo_baro_section  => section_vals_get_subs_vals(md_section,"BAROSTAT%THERMOSTAT",error=error)
    barostat_section     => section_vals_get_subs_vals(md_section,"BAROSTAT",error=error)
    print_section => section_vals_get_subs_vals(md_section,"PRINT",error=error)

    CALL force_env_get(force_env, qmmm_env=qmmm_env, subsys=subsys, cell=cell, error=error )
    CALL section_vals_get(barostat_section,  explicit=explicit_barostat_section,  error=error)
    CALL section_vals_val_get(global_section,"SAVE_MEM",l_val=save_mem,  error=error)
    CALL section_vals_get(thermo_part_section,  explicit=explicit_part,  error=error)
    CALL section_vals_get(thermo_shell_section, explicit=explicit_shell, error=error)
    CALL section_vals_get(thermo_baro_section,  explicit=explicit_baro,  error=error)

    apply_thermo_baro = (simpar%ensemble == npt_f_ensemble).OR.&
                        (simpar%ensemble == npt_i_ensemble)

    apply_general_thermo = apply_thermo_baro .OR. (simpar%ensemble == nvt_ensemble)

    apply_thermo_shell = (simpar%ensemble == nve_ensemble)  .OR.&
                         (simpar%ensemble == nvt_ensemble)  .OR.&
                         (simpar%ensemble == npt_f_ensemble).OR.&
                         (simpar%ensemble == npt_i_ensemble).OR.&
                         (simpar%ensemble == npe_i_ensemble).OR.&
                         (simpar%ensemble == npe_f_ensemble)

    binary_restart_file_name = ""
    CALL section_vals_val_get(force_env%root_section,"EXT_RESTART%BINARY_RESTART_FILE_NAME",&
                              c_val=binary_restart_file_name,error=error)

    ! Compute Degrees of Freedom
    region = do_region_global
    region_sections => section_vals_get_subs_vals(thermo_part_section,"DEFINE_REGION",error=error)
    IF (explicit_part) CALL section_vals_val_get(thermo_part_section,"REGION",i_val=region,error=error)
    CALL cp_subsys_get(subsys, molecule_kinds_new=molecule_kinds_new, local_molecules_new=local_molecules,&
            molecules_new=molecules_new, gci=gci, particles=particles, error=error)
    CALL compute_degrees_of_freedom(thermostats, cell, simpar, molecule_kinds_new%els,&
         local_molecules, molecules_new, particles, print_section, region_sections, gci,&
         region, qmmm_env, error)

    ! Particles
    ! For constant temperature ensembles the thermostat is activated by default
    IF (explicit_part) THEN
       IF (apply_general_thermo) THEN
          CALL create_thermostat_type(thermostats%thermostat_part, simpar, thermo_part_section,&
               label="PARTICLES", error=error)
          ! Initialize thermostat
          IF (thermostats%thermostat_part%type_of_thermostat == do_thermo_nose) THEN
             ! Initialize or possibly restart Nose on Particles
             work_section => section_vals_get_subs_vals(thermo_part_section,"NOSE",error=error)
             CALL initialize_nhc_part(thermostats%thermostat_info_part, simpar, local_molecules,&
                  molecules_new%els, molecule_kinds_new%els, para_env, globenv,&
                  thermostats%thermostat_part%nhc, nose_section=work_section, gci=gci,&
                  save_mem=save_mem, binary_restart_file_name=binary_restart_file_name,&
                  error=error)
          ELSE IF (thermostats%thermostat_part%type_of_thermostat==do_thermo_csvr) THEN
             ! Initialize or possibly restart CSVR thermostat on Particles
             work_section => section_vals_get_subs_vals(thermo_part_section,"CSVR",error=error)
             CALL initialize_csvr_part(thermostats%thermostat_info_part,simpar,local_molecules,&
                  molecules_new%els, molecule_kinds_new%els, para_env, &
                  thermostats%thermostat_part%csvr, csvr_section=work_section,&
                  gci=gci, error=error)
          ELSE IF (thermostats%thermostat_part%type_of_thermostat==do_thermo_gle) THEN
             ! Initialize or possibly restart GLE thermostat on Particles
             work_section => section_vals_get_subs_vals(thermo_part_section,"GLE",error=error)
             CALL initialize_gle_part(thermostats%thermostat_info_part,simpar,local_molecules,&
                  molecules_new%els, molecule_kinds_new%els, particles%els, para_env, &
                  thermostats%thermostat_part%gle, gle_section=work_section,&
                  gci=gci,save_mem=save_mem, error=error)
          END IF
          CALL thermostat_info(thermostats%thermostat_part, "PARTICLES", thermo_part_section, &
                               simpar, para_env, error)
       ELSE
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
               "Thermostat for Particles has been defined but the ensemble provided "//&
               "does not support thermostat for Particles! Ignoring thermostat input."//&
CPSourceFileRef,&
               only_ionode=.TRUE.)          
       END IF
    ELSE IF (apply_general_thermo) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "One constant temperature ensemble has been required, but no thermostat for the "//&
            "particles has been defined. You may want to change your input and add a "//&
            "THERMOSTAT section in the MD section."//&
CPSourceFileRef,&
            only_ionode=.TRUE.) 
    END IF

    ! Core-Shell Model
    CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, error=error)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els,shell_present=shell_present)
    IF (shell_present) THEN
       IF (explicit_shell) THEN
          ! The thermostat is activated only if explicitely required
          ! It can be used to thermalize the shell-core motion when the temperature is not constant (nve, npe) 
          IF (apply_thermo_shell) THEN
             CALL create_thermostat_type(thermostats%thermostat_shell, simpar, thermo_shell_section,&
                  label="SHELL", error=error)
             CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els, shell_adiabatic=shell_adiabatic)
             region_sections => section_vals_get_subs_vals(thermo_shell_section,"DEFINE_REGION",error=error)
             CALL section_vals_val_get(thermo_shell_section,"REGION",i_val=region,error=error)
             CALL setup_thermostat_info(thermostats%thermostat_info_shell, molecule_kinds_new%els,&
                  local_molecules, molecules_new, particles, region, simpar%ensemble, shell=shell_adiabatic,&
                  region_sections=region_sections, qmmm_env=qmmm_env, error=error)
             IF (shell_adiabatic) THEN
                ! Initialize thermostat
                IF (thermostats%thermostat_shell%type_of_thermostat == do_thermo_nose) THEN
                   ! Initialize or possibly restart Nose on Shells
                   work_section => section_vals_get_subs_vals(thermo_shell_section,"NOSE",error=error)
                   CALL initialize_nhc_shell(thermostats%thermostat_info_shell, simpar, local_molecules,&
                        molecules_new%els, molecule_kinds_new%els, para_env, globenv,&
                        thermostats%thermostat_shell%nhc, nose_section=work_section,gci=gci,&
                        save_mem=save_mem, binary_restart_file_name=binary_restart_file_name,&
                        error=error)
                ELSE IF (thermostats%thermostat_shell%type_of_thermostat==do_thermo_csvr) THEN
                   ! Initialize or possibly restart CSVR thermostat on Shells
                   work_section => section_vals_get_subs_vals(thermo_shell_section,"CSVR",error=error)
                   CALL initialize_csvr_shell(thermostats%thermostat_info_shell, simpar, local_molecules,&
                        molecules_new%els, molecule_kinds_new%els, para_env, &
                        thermostats%thermostat_shell%csvr, csvr_section=work_section, gci=gci,&
                        error=error)               
                END IF
                CALL thermostat_info(thermostats%thermostat_shell, "CORE-SHELL", thermo_shell_section, &
                                     simpar, para_env, error)
             ELSE
                CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                     "Thermostat for Core-Shell motion only with adiabatic shell-model. "//&
                     "Continuing calculation ignoring the thermostat info! No Thermostat "//&
                     "applied to Shells!"//&
CPSourceFileRef,&
                     only_ionode=.TRUE.) 
                CALL release_thermostat_type(thermostats%thermostat_shell, error)
                CALL release_thermostat_info(thermostats%thermostat_info_shell, error)
             END IF
          ELSE
             CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
             "Thermostat for Shells has been defined but for the selected ensemble the adiabatic  "//&
             " shell model has not been implemented! Ignoring thermostat input."//&
CPSourceFileRef,&
             only_ionode=.TRUE.)          
          END IF
       END IF
    ELSE IF (explicit_shell) THEN
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
            "Thermostat for Shells has been defined but the system provided "//&
            "does not contain any Shells! Ignoring thermostat input."//&
CPSourceFileRef,&
            only_ionode=.TRUE.) 
    END IF 

    ! Barostat Temperature (not necessarily to be controlled by a thermostat)
    IF (explicit_barostat_section) THEN
       simpar%temp_baro_ext = simpar%temp_ext 
       CALL section_vals_val_get(md_section,"BAROSTAT%TEMPERATURE",n_rep_val=n_rep,error=error)
       IF (n_rep/=0) THEN
          CALL section_vals_val_get(md_section,"BAROSTAT%TEMPERATURE",r_val=simpar%temp_baro_ext,error=error)
          CPPostcondition(simpar%temp_baro_ext>=0.0_dp,cp_failure_level,routineP,error,failure)
       END IF

       ! Setup Barostat Thermostat 
       IF (apply_thermo_baro) THEN
          ! Check if we use the same thermostat as particles
          CALL section_vals_val_get(thermo_baro_section,"TYPE",i_val=thermostat_type,error=error)
          work_section => thermo_baro_section
          IF (thermostat_type==do_thermo_same_as_part) work_section => thermo_part_section
         
          CALL create_thermostat_type(thermostats%thermostat_baro, simpar, work_section, skip_region=.TRUE.,&
               label="BAROSTAT", error=error)
          ! Initialize thermostat
          IF      (thermostats%thermostat_baro%type_of_thermostat==do_thermo_nose) THEN
             ! Initialize or possibly restart Nose on Barostat
             work_section => section_vals_get_subs_vals(thermo_baro_section,"NOSE",error=error)
             CALL initialize_nhc_baro(simpar, para_env, globenv, thermostats%thermostat_baro%nhc,&
                  nose_section=work_section, save_mem=save_mem, error=error)
          ELSE IF (thermostats%thermostat_baro%type_of_thermostat==do_thermo_csvr) THEN
             ! Initialize or possibly restart CSVR thermostat on Barostat
             work_section => section_vals_get_subs_vals(thermo_baro_section,"CSVR",error=error)
             CALL initialize_csvr_baro(simpar, thermostats%thermostat_baro%csvr,&
                  csvr_section=work_section, error=error)
          END IF
          CALL thermostat_info(thermostats%thermostat_baro, "BAROSTAT", thermo_baro_section,&
                               simpar, para_env, error)
          ! If thermostat for barostat uses a diffent kind than the one of the particles
          ! let's update infos in the input structure..
          IF (thermostat_type==do_thermo_same_as_part) THEN
             CALL update_thermo_baro_section(thermostats%thermostat_baro, thermo_baro_section, error)
          END IF
       ELSE
          IF (explicit_baro) THEN
             CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                  "Thermostat for Barostat has been defined but the ensemble provided "//&
                  "does not support thermostat for Barostat! Ignoring thermostat input."//&
CPSourceFileRef,&
                  only_ionode=.TRUE.)
          END IF
          ! Let's remove the section and set it as required.
          CALL section_vals_remove_values(thermo_baro_section, error)
          thermo_baro_section%section%required = .TRUE.
       END IF
    END IF

    ! Release the thermostats info..
    CALL release_thermostat_info(thermostats%thermostat_info_part, error)
    CALL release_thermostat_info(thermostats%thermostat_info_shell, error)

    ! If no thermostats have been allocated deallocate the full structure
    IF(  (.NOT.ASSOCIATED(thermostats%thermostat_part)).AND.&
         (.NOT.ASSOCIATED(thermostats%thermostat_shell)).AND.&
         (.NOT.ASSOCIATED(thermostats%thermostat_baro)) ) THEN
       CALL release_thermostats(thermostats, error)
    END IF

  END SUBROUTINE create_thermostats

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_thermo_baro_section(thermostat, section, error)
    TYPE(thermostat_type), POINTER           :: thermostat
    TYPE(section_vals_type), POINTER         :: section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_thermo_baro_section', &
      routineP = moduleN//':'//routineN

    TYPE(section_vals_type), POINTER         :: work_section

    CALL section_vals_val_set(section,"TYPE",i_val=thermostat%type_of_thermostat,error=error)
    SELECT CASE(thermostat%type_of_thermostat)
    CASE(do_thermo_nose)
       work_section => section_vals_get_subs_vals(section,"NOSE",error=error)
       CALL section_vals_val_set(work_section,"LENGTH",i_val=thermostat%nhc%nhc_len,error=error)
       CALL section_vals_val_set(work_section,"YOSHIDA",i_val=thermostat%nhc%nyosh,error=error)
       CALL section_vals_val_set(work_section,"TIMECON",r_val=thermostat%nhc%tau_nhc,error=error)
       CALL section_vals_val_set(work_section,"MTS",i_val=thermostat%nhc%nc,error=error)
    CASE(do_thermo_csvr)
       work_section => section_vals_get_subs_vals(section,"CSVR",error=error)
       CALL section_vals_val_set(work_section,"TIMECON",r_val=thermostat%csvr%tau_csvr,error=error)
    END SELECT

  END SUBROUTINE update_thermo_baro_section

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE thermostat_info(thermostat, label, section,  simpar, para_env, error)
    TYPE(thermostat_type), POINTER           :: thermostat
    CHARACTER(LEN=*), INTENT(IN)             :: label
    TYPE(section_vals_type), POINTER         :: section
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'thermostat_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iw
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: kin_energy, pot_energy, tmp
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw=cp_print_key_unit_nr(logger,section,"PRINT%THERMOSTAT_INFO",extension=".log",error=error)
    ! Total Tehrmostat Energy
    CALL get_thermostat_energies(thermostat, pot_energy, kin_energy, para_env, error=error)
    IF (iw>0) THEN
       WRITE ( iw, '( /,A )' ) ' THERMOSTAT| Thermostat Info for '//TRIM(label)
       SELECT CASE(thermostat%type_of_thermostat)
       CASE(do_thermo_nose)
          WRITE ( iw, '( A,T63,A )' ) ' THERMOSTAT| Type of thermostat','Nose-Hoover-Chains'
          WRITE ( iw, '( A,T77,I4 )' ) ' THERMOSTAT| Nose-Hoover-Chain length ', thermostat%nhc%nhc_len
          tmp = cp_unit_from_cp2k(thermostat%nhc%tau_nhc,"fs",error=error)
          WRITE ( iw, '( A,A4,A,T71,F10.2 )' ) &
               ' THERMOSTAT| Nose-Hoover-Chain time constant [', 'fs' ,'] ',tmp
          WRITE ( iw, '( A,T77,I4 )' ) ' THERMOSTAT| Order of Yoshida integrator ', &
               thermostat%nhc%nyosh
          WRITE ( iw, '( A,T77,I4 )' ) ' THERMOSTAT| Number of multiple time steps', &
               thermostat%nhc%nc
          WRITE ( iw, '( A,T69,F12.6 )' ) ' THERMOSTAT| Initial Potential Energy', &
               pot_energy
          WRITE ( iw, '( A,T69,F12.6 )' ) ' THERMOSTAT| Initial Kinetic Energy', &
               kin_energy
       CASE(do_thermo_csvr)
          WRITE ( iw, '( A,T44,A )' ) ' THERMOSTAT| Type of thermostat','Canonical Sampling/Velocity Rescaling'
          tmp = cp_unit_from_cp2k(thermostat%csvr%tau_csvr,"fs",error=error)*0.5_dp*simpar%dt
          WRITE ( iw, '( A,A4,A,T71,F10.2 )' ) &
               ' THERMOSTAT| CSVR time constant [', 'fs' ,'] ',tmp
          WRITE ( iw, '( A,T69,F12.6 )' ) &
               ' THERMOSTAT| Initial Kinetic Energy ',kin_energy
       END SELECT
       WRITE ( iw, '( A,/ )' ) ' THERMOSTAT| End of Thermostat Info for '//TRIM(label)
    END IF
    CALL cp_print_key_finished_output(iw,logger,section,"PRINT%THERMOSTAT_INFO",error=error)

  END SUBROUTINE thermostat_info

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE apply_thermostat_baro(thermostat, npt, group, error)
    TYPE(thermostat_type), POINTER           :: thermostat
    TYPE(npt_info_type), DIMENSION(:, :), &
      POINTER                                :: npt
    INTEGER, INTENT(IN)                      :: group
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_thermostat_baro', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(thermostat)) THEN
       IF      (thermostat%type_of_thermostat==do_thermo_nose) THEN
          ! Apply Nose-Hoover Thermostat
          CPPostcondition(ASSOCIATED(thermostat%nhc),cp_failure_level,routineP,error,failure)
          CALL lnhc_barostat ( thermostat%nhc, npt, group, error )
       ELSE IF (thermostat%type_of_thermostat==do_thermo_csvr) THEN  
          ! Apply CSVR Thermostat
          CPPostcondition(ASSOCIATED(thermostat%csvr),cp_failure_level,routineP,error,failure)
          CALL csvr_barostat ( thermostat%csvr, npt, group, error )
       END IF
    END IF
  END SUBROUTINE apply_thermostat_baro

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE apply_thermostat_particles(thermostat, molecule_kind_set, molecule_set,&
       particle_set, local_molecules, group, shell_adiabatic, shell_particle_set,&
       core_particle_set, vel, shell_vel, core_vel, error)

    TYPE(thermostat_type), POINTER           :: thermostat
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    INTEGER, INTENT(IN)                      :: group
    LOGICAL, INTENT(IN), OPTIONAL            :: shell_adiabatic
    TYPE(particle_type), OPTIONAL, POINTER   :: shell_particle_set( : ), &
                                                core_particle_set( : )
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: vel(:,:), shell_vel(:,:), &
                                                core_vel(:,:)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_thermostat_particles', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(thermostat)) THEN
       IF      (thermostat%type_of_thermostat==do_thermo_nose) THEN
          ! Apply Nose-Hoover Thermostat
          CPPostcondition(ASSOCIATED(thermostat%nhc),cp_failure_level,routineP,error,failure)
          CALL lnhc_particles(thermostat%nhc,molecule_kind_set,molecule_set,&
               particle_set,local_molecules,group,shell_adiabatic,shell_particle_set,&
               core_particle_set, vel, shell_vel, core_vel, error)
       ELSE IF (thermostat%type_of_thermostat==do_thermo_csvr) THEN  
          ! Apply CSVR Thermostat
          CPPostcondition(ASSOCIATED(thermostat%csvr),cp_failure_level,routineP,error,failure)
          CALL csvr_particles(thermostat%csvr,molecule_kind_set,molecule_set,&
               particle_set,local_molecules,group,shell_adiabatic,shell_particle_set,&
               core_particle_set, vel, shell_vel, core_vel, error)
       ELSE IF (thermostat%type_of_thermostat==do_thermo_gle) THEN
         ! Apply GLE Thermostat
          CPPostcondition(ASSOCIATED(thermostat%gle),cp_failure_level,routineP,error,failure)
          CALL gle_particles(thermostat%gle,molecule_kind_set,molecule_set,&
               particle_set,local_molecules,group,shell_adiabatic,shell_particle_set,&
               core_particle_set, vel, shell_vel, core_vel, error=error)
       END IF
    END IF
  END SUBROUTINE apply_thermostat_particles

! *****************************************************************************
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE apply_thermostat_shells(thermostat,atomic_kind_set, particle_set,&
       local_particles, group, shell_particle_set, core_particle_set, vel, shell_vel,&
       core_vel, error)

    TYPE(thermostat_type), POINTER           :: thermostat
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    INTEGER, INTENT(IN)                      :: group
    TYPE(particle_type), OPTIONAL, POINTER   :: shell_particle_set(:), &
                                                core_particle_set(:)
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: vel(:,:), shell_vel(:,:), &
                                                core_vel(:,:)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_thermostat_shells', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(thermostat)) THEN
       IF      (thermostat%type_of_thermostat==do_thermo_nose) THEN
          ! Apply Nose-Hoover Thermostat
          CPPostcondition(ASSOCIATED(thermostat%nhc),cp_failure_level,routineP,error,failure)
          CALL lnhc_shells(thermostat%nhc, atomic_kind_set, particle_set, local_particles, &
               group, shell_particle_set, core_particle_set, vel, shell_vel, core_vel, error)
       ELSE IF (thermostat%type_of_thermostat==do_thermo_csvr) THEN  
          ! Apply CSVR Thermostat
          CPPostcondition(ASSOCIATED(thermostat%csvr),cp_failure_level,routineP,error,failure)
          CALL csvr_shells(thermostat%csvr, atomic_kind_set, particle_set, local_particles, &
               group, shell_particle_set, core_particle_set, vel, shell_vel, core_vel, error)
       END IF
    END IF
  END SUBROUTINE apply_thermostat_shells

END MODULE thermostat_methods
