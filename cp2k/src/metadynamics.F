!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Performs the metadynamics calculation
!> \par History
!>      01.2005 created [fawzi and ale]
!>      11.2007 Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
MODULE metadynamics
  USE bibliography,                    ONLY: VandenCic2006
  USE colvar_methods,                  ONLY: colvar_eval_glob_f
  USE colvar_types,                    ONLY: colvar_p_type,&
                                             torsion_colvar_id
  USE constraint_fxd,                  ONLY: fix_atom_control
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_constants,                 ONLY: do_wall_m,&
                                             do_wall_p,&
                                             do_wall_reflective
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE metadynamics_types,              ONLY: hills_env_type,&
                                             meta_env_type,&
                                             metavar_type,&
                                             multiple_walkers_type
  USE metadynamics_utils,              ONLY: add_hill_single,&
                                             get_meta_iter_level,&
                                             meta_walls,&
                                             restart_hills,&
                                             synchronize_multiple_walkers
  USE parallel_rng_types,              ONLY: next_random_number
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: boltzmann,&
                                             femtoseconds,&
                                             joule,&
                                             kelvin
  USE reference_manager,               ONLY: cite_reference
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE 

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'metadynamics'
  INTEGER, SAVE, PRIVATE :: last_meta_env_id=0

  PUBLIC :: metadyn_forces, metadyn_integrator
  PUBLIC :: metadyn_velocities_colvar, metadyn_write_colvar

CONTAINS
! *****************************************************************************
!> \brief  General driver for applying metadynamics
!> \author Teodoro Laino 
!> \date   01.2009
!> \par History
!>      01.2009 created
! *****************************************************************************
  SUBROUTINE metadyn_integrator(force_env, vel, rand, error)
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: vel
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: rand
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_integrator', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure = .FALSE.
    ! Apply Metadynamics
    IF (ASSOCIATED(force_env%meta_env)) THEN
       IF(force_env%meta_env%langevin) THEN
          IF (.NOT.PRESENT(rand)) THEN
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Langevin on COLVAR not implemented for this MD ensemble!", &
                  error=error, error_level=cp_fatal_level)
          END IF
          !    *** Velocity Verlet for Langevin S(t)->S(t+1) 
          CALL metadyn_position_colvar(force_env, error=error)
          !    *** Forces from Vs and  S(X(t+1))
          CALL metadyn_forces(force_env, error=error)
          !    *** Velocity Verlet for Langeving *** v(t+1/2)--> v(t)
          CALL metadyn_velocities_colvar(force_env, rand,error=error)
       ELSE
          CALL metadyn_forces(force_env, vel,error=error)
       ENDIF
       !    *** Write down COVAR informations
       CALL metadyn_write_colvar(force_env, error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE metadyn_integrator

! *****************************************************************************
!> \brief add forces to the subsys due to the metadynamics run
!>      possibly modifies the velocites (if reflective walls are applied)
!> \par History
!>      04.2004 created
! *****************************************************************************
  SUBROUTINE metadyn_forces(force_env,vel,error)
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: vel
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, i_c, icolvar, ii, &
                                                iwall
    LOGICAL                                  :: explicit, failure
    REAL(kind=dp)                            :: check_val, diff_ss, dt, &
                                                ekin_w, fac_t, fft, norm, &
                                                rval, scal, scalf, ss0_test, &
                                                tol_ekin
    TYPE(colvar_p_type), DIMENSION(:), &
      POINTER                                :: colvar_p
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(metavar_type), POINTER              :: cv
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: ss0_section, vvp_section

    failure=.FALSE.
    NULLIFY(logger,meta_env)
    meta_env => force_env%meta_env
    IF (.NOT.ASSOCIATED(meta_env)) RETURN

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)
    NULLIFY(colvar_p,subsys,cv,ss0_section, vvp_section)
    CALL force_env_get(force_env, subsys=subsys, error=error)

    dt = meta_env%dt
    IF (.NOT.meta_env%restart) meta_env%n_steps=meta_env%n_steps+1

    ! Initialize velocity
    IF (meta_env%restart.AND.meta_env%extended_lagrange) THEN
       meta_env%ekin_s = 0.0_dp
       DO i_c=1,meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          cv%vvp = next_random_number(force_env%globenv%gaussian_rng_stream,error=error)
          meta_env%ekin_s = meta_env%ekin_s + 0.5_dp*cv%mass*cv%vvp**2
       END DO
       ekin_w = 0.5_dp*meta_env%temp_wanted*REAL(meta_env%n_colvar,KIND=dp)
       fac_t  = SQRT(ekin_w/MAX(meta_env%ekin_s,1.0E-8_dp))
       DO i_c = 1, meta_env%n_colvar
          cv     => meta_env%metavar(i_c)
          cv%vvp = cv%vvp*fac_t
       ENDDO
       meta_env%ekin_s = 0.0_dp
    END IF

          !    *** Velocity Verlet for Langevin S(t)->S(t+1)
    ! compute ss and the derivative of ss with respect to the atomic positions
    DO i_c=1,meta_env%n_colvar
       cv => meta_env%metavar(i_c)
       icolvar = cv%icolvar
       CALL colvar_eval_glob_f(icolvar,force_env,error=error)
       cv%ss       = subsys%colvar_p(icolvar)%colvar%ss

       ! Setup the periodic flag if the COLVAR is (-pi,pi] periodic
       cv%periodic = (subsys%colvar_p(icolvar)%colvar%type_id==torsion_colvar_id)

       ! Restart for Extended Lagrangian Metadynamics
       IF (meta_env%restart) THEN
          ! Initialize the position of the collective variable in the extended lagrange
          ss0_section => section_vals_get_subs_vals(meta_env%metadyn_section,"EXT_LAGRANGE_SS0",error=error)
          CALL section_vals_get(ss0_section, explicit=explicit, error=error)
          IF (explicit) THEN
             CALL section_vals_val_get(ss0_section,"_DEFAULT_KEYWORD_",&
                  i_rep_val=i_c, r_val=rval, error=error)
             cv%ss0 = rval
          ELSE
             cv%ss0 = cv%ss
          END IF
          vvp_section => section_vals_get_subs_vals(meta_env%metadyn_section,"EXT_LAGRANGE_VVP",error=error)
          CALL section_vals_get(vvp_section, explicit=explicit, error=error)
          IF (explicit) THEN
             CALL section_vals_val_get(vvp_section,"_DEFAULT_KEYWORD_",&
                  i_rep_val=i_c, r_val=rval, error=error)
             cv%vvp = rval
          END IF
       END IF
       !
       IF(.NOT.meta_env%extended_lagrange) THEN 
          cv%ss0 = cv%ss
          cv%vvp = 0.0_dp
       END IF
    ENDDO
    ! History dependent forces (evaluated at s0)
    IF(meta_env%do_hills) CALL  hills(meta_env,error)
    
    ! Apply walls to the colvars
    CALL meta_walls(meta_env, error)

    meta_env%restart = .FALSE.
    IF(.NOT.meta_env%extended_lagrange)THEN
       meta_env%ekin_s = 0.0_dp
       meta_env%epot_s = 0.0_dp
       meta_env%epot_walls = 0.0_dp
       DO i_c = 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          cv%epot_s = 0.0_dp
          cv%ff_s   = 0.0_dp
          meta_env%epot_walls = meta_env%epot_walls + cv%epot_walls
          icolvar=cv%icolvar
          NULLIFY(particles)
          CALL cp_subsys_get(subsys, colvar_p=colvar_p, &
               particles=particles,error=error)
          DO ii=1,colvar_p(icolvar)%colvar%n_atom_s
             i=colvar_p(icolvar)%colvar%i_atom(ii)
             fft = cv%ff_hills + cv%ff_walls
             particles%els(i)%f=particles%els(i)%f+fft*colvar_p(icolvar)%colvar%dsdr(:,ii)
          ENDDO
       ENDDO
    ELSE
       meta_env%ekin_s = 0.0_dp
       meta_env%epot_s = 0.0_dp
       meta_env%epot_walls = 0.0_dp
       DO i_c= 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          diff_ss = cv%ss-cv%ss0
          IF (cv%periodic) THEN
             ! The difference of a periodic COLVAR is always within [-pi,pi]
             diff_ss = SIGN(1.0_dp,ASIN(SIN(diff_ss)))*ACOS(COS(diff_ss))
          END IF
          cv%epot_s = 0.5_dp*cv%lambda*(diff_ss)**2.0_dp
          cv%ff_s   = cv%lambda*(diff_ss)
          icolvar   = cv%icolvar
          ! forces on the atoms
          NULLIFY(particles)
          CALL cp_subsys_get(subsys, colvar_p=colvar_p, &
               particles=particles,error=error)
          DO ii=1,colvar_p(icolvar)%colvar%n_atom_s
             i=colvar_p(icolvar)%colvar%i_atom(ii)
             particles%els(i)%f=particles%els(i)%f- cv%ff_s*colvar_p(icolvar)%colvar%dsdr(:,ii)
          ENDDO
          !  velocity verlet on the s0 if NOT langevin
          IF (.NOT.meta_env%langevin) THEN
             fft    = cv%ff_s + cv%ff_hills + cv%ff_walls
             cv%vvp = cv%vvp+dt*fft/cv%mass
             meta_env%ekin_s = meta_env%ekin_s + 0.5_dp*cv%mass*cv%vvp**2
             meta_env%epot_s = meta_env%epot_s + cv%epot_s
             meta_env%epot_walls = meta_env%epot_walls + cv%epot_walls
          END IF
       ENDDO
       !  velocity rescaling on the s0
       IF(meta_env%tempcontrol.AND.(.NOT.meta_env%langevin))THEN
          ekin_w   = 0.5_dp*meta_env%temp_wanted*REAL(meta_env%n_colvar,KIND=dp)
          tol_ekin = 0.5_dp*meta_env%toll_temp*REAL(meta_env%n_colvar,KIND=dp)
          IF(ABS(ekin_w-meta_env%ekin_s)> tol_ekin)THEN
             fac_t=SQRT(ekin_w/MAX(meta_env%ekin_s,1.0E-8_dp))
             DO i_c = 1, meta_env%n_colvar
                cv     => meta_env%metavar(i_c)
                cv%vvp = cv%vvp*fac_t
             ENDDO
             meta_env%ekin_s = ekin_w
          ENDIF
       ENDIF
       ! Reflective Wall only for s0
       DO i_c = 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          IF (cv%do_wall) THEN
             DO iwall = 1, SIZE(cv%walls)
                SELECT CASE(cv%walls(iwall)%id_type)
                CASE(do_wall_reflective)
                   ss0_test=cv%ss0+dt*cv%vvp
                   IF (cv%periodic) THEN
                      ! A periodic COLVAR is always within [-pi,pi]
                      ss0_test = SIGN(1.0_dp,ASIN(SIN(ss0_test)))*ACOS(COS(ss0_test))
                   END IF
                   SELECT CASE(cv%walls(iwall)%id_direction)
                   CASE(do_wall_p)
                      IF ((ss0_test > cv%walls(iwall)%pos).AND.(cv%vvp>0)) cv%vvp = -cv%vvp
                   CASE(do_wall_m)
                      IF ((ss0_test < cv%walls(iwall)%pos).AND.(cv%vvp<0)) cv%vvp = -cv%vvp
                   END SELECT
                END SELECT
             END DO
          ENDIF
       ENDDO
       ! Update of ss0 if NOT langevin
       IF (.NOT.meta_env%langevin) THEN
          DO i_c = 1, meta_env%n_colvar
             cv => meta_env%metavar(i_c)
             cv%ss0 = cv%ss0 + dt*cv%vvp
             IF (cv%periodic) THEN
                ! A periodic COLVAR is always within [-pi,pi]
                cv%ss0 = SIGN(1.0_dp,ASIN(SIN(cv%ss0)))*ACOS(COS(cv%ss0))
             END IF
          ENDDO
       END IF
    ENDIF
    ! Constraints ONLY of Fixed Atom type
    CALL fix_atom_control(force_env, error=error)

    ! Reflective Wall only for ss
    DO i_c = 1, meta_env%n_colvar
       cv => meta_env%metavar(i_c)
       IF (cv%do_wall) THEN
          DO iwall = 1, SIZE(cv%walls)
             SELECT CASE(cv%walls(iwall)%id_type)
             CASE(do_wall_reflective)
                SELECT CASE(cv%walls(iwall)%id_direction)
                CASE(do_wall_p)
                   IF (cv%ss < cv%walls(iwall)%pos) CYCLE
                   check_val = -1.0_dp
                CASE(do_wall_m)
                   IF (cv%ss > cv%walls(iwall)%pos) CYCLE
                   check_val =  1.0_dp
                END SELECT
                NULLIFY(particles)
                icolvar=cv%icolvar
                CALL cp_subsys_get(subsys, colvar_p=colvar_p, particles=particles,error=error)
                scal=0.0_dp
                scalf=0.0_dp
                norm=0.0_dp
                DO ii=1,colvar_p(icolvar)%colvar%n_atom_s
                   i=colvar_p(icolvar)%colvar%i_atom(ii)
                   IF (PRESENT(vel)) THEN
                      scal=scal+vel(1,i)*colvar_p(icolvar)%colvar%dsdr(1,ii)
                      scal=scal+vel(2,i)*colvar_p(icolvar)%colvar%dsdr(2,ii)
                      scal=scal+vel(3,i)*colvar_p(icolvar)%colvar%dsdr(3,ii)                      
                   ELSE
                      scal=scal+particles%els(i)%v(1)*colvar_p(icolvar)%colvar%dsdr(1,ii)
                      scal=scal+particles%els(i)%v(2)*colvar_p(icolvar)%colvar%dsdr(2,ii)
                      scal=scal+particles%els(i)%v(3)*colvar_p(icolvar)%colvar%dsdr(3,ii)
                   END IF
                   scalf=scalf+particles%els(i)%f(1)*colvar_p(icolvar)%colvar%dsdr(1,ii)
                   scalf=scalf+particles%els(i)%f(2)*colvar_p(icolvar)%colvar%dsdr(2,ii)
                   scalf=scalf+particles%els(i)%f(3)*colvar_p(icolvar)%colvar%dsdr(3,ii)
                   norm=norm+colvar_p(icolvar)%colvar%dsdr(1,ii)**2
                   norm=norm+colvar_p(icolvar)%colvar%dsdr(2,ii)**2
                   norm=norm+colvar_p(icolvar)%colvar%dsdr(3,ii)**2
                ENDDO
                IF (norm /= 0.0_dp) scal = scal/norm
                IF (norm /= 0.0_dp) scalf= scalf/norm

                IF (scal*check_val>0.0_dp) CYCLE
                DO ii=1,colvar_p(icolvar)%colvar%n_atom_s
                   i=colvar_p(icolvar)%colvar%i_atom(ii)
                   IF (PRESENT(vel)) THEN
                      vel(:,i)=vel(:,i)-2.0_dp*colvar_p(icolvar)%colvar%dsdr(:,ii)*scal
                   ELSE
                      particles%els(i)%v(:)=particles%els(i)%v(:)-2.0_dp*colvar_p(icolvar)%colvar%dsdr(:,ii)*scal
                   END IF
                   ! Nullify forces along the colvar (this avoids the weird behaviors of the reflective wall)
                    particles%els(i)%f(:)=particles%els(i)%f(:)-colvar_p(icolvar)%colvar%dsdr(:,ii)*scalf
                ENDDO
             END SELECT
          END DO
       END IF
    ENDDO

    CALL timestop(handle)
  END SUBROUTINE metadyn_forces

! *****************************************************************************
!> \brief Evolves velocities COLVAR according to 
!>        Vanden-Eijnden Ciccotti C.Phys.Letter 429 (2006) 310-316
!> \author Fabio Sterpone and Teodoro Laino
!> \date  01.2009
! *****************************************************************************
  SUBROUTINE metadyn_velocities_colvar(force_env,rand,error)
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT), OPTIONAL                :: rand
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_velocities_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i_c
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: diff_ss, dt, fft, sigma
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(metavar_type), POINTER              :: cv

    failure=.FALSE.
    NULLIFY(logger,meta_env,cv)
    meta_env => force_env%meta_env
    IF (.NOT.ASSOCIATED(meta_env)) RETURN

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)
    ! Add citation
    IF (meta_env%langevin) CALL cite_reference(VandenCic2006)

    dt = meta_env%dt
    ! History dependent forces (evaluated at s0)
    IF(meta_env%do_hills) CALL  hills(meta_env,error)

    ! Evolve Velocities
    meta_env%ekin_s = 0.0_dp
    meta_env%epot_walls = 0.0_dp
    DO i_c= 1, meta_env%n_colvar
       cv => meta_env%metavar(i_c)
       diff_ss = cv%ss-cv%ss0
       IF (cv%periodic) THEN
          ! The difference of a periodic COLVAR is always within [-pi,pi]
          diff_ss = SIGN(1.0_dp,ASIN(SIN(diff_ss)))*ACOS(COS(diff_ss))
       END IF
       cv%epot_s = 0.5_dp*cv%lambda*(diff_ss)**2.0_dp
       cv%ff_s   = cv%lambda*(diff_ss)
       
       fft    = cv%ff_s+cv%ff_hills
       sigma  = SQRT((meta_env%temp_wanted*kelvin)*2.0_dp*(boltzmann/joule)*cv%gamma/cv%mass)
       cv%vvp = cv%vvp+0.5_dp*dt*fft/cv%mass-0.5_dp*dt*cv%gamma*cv%vvp+&
                0.5_dp*SQRT(dt)*sigma*rand(i_c)
       meta_env%ekin_s = meta_env%ekin_s + 0.5_dp*cv%mass*cv%vvp**2
       meta_env%epot_walls = meta_env%epot_walls + cv%epot_walls
    ENDDO
    CALL timestop(handle)

  END SUBROUTINE metadyn_velocities_colvar

! *****************************************************************************
!> \brief Evolves COLVAR position
!>        Vanden-Eijnden Ciccotti C.Phys.Letter 429 (2006) 310-316
!> \author Fabio Sterpone and Teodoro Laino
!> \date  01.2009
! *****************************************************************************
  SUBROUTINE metadyn_position_colvar(force_env,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_position_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i_c
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: dt
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(metavar_type), POINTER              :: cv

    failure=.FALSE.
    NULLIFY(logger,meta_env,cv)
    meta_env => force_env%meta_env
    IF (.NOT.ASSOCIATED(meta_env)) RETURN

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)

    ! Add citation
    IF (meta_env%langevin) CALL cite_reference(VandenCic2006)
    dt = meta_env%dt

    ! Update of ss0
    DO i_c = 1, meta_env%n_colvar
       cv => meta_env%metavar(i_c)             
       cv%ss0 =cv%ss0 +dt*cv%vvp
       IF (cv%periodic) THEN
          ! A periodic COLVAR is always within [-pi,pi]
          cv%ss0 = SIGN(1.0_dp,ASIN(SIN(cv%ss0)))*ACOS(COS(cv%ss0))
       END IF
    ENDDO
    CALL timestop(handle)

  END SUBROUTINE metadyn_position_colvar

! *****************************************************************************
!> \brief Write down COLVAR information evolved  according to 
!>        Vanden-Eijnden Ciccotti C.Phys.Letter 429 (2006) 310-316
!> \author Fabio Sterpone and Teodoro Laino
!> \date  01.2009
! *****************************************************************************
  SUBROUTINE metadyn_write_colvar(force_env,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_write_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, i_c, iw
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: diff_ss, temp
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(metavar_type), POINTER              :: cv

    failure=.FALSE.
    NULLIFY(logger,meta_env,cv)
    meta_env => force_env%meta_env
    IF (.NOT.ASSOCIATED(meta_env)) RETURN

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)

    ! If Langevin we need to recompute few quantities
    ! This does not apply to the standard lagrangian scheme since it is
    ! implemented with a plain Newton integration scheme.. while Langevin
    ! follows the correct Verlet integration.. This will have to be made
    ! uniform in the future (Teodoro Laino - 01.2009)
    IF (meta_env%langevin) THEN
       meta_env%ekin_s = 0.0_dp
       meta_env%epot_s = 0.0_dp
       DO i_c= 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          diff_ss = cv%ss-cv%ss0
          IF (cv%periodic) THEN
             ! The difference of a periodic COLVAR is always within [-pi,pi]
             diff_ss = SIGN(1.0_dp,ASIN(SIN(diff_ss)))*ACOS(COS(diff_ss))
          END IF
          cv%epot_s = 0.5_dp*cv%lambda*(diff_ss)**2.0_dp
          cv%ff_s   = cv%lambda*(diff_ss)
          
          meta_env%epot_s = meta_env%epot_s + cv%epot_s
          meta_env%ekin_s = meta_env%ekin_s + 0.5_dp*cv%mass*cv%vvp**2
       ENDDO
    END IF

    ! write COLVAR file
    iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
         "PRINT%COLVAR",extension=".metadynLog",error=error)
    IF (iw>0) THEN
       IF (meta_env%extended_lagrange) THEN
          WRITE(iw,'(f16.8,70f15.8)')meta_env%time*femtoseconds, &
               (meta_env%metavar(i)%ss0,i=1,meta_env%n_colvar), &
               (meta_env%metavar(i)%ss,i=1,meta_env%n_colvar), &
               (meta_env%metavar(i)%ff_s,i=1,meta_env%n_colvar), &
               (meta_env%metavar(i)%ff_hills,i=1,meta_env%n_colvar), &
               (meta_env%metavar(i)%ff_walls,i=1,meta_env%n_colvar), &
               (meta_env%metavar(i)%vvp,i=1,meta_env%n_colvar), &
                meta_env%epot_s, &
                meta_env%hills_env%energy, &
                meta_env%epot_walls, &
               (meta_env%ekin_s)*2.0_dp/(REAL(meta_env%n_colvar,KIND=dp))*kelvin
       ELSE
          WRITE(iw,'(f16.8,40f13.5)') meta_env%time*femtoseconds,&
               (meta_env%metavar(i)%ss0,i=1,meta_env%n_colvar),&
               (meta_env%metavar(i)%ff_hills,i=1,meta_env%n_colvar),&
               (meta_env%metavar(i)%ff_walls,i=1,meta_env%n_colvar),&
                meta_env%hills_env%energy,&
                meta_env%epot_walls
       END IF
    END IF
    CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
         "PRINT%COLVAR", error=error)
    
    ! Temperature for COLVAR
    IF (meta_env%extended_lagrange) THEN
       temp = meta_env%ekin_s*2.0_dp/(REAL(meta_env%n_colvar,KIND=dp))*kelvin
       meta_env%avg_temp = (meta_env%avg_temp*REAL(meta_env%n_steps,KIND=dp)+&
            temp)/REAL(meta_env%n_steps+1,KIND=dp)
       iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
            "PRINT%TEMPERATURE_COLVAR",extension=".metadynLog",error=error)
       IF (iw > 0) THEN
          WRITE (iw, '(T2,79("-"))')
          WRITE (iw,'( A,T51,f10.2,T71,f10.2)' )' COLVARS INSTANTANEOUS/AVERAGE TEMPERATURE ',&
               temp, meta_env%avg_temp
          WRITE (iw, '(T2,79("-"))')
       ENDIF
       CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
            "PRINT%TEMPERATURE_COLVAR", error=error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE metadyn_write_colvar

! *****************************************************************************
!> \brief Major driver for adding hills and computing forces due to the history
!>        dependent term
!> \par History
!>      04.2004 created
!>      10.2008 Teodoro Laino [tlaino] - University of Zurich
!>              Major rewriting and addition of multiple walkers
! *****************************************************************************
  SUBROUTINE hills(meta_env,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'hills', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, i_hills, ih, &
                                                intermeta_steps, iter_nr, iw, &
                                                n_colvar, n_hills_start, &
                                                n_step, stat
    LOGICAL                                  :: failure, force_gauss
    REAL(KIND=dp)                            :: dfunc, diff_ss, dp2, &
                                                V_now_here, V_to_fes, wtww, ww
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ddp, local_last_hills
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(hills_env_type), POINTER            :: hills_env
    TYPE(metavar_type), DIMENSION(:), &
      POINTER                                :: colvars
    TYPE(multiple_walkers_type), POINTER     :: multiple_walkers

    CALL timeset(routineN,handle)

    failure   = .FALSE.
    NULLIFY(hills_env, multiple_walkers, logger, colvars, ddp, local_last_hills)
    hills_env => meta_env%hills_env
    logger    => cp_error_get_logger(error)
    colvars   => meta_env%metavar
    n_colvar  =  meta_env%n_colvar
    n_step    =  meta_env%n_steps

    ! Create a temporary logger level specific for metadynamics
    CALL cp_add_iter_level(logger%iter_info,"METADYNAMICS",error=error)
    CALL get_meta_iter_level(meta_env, iter_nr, error)
    CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=iter_nr,error=error)   
 
    ! Set-up restart if any
    IF (meta_env%hills_env%restart) THEN
       meta_env%hills_env%restart = .FALSE.
       IF (meta_env%well_tempered) THEN
          CALL restart_hills(hills_env%ss_history, hills_env%delta_s_history, hills_env%ww_history,&
               hills_env%ww, hills_env%n_hills, n_colvar, colvars, meta_env%metadyn_section, error,&
               invdt_history=hills_env%invdt_history)
       ELSE
          CALL restart_hills(hills_env%ss_history, hills_env%delta_s_history, hills_env%ww_history,&
               hills_env%ww, hills_env%n_hills, n_colvar, colvars, meta_env%metadyn_section, error)
       END IF
    END IF

 
    ! Proceed with normal calculation
    intermeta_steps=n_step-hills_env%old_hill_step
    force_gauss=.FALSE.
    IF ((hills_env%min_disp > 0.0_dp).AND.(hills_env%old_hill_number> 0).AND.&
        (MODULO(intermeta_steps,hills_env%min_nt_hills)==0)) THEN
       ALLOCATE(ddp(meta_env%n_colvar), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(local_last_hills(meta_env%n_colvar), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)

       local_last_hills(1:n_colvar)=hills_env%ss_history(1:n_colvar,hills_env%old_hill_number)
     
       !RG Calculate the displacement 
       dp2=0.0_dp
       DO i = 1, n_colvar
          ddp(i) = colvars(i)%ss0 - local_last_hills(i)
          IF (colvars(i)%periodic) THEN
             ! The difference of a periodic COLVAR is always within [-pi,pi]
             ddp(i) = SIGN(1.0_dp,ASIN(SIN(ddp(i))))*ACOS(COS(ddp(i)))
          END IF
          dp2 = dp2 + ddp(i)**2
       ENDDO
       dp2=SQRT(dp2)
 
      IF(dp2 > hills_env%min_disp) force_gauss=.TRUE.
      DEALLOCATE(ddp, stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
      DEALLOCATE(local_last_hills, stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    END IF

    !RG keep into account adaptive hills
    IF (((MODULO(intermeta_steps,hills_env%nt_hills)==0).OR.force_gauss)&
        .AND.(.NOT.meta_env%restart).AND.(hills_env%nt_hills>0)) THEN
       IF (meta_env%do_multiple_walkers) multiple_walkers => meta_env%multiple_walkers

       n_hills_start = hills_env%n_hills
       ! Add the hill corresponding to this location
       IF (meta_env%well_tempered) THEN
          ! Well-Tempered scaling of hills height
          V_now_here = 0._dp
          DO ih=1,hills_env%n_hills
             dp2 = 0._dp
             DO i = 1, n_colvar
               diff_ss = colvars(i)%ss0-hills_env%ss_history(i,ih)
               IF (colvars(i)%periodic) THEN
                  ! The difference of a periodic COLVAR is always within [-pi,pi]
                  diff_ss = SIGN(1.0_dp,ASIN(SIN(diff_ss)))*ACOS(COS(diff_ss))
               END IF
               diff_ss = (diff_ss)/hills_env%delta_s_history(i,ih)
               dp2 = dp2 + diff_ss**2
             ENDDO
             V_to_fes   = 1.0_dp + meta_env%wttemperature*hills_env%invdt_history(ih)
             V_now_here = V_now_here+hills_env%ww_history(ih)/V_to_fes*EXP(-0.5_dp*dp2)
          ENDDO
          wtww=hills_env%ww*EXP(-V_now_here*meta_env%invdt)
          ww=wtww*(1.0_dp + meta_env%wttemperature*meta_env%invdt)
          CALL add_hill_single(hills_env, colvars,           ww, hills_env%n_hills, n_colvar, error, meta_env%invdt)
       ELSE
          CALL add_hill_single(hills_env, colvars, hills_env%ww, hills_env%n_hills, n_colvar, error)
       END IF
       ! Update local n_hills counter
       IF (meta_env%do_multiple_walkers) multiple_walkers%n_hills_local = multiple_walkers%n_hills_local+1

       hills_env%old_hill_number=hills_env%n_hills
       hills_env%old_hill_step=n_step

       ! Update iteration level for printing
       CALL get_meta_iter_level(meta_env, iter_nr, error)
       CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=iter_nr,error=error)

       ! Print just program_run_info
       iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
            "PRINT%PROGRAM_RUN_INFO",extension=".metadynLog",error=error)
       IF (iw>0) THEN
          IF (meta_env%do_multiple_walkers) THEN
             WRITE(iw,'(/,1X,"METADYN|",A,I0,A,I0,A,/)')&
                  ' Global/Local Hills number (',hills_env%n_hills,'/',multiple_walkers%n_hills_local,&
                  ') added.'
          ELSE
             WRITE(iw,'(/,1X,"METADYN|",A,I0,A,/)')' Hills number  (',hills_env%n_hills,') added.'
          END IF
       END IF
       CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
            "PRINT%PROGRAM_RUN_INFO", error=error)

       ! Handle Multiple Walkers
       IF (meta_env%do_multiple_walkers) THEN
          ! Print Local Hills file if requested
          iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
               "PRINT%HILLS",middle_name="LOCAL",extension=".metadynLog",error=error)
          IF (iw>0) THEN
             WRITE(iw,'(f12.1,30f13.5)')meta_env%time*femtoseconds,&
                  (hills_env%ss_history(ih,hills_env%n_hills),ih=1,n_colvar),&
                  (hills_env%delta_s_history(ih,hills_env%n_hills),ih=1,n_colvar),&
                   hills_env%ww_history(hills_env%n_hills)
          END IF
          CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
               "PRINT%HILLS", error=error)

          ! Check the communication buffer of the other walkers
          CALL synchronize_multiple_walkers(multiple_walkers, hills_env, colvars,&
               n_colvar, meta_env%metadyn_section, error)
       END IF

       ! Print Hills file if requested (for multiple walkers this file includes
       ! the Hills coming from all the walkers).
       iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
            "PRINT%HILLS",extension=".metadynLog",error=error)
       IF (iw>0) THEN
          DO i_hills = n_hills_start+1, hills_env%n_hills
             WRITE(iw,'(f12.1,30f13.5)')meta_env%time*femtoseconds,&
                  (hills_env%ss_history(ih,i_hills),ih=1,n_colvar),&
                  (hills_env%delta_s_history(ih,i_hills),ih=1,n_colvar),&
                   hills_env%ww_history(i_hills)
          END DO
       END IF
       CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
            "PRINT%HILLS", error=error)
    END IF

    ! Computes forces due to the hills: history dependent term
    ALLOCATE(ddp(meta_env%n_colvar), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    hills_env%energy = 0.0_dp
    DO ih=1,n_colvar
       colvars(ih)%ff_hills = 0.0_dp
    ENDDO
    DO ih=1,hills_env%n_hills
       dp2 = 0._dp
       DO i = 1, n_colvar
          diff_ss = colvars(i)%ss0-hills_env%ss_history(i,ih)
          IF (colvars(i)%periodic) THEN
             ! The difference of a periodic COLVAR is always within [-pi,pi]
             diff_ss = SIGN(1.0_dp,ASIN(SIN(diff_ss)))*ACOS(COS(diff_ss))
          END IF
          IF (hills_env%delta_s_history(i,ih) == 0.0_dp) THEN
            ! trick: scale = 0 is interpreted as infinitely wide Gaussian hill
            ! instead of infinitely narrow. This way one can combine several
            ! one-dimensional bias potentials in a multi-dimensional metadyn
            ! simulation.
            ddp(i) = 0.0_dp
          ELSE
            ddp(i) = (diff_ss)/hills_env%delta_s_history(i,ih)
          END IF
          dp2 = dp2 + ddp(i)**2
       ENDDO
       ! ff_hills contains the "force" due to the hills
       dfunc = hills_env%ww_history(ih)*EXP(-0.5_dp*dp2)
       IF (meta_env%well_tempered) dfunc=dfunc/(1.0_dp + meta_env%wttemperature*hills_env%invdt_history(ih))
       hills_env%energy = hills_env%energy + dfunc
       DO i=1,n_colvar
          IF (hills_env%delta_s_history(i,ih) /= 0.0_dp) THEN
            ! only apply a force when the Gaussian hill has a finite width in
            ! this direction
            colvars(i)%ff_hills = colvars(i)%ff_hills + ddp(i)/hills_env%delta_s_history(i,ih)*dfunc
          END IF
       ENDDO
    ENDDO
    DEALLOCATE(ddp, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CALL cp_rm_iter_level(logger%iter_info,"METADYNAMICS",error=error)    

    CALL timestop(handle)

  END SUBROUTINE hills

END MODULE metadynamics
