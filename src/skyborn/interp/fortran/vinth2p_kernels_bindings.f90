! =============================================================================
! File    : vinth2p_kernels_bindings.f90
! Purpose : C interoperability bindings split from vinth2p_kernels.f90.
! Notes   : Keep bind(C) entry points separate from the scientific kernels.
! =============================================================================

subroutine dvinth2p_nodes_pa_c(dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
                               nlevi, ncol, nlevo) bind(C, name="dvinth2p_nodes_pa_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dati, dato, hbcofa, hbcofb, plevo, psfc
    real(c_double), value :: p0, spvl
    integer(c_int), value :: intyp, kxtrp, nlevi, ncol, nlevo
    real(real64), pointer :: dati_f(:, :), dato_f(:, :), hbcofa_f(:), hbcofb_f(:), plevo_f(:), psfc_f(:)

    call c_f_pointer(dati, dati_f, [nlevi, ncol])
    call c_f_pointer(dato, dato_f, [nlevo, ncol])
    call c_f_pointer(hbcofa, hbcofa_f, [nlevi])
    call c_f_pointer(hbcofb, hbcofb_f, [nlevi])
    call c_f_pointer(plevo, plevo_f, [nlevo])
    call c_f_pointer(psfc, psfc_f, [ncol])

    call dvinth2p_nodes_pa(dati_f, dato_f, hbcofa_f, hbcofb_f, p0, plevo_f, intyp, psfc_f, spvl, kxtrp, &
                           nlevi, ncol, nlevo)
end subroutine dvinth2p_nodes_pa_c


subroutine dvinth2p_nodes_pa_into_c(dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
                                    nlevi, ncol, nlevo) bind(C, name="dvinth2p_nodes_pa_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dati, dato, hbcofa, hbcofb, plevo, psfc
    real(c_double), value :: p0, spvl
    integer(c_int), value :: intyp, kxtrp, nlevi, ncol, nlevo
    real(real64), pointer :: dati_f(:, :), dato_f(:, :), hbcofa_f(:), hbcofb_f(:), plevo_f(:), psfc_f(:)

    call c_f_pointer(dati, dati_f, [nlevi, ncol])
    call c_f_pointer(dato, dato_f, [nlevo, ncol])
    call c_f_pointer(hbcofa, hbcofa_f, [nlevi])
    call c_f_pointer(hbcofb, hbcofb_f, [nlevi])
    call c_f_pointer(plevo, plevo_f, [nlevo])
    call c_f_pointer(psfc, psfc_f, [ncol])

    call dvinth2p_nodes_pa_into(dati_f, dato_f, hbcofa_f, hbcofb_f, p0, plevo_f, intyp, psfc_f, spvl, kxtrp, &
                                nlevi, ncol, nlevo)
end subroutine dvinth2p_nodes_pa_into_c


subroutine dvinth2p_ecmwf_nodes_pa_c(dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
                                     nlevi, ncol, nlevo, varflg, tbot, phis) &
        bind(C, name="dvinth2p_ecmwf_nodes_pa_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dati, dato, hbcofa, hbcofb, plevo, psfc, tbot, phis
    real(c_double), value :: p0, spvl
    integer(c_int), value :: intyp, kxtrp, nlevi, ncol, nlevo, varflg
    real(real64), pointer :: dati_f(:, :), dato_f(:, :), hbcofa_f(:), hbcofb_f(:), plevo_f(:), psfc_f(:)
    real(real64), pointer :: tbot_f(:), phis_f(:)

    call c_f_pointer(dati, dati_f, [nlevi, ncol])
    call c_f_pointer(dato, dato_f, [nlevo, ncol])
    call c_f_pointer(hbcofa, hbcofa_f, [nlevi])
    call c_f_pointer(hbcofb, hbcofb_f, [nlevi])
    call c_f_pointer(plevo, plevo_f, [nlevo])
    call c_f_pointer(psfc, psfc_f, [ncol])
    call c_f_pointer(tbot, tbot_f, [ncol])
    call c_f_pointer(phis, phis_f, [ncol])

    call dvinth2p_ecmwf_nodes_pa(dati_f, dato_f, hbcofa_f, hbcofb_f, p0, plevo_f, intyp, psfc_f, spvl, kxtrp, &
                                 nlevi, ncol, nlevo, varflg, tbot_f, phis_f)
end subroutine dvinth2p_ecmwf_nodes_pa_c


subroutine dvinth2p_ecmwf_nodes_pa_into_c(dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
                                          nlevi, ncol, nlevo, varflg, tbot, phis) &
        bind(C, name="dvinth2p_ecmwf_nodes_pa_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dati, dato, hbcofa, hbcofb, plevo, psfc, tbot, phis
    real(c_double), value :: p0, spvl
    integer(c_int), value :: intyp, kxtrp, nlevi, ncol, nlevo, varflg
    real(real64), pointer :: dati_f(:, :), dato_f(:, :), hbcofa_f(:), hbcofb_f(:), plevo_f(:), psfc_f(:)
    real(real64), pointer :: tbot_f(:), phis_f(:)

    call c_f_pointer(dati, dati_f, [nlevi, ncol])
    call c_f_pointer(dato, dato_f, [nlevo, ncol])
    call c_f_pointer(hbcofa, hbcofa_f, [nlevi])
    call c_f_pointer(hbcofb, hbcofb_f, [nlevi])
    call c_f_pointer(plevo, plevo_f, [nlevo])
    call c_f_pointer(psfc, psfc_f, [ncol])
    call c_f_pointer(tbot, tbot_f, [ncol])
    call c_f_pointer(phis, phis_f, [ncol])

    call dvinth2p_ecmwf_nodes_pa_into(dati_f, dato_f, hbcofa_f, hbcofb_f, p0, plevo_f, intyp, psfc_f, spvl, &
                                      kxtrp, nlevi, ncol, nlevo, varflg, tbot_f, phis_f)
end subroutine dvinth2p_ecmwf_nodes_pa_into_c


subroutine ddelta_pressure_hybrid_pa_c(psfc, dph, hbcofa, hbcofb, p0, ncol, nlev, nlevo) &
        bind(C, name="ddelta_pressure_hybrid_pa_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: psfc, dph, hbcofa, hbcofb
    real(c_double), value :: p0
    integer(c_int), value :: ncol, nlev, nlevo
    real(real64), pointer :: psfc_f(:), dph_f(:, :), hbcofa_f(:), hbcofb_f(:)

    call c_f_pointer(psfc, psfc_f, [ncol])
    call c_f_pointer(dph, dph_f, [nlevo, ncol])
    call c_f_pointer(hbcofa, hbcofa_f, [nlev])
    call c_f_pointer(hbcofb, hbcofb_f, [nlev])

    call ddelta_pressure_hybrid_pa(psfc_f, dph_f, hbcofa_f, hbcofb_f, p0, ncol, nlev, nlevo)
end subroutine ddelta_pressure_hybrid_pa_c


subroutine ddelta_pressure_hybrid_pa_into_c(psfc, dph, hbcofa, hbcofb, p0, ncol, nlev, nlevo) &
        bind(C, name="ddelta_pressure_hybrid_pa_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: psfc, dph, hbcofa, hbcofb
    real(c_double), value :: p0
    integer(c_int), value :: ncol, nlev, nlevo
    real(real64), pointer :: psfc_f(:), dph_f(:, :), hbcofa_f(:), hbcofb_f(:)

    call c_f_pointer(psfc, psfc_f, [ncol])
    call c_f_pointer(dph, dph_f, [nlevo, ncol])
    call c_f_pointer(hbcofa, hbcofa_f, [nlev])
    call c_f_pointer(hbcofb, hbcofb_f, [nlev])

    call ddelta_pressure_hybrid_pa_into(psfc_f, dph_f, hbcofa_f, hbcofb_f, p0, ncol, nlev, nlevo)
end subroutine ddelta_pressure_hybrid_pa_into_c


subroutine dpressure_at_hybrid_levels_pa_c(psfc, pressure, hbcofa, hbcofb, p0, ncol, nlev) &
        bind(C, name="dpressure_at_hybrid_levels_pa_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: psfc, pressure, hbcofa, hbcofb
    real(c_double), value :: p0
    integer(c_int), value :: ncol, nlev
    real(real64), pointer :: psfc_f(:), pressure_f(:, :), hbcofa_f(:), hbcofb_f(:)

    call c_f_pointer(psfc, psfc_f, [ncol])
    call c_f_pointer(pressure, pressure_f, [nlev, ncol])
    call c_f_pointer(hbcofa, hbcofa_f, [nlev])
    call c_f_pointer(hbcofb, hbcofb_f, [nlev])

    call dpressure_at_hybrid_levels_pa(psfc_f, pressure_f, hbcofa_f, hbcofb_f, p0, ncol, nlev)
end subroutine dpressure_at_hybrid_levels_pa_c


subroutine dpressure_at_hybrid_levels_pa_into_c(psfc, pressure, hbcofa, hbcofb, p0, ncol, nlev) &
        bind(C, name="dpressure_at_hybrid_levels_pa_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: psfc, pressure, hbcofa, hbcofb
    real(c_double), value :: p0
    integer(c_int), value :: ncol, nlev
    real(real64), pointer :: psfc_f(:), pressure_f(:, :), hbcofa_f(:), hbcofb_f(:)

    call c_f_pointer(psfc, psfc_f, [ncol])
    call c_f_pointer(pressure, pressure_f, [nlev, ncol])
    call c_f_pointer(hbcofa, hbcofa_f, [nlev])
    call c_f_pointer(hbcofb, hbcofb_f, [nlev])

    call dpressure_at_hybrid_levels_pa_into(psfc_f, pressure_f, hbcofa_f, hbcofb_f, p0, ncol, nlev)
end subroutine dpressure_at_hybrid_levels_pa_into_c


subroutine dgeopotential_height_hybrid_corder_pa_into_c(temp_flat, q_flat, z3_flat, psfc, phis, hyai, hybi, p0, &
                                                        nouter, nlev, ninner) &
        bind(C, name="dgeopotential_height_hybrid_corder_pa_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: temp_flat, q_flat, z3_flat, psfc, phis, hyai, hybi
    real(c_double), value :: p0
    integer(c_int), value :: nouter, nlev, ninner
    real(real64), pointer :: temp_flat_f(:), q_flat_f(:), z3_flat_f(:), psfc_f(:), phis_f(:), hyai_f(:), hybi_f(:)

    call c_f_pointer(temp_flat, temp_flat_f, [nouter * nlev * ninner])
    call c_f_pointer(q_flat, q_flat_f, [nouter * nlev * ninner])
    call c_f_pointer(z3_flat, z3_flat_f, [nouter * nlev * ninner])
    call c_f_pointer(psfc, psfc_f, [nouter * ninner])
    call c_f_pointer(phis, phis_f, [nouter * ninner])
    call c_f_pointer(hyai, hyai_f, [nlev + 1])
    call c_f_pointer(hybi, hybi_f, [nlev + 1])

    call dgeopotential_height_hybrid_corder_pa_into(temp_flat_f, q_flat_f, z3_flat_f, psfc_f, phis_f, hyai_f, hybi_f, &
                                                    p0, nouter, nlev, ninner)
end subroutine dgeopotential_height_hybrid_corder_pa_into_c


subroutine dsigma2hybrid_nodes_c(dsigmai, dato, sigmai, sigmao, intyp, spvl, nlevi, ncol, nlevo) &
        bind(C, name="dsigma2hybrid_nodes_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dsigmai, dato, sigmai, sigmao
    real(c_double), value :: spvl
    integer(c_int), value :: intyp, nlevi, ncol, nlevo
    real(real64), pointer :: dati_f(:, :), dato_f(:, :), sigmai_f(:), sigmao_f(:, :)

    call c_f_pointer(dsigmai, dati_f, [nlevi, ncol])
    call c_f_pointer(dato, dato_f, [nlevo, ncol])
    call c_f_pointer(sigmai, sigmai_f, [nlevi])
    call c_f_pointer(sigmao, sigmao_f, [nlevo, ncol])

    call dsigma2hybrid_nodes(dati_f, dato_f, sigmai_f, sigmao_f, intyp, spvl, nlevi, ncol, nlevo)
end subroutine dsigma2hybrid_nodes_c


subroutine dsigma2hybrid_nodes_into_c(dsigmai, dato, sigmai, sigmao, intyp, spvl, nlevi, ncol, nlevo) &
        bind(C, name="dsigma2hybrid_nodes_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dsigmai, dato, sigmai, sigmao
    real(c_double), value :: spvl
    integer(c_int), value :: intyp, nlevi, ncol, nlevo
    real(real64), pointer :: dati_f(:, :), dato_f(:, :), sigmai_f(:), sigmao_f(:, :)

    call c_f_pointer(dsigmai, dati_f, [nlevi, ncol])
    call c_f_pointer(dato, dato_f, [nlevo, ncol])
    call c_f_pointer(sigmai, sigmai_f, [nlevi])
    call c_f_pointer(sigmao, sigmao_f, [nlevo, ncol])

    call dsigma2hybrid_nodes_into(dati_f, dato_f, sigmai_f, sigmao_f, intyp, spvl, nlevi, ncol, nlevo)
end subroutine dsigma2hybrid_nodes_into_c


subroutine dsigma2hybrid_nodes_corder_into_c(dati_flat, dato_flat, sigmai, hyam, hybm, p0, psfc, intyp, spvl, &
                                             nouter, nlevi, ninner, nlevo) &
        bind(C, name="dsigma2hybrid_nodes_corder_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dati_flat, dato_flat, sigmai, hyam, hybm, psfc
    real(c_double), value :: p0, spvl
    integer(c_int), value :: intyp, nouter, nlevi, ninner, nlevo
    real(real64), pointer :: dati_flat_f(:), dato_flat_f(:), sigmai_f(:), hyam_f(:), hybm_f(:), psfc_f(:)

    call c_f_pointer(dati_flat, dati_flat_f, [nouter * nlevi * ninner])
    call c_f_pointer(dato_flat, dato_flat_f, [nouter * nlevo * ninner])
    call c_f_pointer(sigmai, sigmai_f, [nlevi])
    call c_f_pointer(hyam, hyam_f, [nlevo])
    call c_f_pointer(hybm, hybm_f, [nlevo])
    call c_f_pointer(psfc, psfc_f, [nouter * ninner])

    call dsigma2hybrid_nodes_corder_into(dati_flat_f, dato_flat_f, sigmai_f, hyam_f, hybm_f, p0, psfc_f, intyp, spvl, &
                                         nouter, nlevi, ninner, nlevo)
end subroutine dsigma2hybrid_nodes_corder_into_c


subroutine dvinth2p_nodes_corder_pa_into_c(dati_flat, dato_flat, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, &
                                           kxtrp, nouter, nlevi, ninner, nlevo) &
        bind(C, name="dvinth2p_nodes_corder_pa_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dati_flat, dato_flat, hbcofa, hbcofb, plevo, psfc
    real(c_double), value :: p0, spvl
    integer(c_int), value :: intyp, kxtrp, nouter, nlevi, ninner, nlevo
    real(real64), pointer :: dati_flat_f(:), dato_flat_f(:), hbcofa_f(:), hbcofb_f(:), plevo_f(:), psfc_f(:)

    call c_f_pointer(dati_flat, dati_flat_f, [nouter * nlevi * ninner])
    call c_f_pointer(dato_flat, dato_flat_f, [nouter * nlevo * ninner])
    call c_f_pointer(hbcofa, hbcofa_f, [nlevi])
    call c_f_pointer(hbcofb, hbcofb_f, [nlevi])
    call c_f_pointer(plevo, plevo_f, [nlevo])
    call c_f_pointer(psfc, psfc_f, [nouter * ninner])

    call dvinth2p_nodes_corder_pa_into(dati_flat_f, dato_flat_f, hbcofa_f, hbcofb_f, p0, plevo_f, intyp, psfc_f, spvl, &
                                       kxtrp, nouter, nlevi, ninner, nlevo)
end subroutine dvinth2p_nodes_corder_pa_into_c


subroutine dvinth2p_ecmwf_nodes_corder_pa_into_c(dati_flat, dato_flat, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, &
                                                 kxtrp, nouter, nlevi, ninner, nlevo, varflg, tbot, phis) &
        bind(C, name="dvinth2p_ecmwf_nodes_corder_pa_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dati_flat, dato_flat, hbcofa, hbcofb, plevo, psfc, tbot, phis
    real(c_double), value :: p0, spvl
    integer(c_int), value :: intyp, kxtrp, nouter, nlevi, ninner, nlevo, varflg
    real(real64), pointer :: dati_flat_f(:), dato_flat_f(:), hbcofa_f(:), hbcofb_f(:), plevo_f(:), psfc_f(:)
    real(real64), pointer :: tbot_f(:), phis_f(:)

    call c_f_pointer(dati_flat, dati_flat_f, [nouter * nlevi * ninner])
    call c_f_pointer(dato_flat, dato_flat_f, [nouter * nlevo * ninner])
    call c_f_pointer(hbcofa, hbcofa_f, [nlevi])
    call c_f_pointer(hbcofb, hbcofb_f, [nlevi])
    call c_f_pointer(plevo, plevo_f, [nlevo])
    call c_f_pointer(psfc, psfc_f, [nouter * ninner])
    call c_f_pointer(tbot, tbot_f, [nouter * ninner])
    call c_f_pointer(phis, phis_f, [nouter * ninner])

    call dvinth2p_ecmwf_nodes_corder_pa_into(dati_flat_f, dato_flat_f, hbcofa_f, hbcofb_f, p0, plevo_f, intyp, psfc_f, &
                                             spvl, kxtrp, nouter, nlevi, ninner, nlevo, varflg, tbot_f, phis_f)
end subroutine dvinth2p_ecmwf_nodes_corder_pa_into_c
