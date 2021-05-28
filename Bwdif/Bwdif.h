#pragma once

#include <type_traits>

#ifdef BWDIF_X86
#include "VCL2/vectorclass.h"
#endif

/*
 * Filter coefficients coef_lf and coef_hf taken from BBC PH-2071 (Weston 3 Field Deinterlacer).
 * Used when there is spatial and temporal interpolation.
 * Filter coefficients coef_sp are used when there is spatial interpolation only.
 * Adjusted for matching visual sharpness impression of spatial and temporal interpolation.
 */
static constexpr uint16_t coef_lf[2] = { 4309, 213 };
static constexpr uint16_t coef_hf[3] = { 5570, 3801, 1016 };
static constexpr uint16_t coef_sp[2] = { 5077, 981 };
static constexpr float coef_lf_f[2] = { 4309 / 8192.0f, 213 / 8192.0f };
static constexpr float coef_hf_f[3] = { 5570 / 8192.0f, 3801 / 8192.0f, 1016 / 8192.0f };
static constexpr float coef_sp_f[2] = { 5077 / 8192.0f, 981 / 8192.0f };

template<typename T>
static constexpr auto coefLF() noexcept {
    if constexpr (std::is_integral_v<T>)
        return coef_lf;
    else
        return coef_lf_f;
}

template<typename T>
static constexpr auto coefHF() noexcept {
    if constexpr (std::is_integral_v<T>)
        return coef_hf;
    else
        return coef_hf_f;
}

template<typename T>
static constexpr auto coefSP() noexcept {
    if constexpr (std::is_integral_v<T>)
        return coef_sp;
    else
        return coef_sp_f;
}

template<typename T, typename V>
static constexpr auto div2(const V x) noexcept {
    if constexpr (std::is_integral_v<T>)
        return x >> 1;
    else
        return x * (1 / 2.0f);
}

template<typename T, typename V>
static constexpr auto div4(const V x) noexcept {
    if constexpr (std::is_integral_v<T>)
        return x >> 2;
    else
        return x * (1 / 4.0f);
}
