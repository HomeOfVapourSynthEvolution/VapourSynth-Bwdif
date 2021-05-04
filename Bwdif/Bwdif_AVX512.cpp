#ifdef BWDIF_X86
#include "Bwdif.h"

template<typename pixel_t, bool spat, bool hasEdeint>
void filterEdge_avx512(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst,
                       const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept {
    const pixel_t* prev2 = reinterpret_cast<const pixel_t*>(_prev2);
    const pixel_t* prev = reinterpret_cast<const pixel_t*>(_prev);
    const pixel_t* cur = reinterpret_cast<const pixel_t*>(_cur);
    const pixel_t* next = reinterpret_cast<const pixel_t*>(_next);
    const pixel_t* next2 = reinterpret_cast<const pixel_t*>(_next2);
    const pixel_t* edeint = reinterpret_cast<const pixel_t*>(_edeint);
    pixel_t* dst = reinterpret_cast<pixel_t*>(_dst);

    const pixel_t* prev2Above2 = prev2 - stride2;
    const pixel_t* prev2Below2 = prev2 + stride2;
    const pixel_t* prevAbove = prev + negativeStride;
    const pixel_t* prevBelow = prev + positiveStride;
    const pixel_t* curAbove = cur + negativeStride;
    const pixel_t* curBelow = cur + positiveStride;
    const pixel_t* nextAbove = next + negativeStride;
    const pixel_t* nextBelow = next + positiveStride;
    const pixel_t* next2Above2 = next2 - stride2;
    const pixel_t* next2Below2 = next2 + stride2;

    for (int x = 0; x < width; x += step) {
        if constexpr (std::is_same_v<pixel_t, uint8_t>) {
            const Vec32s c = Vec32s().load_32uc(curAbove + x);
            const Vec32s d = (Vec32s().load_32uc(prev2 + x) + Vec32s().load_32uc(next2 + x)) >> 1;
            const Vec32s e = Vec32s().load_32uc(curBelow + x);
            const Vec32s temporal_diff0 = abs(Vec32s().load_32uc(prev2 + x) - Vec32s().load_32uc(next2 + x));
            const Vec32s temporal_diff1 = (abs(Vec32s().load_32uc(prevAbove + x) - c) + abs(Vec32s().load_32uc(prevBelow + x) - e)) >> 1;
            const Vec32s temporal_diff2 = (abs(Vec32s().load_32uc(nextAbove + x) - c) + abs(Vec32s().load_32uc(nextBelow + x) - e)) >> 1;
            const Vec32s temporal_diff = max(max(temporal_diff0 >> 1, temporal_diff1), temporal_diff2);

            Vec32s diff = temporal_diff;
            if constexpr (spat) {
                const Vec32s b = ((Vec32s().load_32uc(prev2Above2 + x) + Vec32s().load_32uc(next2Above2 + x)) >> 1) - c;
                const Vec32s f = ((Vec32s().load_32uc(prev2Below2 + x) + Vec32s().load_32uc(next2Below2 + x)) >> 1) - e;
                const Vec32s dc = d - c;
                const Vec32s de = d - e;
                const Vec32s maximum = max(max(de, dc), min(b, f));
                const Vec32s minimum = min(min(de, dc), max(b, f));
                diff = max(max(temporal_diff, minimum), -maximum);
            }

            Vec32s interpol;
            if constexpr (hasEdeint)
                interpol = Vec32s().load_32uc(edeint + x);
            else
                interpol = (c + e) >> 1;
            interpol = min(max(interpol, d - diff), d + diff);
            interpol = select(temporal_diff == 0, d, interpol);

            const auto result = compress_saturated_s2u(interpol, zero_si512()).get_low();
            result.store_nt(dst + x);
        } else if constexpr (std::is_same_v<pixel_t, uint16_t>) {
            const Vec16i c = Vec16i().load_16us(curAbove + x);
            const Vec16i d = (Vec16i().load_16us(prev2 + x) + Vec16i().load_16us(next2 + x)) >> 1;
            const Vec16i e = Vec16i().load_16us(curBelow + x);
            const Vec16i temporal_diff0 = abs(Vec16i().load_16us(prev2 + x) - Vec16i().load_16us(next2 + x));
            const Vec16i temporal_diff1 = (abs(Vec16i().load_16us(prevAbove + x) - c) + abs(Vec16i().load_16us(prevBelow + x) - e)) >> 1;
            const Vec16i temporal_diff2 = (abs(Vec16i().load_16us(nextAbove + x) - c) + abs(Vec16i().load_16us(nextBelow + x) - e)) >> 1;
            const Vec16i temporal_diff = max(max(temporal_diff0 >> 1, temporal_diff1), temporal_diff2);

            Vec16i diff = temporal_diff;
            if constexpr (spat) {
                const Vec16i b = ((Vec16i().load_16us(prev2Above2 + x) + Vec16i().load_16us(next2Above2 + x)) >> 1) - c;
                const Vec16i f = ((Vec16i().load_16us(prev2Below2 + x) + Vec16i().load_16us(next2Below2 + x)) >> 1) - e;
                const Vec16i dc = d - c;
                const Vec16i de = d - e;
                const Vec16i maximum = max(max(de, dc), min(b, f));
                const Vec16i minimum = min(min(de, dc), max(b, f));
                diff = max(max(temporal_diff, minimum), -maximum);
            }

            Vec16i interpol;
            if constexpr (hasEdeint)
                interpol = Vec16i().load_16us(edeint + x);
            else
                interpol = (c + e) >> 1;
            interpol = min(max(interpol, d - diff), d + diff);
            interpol = select(temporal_diff == 0, d, interpol);

            const auto result = compress_saturated_s2u(interpol, zero_si512()).get_low();
            result.store_nt(dst + x);
        } else {
            const Vec16f c = Vec16f().load_a(curAbove + x);
            const Vec16f d = (Vec16f().load_a(prev2 + x) + Vec16f().load_a(next2 + x)) * 0.5f;
            const Vec16f e = Vec16f().load_a(curBelow + x);
            const Vec16f temporal_diff0 = abs(Vec16f().load_a(prev2 + x) - Vec16f().load_a(next2 + x));
            const Vec16f temporal_diff1 = (abs(Vec16f().load_a(prevAbove + x) - c) + abs(Vec16f().load_a(prevBelow + x) - e)) * 0.5f;
            const Vec16f temporal_diff2 = (abs(Vec16f().load_a(nextAbove + x) - c) + abs(Vec16f().load_a(nextBelow + x) - e)) * 0.5f;
            const Vec16f temporal_diff = max(max(temporal_diff0 * 0.5f, temporal_diff1), temporal_diff2);

            Vec16f diff = temporal_diff;
            if constexpr (spat) {
                const Vec16f b = ((Vec16f().load_a(prev2Above2 + x) + Vec16f().load_a(next2Above2 + x)) * 0.5f) - c;
                const Vec16f f = ((Vec16f().load_a(prev2Below2 + x) + Vec16f().load_a(next2Below2 + x)) * 0.5f) - e;
                const Vec16f dc = d - c;
                const Vec16f de = d - e;
                const Vec16f maximum = max(max(de, dc), min(b, f));
                const Vec16f minimum = min(min(de, dc), max(b, f));
                diff = max(max(temporal_diff, minimum), -maximum);
            }

            Vec16f interpol;
            if constexpr (hasEdeint)
                interpol = Vec16f().load_a(edeint + x);
            else
                interpol = (c + e) * 0.5f;
            interpol = min(max(interpol, d - diff), d + diff);
            interpol = select(temporal_diff == 0, d, interpol);

            interpol.store_nt(dst + x);
        }
    }
}

template<typename pixel_t, bool hasEdeint>
void filterLine_avx512(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst,
                       const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step) noexcept {
    const pixel_t* prev2 = reinterpret_cast<const pixel_t*>(_prev2);
    const pixel_t* prev = reinterpret_cast<const pixel_t*>(_prev);
    const pixel_t* cur = reinterpret_cast<const pixel_t*>(_cur);
    const pixel_t* next = reinterpret_cast<const pixel_t*>(_next);
    const pixel_t* next2 = reinterpret_cast<const pixel_t*>(_next2);
    const pixel_t* edeint = reinterpret_cast<const pixel_t*>(_edeint);
    pixel_t* dst = reinterpret_cast<pixel_t*>(_dst);

    const pixel_t* prev2Above4 = prev2 - stride4;
    const pixel_t* prev2Above2 = prev2 - stride2;
    const pixel_t* prev2Below2 = prev2 + stride2;
    const pixel_t* prev2Below4 = prev2 + stride4;
    const pixel_t* prevAbove = prev - stride;
    const pixel_t* prevBelow = prev + stride;
    const pixel_t* curAbove3 = cur - stride3;
    const pixel_t* curAbove = cur - stride;
    const pixel_t* curBelow = cur + stride;
    const pixel_t* curBelow3 = cur + stride3;
    const pixel_t* nextAbove = next - stride;
    const pixel_t* nextBelow = next + stride;
    const pixel_t* next2Above4 = next2 - stride4;
    const pixel_t* next2Above2 = next2 - stride2;
    const pixel_t* next2Below2 = next2 + stride2;
    const pixel_t* next2Below4 = next2 + stride4;

    for (int x = 0; x < width; x += step) {
        if constexpr (std::is_same_v<pixel_t, uint8_t>) {
            const Vec16i c = Vec16i().load_16uc(curAbove + x);
            const Vec16i d = (Vec16i().load_16uc(prev2 + x) + Vec16i().load_16uc(next2 + x)) >> 1;
            const Vec16i e = Vec16i().load_16uc(curBelow + x);
            const Vec16i temporal_diff0 = abs(Vec16i().load_16uc(prev2 + x) - Vec16i().load_16uc(next2 + x));
            const Vec16i temporal_diff1 = (abs(Vec16i().load_16uc(prevAbove + x) - c) + abs(Vec16i().load_16uc(prevBelow + x) - e)) >> 1;
            const Vec16i temporal_diff2 = (abs(Vec16i().load_16uc(nextAbove + x) - c) + abs(Vec16i().load_16uc(nextBelow + x) - e)) >> 1;
            const Vec16i temporal_diff = max(max(temporal_diff0 >> 1, temporal_diff1), temporal_diff2);

            const Vec16i b = ((Vec16i().load_16uc(prev2Above2 + x) + Vec16i().load_16uc(next2Above2 + x)) >> 1) - c;
            const Vec16i f = ((Vec16i().load_16uc(prev2Below2 + x) + Vec16i().load_16uc(next2Below2 + x)) >> 1) - e;
            const Vec16i dc = d - c;
            const Vec16i de = d - e;
            const Vec16i maximum = max(max(de, dc), min(b, f));
            const Vec16i minimum = min(min(de, dc), max(b, f));
            const Vec16i diff = max(max(temporal_diff, minimum), -maximum);

            const Vec16i interpol1 = (((coef_hf[0] * (Vec16i().load_16uc(prev2 + x) + Vec16i().load_16uc(next2 + x))
                                        - coef_hf[1] * (Vec16i().load_16uc(prev2Above2 + x) + Vec16i().load_16uc(next2Above2 + x) + Vec16i().load_16uc(prev2Below2 + x) + Vec16i().load_16uc(next2Below2 + x))
                                        + coef_hf[2] * (Vec16i().load_16uc(prev2Above4 + x) + Vec16i().load_16uc(next2Above4 + x) + Vec16i().load_16uc(prev2Below4 + x) + Vec16i().load_16uc(next2Below4 + x))) >> 2)
                                      + coef_lf[0] * (c + e) - coef_lf[1] * (Vec16i().load_16uc(curAbove3 + x) + Vec16i().load_16uc(curBelow3 + x))) >> 13;

            Vec16i interpol2;
            if constexpr (hasEdeint)
                interpol2 = Vec16i().load_16uc(edeint + x);
            else
                interpol2 = (coef_sp[0] * (c + e) - coef_sp[1] * (Vec16i().load_16uc(curAbove3 + x) + Vec16i().load_16uc(curBelow3 + x))) >> 13;

            Vec16i interpol = select(abs(c - e) > temporal_diff0, interpol1, interpol2);
            interpol = min(max(interpol, d - diff), d + diff);
            interpol = select(temporal_diff == 0, d, interpol);

            const auto result = compress_saturated_s2u(compress_saturated(interpol, zero_si512()), zero_si512()).get_low().get_low();
            result.store_nt(dst + x);
        } else if constexpr (std::is_same_v<pixel_t, uint16_t>) {
            const Vec16i c = Vec16i().load_16us(curAbove + x);
            const Vec16i d = (Vec16i().load_16us(prev2 + x) + Vec16i().load_16us(next2 + x)) >> 1;
            const Vec16i e = Vec16i().load_16us(curBelow + x);
            const Vec16i temporal_diff0 = abs(Vec16i().load_16us(prev2 + x) - Vec16i().load_16us(next2 + x));
            const Vec16i temporal_diff1 = (abs(Vec16i().load_16us(prevAbove + x) - c) + abs(Vec16i().load_16us(prevBelow + x) - e)) >> 1;
            const Vec16i temporal_diff2 = (abs(Vec16i().load_16us(nextAbove + x) - c) + abs(Vec16i().load_16us(nextBelow + x) - e)) >> 1;
            const Vec16i temporal_diff = max(max(temporal_diff0 >> 1, temporal_diff1), temporal_diff2);

            const Vec16i b = ((Vec16i().load_16us(prev2Above2 + x) + Vec16i().load_16us(next2Above2 + x)) >> 1) - c;
            const Vec16i f = ((Vec16i().load_16us(prev2Below2 + x) + Vec16i().load_16us(next2Below2 + x)) >> 1) - e;
            const Vec16i dc = d - c;
            const Vec16i de = d - e;
            const Vec16i maximum = max(max(de, dc), min(b, f));
            const Vec16i minimum = min(min(de, dc), max(b, f));
            const Vec16i diff = max(max(temporal_diff, minimum), -maximum);

            const Vec16i interpol1 = (((coef_hf[0] * (Vec16i().load_16us(prev2 + x) + Vec16i().load_16us(next2 + x))
                                        - coef_hf[1] * (Vec16i().load_16us(prev2Above2 + x) + Vec16i().load_16us(next2Above2 + x) + Vec16i().load_16us(prev2Below2 + x) + Vec16i().load_16us(next2Below2 + x))
                                        + coef_hf[2] * (Vec16i().load_16us(prev2Above4 + x) + Vec16i().load_16us(next2Above4 + x) + Vec16i().load_16us(prev2Below4 + x) + Vec16i().load_16us(next2Below4 + x))) >> 2)
                                      + coef_lf[0] * (c + e) - coef_lf[1] * (Vec16i().load_16us(curAbove3 + x) + Vec16i().load_16us(curBelow3 + x))) >> 13;

            Vec16i interpol2;
            if constexpr (hasEdeint)
                interpol2 = Vec16i().load_16us(edeint + x);
            else
                interpol2 = (coef_sp[0] * (c + e) - coef_sp[1] * (Vec16i().load_16us(curAbove3 + x) + Vec16i().load_16us(curBelow3 + x))) >> 13;

            Vec16i interpol = select(abs(c - e) > temporal_diff0, interpol1, interpol2);
            interpol = min(max(interpol, d - diff), d + diff);
            interpol = select(temporal_diff == 0, d, interpol);

            const auto result = compress_saturated_s2u(interpol, zero_si512()).get_low();
            result.store_nt(dst + x);
        } else {
            const Vec16f c = Vec16f().load_a(curAbove + x);
            const Vec16f d = (Vec16f().load_a(prev2 + x) + Vec16f().load_a(next2 + x)) * 0.5f;
            const Vec16f e = Vec16f().load_a(curBelow + x);
            const Vec16f temporal_diff0 = abs(Vec16f().load_a(prev2 + x) - Vec16f().load_a(next2 + x));
            const Vec16f temporal_diff1 = (abs(Vec16f().load_a(prevAbove + x) - c) + abs(Vec16f().load_a(prevBelow + x) - e)) * 0.5f;
            const Vec16f temporal_diff2 = (abs(Vec16f().load_a(nextAbove + x) - c) + abs(Vec16f().load_a(nextBelow + x) - e)) * 0.5f;
            const Vec16f temporal_diff = max(max(temporal_diff0 * 0.5f, temporal_diff1), temporal_diff2);

            const Vec16f b = ((Vec16f().load_a(prev2Above2 + x) + Vec16f().load_a(next2Above2 + x)) * 0.5f) - c;
            const Vec16f f = ((Vec16f().load_a(prev2Below2 + x) + Vec16f().load_a(next2Below2 + x)) * 0.5f) - e;
            const Vec16f dc = d - c;
            const Vec16f de = d - e;
            const Vec16f maximum = max(max(de, dc), min(b, f));
            const Vec16f minimum = min(min(de, dc), max(b, f));
            const Vec16f diff = max(max(temporal_diff, minimum), -maximum);

            const Vec16f interpol1 = ((coef_hf_f[0] * (Vec16f().load_a(prev2 + x) + Vec16f().load_a(next2 + x))
                                       - coef_hf_f[1] * (Vec16f().load_a(prev2Above2 + x) + Vec16f().load_a(next2Above2 + x) + Vec16f().load_a(prev2Below2 + x) + Vec16f().load_a(next2Below2 + x))
                                       + coef_hf_f[2] * (Vec16f().load_a(prev2Above4 + x) + Vec16f().load_a(next2Above4 + x) + Vec16f().load_a(prev2Below4 + x) + Vec16f().load_a(next2Below4 + x))) * 0.25f
                                      + coef_lf_f[0] * (c + e) - coef_lf_f[1] * (Vec16f().load_a(curAbove3 + x) + Vec16f().load_a(curBelow3 + x)));

            Vec16f interpol2;
            if constexpr (hasEdeint)
                interpol2 = Vec16f().load_a(edeint + x);
            else
                interpol2 = coef_sp_f[0] * (c + e) - coef_sp_f[1] * (Vec16f().load_a(curAbove3 + x) + Vec16f().load_a(curBelow3 + x));

            Vec16f interpol = select(abs(c - e) > temporal_diff0, interpol1, interpol2);
            interpol = min(max(interpol, d - diff), d + diff);
            interpol = select(temporal_diff == 0, d, interpol);

            interpol.store_nt(dst + x);
        }
    }
}

template void filterEdge_avx512<uint8_t, true, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx512<uint8_t, true, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx512<uint8_t, false, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx512<uint8_t, false, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx512<uint16_t, true, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx512<uint16_t, true, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx512<uint16_t, false, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx512<uint16_t, false, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx512<float, true, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx512<float, true, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx512<float, false, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx512<float, false, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;

template void filterLine_avx512<uint8_t, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step) noexcept;
template void filterLine_avx512<uint8_t, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step) noexcept;
template void filterLine_avx512<uint16_t, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step) noexcept;
template void filterLine_avx512<uint16_t, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step) noexcept;
template void filterLine_avx512<float, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step) noexcept;
template void filterLine_avx512<float, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step) noexcept;
#endif
