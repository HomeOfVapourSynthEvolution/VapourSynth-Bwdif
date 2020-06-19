#ifdef BWDIF_X86
#include "Bwdif.h"

template<typename pixel_t, bool spat>
void filterEdge_avx2(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width,
                     const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept {
    const pixel_t * prev2 = reinterpret_cast<const pixel_t *>(_prev2);
    const pixel_t * prev = reinterpret_cast<const pixel_t *>(_prev);
    const pixel_t * cur = reinterpret_cast<const pixel_t *>(_cur);
    const pixel_t * next = reinterpret_cast<const pixel_t *>(_next);
    const pixel_t * next2 = reinterpret_cast<const pixel_t *>(_next2);
    pixel_t * dst = reinterpret_cast<pixel_t *>(_dst);

    const pixel_t * prev2Above2 = prev2 - stride2;
    const pixel_t * prev2Below2 = prev2 + stride2;
    const pixel_t * prevAbove = prev + negativeStride;
    const pixel_t * prevBelow = prev + positiveStride;
    const pixel_t * curAbove = cur + negativeStride;
    const pixel_t * curBelow = cur + positiveStride;
    const pixel_t * nextAbove = next + negativeStride;
    const pixel_t * nextBelow = next + positiveStride;
    const pixel_t * next2Above2 = next2 - stride2;
    const pixel_t * next2Below2 = next2 + stride2;

    for (int x = 0; x < width; x += step) {
        if constexpr (std::is_same_v<pixel_t, uint8_t>) {
            const Vec16s c = Vec16s().load_16uc(curAbove + x);
            const Vec16s d = (Vec16s().load_16uc(prev2 + x) + Vec16s().load_16uc(next2 + x)) >> 1;
            const Vec16s e = Vec16s().load_16uc(curBelow + x);
            const Vec16s temporal_diff0 = abs(Vec16s().load_16uc(prev2 + x) - Vec16s().load_16uc(next2 + x));
            const Vec16s temporal_diff1 = (abs(Vec16s().load_16uc(prevAbove + x) - c) + abs(Vec16s().load_16uc(prevBelow + x) - e)) >> 1;
            const Vec16s temporal_diff2 = (abs(Vec16s().load_16uc(nextAbove + x) - c) + abs(Vec16s().load_16uc(nextBelow + x) - e)) >> 1;
            const Vec16s temporal_diff = max(max(temporal_diff0 >> 1, temporal_diff1), temporal_diff2);

            Vec16s diff;
            if constexpr (spat) {
                const Vec16s b = ((Vec16s().load_16uc(prev2Above2 + x) + Vec16s().load_16uc(next2Above2 + x)) >> 1) - c;
                const Vec16s f = ((Vec16s().load_16uc(prev2Below2 + x) + Vec16s().load_16uc(next2Below2 + x)) >> 1) - e;
                const Vec16s dc = d - c;
                const Vec16s de = d - e;
                const Vec16s maximum = max(max(de, dc), min(b, f));
                const Vec16s minimum = min(min(de, dc), max(b, f));
                diff = max(max(temporal_diff, minimum), -maximum);
            }

            Vec16s interpol = min(max((c + e) >> 1, d - diff), d + diff);
            interpol = select(temporal_diff == 0, d, interpol);
            const auto result = compress_saturated_s2u(interpol, zero_si256()).get_low();
            result.store_nt(dst + x);
        } else if constexpr (std::is_same_v<pixel_t, uint16_t>) {
            const Vec8i c = Vec8i().load_8us(curAbove + x);
            const Vec8i d = (Vec8i().load_8us(prev2 + x) + Vec8i().load_8us(next2 + x)) >> 1;
            const Vec8i e = Vec8i().load_8us(curBelow + x);
            const Vec8i temporal_diff0 = abs(Vec8i().load_8us(prev2 + x) - Vec8i().load_8us(next2 + x));
            const Vec8i temporal_diff1 = (abs(Vec8i().load_8us(prevAbove + x) - c) + abs(Vec8i().load_8us(prevBelow + x) - e)) >> 1;
            const Vec8i temporal_diff2 = (abs(Vec8i().load_8us(nextAbove + x) - c) + abs(Vec8i().load_8us(nextBelow + x) - e)) >> 1;
            const Vec8i temporal_diff = max(max(temporal_diff0 >> 1, temporal_diff1), temporal_diff2);

            Vec8i diff;
            if constexpr (spat) {
                const Vec8i b = ((Vec8i().load_8us(prev2Above2 + x) + Vec8i().load_8us(next2Above2 + x)) >> 1) - c;
                const Vec8i f = ((Vec8i().load_8us(prev2Below2 + x) + Vec8i().load_8us(next2Below2 + x)) >> 1) - e;
                const Vec8i dc = d - c;
                const Vec8i de = d - e;
                const Vec8i maximum = max(max(de, dc), min(b, f));
                const Vec8i minimum = min(min(de, dc), max(b, f));
                diff = max(max(temporal_diff, minimum), -maximum);
            }

            Vec8i interpol = min(max((c + e) >> 1, d - diff), d + diff);
            interpol = select(temporal_diff == 0, d, interpol);
            const auto result = compress_saturated_s2u(interpol, zero_si256()).get_low();
            min(result, peak).store_nt(dst + x);
        } else {
            const Vec8f c = Vec8f().load_a(curAbove + x);
            const Vec8f d = (Vec8f().load_a(prev2 + x) + Vec8f().load_a(next2 + x)) * 0.5f;
            const Vec8f e = Vec8f().load_a(curBelow + x);
            const Vec8f temporal_diff0 = abs(Vec8f().load_a(prev2 + x) - Vec8f().load_a(next2 + x));
            const Vec8f temporal_diff1 = (abs(Vec8f().load_a(prevAbove + x) - c) + abs(Vec8f().load_a(prevBelow + x) - e)) * 0.5f;
            const Vec8f temporal_diff2 = (abs(Vec8f().load_a(nextAbove + x) - c) + abs(Vec8f().load_a(nextBelow + x) - e)) * 0.5f;
            const Vec8f temporal_diff = max(max(temporal_diff0 * 0.5f, temporal_diff1), temporal_diff2);

            Vec8f diff;
            if constexpr (spat) {
                const Vec8f b = ((Vec8f().load_a(prev2Above2 + x) + Vec8f().load_a(next2Above2 + x)) * 0.5f) - c;
                const Vec8f f = ((Vec8f().load_a(prev2Below2 + x) + Vec8f().load_a(next2Below2 + x)) * 0.5f) - e;
                const Vec8f dc = d - c;
                const Vec8f de = d - e;
                const Vec8f maximum = max(max(de, dc), min(b, f));
                const Vec8f minimum = min(min(de, dc), max(b, f));
                diff = max(max(temporal_diff, minimum), -maximum);
            }

            const Vec8f interpol = min(max((c + e) * 0.5f, d - diff), d + diff);
            select(temporal_diff == 0, d, interpol).store_nt(dst + x);
        }
    }
}

template<typename pixel_t>
void filterLine_avx2(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width,
                     const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept {
    const pixel_t * prev2 = reinterpret_cast<const pixel_t *>(_prev2);
    const pixel_t * prev = reinterpret_cast<const pixel_t *>(_prev);
    const pixel_t * cur = reinterpret_cast<const pixel_t *>(_cur);
    const pixel_t * next = reinterpret_cast<const pixel_t *>(_next);
    const pixel_t * next2 = reinterpret_cast<const pixel_t *>(_next2);
    pixel_t * dst = reinterpret_cast<pixel_t *>(_dst);

    const pixel_t * prev2Above4 = prev2 - stride4;
    const pixel_t * prev2Above2 = prev2 - stride2;
    const pixel_t * prev2Below2 = prev2 + stride2;
    const pixel_t * prev2Below4 = prev2 + stride4;
    const pixel_t * prevAbove = prev - stride;
    const pixel_t * prevBelow = prev + stride;
    const pixel_t * curAbove3 = cur - stride3;
    const pixel_t * curAbove = cur - stride;
    const pixel_t * curBelow = cur + stride;
    const pixel_t * curBelow3 = cur + stride3;
    const pixel_t * nextAbove = next - stride;
    const pixel_t * nextBelow = next + stride;
    const pixel_t * next2Above4 = next2 - stride4;
    const pixel_t * next2Above2 = next2 - stride2;
    const pixel_t * next2Below2 = next2 + stride2;
    const pixel_t * next2Below4 = next2 + stride4;

    for (int x = 0; x < width; x += step) {
        if constexpr (std::is_same_v<pixel_t, uint8_t>) {
            const Vec8i c = Vec8i().load_8uc(curAbove + x);
            const Vec8i d = (Vec8i().load_8uc(prev2 + x) + Vec8i().load_8uc(next2 + x)) >> 1;
            const Vec8i e = Vec8i().load_8uc(curBelow + x);
            const Vec8i temporal_diff0 = abs(Vec8i().load_8uc(prev2 + x) - Vec8i().load_8uc(next2 + x));
            const Vec8i temporal_diff1 = (abs(Vec8i().load_8uc(prevAbove + x) - c) + abs(Vec8i().load_8uc(prevBelow + x) - e)) >> 1;
            const Vec8i temporal_diff2 = (abs(Vec8i().load_8uc(nextAbove + x) - c) + abs(Vec8i().load_8uc(nextBelow + x) - e)) >> 1;
            const Vec8i temporal_diff = max(max(temporal_diff0 >> 1, temporal_diff1), temporal_diff2);

            const Vec8i b = ((Vec8i().load_8uc(prev2Above2 + x) + Vec8i().load_8uc(next2Above2 + x)) >> 1) - c;
            const Vec8i f = ((Vec8i().load_8uc(prev2Below2 + x) + Vec8i().load_8uc(next2Below2 + x)) >> 1) - e;
            const Vec8i dc = d - c;
            const Vec8i de = d - e;
            const Vec8i maximum = max(max(de, dc), min(b, f));
            const Vec8i minimum = min(min(de, dc), max(b, f));
            const Vec8i diff = max(max(temporal_diff, minimum), -maximum);

            const Vec8i interpol1 = (((coef_hf[0] * (Vec8i().load_8uc(prev2 + x) + Vec8i().load_8uc(next2 + x))
                                        - coef_hf[1] * (Vec8i().load_8uc(prev2Above2 + x) + Vec8i().load_8uc(next2Above2 + x) + Vec8i().load_8uc(prev2Below2 + x) + Vec8i().load_8uc(next2Below2 + x))
                                        + coef_hf[2] * (Vec8i().load_8uc(prev2Above4 + x) + Vec8i().load_8uc(next2Above4 + x) + Vec8i().load_8uc(prev2Below4 + x) + Vec8i().load_8uc(next2Below4 + x))) >> 2)
                                      + coef_lf[0] * (c + e) - coef_lf[1] * (Vec8i().load_8uc(curAbove3 + x) + Vec8i().load_8uc(curBelow3 + x))) >> 13;
            const Vec8i interpol2 = (coef_sp[0] * (c + e) - coef_sp[1] * (Vec8i().load_8uc(curAbove3 + x) + Vec8i().load_8uc(curBelow3 + x))) >> 13;

            Vec8i interpol = select(abs(c - e) > temporal_diff0, interpol1, interpol2);
            interpol = min(max(interpol, d - diff), d + diff);
            interpol = select(temporal_diff == 0, d, interpol);
            const auto result = compress_saturated_s2u(compress_saturated(interpol, zero_si256()), zero_si256()).get_low();
            result.storel(dst + x);
        } else if constexpr (std::is_same_v<pixel_t, uint16_t>) {
            const Vec8i c = Vec8i().load_8us(curAbove + x);
            const Vec8i d = (Vec8i().load_8us(prev2 + x) + Vec8i().load_8us(next2 + x)) >> 1;
            const Vec8i e = Vec8i().load_8us(curBelow + x);
            const Vec8i temporal_diff0 = abs(Vec8i().load_8us(prev2 + x) - Vec8i().load_8us(next2 + x));
            const Vec8i temporal_diff1 = (abs(Vec8i().load_8us(prevAbove + x) - c) + abs(Vec8i().load_8us(prevBelow + x) - e)) >> 1;
            const Vec8i temporal_diff2 = (abs(Vec8i().load_8us(nextAbove + x) - c) + abs(Vec8i().load_8us(nextBelow + x) - e)) >> 1;
            const Vec8i temporal_diff = max(max(temporal_diff0 >> 1, temporal_diff1), temporal_diff2);

            const Vec8i b = ((Vec8i().load_8us(prev2Above2 + x) + Vec8i().load_8us(next2Above2 + x)) >> 1) - c;
            const Vec8i f = ((Vec8i().load_8us(prev2Below2 + x) + Vec8i().load_8us(next2Below2 + x)) >> 1) - e;
            const Vec8i dc = d - c;
            const Vec8i de = d - e;
            const Vec8i maximum = max(max(de, dc), min(b, f));
            const Vec8i minimum = min(min(de, dc), max(b, f));
            const Vec8i diff = max(max(temporal_diff, minimum), -maximum);

            const Vec8i interpol1 = (((coef_hf[0] * (Vec8i().load_8us(prev2 + x) + Vec8i().load_8us(next2 + x))
                                       - coef_hf[1] * (Vec8i().load_8us(prev2Above2 + x) + Vec8i().load_8us(next2Above2 + x) + Vec8i().load_8us(prev2Below2 + x) + Vec8i().load_8us(next2Below2 + x))
                                       + coef_hf[2] * (Vec8i().load_8us(prev2Above4 + x) + Vec8i().load_8us(next2Above4 + x) + Vec8i().load_8us(prev2Below4 + x) + Vec8i().load_8us(next2Below4 + x))) >> 2)
                                     + coef_lf[0] * (c + e) - coef_lf[1] * (Vec8i().load_8us(curAbove3 + x) + Vec8i().load_8us(curBelow3 + x))) >> 13;
            const Vec8i interpol2 = (coef_sp[0] * (c + e) - coef_sp[1] * (Vec8i().load_8us(curAbove3 + x) + Vec8i().load_8us(curBelow3 + x))) >> 13;

            Vec8i interpol = select(abs(c - e) > temporal_diff0, interpol1, interpol2);
            interpol = min(max(interpol, d - diff), d + diff);
            interpol = select(temporal_diff == 0, d, interpol);
            const auto result = compress_saturated_s2u(interpol, zero_si256()).get_low();
            min(result, peak).store_nt(dst + x);
        } else {
            const Vec8f c = Vec8f().load_a(curAbove + x);
            const Vec8f d = (Vec8f().load_a(prev2 + x) + Vec8f().load_a(next2 + x)) * 0.5f;
            const Vec8f e = Vec8f().load_a(curBelow + x);
            const Vec8f temporal_diff0 = abs(Vec8f().load_a(prev2 + x) - Vec8f().load_a(next2 + x));
            const Vec8f temporal_diff1 = (abs(Vec8f().load_a(prevAbove + x) - c) + abs(Vec8f().load_a(prevBelow + x) - e)) * 0.5f;
            const Vec8f temporal_diff2 = (abs(Vec8f().load_a(nextAbove + x) - c) + abs(Vec8f().load_a(nextBelow + x) - e)) * 0.5f;
            const Vec8f temporal_diff = max(max(temporal_diff0 * 0.5f, temporal_diff1), temporal_diff2);

            const Vec8f b = ((Vec8f().load_a(prev2Above2 + x) + Vec8f().load_a(next2Above2 + x)) * 0.5f) - c;
            const Vec8f f = ((Vec8f().load_a(prev2Below2 + x) + Vec8f().load_a(next2Below2 + x)) * 0.5f) - e;
            const Vec8f dc = d - c;
            const Vec8f de = d - e;
            const Vec8f maximum = max(max(de, dc), min(b, f));
            const Vec8f minimum = min(min(de, dc), max(b, f));
            const Vec8f diff = max(max(temporal_diff, minimum), -maximum);

            const Vec8f interpol1 = ((coef_hf_f[0] * (Vec8f().load_a(prev2 + x) + Vec8f().load_a(next2 + x))
                                      - coef_hf_f[1] * (Vec8f().load_a(prev2Above2 + x) + Vec8f().load_a(next2Above2 + x) + Vec8f().load_a(prev2Below2 + x) + Vec8f().load_a(next2Below2 + x))
                                      + coef_hf_f[2] * (Vec8f().load_a(prev2Above4 + x) + Vec8f().load_a(next2Above4 + x) + Vec8f().load_a(prev2Below4 + x) + Vec8f().load_a(next2Below4 + x))) * 0.25f
                                     + coef_lf_f[0] * (c + e) - coef_lf_f[1] * (Vec8f().load_a(curAbove3 + x) + Vec8f().load_a(curBelow3 + x)));
            const Vec8f interpol2 = coef_sp_f[0] * (c + e) - coef_sp_f[1] * (Vec8f().load_a(curAbove3 + x) + Vec8f().load_a(curBelow3 + x));

            Vec8f interpol = select(abs(c - e) > temporal_diff0, interpol1, interpol2);
            interpol = min(max(interpol, d - diff), d + diff);
            select(temporal_diff == 0, d, interpol).store_nt(dst + x);
        }
    }
}

template void filterEdge_avx2<uint8_t, true>(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
template void filterEdge_avx2<uint8_t, false>(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
template void filterEdge_avx2<uint16_t, true>(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
template void filterEdge_avx2<uint16_t, false>(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
template void filterEdge_avx2<float, true>(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
template void filterEdge_avx2<float, false>(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;

template void filterLine_avx2<uint8_t>(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
template void filterLine_avx2<uint16_t>(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
template void filterLine_avx2<float>(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
#endif
