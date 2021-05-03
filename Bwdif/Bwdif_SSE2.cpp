#ifdef BWDIF_X86
#include "Bwdif.h"

template<typename pixel_t, bool spat>
void filterEdge_sse2(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const int width,
                     const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept {
    const pixel_t* prev2 = reinterpret_cast<const pixel_t*>(_prev2);
    const pixel_t* prev = reinterpret_cast<const pixel_t*>(_prev);
    const pixel_t* cur = reinterpret_cast<const pixel_t*>(_cur);
    const pixel_t* next = reinterpret_cast<const pixel_t*>(_next);
    const pixel_t* next2 = reinterpret_cast<const pixel_t*>(_next2);
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
            const Vec8s c = Vec8s().load_8uc(curAbove + x);
            const Vec8s d = (Vec8s().load_8uc(prev2 + x) + Vec8s().load_8uc(next2 + x)) >> 1;
            const Vec8s e = Vec8s().load_8uc(curBelow + x);
            const Vec8s temporal_diff0 = abs(Vec8s().load_8uc(prev2 + x) - Vec8s().load_8uc(next2 + x));
            const Vec8s temporal_diff1 = (abs(Vec8s().load_8uc(prevAbove + x) - c) + abs(Vec8s().load_8uc(prevBelow + x) - e)) >> 1;
            const Vec8s temporal_diff2 = (abs(Vec8s().load_8uc(nextAbove + x) - c) + abs(Vec8s().load_8uc(nextBelow + x) - e)) >> 1;
            const Vec8s temporal_diff = max(max(temporal_diff0 >> 1, temporal_diff1), temporal_diff2);

            Vec8s diff;
            if constexpr (spat) {
                const Vec8s b = ((Vec8s().load_8uc(prev2Above2 + x) + Vec8s().load_8uc(next2Above2 + x)) >> 1) - c;
                const Vec8s f = ((Vec8s().load_8uc(prev2Below2 + x) + Vec8s().load_8uc(next2Below2 + x)) >> 1) - e;
                const Vec8s dc = d - c;
                const Vec8s de = d - e;
                const Vec8s maximum = max(max(de, dc), min(b, f));
                const Vec8s minimum = min(min(de, dc), max(b, f));
                diff = max(max(temporal_diff, minimum), -maximum);
            }

            Vec8s interpol = min(max((c + e) >> 1, d - diff), d + diff);
            interpol = select(temporal_diff == 0, d, interpol);
            const auto result = compress_saturated_s2u(interpol, zero_si128());
            result.storel(dst + x);
        } else if constexpr (std::is_same_v<pixel_t, uint16_t>) {
            const Vec4i c = Vec4i().load_4us(curAbove + x);
            const Vec4i d = (Vec4i().load_4us(prev2 + x) + Vec4i().load_4us(next2 + x)) >> 1;
            const Vec4i e = Vec4i().load_4us(curBelow + x);
            const Vec4i temporal_diff0 = abs(Vec4i().load_4us(prev2 + x) - Vec4i().load_4us(next2 + x));
            const Vec4i temporal_diff1 = (abs(Vec4i().load_4us(prevAbove + x) - c) + abs(Vec4i().load_4us(prevBelow + x) - e)) >> 1;
            const Vec4i temporal_diff2 = (abs(Vec4i().load_4us(nextAbove + x) - c) + abs(Vec4i().load_4us(nextBelow + x) - e)) >> 1;
            const Vec4i temporal_diff = max(max(temporal_diff0 >> 1, temporal_diff1), temporal_diff2);

            Vec4i diff;
            if constexpr (spat) {
                const Vec4i b = ((Vec4i().load_4us(prev2Above2 + x) + Vec4i().load_4us(next2Above2 + x)) >> 1) - c;
                const Vec4i f = ((Vec4i().load_4us(prev2Below2 + x) + Vec4i().load_4us(next2Below2 + x)) >> 1) - e;
                const Vec4i dc = d - c;
                const Vec4i de = d - e;
                const Vec4i maximum = max(max(de, dc), min(b, f));
                const Vec4i minimum = min(min(de, dc), max(b, f));
                diff = max(max(temporal_diff, minimum), -maximum);
            }

            Vec4i interpol = min(max((c + e) >> 1, d - diff), d + diff);
            interpol = select(temporal_diff == 0, d, interpol);
            const auto result = compress_saturated_s2u(interpol, zero_si128());
            min(result, peak).storel(dst + x);
        } else {
            const Vec4f c = Vec4f().load_a(curAbove + x);
            const Vec4f d = (Vec4f().load_a(prev2 + x) + Vec4f().load_a(next2 + x)) * 0.5f;
            const Vec4f e = Vec4f().load_a(curBelow + x);
            const Vec4f temporal_diff0 = abs(Vec4f().load_a(prev2 + x) - Vec4f().load_a(next2 + x));
            const Vec4f temporal_diff1 = (abs(Vec4f().load_a(prevAbove + x) - c) + abs(Vec4f().load_a(prevBelow + x) - e)) * 0.5f;
            const Vec4f temporal_diff2 = (abs(Vec4f().load_a(nextAbove + x) - c) + abs(Vec4f().load_a(nextBelow + x) - e)) * 0.5f;
            const Vec4f temporal_diff = max(max(temporal_diff0 * 0.5f, temporal_diff1), temporal_diff2);

            Vec4f diff;
            if constexpr (spat) {
                const Vec4f b = ((Vec4f().load_a(prev2Above2 + x) + Vec4f().load_a(next2Above2 + x)) * 0.5f) - c;
                const Vec4f f = ((Vec4f().load_a(prev2Below2 + x) + Vec4f().load_a(next2Below2 + x)) * 0.5f) - e;
                const Vec4f dc = d - c;
                const Vec4f de = d - e;
                const Vec4f maximum = max(max(de, dc), min(b, f));
                const Vec4f minimum = min(min(de, dc), max(b, f));
                diff = max(max(temporal_diff, minimum), -maximum);
            }

            const Vec4f interpol = min(max((c + e) * 0.5f, d - diff), d + diff);
            select(temporal_diff == 0, d, interpol).store_nt(dst + x);
        }
    }
}

template<typename pixel_t>
void filterLine_sse2(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const int width,
                     const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept {
    const pixel_t* prev2 = reinterpret_cast<const pixel_t*>(_prev2);
    const pixel_t* prev = reinterpret_cast<const pixel_t*>(_prev);
    const pixel_t* cur = reinterpret_cast<const pixel_t*>(_cur);
    const pixel_t* next = reinterpret_cast<const pixel_t*>(_next);
    const pixel_t* next2 = reinterpret_cast<const pixel_t*>(_next2);
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
            const Vec4i c = Vec4i().load_4uc(curAbove + x);
            const Vec4i d = (Vec4i().load_4uc(prev2 + x) + Vec4i().load_4uc(next2 + x)) >> 1;
            const Vec4i e = Vec4i().load_4uc(curBelow + x);
            const Vec4i temporal_diff0 = abs(Vec4i().load_4uc(prev2 + x) - Vec4i().load_4uc(next2 + x));
            const Vec4i temporal_diff1 = (abs(Vec4i().load_4uc(prevAbove + x) - c) + abs(Vec4i().load_4uc(prevBelow + x) - e)) >> 1;
            const Vec4i temporal_diff2 = (abs(Vec4i().load_4uc(nextAbove + x) - c) + abs(Vec4i().load_4uc(nextBelow + x) - e)) >> 1;
            const Vec4i temporal_diff = max(max(temporal_diff0 >> 1, temporal_diff1), temporal_diff2);

            const Vec4i b = ((Vec4i().load_4uc(prev2Above2 + x) + Vec4i().load_4uc(next2Above2 + x)) >> 1) - c;
            const Vec4i f = ((Vec4i().load_4uc(prev2Below2 + x) + Vec4i().load_4uc(next2Below2 + x)) >> 1) - e;
            const Vec4i dc = d - c;
            const Vec4i de = d - e;
            const Vec4i maximum = max(max(de, dc), min(b, f));
            const Vec4i minimum = min(min(de, dc), max(b, f));
            const Vec4i diff = max(max(temporal_diff, minimum), -maximum);

            const Vec4i interpol1 = (((coef_hf[0] * (Vec4i().load_4uc(prev2 + x) + Vec4i().load_4uc(next2 + x))
                                       - coef_hf[1] * (Vec4i().load_4uc(prev2Above2 + x) + Vec4i().load_4uc(next2Above2 + x) + Vec4i().load_4uc(prev2Below2 + x) + Vec4i().load_4uc(next2Below2 + x))
                                       + coef_hf[2] * (Vec4i().load_4uc(prev2Above4 + x) + Vec4i().load_4uc(next2Above4 + x) + Vec4i().load_4uc(prev2Below4 + x) + Vec4i().load_4uc(next2Below4 + x))) >> 2)
                                     + coef_lf[0] * (c + e) - coef_lf[1] * (Vec4i().load_4uc(curAbove3 + x) + Vec4i().load_4uc(curBelow3 + x))) >> 13;
            const Vec4i interpol2 = (coef_sp[0] * (c + e) - coef_sp[1] * (Vec4i().load_4uc(curAbove3 + x) + Vec4i().load_4uc(curBelow3 + x))) >> 13;

            Vec4i interpol = select(abs(c - e) > temporal_diff0, interpol1, interpol2);
            interpol = min(max(interpol, d - diff), d + diff);
            interpol = select(temporal_diff == 0, d, interpol);
            const auto result = compress_saturated_s2u(compress_saturated(interpol, zero_si128()), zero_si128());
            result.store_si32(dst + x);
        } else if constexpr (std::is_same_v<pixel_t, uint16_t>) {
            const Vec4i c = Vec4i().load_4us(curAbove + x);
            const Vec4i d = (Vec4i().load_4us(prev2 + x) + Vec4i().load_4us(next2 + x)) >> 1;
            const Vec4i e = Vec4i().load_4us(curBelow + x);
            const Vec4i temporal_diff0 = abs(Vec4i().load_4us(prev2 + x) - Vec4i().load_4us(next2 + x));
            const Vec4i temporal_diff1 = (abs(Vec4i().load_4us(prevAbove + x) - c) + abs(Vec4i().load_4us(prevBelow + x) - e)) >> 1;
            const Vec4i temporal_diff2 = (abs(Vec4i().load_4us(nextAbove + x) - c) + abs(Vec4i().load_4us(nextBelow + x) - e)) >> 1;
            const Vec4i temporal_diff = max(max(temporal_diff0 >> 1, temporal_diff1), temporal_diff2);

            const Vec4i b = ((Vec4i().load_4us(prev2Above2 + x) + Vec4i().load_4us(next2Above2 + x)) >> 1) - c;
            const Vec4i f = ((Vec4i().load_4us(prev2Below2 + x) + Vec4i().load_4us(next2Below2 + x)) >> 1) - e;
            const Vec4i dc = d - c;
            const Vec4i de = d - e;
            const Vec4i maximum = max(max(de, dc), min(b, f));
            const Vec4i minimum = min(min(de, dc), max(b, f));
            const Vec4i diff = max(max(temporal_diff, minimum), -maximum);

            const Vec4i interpol1 = (((coef_hf[0] * (Vec4i().load_4us(prev2 + x) + Vec4i().load_4us(next2 + x))
                                       - coef_hf[1] * (Vec4i().load_4us(prev2Above2 + x) + Vec4i().load_4us(next2Above2 + x) + Vec4i().load_4us(prev2Below2 + x) + Vec4i().load_4us(next2Below2 + x))
                                       + coef_hf[2] * (Vec4i().load_4us(prev2Above4 + x) + Vec4i().load_4us(next2Above4 + x) + Vec4i().load_4us(prev2Below4 + x) + Vec4i().load_4us(next2Below4 + x))) >> 2)
                                     + coef_lf[0] * (c + e) - coef_lf[1] * (Vec4i().load_4us(curAbove3 + x) + Vec4i().load_4us(curBelow3 + x))) >> 13;
            const Vec4i interpol2 = (coef_sp[0] * (c + e) - coef_sp[1] * (Vec4i().load_4us(curAbove3 + x) + Vec4i().load_4us(curBelow3 + x))) >> 13;

            Vec4i interpol = select(abs(c - e) > temporal_diff0, interpol1, interpol2);
            interpol = min(max(interpol, d - diff), d + diff);
            interpol = select(temporal_diff == 0, d, interpol);
            const auto result = compress_saturated_s2u(interpol, zero_si128());
            min(result, peak).storel(dst + x);
        } else {
            const Vec4f c = Vec4f().load_a(curAbove + x);
            const Vec4f d = (Vec4f().load_a(prev2 + x) + Vec4f().load_a(next2 + x)) * 0.5f;
            const Vec4f e = Vec4f().load_a(curBelow + x);
            const Vec4f temporal_diff0 = abs(Vec4f().load_a(prev2 + x) - Vec4f().load_a(next2 + x));
            const Vec4f temporal_diff1 = (abs(Vec4f().load_a(prevAbove + x) - c) + abs(Vec4f().load_a(prevBelow + x) - e)) * 0.5f;
            const Vec4f temporal_diff2 = (abs(Vec4f().load_a(nextAbove + x) - c) + abs(Vec4f().load_a(nextBelow + x) - e)) * 0.5f;
            const Vec4f temporal_diff = max(max(temporal_diff0 * 0.5f, temporal_diff1), temporal_diff2);

            const Vec4f b = ((Vec4f().load_a(prev2Above2 + x) + Vec4f().load_a(next2Above2 + x)) * 0.5f) - c;
            const Vec4f f = ((Vec4f().load_a(prev2Below2 + x) + Vec4f().load_a(next2Below2 + x)) * 0.5f) - e;
            const Vec4f dc = d - c;
            const Vec4f de = d - e;
            const Vec4f maximum = max(max(de, dc), min(b, f));
            const Vec4f minimum = min(min(de, dc), max(b, f));
            const Vec4f diff = max(max(temporal_diff, minimum), -maximum);

            const Vec4f interpol1 = ((coef_hf_f[0] * (Vec4f().load_a(prev2 + x) + Vec4f().load_a(next2 + x))
                                      - coef_hf_f[1] * (Vec4f().load_a(prev2Above2 + x) + Vec4f().load_a(next2Above2 + x) + Vec4f().load_a(prev2Below2 + x) + Vec4f().load_a(next2Below2 + x))
                                      + coef_hf_f[2] * (Vec4f().load_a(prev2Above4 + x) + Vec4f().load_a(next2Above4 + x) + Vec4f().load_a(prev2Below4 + x) + Vec4f().load_a(next2Below4 + x))) * 0.25f
                                     + coef_lf_f[0] * (c + e) - coef_lf_f[1] * (Vec4f().load_a(curAbove3 + x) + Vec4f().load_a(curBelow3 + x)));
            const Vec4f interpol2 = coef_sp_f[0] * (c + e) - coef_sp_f[1] * (Vec4f().load_a(curAbove3 + x) + Vec4f().load_a(curBelow3 + x));

            Vec4f interpol = select(abs(c - e) > temporal_diff0, interpol1, interpol2);
            interpol = min(max(interpol, d - diff), d + diff);
            select(temporal_diff == 0, d, interpol).store_nt(dst + x);
        }
    }
}

template void filterEdge_sse2<uint8_t, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
template void filterEdge_sse2<uint8_t, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
template void filterEdge_sse2<uint16_t, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
template void filterEdge_sse2<uint16_t, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
template void filterEdge_sse2<float, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
template void filterEdge_sse2<float, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;

template void filterLine_sse2<uint8_t>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
template void filterLine_sse2<uint16_t>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
template void filterLine_sse2<float>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
#endif
