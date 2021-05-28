#ifdef BWDIF_X86
#include "Bwdif.h"

template<typename pixel_t, bool spat, bool hasEdeint>
void filterEdge_avx2(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst,
                     const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept {
    static constexpr auto load = [](auto srcp) noexcept {
        if constexpr (std::is_same_v<pixel_t, uint8_t>)
            return Vec16s().load_16uc(srcp);
        else if constexpr (std::is_same_v<pixel_t, uint16_t>)
            return Vec8i().load_8us(srcp);
        else
            return Vec8f().load_a(srcp);
    };

    static constexpr auto store = [](auto interpol, auto dstp) noexcept {
        if constexpr (std::is_integral_v<pixel_t>) {
            auto result{ compress_saturated_s2u(interpol, zero_si256()).get_low() };
            result.store_nt(dstp);
        } else {
            interpol.store_nt(dstp);
        }
    };

    auto prev2{ reinterpret_cast<const pixel_t*>(_prev2) };
    auto prev{ reinterpret_cast<const pixel_t*>(_prev) };
    auto cur{ reinterpret_cast<const pixel_t*>(_cur) };
    auto next{ reinterpret_cast<const pixel_t*>(_next) };
    auto next2{ reinterpret_cast<const pixel_t*>(_next2) };
    auto edeint{ reinterpret_cast<const pixel_t*>(_edeint) };
    auto dst{ reinterpret_cast<pixel_t*>(_dst) };

    auto prev2Above2{ prev2 - stride2 };
    auto prev2Below2{ prev2 + stride2 };
    auto prevAbove{ prev + negativeStride };
    auto prevBelow{ prev + positiveStride };
    auto curAbove{ cur + negativeStride };
    auto curBelow{ cur + positiveStride };
    auto nextAbove{ next + negativeStride };
    auto nextBelow{ next + positiveStride };
    auto next2Above2{ next2 - stride2 };
    auto next2Below2{ next2 + stride2 };

    for (auto x{ 0 }; x < width; x += step) {
        auto c{ load(curAbove + x) };
        auto d{ div2<pixel_t>(load(prev2 + x) + load(next2 + x)) };
        auto e{ load(curBelow + x) };
        auto temporal_diff0{ abs(load(prev2 + x) - load(next2 + x)) };
        auto temporal_diff1{ div2<pixel_t>(abs(load(prevAbove + x) - c) + abs(load(prevBelow + x) - e)) };
        auto temporal_diff2{ div2<pixel_t>(abs(load(nextAbove + x) - c) + abs(load(nextBelow + x) - e)) };
        auto temporal_diff{ max(max(div2<pixel_t>(temporal_diff0), temporal_diff1), temporal_diff2) };

        auto diff{ temporal_diff };
        if constexpr (spat) {
            auto b{ div2<pixel_t>(load(prev2Above2 + x) + load(next2Above2 + x)) - c };
            auto f{ div2<pixel_t>(load(prev2Below2 + x) + load(next2Below2 + x)) - e };
            auto dc{ d - c };
            auto de{ d - e };
            auto maximum{ max(max(de, dc), min(b, f)) };
            auto minimum{ min(min(de, dc), max(b, f)) };
            diff = max(max(diff, minimum), -maximum);
        }

        decltype(d) interpol;
        if constexpr (hasEdeint)
            interpol = load(edeint + x);
        else
            interpol = div2<pixel_t>(c + e);
        interpol = min(max(interpol, d - diff), d + diff);
        interpol = select(temporal_diff == 0, d, interpol);

        store(interpol, dst + x);
    }
}

template<typename pixel_t, bool hasEdeint>
void filterLine_avx2(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst,
                     const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept {
    static constexpr auto load = [](auto srcp) noexcept {
        if constexpr (std::is_same_v<pixel_t, uint8_t>)
            return Vec8i().load_8uc(srcp);
        else if constexpr (std::is_same_v<pixel_t, uint16_t>)
            return Vec8i().load_8us(srcp);
        else
            return Vec8f().load_a(srcp);
    };

    static constexpr auto store = [](auto interpol, auto dstp) noexcept {
        if constexpr (std::is_same_v<pixel_t, uint8_t>) {
            auto result{ compress_saturated_s2u(compress_saturated(interpol, zero_si256()), zero_si256()).get_low() };
            result.storel(dstp);
        } else if constexpr (std::is_same_v<pixel_t, uint16_t>) {
            auto result{ compress_saturated_s2u(interpol, zero_si256()).get_low() };
            result.store_nt(dstp);
        } else {
            interpol.store_nt(dstp);
        }
    };

    auto prev2{ reinterpret_cast<const pixel_t*>(_prev2) };
    auto prev{ reinterpret_cast<const pixel_t*>(_prev) };
    auto cur{ reinterpret_cast<const pixel_t*>(_cur) };
    auto next{ reinterpret_cast<const pixel_t*>(_next) };
    auto next2{ reinterpret_cast<const pixel_t*>(_next2) };
    auto edeint{ reinterpret_cast<const pixel_t*>(_edeint) };
    auto dst{ reinterpret_cast<pixel_t*>(_dst) };

    auto prev2Above4{ prev2 - stride4 };
    auto prev2Above2{ prev2 - stride2 };
    auto prev2Below2{ prev2 + stride2 };
    auto prev2Below4{ prev2 + stride4 };
    auto prevAbove{ prev - stride };
    auto prevBelow{ prev + stride };
    auto curAbove3{ cur - stride3 };
    auto curAbove{ cur - stride };
    auto curBelow{ cur + stride };
    auto curBelow3{ cur + stride3 };
    auto nextAbove{ next - stride };
    auto nextBelow{ next + stride };
    auto next2Above4{ next2 - stride4 };
    auto next2Above2{ next2 - stride2 };
    auto next2Below2{ next2 + stride2 };
    auto next2Below4{ next2 + stride4 };

    for (auto x{ 0 }; x < width; x += step) {
        auto c{ load(curAbove + x) };
        auto d{ div2<pixel_t>(load(prev2 + x) + load(next2 + x)) };
        auto e{ load(curBelow + x) };
        auto temporal_diff0{ abs(load(prev2 + x) - load(next2 + x)) };
        auto temporal_diff1{ div2<pixel_t>(abs(load(prevAbove + x) - c) + abs(load(prevBelow + x) - e)) };
        auto temporal_diff2{ div2<pixel_t>(abs(load(nextAbove + x) - c) + abs(load(nextBelow + x) - e)) };
        auto temporal_diff{ max(max(div2<pixel_t>(temporal_diff0), temporal_diff1), temporal_diff2) };

        auto b{ div2<pixel_t>(load(prev2Above2 + x) + load(next2Above2 + x)) - c };
        auto f{ div2<pixel_t>(load(prev2Below2 + x) + load(next2Below2 + x)) - e };
        auto dc{ d - c };
        auto de{ d - e };
        auto maximum{ max(max(de, dc), min(b, f)) };
        auto minimum{ min(min(de, dc), max(b, f)) };
        auto diff{ max(max(temporal_diff, minimum), -maximum) };

        auto interpol1{ div4<pixel_t>(coefHF<pixel_t>()[0] * (load(prev2 + x) + load(next2 + x))
                                       - coefHF<pixel_t>()[1] * (load(prev2Above2 + x) + load(next2Above2 + x) + load(prev2Below2 + x) + load(next2Below2 + x))
                                       + coefHF<pixel_t>()[2] * (load(prev2Above4 + x) + load(next2Above4 + x) + load(prev2Below4 + x) + load(next2Below4 + x)))
            + coefLF<pixel_t>()[0] * (c + e) - coefLF<pixel_t>()[1] * (load(curAbove3 + x) + load(curBelow3 + x)) };
        if constexpr (std::is_integral_v<pixel_t>)
            interpol1 >>= 13;

        decltype(d) interpol2;
        if constexpr (hasEdeint) {
            interpol2 = load(edeint + x);
        } else {
            interpol2 = coefSP<pixel_t>()[0] * (c + e) - coefSP<pixel_t>()[1] * (load(curAbove3 + x) + load(curBelow3 + x));
            if constexpr (std::is_integral_v<pixel_t>)
                interpol2 >>= 13;
        }

        auto interpol{ select(abs(c - e) > temporal_diff0, interpol1, interpol2) };
        interpol = min(max(interpol, d - diff), d + diff);
        if constexpr (std::is_integral_v<pixel_t>)
            interpol = min(max(interpol, 0), peak);
        interpol = select(temporal_diff == 0, d, interpol);

        store(interpol, dst + x);
    }
}

template void filterEdge_avx2<uint8_t, true, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx2<uint8_t, true, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx2<uint8_t, false, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx2<uint8_t, false, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx2<uint16_t, true, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx2<uint16_t, true, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx2<uint16_t, false, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx2<uint16_t, false, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx2<float, true, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx2<float, true, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx2<float, false, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;
template void filterEdge_avx2<float, false, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step) noexcept;

template void filterLine_avx2<uint8_t, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
template void filterLine_avx2<uint8_t, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
template void filterLine_avx2<uint16_t, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
template void filterLine_avx2<uint16_t, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
template void filterLine_avx2<float, true>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
template void filterLine_avx2<float, false>(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
#endif
