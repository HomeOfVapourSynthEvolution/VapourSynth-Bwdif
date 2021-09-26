/*
 * BobWeaver Deinterlacing Filter
 * Copyright (C) 2016      Thomas Mundt <loudmax@yahoo.de>
 * Copyright (C) 2020-2021 HolyWu
 *
 * Based on YADIF (Yet Another Deinterlacing Filter)
 * Copyright (C) 2006-2011 Michael Niedermayer <michaelni@gmx.at>
 *               2010      James Darnley <james.darnley@gmail.com>
 *
 * With use of Weston 3 Field Deinterlacing Filter algorithm
 * Copyright (C) 2012 British Broadcasting Corporation, All Rights Reserved
 * Author of de-interlace algorithm: Jim Easterbrook for BBC R&D
 * Based on the process described by Martin Weston for BBC R&D
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <VapourSynth4.h>
#include <VSHelper4.h>
#include <VSConstants4.h>

#include "Bwdif.h"

using namespace std::literals;

#ifdef BWDIF_X86
template<typename pixel_t, bool spat, bool hasEdeint> extern void filterEdge_sse2(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const ptrdiff_t positiveStride, const ptrdiff_t negativeStride, const ptrdiff_t stride2, const int step) noexcept;
template<typename pixel_t, bool spat, bool hasEdeint> extern void filterEdge_avx2(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const ptrdiff_t positiveStride, const ptrdiff_t negativeStride, const ptrdiff_t stride2, const int step) noexcept;
template<typename pixel_t, bool spat, bool hasEdeint> extern void filterEdge_avx512(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const ptrdiff_t positiveStride, const ptrdiff_t negativeStride, const ptrdiff_t stride2, const int step) noexcept;

template<typename pixel_t, bool hasEdeint> extern void filterLine_sse2(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const ptrdiff_t stride, const ptrdiff_t stride2, const ptrdiff_t stride3, const ptrdiff_t stride4, const int step, const int peak) noexcept;
template<typename pixel_t, bool hasEdeint> extern void filterLine_avx2(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const ptrdiff_t stride, const ptrdiff_t stride2, const ptrdiff_t stride3, const ptrdiff_t stride4, const int step, const int peak) noexcept;
template<typename pixel_t, bool hasEdeint> extern void filterLine_avx512(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const ptrdiff_t stride, const ptrdiff_t stride2, const ptrdiff_t stride3, const ptrdiff_t stride4, const int step, const int peak) noexcept;
#endif

struct BwdifData final {
    VSNode* node;
    VSNode* edeint;
    VSVideoInfo vi;
    const VSVideoInfo* viSaved;
    int field;
    int edgeStep;
    int lineStep;
    int peak;
    void (*filter)(const VSFrame* prevFrame, const VSFrame* curFrame, const VSFrame* nextFrame, const VSFrame* edeintFrame, VSFrame* dstFrame, const int field, const BwdifData* const VS_RESTRICT d, const VSAPI* vsapi) noexcept;
    void (*filterEdgeWithSpat)(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const ptrdiff_t positiveStride, const ptrdiff_t negativeStride, const ptrdiff_t stride2, const int step) noexcept;
    void (*filterEdgeWithoutSpat)(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const ptrdiff_t positiveStride, const ptrdiff_t negativeStride, const ptrdiff_t stride2, const int step) noexcept;
    void (*filterLine)(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst, const int width, const ptrdiff_t stride, const ptrdiff_t stride2, const ptrdiff_t stride3, const ptrdiff_t stride4, const int step, const int peak) noexcept;
};

template<typename pixel_t, bool spat, bool hasEdeint>
static void filterEdge_c(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst,
                         const int width, const ptrdiff_t positiveStride, const ptrdiff_t negativeStride, const ptrdiff_t stride2, [[maybe_unused]] const int step) noexcept {
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

    for (auto x{ 0 }; x < width; x++) {
        auto c{ curAbove[x] };
        auto d{ div2<pixel_t>(prev2[x] + next2[x]) };
        auto e{ curBelow[x] };
        auto temporal_diff0{ std::abs(prev2[x] - next2[x]) };
        auto temporal_diff1{ div2<pixel_t>(std::abs(prevAbove[x] - c) + std::abs(prevBelow[x] - e)) };
        auto temporal_diff2{ div2<pixel_t>(std::abs(nextAbove[x] - c) + std::abs(nextBelow[x] - e)) };
        auto diff{ std::max({ div2<pixel_t>(temporal_diff0), temporal_diff1, temporal_diff2 }) };

        if (!diff) {
            dst[x] = static_cast<pixel_t>(d);
        } else {
            if constexpr (spat) {
                auto b{ div2<pixel_t>(prev2Above2[x] + next2Above2[x]) - c };
                auto f{ div2<pixel_t>(prev2Below2[x] + next2Below2[x]) - e };
                auto dc{ d - c };
                auto de{ d - e };
                auto maximum{ std::max({ de, dc, std::min(b, f) }) };
                auto minimum{ std::min({ de, dc, std::max(b, f) }) };
                diff = std::max({ diff, minimum, -maximum });
            }

            decltype(d) interpol;
            if constexpr (hasEdeint)
                interpol = edeint[x];
            else
                interpol = div2<pixel_t>(c + e);

            dst[x] = static_cast<pixel_t>(std::clamp(interpol, d - diff, d + diff));
        }
    }
}

template<typename pixel_t, bool hasEdeint>
static void filterLine_c(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, const void* _edeint, void* _dst,
                         const int width, const ptrdiff_t stride, const ptrdiff_t stride2, const ptrdiff_t stride3, const ptrdiff_t stride4,
                         [[maybe_unused]] const int step, const int peak) noexcept {
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

    for (auto x{ 0 }; x < width; x++) {
        auto c{ curAbove[x] };
        auto d{ div2<pixel_t>(prev2[x] + next2[x]) };
        auto e{ curBelow[x] };
        auto temporal_diff0{ std::abs(prev2[x] - next2[x]) };
        auto temporal_diff1{ div2<pixel_t>(std::abs(prevAbove[x] - c) + std::abs(prevBelow[x] - e)) };
        auto temporal_diff2{ div2<pixel_t>(std::abs(nextAbove[x] - c) + std::abs(nextBelow[x] - e)) };
        auto diff{ std::max({ div2<pixel_t>(temporal_diff0), temporal_diff1, temporal_diff2 }) };

        if (!diff) {
            dst[x] = static_cast<pixel_t>(d);
        } else {
            auto b{ div2<pixel_t>(prev2Above2[x] + next2Above2[x]) - c };
            auto f{ div2<pixel_t>(prev2Below2[x] + next2Below2[x]) - e };
            auto dc{ d - c };
            auto de{ d - e };
            auto maximum{ std::max({ de, dc, std::min(b, f) }) };
            auto minimum{ std::min({ de, dc, std::max(b, f) }) };
            diff = std::max({ diff, minimum, -maximum });

            decltype(d) interpol;
            if (std::abs(c - e) > temporal_diff0) {
                interpol = div4<pixel_t>(coefHF<pixel_t>()[0] * (prev2[x] + next2[x])
                                         - coefHF<pixel_t>()[1] * (prev2Above2[x] + next2Above2[x] + prev2Below2[x] + next2Below2[x])
                                         + coefHF<pixel_t>()[2] * (prev2Above4[x] + next2Above4[x] + prev2Below4[x] + next2Below4[x]))
                    + coefLF<pixel_t>()[0] * (c + e) - coefLF<pixel_t>()[1] * (curAbove3[x] + curBelow3[x]);
                if constexpr (std::is_integral_v<pixel_t>)
                    interpol >>= 13;
            } else {
                if constexpr (hasEdeint) {
                    interpol = edeint[x];
                } else {
                    interpol = coefSP<pixel_t>()[0] * (c + e) - coefSP<pixel_t>()[1] * (curAbove3[x] + curBelow3[x]);
                    if constexpr (std::is_integral_v<pixel_t>)
                        interpol >>= 13;
                }
            }

            interpol = std::clamp(interpol, d - diff, d + diff);
            if constexpr (std::is_integral_v<pixel_t>)
                interpol = std::clamp(interpol, 0, peak);

            dst[x] = static_cast<pixel_t>(interpol);
        }
    }
}

template<typename pixel_t>
static void filter(const VSFrame* prevFrame, const VSFrame* curFrame, const VSFrame* nextFrame, const VSFrame* edeintFrame, VSFrame* dstFrame,
                   const int field, const BwdifData* const VS_RESTRICT d, const VSAPI* vsapi) noexcept {
    for (auto plane{ 0 }; plane < d->vi.format.numPlanes; plane++) {
        const auto width{ vsapi->getFrameWidth(curFrame, plane) };
        const auto height{ vsapi->getFrameHeight(curFrame, plane) };
        const auto stride{ vsapi->getStride(curFrame, plane) / d->vi.format.bytesPerSample };
        auto prev{ reinterpret_cast<const pixel_t*>(vsapi->getReadPtr(prevFrame, plane)) };
        auto cur{ reinterpret_cast<const pixel_t*>(vsapi->getReadPtr(curFrame, plane)) };
        auto next{ reinterpret_cast<const pixel_t*>(vsapi->getReadPtr(nextFrame, plane)) };
        auto dst{ reinterpret_cast<pixel_t*>(vsapi->getWritePtr(dstFrame, plane)) };

        decltype(cur) edeint{ nullptr };
        if (edeintFrame)
            edeint = reinterpret_cast<const pixel_t*>(vsapi->getReadPtr(edeintFrame, plane));

        vsh::bitblt(dst + stride * (1 - field),
                    vsapi->getStride(dstFrame, plane) * 2,
                    cur + stride * (1 - field),
                    vsapi->getStride(curFrame, plane) * 2,
                    width * sizeof(pixel_t),
                    height / 2);

        prev += stride * field;
        cur += stride * field;
        next += stride * field;
        if (edeint)
            edeint += stride * field;
        dst += stride * field;

        auto prev2{ field ? prev : cur };
        auto next2{ field ? cur : next };

        for (auto y{ field }; y < height; y += 2) {
            if ((y < 4) || (y + 5 > height)) {
                if ((y < 2) || (y + 3 > height))
                    d->filterEdgeWithoutSpat(prev2, prev, cur, next, next2, edeint, dst,
                                             width,
                                             y + 1 < height ? stride : -stride,
                                             y > 0 ? -stride : stride,
                                             stride * 2,
                                             d->edgeStep);
                else
                    d->filterEdgeWithSpat(prev2, prev, cur, next, next2, edeint, dst,
                                          width,
                                          y + 1 < height ? stride : -stride,
                                          y > 0 ? -stride : stride,
                                          stride * 2,
                                          d->edgeStep);
            } else {
                d->filterLine(prev2, prev, cur, next, next2, edeint, dst,
                              width,
                              stride, stride * 2, stride * 3, stride * 4,
                              d->lineStep, d->peak);
            }

            prev2 += stride * 2;
            prev += stride * 2;
            cur += stride * 2;
            next += stride * 2;
            next2 += stride * 2;
            if (edeint)
                edeint += stride * 2;
            dst += stride * 2;
        }
    }
}

static const VSFrame* VS_CC bwdifGetFrame(int n, int activationReason, void* instanceData, [[maybe_unused]] void** frameData, VSFrameContext* frameCtx, VSCore* core, const VSAPI* vsapi) {
    auto d{ static_cast<const BwdifData*>(instanceData) };

    if (activationReason == arInitial) {
        auto nSaved{ n };
        if (d->field > 1)
            n /= 2;

        if (n > 0)
            vsapi->requestFrameFilter(n - 1, d->node, frameCtx);
        vsapi->requestFrameFilter(n, d->node, frameCtx);
        if (n < d->viSaved->numFrames - 1)
            vsapi->requestFrameFilter(n + 1, d->node, frameCtx);

        if (d->edeint)
            vsapi->requestFrameFilter(nSaved, d->edeint, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        auto nSaved{ n };
        auto field{ d->field };
        if (d->field > 1) {
            n /= 2;
            field -= 2;
        }

        auto prev{ vsapi->getFrameFilter(std::max(n - 1, 0), d->node, frameCtx) };
        auto cur{ vsapi->getFrameFilter(n, d->node, frameCtx) };
        auto next{ vsapi->getFrameFilter(std::min(n + 1, d->viSaved->numFrames - 1), d->node, frameCtx) };
        auto dst{ vsapi->newVideoFrame(&d->vi.format, d->vi.width, d->vi.height, cur, core) };

        decltype(cur) edeint{ nullptr };
        if (d->edeint)
            edeint = vsapi->getFrameFilter(nSaved, d->edeint, frameCtx);

        auto err{ 0 };
        auto fieldBased{ vsapi->mapGetIntSaturated(vsapi->getFramePropertiesRO(cur), "_FieldBased", 0, &err) };
        if (fieldBased == VSC_FIELD_BOTTOM)
            field = 0;
        else if (fieldBased == VSC_FIELD_TOP)
            field = 1;

        if (d->field > 1) {
            if (nSaved & 1)
                field = field == 0;
            else
                field = field == 1;
        }

        d->filter(prev, cur, next, edeint, dst, field, d, vsapi);

        auto props{ vsapi->getFramePropertiesRW(dst) };
        vsapi->mapSetInt(props, "_FieldBased", VSC_FIELD_PROGRESSIVE, maReplace);

        if (d->field > 1) {
            auto errNum{ 0 }, errDen{ 0 };
            auto durationNum{ vsapi->mapGetInt(props, "_DurationNum", 0, &errNum) };
            auto durationDen{ vsapi->mapGetInt(props, "_DurationDen", 0, &errDen) };
            if (!errNum && !errDen) {
                vsh::muldivRational(&durationNum, &durationDen, 1, 2);
                vsapi->mapSetInt(props, "_DurationNum", durationNum, maReplace);
                vsapi->mapSetInt(props, "_DurationDen", durationDen, maReplace);
            }
        }

        vsapi->freeFrame(prev);
        vsapi->freeFrame(cur);
        vsapi->freeFrame(next);
        vsapi->freeFrame(edeint);
        return dst;
    }

    return nullptr;
}

static void VS_CC bwdifFree(void* instanceData, [[maybe_unused]] VSCore* core, const VSAPI* vsapi) {
    auto d{ static_cast<BwdifData*>(instanceData) };
    vsapi->freeNode(d->node);
    vsapi->freeNode(d->edeint);
    delete d;
}

static void VS_CC bwdifCreate(const VSMap* in, VSMap* out, [[maybe_unused]] void* userData, VSCore* core, const VSAPI* vsapi) {
    auto d{ std::make_unique<BwdifData>() };

    try {
        d->node = vsapi->mapGetNode(in, "clip", 0, nullptr);
        d->vi = *vsapi->getVideoInfo(d->node);
        d->viSaved = vsapi->getVideoInfo(d->node);
        auto err{ 0 };

        if (!vsh::isConstantVideoFormat(&d->vi) ||
            (d->vi.format.sampleType == stInteger && d->vi.format.bitsPerSample > 16) ||
            (d->vi.format.sampleType == stFloat && d->vi.format.bitsPerSample != 32))
            throw "only constant format 8-16 bit integer and 32 bit float input supported";

        if (d->vi.height < 4)
            throw "height must be at least 4";

        d->field = vsapi->mapGetIntSaturated(in, "field", 0, nullptr);
        d->edeint = vsapi->mapGetNode(in, "edeint", 0, &err);
        auto opt{ vsapi->mapGetIntSaturated(in, "opt", 0, &err) };

        if (d->field < 0 || d->field > 3)
            throw "field must be 0, 1, 2, or 3";

        if (opt < 0 || opt > 4)
            throw "opt must be 0, 1, 2, 3, or 4";

        {
            if (d->vi.format.bytesPerSample == 1) {
                d->filter = filter<uint8_t>;
                if (d->edeint) {
                    d->filterEdgeWithSpat = filterEdge_c<uint8_t, true, true>;
                    d->filterEdgeWithoutSpat = filterEdge_c<uint8_t, false, true>;
                    d->filterLine = filterLine_c<uint8_t, true>;
                } else {
                    d->filterEdgeWithSpat = filterEdge_c<uint8_t, true, false>;
                    d->filterEdgeWithoutSpat = filterEdge_c<uint8_t, false, false>;
                    d->filterLine = filterLine_c<uint8_t, false>;
                }
            } else if (d->vi.format.bytesPerSample == 2) {
                d->filter = filter<uint16_t>;
                if (d->edeint) {
                    d->filterEdgeWithSpat = filterEdge_c<uint16_t, true, true>;
                    d->filterEdgeWithoutSpat = filterEdge_c<uint16_t, false, true>;
                    d->filterLine = filterLine_c<uint16_t, true>;
                } else {
                    d->filterEdgeWithSpat = filterEdge_c<uint16_t, true, false>;
                    d->filterEdgeWithoutSpat = filterEdge_c<uint16_t, false, false>;
                    d->filterLine = filterLine_c<uint16_t, false>;
                }
            } else {
                d->filter = filter<float>;
                if (d->edeint) {
                    d->filterEdgeWithSpat = filterEdge_c<float, true, true>;
                    d->filterEdgeWithoutSpat = filterEdge_c<float, false, true>;
                    d->filterLine = filterLine_c<float, true>;
                } else {
                    d->filterEdgeWithSpat = filterEdge_c<float, true, false>;
                    d->filterEdgeWithoutSpat = filterEdge_c<float, false, false>;
                    d->filterLine = filterLine_c<float, false>;
                }
            }

#ifdef BWDIF_X86
            auto iset{ instrset_detect() };
            if ((opt == 0 && iset >= 10) || opt == 4) {
                if (d->vi.format.bytesPerSample == 1) {
                    if (d->edeint) {
                        d->filterEdgeWithSpat = filterEdge_avx512<uint8_t, true, true>;
                        d->filterEdgeWithoutSpat = filterEdge_avx512<uint8_t, false, true>;
                        d->filterLine = filterLine_avx512<uint8_t, true>;
                    } else {
                        d->filterEdgeWithSpat = filterEdge_avx512<uint8_t, true, false>;
                        d->filterEdgeWithoutSpat = filterEdge_avx512<uint8_t, false, false>;
                        d->filterLine = filterLine_avx512<uint8_t, false>;
                    }
                    d->edgeStep = 32;
                } else if (d->vi.format.bytesPerSample == 2) {
                    if (d->edeint) {
                        d->filterEdgeWithSpat = filterEdge_avx512<uint16_t, true, true>;
                        d->filterEdgeWithoutSpat = filterEdge_avx512<uint16_t, false, true>;
                        d->filterLine = filterLine_avx512<uint16_t, true>;
                    } else {
                        d->filterEdgeWithSpat = filterEdge_avx512<uint16_t, true, false>;
                        d->filterEdgeWithoutSpat = filterEdge_avx512<uint16_t, false, false>;
                        d->filterLine = filterLine_avx512<uint16_t, false>;
                    }
                    d->edgeStep = 16;
                } else {
                    if (d->edeint) {
                        d->filterEdgeWithSpat = filterEdge_avx512<float, true, true>;
                        d->filterEdgeWithoutSpat = filterEdge_avx512<float, false, true>;
                        d->filterLine = filterLine_avx512<float, true>;
                    } else {
                        d->filterEdgeWithSpat = filterEdge_avx512<float, true, false>;
                        d->filterEdgeWithoutSpat = filterEdge_avx512<float, false, false>;
                        d->filterLine = filterLine_avx512<float, false>;
                    }
                    d->edgeStep = 16;
                }
                d->lineStep = 16;
            } else if ((opt == 0 && iset >= 8) || opt == 3) {
                if (d->vi.format.bytesPerSample == 1) {
                    if (d->edeint) {
                        d->filterEdgeWithSpat = filterEdge_avx2<uint8_t, true, true>;
                        d->filterEdgeWithoutSpat = filterEdge_avx2<uint8_t, false, true>;
                        d->filterLine = filterLine_avx2<uint8_t, true>;
                    } else {
                        d->filterEdgeWithSpat = filterEdge_avx2<uint8_t, true, false>;
                        d->filterEdgeWithoutSpat = filterEdge_avx2<uint8_t, false, false>;
                        d->filterLine = filterLine_avx2<uint8_t, false>;
                    }
                    d->edgeStep = 16;
                } else if (d->vi.format.bytesPerSample == 2) {
                    if (d->edeint) {
                        d->filterEdgeWithSpat = filterEdge_avx2<uint16_t, true, true>;
                        d->filterEdgeWithoutSpat = filterEdge_avx2<uint16_t, false, true>;
                        d->filterLine = filterLine_avx2<uint16_t, true>;
                    } else {
                        d->filterEdgeWithSpat = filterEdge_avx2<uint16_t, true, false>;
                        d->filterEdgeWithoutSpat = filterEdge_avx2<uint16_t, false, false>;
                        d->filterLine = filterLine_avx2<uint16_t, false>;
                    }
                    d->edgeStep = 8;
                } else {
                    if (d->edeint) {
                        d->filterEdgeWithSpat = filterEdge_avx2<float, true, true>;
                        d->filterEdgeWithoutSpat = filterEdge_avx2<float, false, true>;
                        d->filterLine = filterLine_avx2<float, true>;
                    } else {
                        d->filterEdgeWithSpat = filterEdge_avx2<float, true, false>;
                        d->filterEdgeWithoutSpat = filterEdge_avx2<float, false, false>;
                        d->filterLine = filterLine_avx2<float, false>;
                    }
                    d->edgeStep = 8;
                }
                d->lineStep = 8;
            } else if ((opt == 0 && iset >= 2) || opt == 2) {
                if (d->vi.format.bytesPerSample == 1) {
                    if (d->edeint) {
                        d->filterEdgeWithSpat = filterEdge_sse2<uint8_t, true, true>;
                        d->filterEdgeWithoutSpat = filterEdge_sse2<uint8_t, false, true>;
                        d->filterLine = filterLine_sse2<uint8_t, true>;
                    } else {
                        d->filterEdgeWithSpat = filterEdge_sse2<uint8_t, true, false>;
                        d->filterEdgeWithoutSpat = filterEdge_sse2<uint8_t, false, false>;
                        d->filterLine = filterLine_sse2<uint8_t, false>;
                    }
                    d->edgeStep = 8;
                } else if (d->vi.format.bytesPerSample == 2) {
                    if (d->edeint) {
                        d->filterEdgeWithSpat = filterEdge_sse2<uint16_t, true, true>;
                        d->filterEdgeWithoutSpat = filterEdge_sse2<uint16_t, false, true>;
                        d->filterLine = filterLine_sse2<uint16_t, true>;
                    } else {
                        d->filterEdgeWithSpat = filterEdge_sse2<uint16_t, true, false>;
                        d->filterEdgeWithoutSpat = filterEdge_sse2<uint16_t, false, false>;
                        d->filterLine = filterLine_sse2<uint16_t, false>;
                    }
                    d->edgeStep = 4;
                } else {
                    if (d->edeint) {
                        d->filterEdgeWithSpat = filterEdge_sse2<float, true, true>;
                        d->filterEdgeWithoutSpat = filterEdge_sse2<float, false, true>;
                        d->filterLine = filterLine_sse2<float, true>;
                    } else {
                        d->filterEdgeWithSpat = filterEdge_sse2<float, true, false>;
                        d->filterEdgeWithoutSpat = filterEdge_sse2<float, false, false>;
                        d->filterLine = filterLine_sse2<float, false>;
                    }
                    d->edgeStep = 4;
                }
                d->lineStep = 4;
            }
#endif
        }

        if (d->field > 1) {
            if (d->vi.numFrames > INT_MAX / 2)
                throw "resulting clip is too long";
            d->vi.numFrames *= 2;

            vsh::muldivRational(&d->vi.fpsNum, &d->vi.fpsDen, 2, 1);
        }

        if (d->edeint) {
            if (!vsh::isSameVideoInfo(vsapi->getVideoInfo(d->edeint), &d->vi))
                throw "edeint clip must have the same format and dimensions as main clip";

            if (vsapi->getVideoInfo(d->edeint)->numFrames != d->vi.numFrames)
                throw "edeint clip's number of frames does not match";
        }

        if (d->vi.format.sampleType == stInteger)
            d->peak = (1 << d->vi.format.bitsPerSample) - 1;
    } catch (const char* error) {
        vsapi->mapSetError(out, ("Bwdif: "s + error).c_str());
        vsapi->freeNode(d->node);
        vsapi->freeNode(d->edeint);
        return;
    }

    std::vector<VSFilterDependency> deps;
    deps.push_back({ d->node, rpGeneral });
    if (d->edeint)
        deps.push_back({ d->edeint, rpStrictSpatial });

    vsapi->createVideoFilter(out, "Bwdif", &d->vi, bwdifGetFrame, bwdifFree, fmParallel, deps.data(), deps.size(), d.get(), core);
    d.release();
}

//////////////////////////////////////////
// Init

VS_EXTERNAL_API(void) VapourSynthPluginInit2(VSPlugin* plugin, const VSPLUGINAPI* vspapi) {
    vspapi->configPlugin("com.holywu.bwdif", "bwdif", "BobWeaver Deinterlacing Filter", VS_MAKE_VERSION(4, 1), VAPOURSYNTH_API_VERSION, 0, plugin);
    vspapi->registerFunction("Bwdif",
                             "clip:vnode;"
                             "field:int;"
                             "edeint:vnode:opt;"
                             "opt:int:opt;",
                             "clip:vnode;",
                             bwdifCreate, nullptr, plugin);
}
