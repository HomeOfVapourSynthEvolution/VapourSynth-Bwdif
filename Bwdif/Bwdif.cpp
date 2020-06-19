/*
 * BobWeaver Deinterlacing Filter
 * Copyright (C) 2016 Thomas Mundt <loudmax@yahoo.de>
 * Copyright (C) 2020 HolyWu
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

#include <VapourSynth.h>
#include <VSHelper.h>

#include "Bwdif.h"

#ifdef BWDIF_X86
template<typename pixel_t, bool spat> extern void filterEdge_sse2(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
template<typename pixel_t, bool spat> extern void filterEdge_avx2(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
template<typename pixel_t, bool spat> extern void filterEdge_avx512(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;

template<typename pixel_t> extern void filterLine_sse2(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
template<typename pixel_t> extern void filterLine_avx2(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
template<typename pixel_t> extern void filterLine_avx512(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
#endif

struct BwdifData final {
    VSNodeRef * node;
    VSVideoInfo vi;
    const VSVideoInfo * viSaved;
    int field;
    int edgeStep, lineStep, peak;
    void (*filterEdgeWithSpat)(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
    void (*filterEdgeWithoutSpat)(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
    void (*filterLine)(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
};

template<typename pixel_t, bool spat>
static inline void filterEdge_c(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width,
                                const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept {
    const pixel_t * prev2 = reinterpret_cast<const pixel_t *>(_prev2);
    const pixel_t * prev = reinterpret_cast<const pixel_t *>(_prev);
    const pixel_t * cur = reinterpret_cast<const pixel_t *>(_cur);
    const pixel_t * next = reinterpret_cast<const pixel_t *>(_next);
    const pixel_t * next2 = reinterpret_cast<const pixel_t *>(_next2);
    pixel_t * VS_RESTRICT dst = reinterpret_cast<pixel_t *>(_dst);

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

    for (int x = 0; x < width; x++) {
        if constexpr (std::is_integral_v<pixel_t>) {
            const int c = curAbove[x];
            const int d = (prev2[x] + next2[x]) >> 1;
            const int e = curBelow[x];
            const int temporal_diff0 = std::abs(prev2[x] - next2[x]);
            const int temporal_diff1 = (std::abs(prevAbove[x] - c) + std::abs(prevBelow[x] - e)) >> 1;
            const int temporal_diff2 = (std::abs(nextAbove[x] - c) + std::abs(nextBelow[x] - e)) >> 1;
            int diff = std::max({ temporal_diff0 >> 1, temporal_diff1, temporal_diff2 });

            if (!diff) {
                dst[x] = d;
            } else {
                if constexpr (spat) {
                    const int b = ((prev2Above2[x] + next2Above2[x]) >> 1) - c;
                    const int f = ((prev2Below2[x] + next2Below2[x]) >> 1) - e;
                    const int dc = d - c;
                    const int de = d - e;
                    const int maximum = std::max({ de, dc, std::min(b, f) });
                    const int minimum = std::min({ de, dc, std::max(b, f) });
                    diff = std::max({ diff, minimum, -maximum });
                }

                const int interpol = std::clamp((c + e) >> 1, d - diff, d + diff);
                dst[x] = std::clamp(interpol, 0, peak);
            }
        } else {
            const float c = curAbove[x];
            const float d = (prev2[x] + next2[x]) * 0.5f;
            const float e = curBelow[x];
            const float temporal_diff0 = std::abs(prev2[x] - next2[x]);
            const float temporal_diff1 = (std::abs(prevAbove[x] - c) + std::abs(prevBelow[x] - e)) * 0.5f;
            const float temporal_diff2 = (std::abs(nextAbove[x] - c) + std::abs(nextBelow[x] - e)) * 0.5f;
            float diff = std::max({ temporal_diff0 * 0.5f, temporal_diff1, temporal_diff2 });

            if (!diff) {
                dst[x] = d;
            } else {
                if constexpr (spat) {
                    const float b = ((prev2Above2[x] + next2Above2[x]) * 0.5f) - c;
                    const float f = ((prev2Below2[x] + next2Below2[x]) * 0.5f) - e;
                    const float dc = d - c;
                    const float de = d - e;
                    const float maximum = std::max({ de, dc, std::min(b, f) });
                    const float minimum = std::min({ de, dc, std::max(b, f) });
                    diff = std::max({ diff, minimum, -maximum });
                }

                dst[x] = std::clamp((c + e) * 0.5f, d - diff, d + diff);
            }
        }
    }
}

template<typename pixel_t>
static inline void filterLine_c(const void * _prev2, const void * _prev, const void * _cur, const void * _next, const void * _next2, void * _dst, const int width,
                                const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept {
    const pixel_t * prev2 = reinterpret_cast<const pixel_t *>(_prev2);
    const pixel_t * prev = reinterpret_cast<const pixel_t *>(_prev);
    const pixel_t * cur = reinterpret_cast<const pixel_t *>(_cur);
    const pixel_t * next = reinterpret_cast<const pixel_t *>(_next);
    const pixel_t * next2 = reinterpret_cast<const pixel_t *>(_next2);
    pixel_t * VS_RESTRICT dst = reinterpret_cast<pixel_t *>(_dst);

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

    for (int x = 0; x < width; x++) {
        if constexpr (std::is_integral_v<pixel_t>) {
            const int c = curAbove[x];
            const int d = (prev2[x] + next2[x]) >> 1;
            const int e = curBelow[x];
            const int temporal_diff0 = std::abs(prev2[x] - next2[x]);
            const int temporal_diff1 = (std::abs(prevAbove[x] - c) + std::abs(prevBelow[x] - e)) >> 1;
            const int temporal_diff2 = (std::abs(nextAbove[x] - c) + std::abs(nextBelow[x] - e)) >> 1;
            int diff = std::max({ temporal_diff0 >> 1, temporal_diff1, temporal_diff2 });

            if (!diff) {
                dst[x] = d;
            } else {
                const int b = ((prev2Above2[x] + next2Above2[x]) >> 1) - c;
                const int f = ((prev2Below2[x] + next2Below2[x]) >> 1) - e;
                const int dc = d - c;
                const int de = d - e;
                const int maximum = std::max({ de, dc, std::min(b, f) });
                const int minimum = std::min({ de, dc, std::max(b, f) });
                diff = std::max({ diff, minimum, -maximum });

                int interpol;
                if (std::abs(c - e) > temporal_diff0)
                    interpol = (((coef_hf[0] * (prev2[x] + next2[x])
                                  - coef_hf[1] * (prev2Above2[x] + next2Above2[x] + prev2Below2[x] + next2Below2[x])
                                  + coef_hf[2] * (prev2Above4[x] + next2Above4[x] + prev2Below4[x] + next2Below4[x])) >> 2)
                                + coef_lf[0] * (c + e) - coef_lf[1] * (curAbove3[x] + curBelow3[x])) >> 13;
                else
                    interpol = (coef_sp[0] * (c + e) - coef_sp[1] * (curAbove3[x] + curBelow3[x])) >> 13;

                interpol = std::clamp(interpol, d - diff, d + diff);
                dst[x] = std::clamp(interpol, 0, peak);
            }
        } else {
            const float c = curAbove[x];
            const float d = (prev2[x] + next2[x]) * 0.5f;
            const float e = curBelow[x];
            const float temporal_diff0 = std::abs(prev2[x] - next2[x]);
            const float temporal_diff1 = (std::abs(prevAbove[x] - c) + std::abs(prevBelow[x] - e)) * 0.5f;
            const float temporal_diff2 = (std::abs(nextAbove[x] - c) + std::abs(nextBelow[x] - e)) * 0.5f;
            float diff = std::max({ temporal_diff0 * 0.5f, temporal_diff1, temporal_diff2 });

            if (!diff) {
                dst[x] = d;
            } else {
                const float b = ((prev2Above2[x] + next2Above2[x]) * 0.5f) - c;
                const float f = ((prev2Below2[x] + next2Below2[x]) * 0.5f) - e;
                const float dc = d - c;
                const float de = d - e;
                const float maximum = std::max({ de, dc, std::min(b, f) });
                const float minimum = std::min({ de, dc, std::max(b, f) });
                diff = std::max({ diff, minimum, -maximum });

                float interpol;
                if (std::abs(c - e) > temporal_diff0)
                    interpol = ((coef_hf_f[0] * (prev2[x] + next2[x])
                                 - coef_hf_f[1] * (prev2Above2[x] + next2Above2[x] + prev2Below2[x] + next2Below2[x])
                                 + coef_hf_f[2] * (prev2Above4[x] + next2Above4[x] + prev2Below4[x] + next2Below4[x])) * 0.25f
                                + coef_lf_f[0] * (c + e) - coef_lf_f[1] * (curAbove3[x] + curBelow3[x]));
                else
                    interpol = coef_sp_f[0] * (c + e) - coef_sp_f[1] * (curAbove3[x] + curBelow3[x]);

                dst[x] = std::clamp(interpol, d - diff, d + diff);
            }
        }
    }
}

template<typename pixel_t>
static void filter(const VSFrameRef * prevFrame, const VSFrameRef * curFrame, const VSFrameRef * nextFrame, VSFrameRef * dstFrame,
                   const int field, const BwdifData * const VS_RESTRICT d, const VSAPI * vsapi) noexcept {
    for (int plane = 0; plane < d->vi.format->numPlanes; plane++) {
        const int width = vsapi->getFrameWidth(curFrame, plane);
        const int height = vsapi->getFrameHeight(curFrame, plane);
        const int stride = vsapi->getStride(curFrame, plane) / sizeof(pixel_t);
        const pixel_t * prev = reinterpret_cast<const pixel_t *>(vsapi->getReadPtr(prevFrame, plane));
        const pixel_t * cur = reinterpret_cast<const pixel_t *>(vsapi->getReadPtr(curFrame, plane));
        const pixel_t * next = reinterpret_cast<const pixel_t *>(vsapi->getReadPtr(nextFrame, plane));
        pixel_t * VS_RESTRICT dst = reinterpret_cast<pixel_t *>(vsapi->getWritePtr(dstFrame, plane));

        vs_bitblt(dst + stride * (1 - field),
                  vsapi->getStride(dstFrame, plane) * 2,
                  cur + stride * (1 - field),
                  vsapi->getStride(curFrame, plane) * 2,
                  width * sizeof(pixel_t),
                  height / 2);

        prev += stride * field;
        cur += stride * field;
        next += stride * field;
        dst += stride * field;

        const pixel_t * prev2 = field ? prev : cur;
        const pixel_t * next2 = field ? cur : next;

        for (int y = field; y < height; y += 2) {
            if ((y < 4) || (y + 5 > height)) {
                if ((y < 2) || (y + 3 > height))
                    d->filterEdgeWithoutSpat(prev2, prev, cur, next, next2, dst, width,
                                             y + 1 < height ? stride : -stride,
                                             y > 0 ? -stride : stride,
                                             stride * 2,
                                             d->edgeStep, d->peak);
                else
                    d->filterEdgeWithSpat(prev2, prev, cur, next, next2, dst, width,
                                          y + 1 < height ? stride : -stride,
                                          y > 0 ? -stride : stride,
                                          stride * 2,
                                          d->edgeStep, d->peak);
            } else {
                d->filterLine(prev2, prev, cur, next, next2, dst, width,
                              stride, stride * 2, stride * 3, stride * 4,
                              d->lineStep, d->peak);
            }

            prev2 += stride * 2;
            prev += stride * 2;
            cur += stride * 2;
            next += stride * 2;
            next2 += stride * 2;
            dst += stride * 2;
        }
    }
}

static void VS_CC bwdifInit(VSMap * in, VSMap * out, void ** instanceData, VSNode * node, VSCore * core, const VSAPI * vsapi) {
    BwdifData * d = static_cast<BwdifData *>(*instanceData);
    vsapi->setVideoInfo(&d->vi, 1, node);
}

static const VSFrameRef * VS_CC bwdifGetFrame(int n, int activationReason, void ** instanceData, void ** frameData, VSFrameContext * frameCtx, VSCore * core, const VSAPI * vsapi) {
    const BwdifData * d = static_cast<const BwdifData *>(*instanceData);

    if (activationReason == arInitial) {
        if (d->field > 1)
            n /= 2;

        if (n > 0)
            vsapi->requestFrameFilter(n - 1, d->node, frameCtx);
        vsapi->requestFrameFilter(n, d->node, frameCtx);
        if (n < d->viSaved->numFrames - 1)
            vsapi->requestFrameFilter(n + 1, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const int nSaved = n;
        int field = d->field;
        if (d->field > 1) {
            n /= 2;
            field -= 2;
        }

        const VSFrameRef * prev = vsapi->getFrameFilter(std::max(n - 1, 0), d->node, frameCtx);
        const VSFrameRef * cur = vsapi->getFrameFilter(n, d->node, frameCtx);
        const VSFrameRef * next = vsapi->getFrameFilter(std::min(n + 1, d->viSaved->numFrames - 1), d->node, frameCtx);
        VSFrameRef * dst = vsapi->newVideoFrame(d->vi.format, d->vi.width, d->vi.height, cur, core);

        int err;
        const int fieldBased = int64ToIntS(vsapi->propGetInt(vsapi->getFramePropsRO(cur), "_FieldBased", 0, &err));
        if (fieldBased == 1)
            field = 0;
        else if (fieldBased == 2)
            field = 1;

        if (d->field > 1) {
            if (nSaved & 1)
                field = (field == 0);
            else
                field = (field == 1);
        }

        if (d->vi.format->bytesPerSample == 1)
            filter<uint8_t>(prev, cur, next, dst, field, d, vsapi);
        else if (d->vi.format->bytesPerSample == 2)
            filter<uint16_t>(prev, cur, next, dst, field, d, vsapi);
        else
            filter<float>(prev, cur, next, dst, field, d, vsapi);

        VSMap * props = vsapi->getFramePropsRW(dst);
        vsapi->propSetInt(props, "_FieldBased", 0, paReplace);

        if (d->field > 1) {
            int errNum, errDen;
            int64_t durationNum = vsapi->propGetInt(props, "_DurationNum", 0, &errNum);
            int64_t durationDen = vsapi->propGetInt(props, "_DurationDen", 0, &errDen);
            if (!errNum && !errDen) {
                muldivRational(&durationNum, &durationDen, 1, 2);
                vsapi->propSetInt(props, "_DurationNum", durationNum, paReplace);
                vsapi->propSetInt(props, "_DurationDen", durationDen, paReplace);
            }
        }

        vsapi->freeFrame(prev);
        vsapi->freeFrame(cur);
        vsapi->freeFrame(next);
        return dst;
    }

    return nullptr;
}

static void VS_CC bwdifFree(void * instanceData, VSCore * core, const VSAPI * vsapi) {
    BwdifData * d = static_cast<BwdifData *>(instanceData);
    vsapi->freeNode(d->node);
    delete d;
}

static void VS_CC bwdifCreate(const VSMap * in, VSMap * out, void * userData, VSCore * core, const VSAPI * vsapi) {
    using namespace std::literals;

    std::unique_ptr<BwdifData> d = std::make_unique<BwdifData>();

    try {
        d->node = vsapi->propGetNode(in, "clip", 0, nullptr);
        d->vi = *vsapi->getVideoInfo(d->node);
        d->viSaved = vsapi->getVideoInfo(d->node);
        int err;

        if (!isConstantFormat(&d->vi) ||
            (d->vi.format->sampleType == stInteger && d->vi.format->bitsPerSample > 16) ||
            (d->vi.format->sampleType == stFloat && d->vi.format->bitsPerSample != 32))
            throw "only constant format 8-16 bit integer and 32 bit float input supported"sv;

        if (d->vi.height < 4)
            throw "height must be greater than or equal to 4"sv;

        d->field = int64ToIntS(vsapi->propGetInt(in, "field", 0, nullptr));
        const int opt = int64ToIntS(vsapi->propGetInt(in, "opt", 0, &err));

        if (d->field < 0 || d->field > 3)
            throw "field must be 0, 1, 2, or 3"sv;

        if (opt < 0 || opt > 4)
            throw "opt must be 0, 1, 2, 3, or 4"sv;

        {
            if (d->vi.format->bytesPerSample == 1) {
                d->filterEdgeWithSpat = filterEdge_c<uint8_t, true>;
                d->filterEdgeWithoutSpat = filterEdge_c<uint8_t, false>;
                d->filterLine = filterLine_c<uint8_t>;
            } else if (d->vi.format->bytesPerSample == 2) {
                d->filterEdgeWithSpat = filterEdge_c<uint16_t, true>;
                d->filterEdgeWithoutSpat = filterEdge_c<uint16_t, false>;
                d->filterLine = filterLine_c<uint16_t>;
            } else {
                d->filterEdgeWithSpat = filterEdge_c<float, true>;
                d->filterEdgeWithoutSpat = filterEdge_c<float, false>;
                d->filterLine = filterLine_c<float>;
            }

#ifdef BWDIF_X86
            const int iset = instrset_detect();
            if ((opt == 0 && iset >= 10) || opt == 4) {
                if (d->vi.format->bytesPerSample == 1) {
                    d->filterEdgeWithSpat = filterEdge_avx512<uint8_t, true>;
                    d->filterEdgeWithoutSpat = filterEdge_avx512<uint8_t, false>;
                    d->filterLine = filterLine_avx512<uint8_t>;
                    d->edgeStep = 32;
                } else if (d->vi.format->bytesPerSample == 2) {
                    d->filterEdgeWithSpat = filterEdge_avx512<uint16_t, true>;
                    d->filterEdgeWithoutSpat = filterEdge_avx512<uint16_t, false>;
                    d->filterLine = filterLine_avx512<uint16_t>;
                    d->edgeStep = 16;
                } else {
                    d->filterEdgeWithSpat = filterEdge_avx512<float, true>;
                    d->filterEdgeWithoutSpat = filterEdge_avx512<float, false>;
                    d->filterLine = filterLine_avx512<float>;
                    d->edgeStep = 16;
                }
                d->lineStep = 16;
            } else if ((opt == 0 && iset >= 8) || opt == 3) {
                if (d->vi.format->bytesPerSample == 1) {
                    d->filterEdgeWithSpat = filterEdge_avx2<uint8_t, true>;
                    d->filterEdgeWithoutSpat = filterEdge_avx2<uint8_t, false>;
                    d->filterLine = filterLine_avx2<uint8_t>;
                    d->edgeStep = 16;
                } else if (d->vi.format->bytesPerSample == 2) {
                    d->filterEdgeWithSpat = filterEdge_avx2<uint16_t, true>;
                    d->filterEdgeWithoutSpat = filterEdge_avx2<uint16_t, false>;
                    d->filterLine = filterLine_avx2<uint16_t>;
                    d->edgeStep = 8;
                } else {
                    d->filterEdgeWithSpat = filterEdge_avx2<float, true>;
                    d->filterEdgeWithoutSpat = filterEdge_avx2<float, false>;
                    d->filterLine = filterLine_avx2<float>;
                    d->edgeStep = 8;
                }
                d->lineStep = 8;
            } else if ((opt == 0 && iset >= 2) || opt == 2) {
                if (d->vi.format->bytesPerSample == 1) {
                    d->filterEdgeWithSpat = filterEdge_sse2<uint8_t, true>;
                    d->filterEdgeWithoutSpat = filterEdge_sse2<uint8_t, false>;
                    d->filterLine = filterLine_sse2<uint8_t>;
                    d->edgeStep = 8;
                } else if (d->vi.format->bytesPerSample == 2) {
                    d->filterEdgeWithSpat = filterEdge_sse2<uint16_t, true>;
                    d->filterEdgeWithoutSpat = filterEdge_sse2<uint16_t, false>;
                    d->filterLine = filterLine_sse2<uint16_t>;
                    d->edgeStep = 4;
                } else {
                    d->filterEdgeWithSpat = filterEdge_sse2<float, true>;
                    d->filterEdgeWithoutSpat = filterEdge_sse2<float, false>;
                    d->filterLine = filterLine_sse2<float>;
                    d->edgeStep = 4;
                }
                d->lineStep = 4;
            }
#endif
        }

        if (d->field > 1) {
            if (d->vi.numFrames > INT_MAX / 2)
                throw "resulting clip is too long"sv;
            d->vi.numFrames *= 2;

            muldivRational(&d->vi.fpsNum, &d->vi.fpsDen, 2, 1);
        }

        d->peak = (1 << d->vi.format->bitsPerSample) - 1;
    } catch (const std::string_view & error) {
        vsapi->setError(out, ("Bwdif: "s + error.data()).c_str());
        vsapi->freeNode(d->node);
        return;
    }

    vsapi->createFilter(in, out, "Bwdif", bwdifInit, bwdifGetFrame, bwdifFree, fmParallel, 0, d.release(), core);
}

//////////////////////////////////////////
// Init

VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin * plugin) {
    configFunc("com.holywu.bwdif", "bwdif", "BobWeaver Deinterlacing Filter", VAPOURSYNTH_API_VERSION, 1, plugin);
    registerFunc("Bwdif",
                 "clip:clip;"
                 "field:int;"
                 "opt:int:opt;",
                 bwdifCreate, nullptr, plugin);
}
