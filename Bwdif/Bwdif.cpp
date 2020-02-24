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

#include <cstdlib>

#include <algorithm>
#include <memory>
#include <string>

#include <VapourSynth.h>
#include <VSHelper.h>

/*
 * Filter coefficients coef_lf and coef_hf taken from BBC PH-2071 (Weston 3 Field Deinterlacer).
 * Used when there is spatial and temporal interpolation.
 * Filter coefficients coef_sp are used when there is spatial interpolation only.
 * Adjusted for matching visual sharpness impression of spatial and temporal interpolation.
 */
static constexpr uint16_t coef_lf[2] = { 4309, 213 };
static constexpr uint16_t coef_hf[3] = { 5570, 3801, 1016 };
static constexpr uint16_t coef_sp[2] = { 5077, 981 };

struct BwdifData {
    VSNodeRef * node;
    VSVideoInfo vi;
    const VSVideoInfo * viSaved;
    int field;
    int peak;
};

template<typename T, bool spat>
static inline void filterEdge(const T * prev2, const T * prev, const T * cur, const T * next, const T * next2, T * VS_RESTRICT dst, const int width,
                              const int positiveStride, const int negativeStride, const int positiveStride2, const int negativeStride2,
                              const int peak) noexcept {
    const T * prev2Above2 = prev2 + negativeStride2;
    const T * prev2Below2 = prev2 + positiveStride2;
    const T * prevAbove = prev + negativeStride;
    const T * prevBelow = prev + positiveStride;
    const T * curAbove = cur + negativeStride;
    const T * curBelow = cur + positiveStride;
    const T * nextAbove = next + negativeStride;
    const T * nextBelow = next + positiveStride;
    const T * next2Above2 = next2 + negativeStride2;
    const T * next2Below2 = next2 + positiveStride2;

    for (int x = 0; x < width; x++) {
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
            if (spat) {
                const int b = ((prev2Above2[x] + next2Above2[x]) >> 1) - c;
                const int f = ((prev2Below2[x] + next2Below2[x]) >> 1) - e;
                const int dc = d - c;
                const int de = d - e;
                const int max = std::max({ de, dc, std::min(b, f) });
                const int min = std::min({ de, dc, std::max(b, f) });
                diff = std::max({ diff, min, -max });
            }

            const int interpol = std::min(std::max((c + e) >> 1, d - diff), d + diff);
            dst[x] = std::min(std::max(interpol, 0), peak);
        }
    }
}

template<typename T>
static inline void filterLine(const T * prev2, const T * prev, const T * cur, const T * next, const T * next2, T * VS_RESTRICT dst, const int width,
                              const int positiveStride, const int negativeStride, const int positiveStride2, const int negativeStride2,
                              const int positiveStride3, const int negativeStride3, const int positiveStride4, const int negativeStride4,
                              const int peak) noexcept {
    const T * prev2Above4 = prev2 + negativeStride4;
    const T * prev2Below4 = prev2 + positiveStride4;
    const T * prev2Above2 = prev2 + negativeStride2;
    const T * prev2Below2 = prev2 + positiveStride2;
    const T * prevAbove = prev + negativeStride;
    const T * prevBelow = prev + positiveStride;
    const T * curAbove3 = cur + negativeStride3;
    const T * curAbove = cur + negativeStride;
    const T * curBelow = cur + positiveStride;
    const T * curBelow3 = cur + positiveStride3;
    const T * nextAbove = next + negativeStride;
    const T * nextBelow = next + positiveStride;
    const T * next2Above2 = next2 + negativeStride2;
    const T * next2Below2 = next2 + positiveStride2;
    const T * next2Above4 = next2 + negativeStride4;
    const T * next2Below4 = next2 + positiveStride4;

    for (int x = 0; x < width; x++) {
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
            const int max = std::max({ de, dc, std::min(b, f) });
            const int min = std::min({ de, dc, std::max(b, f) });
            diff = std::max({ diff, min, -max });

            int interpol;
            if (std::abs(c - e) > temporal_diff0)
                interpol = (((coef_hf[0] * (prev2[x] + next2[x])
                              - coef_hf[1] * (prev2Above2[x] + next2Above2[x] + prev2Below2[x] + next2Below2[x])
                              + coef_hf[2] * (prev2Above4[x] + next2Above4[x] + prev2Below4[x] + next2Below4[x])) >> 2)
                            + coef_lf[0] * (c + e) - coef_lf[1] * (curAbove3[x] + curBelow3[x])) >> 13;
            else
                interpol = (coef_sp[0] * (c + e) - coef_sp[1] * (curAbove3[x] + curBelow3[x])) >> 13;

            if (interpol > d + diff)
                interpol = d + diff;
            else if (interpol < d - diff)
                interpol = d - diff;

            dst[x] = std::min(std::max(interpol, 0), peak);
        }
    }
}

template<typename T>
static void filter(const VSFrameRef * prevFrame, const VSFrameRef * curFrame, const VSFrameRef * nextFrame, VSFrameRef * dstFrame,
                   const int field, const BwdifData * const VS_RESTRICT d, const VSAPI * vsapi) noexcept {
    for (int plane = 0; plane < d->vi.format->numPlanes; plane++) {
        const int width = vsapi->getFrameWidth(curFrame, plane);
        const int height = vsapi->getFrameHeight(curFrame, plane);
        const int stride = vsapi->getStride(curFrame, plane) / sizeof(T);
        const T * prev = reinterpret_cast<const T *>(vsapi->getReadPtr(prevFrame, plane));
        const T * cur = reinterpret_cast<const T *>(vsapi->getReadPtr(curFrame, plane));
        const T * next = reinterpret_cast<const T *>(vsapi->getReadPtr(nextFrame, plane));
        T * VS_RESTRICT dst = reinterpret_cast<T *>(vsapi->getWritePtr(dstFrame, plane));

        vs_bitblt(dst + stride * (1 - field), vsapi->getStride(dstFrame, plane) * 2, cur + stride * (1 - field), vsapi->getStride(curFrame, plane) * 2, width * sizeof(T), height / 2);

        prev += stride * field;
        cur += stride * field;
        next += stride * field;
        dst += stride * field;

        const T * prev2 = field ? prev : cur;
        const T * next2 = field ? cur : next;

        for (int y = field; y < height; y += 2) {
            if ((y < 4) || (y + 5 > height)) {
                if ((y < 2) || (y + 3 > height))
                    filterEdge<T, false>(prev2, prev, cur, next, next2, dst, width,
                                         y + 1 < height ? stride : -stride,
                                         y > 0 ? -stride : stride,
                                         stride * 2, stride * -2,
                                         d->peak);
                else
                    filterEdge<T, true>(prev2, prev, cur, next, next2, dst, width,
                                        y + 1 < height ? stride : -stride,
                                        y > 0 ? -stride : stride,
                                        stride * 2, stride * -2,
                                        d->peak);
            } else {
                filterLine(prev2, prev, cur, next, next2, dst, width,
                           stride, -stride, stride * 2, stride * -2,
                           stride * 3, stride * -3, stride * 4, stride * -4,
                           d->peak);
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
        else
            filter<uint16_t>(prev, cur, next, dst, field, d, vsapi);

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
    std::unique_ptr<BwdifData> d = std::make_unique<BwdifData>();

    d->node = vsapi->propGetNode(in, "clip", 0, nullptr);
    d->vi = *vsapi->getVideoInfo(d->node);
    d->viSaved = vsapi->getVideoInfo(d->node);

    try {
        if (!isConstantFormat(&d->vi) || d->vi.format->sampleType != stInteger || d->vi.format->bitsPerSample > 16)
            throw std::string{ "only constant format 8-16 bit integer input supported" };

        if (d->vi.height < 4)
            throw std::string{ "height must be greater than or equal to 4" };

        d->field = int64ToIntS(vsapi->propGetInt(in, "field", 0, nullptr));

        if (d->field < 0 || d->field > 3)
            throw std::string{ "field must be 0, 1, 2, or 3" };

        if (d->field > 1) {
            if (d->vi.numFrames > INT_MAX / 2)
                throw std::string{ "resulting clip is too long" };
            d->vi.numFrames *= 2;

            muldivRational(&d->vi.fpsNum, &d->vi.fpsDen, 2, 1);
        }

        d->peak = (1 << d->vi.format->bitsPerSample) - 1;
    } catch (const std::string & error) {
        vsapi->setError(out, ("Bwdif: " + error).c_str());
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
                 "field:int;",
                 bwdifCreate, nullptr, plugin);
}
