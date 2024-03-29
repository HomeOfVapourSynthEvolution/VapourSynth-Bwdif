# Bwdif
Motion adaptive deinterlacing based on yadif with the use of w3fdif and cubic interpolation algorithms.

Ported from FFmpeg's libavfilter https://www.ffmpeg.org/ffmpeg-filters.html#bwdif


## Usage
    bwdif.Bwdif(vnode clip, int field[, vnode edeint=None, int opt=0])

- clip: Clip to process. Any format with either integer sample type of 8-16 bit depth or float sample type of 32 bit depth is supported.

- field: Controls the mode of operation (double vs same rate) and which field is kept.
  - 0 = same rate, keep bottom field
  - 1 = same rate, keep top field
  - 2 = double rate (alternates each frame), starts with bottom
  - 3 = double rate (alternates each frame), starts with top

- edeint: Allows the specification of an external clip from which to take spatial predictions instead of having Bwdif use cubic interpolation. This clip must be the same width, height, and colorspace as the input clip. If using same rate output, this clip should have the same number of frames as the input. If using double rate output, this clip should have twice as many frames as the input.

- opt: Sets which cpu optimizations to use.
  - 0 = auto detect
  - 1 = use c
  - 2 = use sse2
  - 3 = use avx2
  - 4 = use avx512


## Compilation
```
meson build
ninja -C build
ninja -C build install
```
