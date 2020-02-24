Description
===========

Motion adaptive deinterlacing based on yadif with the use of w3fdif and cubic interpolation algorithms.

Ported from FFmpeg's libavfilter https://www.ffmpeg.org/ffmpeg-filters.html#bwdif


Usage
=====

    bwdif.Bwdif(clip clip, int field)

* clip: Clip to process. Any planar format with integer sample type of 8-16 bit depth is supported.

* field: Controls the mode of operation (double vs same rate) and which field is kept.
  * 0 = same rate, keep bottom field
  * 1 = same rate, keep top field
  * 2 = double rate (alternates each frame), starts with bottom
  * 3 = double rate (alternates each frame), starts with top


Compilation
===========

```
meson build
ninja -C build
ninja -C build install
```
