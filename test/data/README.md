TIFF files written by ImageJ 1.54g (batch macro `check3.ijm` in
`Phase/plots/SampledDomains/ImageJ_check` of the private Phase repository): 8 × 6 pixels,
32-bit, value `x + 10y` of the 0-based pixel `(x, y)` of ImageJ (y down), so that every
flip or transposition is visible; pixel 5.5 × 4.0 micron, origin (4, 3) in pixels.
`cal_inverty.tif` with "Invert Y coordinates" (y up), `cal_plain.tif` without.

ImageJ measured the centres of two pixels:

| file | pixel (3, 1), value 13 | pixel (6, 4), value 46 |
|---|---|---|
| cal_inverty.tif | (-2.75, 6) | (13.75, -6) |
| cal_plain.tif | (-2.75, -6) | (13.75, 6) |

In Julia these are `img[2, 4]` and `img[5, 7]`. ImageJ writes the pixel width as the rational
181818/1000000 pixels per micron, i.e. 5.5000055 micron.
