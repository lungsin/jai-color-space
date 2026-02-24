// sRGB color space.
// Note: each channel has value in the range of 0-1.
sRGB :: struct {
    r, g, b: float64;

    #overlay (r) components: [3]float64;
}

// CIEXYZ color space.
CIEXYZ :: struct {
    x, y, z: float64;

    #overlay (x) components: [3]float64;
}

// CIELAB color space.
CIELAB :: struct {
    l, a, b: float64;

    #overlay (l) components: [3]float64;
}

// Oklab color space.
Oklab :: struct {
    l, a, b: float64;

    #overlay (l) components: [3]float64;
}

// D65 reference point for white light using standard 2° degree observer, normalized to Y = 1.0
// https://en.wikipedia.org/wiki/Standard_illuminant#D65_values
D65 :: CIEXYZ.{x = 0.95047, y = 1.00000, z = 1.08883};

/* Color conversion */

/*
    Linear sRGB color conversion
*/

srgb_to_linear_srgb :: (c: sRGB) -> sRGB {
    linearize :: (v: float64) -> float64 {
        if v <= 0.04045 {
            return v / 12.92;
        }

        return pow((v + 0.055) / 1.055, 2.4);
    }
    return sRGB.{linearize(c.r), linearize(c.g), linearize(c.b)};
}

linear_srgb_to_srgb :: (c: sRGB) -> sRGB {
    delinearize :: (v: float64) -> float64 {
        if v <= 0.0031308 {
            return 12.92 * v;
        }

        return 1.055 * pow(v, 1.0 / 2.4) - 0.055;
    }
    return sRGB.{delinearize(c.r), delinearize(c.g), delinearize(c.b)};
}

/*
   CIEXYZ conversion
*/

linear_srgb_to_ciexyz :: (c: sRGB) -> CIEXYZ {
    return .{x = 0.41239079926595948 * c.r + 0.35758433938387796 * c.g + 0.18048078840183429 * c.b,
             y = 0.21263900587151036 * c.r + 0.71516867876775593 * c.g + 0.072192315360733715 * c.b,
             z = 0.019330818715591851 * c.r + 0.11919477979462599 * c.g + 0.95053215224966058 * c.b};
}

srgb_to_ciexyz :: (c: sRGB) -> CIEXYZ {
    lin_c := srgb_to_linear_srgb(c);
    return linear_srgb_to_ciexyz(lin_c);
}

ciexyz_to_linear_srgb :: (c: CIEXYZ) -> sRGB {
    return .{r = 3.2409699419045214 * c.x - 1.5373831775700935 * c.y - 0.49861076029300328 * c.z,
             g = -0.96924363628087983 * c.x + 1.8759675015077207 * c.y + 0.041555057407175613 * c.z,
             b = 0.055630079696993609 * c.x - 0.20397695888897657 * c.y + 1.0569715142428786 * c.z};
}

ciexyz_to_srgb :: (c: CIEXYZ) -> sRGB {
    lin_c := ciexyz_to_linear_srgb(c);
    return linear_srgb_to_srgb(lin_c);
}

/*
    CIELAB conversion
*/

// Convert CIEXYZ color space to CIELAB taking into account a given white reference point.
ciexyz_to_cielab_white_ref :: (c: CIEXYZ, wref: CIEXYZ) -> CIELAB {
    f :: inline (t: float64) -> float64 {
        if t > 216.0 / 24389.0 return cbrt(t);
        // KAPPA = 24389 / 27
        // (KAPPA * t + 16) / 116
        //  = (24389 / 27 / 116) * t + (4 / 29)
        //  = (841 / 108) * t + (4 / 29)
        return 841.0 / 108.0 * t + 4.0 / 29.0;
    }

    fx := f(c.x / wref.x);
    fy := f(c.y / wref.y);
    fz := f(c.z / wref.z);

    return .{
             l = 116.0 * fy - 16.0,
             a = 500.0 * (fx - fy),
             b = 200.0 * (fy - fz)};
}

// Convert CIELAB color space to CIEXYZ taking into account a given white reference point.
cielab_to_ciexyz_white_ref :: (c: CIELAB, wref: CIEXYZ) -> CIEXYZ {
    f_inv :: inline (t: float64) -> float64 {
        if t > 6.0 / 29.0 return cube(t);
        // 3 * (6 / 29)^2 * (t - 4 / 29) 
        // = 108 / 841 * t - 432 / 24389
        return 108.0 / 841.0 * t - 432.0 / 24389.0;
    }

    fy := (c.l + 16.0) / 116.0;
    fx := c.a / 500.0 + fy;
    fz := fy - c.b / 200.0;

    xr := f_inv(fx);
    zr := f_inv(fz);
    // special optimization for yr for more accuracy. Instead of using f_inv(yr):
    // if fy > 6 / 29 
    //      <=> (c.l + 16) / 116 > 6 / 29
    //      <=> c.l > 8
    // then cube(fy)
    // else 108 / 841 * (c.l + 16) / 116 - 432 / 24389
    //      = 27 * c.l / 24389
    yr := ifx c.l > 8.0 then cube(fy) else c.l * 27.0 / 24389.0;

    return .{wref.x * xr,
             wref.y * yr,
             wref.z * zr};
}

// Convert CIEXYZ to CIELAB using D65 as white reference point.
ciexyz_to_cielab :: (c: CIEXYZ) -> CIELAB {
    return ciexyz_to_cielab_white_ref(c, D65);
}

// Convert CIELAB to CIEXYZ using D65 as white reference point.
cielab_to_ciexyz :: (c: CIELAB) -> CIEXYZ {
    return cielab_to_ciexyz_white_ref(c, D65);
}

// Convert sRGB to CIELAB
srgb_to_cielab :: (c: sRGB) -> CIELAB {
    return ciexyz_to_cielab(srgb_to_ciexyz(c));
}

// Convert CIELAB to sRGB
cielab_to_srgb :: (c: CIELAB) -> sRGB {
    return ciexyz_to_srgb(cielab_to_ciexyz(c));
}

/*
    Oklab conversion
*/

linear_srgb_to_oklab :: (c: sRGB) -> Oklab {
    l := cbrt(0.4122214708 * c.r + 0.5363325363 * c.g + 0.0514459929 * c.b);
    m := cbrt(0.2119034982 * c.r + 0.6806995451 * c.g + 0.1073969566 * c.b);
    s := cbrt(0.0883024619 * c.r + 0.2817188376 * c.g + 0.6299787005 * c.b);

    return .{
             l = 0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
             a = 1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
             b = 0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s};
}

oklab_to_linear_srgb :: (c: Oklab) -> sRGB {
    l_ := c.l + 0.3963377774 * c.a + 0.2158037573 * c.b;
    m_ := c.l - 0.1055613458 * c.a - 0.0638541728 * c.b;
    s_ := c.l - 0.0894841775 * c.a - 1.2914855480 * c.b;

    l := cube(l_);
    m := cube(m_);
    s := cube(s_);

    return .{
             r = +4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
             g = -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
             b = -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s};
}

srgb_to_oklab :: (c: sRGB) -> Oklab {
    return linear_srgb_to_oklab(srgb_to_linear_srgb(c));
}

oklab_to_srgb :: (c: Oklab) -> sRGB {
    return linear_srgb_to_srgb(oklab_to_linear_srgb(c));
}

#scope_file

#import "Basic";
#import "Math";
#import "jai-cbrt";

cube :: inline (x: float64) -> float64 {
    return x * x * x;
}
