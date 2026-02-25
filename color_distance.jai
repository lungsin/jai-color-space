// Euclidean distance
distance_euclidean :: (a: $T / .[sRGB, Oklab], b: T) -> float64 {
    return distance_euclidean(a.components, b.components);
}

/* CIE76 */

// Just Noticable Difference constant for CIE76.
JND_CIE76 :: 2.3;

// CIE76 color distance algorithm. Basically euclidean distance.
// Check out JND_CIE76 for the Just Noticable Difference constant.
distance_CIE76 :: (a: CIELAB, b: CIELAB) -> float64 {
    return distance_euclidean(a.components, b.components);
}

/* CIE94 */

// Parameters for CIE94 color distance algorithm.
CIE94_Params :: struct {
    kl: float64;
    kc := 1.0; // default
    kh := 1.0; // default
    k1: float64;
    k2: float64;
}

// Parameters for common domains for CIE94.
CIE94_GRAPHIC_ARTS_PARAMS :: CIE94_Params.{kl = 1.0, k1 = 0.045, k2 = 0.015};
CIE94_TEXTILE_PARAMS :: CIE94_Params.{kl = 2.0, k1 = 0.048, k2 = 0.014};

// CIE94 color distance algorithm.
// There are some predefined params that can be used depending on the application domain:
//  - CIE94_GRAPHIC_ARTS_PARAMS
//  - CIE94_TEXTILE_PARAMS
// Your app probably only need one version of this function, I recommend to create an alias that uses the specific
// params.
distance_CIE94 :: (a: CIELAB, b: CIELAB, p: CIE94_Params) -> float64 {
    delta_l := a.l - b.l;

    c1 := sqrt(sq(a.a) + sq(a.b));
    c2 := sqrt(sq(b.a) + sq(b.b));
    delta_c := c1 - c2;

    delta_a := a.a - b.a;
    delta_b := a.b - b.b;
    delta_h_sq := sq(delta_a) + sq(delta_b) - sq(delta_c);

    sl := 1.0;
    sc := 1.0 + p.k1 * c1;
    sh := 1.0 + p.k2 * c1;

    l_term := sq(delta_l / (p.kl * sl));
    c_term := sq(delta_c / (p.kc * sc));
    h_term := delta_h_sq / sq(p.kh * sh);

    return sqrt(l_term + c_term + h_term);
}

/* CIEDE2000 */

// Parameters for CIEDE2000 color distance algorithm
CIEDE2000_Params :: struct {
    kl := 1.0;
    kc := 1.0;
    kh := 1.0;
}

// CIEDE2000 color distance algorithm.
distance_CIEDE2000 :: (a: CIELAB, b: CIELAB, p: CIEDE2000_Params = .{}) -> float64 {
    POW25_7 :: 6103515625.0; // pow(25, 7)
    DEG_360_IN_RAD :: PI64 * 2;
    DEG_180_IN_RAD :: PI64;

    pow7 :: inline (x: float64) -> float64 {
        x_2 := sq(x);
        return x_2 * x_2 * x_2 * x;
    }

    atan2_non_negative :: inline (y: float64, x: float64) -> float64 {
        h := atan2(y, x);
        if h < 0.0 h += PI64 * 2;
        return h;
    }

    deg_to_rad :: inline (deg: float64) -> float64 {
        DEG_TO_RAD :: PI64 / 180.0;
        return deg * DEG_TO_RAD;
    }

    c1 := sqrt(sq(a.a) + sq(a.b));
    c2 := sqrt(sq(b.a) + sq(b.b));
    c_avg := (c1 + c2) * 0.5;

    c_avg_7 := pow7(c_avg);
    g := 0.5 * (1.0 - sqrt(c_avg_7 / (c_avg_7 + POW25_7)));

    a1_p := (1.0 + g) * a.a;
    a2_p := (1.0 + g) * b.a;
    c1_p := sqrt(sq(a1_p) + sq(a.b));
    c2_p := sqrt(sq(a2_p) + sq(b.b));

    h1_p := ifx c1_p == 0.0 then 0.0 else atan2_non_negative(a.b, a1_p);
    h2_p := ifx c2_p == 0.0 then 0.0 else atan2_non_negative(b.b, a2_p);

    delta_l_p := b.l - a.l;
    delta_c_p := c2_p - c1_p;

    delta_h_p: float64 = 0.0;
    c_p_product := c1_p * c2_p;
    if c_p_product != 0.0 {
        delta_h_p = h2_p - h1_p;
        if delta_h_p > DEG_180_IN_RAD delta_h_p -= DEG_360_IN_RAD;
        if delta_h_p < -DEG_180_IN_RAD delta_h_p += DEG_360_IN_RAD;
    }

    delta_h_term := 2.0 * sqrt(c_p_product) * sin(delta_h_p * 0.5);

    l_avg_p := (a.l + b.l) * 0.5;
    c_avg_p := (c1_p + c2_p) * 0.5;

    h_p_sum := h1_p + h2_p;
    h_avg_p := ifx c_p_product == 0.0 {
        h_p_sum;
    } else {
        ifx abs(h1_p - h2_p) <= DEG_180_IN_RAD then
            h_p_sum * 0.5
        else {
            delta := ifx h_p_sum < DEG_360_IN_RAD then
                         +DEG_360_IN_RAD
                     else
                         -DEG_360_IN_RAD;
            (h_p_sum + delta) * 0.5;
        }
    };

    t := 1.0 -
         0.17 * cos(h_avg_p - deg_to_rad(30.0)) +
         0.24 * cos(2.0 * h_avg_p) +
         0.32 * cos(3.0 * h_avg_p + deg_to_rad(6.0)) -
         0.20 * cos(4.0 * h_avg_p - deg_to_rad(63.0));

    delta_theta := deg_to_rad(30.0) * exp(-sq((h_avg_p - deg_to_rad(275.0)) / deg_to_rad(25.0)));

    c_avg_p_7 := pow7(c_avg_p);
    rc := 2.0 * sqrt(c_avg_p_7 / (c_avg_p_7 + POW25_7));
    rt := -sin(2.0 * delta_theta) * rc;

    l_shift_squared := sq(l_avg_p - 50.0);
    sl := 1.0 + (0.015 * l_shift_squared) / sqrt(20.0 + l_shift_squared);
    sc := 1.0 + 0.045 * c_avg_p;
    sh := 1.0 + 0.015 * c_avg_p * t;

    l_term := delta_l_p / (sl * p.kl);
    c_term := delta_c_p / (sc * p.kc);
    h_term := delta_h_term / (sh * p.kh);

    return sqrt(sq(l_term) + sq(c_term) + sq(h_term) + rt * c_term * h_term);
}

#scope_file
#import "Math";

sq :: inline (x: float64) -> float64 {
    return x * x;
}

distance_euclidean :: inline (a: [3]float64, b: [3]float64) -> float64 {
    return sqrt(distance_squared(a, b));
}

distance_squared :: inline (a: [3]float64, b: [3]float64) -> float64 {
    return sq(a[0] - b[0]) + sq(a[1] - b[1]) + sq(a[2] - b[2]);
}
