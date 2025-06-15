// hapke.js

// pull Math functions into local scope
const { PI, sqrt, sin, cos, tan, log, exp, acos, pow, min, max } = Math;

// degrees → radians factor
const deg = PI / 180.0;

// global Hapke parameters (must be set before calling init_hapke)
let B0, minh, ming, bartheta;
let B, P, tanbartheta;
let mu_i_, mu_e_;

/**
 * Hapke bidirectional scattering function
 * @param {number} f_L  – Lambertian term
 * @param {number} mu_i – cos(incidence angle)
 * @param {number} mu_e – cos(emergence angle)
 * @param {number} alpha – phase angle (radians)
 * @returns {number}
 */
function f_hapke(f_L, mu_i, mu_e, alpha) {
  let f;
  if (mu_i > 0.0 && mu_e > 0.0) {
    const A_w = 4.0 * PI * f_L;
    const Stmp = Sr(mu_i, mu_e, alpha);
    f = f_L / (mu_i_ + mu_e_) *
        ((1.0 + B) * P + H(mu_i_, A_w) * H(mu_e_, A_w) - 1.0) *
        Stmp;
    f *= mu_e_ / mu_e;
  } else {
    f = 0.0;
  }
  return f;
}

/**
 * Initialize the Hapke model for a given phase angle
 * @param {number} alpha – phase angle (radians)
 */
function init_hapke(alpha) {
  B = B0 / (1.0 + 1.0 / minh * tan(alpha / 2.0));  
  P = (1.0 - ming * ming) /
      pow(1.0 + 2.0 * ming * cos(alpha) + ming * ming, 1.5);
  tanbartheta = tan(bartheta);
}

/** Chandrasekhar H-function */
function H(mu, A_w) {
  const gamma = sqrt(1.0 - A_w);
  const r0 = (1.0 - gamma) / (1.0 + gamma);
  return pow(
    1.0 - A_w * mu *
      ( r0 + (1.0 - 2.0 * r0 * mu) / 2.0 * log((1.0 + mu) / mu) ),
    -1.0
  );
}

/** Surface‐roughness correction */
function Sr(mu_i, mu_e, alpha) {
  const cosi = mu_i, cose = mu_e;
  const sini = sqrt(1.0 - cosi * cosi);
  const sine = sqrt(1.0 - cose * cose);
  const tani = sini / cosi;
  const tane = sine / cose;

  const tmp = sini * sine;
  let cospsi = 0.0;
  if (tmp !== 0.0) {
    cospsi = (cos(alpha) - cosi * cose) / tmp;
  }
  const psi = acos(min(max(cospsi, -1.0), 1.0));
  const sinpsihalfsq = sin(psi / 2.0) ** 2;

  const xi = 1.0 / sqrt(1.0 + PI * tanbartheta * tanbartheta);
  const f_val = exp(-2.0 * tan(psi / 2.0));

  const E1i = exp(-2.0 / (PI * tanbartheta * tani));
  const E1e = exp(-2.0 / (PI * tanbartheta * tane));
  const E2i = exp(-1.0 / (PI * (tanbartheta * tani) ** 2));
  const E2e = exp(-1.0 / (PI * (tanbartheta * tane) ** 2));

  const eta_i = xi * (cosi + sini * tanbartheta * E2i / (2.0 - E1i));
  const eta_e = xi * (cose + sine * tanbartheta * E2e / (2.0 - E1e));

  let Sr_val;
  if (sini <= sine) {
    const K1 = cospsi * E2e + sinpsihalfsq * E2i;
    const K2 = 2.0 - E1e - (psi / PI) * E1i;
    const K3 = E2e - sinpsihalfsq * E2i;
    mu_i_ = xi * (cosi + sini * tanbartheta * K1 / K2);
    mu_e_ = xi * (cose + sine * tanbartheta * K3 / K2);
    Sr_val = (mu_i / eta_i) *
             (mu_e_ / eta_e) *
             (xi / (1.0 - f_val + f_val * xi * mu_i / eta_i));
  } else {
    const K1 = cospsi * E2i + sinpsihalfsq * E2e;
    const K2 = 2.0 - E1i - (psi / PI) * E1e;
    const K3 = E2i - sinpsihalfsq * E2e;
    mu_i_ = xi * (cosi + sini * tanbartheta * K3 / K2);
    mu_e_ = xi * (cose + sine * tanbartheta * K1 / K2);
    Sr_val = (mu_i / eta_i) *
             (mu_e_ / eta_e) *
             (xi / (1.0 - f_val + f_val * xi * mu_e / eta_e));
  }

  return Sr_val;
}

// expose for HTML usage
window.deg        = deg;
window.init_hapke = init_hapke;
window.f_hapke    = f_hapke;
