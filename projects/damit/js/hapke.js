#!/usr/bin/env node

const { PI, sqrt, sin, cos, tan, log, exp, acos, pow, min, max } = Math;
const deg = PI / 180.0;

// Global parameters
let B0, minh, ming, bartheta;
let B, P, tanbartheta;
let mu_i_, mu_e_;

/**
 * Hapke law.
 * Miroslav Broz (miroslav.broz@email.cz), Jan 14th 2023
 *
 * References:
 *  - Spjuth (2009)
 *  - Hapke (1984)
 *  - Kuzminykh (2021)
 */
function f_hapke(f_L, mu_i, mu_e, alpha) {
  let f;
  if (mu_i > 0.0 && mu_e > 0.0) {
    const A_w = 4.0 * PI * f_L;
    const Stmp = Sr(mu_i, mu_e, alpha);
    f = f_L / (mu_i_ + mu_e_) * ((1.0 + B) * P + H(mu_i_, A_w) * H(mu_e_, A_w) - 1.0) * Stmp;
    f *= mu_e_ / mu_e;
  } else {
    f = 0.0;
  }
  return f;
}

function init_hapke(alpha) {
  // Eq. (2.26), (2.13)
  B = B0 / (1.0 + 1.0 / minh * tan(alpha / 2.0));
  P = (1.0 - ming * ming) / pow(1.0 + 2.0 * ming * cos(alpha) + ming * ming, 1.5);
  tanbartheta = tan(bartheta);
}

/**
 * Chandrasekhar function H(mu, A_w)
 * Eq. (2.17)–(2.19)
 */
function H(mu, A_w) {
  const gamma = sqrt(1.0 - A_w);  
  const r0    = (1.0 - gamma) / (1.0 + gamma);
  return pow(
    1.0 - A_w * mu * (r0 + (1.0 - 2.0 * r0 * mu) / 2.0 * log((1.0 + mu) / mu)),
    -1.0
  );
}

/**
 * Surface roughness function Sr(mu_i, mu_e, alpha)
 * Eqs. (2.36)–(2.41), (12.48)–(12.51)
 */
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
  const f  = exp(-2.0 * tan(psi / 2.0));

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

    Sr_val = (mu_i / eta_i) * (mu_e_ / eta_e) * (xi / (1.0 - f + f * xi * mu_i / eta_i));
  } else {
    const K1 = cospsi * E2i + sinpsihalfsq * E2e;
    const K2 = 2.0 - E1i - (psi / PI) * E1e;
    const K3 = E2i - sinpsihalfsq * E2e;

    mu_i_ = xi * (cosi + sini * tanbartheta * K3 / K2);
    mu_e_ = xi * (cose + sine * tanbartheta * K1 / K2);

    Sr_val = (mu_i / eta_i) * (mu_e_ / eta_e) * (xi / (1.0 - f + f * xi * mu_e / eta_e));
  }

  return Sr_val;
}

// Example usage, matching your Python `main()`
function main() {
  const Aw = 0.23;
  B0 = 1.32;
  minh = 0.20;
  ming = -0.35;
  bartheta = 10.0 * deg;

  const mu_i = 0.2;
  const mu_e = 0.4;
  const f_L  = Aw / (4.0 * PI);
  const alpha = 30.0 * deg;

  init_hapke(alpha);
  const f = f_hapke(f_L, mu_i, mu_e, alpha);
  console.log("f =", f);
}

if (require.main === module) {
  main();
}
