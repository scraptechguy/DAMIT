function f_lambert(f_L, mu_i, mu_e, alpha) {
    return f_L;
}

function rotate_matrix_inv(rotations) {
	//[rot_z, rot_y, rot_z, rot_x]
	var model_m = Matrix.RotationX(rotations[3])
		.multiply(Matrix.RotationZ(rotations[2]))
		.multiply(Matrix.RotationY(rotations[1]))
		.multiply(Matrix.RotationZ(rotations[0]));
    model_m = model_m.inverse();

	return model_m;
};

function rotate_vector_inv(vector, model_m) {
    var new_vector = [];
	new_vector = model_m.multiply($V(vector));

    return new_vector.elements;
}

/**
 * Intersection of ray A + x·B with triangle t in 3D.
 *
 * @param {number[]} A  – [Ax,Ay,Az]
 * @param {number[]} B  – [Bx,By,Bz]
 * @param {number[][]} t – [ [x0,y0,z0], [x1,y1,z1], [x2,y2,z2] ]
 * @returns {{ C: number[], hasSolution: boolean }}
 */
function intersectABT(A, B, t) {
    const EPS = 0.0;

    // helper: vector subtraction u - v
    function sub(u, v) {
        return [ u[0]-v[0], u[1]-v[1], u[2]-v[2] ];
    }

    // helper: cross product u × v
    function cross(u, v) {
        return [
        u[1]*v[2] - u[2]*v[1],
        u[2]*v[0] - u[0]*v[2],
        u[0]*v[1] - u[1]*v[0]
        ];
    }

    // helper: dot product u ⋅ v
    function dot(u, v) {
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
    }

    // 1) compute the “area” terms for each edge
    const area = [0, 0, 0];
    for (let i = 0; i < 3; i++) {
        // cross of (t[i] – A) and (t[i+1] – A)
        const Ccross = cross(
        sub(t[i % 3], A),
        sub(t[(i + 1) % 3], A)
        );
        area[i] = 0.5 * dot(B, Ccross);
    }

    // 2) check if all areas have the same sign (≥0 or ≤0)
    let hasSolution;
    if (
        area[0] >= -EPS && area[1] >= -EPS && area[2] >= -EPS
        || area[0] <= +EPS && area[1] <= +EPS && area[2] <= +EPS
    ) {
        hasSolution = true;
    } else {
        hasSolution = false;
    }

    // 3) if there is a solution, compute the interpolation point C
    let C = [0, 0, 0];
    if (hasSolution) {
        // normalize barycentrics: area /= sum(area)
        const sumA = area[0] + area[1] + area[2];
        const w = area.map(v => v / sumA);

        // C = w0*t0 + w1*t1 + w2*t2
        for (let i = 0; i < 3; i++) {
        C[0] += w[i] * t[i][0];
        C[1] += w[i] * t[i][1];
        C[2] += w[i] * t[i][2];
        }
    }

    return { C, hasSolution };
}

/**
 * Build non-illumination / non-visibility masks.
 *
 * @param {number[]} mu_i
 * @param {number[]} mu_e
 * @returns {{ nu_i: number[], nu_e: number[] }}
 */
function non(mu_i, mu_e) {
    const n = mu_i.length;
    const nu_i = new Array(n).fill(1);
    const nu_e = new Array(n).fill(1);

    for (let i = 0; i < n; i++) {
        if (mu_i[i] === 0.0 || mu_e[i] === 0.0) {
        nu_i[i] = 0.0;
        nu_e[i] = 0.0;
        }
    }

    return { nu_i, nu_e };
}

/**
 * Shadowing mask for a non‐convex mesh.
 *
 * @param {number[][]} faces    – array of [i0,i1,i2] vertex‐indices
 * @param {number[][]} nodes    – array of [x,y,z] vertex positions
 * @param {number[][]} normals  – array of [nx,ny,nz] triangle normals
 * @param {number[][]} centres  – array of [cx,cy,cz] triangle centroids
 * @param {number[]}   s        – direction vector for the “ray” (length 3)
 * @param {number[]}   nu_i     – per‐triangle visibility flags (1 = visible)
 * @returns {number[]}          – updated nu_i with occluded ones set to 0
 */
function nu(faces, nodes, normals, centres, s, nu_i) {
    // temporary 3×3 triangle vertices buffer
    const t = [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
    ];

    // for each triangle i
    for (let i = 0; i < faces.length; i++) {
        // only test those still marked visible
        if (nu_i[i] > 0.0) {
        const A = centres[i]; // ray origin
        const B = s;          // ray direction

        // test against every other triangle j
        let j = 0;
        while (j < faces.length && nu_i[i] > 0.0) {
            if (j !== i && nu_i[j] > 0.0) {
            // vector from triangle‐i center to triangle‐j center
            const Csub = [
                centres[j][0] - A[0],
                centres[j][1] - A[1],
                centres[j][2] - A[2],
            ];
            // dot with triangle‐i normal
            const tmp =
                Csub[0] * normals[i][0] +
                Csub[1] * normals[i][1] +
                Csub[2] * normals[i][2];

            if (tmp > 0.0) {
                // build the 3×3 vertex matrix for triangle j
                for (let k = 0; k < 3; k++) {
                const vidx = faces[j][k];
                t[k][0] = nodes[vidx][0];
                t[k][1] = nodes[vidx][1];
                t[k][2] = nodes[vidx][2];
                }

                // test intersection of ray A + B·x with triangle j
                const { C: _, hasSolution } = intersectABT(A, B, t);
                if (hasSolution) {
                // if it hits, triangle i is shadowed
                nu_i[i] = 0.0;
                }
            }
            }
            j++;
        }
        }
    }

    return nu_i;
}

/**
 * Compute incidence/emergence cosines and initialize masks.
 *
 * @param {number[]} [s=[1,0,0]] – sun‐direction vector
 * @param {number[]} [o=[0,0,1]] – observer‐direction vector
 */
function getCosinesAndFluxes(s = [1, 0, 0], o = [0, 0, 1], surfaces) {
    // helper: dot product of two 3‐vectors
    // const dot = (u, v) => u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
    // helper: dot product u ⋅ v
    function dot(u, v) {
        let tmp = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
        return tmp;
    }

    // normalise vectors 
    let absS = Math.sqrt(dot(s, s));
    s = [s[0] / absS, s[1] / absS, s[2] / absS];
    let absO = Math.sqrt(dot(o, o));
    o = [o[0] / absO, o[1] / absO, o[2] / absO];

    // phase angle α between s and o
    let alpha = Math.acos(dot(s, o));
    console.log("alpha: " + alpha)

    // compute raw cosines and clamp negatives to zero
    let mu_i = face_normals.map(normals => Math.max(dot(s, normals.elements), 0));
    let mu_e = face_normals.map(normals => Math.max(dot(o, normals.elements), 0));
    console.log("mu_i: " + mu_i)
    console.log("mu_e: " + mu_e)

    // initialize shadow/visibility masks to zero arrays
    const N = faces.length;
    let nu_i = new Array(N).fill(1);
    let nu_e = new Array(N).fill(1);

    // solar constant [W/m²]
    const phi_s = 1361;

    // φ_i = φ_s * μ_i * ν_i  (element-wise)
    let phi_i = mu_i.map((val, i) => phi_s * val * nu_i[i]);

    // Lambertian term
    const A_w = 0.23;
    let f_L = A_w / (4 * Math.PI);

    // f[i] = f_func(f_L, μ_i[i], μ_e[i], α)
    let f = mu_e.map((val, i) => f_lambert(f_L, mu_i[i], val, alpha));

    // I = f * φ_i
    let I = f.map((val, i) => val * phi_i[i]);
    console.log("I: " + I)

    // φ_e = I * μ_e * ν_e
    let phi_e = I.map((val, i) => val * mu_e[i] * nu_e[i]);

    // total = sum(φ_e)
    let phi_eSurf = phi_e.map((val, i) => val * surfaces[i]);
    let total = phi_eSurf.reduce((sum, val) => sum + val, 0);

    return total;
}

function getLightCurve(n = 100, s = [1, 0, 0], o = [0, 0, 1]) {
    var gammas = [];
    var fluxes = [];
    var sv = getSurfVol(vertexes, faces);
    
    // JD at Earth -- the time for which we want to draw the model.
    var jd = Number(document.getElementById("jd").value);
    // Retarded JD.
    var jd_ast_ret = getRetardedJD(jd, o);
    var epsilon = getEpsilon(jd); // in degrees

    // rotate observer and sun vectors
    var rotations = [
        (phi0 + (jd_ast_ret - jd0) / period * 24 * 360) * DEG_TO_RAD,
        (90. - beta) * DEG_TO_RAD,
        lambda * DEG_TO_RAD,
        epsilon * DEG_TO_RAD
    ];
    var model_m = rotate_matrix_inv(rotations);
    s = rotate_vector_inv(s, model_m);
    o = rotate_vector_inv(o, model_m);
    console.log(s)
    console.log(o)

    for (let i = 0; i <= n; i++) {
        // evenly sample gamma from 0 to 2π
        const gamma = 2 * Math.PI * i / n;

        // rotate s about Z by gamma
        var x_ =  s[0] * Math.cos(gamma) + s[1] * Math.sin(gamma);
        var y_ = -s[0] * Math.sin(gamma) + s[1] * Math.cos(gamma);
        var z_ =  s[2];  // still 0
        const s_ = [x_, y_, z_];

        x_ =  o[0] * Math.cos(gamma) + o[1] * Math.sin(gamma);
        y_ = -o[0] * Math.sin(gamma) + o[1] * Math.cos(gamma);
        z_ =  o[2];  // still 0
        const o_ = [x_, y_, z_];

        // compute cosines & shadow masks based on s_ and up
        let flux = getCosinesAndFluxes(s_, o_, sv.surfaces);

        // push a [gamma, flux] pair
        gammas.push(gamma);
        fluxes.push(flux);
    }

    // now normalize columns just like before
    const min0 = Math.min(...gammas), max0 = Math.max(...gammas);
    const min1 = Math.min(...fluxes), max1 = Math.max(...fluxes);
    const range0 = max0 - min0, range1 = max1 - min1;

    gammas.forEach(r => {
        r = (r - min0) / range0;
    });
    fluxes.forEach(r => {
        r = -(r - min1) / range1;

    });

    console.log(gammas);
    console.log(fluxes);
    return [gammas, fluxes];
}

function drawLightCurve(n = 100, s = [1, 0, 0], o = [0, 0, 1]) {
    // get back the pair of arrays
    s = [s[0], s[1], s[2]];
    o = [o[0], o[1], o[2]];

    const curveData = getLightCurve(n, s, o);
    const rawGammas = curveData[0];
    const rawFluxes = curveData[1];

    // compute mins and maxs
    const minG = d3.min(rawGammas), maxG = d3.max(rawGammas);
    const meanF = d3.sum(rawFluxes) / rawFluxes.length;

    // normalize into the desired ranges
    //    phase: [0…1]  => (g - minG)/(maxG - minG)
    //    flux:  [0…2]  => (f - minF)/(maxF - minF) * 2
    const gammas = rawGammas.map(g => (g - minG) / (maxG - minG));
    const fluxes = rawFluxes.map(f => (f / meanF));

    // zip into objects
    const data = gammas.map((g, i) => ({ phase: g, flux: fluxes[i] }));

    // grab the same image-size px that you use for the asteroid canvases:
    const imgSize = +document.getElementById('image_size_px_input').value;

    // resize the SVG to match:
    const svg = d3.select('#fluxChart').attr('width',  imgSize * 1.333333).attr('height', imgSize);

    // set up SVG
    const width  = +svg.attr("width");
    const height = +svg.attr("height");
    const margin = { top: 20, right: 20, bottom: 40, left: 50 };
    const innerW = width  - margin.left - margin.right;
    const innerH = height - margin.top  - margin.bottom;

    svg.selectAll("*").remove();

    // scales — now fixed domains [0,1] and [0,2] and implement padding
    const xExtent = d3.extent(data, d => d.phase);
    const yExtent = [0, 2];

    // pad by 5% of the range
    const xPad = (xExtent[1] - xExtent[0]) * 0.05;
    const yPad = (yExtent[1] - yExtent[0]) * 0.05;

    // build your scales
    const x = d3.scaleLinear()
        .domain([ xExtent[0] - xPad, xExtent[1] + xPad ])
        .nice()
        .range([0, innerW]);

    const y = d3.scaleLinear()
        .domain([ yExtent[0] - yPad, yExtent[1] + yPad ])
        .nice()
        .range([innerH, 0]);

    // container group
    const g = svg.append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);

    // axes – both bottom & top for X, both left & right for Y
    // set a font‐size for tick labels and axis titles
    const tickFontSize = "12px";
    const labelFontSize = "14px";

    // — X bottom
    g.append("g")
    .attr("transform", `translate(0,${innerH})`)
    .call(d3.axisBottom(x).ticks(6).tickSizeInner(6).tickSizeOuter(0))
    .selectAll("text")
        .style("font-size", tickFontSize);

    // — X top
    g.append("g")
    .call(d3.axisTop(x).ticks(6).tickFormat("").tickSizeInner(6).tickSizeOuter(0))
    .selectAll("text")
        .style("font-size", tickFontSize);

    // X label (centered below bottom axis)
    g.append("text")
    .attr("class","axis-label")
    .attr("x", innerW/2)
    .attr("y", innerH + margin.bottom - 5)
    .attr("text-anchor","middle")
    .style("font-size", labelFontSize)
    .text("Phase of rotation");

    // — Y left
    g.append("g")
    .call(d3.axisLeft(y).ticks(6).tickSizeInner(6).tickSizeOuter(0))
    .selectAll("text")
        .style("font-size", tickFontSize);

    // — Y right
    g.append("g")
    .attr("transform", `translate(${innerW},0)`)
    .call(d3.axisRight(y).ticks(6).tickFormat("").tickSizeInner(6).tickSizeOuter(0))
    .selectAll("text")
        .style("font-size", tickFontSize);

    // Y label (rotated, centered alongside left axis)
    g.append("text")
    .attr("class","axis-label")
    .attr("transform", `translate(${-margin.left + 15},${innerH/2}) rotate(-90)`)
    .attr("text-anchor","middle")
    .style("font-size", labelFontSize)
    .text("Relative flux");

    // line generator
    const line = d3.line()
        .x(d => x(d.phase))
        .y(d => y(d.flux))
        .curve(d3.curveMonotoneX);

    // draw path
    g.append("path")
    .datum(data)
    .attr("fill", "none")
    .attr("stroke", "#FF0000")
    .attr("stroke-width", 1.5)
    .attr("d", line);
}

function drawAccurateLightCurve() {
    
}

function saveLightCurveAsSvg() {
    const svg = document.getElementById('fluxChart');
    const serializer = new XMLSerializer();
    const source = serializer.serializeToString(svg);

    // add name‐spaces
    const svgBlob = new Blob([`<?xml version="1.0" standalone="no"?>\n${source}`], {
        type: 'image/svg+xml;charset=utf-8'
    });
    const url = URL.createObjectURL(svgBlob);

    // create a temporary download link and click it
    const a = document.createElement('a');
    a.href = url;
    a.download = 'light_curve.svg';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);

    // cleanup
    URL.revokeObjectURL(url);
}

function saveLightCurveAsTxt() {
    // 1) pull out your data
    const gammas = getLightCurve()[0];
    const fluxes = getLightCurve()[1];
    const pairs = gammas.map((g,i) => [g, fluxes[i]]);
    const text = pairs
    .map(([g, f]) => `${g} ${f}`)
    .join('\n');
  
    // 2) turn it into a Blob
    const blob = new Blob([text], { type: 'text/plain' });
  
    // 3) create an object URL, and a temporary <a> to “click”
    const url = URL.createObjectURL(blob);
    const a   = document.createElement('a');
    a.href        = url;
    a.download    = 'light_curve.txt';   // the suggested filename
    document.body.appendChild(a);
    a.click();
  
    // 4) clean up
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
}
