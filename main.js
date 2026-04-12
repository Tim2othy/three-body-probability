'use strict';

// ════════════════════════════════════════════════════════════════
// CONSTANTS
// ════════════════════════════════════════════════════════════════

const G = 1.0;
const STRIDE = 12; // state vars per member: [x0,y0,x1,y1,x2,y2, vx0,vy0,vx1,vy1,vx2,vy2]
const DT = 0.005;  // integration timestep (fixed)

// Body colors as [r,g,b] in [0,1]
const BODY_RGB = [
    [1.0, 0.20, 0.20],  // red
    [0.15, 0.80, 1.0],  // cyan
    [0.45, 1.0, 0.15],  // lime green
];

// ════════════════════════════════════════════════════════════════
// PRESETS
// ════════════════════════════════════════════════════════════════

function makePresets() {
    return {
        figure8: {
            masses: [1, 1, 1],
            // Chenciner-Montgomery figure-8 orbit (G=1, m=1)
            state: [
                0.97000436, -0.24308753,   // body 0 pos
                0.0, 0.0,          // body 1 pos
                -0.97000436, 0.24308753,   // body 2 pos
                0.46620368, 0.43236573,   // body 0 vel
                -0.93240737, -0.86473146,   // body 1 vel
                0.46620368, 0.43236573,   // body 2 vel
            ],
            scale: 290,
        },

        butterfly: {
            masses: [1, 1, 1],
            // Butterfly orbit initial conditions (Suvakov & Dmitrasinovic 2013)
            state: [
                -1.0, 0.0,
                1.0, 0.0,
                0.0, 0.0,
                0.30689, 0.12551,
                0.30689, 0.12551,
                -0.61378, -0.25102,
            ],
            scale: 260,
        },

        lagrange: {
            masses: [1, 1, 1],
            // Equilateral triangle with angular velocity ω for circular orbit
            // For G=1, m=1, side length L: ω = sqrt(G*m_total/L^3 * L) = sqrt(3G*m/L^3)...
            // Actually for equal masses at equilateral triangle distance L:
            // force on each = G*m^2/L^2 * 2*cos(30°) ... use known result
            state: (() => {
                const L = 1.5;
                const ω = Math.sqrt(3 * G / (L * L * L)); // circular orbit angular vel
                const r = L / Math.sqrt(3); // circumradius
                const s = [];
                for (let i = 0; i < 3; i++) {
                    const θ = (i * 2 * Math.PI) / 3;
                    s.push(r * Math.cos(θ), r * Math.sin(θ));
                }
                for (let i = 0; i < 3; i++) {
                    const θ = (i * 2 * Math.PI) / 3;
                    // velocity = ω × r, perpendicular to position
                    s.push(-ω * r * Math.sin(θ), ω * r * Math.cos(θ));
                }
                return s;
            })(),
            scale: 260,
        },

        hierarchical: {
            masses: [8, 1, 1],
            state: (() => {
                // Heavy body at center, two lighter bodies in different orbits
                const M = 8, m = 1;
                const r1 = 0.6, r2 = 2.2;
                const v1 = Math.sqrt(G * M / r1);
                const v2 = Math.sqrt(G * (M + m) / r2) * 0.88;
                return [
                    0, 0,           // heavy at center
                    r1, 0,           // inner body
                    r2, 0,           // outer body
                    0, 0,           // heavy vel
                    0, v1,          // inner vel
                    0, v2,          // outer vel
                ];
            })(),
            scale: 160,
        },

        random: null, // generated dynamically
    };
}

function generateRandom() {
    const mt = [
        0.5 + Math.random() * 2,
        0.5 + Math.random() * 2,
        0.5 + Math.random() * 2,
    ];
    // Random positions in a disk of radius 1.2
    const pos = [];
    for (let i = 0; i < 3; i++) {
        const r = 0.3 + Math.random() * 0.9;
        const θ = Math.random() * 2 * Math.PI;
        pos.push(r * Math.cos(θ), r * Math.sin(θ));
    }
    // Random velocities, then zero total momentum
    const vel = [];
    let pvx = 0, pvy = 0;
    for (let i = 0; i < 3; i++) {
        const vx = (Math.random() - 0.5) * 1.2;
        const vy = (Math.random() - 0.5) * 1.2;
        vel.push(vx, vy);
        pvx += mt[i] * vx;
        pvy += mt[i] * vy;
    }
    const totalM = mt[0] + mt[1] + mt[2];
    const state = [...pos];
    for (let i = 0; i < 3; i++) {
        state.push(vel[i * 2] - pvx / totalM, vel[i * 2 + 1] - pvy / totalM);
    }
    return { masses: mt, state, scale: 210 };
}

// ════════════════════════════════════════════════════════════════
// PHYSICS
// ════════════════════════════════════════════════════════════════

function deriv(s, m, eps) {
    const x0 = s[0], y0 = s[1], x1 = s[2], y1 = s[3], x2 = s[4], y2 = s[5];
    const vx0 = s[6], vy0 = s[7], vx1 = s[8], vy1 = s[9], vx2 = s[10], vy2 = s[11];

    const dx01 = x1 - x0, dy01 = y1 - y0;
    const dx02 = x2 - x0, dy02 = y2 - y0;
    const dx12 = x2 - x1, dy12 = y2 - y1;

    const eps2 = eps * eps;
    const r01 = Math.sqrt(dx01 * dx01 + dy01 * dy01 + eps2);
    const r02 = Math.sqrt(dx02 * dx02 + dy02 * dy02 + eps2);
    const r12 = Math.sqrt(dx12 * dx12 + dy12 * dy12 + eps2);

    const c01 = G / (r01 * r01 * r01);
    const c02 = G / (r02 * r02 * r02);
    const c12 = G / (r12 * r12 * r12);

    return [
        vx0, vy0, vx1, vy1, vx2, vy2,
        m[1] * c01 * dx01 + m[2] * c02 * dx02,
        m[1] * c01 * dy01 + m[2] * c02 * dy02,
        m[0] * c01 * (-dx01) + m[2] * c12 * dx12,
        m[0] * c01 * (-dy01) + m[2] * c12 * dy12,
        m[0] * c02 * (-dx02) + m[1] * c12 * (-dx12),
        m[0] * c02 * (-dy02) + m[1] * c12 * (-dy12),
    ];
}

// Reusable temp arrays to reduce GC pressure
const _k1 = new Array(STRIDE), _k2 = new Array(STRIDE),
    _k3 = new Array(STRIDE), _k4 = new Array(STRIDE),
    _tmp = new Array(STRIDE);

function rk4Step(s_in, s_out, m, dt, eps) {
    const k1 = _k1, k2 = _k2, k3 = _k3, k4 = _k4, tmp = _tmp;
    const d1 = deriv(s_in, m, eps);
    for (let i = 0; i < STRIDE; i++) { k1[i] = d1[i]; tmp[i] = s_in[i] + 0.5 * dt * k1[i]; }
    const d2 = deriv(tmp, m, eps);
    for (let i = 0; i < STRIDE; i++) { k2[i] = d2[i]; tmp[i] = s_in[i] + 0.5 * dt * k2[i]; }
    const d3 = deriv(tmp, m, eps);
    for (let i = 0; i < STRIDE; i++) { k3[i] = d3[i]; tmp[i] = s_in[i] + dt * k3[i]; }
    const d4 = deriv(tmp, m, eps);
    for (let i = 0; i < STRIDE; i++) {
        k4[i] = d4[i];
        s_out[i] = s_in[i] + (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
}

// ════════════════════════════════════════════════════════════════
// SIMULATION STATE
// ════════════════════════════════════════════════════════════════

let N = 500;
let masses = [1, 1, 1];
let ensemble;    // Float64Array [N * STRIDE]
let ensembleTmp; // Float64Array [N * STRIDE] — swap buffer
let nominal;     // Float64Array [STRIDE]
let nominalTmp;  // Float64Array [STRIDE]
let simTime = 0;
let perturbScale = 1e-4;
let softeningEps = 0.35;

// ════════════════════════════════════════════════════════════════
// INITIALIZATION
// ════════════════════════════════════════════════════════════════

let currentPreset = 'random';
const PRESETS = makePresets();

function gaussianRand() {
    // Box-Muller transform
    const u1 = Math.random() + 1e-15;
    const u2 = Math.random();
    return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}

function initSim(presetName) {
    let preset;
    if (presetName === 'random') {
        preset = generateRandom();
    } else {
        preset = PRESETS[presetName];
    }

    masses = [...preset.masses];
    viewScale = preset.scale;

    nominal = new Float64Array(preset.state);
    nominalTmp = new Float64Array(STRIDE);
    ensemble = new Float64Array(N * STRIDE);
    ensembleTmp = new Float64Array(N * STRIDE);

    const s_in = Array.from(nominal);
    const s_out = new Array(STRIDE);

    for (let e = 0; e < N; e++) {
        const base = e * STRIDE;
        for (let j = 0; j < STRIDE; j++) {
            // Perturb position components (indices 0-5) with Gaussian noise
            const noise = (j < 6) ? gaussianRand() * perturbScale : 0;
            ensemble[base + j] = nominal[j] + noise;
        }
    }

    simTime = 0;
    spreadHistory = [];
    lyapunovEst = null;
    clearDensity();
}

// ════════════════════════════════════════════════════════════════
// STEP
// ════════════════════════════════════════════════════════════════

const _s_in = new Array(STRIDE);
const _s_out = new Array(STRIDE);

function stepAll(dt) {
    const eps = softeningEps;

    // Step each ensemble member
    for (let e = 0; e < N; e++) {
        const base = e * STRIDE;
        for (let j = 0; j < STRIDE; j++) _s_in[j] = ensemble[base + j];
        rk4Step(_s_in, _s_out, masses, dt, eps);
        for (let j = 0; j < STRIDE; j++) ensembleTmp[base + j] = _s_out[j];
    }
    ensemble.set(ensembleTmp.subarray(0, N * STRIDE));

    // Step nominal
    for (let j = 0; j < STRIDE; j++) _s_in[j] = nominal[j];
    rk4Step(_s_in, _s_out, masses, dt, eps);
    for (let j = 0; j < STRIDE; j++) nominal[j] = _s_out[j];
}

// ════════════════════════════════════════════════════════════════
// RENDERING
// ════════════════════════════════════════════════════════════════

const canvas = document.getElementById('main');
const ctx = canvas.getContext('2d');

let W = 0, H = 0;
let imgData;
let densityBuf; // Float32Array [W * H * 3], stores RGB floats

let viewScale = 290;
let viewCX = 0, viewCY = 0; // world-space center
let fadeFactor = 0.97;
let sensitivity = 2.0;
let showNominal = false;
let showTrails = false; // false = snapshot of current positions; true = accumulate history
let showBBox = false;

function resize() {
    W = canvas.width = canvas.offsetWidth;
    H = canvas.height = canvas.offsetHeight;
    imgData = ctx.createImageData(W, H);
    densityBuf = new Float32Array(W * H * 3);
    // Fill alpha channel = 255 permanently
    for (let i = 3; i < imgData.data.length; i += 4) imgData.data[i] = 255;
}

function clearDensity() {
    if (densityBuf) densityBuf.fill(0);
    if (imgData) {
        const d = imgData.data;
        for (let i = 0; i < d.length; i += 4) {
            d[i] = d[i + 1] = d[i + 2] = 0;
            d[i + 3] = 255;
        }
    }
}

function worldToScreen(wx, wy) {
    return [
        (wx - viewCX) * viewScale + W * 0.5,
        -(wy - viewCY) * viewScale + H * 0.5,
    ];
}

// 11×11 Gaussian splat kernel (σ=2.5 px) — [dx, dy, weight]
// Larger kernel means each particle contributes to ~100 pixels instead of 9,
// making the density visible even when the ensemble is sparse or tightly clustered.
const KERNEL = (() => {
    const entries = [];
    const inv2s2 = 1 / (2 * 2.5 * 2.5); // 1/(2σ²), σ=2.5
    for (let dy = -5; dy <= 5; dy++) {
        for (let dx = -5; dx <= 5; dx++) {
            const w = Math.exp(-(dx * dx + dy * dy) * inv2s2);
            if (w > 0.01) entries.push([dx, dy, w]);
        }
    }
    return entries;
})();

function splat(px, py, cr, cg, cb, amount) {
    const ipx = px | 0, ipy = py | 0;
    for (let ki = 0; ki < KERNEL.length; ki++) {
        const kp = KERNEL[ki];
        const x = ipx + kp[0], y = ipy + kp[1];
        if (x < 0 || x >= W || y < 0 || y >= H) continue;
        const idx = (y * W + x) * 3;
        const a = amount * kp[2];
        densityBuf[idx] += cr * a;
        densityBuf[idx + 1] += cg * a;
        densityBuf[idx + 2] += cb * a;
    }
}

function renderDensity() {
    const buf = densityBuf;
    const len = buf.length;

    // Zoom compensation: a pixel covers (1/viewScale)² world area.
    // Zooming in spreads particles across more pixels → fewer hits per pixel → dimmer.
    // Scale perMember by viewScale² so probability density looks the same at any zoom.
    // Reference scale 290 keeps the default zoom calibrated to sensitivity=1 meaning
    // "all N particles overlapping = saturated".
    const zoomComp = (viewScale / 290) ** 2;

    let perMember;
    if (showTrails) {
        // Trail/history mode: decay old density and accumulate current frame on top.
        const fade = fadeFactor;
        for (let i = 0; i < len; i++) buf[i] *= fade;
        perMember = sensitivity * (1 - fade) * 255 / N * zoomComp;
    } else {
        // Snapshot mode: clear buffer and show only WHERE ensemble members ARE RIGHT NOW.
        buf.fill(0);
        perMember = sensitivity * 255 / N * zoomComp;
    }

    for (let e = 0; e < N; e++) {
        const base = e * STRIDE;
        for (let b = 0; b < 3; b++) {
            const wx = ensemble[base + b * 2];
            const wy = ensemble[base + b * 2 + 1];
            const [sx, sy] = worldToScreen(wx, wy);
            if (sx < -6 || sx >= W + 6 || sy < -6 || sy >= H + 6) continue;
            const col = BODY_RGB[b];
            splat(sx, sy, col[0], col[1], col[2], perMember);
        }
    }

    // Gamma tone-map: lower exponent = more aggressive boost of dim regions.
    // 0.35 makes even 1 particle out of 500 clearly visible against the black background.
    const data = imgData.data;
    for (let i = 0; i < W * H; i++) {
        const bi = i * 3, pi = i * 4;
        const r = buf[bi], g = buf[bi + 1], b = buf[bi + 2];
        const scale = 1 / 255;
        data[pi] = Math.min(255, Math.pow(Math.min(1, r * scale), 0.35) * 255) | 0;
        data[pi + 1] = Math.min(255, Math.pow(Math.min(1, g * scale), 0.35) * 255) | 0;
        data[pi + 2] = Math.min(255, Math.pow(Math.min(1, b * scale), 0.35) * 255) | 0;
    }

    ctx.putImageData(imgData, 0, 0);
}

function renderBBoxes() {
    const pct = 0.05; // clip bottom/top 5% → 90% box
    for (let b = 0; b < 3; b++) {
        const xs = [], ys = [];
        for (let e = 0; e < N; e++) {
            xs.push(ensemble[e * STRIDE + b * 2]);
            ys.push(ensemble[e * STRIDE + b * 2 + 1]);
        }
        xs.sort((a, z) => a - z);
        ys.sort((a, z) => a - z);
        const lo = Math.floor(pct * N);
        const hi = Math.ceil((1 - pct) * N) - 1;
        const [x0, y0] = worldToScreen(xs[lo], ys[hi]); // top-left (ys[hi] = maxY = top in screen)
        const [x1, y1] = worldToScreen(xs[hi], ys[lo]); // bottom-right
        const [cr, cg, cb] = BODY_RGB[b];
        ctx.strokeStyle = `rgba(${(cr * 255) | 0},${(cg * 255) | 0},${(cb * 255) | 0},0.85)`;
        ctx.lineWidth = 2;
        ctx.setLineDash([6, 4]);
        ctx.strokeRect(x0, y0, x1 - x0, y1 - y0);
    }
    ctx.setLineDash([]);
}

function renderNominal() {
    const s = nominal;
    for (let b = 0; b < 3; b++) {
        const [sx, sy] = worldToScreen(s[b * 2], s[b * 2 + 1]);
        const [cr, cg, cb] = BODY_RGB[b];
        const hex = `rgb(${(cr * 255) | 0},${(cg * 255) | 0},${(cb * 255) | 0})`;
        const sz = 7;
        ctx.lineWidth = 1.5;
        ctx.strokeStyle = hex;
        ctx.beginPath();
        ctx.moveTo(sx - sz, sy); ctx.lineTo(sx + sz, sy);
        ctx.moveTo(sx, sy - sz); ctx.lineTo(sx, sy + sz);
        ctx.stroke();
        // Bright center dot
        ctx.fillStyle = '#fff';
        ctx.beginPath();
        ctx.arc(sx, sy, 2, 0, Math.PI * 2);
        ctx.fill();
    }
}

function clearCanvas() {
    ctx.fillStyle = '#080810';
    ctx.fillRect(0, 0, W, H);
}

// ════════════════════════════════════════════════════════════════
// STATISTICS
// ════════════════════════════════════════════════════════════════

let spreadHistory = [];
let lyapunovEst = null;

function computeStats() {
    // Compute max spread across all 3 bodies
    let maxSpread = 0;
    for (let b = 0; b < 3; b++) {
        let cx = 0, cy = 0;
        for (let e = 0; e < N; e++) {
            cx += ensemble[e * STRIDE + b * 2];
            cy += ensemble[e * STRIDE + b * 2 + 1];
        }
        cx /= N; cy /= N;
        for (let e = 0; e < N; e++) {
            const dx = ensemble[e * STRIDE + b * 2] - cx;
            const dy = ensemble[e * STRIDE + b * 2 + 1] - cy;
            const d = Math.sqrt(dx * dx + dy * dy);
            if (d > maxSpread) maxSpread = d;
        }
    }

    // Lyapunov estimate using log growth of spread
    if (spreadHistory.length > 0 && maxSpread > 0 && simTime > 0.2) {
        const old = spreadHistory[0];
        if (old.spread > 1e-12) {
            const dt_hist = simTime - old.t;
            if (dt_hist > 0.5) {
                lyapunovEst = Math.log(maxSpread / old.spread) / dt_hist;
                spreadHistory.shift();
            }
        }
    }
    spreadHistory.push({ t: simTime, spread: maxSpread });
    if (spreadHistory.length > 80) spreadHistory.shift();

    // Update DOM
    document.getElementById('stat-time').textContent = simTime.toFixed(2);
    document.getElementById('stat-spread').textContent = maxSpread.toFixed(4);
    document.getElementById('stat-lyap').textContent =
        (lyapunovEst !== null) ? lyapunovEst.toFixed(3) : '—';
    document.getElementById('stat-active').textContent = N;

    const chaosNorm = Math.min(1, maxSpread / 4);
    document.getElementById('chaos-fill').style.width = (chaosNorm * 100).toFixed(1) + '%';
}

// ════════════════════════════════════════════════════════════════
// ANIMATION LOOP
// ════════════════════════════════════════════════════════════════

let running = true;
let stepsPerFrame = 8;
let frameCount = 0;

function frame() {
    if (running && ensemble) {
        for (let i = 0; i < stepsPerFrame; i++) {
            stepAll(DT);
            simTime += DT;
        }

        renderDensity(); // always render current probability mass (snapshot or trail mode)
        if (showBBox) renderBBoxes();
        if (showNominal) renderNominal();

        frameCount++;
        if (frameCount % 12 === 0) computeStats();
    }

    requestAnimationFrame(frame);
}

// ════════════════════════════════════════════════════════════════
// CONTROLS
// ════════════════════════════════════════════════════════════════

function slider(id, valId, onChange, fmt) {
    const sl = document.getElementById(id);
    const vl = document.getElementById(valId);
    const update = () => {
        const v = parseFloat(sl.value);
        vl.textContent = fmt ? fmt(v) : v;
        onChange(v);
    };
    sl.addEventListener('input', update);
    // Initialize display
    vl.textContent = fmt ? fmt(parseFloat(sl.value)) : sl.value;
}

function resetSim() {
    spreadHistory = [];
    lyapunovEst = null;
    initSim(currentPreset);
}

function setup() {
    resize();
    window.addEventListener('resize', () => { resize(); clearDensity(); });

    initSim(currentPreset);

    // ── Play / Pause
    const btnPlay = document.getElementById('btn-play');
    btnPlay.addEventListener('click', function () {
        running = !running;
        this.textContent = running ? '⏸ Pause' : '▶ Play';
        this.classList.toggle('active', running);
    });

    // ── Reset
    document.getElementById('btn-reset').addEventListener('click', resetSim);

    // ── Preset
    document.getElementById('preset-select').addEventListener('change', function () {
        currentPreset = this.value;
        resetSim();
    });

    // ── Ensemble N
    slider('sl-n', 'val-n', v => {
        N = Math.round(v);
        // Re-allocate buffers and re-init
        resetSim();
    });

    // ── Perturbation (log scale)
    slider('sl-perturb', 'val-perturb', v => {
        perturbScale = Math.pow(10, v);
        resetSim();
    }, v => {
        const exp = Math.round(v);
        const man = Math.pow(10, v - exp);
        return man.toFixed(1) + 'e' + exp;
    });

    // ── Speed
    slider('sl-speed', 'val-speed', v => {
        stepsPerFrame = Math.round(v);
    }, v => Math.round(v) + '×');

    // ── Fade
    slider('sl-fade', 'val-fade', v => {
        fadeFactor = v;
    }, v => v.toFixed(3));

    // ── Sensitivity
    slider('sl-bright', 'val-bright', v => {
        sensitivity = v;
    }, v => v.toFixed(1));

    // ── Softening
    slider('sl-soft', 'val-soft', v => {
        softeningEps = v;
    }, v => v.toFixed(3));

    // ── Zoom
    document.getElementById('btn-zoom-in').addEventListener('click', () => {
        viewScale *= 1.4; clearDensity();
    });
    document.getElementById('btn-zoom-out').addEventListener('click', () => {
        viewScale /= 1.4; clearDensity();
    });
    document.getElementById('btn-center').addEventListener('click', () => {
        viewCX = 0; viewCY = 0; clearDensity();
    });

    // ── Checkboxes
    document.getElementById('cb-bbox').addEventListener('change', function () {
        showBBox = this.checked;
    });
    document.getElementById('cb-nominal').addEventListener('change', function () {
        showNominal = this.checked;
    });
    document.getElementById('cb-trails').addEventListener('change', function () {
        showTrails = this.checked;
        if (!showTrails) clearDensity();
    });

    // ── Pan
    let dragging = false, lastMX = 0, lastMY = 0;
    canvas.addEventListener('mousedown', e => {
        dragging = true; lastMX = e.clientX; lastMY = e.clientY;
    });
    window.addEventListener('mouseup', () => { dragging = false; });
    canvas.addEventListener('mousemove', e => {
        if (!dragging) return;
        viewCX -= (e.clientX - lastMX) / viewScale;
        viewCY += (e.clientY - lastMY) / viewScale;
        lastMX = e.clientX; lastMY = e.clientY;
        clearDensity();
    });

    // ── Scroll zoom
    canvas.addEventListener('wheel', e => {
        e.preventDefault();
        const factor = e.deltaY > 0 ? 1 / 1.12 : 1.12;
        viewScale *= factor;
        clearDensity();
    }, { passive: false });

    // ── Keyboard shortcuts
    window.addEventListener('keydown', e => {
        if (e.key === ' ') {
            e.preventDefault();
            document.getElementById('btn-play').click();
        }
        if (e.key === 'r' || e.key === 'R') resetSim();
        if (e.key === '+' || e.key === '=') { viewScale *= 1.2; clearDensity(); }
        if (e.key === '-') { viewScale /= 1.2; clearDensity(); }
        if (e.key === 'c' || e.key === 'C') { viewCX = 0; viewCY = 0; clearDensity(); }
    });

    frame();
}

setup();
