// generator-core.js (ES module)
// 10x10, K=2 “🌍 Battle” generator + solver.
//
// HARD REQUIREMENT:
// - Accepted levels must be human solvable WITHOUT guessing.
//   => logic propagation ONLY (no lookahead), must fully solve.
// - Must be UNIQUE.
// - Must have an "easy start": at least 1 forced globe on the FIRST propagation pass.
//
// Exports:
// - generateLevel({ requestedDifficulty, levelNumber, onProgress }) -> { level, stats }
// - validateLevel(level) -> { solvable, unique, solutions, nodes }
// - stringifyLevel(level)

const N = 10;
const K = 2;

// ---------- RNG ----------
function mulberry32(seed) {
  let t = seed >>> 0;
  return function () {
    t += 0x6D2B79F5;
    let x = t;
    x = Math.imul(x ^ (x >>> 15), x | 1);
    x ^= x + Math.imul(x ^ (x >>> 7), x | 61);
    return ((x ^ (x >>> 14)) >>> 0) / 4294967296;
  };
}

function shuffleInPlace(arr, rng) {
  for (let i = arr.length - 1; i > 0; i--) {
    const j = Math.floor(rng() * (i + 1));
    [arr[i], arr[j]] = [arr[j], arr[i]];
  }
}

function clamp(n, a, b) {
  return Math.max(a, Math.min(b, n));
}

function randChoiceWeighted(items, weights, rng) {
  let sum = 0;
  for (const w of weights) sum += w;
  if (sum <= 0) return items[Math.floor(rng() * items.length)];
  let r = rng() * sum;
  for (let i = 0; i < items.length; i++) {
    r -= weights[i];
    if (r <= 0) return items[i];
  }
  return items[items.length - 1];
}

// ---------- Row patterns: 2 bits, no horizontal adjacency ----------
function maskToCols(mask) {
  const cols = [];
  for (let c = 0; c < N; c++) {
    if (mask & (1 << c)) cols.push(c);
  }
  return cols;
}

function buildRowPatterns() {
  const patterns = [];
  for (let a = 0; a < N; a++) {
    for (let b = a + 1; b < N; b++) {
      if (b - a === 1) continue;
      let m = 0;
      m |= (1 << a);
      m |= (1 << b);
      patterns.push(m);
    }
  }
  return patterns;
}

const ROW_PATTERNS = buildRowPatterns();
const ROW_PATTERN_COLS = new Map(ROW_PATTERNS.map(m => [m, maskToCols(m)]));

// ---------- Generate a valid full solution (row-by-row) ----------
function generateSolutionMasks(rng) {
  const patterns = ROW_PATTERNS.slice();
  const colCounts = Array(N).fill(0);
  const out = [];

  function rec(row, prevMask) {
    if (row === N) return colCounts.every(v => v === K);

    const rowsLeft = N - row;
    for (let c = 0; c < N; c++) {
      if (colCounts[c] > K) return false;
      if (K - colCounts[c] > rowsLeft) return false;
    }

    shuffleInPlace(patterns, rng);

    for (const m of patterns) {
      if (prevMask != null) {
        if (m & prevMask) continue;
        if ((m << 1) & prevMask) continue;
        if ((m >>> 1) & prevMask) continue;
      }

      const cols = ROW_PATTERN_COLS.get(m);
      let ok = true;
      for (const c of cols) {
        if (colCounts[c] + 1 > K) { ok = false; break; }
      }
      if (!ok) continue;

      for (const c of cols) colCounts[c] += 1;
      out.push(m);

      if (rec(row + 1, m)) return true;

      out.pop();
      for (const c of cols) colCounts[c] -= 1;
    }
    return false;
  }

  return rec(0, null) ? out : null;
}

function solutionMasksToStars(masks) {
  const stars = [];
  for (let r = 0; r < N; r++) {
    const m = masks[r];
    for (let c = 0; c < N; c++) {
      if (m & (1 << c)) stars.push([r, c]);
    }
  }
  return stars;
}

function solutionMasksToGrid(masks) {
  const grid = Array.from({ length: N }, () => Array(N).fill(0));
  for (let r = 0; r < N; r++) {
    const m = masks[r];
    for (let c = 0; c < N; c++) {
      if (m & (1 << c)) grid[r][c] = 1;
    }
  }
  return grid;
}

// ---------- Regions (compact-ish, 10 regions) ----------
const DIR4 = [[1,0],[-1,0],[0,1],[0,-1]];
const DIR8 = [
  [-1,-1], [-1,0], [-1,1],
  [ 0,-1],         [ 0,1],
  [ 1,-1], [ 1,0], [ 1,1],
];

function inBounds(r, c) {
  return r >= 0 && r < N && c >= 0 && c < N;
}

function manhattan(a, b) {
  return Math.abs(a[0] - b[0]) + Math.abs(a[1] - b[1]);
}

function keyOf(rc) {
  return rc[0] + "," + rc[1];
}

function scorePairs(stars, rng) {
  const pairs = [];
  for (let i = 0; i < stars.length; i++) {
    for (let j = i + 1; j < stars.length; j++) {
      const d = manhattan(stars[i], stars[j]);
      pairs.push([d + rng() * 0.25, i, j]);
    }
  }
  pairs.sort((a, b) => a[0] - b[0]);
  return pairs;
}

function pickPairs(stars, rng) {
  const pairsSorted = scorePairs(stars, rng);
  const used = new Set();
  const out = [];
  for (const [, i, j] of pairsSorted) {
    if (used.has(i) || used.has(j)) continue;
    used.add(i); used.add(j);
    out.push([stars[i], stars[j]]);
    if (out.length === N) break;
  }
  return out.length === N ? out : null;
}

function bfsRandomShortest(a, b, blockedSet, occupiedSet, rng) {
  const q = [];
  let qi = 0;
  const prev = new Map();
  prev.set(keyOf(a), null);
  q.push(a);

  while (qi < q.length) {
    const [r, c] = q[qi++];
    if (r === b[0] && c === b[1]) break;

    const dirs = DIR4.slice();
    shuffleInPlace(dirs, rng);

    for (const [dr, dc] of dirs) {
      const nr = r + dr, nc = c + dc;
      if (!inBounds(nr, nc)) continue;
      const nk = nr + "," + nc;
      if (prev.has(nk)) continue;

      const isB = (nr === b[0] && nc === b[1]);
      if (!isB) {
        if (blockedSet.has(nk)) continue;
        if (occupiedSet.has(nk)) continue;
      }

      prev.set(nk, r + "," + c);
      q.push([nr, nc]);
    }
  }

  const endK = keyOf(b);
  if (!prev.has(endK)) return null;

  const path = [];
  let cur = endK;
  while (cur != null) {
    const [rs, cs] = cur.split(",");
    path.push([Number(rs), Number(cs)]);
    cur = prev.get(cur);
  }
  path.reverse();
  return path;
}

function buildRegionTargets(requestedDifficulty, rng) {
  const d = clamp(requestedDifficulty, 1, 10);
  const targets = Array(N).fill(10);

  if (d >= 8) {
    const ids = [...Array(N).keys()];
    shuffleInPlace(ids, rng);
    targets[ids[0]] = 12;
    targets[ids[1]] = 12;
    targets[ids[2]] = 8;
    targets[ids[3]] = 8;
  } else if (d >= 5) {
    const ids = [...Array(N).keys()];
    shuffleInPlace(ids, rng);
    targets[ids[0]] = 11;
    targets[ids[1]] = 11;
    targets[ids[2]] = 11;
    targets[ids[3]] = 9;
    targets[ids[4]] = 9;
    targets[ids[5]] = 9;
  }

  return targets;
}

function fillRegions(regions, targets, maxSize, rng, jaggedness) {
  const sizes = Array(N).fill(0);
  const inFrontier = new Set();
  const frontier = [];

  function addFrontier(r, c) {
    const k = r + "," + c;
    if (inFrontier.has(k)) return;
    inFrontier.add(k);
    frontier.push([r, c]);
  }

  for (let r = 0; r < N; r++) {
    for (let c = 0; c < N; c++) {
      const rid = regions[r][c];
      if (rid !== -1) {
        sizes[rid] += 1;
        for (const [dr, dc] of DIR4) {
          const nr = r + dr, nc = c + dc;
          if (inBounds(nr, nc) && regions[nr][nc] === -1) addFrontier(nr, nc);
        }
      }
    }
  }

  let guard = 0;
  while (frontier.length > 0 && guard < N * N * 140) {
    guard++;

    const pickIndex = Math.floor(rng() * frontier.length);
    const [r, c] = frontier[pickIndex];
    frontier[pickIndex] = frontier[frontier.length - 1];
    frontier.pop();
    inFrontier.delete(r + "," + c);

    if (regions[r][c] !== -1) continue;

    const candidates = [];
    for (const [dr, dc] of DIR4) {
      const nr = r + dr, nc = c + dc;
      if (!inBounds(nr, nc)) continue;
      const rid = regions[nr][nc];
      if (rid !== -1 && !candidates.includes(rid)) candidates.push(rid);
    }
    if (candidates.length === 0) continue;

    const weights = candidates.map(rid => {
      if (sizes[rid] >= maxSize) return 0;
      const deficit = targets[rid] - sizes[rid];
      const base = 0.25 + (1 - jaggedness) * 0.35;
      const gain = Math.max(-2, Math.min(6, deficit));
      return base + (gain + 2) * (0.15 + (1 - jaggedness) * 0.12);
    });

    const chosen = weights.every(w => w <= 0)
      ? candidates[Math.floor(rng() * candidates.length)]
      : randChoiceWeighted(candidates, weights, rng);

    regions[r][c] = chosen;
    sizes[chosen] += 1;

    for (const [dr, dc] of DIR4) {
      const nr = r + dr, nc = c + dc;
      if (inBounds(nr, nc) && regions[nr][nc] === -1) addFrontier(nr, nc);
    }
  }

  for (let r = 0; r < N; r++) for (let c = 0; c < N; c++) {
    if (regions[r][c] === -1) return null;
  }
  return sizes;
}

function buildRegionsFromSolution(stars, rng, requestedDifficulty, maxTries = 600) {
  const d = clamp(requestedDifficulty, 1, 10);
  const jaggedness = clamp((d - 1) / 9, 0, 1);

  const minSize = (d >= 8) ? 7 : (d >= 5) ? 8 : 9;
  const maxSize = (d >= 8) ? 13 : (d >= 5) ? 12 : 11;

  const starSet = new Set(stars.map(keyOf));

  for (let t = 0; t < maxTries; t++) {
    const pairs = pickPairs(stars, rng);
    if (!pairs) continue;

    const regions = Array.from({ length: N }, () => Array(N).fill(-1));
    const occupied = new Set();
    let ok = true;

    for (let rid = 0; rid < pairs.length; rid++) {
      const [a, b] = pairs[rid];

      const blocked = new Set(starSet);
      blocked.delete(keyOf(a));
      blocked.delete(keyOf(b));

      const path = bfsRandomShortest(a, b, blocked, occupied, rng);
      if (!path) { ok = false; break; }
      if (path.length > maxSize) { ok = false; break; }

      for (const [r, c] of path) {
        if (regions[r][c] !== -1 && regions[r][c] !== rid) { ok = false; break; }
        regions[r][c] = rid;
        occupied.add(r + "," + c);
      }
      if (!ok) break;
    }
    if (!ok) continue;

    const targets = buildRegionTargets(requestedDifficulty, rng);
    const sizes = fillRegions(regions, targets, maxSize, rng, jaggedness);
    if (!sizes) continue;

    if (sizes.some(s => s < minSize || s > maxSize)) continue;

    const present = new Set();
    for (let r = 0; r < N; r++) for (let c = 0; c < N; c++) present.add(regions[r][c]);
    if (present.size !== N) continue;

    return regions;
  }

  return null;
}

// ---------- Brute-force solver: count solutions up to maxSolutions ----------
function precomputeRemainingCellsPerRegion(regions) {
  const rem = Array.from({ length: N + 1 }, () => Array(N).fill(0));
  for (let row = N - 1; row >= 0; row--) {
    for (let rid = 0; rid < N; rid++) rem[row][rid] = rem[row + 1][rid];
    for (let c = 0; c < N; c++) rem[row][regions[row][c]] += 1;
  }
  return rem;
}

function countSolutions(regions, maxSolutions = 2) {
  const colCounts = Array(N).fill(0);
  const regCounts = Array(N).fill(0);
  const rem = precomputeRemainingCellsPerRegion(regions);

  let nodes = 0;
  let solutions = 0;

  function rec(row, prevMask) {
    if (solutions >= maxSolutions) return;

    if (row === N) {
      if (colCounts.every(v => v === K) && regCounts.every(v => v === K)) solutions += 1;
      return;
    }

    const rowsLeft = N - row;

    for (let c = 0; c < N; c++) {
      const need = K - colCounts[c];
      if (need < 0 || need > rowsLeft) return;
    }

    for (let rid = 0; rid < N; rid++) {
      const need = K - regCounts[rid];
      if (need < 0 || need > rem[row][rid]) return;
    }

    for (const m of ROW_PATTERNS) {
      nodes += 1;

      if (prevMask != null) {
        if (m & prevMask) continue;
        if ((m << 1) & prevMask) continue;
        if ((m >>> 1) & prevMask) continue;
      }

      const cols = ROW_PATTERN_COLS.get(m);

      let ok = true;
      for (const c of cols) {
        if (colCounts[c] + 1 > K) { ok = false; break; }
      }
      if (!ok) continue;

      const rids = cols.map(c => regions[row][c]);
      for (const rid of rids) {
        if (rid < 0 || rid >= N) { ok = false; break; }
        if (regCounts[rid] + 1 > K) { ok = false; break; }
      }
      if (!ok) continue;

      for (const c of cols) colCounts[c] += 1;
      for (const rid of rids) regCounts[rid] += 1;

      rec(row + 1, m);

      for (const c of cols) colCounts[c] -= 1;
      for (const rid of rids) regCounts[rid] -= 1;

      if (solutions >= maxSolutions) return;
    }
  }

  rec(0, null);
  return { count: solutions, nodes };
}

function difficultyFromNodes(nodes) {
  const s = Math.log10(nodes + 1);
  const d = Math.round((s - 2.2) * 3.0 + 4);
  return clamp(d, 1, 10);
}

// ---------- Logic solver (NO guessing) ----------
function buildUnits(level) {
  const regions = level.regions;

  const rows = [];
  const cols = [];
  const regs = Array.from({ length: N }, () => []);

  for (let r = 0; r < N; r++) {
    const row = [];
    for (let c = 0; c < N; c++) row.push([r, c]);
    rows.push(row);
  }
  for (let c = 0; c < N; c++) {
    const col = [];
    for (let r = 0; r < N; r++) col.push([r, c]);
    cols.push(col);
  }
  for (let r = 0; r < N; r++) {
    for (let c = 0; c < N; c++) {
      regs[regions[r][c]].push([r, c]);
    }
  }
  return { rows, cols, regs };
}

function enumerateCombos(cells, need, canPlaceFn, touchFn) {
  if (need === 0) return [[]];
  const valid = cells.filter(canPlaceFn);
  if (valid.length < need) return [];

  if (need === 1) return valid.map(c => [c]);

  // need === 2
  const combos = [];
  for (let i = 0; i < valid.length; i++) {
    for (let j = i + 1; j < valid.length; j++) {
      const a = valid[i], b = valid[j];
      if (touchFn(a, b)) continue;
      combos.push([a, b]);
    }
  }
  return combos;
}

function makeLogicState(level) {
  const { rows, cols, regs } = buildUnits(level);
  const globe = Array.from({ length: N }, () => Array(N).fill(false));
  const empty = Array.from({ length: N }, () => Array(N).fill(false));

  function setEmpty(r, c) {
    if (globe[r][c]) return false;
    if (empty[r][c]) return true;
    empty[r][c] = true;
    return true;
  }

  function setGlobe(r, c) {
    if (empty[r][c]) return false;
    if (globe[r][c]) return true;
    globe[r][c] = true;
    for (const [dr, dc] of DIR8) {
      const nr = r + dr, nc = c + dc;
      if (inBounds(nr, nc)) {
        if (!setEmpty(nr, nc)) return false;
      }
    }
    return true;
  }

  function countIn(unit, arr) {
    let cnt = 0;
    for (const [r, c] of unit) if (arr[r][c]) cnt++;
    return cnt;
  }

  function canPlace([r, c]) {
    if (empty[r][c]) return false;
    if (globe[r][c]) return true;
    for (const [dr, dc] of DIR8) {
      const nr = r + dr, nc = c + dc;
      if (inBounds(nr, nc) && globe[nr][nc]) return false;
    }
    return true;
  }

  function applyUnit(unit, touchFn) {
    const placed = countIn(unit, globe);
    const need = K - placed;
    if (need < 0) return { ok: false, forcedGlobes: 0 };

    let forcedGlobes = 0;

    if (need === 0) {
      for (const [r, c] of unit) {
        if (!globe[r][c]) {
          if (!setEmpty(r, c)) return { ok: false, forcedGlobes };
        }
      }
      return { ok: true, forcedGlobes };
    }

    const combos = enumerateCombos(unit, need, canPlace, touchFn);
    if (combos.length === 0) return { ok: false, forcedGlobes };

    const union = new Set();
    const inter = new Set(combos[0].map(keyOf));

    for (const combo of combos) {
      const s = new Set(combo.map(keyOf));
      for (const k of s) union.add(k);
      for (const k of [...inter]) if (!s.has(k)) inter.delete(k);
    }

    for (const [r, c] of unit) {
      const k = `${r},${c}`;
      if (globe[r][c]) continue;
      if (!union.has(k)) {
        if (!setEmpty(r, c)) return { ok: false, forcedGlobes };
      }
    }

    for (const k of inter) {
      const [rs, cs] = k.split(",");
      const r = Number(rs), c = Number(cs);
      if (!globe[r][c]) forcedGlobes++;
      if (!setGlobe(r, c)) return { ok: false, forcedGlobes };
    }

    return { ok: true, forcedGlobes };
  }

  function isSolved() {
    for (const u of rows) if (countIn(u, globe) !== K) return false;
    for (const u of cols) if (countIn(u, globe) !== K) return false;
    for (const u of regs) if (countIn(u, globe) !== K) return false;
    return true;
  }

  function touchRow(a, b) { return a[0] === b[0] && Math.abs(a[1] - b[1]) === 1; }
  function touchCol(a, b) { return a[1] === b[1] && Math.abs(a[0] - b[0]) === 1; }
  function touchKing(a, b) { return Math.abs(a[0]-b[0]) <= 1 && Math.abs(a[1]-b[1]) <= 1; }

  function propagateOnce() {
    let forcedGlobes = 0;

    for (const u of rows) {
      const r = applyUnit(u, touchRow);
      if (!r.ok) return { ok: false, forcedGlobes };
      forcedGlobes += r.forcedGlobes;
    }
    for (const u of cols) {
      const r = applyUnit(u, touchCol);
      if (!r.ok) return { ok: false, forcedGlobes };
      forcedGlobes += r.forcedGlobes;
    }
    for (const u of regs) {
      const r = applyUnit(u, touchKing);
      if (!r.ok) return { ok: false, forcedGlobes };
      forcedGlobes += r.forcedGlobes;
    }

    return { ok: true, forcedGlobes };
  }

  return { globe, empty, setGlobe, setEmpty, isSolved, propagateOnce };
}

// NO GUESS: allowOneGuess is gone. Pure propagation only.
function logicSolveNoGuess(level) {
  const st = makeLogicState(level);

  let forcedFirst = 0;
  let changed = true;

  for (let i = 0; i < 300 && changed; i++) {
    const snap = JSON.stringify({ g: st.globe, e: st.empty });
    const res = st.propagateOnce();
    if (!res.ok) return { solved: false, forcedFirst };

    if (i === 0) forcedFirst = res.forcedGlobes;

    if (st.isSolved()) return { solved: true, forcedFirst };

    const now = JSON.stringify({ g: st.globe, e: st.empty });
    changed = (now !== snap);
  }

  return { solved: st.isSolved(), forcedFirst };
}

// ---------- Public API ----------
export function stringifyLevel(levelObj) {
  return JSON.stringify(levelObj, null, 2);
}

export function validateLevel(level) {
  if (!level || level.n !== 10 || level.k !== 2 || !Array.isArray(level.regions) || level.regions.length !== 10) {
    return { solvable: false, unique: false, solutions: 0, nodes: 0 };
  }
  const res = countSolutions(level.regions, 2);
  return { solvable: res.count >= 1, unique: res.count === 1, solutions: res.count, nodes: res.nodes };
}

export async function generateLevel({ requestedDifficulty, levelNumber = 1, onProgress }) {
  // With "no guessing" + uniqueness + easy-start + full-logic-solve, you need more attempts.
  // Humans: "Make it strict but instant." Reality: no.
  const maxAttempts = 120000;

  const seed = (Date.now() ^ (requestedDifficulty * 0x9E3779B9) ^ (levelNumber * 0x85EBCA6B)) >>> 0;
  const rng = mulberry32(seed);

  let best = null;
  let bestStats = null;
  let bestDelta = Infinity;
  let bestForcedFirst = 0;

  // Difficulty tolerance: your node-based estimator is coarse.
  // Tight for higher difficulties, looser for low ones.
  const acceptableDelta = (requestedDifficulty <= 2) ? 2
                        : (requestedDifficulty <= 4) ? 1
                        : 1;

  for (let attempt = 1; attempt <= maxAttempts; attempt++) {
    if (attempt % 50 === 0) await Promise.resolve();

    const sol = generateSolutionMasks(rng);
    if (!sol) continue;

    const stars = solutionMasksToStars(sol);
    const solutionGrid = solutionMasksToGrid(sol);

    const regions = buildRegionsFromSolution(stars, rng, requestedDifficulty);
    if (!regions) continue;

    // Uniqueness check
    const brute = countSolutions(regions, 2);
    if (brute.count !== 1) continue;

    const estimatedDifficulty = difficultyFromNodes(brute.nodes);
    const delta = Math.abs(estimatedDifficulty - requestedDifficulty);

    // Must be solvable by pure logic propagation (no guessing)
    const logic = logicSolveNoGuess({ n: N, k: K, regions });
    if (!logic.solved) continue;

    // Must have an easy start (forced globe immediately)
    if (logic.forcedFirst < 1) continue;

    const better =
      (delta < bestDelta) ||
      (delta === bestDelta && logic.forcedFirst > bestForcedFirst);

    if (better) {
      bestDelta = delta;
      bestForcedFirst = logic.forcedFirst;

      best = {
        version: 1,
        n: N,
        k: K,
        level: levelNumber,
        difficulty: estimatedDifficulty,
        unique: true,
        humanSolvable: true,
        easyStart: true,
        regions,
        solution: solutionGrid
      };

      bestStats = {
        attempts: attempt,
        requestedDifficulty,
        estimatedDifficulty,
        solutions: brute.count,
        nodes: brute.nodes,
        forcedFirst: logic.forcedFirst
      };
    }

    if (onProgress && attempt % 200 === 0) {
      onProgress({
        attempt,
        maxAttempts,
        bestEstimatedDifficulty: best?.difficulty ?? null,
        bestDelta: best ? bestDelta : null,
        bestUnique: best ? true : false,
        bestForcedFirst: best ? bestForcedFirst : null
      });
    }

    // If we found something close enough, stop early.
    if (best && bestDelta <= acceptableDelta) {
      return { level: best, stats: bestStats };
    }
  }

  if (!best) {
    throw new Error(
      "Couldn’t generate a level matching constraints in time. " +
      "With UNIQUE + NO-GUESS LOGIC-SOLVE + EASY-START, this can be extremely rare."
    );
  }

  // Return best found, even if not within acceptableDelta after maxAttempts.
  return { level: best, stats: bestStats };
}
