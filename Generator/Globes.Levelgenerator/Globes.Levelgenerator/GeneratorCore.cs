using System.Diagnostics;
using System.Security.Cryptography;

namespace LevelGen;

public static class GeneratorCore
{
    private const int N = 10;
    private const int K = 2;

    private static readonly (int dr, int dc)[] DIR4 =
    [
        (1, 0), (-1, 0), (0, 1), (0, -1)
    ];

    private static readonly (int dr, int dc)[] DIR8 =
    [
        (-1,-1), (-1,0), (-1,1),
        ( 0,-1),         ( 0,1),
        ( 1,-1), ( 1,0), ( 1,1),
    ];

    private static readonly int[] ROW_PATTERNS = BuildRowPatterns();
    private static readonly Dictionary<int, int[]> ROW_PATTERN_COLS =
        ROW_PATTERNS.ToDictionary(m => m, MaskToCols);

    private static int CellId(int r, int c) => r * N + c;
    private static (int r, int c) FromId(int id) => (id / N, id % N);

    public static async Task<LevelObject> GenerateRawLevelAsync(int requestedDifficulty, int levelNumber, Action<string>? progress)
    {
        requestedDifficulty = Utils.Clamp(requestedDifficulty, 1, 10);

        const int maxAttemptsBeforeReportSpam = int.MaxValue;

        uint seed = (uint)RandomNumberGenerator.GetInt32(int.MinValue, int.MaxValue);
        seed ^= (uint)(requestedDifficulty * unchecked((int)0x9E3779B9));
        seed ^= (uint)(levelNumber * unchecked((int)0x85EBCA6B));
        var rng = new Mulberry32(seed);

        int attempts = 0;
        int madeSolution = 0;
        int madeRegions = 0;
        string lastFail = "init";

        var sw = Stopwatch.StartNew();
        long lastReport = -9999;

        while (attempts < maxAttemptsBeforeReportSpam)
        {
            attempts++;

            if (attempts % 100 == 0)
                await Task.Yield();

            var solMasks = GenerateSolutionMasks(rng);
            if (solMasks is null)
            {
                lastFail = "solution";
                Report();
                continue;
            }
            madeSolution++;

            var stars = SolutionMasksToStars(solMasks);
            var solGrid = SolutionMasksToGrid(solMasks);

            var regions = BuildRegionsFromSolution(stars, rng, requestedDifficulty, maxTries: 600);
            if (regions is null)
            {
                lastFail = "regions";
                Report();
                continue;
            }
            madeRegions++;

            Report(force: true);

            return new LevelObject
            {
                Version = 1,
                N = N,
                K = K,
                Level = levelNumber,
                Difficulty = 0,
                Unique = false,
                HumanSolvable = false,
                EasyStart = false,
                Regions = regions,
                Solution = solGrid
            };
        }

        throw new InvalidOperationException("Raw generation loop ended unexpectedly.");

        void Report(bool force = false)
        {
            if (progress is null) return;

            long ms = sw.ElapsedMilliseconds;
            if (!force && ms - lastReport < 500) return;
            lastReport = ms;

            progress($"[GEN] Level-Number={levelNumber} Solution-Ok={madeSolution} Regions-Ok={madeRegions}");
        }
    }

    public static SolvabilityInfo CheckHumanSolvable(LevelObject level, bool requireEasyStart)
    {
        var logic = LogicSolveNoGuess(level);

        if (!logic.Solved)
        {
            return new SolvabilityInfo
            {
                Solved = false,
                EasyStart = false,
                ForcedFirst = logic.ForcedFirst,
                Solutions = 0,
                Nodes = 0,
                FailReason = "logic"
            };
        }

        bool easy = logic.ForcedFirst >= 1;

        if (requireEasyStart && !easy)
        {
            return new SolvabilityInfo
            {
                Solved = true,
                EasyStart = false,
                ForcedFirst = logic.ForcedFirst,
                Solutions = 0,
                Nodes = 0,
                FailReason = "easyStart"
            };
        }

        var brute = CountSolutions(level.Regions, maxSolutions: 2);

        level.HumanSolvable = true;
        level.EasyStart = easy;
        level.Unique = (brute.Count == 1);

        return new SolvabilityInfo
        {
            Solved = true,
            EasyStart = easy,
            ForcedFirst = logic.ForcedFirst,
            Solutions = brute.Count,
            Nodes = brute.Nodes,
            FailReason = "PASS"
        };
    }

    public static int DifficultyFromNodes(int nodes)
    {
        double s = Math.Log10(nodes + 1);
        int d = (int)Math.Round((s - 2.2) * 3.0 + 4);
        return Utils.Clamp(d, 1, 10);
    }

    // ----------------------------
    // RNG
    // ----------------------------
    private sealed class Mulberry32
    {
        private uint _t;
        public Mulberry32(uint seed) => _t = seed;

        public double NextDouble()
        {
            _t += 0x6D2B79F5;
            uint x = _t;
            x = (uint)(unchecked((int)(x ^ (x >> 15))) * unchecked((int)(x | 1)));
            x ^= x + (uint)(unchecked((int)(x ^ (x >> 7))) * unchecked((int)(x | 61)));
            return ((x ^ (x >> 14)) & 0xFFFFFFFFu) / 4294967296.0;
        }

        public void Shuffle<T>(IList<T> arr)
        {
            for (int i = arr.Count - 1; i > 0; i--)
            {
                int j = (int)Math.Floor(NextDouble() * (i + 1));
                (arr[i], arr[j]) = (arr[j], arr[i]);
            }
        }
    }

    // ----------------------------
    // Row patterns / solution gen
    // ----------------------------
    private static int[] MaskToCols(int mask)
    {
        var cols = new List<int>(K);
        for (int c = 0; c < N; c++)
            if ((mask & (1 << c)) != 0) cols.Add(c);
        return cols.ToArray();
    }

    private static int[] BuildRowPatterns()
    {
        var patterns = new List<int>();
        for (int a = 0; a < N; a++)
            for (int b = a + 1; b < N; b++)
            {
                if (b - a == 1) continue;
                patterns.Add((1 << a) | (1 << b));
            }
        return patterns.ToArray();
    }

    private static int[]? GenerateSolutionMasks(Mulberry32 rng)
    {
        var colCounts = new int[N];
        var outMasks = new List<int>(N);

        bool Rec(int row, int? prevMask)
        {
            if (row == N) return colCounts.All(v => v == K);

            int rowsLeft = N - row;

            for (int c = 0; c < N; c++)
            {
                if (colCounts[c] > K) return false;
                if (K - colCounts[c] > rowsLeft) return false;
            }

            // IMPORTANT: shuffle a copy, not a shared list (prevents "collection modified" hell)
            var order = (int[])ROW_PATTERNS.Clone();
            rng.Shuffle(order);

            foreach (var m in order)
            {
                if (prevMask is not null)
                {
                    int pm = prevMask.Value;
                    if ((m & pm) != 0) continue;
                    if (((m << 1) & pm) != 0) continue;
                    if (((m >> 1) & pm) != 0) continue;
                }

                var cols = ROW_PATTERN_COLS[m];

                bool ok = true;
                foreach (var c in cols)
                    if (colCounts[c] + 1 > K) { ok = false; break; }
                if (!ok) continue;

                foreach (var c in cols) colCounts[c]++;
                outMasks.Add(m);

                if (Rec(row + 1, m)) return true;

                outMasks.RemoveAt(outMasks.Count - 1);
                foreach (var c in cols) colCounts[c]--;
            }

            return false;
        }

        return Rec(0, null) ? outMasks.ToArray() : null;
    }

    private static List<int> SolutionMasksToStars(int[] masks)
    {
        var stars = new List<int>(N * K);
        for (int r = 0; r < N; r++)
        {
            int m = masks[r];
            for (int c = 0; c < N; c++)
                if ((m & (1 << c)) != 0) stars.Add(CellId(r, c));
        }
        return stars;
    }

    private static int[][] SolutionMasksToGrid(int[] masks)
    {
        var grid = Enumerable.Range(0, N).Select(_ => new int[N]).ToArray();
        for (int r = 0; r < N; r++)
        {
            int m = masks[r];
            for (int c = 0; c < N; c++)
                if ((m & (1 << c)) != 0) grid[r][c] = 1;
        }
        return grid;
    }

    // ----------------------------
    // Regions
    // ----------------------------
    private static bool InBounds(int r, int c) => r >= 0 && r < N && c >= 0 && c < N;

    private static int ManhattanId(int aId, int bId)
    {
        var (ar, ac) = FromId(aId);
        var (br, bc) = FromId(bId);
        return Math.Abs(ar - br) + Math.Abs(ac - bc);
    }

    private static List<(double score, int i, int j)> ScorePairs(List<int> stars, Mulberry32 rng)
    {
        var pairs = new List<(double score, int i, int j)>();
        for (int i = 0; i < stars.Count; i++)
            for (int j = i + 1; j < stars.Count; j++)
            {
                double d = ManhattanId(stars[i], stars[j]);
                pairs.Add((d + rng.NextDouble() * 0.25, i, j));
            }

        pairs.Sort((a, b) => a.score.CompareTo(b.score));
        return pairs;
    }

    private static List<(int aId, int bId)>? PickPairs(List<int> stars, Mulberry32 rng)
    {
        var sorted = ScorePairs(stars, rng);
        var used = new HashSet<int>();
        var outPairs = new List<(int, int)>(N);

        foreach (var (_, i, j) in sorted)
        {
            if (used.Contains(i) || used.Contains(j)) continue;
            used.Add(i); used.Add(j);
            outPairs.Add((stars[i], stars[j]));
            if (outPairs.Count == N) break;
        }

        return outPairs.Count == N ? outPairs : null;
    }

    private static List<int>? BfsRandomShortest(int aId, int bId, HashSet<int> blocked, HashSet<int> occupied, Mulberry32 rng)
    {
        var q = new List<int>();
        int qi = 0;

        var prev = new Dictionary<int, int?>(256);
        prev[aId] = null;
        q.Add(aId);

        while (qi < q.Count)
        {
            int cur = q[qi++];
            if (cur == bId) break;

            var (r, c) = FromId(cur);

            var dirs = (ValueTuple<int, int>[])DIR4.Clone();
            for (int i = dirs.Length - 1; i > 0; i--)
            {
                int j = (int)Math.Floor(rng.NextDouble() * (i + 1));
                (dirs[i], dirs[j]) = (dirs[j], dirs[i]);
            }

            foreach (var (dr, dc) in dirs)
            {
                int nr = r + dr, nc = c + dc;
                if (!InBounds(nr, nc)) continue;

                int nid = CellId(nr, nc);
                if (prev.ContainsKey(nid)) continue;

                bool isB = (nid == bId);
                if (!isB)
                {
                    if (blocked.Contains(nid)) continue;
                    if (occupied.Contains(nid)) continue;
                }

                prev[nid] = cur;
                q.Add(nid);
            }
        }

        if (!prev.ContainsKey(bId)) return null;

        var path = new List<int>();
        int? p = bId;
        while (p is not null)
        {
            path.Add(p.Value);
            p = prev[p.Value];
        }
        path.Reverse();
        return path;
    }

    private static int[] BuildRegionTargets(int requestedDifficulty, Mulberry32 rng)
    {
        int d = Utils.Clamp(requestedDifficulty, 1, 10);
        var targets = Enumerable.Repeat(10, N).ToArray();

        if (d >= 8)
        {
            var ids = Enumerable.Range(0, N).ToList();
            rng.Shuffle(ids);
            targets[ids[0]] = 12;
            targets[ids[1]] = 12;
            targets[ids[2]] = 8;
            targets[ids[3]] = 8;
        }
        else if (d >= 5)
        {
            var ids = Enumerable.Range(0, N).ToList();
            rng.Shuffle(ids);
            targets[ids[0]] = 11;
            targets[ids[1]] = 11;
            targets[ids[2]] = 11;
            targets[ids[3]] = 9;
            targets[ids[4]] = 9;
            targets[ids[5]] = 9;
        }

        return targets;
    }

    private static int[]? FillRegions(int[][] regions, int[] targets, int maxSize, Mulberry32 rng, double jaggedness)
    {
        var sizes = new int[N];
        var inFrontier = new HashSet<int>();
        var frontier = new List<int>();

        void AddFrontier(int r, int c)
        {
            int id = CellId(r, c);
            if (inFrontier.Add(id))
                frontier.Add(id);
        }

        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++)
            {
                int rid = regions[r][c];
                if (rid != -1)
                {
                    sizes[rid]++;
                    foreach (var (dr, dc) in DIR4)
                    {
                        int nr = r + dr, nc = c + dc;
                        if (InBounds(nr, nc) && regions[nr][nc] == -1) AddFrontier(nr, nc);
                    }
                }
            }

        int guard = 0;
        while (frontier.Count > 0 && guard < N * N * 140)
        {
            guard++;

            int pickIndex = (int)Math.Floor(rng.NextDouble() * frontier.Count);
            int id = frontier[pickIndex];
            frontier[pickIndex] = frontier[^1];
            frontier.RemoveAt(frontier.Count - 1);
            inFrontier.Remove(id);

            var (r, c) = FromId(id);

            if (regions[r][c] != -1) continue;

            var candidates = new List<int>();
            foreach (var (dr, dc) in DIR4)
            {
                int nr = r + dr, nc = c + dc;
                if (!InBounds(nr, nc)) continue;

                int rid = regions[nr][nc];
                if (rid != -1 && !candidates.Contains(rid)) candidates.Add(rid);
            }
            if (candidates.Count == 0) continue;

            var weights = candidates.Select(rid =>
            {
                if (sizes[rid] >= maxSize) return 0.0;
                int deficit = targets[rid] - sizes[rid];
                double @base = 0.25 + (1 - jaggedness) * 0.35;
                int gain = Math.Max(-2, Math.Min(6, deficit));
                return @base + (gain + 2) * (0.15 + (1 - jaggedness) * 0.12);
            }).ToArray();

            int chosen = weights.All(w => w <= 0)
                ? candidates[(int)Math.Floor(rng.NextDouble() * candidates.Count)]
                : RandChoiceWeighted(candidates, weights, rng);

            regions[r][c] = chosen;
            sizes[chosen]++;

            foreach (var (dr, dc) in DIR4)
            {
                int nr = r + dr, nc = c + dc;
                if (InBounds(nr, nc) && regions[nr][nc] == -1) AddFrontier(nr, nc);
            }
        }

        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++)
                if (regions[r][c] == -1) return null;

        return sizes;
    }

    private static int RandChoiceWeighted(List<int> items, double[] weights, Mulberry32 rng)
    {
        double sum = 0;
        for (int i = 0; i < weights.Length; i++) sum += weights[i];

        if (sum <= 0)
            return items[(int)Math.Floor(rng.NextDouble() * items.Count)];

        double r = rng.NextDouble() * sum;
        for (int i = 0; i < items.Count; i++)
        {
            r -= weights[i];
            if (r <= 0) return items[i];
        }
        return items[^1];
    }

    private static int[][]? BuildRegionsFromSolution(List<int> stars, Mulberry32 rng, int requestedDifficulty, int maxTries = 600)
    {
        int d = Utils.Clamp(requestedDifficulty, 1, 10);
        double jaggedness = Utils.Clamp((d - 1) / 9.0, 0.0, 1.0);

        int minSize = (d >= 8) ? 7 : (d >= 5) ? 8 : 9;
        int maxSize = (d >= 8) ? 13 : (d >= 5) ? 12 : 11;

        var starSet = new HashSet<int>(stars);

        for (int t = 0; t < maxTries; t++)
        {
            var pairs = PickPairs(stars, rng);
            if (pairs is null) continue;

            var regions = Enumerable.Range(0, N).Select(_ => Enumerable.Repeat(-1, N).ToArray()).ToArray();
            var occupied = new HashSet<int>();
            bool ok = true;

            for (int rid = 0; rid < pairs.Count; rid++)
            {
                var (aId, bId) = pairs[rid];

                var blocked = new HashSet<int>(starSet);
                blocked.Remove(aId);
                blocked.Remove(bId);

                var path = BfsRandomShortest(aId, bId, blocked, occupied, rng);
                if (path is null) { ok = false; break; }
                if (path.Count > maxSize) { ok = false; break; }

                foreach (var pid in path)
                {
                    var (r, c) = FromId(pid);

                    if (regions[r][c] != -1 && regions[r][c] != rid) { ok = false; break; }
                    regions[r][c] = rid;
                    occupied.Add(pid);
                }

                if (!ok) break;
            }

            if (!ok) continue;

            var targets = BuildRegionTargets(requestedDifficulty, rng);
            var sizes = FillRegions(regions, targets, maxSize, rng, jaggedness);
            if (sizes is null) continue;

            if (sizes.Any(s => s < minSize || s > maxSize)) continue;

            var present = new HashSet<int>();
            for (int r = 0; r < N; r++)
                for (int c = 0; c < N; c++)
                    present.Add(regions[r][c]);

            if (present.Count != N) continue;

            return regions;
        }

        return null;
    }

    // ----------------------------
    // Brute solver (uniqueness + nodes)
    // ----------------------------
    private static int[][] PrecomputeRemainingCellsPerRegion(int[][] regions)
    {
        var rem = Enumerable.Range(0, N + 1).Select(_ => new int[N]).ToArray();

        for (int row = N - 1; row >= 0; row--)
        {
            for (int rid = 0; rid < N; rid++)
                rem[row][rid] = rem[row + 1][rid];

            for (int c = 0; c < N; c++)
                rem[row][regions[row][c]] += 1;
        }

        return rem;
    }

    private sealed record CountRes(int Count, int Nodes);

    private static CountRes CountSolutions(int[][] regions, int maxSolutions = 2)
    {
        var colCounts = new int[N];
        var regCounts = new int[N];
        var rem = PrecomputeRemainingCellsPerRegion(regions);

        int nodes = 0;
        int solutions = 0;

        void Rec(int row, int? prevMask)
        {
            if (solutions >= maxSolutions) return;

            if (row == N)
            {
                if (colCounts.All(v => v == K) && regCounts.All(v => v == K))
                    solutions++;
                return;
            }

            int rowsLeft = N - row;

            for (int c = 0; c < N; c++)
            {
                int need = K - colCounts[c];
                if (need < 0 || need > rowsLeft) return;
            }

            for (int rid = 0; rid < N; rid++)
            {
                int need = K - regCounts[rid];
                if (need < 0 || need > rem[row][rid]) return;
            }

            foreach (var m in ROW_PATTERNS)
            {
                nodes++;

                if (prevMask is not null)
                {
                    int pm = prevMask.Value;
                    if ((m & pm) != 0) continue;
                    if (((m << 1) & pm) != 0) continue;
                    if (((m >> 1) & pm) != 0) continue;
                }

                var cols = ROW_PATTERN_COLS[m];

                bool ok = true;
                foreach (var c in cols)
                    if (colCounts[c] + 1 > K) { ok = false; break; }
                if (!ok) continue;

                var rids = cols.Select(c => regions[row][c]).ToArray();
                foreach (var rid in rids)
                {
                    if (rid < 0 || rid >= N) { ok = false; break; }
                    if (regCounts[rid] + 1 > K) { ok = false; break; }
                }
                if (!ok) continue;

                foreach (var c in cols) colCounts[c]++;
                foreach (var rid in rids) regCounts[rid]++;

                Rec(row + 1, m);

                foreach (var c in cols) colCounts[c]--;
                foreach (var rid in rids) regCounts[rid]--;

                if (solutions >= maxSolutions) return;
            }
        }

        Rec(0, null);
        return new CountRes(solutions, nodes);
    }

    // ----------------------------
    // Logic solver (NO guessing)
    // ----------------------------
    private sealed class Units
    {
        public required List<int>[] Rows { get; init; }
        public required List<int>[] Cols { get; init; }
        public required List<int>[] Regs { get; init; }
    }

    private static Units BuildUnits(LevelObject level)
    {
        var regions = level.Regions;

        var rows = new List<int>[N];
        var cols = new List<int>[N];
        var regs = Enumerable.Range(0, N).Select(_ => new List<int>()).ToArray();

        for (int r = 0; r < N; r++)
        {
            rows[r] = new List<int>(N);
            for (int c = 0; c < N; c++) rows[r].Add(CellId(r, c));
        }

        for (int c = 0; c < N; c++)
        {
            cols[c] = new List<int>(N);
            for (int r = 0; r < N; r++) cols[c].Add(CellId(r, c));
        }

        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++)
                regs[regions[r][c]].Add(CellId(r, c));

        return new Units { Rows = rows, Cols = cols, Regs = regs };
    }

    private static List<int[]> EnumerateCombos(List<int> cells, int need, Func<int, bool> canPlace, Func<int, int, bool> touchFn)
    {
        if (need == 0) return [Array.Empty<int>()];

        var valid = cells.Where(canPlace).ToList();
        if (valid.Count < need) return [];

        if (need == 1) return valid.Select(v => new[] { v }).ToList();

        var combos = new List<int[]>();
        for (int i = 0; i < valid.Count; i++)
            for (int j = i + 1; j < valid.Count; j++)
            {
                int a = valid[i], b = valid[j];
                if (touchFn(a, b)) continue;
                combos.Add(new[] { a, b });
            }

        return combos;
    }

    private sealed class LogicState
    {
        public bool[][] Globe { get; }
        public bool[][] Empty { get; }

        private readonly Units _units;

        public LogicState(LevelObject level)
        {
            _units = BuildUnits(level);
            Globe = Enumerable.Range(0, N).Select(_ => new bool[N]).ToArray();
            Empty = Enumerable.Range(0, N).Select(_ => new bool[N]).ToArray();
        }

        public bool SetEmpty(int id)
        {
            var (r, c) = FromId(id);

            if (Globe[r][c]) return false;
            if (Empty[r][c]) return true;

            Empty[r][c] = true;
            return true;
        }

        public bool SetGlobe(int id)
        {
            var (r, c) = FromId(id);

            if (Empty[r][c]) return false;
            if (Globe[r][c]) return true;

            Globe[r][c] = true;

            foreach (var (dr, dc) in DIR8)
            {
                int nr = r + dr, nc = c + dc;
                if (InBounds(nr, nc))
                {
                    if (!SetEmpty(CellId(nr, nc))) return false;
                }
            }

            return true;
        }

        private static int CountIn(List<int> unit, bool[][] arr)
        {
            int cnt = 0;
            foreach (var id in unit)
            {
                var (r, c) = FromId(id);
                if (arr[r][c]) cnt++;
            }
            return cnt;
        }

        private bool CanPlace(int id)
        {
            var (r, c) = FromId(id);

            if (Empty[r][c]) return false;
            if (Globe[r][c]) return true;

            foreach (var (dr, dc) in DIR8)
            {
                int nr = r + dr, nc = c + dc;
                if (InBounds(nr, nc) && Globe[nr][nc]) return false;
            }

            return true;
        }

        private (bool Ok, int ForcedGlobes) ApplyUnit(List<int> unit, Func<int, int, bool> touchFn)
        {
            int placed = CountIn(unit, Globe);
            int need = K - placed;
            if (need < 0) return (false, 0);

            int forcedGlobes = 0;

            if (need == 0)
            {
                foreach (var id in unit)
                {
                    var (r, c) = FromId(id);
                    if (!Globe[r][c])
                        if (!SetEmpty(id)) return (false, forcedGlobes);
                }
                return (true, forcedGlobes);
            }

            var combos = EnumerateCombos(unit, need, CanPlace, touchFn);
            if (combos.Count == 0) return (false, forcedGlobes);

            var union = new HashSet<int>();
            var inter = new HashSet<int>(combos[0]);

            foreach (var combo in combos)
            {
                foreach (var id in combo) union.Add(id);

                // intersection without modifying during enumeration of inter
                inter.IntersectWith(combo);
            }

            foreach (var id in unit)
            {
                var (r, c) = FromId(id);
                if (Globe[r][c]) continue;

                if (!union.Contains(id))
                    if (!SetEmpty(id)) return (false, forcedGlobes);
            }

            foreach (var id in inter.ToArray())
            {
                var (r, c) = FromId(id);
                if (!Globe[r][c]) forcedGlobes++;
                if (!SetGlobe(id)) return (false, forcedGlobes);
            }

            return (true, forcedGlobes);
        }

        public bool IsSolved()
        {
            foreach (var u in _units.Rows) if (CountIn(u, Globe) != K) return false;
            foreach (var u in _units.Cols) if (CountIn(u, Globe) != K) return false;
            foreach (var u in _units.Regs) if (CountIn(u, Globe) != K) return false;
            return true;
        }

        private static bool TouchRow(int aId, int bId)
        {
            var (ar, ac) = FromId(aId);
            var (br, bc) = FromId(bId);
            return ar == br && Math.Abs(ac - bc) == 1;
        }

        private static bool TouchCol(int aId, int bId)
        {
            var (ar, ac) = FromId(aId);
            var (br, bc) = FromId(bId);
            return ac == bc && Math.Abs(ar - br) == 1;
        }

        private static bool TouchKing(int aId, int bId)
        {
            var (ar, ac) = FromId(aId);
            var (br, bc) = FromId(bId);
            return Math.Abs(ar - br) <= 1 && Math.Abs(ac - bc) <= 1;
        }

        public (bool Ok, int ForcedGlobes) PropagateOnce()
        {
            int forcedGlobes = 0;

            foreach (var u in _units.Rows)
            {
                var r = ApplyUnit(u, TouchRow);
                if (!r.Ok) return (false, forcedGlobes);
                forcedGlobes += r.ForcedGlobes;
            }

            foreach (var u in _units.Cols)
            {
                var r = ApplyUnit(u, TouchCol);
                if (!r.Ok) return (false, forcedGlobes);
                forcedGlobes += r.ForcedGlobes;
            }

            foreach (var u in _units.Regs)
            {
                var r = ApplyUnit(u, TouchKing);
                if (!r.Ok) return (false, forcedGlobes);
                forcedGlobes += r.ForcedGlobes;
            }

            return (true, forcedGlobes);
        }
    }

    private sealed record LogicRes(bool Solved, int ForcedFirst);

    private static LogicRes LogicSolveNoGuess(LevelObject level)
    {
        var st = new LogicState(level);

        int forcedFirst = 0;
        bool changed = true;

        for (int i = 0; i < 300 && changed; i++)
        {
            ulong snap = Snapshot(st);
            var res = st.PropagateOnce();
            if (!res.Ok) return new LogicRes(false, forcedFirst);

            if (i == 0) forcedFirst = res.ForcedGlobes;
            if (st.IsSolved()) return new LogicRes(true, forcedFirst);

            ulong now = Snapshot(st);
            changed = now != snap;
        }

        return new LogicRes(st.IsSolved(), forcedFirst);

        static ulong Snapshot(LogicState s)
        {
            const ulong FnvOffset = 1469598103934665603UL;
            const ulong FnvPrime = 1099511628211UL;

            ulong h = FnvOffset;

            for (int r = 0; r < N; r++)
                for (int c = 0; c < N; c++)
                {
                    h ^= s.Globe[r][c] ? 1UL : 0UL;
                    h *= FnvPrime;
                }

            for (int r = 0; r < N; r++)
                for (int c = 0; c < N; c++)
                {
                    h ^= s.Empty[r][c] ? 1UL : 0UL;
                    h *= FnvPrime;
                }

            return h;
        }
    }
}
