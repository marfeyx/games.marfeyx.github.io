using System.Diagnostics;
using System.Text;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace LevelGen;

internal static class Program
{
    private static readonly JsonSerializerOptions JsonOpts = new()
    {
        WriteIndented = true,
        PropertyNamingPolicy = JsonNamingPolicy.CamelCase,
        DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull
    };

    public static async Task<int> Main()
    {
        Console.OutputEncoding = Encoding.UTF8;

        Console.WriteLine("=== 🌍 Battle Level Generator ===\n");

        string outDir = AskPath("Enter Folder Path: ");
        Directory.CreateDirectory(outDir);

        int amount = AskInt("Enter Level Amount: ", 1, 1_000_000);

        int centerDifficulty = AskInt("Enter Difficulty (1-10): ", 1, 10);
        int tolerance = AskInt("Enter Difficulty Tolerance (0-9) [example: 1 => 4-6]: ", 0, 9);

        int minDifficulty = Utils.Clamp(centerDifficulty - tolerance, 1, 10);
        int maxDifficulty = Utils.Clamp(centerDifficulty + tolerance, 1, 10);

        bool requireEasyStart = AskYesNo("Require Easy Start (forced move at start)? (y/n): ");

        int startLevel = Utils.GetNextLevelNumber(outDir);
        int endLevel = startLevel + amount - 1;

        Console.WriteLine($"\nNext level number: {startLevel}");
        Console.WriteLine($"Generating {amount} levels ({startLevel}..{endLevel})");
        Console.WriteLine($"Difficulty range: {minDifficulty}..{maxDifficulty} (center {centerDifficulty})");
        Console.WriteLine($"Easy start required: {(requireEasyStart ? "YES" : "NO")}\n");

        var slots = Enumerable.Range(0, amount)
            .Select(i => new Slot(index: i, levelNumber: startLevel + i))
            .ToList();

        // -------------------------
        // STAGE 1: Fill all slots with RAW levels
        // -------------------------
        Console.WriteLine("== Stage 1/3: Generate RAW levels (regions + full solution) ==");
        await FillMissingRawAsync(slots, centerDifficulty);

        // -------------------------
        // STAGE 2: Solvability pass until ALL are solvable (and easy-start if enabled)
        // -------------------------
        Console.WriteLine("\n== Stage 2/3: Verify human solvable (no guessing) ==");
        await EnsureAllSolvableAsync(slots, centerDifficulty, requireEasyStart);

        // -------------------------
        // STAGE 3: Difficulty pass until ALL within range (unless user accepts mismatches)
        // -------------------------
        Console.WriteLine("\n== Stage 3/3: Difficulty check ==");
        await EnsureAllDifficultyAsync(slots, centerDifficulty, minDifficulty, maxDifficulty, requireEasyStart);

        // -------------------------
        // SAVE
        // -------------------------
        Console.WriteLine("\n== Saving ==");
        foreach (var s in slots)
        {
            var lvl = s.Level!;

            string jsonName = $"Level {lvl.Level} {centerDifficulty}.json";
            string pngName = $"Level {lvl.Level} {centerDifficulty}.png";

            string jsonPath = Path.Combine(outDir, jsonName);
            string pngPath = Path.Combine(outDir, pngName);

            var json = JsonSerializer.Serialize(lvl, JsonOpts);
            await File.WriteAllTextAsync(jsonPath, json, Encoding.UTF8);

            SolutionPngRenderer.SaveSolutionPng(lvl, pngPath, cellSize: 32, padding: 16);

            Console.WriteLine($"Saved: {jsonName} ✅ (estDiff={lvl.Difficulty})");
            Console.WriteLine($"Saved: {pngName} ✅");
        }

        Console.WriteLine("\nAll done. 🎯");
        return 0;
    }

    private static async Task FillMissingRawAsync(List<Slot> slots, int centerDifficulty)
    {
        var sw = Stopwatch.StartNew();
        long lastPrint = -9999;

        while (slots.Any(s => s.Level is null))
        {
            for (int i = 0; i < slots.Count; i++)
            {
                var s = slots[i];
                if (s.Level is not null) continue;

                var raw = await GeneratorCore.GenerateRawLevelAsync(
                    requestedDifficulty: centerDifficulty,
                    levelNumber: s.LevelNumber,
                    progress: msg => ThrottledPrint(sw, ref lastPrint, msg)
                );

                s.Level = raw;
                s.Solvability = null;

                Console.WriteLine($"  ✅ Ready: slot {s.Index + 1}/{slots.Count} -> L{s.LevelNumber}");
            }
        }

        static void ThrottledPrint(Stopwatch sw, ref long lastPrint, string msg)
        {
            long ms = sw.ElapsedMilliseconds;
            if (ms - lastPrint < 500) return;
            lastPrint = ms;
            Console.WriteLine("  " + msg);
        }
    }

    private static async Task EnsureAllSolvableAsync(List<Slot> slots, int centerDifficulty, bool requireEasyStart)
    {
        int pass = 0;

        while (true)
        {
            pass++;
            Console.WriteLine($"\nSolvability pass #{pass}");

            var sw = Stopwatch.StartNew();
            long lastPrint = -9999;

            int checkedCount = 0;
            int okCount = 0;
            int failLogic = 0;
            int failEasy = 0;

            var failedIndices = new List<int>();

            for (int i = 0; i < slots.Count; i++)
            {
                var s = slots[i];
                var lvl = s.Level!;
                checkedCount++;

                var info = GeneratorCore.CheckHumanSolvable(lvl, requireEasyStart);
                s.Solvability = info;

                if (info.Solved && (!requireEasyStart || info.EasyStart))
                {
                    okCount++;
                }
                else
                {
                    failedIndices.Add(i);
                    s.Level = null;
                    s.Solvability = null;

                    if (info.FailReason == "logic") failLogic++;
                    else if (info.FailReason == "easyStart") failEasy++;
                }

                ThrottledPrint(sw, ref lastPrint,
                    $"[SOLVE] {checkedCount}/{slots.Count} ok={okCount} lastFailReason={info.FailReason}");
                await Task.Yield();
            }

            Console.WriteLine($"Pass #{pass} result: ok={okCount}/{slots.Count}, regen={failedIndices.Count}");

            if (failedIndices.Count == 0)
            {
                Console.WriteLine("✅ All levels are solvable (no guessing).");
                return;
            }

            Console.WriteLine("Regenerating failed slots...");
            await FillMissingRawAsync(slots, centerDifficulty);
        }

        static void ThrottledPrint(Stopwatch sw, ref long lastPrint, string msg)
        {
            long ms = sw.ElapsedMilliseconds;
            if (ms - lastPrint < 500) return;
            lastPrint = ms;
            Console.WriteLine("  " + msg);
        }
    }

    private static async Task EnsureAllDifficultyAsync(
        List<Slot> slots,
        int centerDifficulty,
        int minDifficulty,
        int maxDifficulty,
        bool requireEasyStart)
    {
        int pass = 0;

        while (true)
        {
            pass++;
            Console.WriteLine($"\n-- Difficulty pass #{pass} (need {minDifficulty}..{maxDifficulty}) --");

            var sw = Stopwatch.StartNew();
            long lastPrint = -9999;

            int checkedCount = 0;
            int inRange = 0;
            int outRange = 0;

            var mismatchingIndices = new List<int>();

            for (int i = 0; i < slots.Count; i++)
            {
                var s = slots[i];
                var lvl = s.Level!;
                var sol = s.Solvability!;

                checkedCount++;

                // We estimate from brute-force nodes gathered in solvability check
                int est = GeneratorCore.DifficultyFromNodes(sol.Nodes);
                lvl.Difficulty = est;

                if (est >= minDifficulty && est <= maxDifficulty)
                {
                    inRange++;
                }
                else
                {
                    outRange++;
                    mismatchingIndices.Add(i);
                }

                ThrottledPrint(sw, ref lastPrint,
                    $"[DIFF] {checkedCount}/{slots.Count} inRange={inRange} outRange={outRange} lastEst={est} need={minDifficulty}..{maxDifficulty}");
                await Task.Yield();
            }

            Console.WriteLine($"Pass #{pass} result: inRange={inRange}/{slots.Count}, outRange={outRange}");

            if (mismatchingIndices.Count == 0)
            {
                Console.WriteLine("✅ All levels within difficulty range.");
                return;
            }

            Console.Write($"Found {mismatchingIndices.Count} out-of-range levels. Save anyway and stop hunting? (y/n): ");
            var ans = (Console.ReadLine() ?? "").Trim().ToLowerInvariant();

            if (ans is "y" or "yes")
            {
                Console.WriteLine("✅ Accepting difficulty mismatches. Saving what we have.");
                return;
            }

            Console.WriteLine("Regenerating out-of-range slots...");
            foreach (var idx in mismatchingIndices)
            {
                slots[idx].Level = null;
                slots[idx].Solvability = null;
            }

            await FillMissingRawAsync(slots, centerDifficulty);

            Console.WriteLine("\nRe-checking solvability after replacements...");
            await EnsureAllSolvableAsync(slots, centerDifficulty, requireEasyStart);
        }

        static void ThrottledPrint(Stopwatch sw, ref long lastPrint, string msg)
        {
            long ms = sw.ElapsedMilliseconds;
            if (ms - lastPrint < 500) return;
            lastPrint = ms;
            Console.WriteLine("  " + msg);
        }
    }

    // -------------------------
    // Console input helpers
    // -------------------------

    private static string AskPath(string prompt)
    {
        while (true)
        {
            Console.Write(prompt);
            var s = (Console.ReadLine() ?? "").Trim().Trim('"');
            if (!string.IsNullOrWhiteSpace(s))
                return Path.GetFullPath(s);

            Console.WriteLine("Please enter a valid folder path.");
        }
    }

    private static int AskInt(string prompt, int min, int max)
    {
        while (true)
        {
            Console.Write(prompt);
            var s = (Console.ReadLine() ?? "").Trim();
            if (int.TryParse(s, out int n) && n >= min && n <= max)
                return n;

            Console.WriteLine($"Please enter an integer in [{min}..{max}].");
        }
    }

    private static bool AskYesNo(string prompt)
    {
        while (true)
        {
            Console.Write(prompt);
            var s = (Console.ReadLine() ?? "").Trim().ToLowerInvariant();
            if (s is "y" or "yes") return true;
            if (s is "n" or "no") return false;
            Console.WriteLine("Please enter y or n.");
        }
    }

    private sealed class Slot
    {
        public int Index { get; }
        public int LevelNumber { get; }
        public LevelObject? Level { get; set; }
        public SolvabilityInfo? Solvability { get; set; }

        public Slot(int index, int levelNumber)
        {
            Index = index;
            LevelNumber = levelNumber;
        }
    }
}
