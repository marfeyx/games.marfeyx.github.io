using System.Text.RegularExpressions;

namespace LevelGen;

public static class Utils
{
    private static readonly Regex ReNew = new(@"^Level\s+(\d+)\s+(\d+)\.json$", RegexOptions.IgnoreCase | RegexOptions.Compiled);
    private static readonly Regex ReOld = new(@"^lvl(\d+)\.json$", RegexOptions.IgnoreCase | RegexOptions.Compiled);

    public static int GetNextLevelNumber(string dirPath)
    {
        if (!Directory.Exists(dirPath)) return 1;

        int maxLevel = 0;

        foreach (var file in Directory.EnumerateFiles(dirPath, "*.json", SearchOption.TopDirectoryOnly))
        {
            var name = Path.GetFileName(file);

            var m1 = ReNew.Match(name);
            if (m1.Success && int.TryParse(m1.Groups[1].Value, out int n1))
            {
                maxLevel = Math.Max(maxLevel, n1);
                continue;
            }

            var m2 = ReOld.Match(name);
            if (m2.Success && int.TryParse(m2.Groups[1].Value, out int n2))
            {
                maxLevel = Math.Max(maxLevel, n2);
            }
        }

        return maxLevel + 1;
    }

    public static int Clamp(int n, int a, int b) => Math.Max(a, Math.Min(b, n));
    public static double Clamp(double v, double a, double b) => Math.Max(a, Math.Min(b, v));
}
